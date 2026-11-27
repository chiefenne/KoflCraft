//! Edits of a block topology file (`block.*` edit tools, BS-04). Each edit
//! takes a file and returns a new one with what changed; the new file is
//! checked before it is returned, so an edit that breaks the topology is
//! refused and nothing is written. Blocks and edges are named by their
//! vertex indices, as in the file.

use std::collections::BTreeSet;

use mesh_algorithms::block as algo;
use mesh_core::block::{BLOCK_FACES, block_edges, corner_slot, face_corner_slot, face_key};
use mesh_core::geometry::Vec3;
use mesh_io::handoff::PatchRole;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::block::{
    BlockFile, CurveFile, DistributionFile, EdgeFile, PatchFile, TopologyFile, check,
};
use crate::{EngineError, ErrorCode};

/// One edit of a topology.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(tag = "op", rename_all = "snake_case", deny_unknown_fields)]
pub enum Edit {
    /// Set the cell count of the chain of opposite edges through `edge`
    /// (BS-05): every block edge in the chain gets `cells`.
    SetCount { edge: [usize; 2], cells: usize },
    /// Set the law of one edge, read from `edge[0]` to `edge[1]`.
    SetDistribution {
        edge: [usize; 2],
        distribution: DistributionFile,
    },
    /// Set the curve of one edge; spline points run from `edge[0]` to
    /// `edge[1]`.
    SetCurve { edge: [usize; 2], curve: CurveFile },
    /// Move a vertex to `to`, or by `by`.
    MoveVertex {
        vertex: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        to: Option<[f64; 3]>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        by: Option<[f64; 3]>,
    },
    /// Append a vertex; its index is the old vertex count.
    AddVertex { point: [f64; 3] },
    /// Remove a vertex that nothing uses; higher indices move down by one.
    RemoveVertex { vertex: usize },
    /// Append a block. Its faces that meet another block's face become
    /// shared (and leave their patches); its other faces go to `patch`.
    AddBlock {
        vertices: [usize; 8],
        cells: [usize; 3],
        #[serde(default, skip_serializing_if = "Option::is_none")]
        domain: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        patch: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        role: Option<PatchRole>,
    },
    /// Remove a block. Faces it shared become boundary faces of `patch`;
    /// its boundary faces leave their patches and surfaces.
    RemoveBlock {
        block: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        patch: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        role: Option<PatchRole>,
    },
    /// Split a block along `axis` after cell `at`; every block of the
    /// chain through that axis is split with it, so faces stay
    /// point-matched (BS-05). New vertices sit on the TFI grid nodes.
    SplitBlock {
        block: usize,
        axis: usize,
        at: usize,
    },
    /// Split a block whose cross-section is a triangle (its corner at
    /// `vertex` at least 120°, nearly flat) into three blocks around a
    /// centre vertex; a new row of `cells` cells (from the triangle's
    /// spacing when absent) runs from two of its sides through the
    /// neighbouring blocks to the boundary (BS-03).
    SplitTriangle {
        block: usize,
        vertex: usize,
        /// The other end of the flat edge (3D); the flattest edge at
        /// `vertex` when absent.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        along: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        cells: Option<usize>,
    },
    /// Rename a patch and/or set its role.
    SetPatch {
        name: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        rename: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        role: Option<PatchRole>,
    },
}

/// What an edit changed.
#[derive(Clone, Debug, Default, Serialize, PartialEq)]
pub struct Change {
    pub blocks: Vec<usize>,
    pub edges: Vec<[usize; 2]>,
    pub vertices: Vec<usize>,
    pub patches: Vec<String>,
    /// The template's analytic volume no longer describes the shape after a
    /// geometric edit, so the file's `reference` is removed.
    pub reference_removed: bool,
    /// Blocks and vertices appended by the edit (indices in the new file).
    pub added_blocks: Vec<usize>,
    pub added_vertices: Vec<usize>,
    /// Patches or surfaces left without faces, removed.
    pub removed_patches: Vec<String>,
    /// Edges whose law or curve was fitted to the grid nodes (split).
    pub fitted_edges: Vec<[usize; 2]>,
}

pub(crate) fn refused(message: impl std::fmt::Display) -> EngineError {
    EngineError::new(ErrorCode::InvalidArguments, message)
}

pub(crate) fn key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

/// The edges of block `b` along each axis, as vertex keys.
fn axis_edges(file: &TopologyFile, b: usize) -> [Vec<(usize, usize)>; 3] {
    let v = &file.blocks[b].vertices;
    let mut out: [Vec<(usize, usize)>; 3] = Default::default();
    for (axis, c) in block_edges() {
        let mut e = c;
        e[axis] = 1;
        let (p, q) = (
            v[corner_slot(c[0], c[1], c[2])],
            v[corner_slot(e[0], e[1], e[2])],
        );
        out[axis].push(key(p, q));
    }
    out
}

/// The blocks and axes that carry `edge`; refused when it is no block edge.
fn carriers(file: &TopologyFile, edge: [usize; 2]) -> Result<Vec<(usize, usize)>, EngineError> {
    let k = key(edge[0], edge[1]);
    let found: Vec<(usize, usize)> = (0..file.blocks.len())
        .flat_map(|b| {
            let edges = axis_edges(file, b);
            (0..3)
                .filter(move |&a| edges[a].contains(&k))
                .map(move |a| (b, a))
                .collect::<Vec<_>>()
        })
        .collect();
    if found.is_empty() {
        return Err(refused(format!(
            "edge {}–{} is not an edge of any block",
            edge[0], edge[1]
        )));
    }
    Ok(found)
}

/// The same law read from the other end.
pub(crate) fn reversed_law(d: &DistributionFile) -> DistributionFile {
    match d {
        DistributionFile::Uniform => DistributionFile::Uniform,
        DistributionFile::Geometric { first, last } => DistributionFile::Geometric {
            first: *last,
            last: *first,
        },
        DistributionFile::Tanh { first, last } => DistributionFile::Tanh {
            first: *last,
            last: *first,
        },
    }
}

pub(crate) fn reversed_curve(c: &CurveFile) -> CurveFile {
    match c {
        CurveFile::Spline { through } => CurveFile::Spline {
            through: through.iter().rev().copied().collect(),
        },
        other => other.clone(),
    }
}

/// The file's entry for `edge`, created as a uniform line when absent,
/// and whether it runs the other way.
fn edge_entry(file: &mut TopologyFile, edge: [usize; 2]) -> (&mut EdgeFile, bool) {
    let k = key(edge[0], edge[1]);
    let i = match file
        .edges
        .iter()
        .position(|e| key(e.vertices[0], e.vertices[1]) == k)
    {
        Some(i) => i,
        None => {
            file.edges.push(EdgeFile {
                vertices: edge,
                curve: CurveFile::Line,
                distribution: DistributionFile::Uniform,
            });
            file.edges.len() - 1
        }
    };
    let e = &mut file.edges[i];
    let flipped = e.vertices != edge;
    (e, flipped)
}

/// The four vertices of block face `f` in local order (00, 10, 11, 01).
pub(crate) fn face_vertices(v: &[usize; 8], f: usize) -> [usize; 4] {
    let face = BLOCK_FACES[f];
    [(0, 0), (1, 0), (1, 1), (0, 1)].map(|(a, b)| v[face_corner_slot(face, a, b)])
}

/// How many blocks have a face with this key.
fn face_users(file: &TopologyFile, k: [usize; 4]) -> usize {
    file.blocks
        .iter()
        .map(|b| {
            (0..6)
                .filter(|&f| face_key(face_vertices(&b.vertices, f)) == k)
                .count()
        })
        .sum()
}

/// Take faces with these keys out of every patch and surface.
fn unlist_faces(file: &mut TopologyFile, keys: &BTreeSet<[usize; 4]>) {
    for p in &mut file.patches {
        p.faces.retain(|f| !keys.contains(&face_key(*f)));
    }
    for s in &mut file.surfaces {
        s.faces.retain(|f| !keys.contains(&face_key(*f)));
    }
}

/// Patches and surfaces left without faces are removed.
fn drop_empty(file: &mut TopologyFile, change: &mut Change) {
    for p in file.patches.iter().filter(|p| p.faces.is_empty()) {
        change.removed_patches.push(p.name.clone());
    }
    for s in file.surfaces.iter().filter(|s| s.faces.is_empty()) {
        change.removed_patches.push(format!("surface {}", s.name));
    }
    file.patches.retain(|p| !p.faces.is_empty());
    file.surfaces.retain(|s| !s.faces.is_empty());
}

/// Add faces to the patch `name`, created with `role` when absent.
fn list_faces(
    file: &mut TopologyFile,
    name: &str,
    role: Option<PatchRole>,
    faces: Vec<[usize; 4]>,
    change: &mut Change,
) {
    if faces.is_empty() {
        return;
    }
    match file.patches.iter_mut().find(|p| p.name == name) {
        Some(p) => {
            p.faces.extend(faces);
            if role.is_some() {
                p.role = role;
            }
        }
        None => file.patches.push(PatchFile {
            name: name.into(),
            role,
            faces,
        }),
    }
    change.patches.push(name.into());
}

/// Edge entries whose vertex pair is no longer a block edge are removed.
fn drop_loose_edges(file: &mut TopologyFile) {
    let keys: BTreeSet<(usize, usize)> = (0..file.blocks.len())
        .flat_map(|b| axis_edges(file, b).into_iter().flatten())
        .collect();
    file.edges
        .retain(|e| keys.contains(&key(e.vertices[0], e.vertices[1])));
}

fn remove_vertex(
    file: &mut TopologyFile,
    v: usize,
    change: &mut Change,
) -> Result<(), EngineError> {
    if v >= file.vertices.len() {
        return Err(refused(format!("vertex {v} does not exist")));
    }
    let used = file.blocks.iter().any(|b| b.vertices.contains(&v))
        || file.edges.iter().any(|e| e.vertices.contains(&v))
        || file
            .patches
            .iter()
            .any(|p| p.faces.iter().any(|f| f.contains(&v)))
        || file
            .surfaces
            .iter()
            .any(|s| s.faces.iter().any(|f| f.contains(&v)));
    if used {
        return Err(refused(format!(
            "vertex {v} is still used by a block, edge, patch or surface; remove those first"
        )));
    }
    file.vertices.remove(v);
    let shift = |i: &mut usize| {
        if *i > v {
            *i -= 1;
        }
    };
    for b in &mut file.blocks {
        b.vertices.iter_mut().for_each(shift);
    }
    for e in &mut file.edges {
        e.vertices.iter_mut().for_each(shift);
    }
    for p in &mut file.patches {
        p.faces.iter_mut().flatten().for_each(shift);
    }
    for s in &mut file.surfaces {
        s.faces.iter_mut().flatten().for_each(shift);
    }
    change.vertices.push(v);
    Ok(())
}

fn add_block(
    file: &mut TopologyFile,
    block: BlockFile,
    patch: Option<&str>,
    role: Option<PatchRole>,
    change: &mut Change,
) -> Result<(), EngineError> {
    let mut shared = BTreeSet::new();
    let mut open = Vec::new();
    for f in 0..6 {
        let fv = face_vertices(&block.vertices, f);
        if face_users(file, face_key(fv)) > 0 {
            shared.insert(face_key(fv));
        } else {
            open.push(fv);
        }
    }
    if !open.is_empty() && patch.is_none() {
        return Err(refused(format!(
            "the new block has {} boundary faces; name the `patch` they go to",
            open.len()
        )));
    }
    file.blocks.push(block);
    change.added_blocks.push(file.blocks.len() - 1);
    unlist_faces(file, &shared);
    if let Some(name) = patch {
        list_faces(file, name, role, open, change);
    }
    drop_empty(file, change);
    Ok(())
}

fn remove_block(
    file: &mut TopologyFile,
    b: usize,
    patch: Option<&str>,
    role: Option<PatchRole>,
    change: &mut Change,
) -> Result<(), EngineError> {
    if b >= file.blocks.len() {
        return Err(refused(format!("block {b} does not exist")));
    }
    if file.blocks.len() == 1 {
        return Err(refused("the last block cannot be removed"));
    }
    let removed = file.blocks.remove(b);
    let mut boundary = BTreeSet::new();
    let mut opened = Vec::new();
    for f in 0..6 {
        let fv = face_vertices(&removed.vertices, f);
        let k = face_key(fv);
        if face_users(file, k) > 0 {
            // the neighbour's face, which is wound for the neighbour
            let n = file
                .blocks
                .iter()
                .flat_map(|o| (0..6).map(move |g| face_vertices(&o.vertices, g)))
                .find(|g| face_key(*g) == k)
                .unwrap_or(fv);
            opened.push(n);
        } else {
            boundary.insert(k);
        }
    }
    if !opened.is_empty() && patch.is_none() {
        return Err(refused(format!(
            "removing block {b} opens {} faces of its neighbours; name the `patch` they go to",
            opened.len()
        )));
    }
    unlist_faces(file, &boundary);
    if let Some(name) = patch {
        list_faces(file, name, role, opened, change);
    }
    drop_loose_edges(file);
    drop_empty(file, change);
    change.blocks.push(b);
    Ok(())
}

pub(crate) fn a3(p: Vec3) -> [f64; 3] {
    [p.x, p.y, p.z]
}

/// A law for nodes along an edge: uniform when the spacings agree,
/// otherwise the `like` law's kind with the spacings at the ends.
pub(crate) fn fitted_law(
    nodes: &[Vec3],
    curve: &CurveFile,
    like: &DistributionFile,
) -> DistributionFile {
    let s = arc_spacings(nodes, curve);
    let (lo, hi) = s
        .iter()
        .fold((f64::INFINITY, 0.0f64), |(a, b), &x| (a.min(x), b.max(x)));
    if s.len() < 2 || hi - lo <= 1e-9 * hi {
        return DistributionFile::Uniform;
    }
    let (first, last) = (s[0], s[s.len() - 1]);
    if s.len() == 2 {
        // two cells: the first spacing fixes both
        return DistributionFile::Geometric {
            first: Some(first),
            last: None,
        };
    }
    match like {
        DistributionFile::Geometric { first: Some(_), .. } => DistributionFile::Geometric {
            first: Some(first),
            last: None,
        },
        DistributionFile::Geometric { last: Some(_), .. } => DistributionFile::Geometric {
            first: None,
            last: Some(last),
        },
        _ => {
            // Vinokur's end spacings are slopes, not cell lengths: correct
            // them until the cells at the ends have the measured lengths
            let total: f64 = s.iter().sum();
            let n = s.len();
            let (mut f, mut l) = (first, last);
            for _ in 0..100 {
                let Ok(x) = algo::distribution::fractions(
                    mesh_core::block::Distribution::Tanh { first: f, last: l },
                    n,
                    total,
                ) else {
                    break;
                };
                let (af, al) = (x[1] * total, (1.0 - x[n - 1]) * total);
                let (rf, rl) = (first / af, last / al);
                if (rf - 1.0).abs() < 1e-12 && (rl - 1.0).abs() < 1e-12 {
                    break;
                }
                f *= rf;
                l *= rl;
            }
            DistributionFile::Tanh { first: f, last: l }
        }
    }
}

/// A curve through the nodes: a line when they are straight, an arc when
/// they lie on one circle, else a spline through the inner nodes.
pub(crate) fn fitted_curve(nodes: &[Vec3]) -> CurveFile {
    let (a, b) = (nodes[0], nodes[nodes.len() - 1]);
    let len = (b - a).length();
    let inner = &nodes[1..nodes.len() - 1];
    let off_line = inner
        .iter()
        .map(|&p| {
            let t = (p - a).dot(b - a) / (len * len);
            (p - (a + (b - a) * t)).length()
        })
        .fold(0.0f64, f64::max);
    if inner.is_empty() || off_line <= 1e-12 * len {
        return CurveFile::Line;
    }
    let mid = inner[inner.len() / 2];
    if let Some((centre, radius)) = circle(a, mid, b) {
        let on = inner
            .iter()
            .all(|&p| ((p - centre).length() - radius).abs() <= 1e-9 * radius);
        let w = (mid - a).cross(b - a);
        let plane = w / w.length();
        let flat = inner
            .iter()
            .all(|&p| (p - a).dot(plane).abs() <= 1e-9 * len);
        if on && flat {
            return CurveFile::Arc { through: a3(mid) };
        }
    }
    CurveFile::Spline {
        through: inner.iter().map(|&p| a3(p)).collect(),
    }
}

pub(crate) fn v3(p: [f64; 3]) -> Vec3 {
    Vec3::new(p[0], p[1], p[2])
}

/// The through points of a spline from `a` to `b` on either side of the
/// point `m` on it: the polyline segment nearest to `m` decides.
pub(crate) fn cut_spline(a: Vec3, through: &[Vec3], b: Vec3, m: Vec3) -> (Vec<Vec3>, Vec<Vec3>) {
    let mut poly = vec![a];
    poly.extend_from_slice(through);
    poly.push(b);
    let mut best = (f64::INFINITY, 0);
    for (i, w) in poly.windows(2).enumerate() {
        let d = w[1] - w[0];
        let t = ((m - w[0]).dot(d) / d.dot(d)).clamp(0.0, 1.0);
        let dist = (m - (w[0] + d * t)).length();
        if dist < best.0 {
            best = (dist, i);
        }
    }
    // segment i runs from poly[i] to poly[i + 1]: through[..i] before m
    let i = best.1;
    let len = (b - a).length();
    let apart = |p: &Vec3| (*p - m).length() > 1e-9 * len;
    let first = through[..i].iter().copied().filter(apart).collect();
    let second = through[i..].iter().copied().filter(apart).collect();
    (first, second)
}

/// The spacings of the nodes measured along `curve` from the first node to
/// the last, which is how a law places them (chords where it fails).
fn arc_spacings(nodes: &[Vec3], curve: &CurveFile) -> Vec<f64> {
    let chords = || nodes.windows(2).map(|w| (w[1] - w[0]).length()).collect();
    let c = match curve {
        CurveFile::Line => return chords(),
        CurveFile::Arc { through } => mesh_core::block::Curve::Arc {
            through: v3(*through),
        },
        CurveFile::Spline { through } => mesh_core::block::Curve::Spline {
            through: through.iter().map(|&p| v3(p)).collect(),
        },
    };
    let (a, b) = (nodes[0], nodes[nodes.len() - 1]);
    let Ok(ev) = algo::distribution::Evaluator::new(a, b, &c) else {
        return chords();
    };
    const SAMPLES: usize = 8192;
    let pts: Vec<Vec3> = (0..=SAMPLES)
        .map(|i| ev.at(i as f64 / SAMPLES as f64))
        .collect();
    // fraction of each node: projection onto the nearest sample segment
    let mut start = 0;
    let mut fr = Vec::with_capacity(nodes.len());
    for (k, &p) in nodes.iter().enumerate() {
        if k == 0 {
            fr.push(0.0);
            continue;
        }
        if k == nodes.len() - 1 {
            fr.push(1.0);
            continue;
        }
        let mut best = (f64::INFINITY, 0.0, start);
        for (i, w) in pts.windows(2).enumerate().skip(start) {
            let d = w[1] - w[0];
            let t = ((p - w[0]).dot(d) / d.dot(d)).clamp(0.0, 1.0);
            let dist = (p - (w[0] + d * t)).length();
            if dist < best.0 {
                best = (dist, (i as f64 + t) / SAMPLES as f64, i);
            }
        }
        start = best.2;
        fr.push(best.1);
    }
    let l = ev.length();
    fr.windows(2).map(|w| (w[1] - w[0]) * l).collect()
}

/// The circle through three points.
fn circle(a: Vec3, b: Vec3, c: Vec3) -> Option<(Vec3, f64)> {
    let (u, v) = (b - a, c - a);
    let w = u.cross(v);
    let d = 2.0 * w.dot(w);
    if d <= f64::EPSILON * u.dot(u) * v.dot(v) {
        return None;
    }
    let centre = a + (w.cross(u) * v.dot(v) + v.cross(w) * u.dot(u)) / d;
    Some((centre, (a - centre).length()))
}

fn split_block(
    file: &mut TopologyFile,
    b0: usize,
    axis0: usize,
    at: usize,
    change: &mut Change,
) -> Result<(), EngineError> {
    if b0 >= file.blocks.len() || axis0 > 2 {
        return Err(refused(format!(
            "block {b0} axis {axis0} does not exist (axes are 0, 1, 2)"
        )));
    }
    let n0 = file.blocks[b0].cells[axis0];
    if at == 0 || at >= n0 {
        return Err(refused(format!(
            "split after cell {at}: block {b0} has {n0} cells along axis {axis0}, so `at` must be 1 to {}",
            n0.saturating_sub(1)
        )));
    }
    let t = file.topology()?;
    let grid = algo::build(&t).map_err(|issues| {
        refused(
            issues
                .iter()
                .map(|i| i.message.clone())
                .collect::<Vec<_>>()
                .join("; "),
        )
    })?;
    // oriented split edges: key → (start vertex, split index from start)
    let mut split: std::collections::BTreeMap<(usize, usize), (usize, usize)> =
        std::collections::BTreeMap::new();
    // chain blocks: block → (axis, split index from the local start)
    let mut chain: std::collections::BTreeMap<usize, (usize, usize)> =
        std::collections::BTreeMap::new();
    let ends = |file: &TopologyFile, b: usize, a: usize| -> Vec<(usize, usize)> {
        let v = &file.blocks[b].vertices;
        block_edges()
            .into_iter()
            .filter(|(ax, _)| *ax == a)
            .map(|(_, c)| {
                let mut e = c;
                e[a] = 1;
                (
                    v[corner_slot(c[0], c[1], c[2])],
                    v[corner_slot(e[0], e[1], e[2])],
                )
            })
            .collect()
    };
    let mut queue = vec![(b0, axis0, at)];
    while let Some((b, a, k)) = queue.pop() {
        if let Some(&(a2, k2)) = chain.get(&b) {
            if a2 != a || k2 != k {
                return Err(refused(format!(
                    "the split reaches block {b} twice in different ways; split it by hand"
                )));
            }
            continue;
        }
        chain.insert(b, (a, k));
        let n = file.blocks[b].cells[a];
        for (p, q) in ends(file, b, a) {
            let mine = (p, k);
            match split.get(&key(p, q)) {
                Some(&(s, ks)) => {
                    let same = if s == p { ks == k } else { ks == n - k };
                    if !same {
                        return Err(refused(format!(
                            "edge {p}–{q} would be split at two places; split by hand"
                        )));
                    }
                }
                None => {
                    split.insert(key(p, q), mine);
                }
            }
        }
        // neighbours that share one of these edges along some axis
        for o in 0..file.blocks.len() {
            if chain.contains_key(&o) {
                continue;
            }
            for a2 in 0..3 {
                for (p, q) in ends(file, o, a2) {
                    if let Some(&(s, ks)) = split.get(&key(p, q)) {
                        let n2 = file.blocks[o].cells[a2];
                        let k2 = if s == p { ks } else { n2 - ks };
                        queue.push((o, a2, k2));
                    }
                }
            }
        }
    }
    // the nodes along a grid line of block b: axis a, from local corner c
    let line = |b: usize, a: usize, c: [usize; 3]| -> Vec<Vec3> {
        let n = file.blocks[b].cells;
        let mut idx = [c[0] * n[0], c[1] * n[1], c[2] * n[2]];
        (0..=n[a])
            .map(|s| {
                idx[a] = s;
                grid.points[grid.block_corner(b, idx[0], idx[1], idx[2])]
            })
            .collect()
    };
    // new vertices and the two halves of each split edge
    let mut new_vertex: std::collections::BTreeMap<(usize, usize), usize> =
        std::collections::BTreeMap::new();
    let mut new_edges: Vec<EdgeFile> = Vec::new();
    let old_edges = file.edges.clone();
    for (&b, &(a, _)) in &chain {
        for (_, c) in block_edges().into_iter().filter(|(ax, _)| *ax == a) {
            let v = &file.blocks[b].vertices;
            let mut e = c;
            e[a] = 1;
            let (p, q) = (
                v[corner_slot(c[0], c[1], c[2])],
                v[corner_slot(e[0], e[1], e[2])],
            );
            if new_vertex.contains_key(&key(p, q)) {
                continue;
            }
            let (s, ks) = split[&key(p, q)];
            let mut nodes = line(b, a, c);
            if s != p {
                nodes.reverse();
            }
            let other = if s == p { q } else { p };
            file.vertices.push(a3(nodes[ks]));
            let m = file.vertices.len() - 1;
            new_vertex.insert(key(p, q), m);
            change.added_vertices.push(m);
            let old = old_edges
                .iter()
                .find(|x| key(x.vertices[0], x.vertices[1]) == key(p, q));
            if let Some(old) = old {
                let like = if old.vertices[0] == s {
                    old.distribution.clone()
                } else {
                    reversed_law(&old.distribution)
                };
                // a spline keeps its own points, cut at the new vertex
                let halves = match &old.curve {
                    CurveFile::Spline { through } => {
                        let mut pts: Vec<Vec3> = through
                            .iter()
                            .map(|&p| Vec3::new(p[0], p[1], p[2]))
                            .collect();
                        if old.vertices[0] != s {
                            pts.reverse();
                        }
                        Some(cut_spline(
                            v3(file.vertices[s]),
                            &pts,
                            v3(file.vertices[other]),
                            nodes[ks],
                        ))
                    }
                    _ => None,
                };
                for (h, (from, to, part)) in [(s, m, &nodes[..=ks]), (m, other, &nodes[ks..])]
                    .into_iter()
                    .enumerate()
                {
                    let curve = match (&old.curve, &halves) {
                        (CurveFile::Line, _) => CurveFile::Line,
                        (_, Some((a, b))) => {
                            let pts = if h == 0 { a } else { b };
                            if pts.is_empty() {
                                fitted_curve(part)
                            } else {
                                CurveFile::Spline {
                                    through: pts.iter().map(|&p| a3(p)).collect(),
                                }
                            }
                        }
                        _ => fitted_curve(part),
                    };
                    let distribution = if like == DistributionFile::Uniform {
                        DistributionFile::Uniform
                    } else {
                        fitted_law(part, &curve, &like)
                    };
                    if curve != CurveFile::Line || distribution != DistributionFile::Uniform {
                        new_edges.push(EdgeFile {
                            vertices: [from, to],
                            curve,
                            distribution,
                        });
                        change.fitted_edges.push([from, to]);
                    }
                }
            }
        }
    }
    file.edges
        .retain(|x| !new_vertex.contains_key(&key(x.vertices[0], x.vertices[1])));
    // the blocks: the low half keeps the index, the high half is appended
    let mut faces_split: std::collections::BTreeMap<[usize; 4], [[usize; 4]; 2]> =
        std::collections::BTreeMap::new();
    let mut cut_edges = BTreeSet::new();
    for (&b, &(a, k)) in &chain {
        let old = file.blocks[b].clone();
        let mut low = old.clone();
        let mut high = old.clone();
        for (_, c) in block_edges().into_iter().filter(|(ax, _)| *ax == a) {
            let mut e = c;
            e[a] = 1;
            let (sa, sb) = (corner_slot(c[0], c[1], c[2]), corner_slot(e[0], e[1], e[2]));
            let m = new_vertex[&key(old.vertices[sa], old.vertices[sb])];
            low.vertices[sb] = m;
            high.vertices[sa] = m;
        }
        low.cells[a] = k;
        high.cells[a] = old.cells[a] - k;
        for (f, face) in BLOCK_FACES.iter().enumerate() {
            if face.0 != a {
                faces_split.insert(
                    face_key(face_vertices(&old.vertices, f)),
                    [
                        face_vertices(&low.vertices, f),
                        face_vertices(&high.vertices, f),
                    ],
                );
            }
        }
        // the edges of the cut face, fitted to the grid line at index k
        for (a2, c) in block_edges() {
            if a2 == a || c[a] != 0 {
                continue;
            }
            let mut e = c;
            e[a2] = 1;
            let (sa, sb) = (corner_slot(c[0], c[1], c[2]), corner_slot(e[0], e[1], e[2]));
            // on the cut face: the high block's low side
            let (p, q) = (high.vertices[sa], high.vertices[sb]);
            if !cut_edges.insert(key(p, q)) {
                continue;
            }
            let n = old.cells;
            let mut idx = [c[0] * n[0], c[1] * n[1], c[2] * n[2]];
            idx[a] = k;
            let nodes: Vec<Vec3> = (0..=n[a2])
                .map(|s| {
                    idx[a2] = s;
                    grid.points[grid.block_corner(b, idx[0], idx[1], idx[2])]
                })
                .collect();
            let curve = fitted_curve(&nodes);
            // the TFI line inside a block follows no law exactly; a tanh
            // law with the measured end cells comes closest
            let distribution = fitted_law(
                &nodes,
                &curve,
                &DistributionFile::Tanh {
                    first: 0.0,
                    last: 0.0,
                },
            );
            if curve != CurveFile::Line || distribution != DistributionFile::Uniform {
                new_edges.push(EdgeFile {
                    vertices: [p, q],
                    curve,
                    distribution,
                });
                change.fitted_edges.push([p, q]);
            }
        }
        file.blocks[b] = low;
        file.blocks.push(high);
        change.blocks.push(b);
        change.added_blocks.push(file.blocks.len() - 1);
    }
    file.edges.extend(new_edges);
    let replace = |faces: &mut Vec<[usize; 4]>| {
        *faces = faces
            .iter()
            .flat_map(|f| match faces_split.get(&face_key(*f)) {
                Some(two) => two.to_vec(),
                None => vec![*f],
            })
            .collect();
    };
    for p in &mut file.patches {
        replace(&mut p.faces);
    }
    for s in &mut file.surfaces {
        replace(&mut s.faces);
    }
    Ok(())
}

/// Apply one edit; the result is checked (BS-05) or refused.
pub fn apply(file: &TopologyFile, edit: &Edit) -> Result<(TopologyFile, Change), EngineError> {
    // the input must be a valid topology, so a refusal names this edit
    check(file)?;
    let mut out = file.clone();
    let mut change = Change::default();
    match edit {
        Edit::SetCount { edge, cells } => {
            if *cells == 0 {
                return Err(refused("a cell count must be at least 1"));
            }
            let mut marked: BTreeSet<(usize, usize)> = carriers(file, *edge)?.into_iter().collect();
            let mut keys: BTreeSet<(usize, usize)> = BTreeSet::new();
            // grow the chain: a block axis joins when it shares an edge
            loop {
                for &(b, a) in &marked {
                    keys.extend(axis_edges(file, b)[a].iter().copied());
                }
                let before = marked.len();
                for b in 0..file.blocks.len() {
                    for (a, edges) in axis_edges(file, b).iter().enumerate() {
                        if edges.iter().any(|k| keys.contains(k)) {
                            marked.insert((b, a));
                        }
                    }
                }
                if marked.len() == before {
                    break;
                }
            }
            for &(b, a) in &marked {
                if out.blocks[b].cells[a] != *cells {
                    out.blocks[b].cells[a] = *cells;
                    change.blocks.push(b);
                }
            }
            change.blocks.dedup();
            change.edges = keys.into_iter().map(|(p, q)| [p, q]).collect();
        }
        Edit::SetDistribution { edge, distribution } => {
            carriers(file, *edge)?;
            let (e, flipped) = edge_entry(&mut out, *edge);
            e.distribution = if flipped {
                reversed_law(distribution)
            } else {
                distribution.clone()
            };
            change.edges.push(*edge);
        }
        Edit::SetCurve { edge, curve } => {
            carriers(file, *edge)?;
            let (e, flipped) = edge_entry(&mut out, *edge);
            e.curve = if flipped {
                reversed_curve(curve)
            } else {
                curve.clone()
            };
            change.edges.push(*edge);
            change.reference_removed = out.reference.take().is_some();
        }
        Edit::MoveVertex { vertex, to, by } => {
            let Some(p) = out.vertices.get_mut(*vertex) else {
                return Err(refused(format!(
                    "vertex {vertex} does not exist (the file has {})",
                    file.vertices.len()
                )));
            };
            *p = match (to, by) {
                (Some(t), None) => *t,
                (None, Some(d)) => [p[0] + d[0], p[1] + d[1], p[2] + d[2]],
                _ => return Err(refused("give exactly one of `to` and `by`")),
            };
            if p.iter().any(|c| !c.is_finite()) {
                return Err(refused("the new position is not finite"));
            }
            change.vertices.push(*vertex);
            change.reference_removed = out.reference.take().is_some();
        }
        Edit::AddVertex { point } => {
            if point.iter().any(|c| !c.is_finite()) {
                return Err(refused("the point is not finite"));
            }
            out.vertices.push(*point);
            change.added_vertices.push(out.vertices.len() - 1);
        }
        Edit::RemoveVertex { vertex } => remove_vertex(&mut out, *vertex, &mut change)?,
        Edit::AddBlock {
            vertices,
            cells,
            domain,
            patch,
            role,
        } => add_block(
            &mut out,
            BlockFile {
                vertices: *vertices,
                cells: *cells,
                domain: domain.clone().unwrap_or_else(|| "fluid".into()),
            },
            patch.as_deref(),
            *role,
            &mut change,
        )?,
        Edit::RemoveBlock { block, patch, role } => {
            remove_block(&mut out, *block, patch.as_deref(), *role, &mut change)?
        }
        Edit::SplitBlock { block, axis, at } => {
            split_block(&mut out, *block, *axis, *at, &mut change)?
        }
        Edit::SplitTriangle {
            block,
            vertex,
            along,
            cells,
        } => {
            crate::block_triangle::split_triangle(
                &mut out,
                *block,
                *vertex,
                *along,
                *cells,
                &mut change,
            )?;
        }
        Edit::SetPatch { name, rename, role } => {
            if rename.is_none() && role.is_none() {
                return Err(refused("give `rename`, `role` or both"));
            }
            let Some(i) = out.patches.iter().position(|p| &p.name == name) else {
                let known: Vec<&str> = file.patches.iter().map(|p| p.name.as_str()).collect();
                return Err(EngineError::new(
                    ErrorCode::UnknownPatch,
                    format!("no patch `{name}` (patches: {})", known.join(", ")),
                ));
            };
            if let Some(new) = rename {
                if new != name && out.patches.iter().any(|p| &p.name == new) {
                    return Err(refused(format!("a patch `{new}` exists already")));
                }
                out.patches[i].name = new.clone();
            }
            if let Some(r) = role {
                out.patches[i].role = Some(*r);
            }
            change.patches.push(out.patches[i].name.clone());
        }
    }
    check(&out)?;
    Ok((out, change))
}

/// The contents of a topology for `block.list`: vertices, blocks, edges
/// with their curves and laws, patches with roles, surfaces, the check.
pub fn list(file: &TopologyFile) -> Result<Value, EngineError> {
    let checked = check(file)?;
    let curve = |c: &CurveFile| match c {
        CurveFile::Line => "line",
        CurveFile::Arc { .. } => "arc",
        CurveFile::Spline { .. } => "spline",
    };
    Ok(json!({
        "vertices": file.vertices.iter().enumerate()
            .map(|(i, p)| json!({"index": i, "point": p})).collect::<Vec<_>>(),
        "blocks": file.blocks.iter().enumerate()
            .map(|(i, b)| json!({"index": i, "vertices": b.vertices, "cells": b.cells, "domain": b.domain}))
            .collect::<Vec<_>>(),
        "edges": file.edges.iter()
            .map(|e| json!({"vertices": e.vertices, "curve": curve(&e.curve), "distribution": e.distribution}))
            .collect::<Vec<_>>(),
        "patches": file.patches.iter()
            .map(|p| json!({"name": p.name, "role": p.role, "faces": p.faces.len()}))
            .collect::<Vec<_>>(),
        "surfaces": file.surfaces.iter()
            .map(|s| json!({"name": s.name, "faces": s.faces.len()})).collect::<Vec<_>>(),
        "reference": file.reference.as_ref().map(|r| json!({"template": r.template, "analytic_volume": r.analytic_volume})),
        "check": checked,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::block::{TemplateConfig, template};

    fn pipe() -> TopologyFile {
        template(&TemplateConfig::PipeOgrid {
            radius: 0.05,
            length: 0.2,
            core_fraction: None,
            core_cells: 4,
            radial_cells: 3,
            axial_cells: 6,
            wall_spacing: None,
        })
        .unwrap()
    }

    /// Invariant (BS-04, BS-05): a count set on one edge reaches every
    /// block of its chain, so the result checks, and the cell total follows.
    #[test]
    fn a_count_follows_its_chain() {
        let f = pipe();
        // the axial edge of the core block runs from vertex slot 0 to 4
        let v = f.blocks[0].vertices;
        let (g, c) = apply(
            &f,
            &Edit::SetCount {
                edge: [v[0], v[4]],
                cells: 9,
            },
        )
        .unwrap();
        assert_eq!(c.blocks, [0, 1, 2, 3, 4]);
        assert!(g.blocks.iter().all(|b| b.cells[2] == 9));
        let cells = check(&g).unwrap()["cells"].as_u64().unwrap();
        assert_eq!(cells, 9 * (16 + 4 * 12));
        assert_eq!(g.reference, f.reference, "a count keeps the shape");
    }

    /// Invariant: a law given from the other end of an edge is stored
    /// reversed, so it means the same spacing at the same vertex; a law
    /// on a vertex pair that is no block edge is refused.
    #[test]
    fn a_law_is_stored_in_the_edge_direction() {
        let f = pipe();
        let v = f.blocks[0].vertices;
        let (g, _) = apply(
            &f,
            &Edit::SetDistribution {
                edge: [v[4], v[0]],
                distribution: DistributionFile::Geometric {
                    first: Some(0.01),
                    last: None,
                },
            },
        )
        .unwrap();
        let e = g
            .edges
            .iter()
            .find(|e| key(e.vertices[0], e.vertices[1]) == key(v[0], v[4]))
            .unwrap();
        let expected = if e.vertices == [v[4], v[0]] {
            DistributionFile::Geometric {
                first: Some(0.01),
                last: None,
            }
        } else {
            DistributionFile::Geometric {
                first: None,
                last: Some(0.01),
            }
        };
        assert_eq!(e.distribution, expected);
        let err = apply(
            &f,
            &Edit::SetDistribution {
                edge: [v[0], v[6]],
                distribution: DistributionFile::Uniform,
            },
        )
        .unwrap_err();
        assert!(err.to_string().contains("not an edge"), "{err}");
    }

    /// Invariant: moving a vertex that a projected face needs off its
    /// surface is refused by the check, nothing else changes; a move that
    /// keeps the topology valid drops the template's analytic reference.
    #[test]
    fn moves_are_checked_and_drop_the_reference() {
        let f = template(&TemplateConfig::Elbow {
            radius: 0.5,
            bend_radius: 1.5,
            angle: 90.0,
            inlet_length: 1.0,
            outlet_length: 1.0,
            core_fraction: None,
            core_cells: 4,
            radial_cells: 3,
            axial_cells: 12,
            wall_spacing: None,
        })
        .unwrap();
        // a wall vertex is a corner of a projected face
        let wall = f.surfaces[0].faces[0][0];
        let err = apply(
            &f,
            &Edit::MoveVertex {
                vertex: wall,
                to: None,
                by: Some([0.0, 0.0, 0.01]),
            },
        )
        .unwrap_err();
        assert!(err.to_string().contains("surface"), "{err}");
        // a corner of the pipe's core block is on no surface: a small move
        // inside the core keeps the topology valid
        let p = pipe();
        let core = p.blocks[0].vertices[4];
        let (g, c) = apply(
            &p,
            &Edit::MoveVertex {
                vertex: core,
                to: None,
                by: Some([0.001, 0.0, 0.0]),
            },
        )
        .unwrap();
        assert_eq!(g.vertices[core][0], p.vertices[core][0] + 0.001);
        assert!(c.reference_removed && g.reference.is_none());
        assert!(
            apply(
                &p,
                &Edit::MoveVertex {
                    vertex: 10_000,
                    to: Some([0.0; 3]),
                    by: None
                }
            )
            .is_err()
        );
    }

    /// Invariant: a patch keeps its faces when renamed or given a role; a
    /// taken name and an unknown patch are refused.
    #[test]
    fn patches_are_renamed_and_given_roles() {
        let f = pipe();
        let (g, c) = apply(
            &f,
            &Edit::SetPatch {
                name: "walls".into(),
                rename: Some("pipe_wall".into()),
                role: Some(PatchRole::Wall),
            },
        )
        .unwrap();
        assert_eq!(c.patches, ["pipe_wall"]);
        let p = g.patches.iter().find(|p| p.name == "pipe_wall").unwrap();
        assert_eq!(p.role, Some(PatchRole::Wall));
        assert_eq!(
            p.faces,
            f.patches.iter().find(|p| p.name == "walls").unwrap().faces
        );
        let taken = Edit::SetPatch {
            name: "walls".into(),
            rename: Some("inlet".into()),
            role: None,
        };
        assert!(apply(&f, &taken).is_err());
        let unknown = Edit::SetPatch {
            name: "nope".into(),
            rename: None,
            role: Some(PatchRole::Wall),
        };
        assert!(apply(&f, &unknown).is_err());
    }

    fn grid_of(f: &TopologyFile) -> algo::BlockGrid {
        algo::build(&f.topology().unwrap()).unwrap()
    }

    /// The largest distance from a node of `a` to the nearest node of `b`.
    fn node_gap(a: &algo::BlockGrid, b: &algo::BlockGrid) -> f64 {
        a.points
            .iter()
            .map(|p| {
                b.points
                    .iter()
                    .map(|q| (*p - *q).length())
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(0.0, f64::max)
    }

    /// Invariant (BS-04, BS-05): splitting a pipe O-grid block along its
    /// axis splits all five blocks of the chain, keeps the cell count, the
    /// wall projection and the volume, and puts the cut-face nodes where
    /// the unsplit grid had them (wall edges become arcs).
    #[test]
    fn a_split_follows_its_chain_and_keeps_the_grid() {
        let f = template(&TemplateConfig::PipeOgrid {
            radius: 0.05,
            length: 0.2,
            core_fraction: None,
            core_cells: 4,
            radial_cells: 5,
            axial_cells: 6,
            wall_spacing: Some(0.002),
        })
        .unwrap();
        let (g, c) = apply(
            &f,
            &Edit::SplitBlock {
                block: 2,
                axis: 2,
                at: 2,
            },
        )
        .unwrap();
        assert_eq!(g.blocks.len(), 10);
        // four core and four wall edges along the axis
        assert_eq!(c.added_vertices.len(), 8);
        let cells = |f: &TopologyFile| check(f).unwrap()["cells"].as_u64().unwrap();
        assert_eq!(cells(&g), cells(&f));
        let (before, after) = (grid_of(&f), grid_of(&g));
        let gap = node_gap(&before, &after);
        assert!(gap < 1e-9, "{gap}");
        let arcs = g
            .edges
            .iter()
            .filter(|e| matches!(e.curve, CurveFile::Arc { .. }))
            .count();
        assert!(arcs >= 4 * 3, "{arcs}");
    }

    /// Invariant: a split along a clustered direction off the middle puts
    /// every node where the unsplit grid had it (the halves' tanh laws are
    /// fitted to the measured end cells), and a split index outside the
    /// block is refused.
    #[test]
    fn a_split_keeps_the_clustering() {
        let f = template(&TemplateConfig::Channel {
            length: 1.0,
            width: 0.2,
            height: 0.1,
            cells: [8, 12, 6],
            wall_spacing: Some(0.002),
        })
        .unwrap();
        let (g, c) = apply(
            &f,
            &Edit::SplitBlock {
                block: 0,
                axis: 1,
                at: 4,
            },
        )
        .unwrap();
        assert_eq!(g.blocks.len(), 2);
        assert!(!c.fitted_edges.is_empty());
        let gap = node_gap(&grid_of(&f), &grid_of(&g));
        // a piece of a Vinokur distribution is again one, so the fitted
        // halves reproduce the nodes (measured 3e-14)
        assert!(gap < 1e-9, "{gap}");
        for at in [0, 12] {
            assert!(
                apply(
                    &f,
                    &Edit::SplitBlock {
                        block: 0,
                        axis: 1,
                        at
                    }
                )
                .is_err()
            );
        }
    }

    /// Invariant: a block added on the outlet takes the outlet face out of
    /// its patch (the empty patch is removed) and puts its open faces into
    /// the named patch; removing it again opens the face into a patch, and
    /// its now unused vertices can be removed, which gives back the cell
    /// count and vertex count of the start.
    #[test]
    fn blocks_are_added_and_removed() {
        let f = template(&TemplateConfig::Channel {
            length: 1.0,
            width: 0.2,
            height: 0.1,
            cells: [4, 3, 2],
            wall_spacing: None,
        })
        .unwrap();
        let outlet = f.patches.iter().find(|p| p.name == "outlet").unwrap().faces[0];
        // the outlet face is x = 1; the new block reaches to x = 1.5
        let v = f.blocks[0].vertices;
        let mut g = f.clone();
        let mut fresh = [0usize; 8];
        for (slot, &old) in v.iter().enumerate() {
            let p = f.vertices[old];
            if outlet.contains(&old) {
                fresh[slot] = old;
            } else {
                g = apply(
                    &g,
                    &Edit::AddVertex {
                        point: [1.5, p[1], p[2]],
                    },
                )
                .unwrap()
                .0;
                fresh[slot] = g.vertices.len() - 1;
            }
        }
        // the old block's high-x corners are the new block's low-x corners
        let mut nb = [0usize; 8];
        for (slot, corner) in nb.iter_mut().enumerate() {
            let (i, j, k) = match slot {
                0 => (0, 0, 0),
                1 => (1, 0, 0),
                2 => (1, 1, 0),
                3 => (0, 1, 0),
                4 => (0, 0, 1),
                5 => (1, 0, 1),
                6 => (1, 1, 1),
                _ => (0, 1, 1),
            };
            *corner = if i == 0 {
                v[corner_slot(1, j, k)]
            } else {
                fresh[corner_slot(0, j, k)]
            };
        }
        let refused_without_patch = apply(
            &g,
            &Edit::AddBlock {
                vertices: nb,
                cells: [2, 3, 2],
                domain: None,
                patch: None,
                role: None,
            },
        );
        assert!(refused_without_patch.is_err());
        let (h, c) = apply(
            &g,
            &Edit::AddBlock {
                vertices: nb,
                cells: [2, 3, 2],
                domain: None,
                patch: Some("extension".into()),
                role: Some(PatchRole::Wall),
            },
        )
        .unwrap();
        assert_eq!(c.removed_patches, ["outlet"]);
        assert_eq!(check(&h).unwrap()["cells"], 4 * 3 * 2 + 2 * 3 * 2);
        let (back, _) = apply(
            &h,
            &Edit::RemoveBlock {
                block: 1,
                patch: Some("outlet".into()),
                role: Some(PatchRole::Outlet),
            },
        )
        .unwrap();
        let mut back = back;
        while back.vertices.len() > f.vertices.len() {
            let last = back.vertices.len() - 1;
            back = apply(&back, &Edit::RemoveVertex { vertex: last })
                .unwrap()
                .0;
        }
        assert_eq!(check(&back).unwrap()["cells"], check(&f).unwrap()["cells"]);
        assert!(
            apply(&back, &Edit::RemoveVertex { vertex: 0 }).is_err(),
            "a used vertex stays"
        );
        assert!(
            apply(
                &back,
                &Edit::RemoveBlock {
                    block: 0,
                    patch: None,
                    role: None
                }
            )
            .is_err()
        );
    }

    /// The largest distance of a node on patch `name` of `a` to the
    /// nearest node of `b`.
    fn patch_gap(f: &TopologyFile, a: &algo::BlockGrid, b: &algo::BlockGrid, name: &str) -> f64 {
        let patch = f.patches.iter().position(|p| p.name == name).unwrap();
        let ids: BTreeSet<usize> = a
            .boundary
            .iter()
            .filter(|(_, p, _)| *p == patch)
            .flat_map(|(q, _, _)| *q)
            .collect();
        ids.iter()
            .map(|&i| {
                b.points
                    .iter()
                    .map(|q| (a.points[i] - *q).length())
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(0.0, f64::max)
    }

    /// Invariant: splitting the airfoil C-grid's front block along the wall
    /// keeps the spline wall's own points in each half and fits the halves'
    /// laws by arc length, so the wall nodes stay (measured 1.5e-8 chords).
    /// Interior nodes follow the TFI of the new blocks, which is not the
    /// TFI of the whole: they move by up to 2.2e-3 chords on the airfoil
    /// and, on the elbow, the projected wall face interiors slide along the
    /// torus by 4.8e-5 (radius 0.5), all measured.
    #[test]
    fn curved_walls_survive_a_split() {
        let airfoil = template(&TemplateConfig::AirfoilCgrid {
            naca: Some("2412".into()),
            coordinates: None,
            chord: 1.0,
            span: 0.1,
            farfield: 5.0,
            wake_length: 5.0,
            surface_cells: 40,
            normal_cells: 12,
            wake_cells: 10,
            span_cells: 1,
            wall_spacing: 2e-3,
            leading_edge_spacing: None,
            trailing_edge_spacing: None,
        })
        .unwrap();
        let (g, _) = apply(
            &airfoil,
            &Edit::SplitBlock {
                block: 0,
                axis: 0,
                at: 11,
            },
        )
        .unwrap();
        assert_eq!(g.blocks.len(), 7);
        let (a, b) = (grid_of(&airfoil), grid_of(&g));
        let wall = patch_gap(&airfoil, &a, &b, "airfoil_wall");
        let all = node_gap(&a, &b);
        assert!(wall < 1e-7, "{wall}");
        assert!(all < 3e-3, "{all}");

        let elbow = template(&TemplateConfig::Elbow {
            radius: 0.5,
            bend_radius: 1.5,
            angle: 90.0,
            inlet_length: 1.0,
            outlet_length: 1.0,
            core_fraction: None,
            core_cells: 4,
            radial_cells: 3,
            axial_cells: 24,
            wall_spacing: None,
        })
        .unwrap();
        let (g, _) = apply(
            &elbow,
            &Edit::SplitBlock {
                block: 0,
                axis: 0,
                at: 3,
            },
        )
        .unwrap();
        let (a, b) = (grid_of(&elbow), grid_of(&g));
        let all = node_gap(&a, &b);
        assert!(all < 1e-4, "{all}");
    }
}
