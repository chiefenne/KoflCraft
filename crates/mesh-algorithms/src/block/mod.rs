//! The block-structured mesher variant (OBJ-21, section 4c): checks a block
//! topology (BS-02, BS-05), distributes nodes on its edges (BS-06), fills
//! faces and blocks by transfinite interpolation (BS-09, TFI only in MB1a)
//! and assembles point-matched hexahedra whose shared faces are built once.
//!
//! Design notes and sources: `research/block-structured-concept.md`.

pub mod counts;
pub mod distribution;
pub mod grid;
pub mod junction;
pub mod layout;
pub mod multi_element;
pub mod overlay;
pub mod quality3d;
pub mod smooth;
pub mod templates;
pub mod tri_center;

use std::collections::{BTreeMap, HashMap};

use mesh_core::block::{
    BLOCK_FACES, BlockTopology, Curve, Distribution, block_edges, corner_slot, face_corner_slot,
    face_key,
};
use mesh_core::geometry::Vec3;

pub use grid::{BlockGrid, BlockQuality, BlockRange, EdgeReport, build};
pub use junction::{Cover, Junctions};

/// One finding of the topology check.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Issue {
    /// `vertex`, `block`, `edge`, `count_chain`, `face`, `patch` or
    /// `surface`.
    pub kind: &'static str,
    pub message: String,
}

fn issue(kind: &'static str, message: impl Into<String>) -> Issue {
    Issue {
        kind,
        message: message.into(),
    }
}

/// A chain of opposite edges that must carry the same cell count (BS-05).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Chain {
    pub cells: usize,
    pub edges: usize,
    pub blocks: Vec<usize>,
}

/// What the check derived from a consistent topology.
#[derive(Clone, Debug, PartialEq)]
pub struct Checked {
    pub cells: usize,
    pub points: usize,
    pub block_faces_shared: usize,
    pub boundary_faces: BTreeMap<String, usize>,
    pub chains: Vec<Chain>,
    /// The T-junctions and the coarse faces they make (BS-21).
    pub junctions: Junctions,
    /// The cell count of every block edge by its key.
    pub edge_counts: HashMap<(usize, usize), usize>,
}

pub(crate) fn edge_key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

/// The vertex slots at the start and end of each block edge.
pub(crate) fn edge_slots() -> [(usize, usize, usize); 12] {
    block_edges().map(|(axis, c)| {
        let mut e = c;
        e[axis] = 1;
        (
            axis,
            corner_slot(c[0], c[1], c[2]),
            corner_slot(e[0], e[1], e[2]),
        )
    })
}

/// The four vertices of block face `f` in local order (00, 10, 11, 01).
pub(crate) fn face_vertices(verts: &[usize; 8], f: usize) -> [usize; 4] {
    let face = BLOCK_FACES[f];
    [(0, 0), (1, 0), (1, 1), (0, 1)].map(|(a, b)| verts[face_corner_slot(face, a, b)])
}

struct UnionFind(Vec<usize>);

impl UnionFind {
    fn find(&mut self, mut x: usize) -> usize {
        while self.0[x] != x {
            self.0[x] = self.0[self.0[x]];
            x = self.0[x];
        }
        x
    }
    fn union(&mut self, a: usize, b: usize) {
        let (a, b) = (self.find(a), self.find(b));
        self.0[a] = b;
    }
}

/// Check a topology before any grid is built: indices, right-handed
/// blocks, edge curves and laws, consistent counts along every chain of
/// opposite edges, faces shared by at most two blocks, every boundary face
/// in exactly one patch.
pub fn check(t: &BlockTopology) -> Result<Checked, Vec<Issue>> {
    let mut issues = Vec::new();
    let nv = t.vertices.len();
    for (v, p) in t.vertices.iter().enumerate() {
        if !p.is_finite() {
            issues.push(issue("vertex", format!("vertex {v} is not finite")));
        }
    }
    if t.blocks.is_empty() {
        issues.push(issue("block", "the topology has no blocks"));
    }
    let mut cells: usize = 0;
    let mut valid_blocks = Vec::new();
    for (b, block) in t.blocks.iter().enumerate() {
        let v = block.vertices;
        if v.iter().any(|&i| i >= nv) {
            issues.push(issue(
                "block",
                format!("block {b} names a vertex that does not exist"),
            ));
            continue;
        }
        if (0..8).any(|i| v[..i].contains(&v[i])) {
            issues.push(issue(
                "block",
                format!("block {b} repeats a vertex (collapsed blocks come later)"),
            ));
            continue;
        }
        if block.cells.contains(&0) {
            issues.push(issue("block", format!("block {b} has a zero cell count")));
            continue;
        }
        if block.domain.is_empty() {
            issues.push(issue("block", format!("block {b} has no domain name")));
        }
        if let Some(c) = (0..8).find(|&s| curved_corner_det(t, &v, s) <= 0.0) {
            issues.push(issue(
                "block",
                format!(
                    "block {b} is not a right-handed hexahedron at its corner {c} (vertex {}); list the vertices counter-clockwise about the k axis",
                    v[c]
                ),
            ));
            continue;
        }
        match block
            .cells
            .iter()
            .try_fold(1usize, |a, &n| a.checked_mul(n))
        {
            Some(n) => cells = cells.saturating_add(n),
            None => issues.push(issue("block", format!("block {b} cell count overflows"))),
        }
        valid_blocks.push(b);
    }

    // Edges of the blocks, with their counts and chains.
    let mut edges: HashMap<(usize, usize), (usize, usize, usize)> = HashMap::new(); // key → (count, block, edge id)
    let mut edge_ids: HashMap<(usize, usize), usize> = HashMap::new();
    let mut uf = UnionFind(Vec::new());
    let mut chain_blocks: Vec<Vec<usize>> = Vec::new();
    for &b in &valid_blocks {
        let block = &t.blocks[b];
        let mut first_of_axis = [usize::MAX; 3];
        for (axis, s, e) in edge_slots() {
            let key = edge_key(block.vertices[s], block.vertices[e]);
            let n = block.cells[axis];
            let id = *edge_ids.entry(key).or_insert_with(|| {
                uf.0.push(uf.0.len());
                chain_blocks.push(Vec::new());
                uf.0.len() - 1
            });
            chain_blocks[id].push(b);
            if first_of_axis[axis] == usize::MAX {
                first_of_axis[axis] = id;
            } else {
                uf.union(id, first_of_axis[axis]);
            }
            match edges.get(&key) {
                Some(&(m, other, _)) if m != n => issues.push(issue(
                    "count_chain",
                    format!(
                        "edge {}–{} has {m} cells in block {other} and {n} in block {b}; every chain of opposite edges needs one count",
                        key.0, key.1
                    ),
                )),
                Some(_) => {}
                None => {
                    edges.insert(key, (n, b, id));
                }
            }
        }
    }
    let mut chains: BTreeMap<usize, Chain> = BTreeMap::new();
    for (&_key, &(n, _, id)) in &edges {
        let root = uf.find(id);
        let c = chains.entry(root).or_insert(Chain {
            cells: n,
            edges: 0,
            blocks: Vec::new(),
        });
        c.edges += 1;
        c.blocks.extend(&chain_blocks[id]);
    }
    let mut chains: Vec<Chain> = chains
        .into_values()
        .map(|mut c| {
            c.blocks.sort_unstable();
            c.blocks.dedup();
            c
        })
        .collect();
    chains.sort_by(|a, b| a.blocks.cmp(&b.blocks).then(a.cells.cmp(&b.cells)));

    // T-junctions: each piece an edge of a block, the counts adding up.
    let mut junctions = Junctions::read(t, &mut issues);
    let mut jkeys: Vec<_> = junctions.through.keys().copied().collect();
    jkeys.sort_unstable();
    for key in jkeys {
        let chain = junctions.through[&key].clone();
        let Some(&(whole, _, _)) = edges.get(&key) else {
            issues.push(issue(
                "junction",
                format!(
                    "the junction edge {}–{} is no edge of any block",
                    key.0, key.1
                ),
            ));
            continue;
        };
        let mut sum = 0;
        let mut complete = true;
        for w in chain.windows(2) {
            match edges.get(&edge_key(w[0], w[1])) {
                Some(&(n, _, _)) => sum += n,
                None => {
                    complete = false;
                    issues.push(issue(
                        "junction",
                        format!(
                            "the piece {}–{} of junction edge {}–{} is no edge of any block",
                            w[0], w[1], key.0, key.1
                        ),
                    ));
                }
            }
        }
        if complete && sum != whole {
            issues.push(issue(
                "count_chain",
                format!(
                    "edge {}–{} has {whole} cells but its pieces through the T-junctions have {sum}; the counts must add up",
                    key.0, key.1
                ),
            ));
        }
    }

    // Explicit edges.
    let mut seen = HashMap::new();
    for (e, spec) in t.edges.iter().enumerate() {
        let [a, b] = spec.vertices;
        if a >= nv || b >= nv || a == b {
            issues.push(issue("edge", format!("edge spec {e} has invalid vertices")));
            continue;
        }
        let key = edge_key(a, b);
        if seen.insert(key, e).is_some() {
            issues.push(issue("edge", format!("edge {a}–{b} is given twice")));
            continue;
        }
        if junctions.through.contains_key(&key) {
            issues.push(issue(
                "edge",
                format!("edge {a}–{b} passes through T-junctions; give the curves of its pieces"),
            ));
            continue;
        }
        let Some(&(n, _, _)) = edges.get(&key) else {
            issues.push(issue(
                "edge",
                format!("edge {a}–{b} is no edge of any block"),
            ));
            continue;
        };
        match distribution::Evaluator::new(t.vertices[a], t.vertices[b], &spec.curve) {
            Ok(ev) => {
                if let Err(m) = distribution::fractions(spec.distribution, n, ev.length()) {
                    issues.push(issue("edge", format!("edge {a}–{b}: {m}")));
                }
            }
            Err(m) => issues.push(issue("edge", format!("edge {a}–{b}: {m}"))),
        }
    }

    // Faces and patches.
    let mut faces: HashMap<[usize; 4], Vec<usize>> = HashMap::new();
    for &b in &valid_blocks {
        for f in 0..6 {
            faces
                .entry(face_key(face_vertices(&t.blocks[b].vertices, f)))
                .or_default()
                .push(b);
        }
    }
    let mut shared = 0;
    for (key, blocks) in &faces {
        if blocks.len() > 2 {
            issues.push(issue(
                "face",
                format!("face {key:?} is shared by more than two blocks: {blocks:?}"),
            ));
        } else if blocks.len() == 2 {
            shared += 1;
        }
    }
    let mut patch_of: HashMap<[usize; 4], usize> = HashMap::new();
    let mut boundary_faces = BTreeMap::new();
    for (p, patch) in t.patches.iter().enumerate() {
        if patch.name.is_empty() || t.patches[..p].iter().any(|q| q.name == patch.name) {
            issues.push(issue(
                "patch",
                format!("patch {p} has an empty or repeated name"),
            ));
        }
        if patch.faces.is_empty() {
            issues.push(issue("patch", format!("patch {} has no faces", patch.name)));
        }
        let mut count = 0;
        for &face in &patch.faces {
            let key = face_key(face);
            match faces.get(&key).map(Vec::len) {
                Some(1) => {}
                Some(_) => {
                    issues.push(issue(
                        "patch",
                        format!(
                            "patch {}: face {face:?} lies between two blocks",
                            patch.name
                        ),
                    ));
                    continue;
                }
                None => {
                    issues.push(issue(
                        "patch",
                        format!("patch {}: {face:?} is no block face", patch.name),
                    ));
                    continue;
                }
            }
            if let Some(q) = patch_of.insert(key, p) {
                issues.push(issue(
                    "patch",
                    format!(
                        "face {face:?} is in patch {} and patch {}",
                        t.patches[q].name, patch.name
                    ),
                ));
            }
            count += 1;
        }
        boundary_faces.insert(patch.name.clone(), count);
    }
    let in_patch: std::collections::HashSet<[usize; 4]> = t
        .patches
        .iter()
        .flat_map(|p| p.faces.iter().map(|&f| face_key(f)))
        .collect();
    junctions.find_covers(t, &valid_blocks, &faces, &in_patch, &mut issues);
    for k in junctions.sub.keys() {
        if in_patch.contains(k) {
            issues.push(issue(
                "patch",
                format!("face {k:?} lies inside a coarse face and cannot be in a patch"),
            ));
        }
    }
    for cover in &junctions.covers {
        let m = junctions.coarse_dims(t, cover);
        let count = |p: usize, q: usize| edges.get(&edge_key(p, q)).map_or(0, |e| e.0);
        let (run, across) = if cover.along_first { (0, 1) } else { (1, 0) };
        let mut sum = 0;
        for s in &cover.subs {
            let side = [count(s[0], s[1]), count(s[0], s[3])];
            sum += side[run];
            if side[across] != m[across] {
                issues.push(issue(
                    "count_chain",
                    format!(
                        "face {:?} has {} cells across but the coarse face {:?} of block {} has {}",
                        face_key(*s),
                        side[across],
                        face_key(cover.corners),
                        cover.block,
                        m[across]
                    ),
                ));
            }
        }
        if sum != m[run] {
            issues.push(issue(
                "count_chain",
                format!(
                    "the faces covering face {:?} of block {} have {sum} cells along it, the face has {}",
                    face_key(cover.corners),
                    cover.block,
                    m[run]
                ),
            ));
        }
    }
    for (e, spec) in t.edges.iter().enumerate() {
        let key = edge_key(spec.vertices[0], spec.vertices[1]);
        if junctions.cuts.contains_key(&key) {
            issues.push(issue(
                "edge",
                format!(
                    "edge spec {e} ({}–{}) lies inside a coarse face, whose grid places its nodes",
                    key.0, key.1
                ),
            ));
        }
    }
    let mut open: Vec<_> = faces
        .iter()
        .filter(|(k, b)| {
            b.len() == 1
                && !patch_of.contains_key(*k)
                && !junctions.coarse.contains_key(*k)
                && !junctions.sub.contains_key(*k)
        })
        .map(|(k, b)| (b[0], *k))
        .collect();
    open.sort_unstable();
    for (b, key) in open {
        issues.push(issue(
            "patch",
            format!("boundary face {key:?} of block {b} is in no patch"),
        ));
    }

    // Surfaces: valid, on boundary faces only, each face projected once,
    // its corners on the surface.
    let mut projected: HashMap<[usize; 4], usize> = HashMap::new();
    for (sn, spec) in t.surfaces.iter().enumerate() {
        if spec.name.is_empty() || t.surfaces[..sn].iter().any(|q| q.name == spec.name) {
            issues.push(issue(
                "surface",
                format!("surface {sn} has an empty or repeated name"),
            ));
        }
        if let Some(m) = spec.surface.problem() {
            issues.push(issue("surface", format!("surface {}: {m}", spec.name)));
            continue;
        }
        let tol = 1e-6 * spec.surface.scale();
        for &face in &spec.faces {
            let key = face_key(face);
            if faces.get(&key).map(Vec::len) != Some(1)
                || junctions.coarse.contains_key(&key)
                || junctions.sub.contains_key(&key)
            {
                issues.push(issue(
                    "surface",
                    format!("surface {}: {face:?} is no boundary block face", spec.name),
                ));
                continue;
            }
            if let Some(q) = projected.insert(key, sn) {
                issues.push(issue(
                    "surface",
                    format!(
                        "face {face:?} is projected onto {} and {}",
                        t.surfaces[q].name, spec.name
                    ),
                ));
            }
            for v in face {
                let p = t.vertices[v];
                if (spec.surface.project(p) - p).length() > tol {
                    issues.push(issue(
                        "surface",
                        format!(
                            "surface {}: vertex {v} of face {face:?} is not on the surface",
                            spec.name
                        ),
                    ));
                }
            }
        }
    }

    if !issues.is_empty() {
        return Err(issues);
    }
    // points: vertices + edge interiors + face interiors + block interiors
    let mut points = nv;
    for (key, &(n, _, _)) in &edges {
        if !junctions.through.contains_key(key) && !junctions.cuts.contains_key(key) {
            points += n - 1;
        }
    }
    let mut counted = std::collections::HashSet::new();
    for &b in &valid_blocks {
        let n = t.blocks[b].cells;
        for (f, &(_, _, a, c)) in BLOCK_FACES.iter().enumerate() {
            let key = face_key(face_vertices(&t.blocks[b].vertices, f));
            if !junctions.sub.contains_key(&key) && counted.insert(key) {
                points += (n[a] - 1) * (n[c] - 1);
            }
        }
        points += (n[0] - 1) * (n[1] - 1) * (n[2] - 1);
    }
    let edge_counts = edges.iter().map(|(&k, &(n, _, _))| (k, n)).collect();
    Ok(Checked {
        cells,
        points,
        block_faces_shared: shared,
        boundary_faces,
        chains,
        junctions,
        edge_counts,
    })
}

/// The determinant of the three edges leaving a corner slot of a block,
/// each along its curve's tangent where the edge is curved (scaled to its
/// chord), else along its chord; positive for a right-handed corner. An
/// arc or spline may leave a corner at a right angle where the chords are
/// collinear, as at an airfoil's leading edge.
pub(crate) fn curved_corner_det(t: &BlockTopology, v: &[usize; 8], slot: usize) -> f64 {
    let p = v.map(|i| t.vertices[i]);
    let mut e = corner_edges(&p, slot);
    let c = slot_corner(slot);
    for (axis, dir) in e.iter_mut().enumerate() {
        let mut o = c;
        o[axis] = 1 - c[axis];
        let (from, to) = (v[slot], v[corner_slot(o[0], o[1], o[2])]);
        let Some(spec) = t
            .edges
            .iter()
            .find(|s| edge_key(s.vertices[0], s.vertices[1]) == edge_key(from, to))
        else {
            continue;
        };
        if spec.curve == Curve::Line {
            continue;
        }
        let [a, b] = spec.vertices;
        let Ok(ev) = distribution::Evaluator::new(t.vertices[a], t.vertices[b], &spec.curve) else {
            continue;
        };
        let tangent = ev.end_tangent(a == from);
        if tangent.is_finite() {
            let sign = if c[axis] == 0 { 1.0 } else { -1.0 };
            *dir = tangent * (sign * (t.vertices[to] - t.vertices[from]).length());
        }
    }
    e[0].cross(e[1]).dot(e[2])
}

pub(crate) fn corner_edges(p: &[Vec3; 8], slot: usize) -> [Vec3; 3] {
    let c = slot_corner(slot);
    let mut out = [Vec3::ZERO; 3];
    for axis in 0..3 {
        let mut o = c;
        o[axis] = 1 - c[axis];
        let sign = if c[axis] == 0 { 1.0 } else { -1.0 };
        out[axis] = (p[corner_slot(o[0], o[1], o[2])] - p[slot]) * sign;
    }
    out
}

pub(crate) fn slot_corner(slot: usize) -> [usize; 3] {
    let k = slot / 4;
    let (i, j) = match slot % 4 {
        0 => (0, 0),
        1 => (1, 0),
        2 => (1, 1),
        _ => (0, 1),
    };
    [i, j, k]
}

/// The default of an edge that is not given explicitly.
pub(crate) const DEFAULT_EDGE: (Curve, Distribution) = (Curve::Line, Distribution::Uniform);

#[cfg(test)]
mod tests;
