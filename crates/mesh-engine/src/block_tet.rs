//! Tet split of a block (BS-20, step 2 of the 3D tet/prism insertion): a
//! hexahedral block whose shape is a tetrahedron becomes the four
//! hexahedra of the tetrahedron, one at each corner around a centre
//! vertex of valence 4.
//!
//! **The shape.** In the right local frame the corners of the tetrahedron
//! are P0 = v0, P1 = v1, P2 = v3 and P3 = v4. The block's three faces at
//! v0 are the triangles P0P1P2 (k = 0), P0P1P3 (j = 0) and P0P2P3
//! (i = 0); each has one flat corner, v2 on side P1P2, v5 on side P1P3
//! and v7 on side P2P3. The other three faces of the block meet at v6
//! inside the fourth triangle P1P2P3 and are already its three
//! quadrilaterals. A block is a tet when the three flat corners are
//! flatter than 150° across their faces and v6 is nearly flat.
//!
//! **Cell counts.** Give each corner i a count n_i; the hexahedron at
//! corner i has n_j cells along its side toward corner j, so tet edge
//! (i, j) has n_i + n_j cells. The block's counts (a, b, c) already give
//! n1 = a, n2 = b, n3 = c on the face P1P2P3, which stays as it is. The
//! fourth count n0 is new: the edges P0–P1, P0–P2 and P0–P3 get n0 more
//! cells, and each of them gets a new vertex n0 nodes from its far end,
//! on the edge's own curve and law. The three triangles at P0 become
//! three quadrilaterals each, around their centres.
//!
//! This is the classic tet-to-hex split (the "tri-center" map of a
//! triangle lifted to a tetrahedron), as in hex-dominant meshing and in
//! the integer-grid hex layouts of Lyon, Campen and Kobbelt (2021) and
//! Brückler, Bommes and Campen (2022). No code was read.

use std::collections::BTreeSet;

use mesh_algorithms::block as algo;
use mesh_core::block::{corner_slot, face_key};
use mesh_core::geometry::Vec3;
use serde::{Deserialize, Serialize};

use crate::EngineError;
use crate::block::{BlockFile, TopologyFile};
use crate::block_edit::{Change, a3, face_vertices, key, refused, v3};
use crate::block_triangle::{
    CountSolve, RowStop, Start, corner_of, cross_angle, rotated, rotations, split_prisms,
    t_junctions,
};

/// A tet's flat corners must be flatter than this (degrees) for the
/// detection; a split by hand accepts [`MIN_TET_CORNER`].
pub const FLAT_TET_CORNER: f64 = 150.0;
/// The least flat-corner angle a split by hand accepts.
pub const MIN_TET_CORNER: f64 = 120.0;
/// Largest scaled Jacobian at v6 (the point inside the triangle P1P2P3)
/// for the block to count as a tet: 0 when the three faces there are
/// coplanar, 1 at a right-angled hexahedron corner.
const MAX_INNER_JACOBIAN: f64 = 0.5;

/// A block whose shape is a tetrahedron, in a frame that puts its corner
/// P0 at v0.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct TetBlock {
    pub block: usize,
    /// The vertex at P0, the corner of the three triangles that are split.
    pub corner: usize,
    /// The tet's other corners P1, P2, P3.
    pub corners: [usize; 3],
    /// The flat vertices on the sides P1P2, P1P3 and P2P3.
    pub flat: [usize; 3],
    /// The least flat of the three flat corners (degrees).
    pub corner_deg: f64,
}

/// What a tet split did.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct TetSplit {
    pub block: usize,
    pub corner: usize,
    pub corner_deg: f64,
    /// The new count n0 (cells from the tet's centre toward P0's
    /// neighbours' sides).
    pub cells: usize,
    /// The four blocks, at P0, P1, P2 and P3.
    pub blocks: [usize; 4],
    /// The blocks on the triangles at P0 and stacked on them, each split
    /// into three (a triangle split) so the faces stay point-matched.
    pub neighbours: Vec<usize>,
    /// Blocks the new row cut in two (`sheet`), and with `t_junction` what
    /// the count solve did.
    pub row_blocks: Vec<usize>,
    pub stop: RowStop,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub counts: Option<CountSolve>,
    pub cells_before: usize,
    pub cells_after: usize,
    /// How the tet meets its neighbours: `conformal` (prism columns) or
    /// `non_conformal` (interfaces joined into polyhedra, BS-16).
    pub join: &'static str,
    /// The coarse patch of every interface the split made.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub interfaces: Vec<String>,
    /// Why the conformal split was refused, when the automatic join fell
    /// back to non-conformal.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fallback: Option<String>,
}

/// How a split piece meets its neighbours (card A of step 3,
/// 2026-11-28).
#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Join {
    /// Conformal where the neighbours can take it, else non-conformal.
    #[default]
    Auto,
    /// Conformal only; refused where the neighbours cannot take it.
    ConformalOnly,
    /// Non-conformal: the neighbours stay as they are.
    NonConformal,
}

/// The rotations that put a tet's P0 at v0, with the least flat angle of
/// its three flat corners, best first. `angle(slot, axis)` is the
/// cross-section angle in the block's own frame.
fn tet_frames(
    block: &BlockFile,
    angle: impl Fn(usize, usize) -> f64,
    inner: impl Fn(&BlockFile) -> f64,
) -> Vec<(BlockFile, f64)> {
    let mut out = Vec::new();
    for (p, s) in rotations() {
        // new slot → old slot
        let old_slot = |slot: usize| {
            let c = corner_of(slot);
            let mut o = [0usize; 3];
            for d in 0..3 {
                o[p[d]] = if s[d] { 1 - c[d] } else { c[d] };
            }
            corner_slot(o[0], o[1], o[2])
        };
        // the flat corner at new slot 2 lies across new axis 2, 5 across
        // 1 and 7 across 0
        let deg = [(2, 2), (5, 1), (7, 0)]
            .into_iter()
            .map(|(slot, axis)| angle(old_slot(slot), p[axis]))
            .fold(f64::INFINITY, f64::min);
        let r = rotated(block, p, s);
        if inner(&r) <= MAX_INNER_JACOBIAN {
            out.push((r, deg));
        }
    }
    out.sort_by(|x, y| y.1.total_cmp(&x.1));
    out
}

/// |scaled Jacobian| at v6 of a block, on its straight corner edges.
fn inner_jacobian(file: &TopologyFile, b: &BlockFile) -> f64 {
    let p = |s: usize| v3(file.vertices[b.vertices[s]]);
    let e = [p(7) - p(6), p(5) - p(6), p(2) - p(6)];
    let l: f64 = e.iter().map(|x| x.length()).product();
    if l == 0.0 {
        return 0.0;
    }
    (e[0].cross(e[1]).dot(e[2]) / l).abs()
}

fn grid(file: &TopologyFile) -> Result<algo::BlockGrid, EngineError> {
    algo::build(&file.topology()?).map_err(|i| {
        refused(
            i.iter()
                .map(|x| x.message.clone())
                .collect::<Vec<_>>()
                .join("; "),
        )
    })
}

fn add(file: &mut TopologyFile, p: Vec3, change: &mut Change) -> usize {
    file.vertices.push(a3(p));
    let n = file.vertices.len() - 1;
    change.added_vertices.push(n);
    n
}

/// The non-conformal split with the coarsest n0, from half the tet's
/// largest count up, whose four hexahedra reach a worst scaled Jacobian of
/// [`crate::block_insert::MIN_SCALED_JACOBIAN`] (card E of step 3, the
/// closest the piece's cells come to the neighbours' size at that
/// quality); the best one tried when none does.
fn non_conformal_count(
    file: &mut TopologyFile,
    b: usize,
    corner: Option<usize>,
    stop: RowStop,
    change: &mut Change,
) -> Result<TetSplit, EngineError> {
    let Some(block) = file.blocks.get(b) else {
        return Err(refused(format!("there is no block {b}")));
    };
    let largest = block.cells.iter().copied().max().unwrap_or(1);
    let mut best: Option<(f64, TopologyFile, Change, TetSplit)> = None;
    for n in largest.div_ceil(2).max(1)..=largest.max(1) {
        let mut f = file.clone();
        let mut c = change.clone();
        let split = split_tet_with(&mut f, b, corner, Some(n), stop, true, &mut c)?;
        let g = grid(&f)?;
        let (_, per) = algo::quality3d::measure(&g);
        let sj = split
            .blocks
            .iter()
            .map(|&x| per[x].min_scaled_jacobian)
            .fold(f64::INFINITY, f64::min);
        let good = sj >= crate::block_insert::MIN_SCALED_JACOBIAN;
        if best.as_ref().is_none_or(|x| sj > x.0) {
            best = Some((sj, f, c, split));
        }
        if good {
            break;
        }
    }
    let Some((_, f, c, split)) = best else {
        return Err(refused(format!("block {b} has no cells")));
    };
    *file = f;
    *change = c;
    Ok(split)
}

/// The tet-shaped blocks of a topology whose flat corners are flatter
/// than `min_deg`, flattest first.
pub fn tet_blocks(file: &TopologyFile, min_deg: f64) -> Result<Vec<TetBlock>, EngineError> {
    let g = grid(file)?;
    let mut out = Vec::new();
    for (b, block) in file.blocks.iter().enumerate() {
        let frames = tet_frames(
            block,
            |slot, axis| cross_angle(&g, file, b, slot, axis),
            |r| inner_jacobian(file, r),
        );
        if let Some((r, deg)) = frames.into_iter().next()
            && deg >= min_deg
        {
            out.push(described(b, &r, deg));
        }
    }
    out.sort_by(|x, y| y.corner_deg.total_cmp(&x.corner_deg));
    Ok(out)
}

fn described(b: usize, r: &BlockFile, deg: f64) -> TetBlock {
    let v = r.vertices;
    TetBlock {
        block: b,
        corner: v[0],
        corners: [v[1], v[3], v[4]],
        flat: [v[2], v[5], v[7]],
        corner_deg: deg,
    }
}

/// Split the tet-shaped block `b` into four hexahedra (BS-20). `corner`
/// picks P0 (the best frame when absent); `cells` is n0, from the tet's
/// counts when absent (the mean of a, b and c, which puts the new
/// vertices near the middles of the edges at P0).
///
/// For now the three triangles at P0 must be boundary faces: a block that
/// shares one of them, or one of the edges at P0, is refused.
pub(crate) fn split_tet(
    file: &mut TopologyFile,
    b: usize,
    corner: Option<usize>,
    cells: Option<usize>,
    stop: RowStop,
    join: Join,
    change: &mut Change,
) -> Result<TetSplit, EngineError> {
    match join {
        Join::ConformalOnly => split_tet_with(file, b, corner, cells, stop, false, change),
        Join::NonConformal => split_tet_with(file, b, corner, cells, stop, true, change),
        Join::Auto => {
            let mut trial = file.clone();
            let mut trial_change = change.clone();
            match split_tet_with(&mut trial, b, corner, cells, stop, false, &mut trial_change) {
                Ok(done) => {
                    *file = trial;
                    *change = trial_change;
                    Ok(done)
                }
                Err(conformal) => {
                    let mut trial = file.clone();
                    let mut trial_change = change.clone();
                    match split_tet_with(
                        &mut trial,
                        b,
                        corner,
                        cells,
                        stop,
                        true,
                        &mut trial_change,
                    ) {
                        Ok(mut done) if !done.interfaces.is_empty() => {
                            *file = trial;
                            *change = trial_change;
                            done.fallback = Some(conformal.message);
                            change.fallback = done.fallback.clone();
                            Ok(done)
                        }
                        _ => Err(conformal),
                    }
                }
            }
        }
    }
}

/// [`split_tet`] in one join mode.
fn split_tet_with(
    file: &mut TopologyFile,
    b: usize,
    corner: Option<usize>,
    cells: Option<usize>,
    stop: RowStop,
    non_conformal: bool,
    change: &mut Change,
) -> Result<TetSplit, EngineError> {
    if non_conformal && cells.is_none() {
        return non_conformal_count(file, b, corner, stop, change);
    }
    if b >= file.blocks.len() {
        return Err(refused(format!("there is no block {b}")));
    }
    let g = grid(file)?;
    let frames = tet_frames(
        &file.blocks[b],
        |slot, axis| cross_angle(&g, file, b, slot, axis),
        |r| inner_jacobian(file, r),
    );
    let pick = frames
        .into_iter()
        .find(|(r, _)| corner.is_none_or(|c| r.vertices[0] == c));
    let Some((r, deg)) = pick else {
        return Err(refused(match corner {
            Some(c) => format!(
                "block {b} is not a tetrahedron with its corner P0 at vertex {c} (three faces at \
                 the corner, each a triangle with one flat corner)"
            ),
            None => format!("block {b} is not shaped like a tetrahedron"),
        }));
    };
    if deg < MIN_TET_CORNER {
        return Err(refused(format!(
            "block {b} has a flat corner of {deg:.1}°; a tet split needs its three flat corners \
             at {MIN_TET_CORNER}° or more"
        )));
    }
    let v = r.vertices;
    let [a, bb, c] = r.cells;
    let (p0, p1, p2, p3) = (v[0], v[1], v[3], v[4]);
    let (m12, m13, m23, f0) = (v[2], v[5], v[7], v[6]);
    let cells_before: usize = file
        .blocks
        .iter()
        .map(|x| x.cells.iter().product::<usize>())
        .sum();
    let n0 = match cells {
        Some(0) => return Err(refused("the new count needs at least one cell")),
        Some(n) => n,
        None => (((a + bb + c) as f64 / 3.0).round() as usize).max(1),
    };
    // the triangles at P0 (faces k = 0, j = 0, i = 0 of the frame), each
    // with its flat vertex; a block on one is the start of a column of
    // triangle splits that carries the new row away from the tet
    let triangles: [([usize; 4], usize); 3] = [
        (face_vertices(&v, 4), m12),
        (face_vertices(&v, 2), m13),
        (face_vertices(&v, 0), m23),
    ];
    let mut starts = Vec::new();
    let mut started: Vec<Option<usize>> = vec![None; 3];
    for (t, &(face, flat)) in triangles.iter().enumerate() {
        let k = face_key(face);
        let users: Vec<usize> = (0..file.blocks.len())
            .filter(|&o| o != b)
            .filter(|&o| (0..6).any(|g| face_key(face_vertices(&file.blocks[o].vertices, g)) == k))
            .collect();
        let Some(&o) = users.first() else { continue };
        let ov = file.blocks[o].vertices;
        let slot = ov
            .iter()
            .position(|&x| x == flat)
            .expect("the face holds its flat vertex");
        let along = (0..3).find_map(|axis| {
            let mut cc = corner_of(slot);
            cc[axis] = 1 - cc[axis];
            let w = ov[corner_slot(cc[0], cc[1], cc[2])];
            (!face.contains(&w)).then_some(w)
        });
        started[t] = Some(starts.len());
        starts.push(Start {
            block: o,
            vertex: flat,
            along,
        });
    }
    // the faces shared with neighbours, which become interfaces when the
    // join is non-conformal
    let shared: Vec<usize> = (0..3).filter(|&t| started[t].is_some()).collect();
    let exclude: BTreeSet<usize> = [b].into_iter().collect();
    let (prisms, m0) = if non_conformal {
        let mut new_edges = Vec::new();
        let mut m0 = [0usize; 3];
        for (i, (q, n)) in [(p1, a), (p2, bb), (p3, c)].into_iter().enumerate() {
            m0[i] =
                crate::block_triangle::split_edge(file, q, p0, n, n0, &mut new_edges, change)?.0;
        }
        file.edges.extend(new_edges);
        started = vec![None; 3];
        (None, m0)
    } else {
        let prisms = split_prisms(
            file,
            &starts,
            Some(n0),
            stop,
            &exclude,
            &[(p1, p0, a), (p2, p0, bb), (p3, p0, c)],
            change,
        )
        .map_err(|e| {
            refused(format!(
                "the new row of tet block {b} cannot run through its neighbours: {}",
                e.message
            ))
        })?;
        let m0 = [p1, p2, p3].map(|q| prisms.row[&key(p0, q)].vertex);
        (Some(prisms), m0)
    };
    let [m01, m02, m03] = m0;
    // each triangle's centre: the neighbour column's centre, else the mean
    // of its three edge points (the new vertices and the flat one), which
    // lies inside the triangle they span, so the three quadrilaterals are
    // convex; the tet's centre at the mean of the four face centres
    let q = |i: usize| v3(file.vertices[i]);
    let own = [
        (q(m01) + q(m02) + q(m12)) / 3.0,
        (q(m01) + q(m03) + q(m13)) / 3.0,
        (q(m02) + q(m03) + q(m23)) / 3.0,
    ];
    let mut centres = [0usize; 3];
    for t in 0..3 {
        centres[t] = match started[t] {
            Some(i) => match &prisms {
                Some(p) => p.columns[i].centres[0],
                None => add(file, own[t], change),
            },
            None => add(file, own[t], change),
        };
    }
    let [f3, f2, f1] = centres;
    let cc = {
        let q = |i: usize| v3(file.vertices[i]);
        (q(f1) + q(f2) + q(f3) + q(f0)) / 4.0
    };
    let cen = add(file, cc, change);
    let domain = r.domain.clone();
    let hex = |vertices: [usize; 8], cells: [usize; 3]| BlockFile {
        vertices,
        cells,
        domain: domain.clone(),
    };
    let h = [
        hex([p0, m01, f3, m02, m03, f2, cen, f1], [a, bb, c]),
        hex([p1, m01, f2, m13, m12, f3, cen, f0], [n0, c, bb]),
        hex([p2, m02, f3, m12, m23, f1, cen, f0], [n0, a, c]),
        hex([p3, m03, f1, m23, m13, f2, cen, f0], [n0, bb, a]),
    ];
    // every old face → the new faces inside it
    let centre_of = |face: &[usize; 4]| -> Option<usize> {
        if face.contains(&p0) {
            let others: BTreeSet<usize> = face.iter().copied().collect();
            Some(if !others.contains(&p3) {
                f3
            } else if !others.contains(&p2) {
                f2
            } else {
                f1
            })
        } else {
            None
        }
    };
    let mut faces_split: Vec<([usize; 4], Vec<[usize; 4]>)> = Vec::new();
    for f in 0..6 {
        let face = face_vertices(&v, f);
        let mut on: BTreeSet<usize> = face.iter().copied().collect();
        for (q, m) in [(p1, m01), (p2, m02), (p3, m03)] {
            if on.contains(&p0) && on.contains(&q) {
                on.insert(m);
            }
        }
        if let Some(fc) = centre_of(&face) {
            on.insert(fc);
        }
        let parts: Vec<[usize; 4]> = h
            .iter()
            .flat_map(|x| (0..6).map(|g| face_vertices(&x.vertices, g)))
            .filter(|fv| fv.iter().all(|x| on.contains(x)))
            .collect();
        faces_split.push((face_key(face), parts));
    }
    let [h0, h1, h2, h3] = h;
    file.blocks[b] = h0;
    file.blocks.extend([h1, h2, h3]);
    let n = file.blocks.len();
    change.blocks.push(b);
    change.added_blocks.extend([n - 3, n - 2, n - 1]);
    let replace = |faces: &mut Vec<[usize; 4]>| {
        *faces = faces
            .iter()
            .flat_map(|f| {
                let k = face_key(*f);
                match faces_split.iter().find(|(o, _)| *o == k) {
                    Some((_, parts)) => parts.clone(),
                    None => vec![*f],
                }
            })
            .collect();
    };
    for p in &mut file.patches {
        replace(&mut p.faces);
    }
    for s in &mut file.surfaces {
        replace(&mut s.faces);
    }
    // non-conformal: each shared triangle becomes an interface, the
    // neighbour's whole face coarse, the tet's three parts fine
    let mut interfaces = Vec::new();
    if non_conformal {
        for &t in &shared {
            let face = triangles[t].0;
            let parts = faces_split
                .iter()
                .find(|(k, _)| *k == face_key(face))
                .map(|(_, parts)| parts.clone())
                .unwrap_or_default();
            let n = file.interfaces.len();
            let (coarse, fine) = (format!("interface{n}_coarse"), format!("interface{n}_fine"));
            for (name, faces) in [(&coarse, vec![face]), (&fine, parts)] {
                if file.patches.iter().any(|p| p.name == *name) {
                    return Err(refused(format!("a patch is already named {name}")));
                }
                file.patches.push(crate::block::PatchFile {
                    name: name.clone(),
                    role: Some(mesh_io::handoff::PatchRole::Interface),
                    faces,
                });
                change.patches.push(name.clone());
            }
            file.interfaces.push(crate::block::InterfaceFile {
                coarse: coarse.clone(),
                fine: fine.clone(),
            });
            change.interfaces.push(coarse.clone());
            interfaces.push(coarse);
        }
    }
    let counts = match (stop, &prisms) {
        (RowStop::TJunction, Some(p)) => Some(t_junctions(file, &p.row, cells_before)?),
        _ => None,
    };
    Ok(TetSplit {
        block: b,
        corner: p0,
        corner_deg: deg,
        cells: n0,
        blocks: [b, n - 3, n - 2, n - 1],
        neighbours: prisms
            .as_ref()
            .map(|p| {
                p.columns
                    .iter()
                    .flat_map(|c| c.column.iter().copied())
                    .collect()
            })
            .unwrap_or_default(),
        row_blocks: prisms.map(|p| p.row_blocks).unwrap_or_default(),
        stop,
        counts,
        join: if non_conformal {
            "non_conformal"
        } else {
            "conformal"
        },
        interfaces,
        fallback: None,
        cells_before,
        cells_after: file
            .blocks
            .iter()
            .map(|x| x.cells.iter().product::<usize>())
            .sum(),
    })
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::block::{PatchFile, TOPOLOGY_FORMAT, TOPOLOGY_FORMAT_VERSION, check};
    use crate::block_edit::{Edit, apply};
    use mesh_core::block::BLOCK_FACES;

    /// A tetrahedron P0 = origin, P1 = x, P2 = y, P3 = z as one hex block;
    /// the flat vertices and v6 are pushed 0.04 outward in their faces.
    pub(crate) fn tet() -> TopologyFile {
        let s = 0.04 / 2f64.sqrt();
        let t = 0.04 / 3f64.sqrt();
        let vertices = vec![
            [0.0, 0.0, 0.0],                               // 0 P0
            [1.0, 0.0, 0.0],                               // 1 P1
            [0.5 + s, 0.5 + s, 0.0],                       // 2 m12
            [0.0, 1.0, 0.0],                               // 3 P2
            [0.0, 0.0, 1.0],                               // 4 P3
            [0.5 + s, 0.0, 0.5 + s],                       // 5 m13
            [1.0 / 3.0 + t, 1.0 / 3.0 + t, 1.0 / 3.0 + t], // 6 F0
            [0.0, 0.5 + s, 0.5 + s],                       // 7 m23
        ];
        let block = BlockFile {
            vertices: [0, 1, 2, 3, 4, 5, 6, 7],
            cells: [4, 5, 6],
            domain: "fluid".into(),
        };
        let faces: Vec<[usize; 4]> = (0..6).map(|f| face_vertices(&block.vertices, f)).collect();
        TopologyFile {
            format: TOPOLOGY_FORMAT.into(),
            format_version: TOPOLOGY_FORMAT_VERSION.into(),
            vertices,
            edges: Vec::new(),
            blocks: vec![block],
            patches: vec![
                PatchFile {
                    name: "wall".into(),
                    role: None,
                    faces: (0..6)
                        .filter(|&g| BLOCK_FACES[g].1 == 0)
                        .map(|g| faces[g])
                        .collect(),
                },
                PatchFile {
                    name: "outlet".into(),
                    role: None,
                    faces: (0..6)
                        .filter(|&g| BLOCK_FACES[g].1 == 1)
                        .map(|g| faces[g])
                        .collect(),
                },
            ],
            surfaces: Vec::new(),
            junctions: Vec::new(),
            interfaces: Vec::new(),
            reference: None,
        }
    }

    fn volume(file: &TopologyFile) -> f64 {
        let g = algo::build(&file.topology().unwrap()).unwrap();
        g.volumes().iter().sum::<f64>()
    }

    fn quality(file: &TopologyFile) -> algo::quality3d::GridQuality {
        let g = algo::build(&file.topology().unwrap()).unwrap();
        algo::quality3d::measure(&g).0
    }

    /// The tet is found in every one of its 24 frames, a cube is not.
    #[test]
    fn a_tet_shaped_block_is_found() {
        let f = tet();
        let found = tet_blocks(&f, FLAT_TET_CORNER).unwrap();
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!(found[0].corner, 0);
        assert!(found[0].corner_deg > 165.0, "{found:?}");
        for (p, s) in rotations() {
            let mut g = f.clone();
            g.blocks[0] = rotated(&f.blocks[0], p, s);
            let found = tet_blocks(&g, FLAT_TET_CORNER).unwrap();
            assert_eq!(found.len(), 1);
            assert_eq!(found[0].corner, 0);
        }
        let cube = crate::block_triangle::tests::two_layers();
        assert!(tet_blocks(&cube, FLAT_TET_CORNER).unwrap().is_empty());
    }

    /// Invariant (BS-20, BS-05): the tet becomes four hexahedra with the
    /// counts of the corners (n0 new on the edges at P0), the faces on
    /// the patches are replaced by their parts, the volume stays, the
    /// result checks, and the folded flat corners are gone.
    #[test]
    fn a_tet_becomes_four_hexahedra() {
        let f = tet();
        let before = quality(&f);
        assert!(before.min_scaled_jacobian < 0.1, "{before:?}");
        let (g, c) = apply(
            &f,
            &Edit::SplitTet {
                block: 0,
                corner: None,
                cells: Some(3),
                stop: None,
                join: None,
            },
        )
        .unwrap();
        assert_eq!(g.blocks.len(), 4);
        assert_eq!(c.added_blocks, vec![1, 2, 3]);
        let cells: usize = g
            .blocks
            .iter()
            .map(|b| b.cells.iter().product::<usize>())
            .sum();
        assert_eq!(cells, 4 * 5 * 6 + 3 * (6 * 5 + 4 * 6 + 5 * 4));
        let checked = check(&g).unwrap();
        assert_eq!(checked["cells"].as_u64().unwrap() as usize, cells);
        let (vf, vg) = (volume(&f), volume(&g));
        assert!((vf - vg).abs() < 1e-3 * vf, "{vf} {vg}");
        let faces = |file: &TopologyFile, name: &str| {
            file.patches
                .iter()
                .find(|p| p.name == name)
                .unwrap()
                .faces
                .len()
        };
        assert_eq!(faces(&g, "wall"), 9);
        assert_eq!(faces(&g, "outlet"), 3);
        let after = quality(&g);
        assert_eq!(after.folded_corners, 0, "{after:?}");
        assert!(after.min_scaled_jacobian > 0.3, "{after:?}");
        // the default count is the mean of the tet's counts, which puts
        // the new vertices near the middles of the edges at P0
        let (g, _) = apply(
            &f,
            &Edit::SplitTet {
                block: 0,
                corner: None,
                cells: None,
                stop: None,
                join: None,
            },
        )
        .unwrap();
        assert_eq!(g.blocks[1].cells[0], 5);
        let after = quality(&g);
        assert_eq!(after.folded_corners, 0, "{after:?}");
        assert!(after.min_scaled_jacobian > 0.45, "{after:?}");
        assert!(after.max_deg < 55.0, "{after:?}");
    }

    /// The face P1P2P3 stays as it is, so a block on it keeps meeting the
    /// tet node for node.
    #[test]
    fn a_neighbour_on_the_fourth_face_stays_conformal() {
        let mut f = tet();
        // a block on the quadrilateral P1 m12 F0 m13 (old face i = 1),
        // pushed outward along the face normal
        let quad = face_vertices(&f.blocks[0].vertices, 1);
        let n = Vec3::new(1.0, 1.0, 1.0) / 3f64.sqrt() * 0.4;
        let base = f.vertices.len();
        for &v in &quad {
            f.vertices.push(a3(v3(f.vertices[v]) + n));
        }
        // the face's (j, k) corners as the new block's k = 0 face: j × k
        // points outward, so k runs outward too
        let nb = BlockFile {
            vertices: [
                quad[0],
                quad[1],
                quad[2],
                quad[3],
                base,
                base + 1,
                base + 2,
                base + 3,
            ],
            cells: [f.blocks[0].cells[1], f.blocks[0].cells[2], 3],
            domain: "fluid".into(),
        };
        let qk = face_key(quad);
        for p in &mut f.patches {
            p.faces.retain(|x| face_key(*x) != qk);
        }
        let outer: Vec<[usize; 4]> = (0..6)
            .map(|g| face_vertices(&nb.vertices, g))
            .filter(|x| face_key(*x) != qk)
            .collect();
        f.blocks.push(nb);
        f.patches[1].faces.extend(outer);
        let before = quality(&f);
        assert!(before.min_scaled_jacobian < 0.1, "{before:?}");
        check(&f).unwrap();
        let (g, _) = apply(
            &f,
            &Edit::SplitTet {
                block: 0,
                corner: None,
                cells: None,
                stop: None,
                join: None,
            },
        )
        .unwrap();
        check(&g).unwrap();
        assert_eq!(g.blocks.len(), 5);
        let after = quality(&g);
        assert!(after.min_scaled_jacobian > 0.4, "{after:?}");
    }

    /// The tet with a slab under its triangle P0 P1 m12 P2 (z = 0), the
    /// slab a prism along z.
    pub(crate) fn tet_on_slab() -> TopologyFile {
        let mut f = tet();
        let base = f.vertices.len();
        for i in [0, 1, 2, 3] {
            let p = f.vertices[i];
            f.vertices.push([p[0], p[1], -0.5]);
        }
        f.blocks.push(BlockFile {
            vertices: [base, base + 1, base + 2, base + 3, 0, 1, 2, 3],
            cells: [4, 5, 6],
            domain: "fluid".into(),
        });
        let k0 = face_key([0, 1, 2, 3]);
        for p in &mut f.patches {
            p.faces.retain(|x| face_key(*x) != k0);
        }
        let slab: Vec<[usize; 4]> = (0..6)
            .map(|g| face_vertices(&f.blocks[1].vertices, g))
            .filter(|x| face_key(*x) != k0)
            .collect();
        f.patches[1].faces.extend(slab);
        check(&f).unwrap();
        f
    }

    /// Invariant (BS-20, card 2026-11-27 "Conformal, prism columns"): a
    /// block on a triangle at P0 is split into three with the tet (a
    /// triangle split along its edge away from the tet, the new row's n0
    /// cells through it), the two meet node for node, nothing folds, and
    /// both row ends give the same here (the row reaches the boundary in
    /// the slab).
    #[test]
    fn a_neighbour_on_a_triangle_becomes_a_prism_column() {
        let f = tet_on_slab();
        for stop in [RowStop::Sheet, RowStop::TJunction] {
            let mut g = f.clone();
            let split = split_tet(
                &mut g,
                0,
                None,
                None,
                stop,
                Join::ConformalOnly,
                &mut Change::default(),
            )
            .unwrap();
            assert_eq!(split.neighbours, vec![1]);
            assert!(split.row_blocks.is_empty());
            assert_eq!(g.blocks.len(), 7);
            let checked = check(&g).unwrap();
            let cells: usize = g
                .blocks
                .iter()
                .map(|b| b.cells.iter().product::<usize>())
                .sum();
            assert_eq!(checked["cells"].as_u64().unwrap() as usize, cells);
            // the slab: 4 × 5 × 6 becomes the three quads of the triangle
            // (n1 = 4, n2 = 5, n0 = 5) times 6 layers
            let tet_cells = 4 * 5 * 6 + 5 * (6 * 5 + 4 * 6 + 5 * 4);
            assert_eq!(cells, tet_cells + 6 * (4 * 5 + 5 * 5 + 5 * 4));
            let (vf, vg) = (volume(&f), volume(&g));
            assert!((vf - vg).abs() < 1e-3 * vf, "{vf} {vg}");
            let after = quality(&g);
            assert_eq!(after.folded_corners, 0, "{after:?}");
            assert!(after.min_scaled_jacobian > 0.4, "{after:?}");
        }
    }

    /// The tet on a box whose bottom is a real quadrilateral (no flat
    /// corner), so the box cannot take a prism column.
    pub(crate) fn tet_on_box() -> TopologyFile {
        let mut f = tet_on_slab();
        let base = f.vertices.len() - 4;
        for (i, p) in [[-0.2, -0.2], [1.2, -0.2], [1.2, 1.2], [-0.2, 1.2]]
            .iter()
            .enumerate()
        {
            f.vertices[base + i] = [p[0], p[1], -0.5];
        }
        check(&f).unwrap();
        f
    }

    /// Invariant (BS-16, card A of step 3 "Automatic fallback"): when the
    /// neighbour on a triangle cannot take a prism column, the automatic
    /// join leaves it whole and makes the shared face a non-conformal
    /// interface; `conformal_only` refuses. The joined hand-off mesh closes
    /// every cell, has no inverted face, conserves the volume and the
    /// interface area, and only the cells at the interface are polyhedra.
    #[test]
    fn a_box_neighbour_joins_non_conformally() {
        let f = tet_on_box();
        let mut refused_copy = f.clone();
        assert!(
            split_tet(
                &mut refused_copy,
                0,
                None,
                None,
                RowStop::Sheet,
                Join::ConformalOnly,
                &mut Change::default()
            )
            .is_err()
        );
        let mut g = f.clone();
        let split = split_tet(
            &mut g,
            0,
            None,
            None,
            RowStop::Sheet,
            Join::Auto,
            &mut Change::default(),
        )
        .unwrap();
        assert_eq!(split.join, "non_conformal");
        assert!(split.fallback.is_some());
        assert_eq!(split.interfaces.len(), 1);
        // cards C and E: the coarsest n0 from ceil(max(4, 5, 6) / 2) whose
        // tet cells reach a scaled Jacobian of 0.4
        assert_eq!(split.cells, 4);
        assert_eq!(g.blocks.len(), 5);
        assert_eq!(g.blocks[1].cells, [4, 5, 6], "the box stays as it was");
        check(&g).unwrap();
        let (grid, _) = crate::block::grid(&g, &crate::block::VolumeAlgorithm::Tfi).unwrap();
        let (report, handoff, _) = crate::block::validate(&g, &grid).unwrap();
        let joins = report["interfaces"].as_array().unwrap();
        assert_eq!(joins.len(), 1);
        let j = &joins[0];
        assert_eq!(j["coarse_faces"].as_u64(), Some(4 * 5));
        // snapping moves points by at most 2 % of a local edge, so a cell
        // face's overlay pieces differ from it by about that much
        assert!(j["max_area_error"].as_f64().unwrap() < 0.02, "{j}");
        assert!(j["max_gap"].as_f64().unwrap() < 1e-12, "{j}");
        let summary = handoff.summary();
        assert_eq!(summary.inverted_faces, 0);
        let hex_volume: f64 = grid.volumes().iter().sum();
        let poly_volume: f64 = handoff.cell_volumes().iter().map(|x| x.0).sum();
        assert!(
            (poly_volume / hex_volume - 1.0).abs() < 1e-9,
            "{poly_volume} {hex_volume}"
        );
        let poly = report["polyhedra"].as_u64().unwrap() as usize;
        // the box's top layer and the tet cells on the triangle, at most
        assert!(poly > 0 && poly <= 4 * 5 + 3 * 5 * 6 + 3 * 3 * 3, "{poly}");
        assert!(
            handoff
                .patches
                .iter()
                .all(|p| !p.name.starts_with("interface"))
        );
    }

    /// The new row runs on through a block that shares only an edge at P0
    /// with the tet: `sheet` cuts it in two, `t_junction` keeps it whole
    /// with a T-junction and fits the counts.
    #[test]
    fn the_row_runs_on_through_a_block_on_an_edge_at_p0() {
        let mut f = tet();
        // a box at y < 0, z < 0 whose edge at y = z = 0 is P0–P1: it
        // touches the tet only there
        let base = f.vertices.len();
        for p in [
            [0.0, -0.6, -0.5],
            [1.0, -0.6, -0.5],
            [1.0, 0.0, -0.5],
            [0.0, 0.0, -0.5],
            [0.0, -0.6, 0.0],
            [1.0, -0.6, 0.0],
        ] {
            f.vertices.push(p);
        }
        f.blocks.push(BlockFile {
            vertices: [base, base + 1, base + 2, base + 3, base + 4, base + 5, 1, 0],
            cells: [4, 3, 2],
            domain: "fluid".into(),
        });
        let faces: Vec<[usize; 4]> = (0..6)
            .map(|g| face_vertices(&f.blocks[1].vertices, g))
            .collect();
        f.patches[1].faces.extend(faces);
        let sheet = {
            let mut g = f.clone();
            check(&g).unwrap();
            let s = split_tet(
                &mut g,
                0,
                None,
                None,
                RowStop::Sheet,
                Join::ConformalOnly,
                &mut Change::default(),
            )
            .unwrap();
            assert_eq!(s.row_blocks, vec![1]);
            check(&g).unwrap();
            assert_eq!(g.blocks.len(), 4 + 2);
            g
        };
        assert_eq!(quality(&sheet).folded_corners, 0);
        let mut g = f.clone();
        let s = split_tet(
            &mut g,
            0,
            None,
            None,
            RowStop::TJunction,
            Join::ConformalOnly,
            &mut Change::default(),
        )
        .unwrap();
        check(&g).unwrap();
        assert_eq!(g.blocks.len(), 4 + 1);
        assert_eq!(g.junctions.len(), 1, "{:?}", g.junctions);
        let counts = s.counts.unwrap();
        assert_eq!(counts.junctions, 1);
        assert_eq!(quality(&g).folded_corners, 0);
    }
}
