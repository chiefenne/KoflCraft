//! Triangle split of a block (BS-03, MB1h): a block whose cross-section is
//! a triangle (one corner flat, near 180°) becomes three blocks around a
//! centre vertex of valence 3, the "tri-center" map of a triangle by three
//! quadrilaterals. On an extruded topology the three blocks are the three
//! hexahedra of a prism.
//!
//! **Cell counts.** Give each corner i of the triangle a count n_i. The
//! quadrilateral at corner i has n_j cells along its side toward corner j,
//! so side (i, j) has n_i + n_j cells, the three side counts sum to an even
//! number, and each is below the sum of the other two. In a block whose
//! flat corner M lies on side AB, the counts of the two edges at M give
//! n_B (A–M) and n_A (M–B); the opposite edges B–C and C–A keep them. The
//! third count n_C is new: its cells enter as a row that leaves the
//! triangle through side BC next to B and side CA next to A, and runs
//! through the neighbouring blocks to the boundary, the way
//! `block.split_block` runs a split through a chain. Every block on the
//! row's way is cut in two, the new part with n_C cells; the new vertices
//! sit on the edges' own curves and laws with n_C more cells.

use std::collections::{BTreeMap, BTreeSet};

use mesh_algorithms::block as algo;
use mesh_algorithms::block::layout::{
    self, BlockOrthogonality, LayoutOptions, LayoutReport, Orthogonality,
};
use mesh_core::block::{BLOCK_FACES, Curve, Distribution, EdgeEnd, block_edges, corner_slot, face_key};
use mesh_core::geometry::Vec3;
use serde::Serialize;

use crate::EngineError;
use crate::block::{BlockFile, CurveFile, DistributionFile, EdgeFile, TopologyFile};
use crate::block_edit::{
    Change, a3, cut_spline, face_vertices, fitted_curve, fitted_law, key, refused, reversed_curve,
    reversed_law, v3,
};

/// A corner flatter than this (degrees between its two in-plane edges)
/// marks a block as a triangle.
pub const FLAT_CORNER: f64 = 150.0;
/// The least corner angle a split by hand accepts.
const MIN_CORNER: f64 = 120.0;

fn core_curve(c: &CurveFile) -> Curve {
    match c {
        CurveFile::Line => Curve::Line,
        CurveFile::Arc { through } => Curve::Arc {
            through: v3(*through),
        },
        CurveFile::Spline { through } => Curve::Spline {
            through: through.iter().map(|&p| v3(p)).collect(),
        },
    }
}

fn core_law(d: &DistributionFile) -> Distribution {
    match *d {
        DistributionFile::Uniform => Distribution::Uniform,
        DistributionFile::Geometric {
            first: Some(spacing),
            ..
        } => Distribution::Geometric {
            spacing,
            at: EdgeEnd::Start,
        },
        DistributionFile::Geometric {
            last: Some(spacing),
            ..
        } => Distribution::Geometric {
            spacing,
            at: EdgeEnd::End,
        },
        DistributionFile::Geometric { .. } => Distribution::Uniform,
        DistributionFile::Tanh { first, last } => Distribution::Tanh { first, last },
    }
}

/// The curve and law of the edge from `p` to `q`, read in that direction.
fn edge_of(file: &TopologyFile, p: usize, q: usize) -> (CurveFile, DistributionFile) {
    match file
        .edges
        .iter()
        .find(|e| key(e.vertices[0], e.vertices[1]) == key(p, q))
    {
        Some(e) if e.vertices[0] == p => (e.curve.clone(), e.distribution.clone()),
        Some(e) => (reversed_curve(&e.curve), reversed_law(&e.distribution)),
        None => (CurveFile::Line, DistributionFile::Uniform),
    }
}

/// Fractions of the length of the edge `p` → `q` at its `n` nodes.
fn edge_fractions(
    file: &TopologyFile,
    p: usize,
    q: usize,
    n: usize,
) -> Result<(algo::distribution::Evaluator, Vec<f64>), EngineError> {
    let (curve, law) = edge_of(file, p, q);
    let ev = algo::distribution::Evaluator::new(
        v3(file.vertices[p]),
        v3(file.vertices[q]),
        &core_curve(&curve),
    )
    .map_err(refused)?;
    let f = algo::distribution::fractions(core_law(&law), n, ev.length()).map_err(refused)?;
    Ok((ev, f))
}

/// The local corner (i, j, k) of a vertex slot.
fn corner_of(slot: usize) -> [usize; 3] {
    let (i, j) = [(0, 0), (1, 0), (1, 1), (0, 1)][slot % 4];
    [i, j, slot / 4]
}

/// A split edge of the new cell row: the vertex on the side of the new
/// cells, the old cell count, and (once placed) the new vertex with the
/// edge's nodes from the new cells' side.
struct RowEdge {
    from: usize,
    to: usize,
    cells: usize,
    vertex: usize,
    nodes: Vec<Vec3>,
    /// Fraction of the length from `from` at the new vertex, and the old
    /// nodes' fractions from `from`.
    at: f64,
    old: Vec<f64>,
}

/// What a triangle split did.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct TriangleSplit {
    pub block: usize,
    pub vertex: usize,
    /// Degrees between the two edges at the flat corner (first grid cell).
    pub corner_deg: f64,
    /// Cells of the new row (n_C) and the blocks it cut.
    pub row_cells: usize,
    pub row_blocks: Vec<usize>,
    /// The three blocks of the triangle.
    pub blocks: [usize; 3],
}

/// Split block `b0` at its flat corner `vertex` (BS-03); `cells` is the
/// count of the new row, chosen from the triangle's spacing when absent.
pub(crate) fn split_triangle(
    file: &mut TopologyFile,
    b0: usize,
    vertex: usize,
    cells: Option<usize>,
    change: &mut Change,
) -> Result<TriangleSplit, EngineError> {
    if b0 >= file.blocks.len() {
        return Err(refused(format!("there is no block {b0}")));
    }
    let old = file.blocks[b0].clone();
    let Some(slot) = old.vertices.iter().position(|&v| v == vertex) else {
        return Err(refused(format!("vertex {vertex} is no corner of block {b0}")));
    };
    let f = slot % 4;
    for side in [4, 5] {
        let k = face_key(face_vertices(&old.vertices, side));
        let users = file
            .blocks
            .iter()
            .filter(|b| (0..6).any(|g| face_key(face_vertices(&b.vertices, g)) == k))
            .count();
        if users > 1 {
            return Err(refused(format!(
                "block {b0} has a neighbour across its k = {} face; the triangle split works on one \
                 layer of blocks (an extruded topology)",
                side - 4
            )));
        }
    }
    let grid = algo::build(&file.topology()?).map_err(issues)?;
    let corner_deg = corner_angles(&grid, file, b0)[f];
    if corner_deg < MIN_CORNER {
        return Err(refused(format!(
            "the corner of block {b0} at vertex {vertex} is {corner_deg:.1}°; a triangle split \
             needs a corner of at least {MIN_CORNER}° there (a block whose cross-section is a \
             triangle)"
        )));
    }
    // the triangle A, B, C with the flat corner M on side AB, both layers
    let (ia, im, ib, ic) = ((f + 3) % 4, f, (f + 1) % 4, (f + 2) % 4);
    let v = old.vertices;
    let lay = |i: usize, k: usize| v[i + 4 * k];
    let (n_b, n_a) = (old.cells[ia % 2], old.cells[im % 2]);
    let n_c = match cells {
        Some(0) => return Err(refused("the new row needs at least one cell")),
        Some(c) => c,
        None => {
            let p = |i: usize| v3(file.vertices[lay(i, 0)]);
            let centre = (p(ia) + p(ib) + p(ic)) / 3.0;
            let h = ((p(im) - p(ia)).length() / n_b as f64
                + (p(ib) - p(im)).length() / n_a as f64)
                / 2.0;
            (((centre - p(im)).length() / h).round() as usize).max(1)
        }
    };
    // the row: edges B–C (new cells next to B) and C–A (next to A), both
    // layers, then every block that shares a split edge
    let mut row: BTreeMap<(usize, usize), RowEdge> = BTreeMap::new();
    let start = |from: usize, to: usize, cells: usize| RowEdge {
        from,
        to,
        cells,
        vertex: usize::MAX,
        nodes: Vec::new(),
        at: 0.0,
        old: Vec::new(),
    };
    for k in 0..2 {
        row.insert(
            key(lay(ib, k), lay(ic, k)),
            start(lay(ib, k), lay(ic, k), n_b),
        );
        row.insert(
            key(lay(ia, k), lay(ic, k)),
            start(lay(ia, k), lay(ic, k), n_a),
        );
    }
    let own: BTreeSet<(usize, usize)> = row.keys().copied().collect();
    // crossed block → (axis, the new cells at its local start)
    let mut crossed: BTreeMap<usize, (usize, bool)> = BTreeMap::new();
    let along = |b: &BlockFile, a: usize| -> Vec<(usize, usize)> {
        block_edges()
            .into_iter()
            .filter(|(ax, _)| *ax == a)
            .map(|(_, c)| {
                let mut e = c;
                e[a] = 1;
                (
                    b.vertices[corner_slot(c[0], c[1], c[2])],
                    b.vertices[corner_slot(e[0], e[1], e[2])],
                )
            })
            .collect()
    };
    loop {
        let mut grew = false;
        for (o, b) in file.blocks.iter().enumerate() {
            if o == b0 || crossed.contains_key(&o) {
                continue;
            }
            for a in 0..3 {
                let edges = along(b, a);
                let hit: Vec<bool> = edges
                    .iter()
                    .filter_map(|&(p, q)| row.get(&key(p, q)).map(|r| r.from == p))
                    .collect();
                if hit.is_empty() {
                    continue;
                }
                let at_start = hit[0];
                if hit.iter().any(|&h| h != at_start) {
                    return Err(refused(format!(
                        "the new row would enter block {o} from both ends; split by hand"
                    )));
                }
                if let Some(&(a2, _)) = crossed.get(&o) {
                    if a2 != a {
                        return Err(refused(format!(
                            "the new row would cross block {o} twice; split by hand"
                        )));
                    }
                }
                crossed.insert(o, (a, at_start));
                for (p, q) in edges {
                    let (from, to) = if at_start { (p, q) } else { (q, p) };
                    row.entry(key(p, q))
                        .or_insert_with(|| start(from, to, b.cells[a]));
                }
                grew = true;
            }
        }
        if !grew {
            break;
        }
    }
    for (&o, &(a, _)) in &crossed {
        for g in (0..3).filter(|&g| g != a) {
            if along(&file.blocks[o], g)
                .iter()
                .any(|&(p, q)| row.contains_key(&key(p, q)))
            {
                return Err(refused(format!(
                    "the new row would cross block {o} twice; split by hand"
                )));
            }
        }
    }
    for g in 0..3 {
        for (p, q) in along(&old, g) {
            if row.contains_key(&key(p, q)) && !own.contains(&key(p, q)) {
                return Err(refused(format!(
                    "the new row would come back into block {b0}; split by hand"
                )));
            }
        }
    }
    for &o in crossed.keys() {
        let keys: BTreeSet<[usize; 4]> = (0..6)
            .map(|g| face_key(face_vertices(&file.blocks[o].vertices, g)))
            .collect();
        if file
            .surfaces
            .iter()
            .any(|s| s.faces.iter().any(|f| keys.contains(&face_key(*f))))
        {
            return Err(refused(format!(
                "the new row would cut block {o}, which has a face on a projection surface"
            )));
        }
    }
    // place the new vertices on the edges with n_C more cells
    let mut new_edges: Vec<EdgeFile> = Vec::new();
    for r in row.values_mut() {
        r.old = edge_fractions(file, r.from, r.to, r.cells)?.1;
        let (ev, fr) = edge_fractions(file, r.from, r.to, r.cells + n_c)?;
        let mut nodes: Vec<Vec3> = fr.iter().map(|&s| ev.at(s)).collect();
        nodes[0] = v3(file.vertices[r.from]);
        nodes[r.cells + n_c] = v3(file.vertices[r.to]);
        r.at = fr[n_c];
        file.vertices.push(a3(nodes[n_c]));
        r.vertex = file.vertices.len() - 1;
        change.added_vertices.push(r.vertex);
        let (curve, law) = edge_of(file, r.from, r.to);
        let halves = match &curve {
            CurveFile::Spline { through } => {
                let pts: Vec<Vec3> = through.iter().map(|&p| v3(p)).collect();
                Some(cut_spline(nodes[0], &pts, nodes[r.cells + n_c], nodes[n_c]))
            }
            _ => None,
        };
        for (h, (a, b, part)) in [
            (r.from, r.vertex, &nodes[..=n_c]),
            (r.vertex, r.to, &nodes[n_c..]),
        ]
        .into_iter()
        .enumerate()
        {
            let c = match (&curve, &halves) {
                (CurveFile::Line, _) => CurveFile::Line,
                (_, Some((x, y))) => {
                    let pts = if h == 0 { x } else { y };
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
            let d = if law == DistributionFile::Uniform {
                DistributionFile::Uniform
            } else {
                fitted_law(part, &c, &law)
            };
            if c != CurveFile::Line || d != DistributionFile::Uniform {
                new_edges.push(EdgeFile {
                    vertices: [a, b],
                    curve: c,
                    distribution: d,
                });
                change.fitted_edges.push([a, b]);
            }
        }
        r.nodes = nodes;
    }
    file.edges.retain(|e| !row.contains_key(&key(e.vertices[0], e.vertices[1])));
    // cut every crossed block in two along the row
    let mut faces_split: BTreeMap<[usize; 4], Vec<[usize; 4]>> = BTreeMap::new();
    let mut cut_edges: BTreeSet<(usize, usize)> = BTreeSet::new();
    for (&o, &(a, at_start)) in &crossed {
        let ob = file.blocks[o].clone();
        let n = ob.cells;
        let mut low = ob.clone();
        let mut high = ob.clone();
        // fractional grid index along `a` of each new vertex, from local 0
        let mut index: BTreeMap<usize, f64> = BTreeMap::new();
        for (_, c) in block_edges().into_iter().filter(|(ax, _)| *ax == a) {
            let mut e = c;
            e[a] = 1;
            let (sa, sb) = (corner_slot(c[0], c[1], c[2]), corner_slot(e[0], e[1], e[2]));
            let (p, q) = (ob.vertices[sa], ob.vertices[sb]);
            let r = &row[&key(p, q)];
            low.vertices[sb] = r.vertex;
            high.vertices[sa] = r.vertex;
            // the old cell the new vertex falls in, counted from `from`
            let k = r.old.partition_point(|&s| s <= r.at).clamp(1, r.cells) - 1;
            let t = ((r.at - r.old[k]) / (r.old[k + 1] - r.old[k])).clamp(0.0, 1.0);
            let from_start = k as f64 + t;
            let local = if r.from == p {
                from_start
            } else {
                r.cells as f64 - from_start
            };
            index.insert(r.vertex, local);
        }
        if at_start {
            low.cells[a] = n_c;
            high.cells[a] = n[a];
        } else {
            low.cells[a] = n[a];
            high.cells[a] = n_c;
        }
        for (g, face) in BLOCK_FACES.iter().enumerate() {
            if face.0 != a {
                faces_split.insert(
                    face_key(face_vertices(&ob.vertices, g)),
                    vec![face_vertices(&low.vertices, g), face_vertices(&high.vertices, g)],
                );
            }
        }
        // the edges of the cut face, through the old grid at the row
        let sample = |idx: [f64; 3]| -> Vec3 {
            let fl = idx[a].floor().min((n[a] - 1) as f64).max(0.0);
            let t = idx[a] - fl;
            let mut i0 = [idx[0] as usize, idx[1] as usize, idx[2] as usize];
            i0[a] = fl as usize;
            let mut i1 = i0;
            i1[a] += 1;
            let p0 = grid.points[grid.block_corner(o, i0[0], i0[1], i0[2])];
            let p1 = grid.points[grid.block_corner(o, i1[0], i1[1], i1[2])];
            p0 + (p1 - p0) * t
        };
        for (a2, c) in block_edges() {
            if a2 == a || c[a] != 0 {
                continue;
            }
            let mut e = c;
            e[a2] = 1;
            let (sa, sb) = (corner_slot(c[0], c[1], c[2]), corner_slot(e[0], e[1], e[2]));
            let (p, q) = (high.vertices[sa], high.vertices[sb]);
            if !cut_edges.insert(key(p, q)) {
                continue;
            }
            let (ip, iq) = (index[&p], index[&q]);
            let nodes: Vec<Vec3> = (0..=n[a2])
                .map(|s| {
                    let mut idx = [
                        (c[0] * n[0]) as f64,
                        (c[1] * n[1]) as f64,
                        (c[2] * n[2]) as f64,
                    ];
                    idx[a2] = s as f64;
                    idx[a] = ip + (iq - ip) * s as f64 / n[a2] as f64;
                    if s == 0 {
                        v3(file.vertices[p])
                    } else if s == n[a2] {
                        v3(file.vertices[q])
                    } else {
                        sample(idx)
                    }
                })
                .collect();
            let curve = fitted_curve(&nodes);
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
        file.blocks[o] = low;
        file.blocks.push(high);
        change.blocks.push(o);
        change.added_blocks.push(file.blocks.len() - 1);
    }
    // the triangle: three blocks around the centre O
    let m = |i: usize, j: usize, k: usize| row[&key(lay(i, k), lay(j, k))].vertex;
    let mut centre = [0usize; 2];
    for (k, c) in centre.iter_mut().enumerate() {
        let p = |i: usize| v3(file.vertices[lay(i, k)]);
        file.vertices.push(a3((p(ia) + p(ib) + p(ic)) / 3.0));
        *c = file.vertices.len() - 1;
        change.added_vertices.push(*c);
    }
    let quad = |q: [[usize; 2]; 4]| -> [usize; 8] {
        [
            q[0][0], q[1][0], q[2][0], q[3][0], q[0][1], q[1][1], q[2][1], q[3][1],
        ]
    };
    let at = |i: usize| [lay(i, 0), lay(i, 1)];
    let (mbc, mca) = ([m(ib, ic, 0), m(ib, ic, 1)], [m(ia, ic, 0), m(ia, ic, 1)]);
    let nk = old.cells[2];
    let qa = BlockFile {
        vertices: quad([at(ia), at(im), centre, mca]),
        cells: [n_b, n_c, nk],
        domain: old.domain.clone(),
    };
    let qb = BlockFile {
        vertices: quad([at(ib), mbc, centre, at(im)]),
        cells: [n_c, n_a, nk],
        domain: old.domain.clone(),
    };
    let qc = BlockFile {
        vertices: quad([at(ic), mca, centre, mbc]),
        cells: [n_a, n_b, nk],
        domain: old.domain.clone(),
    };
    let new_vertices: BTreeMap<(usize, usize), usize> =
        row.iter().map(|(&k, r)| (k, r.vertex)).collect();
    for g in 0..6 {
        let face = face_vertices(&old.vertices, g);
        let mut on: BTreeSet<usize> = face.iter().copied().collect();
        for w in 0..4 {
            if let Some(&nv) = new_vertices.get(&key(face[w], face[(w + 1) % 4])) {
                on.insert(nv);
            }
        }
        if BLOCK_FACES[g].0 == 2 {
            on.insert(centre[BLOCK_FACES[g].1]);
        }
        let parts: Vec<[usize; 4]> = [&qa, &qb, &qc]
            .iter()
            .flat_map(|q| (0..6).map(|h| face_vertices(&q.vertices, h)))
            .filter(|fv| fv.iter().all(|x| on.contains(x)))
            .collect();
        faces_split.insert(face_key(face), parts);
    }
    file.blocks[b0] = qa;
    file.blocks.push(qb);
    file.blocks.push(qc);
    let blocks = [b0, file.blocks.len() - 2, file.blocks.len() - 1];
    change.blocks.push(b0);
    change.added_blocks.extend([blocks[1], blocks[2]]);
    file.edges.extend(new_edges);
    let replace = |faces: &mut Vec<[usize; 4]>| {
        *faces = faces
            .iter()
            .flat_map(|f| match faces_split.get(&face_key(*f)) {
                Some(parts) => parts.clone(),
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
    Ok(TriangleSplit {
        block: b0,
        vertex,
        corner_deg,
        row_cells: n_c,
        row_blocks: crossed.keys().copied().collect(),
        blocks,
    })
}

fn issues(issues: Vec<algo::Issue>) -> EngineError {
    refused(
        issues
            .iter()
            .map(|i| i.message.clone())
            .collect::<Vec<_>>()
            .join("; "),
    )
}

/// The four in-plane corner angles (degrees) of block `b` at k = 0,
/// between its two edges, measured on the first cells of the TFI grid.
fn corner_angles(grid: &algo::BlockGrid, file: &TopologyFile, b: usize) -> [f64; 4] {
    let [ni, nj, _] = file.blocks[b].cells;
    let at = |i: usize, j: usize| grid.points[grid.block_corner(b, i, j, 0)];
    let mut out = [0.0; 4];
    for (c, o) in out.iter_mut().enumerate() {
        let [ci, cj, _] = corner_of(c);
        let (i, j) = (ci * ni, cj * nj);
        let di = if ci == 0 { i + 1 } else { i - 1 };
        let dj = if cj == 0 { j + 1 } else { j - 1 };
        let p = at(i, j);
        let (u, w) = (at(di, j) - p, at(i, dj) - p);
        *o = (u.dot(w) / (u.length() * w.length()))
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees();
    }
    out
}

/// A block corner flatter than `min_deg`: the block is a triangle.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct FlatCorner {
    pub block: usize,
    pub vertex: usize,
    pub corner_deg: f64,
}

/// The flat corners of the blocks of one layer (extruded topologies),
/// flattest first.
pub fn flat_corners(file: &TopologyFile, min_deg: f64) -> Result<Vec<FlatCorner>, EngineError> {
    let grid = algo::build(&file.topology()?).map_err(issues)?;
    let mut out = Vec::new();
    for (b, block) in file.blocks.iter().enumerate() {
        let angles = corner_angles(&grid, file, b);
        for (c, &deg) in angles.iter().enumerate() {
            if deg >= min_deg {
                out.push(FlatCorner {
                    block: b,
                    vertex: block.vertices[c],
                    corner_deg: deg,
                });
            }
        }
    }
    out.sort_by(|x, y| y.corner_deg.total_cmp(&x.corner_deg));
    Ok(out)
}

/// Most triangle splits the trigger adds.
pub const MAX_SPLITS: usize = 4;
/// Trials of the trigger per round: the flattest corners first.
const TRIALS_PER_ROUND: usize = 4;
/// A split must lower the worst angle of its triangle by this much.
const MIN_GAIN_DEG: f64 = 2.0;
/// The whole grid's worst angle may rise by this much at most.
const MAX_RISE_DEG: f64 = 0.05;
/// Growth (cell size ratio across a face) is fine up to this ...
const GROWTH_LIMIT: f64 = 2.0;
/// ... and may pass it, or the start's largest growth, by this factor.
const GROWTH_SLACK: f64 = 1.05;

/// A split the trigger kept and what it measured, in degrees of
/// non-orthogonality: the triangle block before against its three
/// blocks after, and the whole grid.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct KeptSplit {
    #[serde(flatten)]
    pub split: TriangleSplit,
    pub triangle_before_deg: f64,
    pub triangle_after_deg: f64,
    pub triangle_growth_before: f64,
    pub triangle_growth_after: f64,
    pub max_before_deg: f64,
    pub max_after_deg: f64,
}

/// A flat corner the trigger tried and left as it was.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct SkippedSplit {
    pub block: usize,
    pub vertex: usize,
    pub corner_deg: f64,
    pub reason: String,
}

/// What [`inject_triangle_splits`] did.
#[derive(Clone, Debug)]
pub struct Injection {
    pub file: TopologyFile,
    pub kept: Vec<KeptSplit>,
    pub skipped: Vec<SkippedSplit>,
    pub before: Orthogonality,
    pub after: Orthogonality,
}

/// The topology with its layout optimised ([`layout::optimize`]), the
/// blocks on wall patches held; patches and the reference as they were.
pub fn optimized(
    file: &TopologyFile,
    hold_growth: bool,
) -> Result<(TopologyFile, LayoutReport), EngineError> {
    let mut t = file.topology()?;
    let report = layout::optimize(
        &mut t,
        &LayoutOptions {
            rigid_patches: crate::block::wall_patches(file),
            hold_growth,
            ..LayoutOptions::default()
        },
    )
    .map_err(refused)?;
    let mut out = TopologyFile::from_topology(&t, file.reference.clone());
    out.patches = file.patches.clone();
    Ok((out, report))
}

/// Optimise the layout of a topology, with triangle splits where they help
/// when `splits` (`block.optimize_layout`): the new topology and the
/// layout report, also kept in its reference.
pub fn optimize_layout(
    file: &TopologyFile,
    splits: bool,
) -> Result<(TopologyFile, serde_json::Value), EngineError> {
    let (mut out, report) = optimized(file, false)?;
    let mut layout = crate::block::layout_report(&report);
    if splits {
        let injection = inject_triangle_splits(file, Some(&out), MAX_SPLITS)?;
        crate::block::note_injection(&mut layout, &injection);
        out = injection.file;
    }
    if let Some(r) = out.reference.as_mut() {
        r.layout = Some(layout.clone());
    }
    Ok((out, layout))
}

fn measure(file: &TopologyFile) -> Result<(Orthogonality, Vec<BlockOrthogonality>), EngineError> {
    let t = file.topology()?;
    Ok((
        layout::orthogonality(&t).map_err(refused)?,
        layout::block_orthogonality(&t).map_err(refused)?,
    ))
}

/// Add triangle splits where they measurably help (BS-03, MB1h). The
/// candidates are the block corners of at least [`FLAT_CORNER`] in `raw`,
/// the layout before its optimisation, flattest first. A trial splits one
/// there and optimises the layout ([`optimized`]); it is compared with
/// `optimised`, the layout optimised without the split (computed from
/// `raw` when absent). A trial is kept when the worst angle of the
/// triangle's blocks drops by [`MIN_GAIN_DEG`] or more while the whole
/// grid gets no worse: its worst angle rises by at most [`MAX_RISE_DEG`],
/// no more cells pass 70°, nothing folds and the growth stays within
/// bounds. The best trial of a round is kept and the next round starts
/// from it, up to `max_splits`.
pub fn inject_triangle_splits(
    raw: &TopologyFile,
    optimised: Option<&TopologyFile>,
    max_splits: usize,
) -> Result<Injection, EngineError> {
    let mut cur_raw = raw.clone();
    let mut cur = match optimised {
        Some(f) => f.clone(),
        None => optimized(raw, false)?.0,
    };
    let (mut q, mut per) = measure(&cur)?;
    let before = q;
    let mut kept: Vec<KeptSplit> = Vec::new();
    let mut skipped: Vec<SkippedSplit> = Vec::new();
    // a block keeps its vertices until something cuts it: a corner tried
    // once is not tried again
    let mut tried: BTreeSet<(usize, [usize; 8])> = BTreeSet::new();
    while kept.len() < max_splits {
        let fresh: Vec<FlatCorner> = flat_corners(&cur_raw, FLAT_CORNER)?
            .into_iter()
            .filter(|c| !tried.contains(&(c.vertex, cur_raw.blocks[c.block].vertices)))
            .take(TRIALS_PER_ROUND)
            .collect();
        if fresh.is_empty() {
            break;
        }
        let mut best: Option<(f64, TopologyFile, TopologyFile, KeptSplit)> = None;
        for c in fresh {
            tried.insert((c.vertex, cur_raw.blocks[c.block].vertices));
            match trial(&cur_raw, &c, q, &per) {
                Ok((gain, raw_split, f, k)) => {
                    if best.as_ref().is_none_or(|b| gain > b.0) {
                        best = Some((gain, raw_split, f, k));
                    }
                }
                Err(reason) => skipped.push(SkippedSplit {
                    block: c.block,
                    vertex: c.vertex,
                    corner_deg: c.corner_deg,
                    reason,
                }),
            }
        }
        if let Some((_, raw_split, f, k)) = best {
            (q, per) = measure(&f)?;
            cur_raw = raw_split;
            cur = f;
            kept.push(k);
        }
    }
    Ok(Injection {
        file: cur,
        kept,
        skipped,
        before,
        after: q,
    })
}

/// One trial of [`inject_triangle_splits`]: the gain, the split layout
/// before and after its optimisation and the record, or why it is left.
fn trial(
    raw: &TopologyFile,
    c: &FlatCorner,
    q: Orthogonality,
    per: &[BlockOrthogonality],
) -> Result<(f64, TopologyFile, TopologyFile, KeptSplit), String> {
    let mut split_raw = raw.clone();
    let split = split_triangle(
        &mut split_raw,
        c.block,
        c.vertex,
        None,
        &mut Change::default(),
    )
    .map_err(|e| e.message)?;
    // the split's own layout optimised so that it trades no new size jumps
    // for angles
    let (f, _) = optimized(&split_raw, true).map_err(|e| e.message)?;
    let (tq, tper) = measure(&f).map_err(|e| e.message)?;
    let region = |p: &[BlockOrthogonality], bs: &[usize]| {
        bs.iter().fold((0.0f64, 1.0f64), |(d, g), &b| {
            (d.max(p[b].max_deg), g.max(p[b].max_growth))
        })
    };
    let (deg0, grow0) = region(per, &[c.block]);
    let (deg1, grow1) = region(&tper, &split.blocks);
    if tq.folded_corners > 0 {
        return Err(format!("{} cell corners fold", tq.folded_corners));
    }
    if deg1 > deg0 - MIN_GAIN_DEG {
        return Err(format!(
            "the triangle's worst angle goes from {deg0:.1}° to {deg1:.1}°, less than \
             {MIN_GAIN_DEG}° better"
        ));
    }
    if tq.max_deg > q.max_deg + MAX_RISE_DEG {
        return Err(format!(
            "the grid's worst angle rises from {:.1}° to {:.1}°",
            q.max_deg, tq.max_deg
        ));
    }
    if tq.cells_over_70 > q.cells_over_70 {
        return Err(format!(
            "cells above 70° go from {} to {}",
            q.cells_over_70, tq.cells_over_70
        ));
    }
    if grow1 > grow0.max(GROWTH_LIMIT) * GROWTH_SLACK
        || tq.max_growth > q.max_growth.max(GROWTH_LIMIT) * GROWTH_SLACK
    {
        return Err(format!(
            "the growth rises: triangle {grow0:.2} to {grow1:.2}, grid {:.2} to {:.2}",
            q.max_growth, tq.max_growth
        ));
    }
    let k = KeptSplit {
        split,
        triangle_before_deg: deg0,
        triangle_after_deg: deg1,
        triangle_growth_before: grow0,
        triangle_growth_after: grow1,
        max_before_deg: q.max_deg,
        max_after_deg: tq.max_deg,
    };
    Ok((deg0 - deg1, split_raw, f, k))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::block::{PatchFile, TOPOLOGY_FORMAT, TOPOLOGY_FORMAT_VERSION, check};
    use crate::block_edit::{Edit, apply};

    /// Two blocks extruded along z: a triangle with a flat corner at
    /// vertex 1 (171°) and a neighbour across its side 2–3.
    pub(crate) fn wedge() -> TopologyFile {
        let plane = [
            [0.0, 0.0],
            [1.0, -0.08],
            [2.0, 0.0],
            [1.0, 1.7],
            [3.2, 0.0],
            [2.4, 1.8],
        ];
        let mut vertices = Vec::new();
        for z in [0.0, 0.1] {
            for p in plane {
                vertices.push([p[0], p[1], z]);
            }
        }
        let b = |q: [usize; 4]| [q[0], q[1], q[2], q[3], q[0] + 6, q[1] + 6, q[2] + 6, q[3] + 6];
        let blocks = vec![
            BlockFile {
                vertices: b([0, 1, 2, 3]),
                cells: [4, 4, 1],
                domain: "fluid".into(),
            },
            BlockFile {
                vertices: b([2, 4, 5, 3]),
                cells: [5, 4, 1],
                domain: "fluid".into(),
            },
        ];
        let side = |p: usize, q: usize| [p, q, q + 6, p + 6];
        let patches = vec![
            PatchFile {
                name: "wall".into(),
                role: None,
                faces: vec![side(0, 1), side(1, 2), side(3, 0)],
            },
            PatchFile {
                name: "outlet".into(),
                role: None,
                faces: vec![side(2, 4), side(4, 5), side(5, 3)],
            },
            PatchFile {
                name: "symmetry".into(),
                role: None,
                faces: blocks
                    .iter()
                    .flat_map(|x| {
                        let v = x.vertices;
                        [[v[0], v[1], v[2], v[3]], [v[4], v[5], v[6], v[7]]]
                    })
                    .collect(),
            },
        ];
        TopologyFile {
            format: TOPOLOGY_FORMAT.into(),
            format_version: TOPOLOGY_FORMAT_VERSION.into(),
            vertices,
            edges: Vec::new(),
            blocks,
            patches,
            surfaces: Vec::new(),
            reference: None,
        }
    }

    fn area(file: &TopologyFile) -> f64 {
        let g = algo::build(&file.topology().unwrap()).unwrap();
        g.volumes().iter().sum::<f64>() / 0.1
    }

    /// Invariant (BS-03, BS-05): the triangle becomes three blocks around
    /// a centre with the counts n_B, n_C, n_A of the tri-center map, the
    /// new row cuts the neighbour across side 2–3 in two (n_C cells next
    /// to the triangle's corner 2), the patches keep their faces, the area
    /// stays and the result checks.
    #[test]
    fn a_triangle_becomes_three_blocks_and_a_row() {
        let f = wedge();
        let found = flat_corners(&f, FLAT_CORNER).unwrap();
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!((found[0].block, found[0].vertex), (0, 1));
        let (g, c) = apply(
            &f,
            &Edit::SplitTriangle {
                block: 0,
                vertex: 1,
                cells: Some(3),
            },
        )
        .unwrap();
        assert_eq!(g.blocks.len(), 5);
        assert_eq!(c.added_blocks.len(), 3);
        let cells: usize = g.blocks.iter().map(|b| b.cells.iter().product::<usize>()).sum();
        assert_eq!(cells, (4 * 3 + 3 * 4 + 4 * 4) + 5 * (4 + 3));
        let checked = check(&g).unwrap();
        assert_eq!(checked["cells"].as_u64().unwrap() as usize, cells);
        assert!((area(&g) - area(&f)).abs() < 1e-12 * area(&f), "{} {}", area(&g), area(&f));
        let faces = |file: &TopologyFile, name: &str| {
            file.patches.iter().find(|p| p.name == name).unwrap().faces.len()
        };
        assert_eq!(faces(&g, "wall"), 4, "side 3–0 is cut by the row");
        assert_eq!(faces(&g, "outlet"), 4, "side 4–5 is cut by the row");
        assert_eq!(faces(&g, "symmetry"), 10);
        let grid = algo::build(&g.topology().unwrap()).unwrap();
        let jmin = grid
            .quality()
            .iter()
            .map(|q| q.min_jacobian)
            .fold(f64::INFINITY, f64::min);
        assert!(jmin > 0.0, "{jmin}");
    }

    /// The wedge with open boundaries (no walls): the optimiser may move
    /// every block.
    fn open_wedge() -> TopologyFile {
        let mut f = wedge();
        f.blocks[0].cells = [8, 8, 1];
        f.blocks[1].cells = [10, 8, 1];
        f.patches[0].name = "inlet".into();
        f
    }

    /// Invariant (MB1h): the trigger keeps a triangle split that lowers the
    /// grid's worst angle (here 62° to 26°) without a size jump above the
    /// growth limit, and records the corners it tried and left with the
    /// reason; the result checks and keeps the area. A trigger that keeps
    /// every split, judges before the optimisation or ignores the growth
    /// breaks it; the split test above does not run the trigger.
    #[test]
    fn the_trigger_keeps_a_split_that_helps() {
        let f = open_wedge();
        let (g, layout) = optimize_layout(&f, true).unwrap();
        let kept = layout["triangle_splits"]["kept"].as_array().unwrap();
        assert_eq!(kept.len(), 1, "{layout:#}");
        assert_eq!(kept[0]["block"], 0);
        assert_eq!(kept[0]["vertex"], 1);
        assert_eq!(kept[0]["row_blocks"], serde_json::json!([1]));
        let before = layout["before"]["max_nonorthogonality_deg"].as_f64().unwrap();
        let after = layout["after"]["max_nonorthogonality_deg"].as_f64().unwrap();
        assert!(before > 60.0 && after < 30.0, "{before} {after}");
        assert!(layout["after"]["max_growth"].as_f64().unwrap() <= GROWTH_LIMIT * GROWTH_SLACK);
        assert_eq!(layout["after"]["folded_corners"], 0);
        // the second round tried a corner of the new blocks and left it
        let skipped = layout["triangle_splits"]["skipped"].as_array().unwrap();
        assert!(!skipped.is_empty());
        assert!(skipped.iter().all(|s| !s["reason"].as_str().unwrap().is_empty()));
        assert_eq!(g.blocks.len(), 5);
        check(&g).unwrap();
        assert!((area(&g) - area(&f)).abs() < 1e-9 * area(&f));
    }

    /// Invariant (MB1h): without a corner of [`FLAT_CORNER`] the trigger
    /// changes nothing: the layout is the plain optimisation's and the
    /// report lists no split. A trigger that splits blocks without flat
    /// corners, or re-optimises a layout it keeps, breaks it.
    #[test]
    fn the_trigger_leaves_a_layout_without_triangles_alone() {
        let mut f = open_wedge();
        f.vertices[1] = [1.0, -0.6, 0.0];
        f.vertices[7] = [1.0, -0.6, 0.1];
        let (plain, _) = optimize_layout(&f, false).unwrap();
        let (g, layout) = optimize_layout(&f, true).unwrap();
        assert_eq!(layout["triangle_splits"]["kept"], serde_json::json!([]));
        assert_eq!(layout["triangle_splits"]["skipped"], serde_json::json!([]));
        assert_eq!(g, plain);
    }

    /// Invariant (MB1g, MB1h): with `hold_growth` the layout optimisation
    /// trades no new size jumps for angles; without it the split wedge's
    /// optimum has a cell size ratio above 3 across the row. Dropping the
    /// growth hold, or holding a folded start, breaks it; the trigger test
    /// only sees the held result.
    #[test]
    fn the_growth_hold_keeps_the_size_ratio() {
        let mut s = open_wedge();
        split_triangle(&mut s, 0, 1, None, &mut Change::default()).unwrap();
        let growth = |hold: bool| {
            let (o, _) = optimized(&s, hold).unwrap();
            let (q, _) = measure(&o).unwrap();
            (q.max_deg, q.max_growth)
        };
        let (free_deg, free_growth) = growth(false);
        let (held_deg, held_growth) = growth(true);
        assert!(free_growth > 3.0, "{free_growth}");
        assert!(held_growth <= GROWTH_LIMIT * GROWTH_SLACK, "{held_growth}");
        assert!(held_deg < free_deg + 1.0, "{held_deg} {free_deg}");
    }

    /// Invariant (BS-03): a split by hand is refused, with nothing
    /// changed, at a corner below 120° (a quadrilateral, no triangle) and
    /// in a block with a neighbour across its k faces (more than one layer).
    /// Dropping either check breaks it.
    #[test]
    fn a_split_is_refused_off_a_triangle_or_layer() {
        let f = wedge();
        let e = apply(
            &f,
            &Edit::SplitTriangle {
                block: 1,
                vertex: 4,
                cells: None,
            },
        )
        .unwrap_err();
        assert!(e.message.contains("needs a corner of at least"), "{}", e.message);
        // a second layer on top of the triangle
        let mut g = wedge();
        let top: Vec<[f64; 3]> = (6..12).map(|v| [g.vertices[v][0], g.vertices[v][1], 0.2]).collect();
        g.vertices.extend(top);
        let mut upper = g.blocks[0].clone();
        upper.vertices = [6, 7, 8, 9, 12, 13, 14, 15];
        g.blocks.push(upper);
        let e = crate::block_triangle::split_triangle(&mut g, 0, 1, None, &mut Change::default())
            .unwrap_err();
        assert!(e.message.contains("one layer of blocks"), "{}", e.message);
    }
}
