//! Node placement and hexahedra of a checked block topology.
//!
//! Every topological entity is gridded once: vertices, then each edge (its
//! curve and law), then each block face by 2D transfinite interpolation of
//! its four edges, then each block interior by 3D transfinite interpolation
//! of its six faces (Gordon and Hall, "Construction of curvilinear
//! co-ordinate systems and applications to mesh generation", Int. J. Numer.
//! Meth. Eng. 7, 1973; Thompson, Warsi and Mastin, *Numerical Grid
//! Generation*, 1985, chapter 8). A block that meets a face or edge built
//! by another block reads the stored nodes through the orientation map of
//! the two local frames, so shared faces match node for node without any
//! tolerance welding (BS-05), and boundary nodes are never moved (BS-08).
//!
//! The interior nodes of a face listed under a surface are projected onto
//! it (BS-07) before the block interiors are filled, so the interiors see
//! the exact wall.
//!
//! Interior parameters follow the edge distributions: where opposite edges
//! carry different laws, the parameter of a node blends the chord-length
//! fractions of its edges, so grid lines join matching nodes.

use std::collections::HashMap;

use mesh_core::block::{BLOCK_FACES, BlockTopology, corner_slot, face_key};
use mesh_core::geometry::Vec3;

use super::distribution::{Evaluator, fractions};
use super::{
    Checked, DEFAULT_EDGE, Issue, check, corner_edges, edge_key, edge_slots, face_vertices,
};

/// The cells of one block: cell `first_cell + i + ni (j + nj k)` is the
/// cell (i, j, k).
#[derive(Clone, Debug, PartialEq)]
pub struct BlockRange {
    pub block: usize,
    pub domain: String,
    pub cells: [usize; 3],
    pub first_cell: usize,
}

/// Achieved spacing of one edge (BS-06).
#[derive(Clone, Debug, PartialEq)]
pub struct EdgeReport {
    pub vertices: [usize; 2],
    pub cells: usize,
    pub length: f64,
    pub first_spacing: f64,
    pub last_spacing: f64,
    /// Largest ratio of neighbouring cell lengths along the edge, ≥ 1.
    pub max_growth: f64,
}

/// Per block quality (BS-13).
#[derive(Clone, Debug, PartialEq)]
pub struct BlockQuality {
    pub block: usize,
    /// Smallest corner Jacobian determinant of any hexahedron.
    pub min_jacobian: f64,
    /// Smallest corner determinant of unit edge vectors, in (−1, 1].
    pub min_scaled_jacobian: f64,
    /// Largest ratio of longest to shortest edge of a hexahedron.
    pub max_aspect_ratio: f64,
}

#[derive(Clone, Debug)]
pub struct BlockGrid {
    pub points: Vec<Vec3>,
    /// Hexahedra in VTK order: (i,j,k), (i+1,j,k), (i+1,j+1,k), (i,j+1,k),
    /// then the same at k+1.
    pub hexes: Vec<[usize; 8]>,
    pub blocks: Vec<BlockRange>,
    /// Boundary quads wound outward, with patch index and owner cell.
    pub boundary: Vec<([usize; 4], usize, usize)>,
    pub edges: Vec<EdgeReport>,
    pub checked: Checked,
    /// The patch of each of a block's six faces (order of `BLOCK_FACES`),
    /// `None` for a face shared with another block.
    pub face_patches: Vec<[Option<usize>; 6]>,
}

/// Outward faces of a hexahedron in local node slots.
pub const HEX_FACES: [[usize; 4]; 6] = [
    [0, 4, 7, 3], // i = 0
    [1, 2, 6, 5], // i = 1
    [0, 1, 5, 4], // j = 0
    [3, 7, 6, 2], // j = 1
    [0, 3, 2, 1], // k = 0
    [4, 5, 6, 7], // k = 1
];

struct StoredFace {
    corners: [usize; 4],
    dims: [usize; 2],
    ids: Vec<usize>,
}

struct Builder<'a> {
    t: &'a BlockTopology,
    points: Vec<Vec3>,
    edges: HashMap<(usize, usize), Vec<usize>>,
    faces: HashMap<[usize; 4], StoredFace>,
    reports: Vec<EdgeReport>,
}

impl Builder<'_> {
    /// Node ids from `a` to `b` along the edge, built on first use.
    fn edge(&mut self, a: usize, b: usize, n: usize) -> Result<Vec<usize>, String> {
        let key = edge_key(a, b);
        if !self.edges.contains_key(&key) {
            let spec = self
                .t
                .edges
                .iter()
                .find(|e| edge_key(e.vertices[0], e.vertices[1]) == key);
            let (from, to, curve, law) = match spec {
                Some(s) => (
                    s.vertices[0],
                    s.vertices[1],
                    s.curve.clone(),
                    s.distribution,
                ),
                None => (key.0, key.1, DEFAULT_EDGE.0, DEFAULT_EDGE.1),
            };
            let ev = Evaluator::new(self.t.vertices[from], self.t.vertices[to], &curve)?;
            let s = fractions(law, n, ev.length())?;
            let mut ids = Vec::with_capacity(n + 1);
            ids.push(from);
            for &f in &s[1..n] {
                ids.push(self.points.len());
                self.points.push(ev.at(f));
            }
            ids.push(to);
            let lengths: Vec<f64> = ids
                .windows(2)
                .map(|w| (self.points[w[1]] - self.points[w[0]]).length())
                .collect();
            let max_growth = lengths
                .windows(2)
                .map(|w| (w[1] / w[0]).max(w[0] / w[1]))
                .fold(1.0, f64::max);
            self.reports.push(EdgeReport {
                vertices: [from, to],
                cells: n,
                length: ev.length(),
                first_spacing: lengths[0],
                last_spacing: lengths[n - 1],
                max_growth,
            });
            if from != key.0 {
                ids.reverse();
            }
            self.edges.insert(key, ids);
        }
        let ids = &self.edges[&key];
        if ids.len() != n + 1 {
            return Err(format!("edge {}–{} has two cell counts", key.0, key.1));
        }
        let mut ids = ids.clone();
        if a != key.0 {
            ids.reverse();
        }
        Ok(ids)
    }
}

/// Chord-length fractions of a polyline of points.
fn chord_fractions(p: &[Vec3]) -> Vec<f64> {
    let mut s = vec![0.0; p.len()];
    for k in 1..p.len() {
        s[k] = s[k - 1] + (p[k] - p[k - 1]).length();
    }
    let total = s[p.len() - 1];
    s.iter_mut().for_each(|v| *v /= total);
    s
}

/// Build the grid of a topology: check, then grid edges, faces, blocks.
pub fn build(t: &BlockTopology) -> Result<BlockGrid, Vec<Issue>> {
    let checked = check(t)?;
    let fail = |m: String| {
        vec![Issue {
            kind: "grid",
            message: m,
        }]
    };
    let mut g = Builder {
        t,
        points: t.vertices.clone(),
        edges: HashMap::new(),
        faces: HashMap::new(),
        reports: Vec::new(),
    };
    let mut hexes = Vec::with_capacity(checked.cells);
    let mut blocks = Vec::with_capacity(t.blocks.len());
    let mut boundary = Vec::new();
    let mut face_patches = Vec::with_capacity(t.blocks.len());
    let mut patch_of: HashMap<[usize; 4], usize> = HashMap::new();
    for (p, patch) in t.patches.iter().enumerate() {
        for &f in &patch.faces {
            patch_of.insert(face_key(f), p);
        }
    }
    let mut surface_of: HashMap<[usize; 4], usize> = HashMap::new();
    for (sn, spec) in t.surfaces.iter().enumerate() {
        for &f in &spec.faces {
            surface_of.insert(face_key(f), sn);
        }
    }
    for (b, block) in t.blocks.iter().enumerate() {
        let n = block.cells;
        let idx = |i: usize, j: usize, k: usize| i + (n[0] + 1) * (j + (n[1] + 1) * k);
        let mut local = vec![usize::MAX; (n[0] + 1) * (n[1] + 1) * (n[2] + 1)];
        let at = |c: [usize; 3]| idx(c[0], c[1], c[2]);
        // edges
        for (axis, s, e) in edge_slots() {
            let ids = g
                .edge(block.vertices[s], block.vertices[e], n[axis])
                .map_err(fail)?;
            let sc = super::slot_corner(s);
            let c: [usize; 3] = std::array::from_fn(|d| sc[d] * n[d]);
            for (m, id) in ids.into_iter().enumerate() {
                let mut q = c;
                q[axis] = m;
                local[at(q)] = id;
            }
        }
        // faces
        for (f, &(fixed, side, a, bx)) in BLOCK_FACES.iter().enumerate() {
            let corners = face_vertices(&block.vertices, f);
            let key = face_key(corners);
            let (m1, m2) = (n[a], n[bx]);
            let pos = |ca: usize, cb: usize| {
                let mut q = [0; 3];
                q[fixed] = side * n[fixed];
                q[a] = ca;
                q[bx] = cb;
                at(q)
            };
            if let Some(sf) = g.faces.get(&key) {
                let place = |v: usize| sf.corners.iter().position(|&c| c == v);
                let qc = |p: usize| -> [i64; 2] {
                    let [n1, n2] = sf.dims.map(|d| d as i64);
                    [[0, 0], [n1, 0], [n1, n2], [0, n2]][p]
                };
                let (Some(p00), Some(p10), Some(p01)) =
                    (place(corners[0]), place(corners[1]), place(corners[3]))
                else {
                    return Err(fail(format!(
                        "face {key:?} does not match its stored corners"
                    )));
                };
                let (o, ea, eb) = (qc(p00), qc(p10), qc(p01));
                let da = [ea[0] - o[0], ea[1] - o[1]];
                let db = [eb[0] - o[0], eb[1] - o[1]];
                let len = |d: [i64; 2]| d[0].abs() + d[1].abs();
                if len(da) != m1 as i64 || len(db) != m2 as i64 {
                    return Err(fail(format!("face {key:?} has two sizes")));
                }
                let (ua, ub) = (
                    [da[0] / m1 as i64, da[1] / m1 as i64],
                    [db[0] / m2 as i64, db[1] / m2 as i64],
                );
                let w = sf.dims[0] as i64 + 1;
                for cb in 0..=m2 {
                    for ca in 0..=m1 {
                        let p = o[0] + ua[0] * ca as i64 + ub[0] * cb as i64;
                        let q = o[1] + ua[1] * ca as i64 + ub[1] * cb as i64;
                        let id = sf.ids[(p + w * q) as usize];
                        let slot = pos(ca, cb);
                        if local[slot] != usize::MAX && local[slot] != id {
                            return Err(fail(format!(
                                "face {key:?} does not match the nodes of its edges"
                            )));
                        }
                        local[slot] = id;
                    }
                }
            } else {
                // 2D TFI of the face from its four edges
                let pt = |local: &Vec<usize>, ca, cb| g.points[local[pos(ca, cb)]];
                let bottom: Vec<Vec3> = (0..=m1).map(|i| pt(&local, i, 0)).collect();
                let top: Vec<Vec3> = (0..=m1).map(|i| pt(&local, i, m2)).collect();
                let left: Vec<Vec3> = (0..=m2).map(|j| pt(&local, 0, j)).collect();
                let right: Vec<Vec3> = (0..=m2).map(|j| pt(&local, m1, j)).collect();
                let (sb, st, sl, sr) = (
                    chord_fractions(&bottom),
                    chord_fractions(&top),
                    chord_fractions(&left),
                    chord_fractions(&right),
                );
                let (p00, p10, p11, p01) = (bottom[0], bottom[m1], top[m1], top[0]);
                for j in 1..m2 {
                    for i in 1..m1 {
                        // u = sb + v (st − sb), v = sl + u (sr − sl)
                        let (db, dl) = (st[i] - sb[i], sr[j] - sl[j]);
                        let u = (sb[i] + sl[j] * db) / (1.0 - dl * db);
                        let v = sl[j] + u * dl;
                        let p =
                            bottom[i] * (1.0 - v) + top[i] * v + left[j] * (1.0 - u) + right[j] * u
                                - (p00 * ((1.0 - u) * (1.0 - v))
                                    + p10 * (u * (1.0 - v))
                                    + p11 * (u * v)
                                    + p01 * ((1.0 - u) * v));
                        local[pos(i, j)] = g.points.len();
                        g.points.push(match surface_of.get(&key) {
                            Some(&sn) => t.surfaces[sn].surface.project(p),
                            None => p,
                        });
                    }
                }
                let mut ids = Vec::with_capacity((m1 + 1) * (m2 + 1));
                for j in 0..=m2 {
                    for i in 0..=m1 {
                        ids.push(local[pos(i, j)]);
                    }
                }
                g.faces.insert(
                    key,
                    StoredFace {
                        corners,
                        dims: [m1, m2],
                        ids,
                    },
                );
            }
        }
        // 3D TFI of the interior
        fill_interior(&mut g.points, &mut local, n);
        if local.contains(&usize::MAX) {
            return Err(fail(format!("block {b} has nodes that were not placed")));
        }
        let first_cell = hexes.len();
        for k in 0..n[2] {
            for j in 0..n[1] {
                for i in 0..n[0] {
                    hexes.push([
                        local[idx(i, j, k)],
                        local[idx(i + 1, j, k)],
                        local[idx(i + 1, j + 1, k)],
                        local[idx(i, j + 1, k)],
                        local[idx(i, j, k + 1)],
                        local[idx(i + 1, j, k + 1)],
                        local[idx(i + 1, j + 1, k + 1)],
                        local[idx(i, j + 1, k + 1)],
                    ]);
                }
            }
        }
        let mut patches_here = [None; 6];
        for (f, &(fixed, side, _, _)) in BLOCK_FACES.iter().enumerate() {
            let Some(&patch) = patch_of.get(&face_key(face_vertices(&block.vertices, f))) else {
                continue;
            };
            patches_here[f] = Some(patch);
            for k in 0..n[2] {
                for j in 0..n[1] {
                    for i in 0..n[0] {
                        let c = [i, j, k];
                        if c[fixed] != side * (n[fixed] - 1) {
                            continue;
                        }
                        let cell = first_cell + i + n[0] * (j + n[1] * k);
                        let hex = hexes[cell];
                        boundary.push((HEX_FACES[f].map(|s| hex[s]), patch, cell));
                    }
                }
            }
        }
        face_patches.push(patches_here);
        blocks.push(BlockRange {
            block: b,
            domain: block.domain.clone(),
            cells: n,
            first_cell,
        });
    }
    if g.points.len() != checked.points {
        return Err(fail(format!(
            "{} nodes placed where the topology has {}",
            g.points.len(),
            checked.points
        )));
    }
    Ok(BlockGrid {
        points: g.points,
        hexes,
        blocks,
        boundary,
        edges: g.reports,
        checked,
        face_patches,
    })
}

fn fill_interior(points: &mut Vec<Vec3>, local: &mut [usize], n: [usize; 3]) {
    let idx = |i: usize, j: usize, k: usize| i + (n[0] + 1) * (j + (n[1] + 1) * k);
    let x = |points: &Vec<Vec3>, local: &[usize], i, j, k| points[local[idx(i, j, k)]];
    // fractions of the four edges along each axis, by the edge's corner
    let edge_s = |axis: usize, c1: usize, c2: usize| -> Vec<f64> {
        let p: Vec<Vec3> = (0..=n[axis])
            .map(|m| {
                let mut q = [0; 3];
                q[axis] = m;
                q[(axis + 1) % 3] = c1 * n[(axis + 1) % 3];
                q[(axis + 2) % 3] = c2 * n[(axis + 2) % 3];
                points[local[idx(q[0], q[1], q[2])]]
            })
            .collect();
        chord_fractions(&p)
    };
    let s: [[Vec<f64>; 4]; 3] = std::array::from_fn(|axis| {
        [(0, 0), (1, 0), (0, 1), (1, 1)].map(|(c1, c2)| edge_s(axis, c1, c2))
    });
    let blend = |e: &[Vec<f64>; 4], m: usize, a: f64, b: f64| {
        (1.0 - a) * (1.0 - b) * e[0][m]
            + a * (1.0 - b) * e[1][m]
            + (1.0 - a) * b * e[2][m]
            + a * b * e[3][m]
    };
    let [ni, nj, nk] = n;
    for k in 1..nk {
        for j in 1..nj {
            for i in 1..ni {
                // parameters: u along i blends the i-edges over (v, w), etc.
                let (mut u, mut v, mut w) = (
                    i as f64 / ni as f64,
                    j as f64 / nj as f64,
                    k as f64 / nk as f64,
                );
                for _ in 0..8 {
                    u = blend(&s[0], i, v, w);
                    v = blend(&s[1], j, w, u);
                    w = blend(&s[2], k, u, v);
                }
                let xp = |a, b, c| x(points, local, a, b, c);
                let faces = xp(0, j, k) * (1.0 - u)
                    + xp(ni, j, k) * u
                    + xp(i, 0, k) * (1.0 - v)
                    + xp(i, nj, k) * v
                    + xp(i, j, 0) * (1.0 - w)
                    + xp(i, j, nk) * w;
                let edges = xp(0, 0, k) * ((1.0 - u) * (1.0 - v))
                    + xp(0, nj, k) * ((1.0 - u) * v)
                    + xp(ni, 0, k) * (u * (1.0 - v))
                    + xp(ni, nj, k) * (u * v)
                    + xp(0, j, 0) * ((1.0 - u) * (1.0 - w))
                    + xp(0, j, nk) * ((1.0 - u) * w)
                    + xp(ni, j, 0) * (u * (1.0 - w))
                    + xp(ni, j, nk) * (u * w)
                    + xp(i, 0, 0) * ((1.0 - v) * (1.0 - w))
                    + xp(i, 0, nk) * ((1.0 - v) * w)
                    + xp(i, nj, 0) * (v * (1.0 - w))
                    + xp(i, nj, nk) * (v * w);
                let mut corners = Vec3::ZERO;
                for (ci, cj, ck) in [
                    (0, 0, 0),
                    (1, 0, 0),
                    (0, 1, 0),
                    (1, 1, 0),
                    (0, 0, 1),
                    (1, 0, 1),
                    (0, 1, 1),
                    (1, 1, 1),
                ] {
                    let wgt = (if ci == 1 { u } else { 1.0 - u })
                        * (if cj == 1 { v } else { 1.0 - v })
                        * (if ck == 1 { w } else { 1.0 - w });
                    corners = corners + xp(ci * ni, cj * nj, ck * nk) * wgt;
                }
                let p = faces - edges + corners;
                local[idx(i, j, k)] = points.len();
                points.push(p);
            }
        }
    }
}

impl BlockGrid {
    /// Corner Jacobians and aspect ratios per block (BS-13).
    pub fn quality(&self) -> Vec<BlockQuality> {
        self.blocks
            .iter()
            .map(|r| {
                let count = r.cells.iter().product::<usize>();
                let mut q = BlockQuality {
                    block: r.block,
                    min_jacobian: f64::INFINITY,
                    min_scaled_jacobian: f64::INFINITY,
                    max_aspect_ratio: 1.0,
                };
                for hex in &self.hexes[r.first_cell..r.first_cell + count] {
                    let p = hex.map(|v| self.points[v]);
                    let (mut lo, mut hi) = (f64::INFINITY, 0.0f64);
                    for slot in 0..8 {
                        let e = corner_edges(&p, slot);
                        let det = e[0].cross(e[1]).dot(e[2]);
                        let l = e.map(Vec3::length);
                        q.min_jacobian = q.min_jacobian.min(det);
                        q.min_scaled_jacobian =
                            q.min_scaled_jacobian.min(det / (l[0] * l[1] * l[2]));
                        for v in l {
                            lo = lo.min(v);
                            hi = hi.max(v);
                        }
                    }
                    q.max_aspect_ratio = q.max_aspect_ratio.max(hi / lo);
                }
                q
            })
            .collect()
    }

    /// The node ids of block `b` in its local grid, index
    /// i + (ni + 1) (j + (nj + 1) k).
    pub fn block_nodes(&self, b: usize) -> Vec<usize> {
        let [ni, nj, nk] = self.blocks[b].cells;
        let mut ids = Vec::with_capacity((ni + 1) * (nj + 1) * (nk + 1));
        for k in 0..=nk {
            for j in 0..=nj {
                for i in 0..=ni {
                    ids.push(self.block_corner(b, i, j, k));
                }
            }
        }
        ids
    }

    /// The node ids of the corner (i, j, k) of block `b` in its local grid.
    pub fn block_corner(&self, b: usize, i: usize, j: usize, k: usize) -> usize {
        let r = &self.blocks[b];
        let [ni, nj, _] = r.cells;
        let (ci, cj, ck) = (i.min(ni - 1), j.min(nj - 1), k.min(r.cells[2] - 1));
        let hex = self.hexes[r.first_cell + ci + ni * (cj + nj * ck)];
        hex[corner_slot(i - ci, j - cj, k - ck)]
    }
}

/// Six times the signed volume enclosed by a quad, by the fan of four
/// triangles from its vertex average, relative to `origin`.
fn quad_volume6(p: &[Vec3], q: [usize; 4], origin: Vec3) -> f64 {
    let c = (p[q[0]] + p[q[1]] + p[q[2]] + p[q[3]]) * 0.25 - origin;
    (0..4)
        .map(|k| {
            let (a, b) = (p[q[k]] - origin, p[q[(k + 1) % 4]] - origin);
            c.dot(a.cross(b))
        })
        .sum()
}

impl BlockGrid {
    /// The volume of every hexahedron, from its six faces split as the
    /// FoehnCraft hand-off splits them (fan from the face average).
    pub fn volumes(&self) -> Vec<f64> {
        self.hexes
            .iter()
            .map(|h| {
                let origin = self.points[h[0]];
                HEX_FACES
                    .iter()
                    .map(|f| quad_volume6(&self.points, f.map(|s| h[s]), origin))
                    .sum::<f64>()
                    / 6.0
            })
            .collect()
    }

    /// The volume enclosed by the boundary quads alone. Equal to the sum of
    /// the cell volumes up to rounding exactly when the interior faces
    /// cancel, which makes it an independent conservation check.
    pub fn boundary_volume(&self) -> f64 {
        let origin = self.points[0];
        self.boundary
            .iter()
            .map(|&(q, _, _)| quad_volume6(&self.points, q, origin))
            .sum::<f64>()
            / 6.0
    }
}
