//! Elliptic smoothing of a multi-block grid with wall orthogonality
//! control (BS-09, MB1b; residual freeze, sliding and line solver MB1d).
//!
//! Every node that is not on the domain boundary moves. A node inside a
//! block, or inside a face shared by two blocks, solves the
//! Thompson–Thames–Mastin system in the index space (ξ¹, ξ², ξ³) of its
//! home block:
//!
//! Σᵢⱼ (aⁱ·aʲ) r_ξⁱξʲ + Σₖ (aᵏ·aᵏ) Pₖ r_ξᵏ = 0,  a¹ = r_ξ² × r_ξ³, …
//!
//! which is J² times the inverse Laplace (Winslow) equation with control
//! functions Pₖ (J. F. Thompson, F. C. Thames, C. W. Mastin, J. Comput.
//! Phys. 15, 1974; Thompson, Warsi and Mastin, *Numerical Grid
//! Generation*, 1985, chapter 6). For a node on a shared face the stencil
//! reaches into the neighbouring block: the node one step off the face is
//! the neighbour's node one step off the same face node, so the two blocks
//! form one structured neighbourhood whatever their local frames are.
//! Nodes on block edges inside the domain (where three or more blocks
//! meet, as at the corners of an O-grid core) have no structured stencil;
//! they move to the average of their edge-connected neighbours, restricted
//! to the plane normal to the block edge so the distribution along the
//! edge is kept. Block vertices inside the domain (where the segments of
//! a swept topology meet) move normal to the grid line through them.
//! Boundary nodes never move (BS-08), except on a planar symmetry patch:
//! there a node slides in the plane, its stencil reaches across the plane
//! to mirror images of the interior nodes, and it is projected back onto
//! the plane after every update. The maintainer chose smoothing across
//! block faces (2026-11-24) and sliding on symmetry planes on cards.
//!
//! **Residual freeze.** The TFI grid does not solve the system with only
//! the distribution controls: at the joins of a swept bend the Laplace
//! part straightens the lines and drags the whole core toward the bend
//! centre. So the residual of every node at the TFI grid is stored and
//! subtracted at every update, which makes the TFI grid a fixed point;
//! only the wall controls move nodes. Thomas–Middlecoff controls were
//! measured as the alternative and were not better on both the pipe and
//! the elbow (concept doc, MB1d).
//!
//! **Solver.** Line Gauss–Seidel: the nodes of each grid line inside a
//! block are solved together by the Thomas algorithm, the line axis turns
//! every sweep. Convergence: the largest move of a sweep below the
//! tolerance times the shortest TFI edge at the node, and the first
//! spacing off every wall within `SPACING_TOLERANCE`.
//!
//! Control functions:
//! - **Distribution.** From the TFI grid, Pₖ = −(r_ξᵏ · r_ξᵏξᵏ)/|r_ξᵏ|²
//!   at every node, frozen. A grid of straight, orthogonal, stretched
//!   lines satisfies the system exactly, so the edge laws survive.
//! - **Walls** (Steger and Sorenson, "Automatic mesh-point clustering near
//!   a boundary in grid generation with elliptic partial differential
//!   equations", J. Comput. Phys. 33, 1979; Sorenson, NASA TM-81198, 1980;
//!   the 3D form of 3DGRAPE, NASA TM-102224, 1989). On a wall face the
//!   target off-wall derivative is the inward normal times the first
//!   spacing of the TFI grid; the system evaluated there with that
//!   derivative gives the three controls that make it hold, which are
//!   relaxed into place iteration by iteration and decay into the block as
//!   exp(−a (d − 1)) with the index distance d from the wall.

#![allow(clippy::needless_range_loop)] // loops over the three axes read clearer by index

use std::collections::HashMap;

use mesh_core::block::BLOCK_FACES;
use mesh_core::geometry::Vec3;

use super::grid::BlockGrid;

/// Node index (i, j, k) inside a block.
type Index = [usize; 3];

/// The largest |ln(h / h_target)| of the first spacing off a wall at
/// convergence.
const SPACING_TOLERANCE: f64 = 1e-4;

/// Options of the elliptic smoother.
#[derive(Clone, Debug, PartialEq)]
pub struct SmoothOptions {
    /// Patches whose faces get orthogonality control.
    pub wall_patches: Vec<usize>,
    /// Patches whose nodes slide in their plane (symmetry); a patch that is
    /// not planar keeps its nodes fixed.
    pub sliding_patches: Vec<usize>,
    pub max_iterations: usize,
    /// Stop when the largest node move of a sweep is below this fraction of
    /// the shortest TFI edge at that node.
    pub tolerance: f64,
    /// Decay rate a of the wall controls per index step.
    pub decay: f64,
    /// Under-relaxation of the wall controls per iteration.
    pub control_relaxation: f64,
    /// Bound on |Pₖ| of a wall control.
    pub max_control: f64,
    /// Gain of the first-spacing feedback per iteration.
    pub spacing_gain: f64,
}

impl Default for SmoothOptions {
    fn default() -> Self {
        Self {
            wall_patches: Vec::new(),
            sliding_patches: Vec::new(),
            max_iterations: 10000,
            tolerance: 1e-4,
            decay: 0.45,
            control_relaxation: 0.05,
            max_control: 4.0,
            spacing_gain: 0.1,
        }
    }
}

/// What smoothing did, over all blocks.
#[derive(Clone, Debug, PartialEq)]
pub struct Smoothing {
    pub iterations: usize,
    pub converged: bool,
    /// Largest move of the last sweep over the shortest TFI edge at the
    /// node.
    pub last_change: f64,
    /// Nodes that moved: inside blocks, on shared faces, on inner block
    /// edges.
    pub block_nodes: usize,
    pub face_nodes: usize,
    pub edge_nodes: usize,
    /// Block vertices inside the domain, moved to their neighbours' mean.
    pub vertex_nodes: usize,
    /// Planar sliding patches, and the nodes on them that moved.
    pub sliding_patches: usize,
    pub sliding_nodes: usize,
    pub wall_faces: usize,
    /// Deviation of the first grid line off a wall from the wall normal,
    /// in degrees, over the interior nodes of the wall faces.
    pub wall_angle_max_before: f64,
    pub wall_angle_max_after: f64,
    pub wall_angle_mean_before: f64,
    pub wall_angle_mean_after: f64,
    /// Largest |h_after / h_before − 1| of the first spacing off a wall.
    pub first_spacing_change: f64,
}

struct Block {
    n: [usize; 3],
    ids: Vec<usize>,
    /// For each face: the other block and its face, if shared.
    across: [Option<(usize, usize)>; 6],
    /// For each face: the plane it slides on, if it is on a sliding patch.
    mirror: [Option<usize>; 6],
}

/// A run of structured nodes of one block along one axis, in order.
struct Line {
    block: usize,
    nodes: Vec<(usize, [usize; 3])>,
}

/// A plane of a sliding patch.
#[derive(Clone, Copy)]
struct Plane {
    point: Vec3,
    normal: Vec3,
}

impl Plane {
    fn project(&self, p: Vec3) -> Vec3 {
        p - self.normal * (p - self.point).dot(self.normal)
    }
    fn reflect(&self, p: Vec3) -> Vec3 {
        p - self.normal * (2.0 * (p - self.point).dot(self.normal))
    }
    fn reflect_direction(&self, d: Vec3) -> Vec3 {
        d - self.normal * (2.0 * d.dot(self.normal))
    }
}

impl Block {
    fn at(&self, c: [usize; 3]) -> usize {
        self.ids[c[0] + (self.n[0] + 1) * (c[1] + (self.n[1] + 1) * c[2])]
    }
}

/// One wall face of a block: the off-wall axis and side, and per face node
/// the target off-wall derivative and the control.
struct Wall {
    block: usize,
    axis: usize,
    side: usize,
    /// Tangential axes.
    t: [usize; 2],
    n: [usize; 3],
    target: Vec<Vec3>,
    control: Vec<[f64; 3]>,
    spacing: Vec<f64>,
    /// Integral feedback on the off-wall control that pins the first
    /// spacing (the system alone holds any spacing with a matching ratio).
    spacing_feedback: Vec<f64>,
    /// The off-wall distribution control of the TFI grid one node in.
    base_off: Vec<f64>,
    /// For each tangential axis and side: the plane of a sliding face there,
    /// whose mirror image continues the wall face (so its rim nodes are
    /// controlled like the interior).
    mirror: [[Option<Plane>; 2]; 2],
}

impl Wall {
    /// The controlled face nodes along tangential axis `t`: the interior,
    /// plus a rim on a sliding face.
    fn range(&self, t: usize) -> std::ops::RangeInclusive<usize> {
        let n = self.n[self.t[t]];
        let lo = if self.mirror[t][0].is_some() { 0 } else { 1 };
        let hi = if self.mirror[t][1].is_some() {
            n
        } else {
            n - 1
        };
        lo..=hi
    }

    /// Index (u, v) of the face, mirrored back into range, and the planes
    /// to reflect by.
    fn resolve(&self, u: isize, v: isize) -> (usize, usize, [Option<Plane>; 2]) {
        let mut planes = [None; 2];
        let mut fold = |w: isize, t: usize| -> usize {
            let n = self.n[self.t[t]] as isize;
            if w < 0 {
                planes[t] = self.mirror[t][0];
                (-w) as usize
            } else if w > n {
                planes[t] = self.mirror[t][1];
                (2 * n - w) as usize
            } else {
                w as usize
            }
        };
        let (uu, vv) = (fold(u, 0), fold(v, 1));
        (uu, vv, planes)
    }

    /// The node at (u, v) and off-wall depth `depth`, mirrored where (u, v)
    /// lies beyond a sliding face.
    fn point(&self, x: &[Vec3], b: &Block, u: isize, v: isize, depth: usize) -> Vec3 {
        let (uu, vv, planes) = self.resolve(u, v);
        planes
            .iter()
            .flatten()
            .fold(x[b.at(self.node(uu, vv, depth))], |p, pl| pl.reflect(p))
    }

    /// The target off-wall derivative at (u, v), mirrored likewise.
    fn target_at(&self, u: isize, v: isize) -> Vec3 {
        let (uu, vv, planes) = self.resolve(u, v);
        planes
            .iter()
            .flatten()
            .fold(self.target[self.slot(uu, vv)], |d, pl| {
                pl.reflect_direction(d)
            })
    }

    fn slot(&self, u: usize, v: usize) -> usize {
        u + (self.n[self.t[0]] + 1) * v
    }
    fn node(&self, u: usize, v: usize, depth: usize) -> [usize; 3] {
        let mut c = [0; 3];
        c[self.t[0]] = u;
        c[self.t[1]] = v;
        c[self.axis] = if self.side == 0 {
            depth
        } else {
            self.n[self.axis] - depth
        };
        c
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Role {
    Fixed,
    /// A structured node: its home block and index there.
    Structured {
        block: usize,
        c: [usize; 3],
    },
    /// On an inner block edge, with its two neighbours along the edge.
    Edge {
        prev: usize,
        next: usize,
    },
    /// A block vertex inside the domain.
    Vertex,
}

struct Stencil<'a> {
    blocks: &'a [Block],
    /// (block, face) → (face node → node one step into that block).
    off: &'a HashMap<(usize, usize), HashMap<usize, usize>>,
    planes: &'a [Plane],
}

impl Stencil<'_> {
    /// The point at offset `d` from index `c` of block `b`, crossing at
    /// most one shared face and one sliding face (whose ghost is the mirror
    /// image of the node inside).
    fn point(&self, x: &[Vec3], b: usize, c: [usize; 3], d: [isize; 3]) -> Option<Vec3> {
        let blk = &self.blocks[b];
        let mut q = [0usize; 3];
        let mut out = None;
        let mut reflect = None;
        for a in 0..3 {
            let v = c[a] as isize + d[a];
            if v < 0 || v > blk.n[a] as isize {
                let f = 2 * a + usize::from(v > 0);
                if let (Some(pl), None) = (blk.mirror[f], reflect) {
                    reflect = Some(pl);
                    q[a] = (c[a] as isize - d[a]) as usize;
                    continue;
                }
                if out.is_some() {
                    return None;
                }
                out = Some(f);
                q[a] = c[a];
            } else {
                q[a] = v as usize;
            }
        }
        let id = match out {
            None => blk.at(q),
            Some(f) => {
                let (other, g) = blk.across[f]?;
                *self.off.get(&(other, g))?.get(&blk.at(q))?
            }
        };
        Some(match reflect {
            Some(pl) => self.planes[pl].reflect(x[id]),
            None => x[id],
        })
    }

    fn derivatives(
        &self,
        x: &[Vec3],
        b: usize,
        c: [usize; 3],
    ) -> Option<([Vec3; 3], [[Vec3; 3]; 3])> {
        let p = |d: [isize; 3]| self.point(x, b, c, d);
        let unit = |a: usize, s: isize| {
            let mut d = [0isize; 3];
            d[a] = s;
            d
        };
        let centre = x[self.blocks[b].at(c)];
        let mut first = [Vec3::ZERO; 3];
        let mut second = [[Vec3::ZERO; 3]; 3];
        for a in 0..3 {
            let (plus, minus) = (p(unit(a, 1))?, p(unit(a, -1))?);
            first[a] = (plus - minus) * 0.5;
            second[a][a] = plus + minus - centre * 2.0;
            for e in a + 1..3 {
                let mut pp = [0isize; 3];
                pp[a] = 1;
                pp[e] = 1;
                let mut pm = pp;
                pm[e] = -1;
                let mut mp = pp;
                mp[a] = -1;
                let mut mm = pm;
                mm[a] = -1;
                let v = (p(pp)? - p(pm)? - p(mp)? + p(mm)?) * 0.25;
                second[a][e] = v;
                second[e][a] = v;
            }
        }
        Some((first, second))
    }
}

/// aⁱ·aʲ of the cofactor vectors.
fn metric(r: &[Vec3; 3]) -> [[f64; 3]; 3] {
    let a = [r[1].cross(r[2]), r[2].cross(r[0]), r[0].cross(r[1])];
    std::array::from_fn(|i| std::array::from_fn(|j| a[i].dot(a[j])))
}

fn solve3(m: [[f64; 3]; 3], rhs: [f64; 3]) -> Option<[f64; 3]> {
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(m);
    if !d.is_normal() {
        return None;
    }
    Some(std::array::from_fn(|c| {
        let mut mc = m;
        for row in mc.iter_mut() {
            row[c] = 0.0;
        }
        for r in 0..3 {
            mc[r][c] = rhs[r];
        }
        det(mc) / d
    }))
}

/// The unit wall normal pointing into the block at a wall face node.
fn inward_normal(x: &[Vec3], b: &Block, w: &Wall, u: usize, v: usize) -> Option<Vec3> {
    let face = |u: usize, v: usize| x[b.at(w.node(u, v, 0))];
    let (u0, u1) = (u.saturating_sub(1), (u + 1).min(w.n[w.t[0]]));
    let (v0, v1) = (v.saturating_sub(1), (v + 1).min(w.n[w.t[1]]));
    let normal = (face(u1, v) - face(u0, v)).cross(face(u, v1) - face(u, v0));
    let len = normal.length();
    if len.is_nan() || len <= 0.0 {
        return None;
    }
    let normal = normal / len;
    let off = x[b.at(w.node(u, v, 1))] - face(u, v);
    Some(if normal.dot(off) < 0.0 {
        normal * -1.0
    } else {
        normal
    })
}

/// Deviation in degrees of the first off-wall grid line from the normal,
/// over the interior nodes of the wall faces: (max, mean).
fn wall_angles(x: &[Vec3], blocks: &[Block], walls: &[Wall]) -> (f64, f64) {
    let (mut max, mut sum, mut count) = (0.0f64, 0.0, 0usize);
    for w in walls {
        let b = &blocks[w.block];
        for v in w.range(1) {
            for u in w.range(0) {
                let Some(nrm) = inward_normal(x, b, w, u, v) else {
                    continue;
                };
                let off = x[b.at(w.node(u, v, 1))] - x[b.at(w.node(u, v, 0))];
                let deg = (nrm.dot(off) / off.length())
                    .clamp(-1.0, 1.0)
                    .acos()
                    .to_degrees();
                max = max.max(deg);
                sum += deg;
                count += 1;
            }
        }
    }
    (max, if count > 0 { sum / count as f64 } else { 0.0 })
}

fn face_corner_ids(b: &Block, f: usize) -> [usize; 4] {
    let (axis, side, t0, t1) = BLOCK_FACES[f];
    let mut k = [(0, 0), (1, 0), (1, 1), (0, 1)].map(|(u, v)| {
        let mut c = [0; 3];
        c[axis] = side * b.n[axis];
        c[t0] = u * b.n[t0];
        c[t1] = v * b.n[t1];
        b.at(c)
    });
    k.sort_unstable();
    k
}

/// Smooth every node that is not on the boundary, in place.
pub fn smooth(grid: &mut BlockGrid, options: &SmoothOptions) -> Result<Smoothing, String> {
    let o = options;
    let mut blocks: Vec<Block> = (0..grid.blocks.len())
        .map(|b| Block {
            n: grid.blocks[b].cells,
            ids: grid.block_nodes(b),
            across: [None; 6],
            mirror: [None; 6],
        })
        .collect();
    // shared faces by their corner nodes
    let mut by_corners: HashMap<[usize; 4], Vec<(usize, usize)>> = HashMap::new();
    for (b, blk) in blocks.iter().enumerate() {
        for f in 0..6 {
            if grid.face_patches[b][f].is_none() {
                by_corners
                    .entry(face_corner_ids(blk, f))
                    .or_default()
                    .push((b, f));
            }
        }
    }
    for pair in by_corners.values() {
        if let [(a, f), (b, g)] = pair[..] {
            blocks[a].across[f] = Some((b, g));
            blocks[b].across[g] = Some((a, f));
        } else {
            return Err("a block face is shared by other than two blocks".into());
        }
    }
    let mut off: HashMap<(usize, usize), HashMap<usize, usize>> = HashMap::new();
    for (b, blk) in blocks.iter().enumerate() {
        for f in 0..6 {
            if blk.across[f].is_none() {
                continue;
            }
            let (axis, side, t0, t1) = BLOCK_FACES[f];
            let mut map = HashMap::new();
            for v in 0..=blk.n[t1] {
                for u in 0..=blk.n[t0] {
                    let mut c = [0; 3];
                    c[t0] = u;
                    c[t1] = v;
                    c[axis] = side * blk.n[axis];
                    let mut d = c;
                    d[axis] = if side == 0 { 1 } else { blk.n[axis] - 1 };
                    map.insert(blk.at(c), blk.at(d));
                }
            }
            off.insert((b, f), map);
        }
    }
    // roles
    // sliding planes: a sliding patch whose nodes lie on one plane
    let mut planes: Vec<Plane> = Vec::new();
    let mut plane_of_patch: HashMap<usize, usize> = HashMap::new();
    for &p in &o.sliding_patches {
        let quads: Vec<&[usize; 4]> = grid
            .boundary
            .iter()
            .filter(|(_, q, _)| *q == p)
            .map(|(f, _, _)| f)
            .collect();
        if quads.is_empty() {
            continue;
        }
        let x = &grid.points;
        let mut normal = Vec3::ZERO;
        let mut centre = Vec3::ZERO;
        let mut count = 0.0;
        let (mut lo, mut hi) = (x[quads[0][0]], x[quads[0][0]]);
        for q in &quads {
            normal = normal + (x[q[2]] - x[q[0]]).cross(x[q[3]] - x[q[1]]);
            for &v in q.iter() {
                centre = centre + x[v];
                count += 1.0;
                lo = Vec3::new(lo.x.min(x[v].x), lo.y.min(x[v].y), lo.z.min(x[v].z));
                hi = Vec3::new(hi.x.max(x[v].x), hi.y.max(x[v].y), hi.z.max(x[v].z));
            }
        }
        let length = normal.length();
        if length.is_nan() || length <= 0.0 {
            continue;
        }
        let plane = Plane {
            point: centre / count,
            normal: normal / length,
        };
        let tol = 1e-9 * (hi - lo).length();
        let flat = quads.iter().all(|q| {
            q.iter()
                .all(|&v| (x[v] - plane.project(x[v])).length() <= tol)
        });
        if flat {
            plane_of_patch.insert(p, planes.len());
            planes.push(plane);
        }
    }
    for (b, blk) in blocks.iter_mut().enumerate() {
        for f in 0..6 {
            blk.mirror[f] = grid.face_patches[b][f].and_then(|p| plane_of_patch.get(&p).copied());
        }
    }
    // the planes each boundary node slides on; a node on any other patch
    // stays fixed
    let mut slides_on: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut role = vec![None::<Role>; grid.points.len()];
    for (q, p, _) in &grid.boundary {
        for &v in q {
            match plane_of_patch.get(p) {
                Some(&pl) => {
                    let e = slides_on.entry(v).or_default();
                    if !e.contains(&pl) {
                        e.push(pl);
                    }
                }
                None => role[v] = Some(Role::Fixed),
            }
        }
    }
    slides_on.retain(|v, _| role[*v].is_none());
    for (b, blk) in blocks.iter().enumerate() {
        let n = blk.n;
        for k in 0..=n[2] {
            for j in 0..=n[1] {
                for i in 0..=n[0] {
                    let c = [i, j, k];
                    let id = blk.at(c);
                    if role[id].is_some() {
                        continue;
                    }
                    // limits on sliding faces do not count: the mirror
                    // image continues the block there
                    let at_limit: Vec<usize> = (0..3)
                        .filter(|&a| {
                            (c[a] == 0 && blk.mirror[2 * a].is_none())
                                || (c[a] == n[a] && blk.mirror[2 * a + 1].is_none())
                        })
                        .collect();
                    role[id] = Some(match at_limit.len() {
                        0 | 1 => Role::Structured { block: b, c },
                        2 if {
                            let a = (0..3).find(|a| !at_limit.contains(a)).unwrap_or(0);
                            c[a] > 0 && c[a] < n[a]
                        } =>
                        {
                            let a = (0..3).find(|a| !at_limit.contains(a)).unwrap_or(0);
                            let (mut p, mut q) = (c, c);
                            p[a] -= 1;
                            q[a] += 1;
                            Role::Edge {
                                prev: blk.at(p),
                                next: blk.at(q),
                            }
                        }
                        _ => Role::Vertex,
                    });
                }
            }
        }
    }
    let stencil = Stencil {
        blocks: &blocks,
        off: &off,
        planes: &planes,
    };
    let x = &mut grid.points;
    // graph neighbours of the edge nodes
    let mut neighbours: HashMap<usize, Vec<usize>> = HashMap::new();
    for h in &grid.hexes {
        for (a, b) in [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ] {
            for (p, q) in [(h[a], h[b]), (h[b], h[a])] {
                if matches!(role[p], Some(Role::Edge { .. } | Role::Vertex)) {
                    let e = neighbours.entry(p).or_default();
                    if !e.contains(&q) {
                        e.push(q);
                    }
                }
            }
        }
    }
    // the shortest TFI edge at each node: moves are measured against it,
    // so a strongly stretched grid converges near the wall and far away alike
    let mut h_node = vec![f64::INFINITY; x.len()];
    for h in &grid.hexes {
        for (a, b) in [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ] {
            let l = (x[h[b]] - x[h[a]]).length();
            h_node[h[a]] = h_node[h[a]].min(l);
            h_node[h[b]] = h_node[h[b]].min(l);
        }
    }
    // distribution controls, frozen from the TFI grid
    let structured: Vec<(usize, usize, [usize; 3])> = role
        .iter()
        .enumerate()
        .filter_map(|(id, r)| match r {
            Some(Role::Structured { block, c }) => Some((id, *block, *c)),
            _ => None,
        })
        .collect();
    let edges: Vec<(usize, usize, usize)> = role
        .iter()
        .enumerate()
        .filter_map(|(id, r)| match r {
            Some(Role::Edge { prev, next }) => Some((id, *prev, *next)),
            _ => None,
        })
        .collect();
    let vertices: Vec<usize> = role
        .iter()
        .enumerate()
        .filter_map(|(id, r)| matches!(r, Some(Role::Vertex)).then_some(id))
        .collect();
    let mut base: HashMap<usize, [f64; 3]> = HashMap::new();
    for &(id, b, c) in &structured {
        let (r, rr) = stencil
            .derivatives(x, b, c)
            .ok_or_else(|| format!("node {id} of block {b} has no structured stencil"))?;
        base.insert(
            id,
            std::array::from_fn(|a| -r[a].dot(rr[a][a]) / r[a].dot(r[a])),
        );
    }
    let mut walls: Vec<Wall> = Vec::new();
    for (b, blk) in blocks.iter().enumerate() {
        for (f, &(axis, side, t0, t1)) in BLOCK_FACES.iter().enumerate() {
            if grid.face_patches[b][f].is_some_and(|p| o.wall_patches.contains(&p))
                && blk.n[axis] >= 2
            {
                let size = (blk.n[t0] + 1) * (blk.n[t1] + 1);
                let mut w = Wall {
                    block: b,
                    axis,
                    side,
                    t: [t0, t1],
                    n: blk.n,
                    target: vec![Vec3::ZERO; size],
                    control: vec![[0.0; 3]; size],
                    spacing: vec![0.0; size],
                    spacing_feedback: vec![0.0; size],
                    base_off: vec![0.0; size],
                    mirror: [[t0, 0], [t0, 1], [t1, 0], [t1, 1]]
                        .map(|[t, side]| blk.mirror[2 * t + side].map(|pl| planes[pl]))
                        .chunks(2)
                        .map(|c| [c[0], c[1]])
                        .collect::<Vec<_>>()
                        .try_into()
                        .unwrap_or([[None; 2]; 2]),
                };
                for v in 0..=blk.n[t1] {
                    for u in 0..=blk.n[t0] {
                        let nrm = inward_normal(x, blk, &w, u, v).ok_or_else(|| {
                            format!("block {b}: a wall face has a degenerate node")
                        })?;
                        let h = (x[blk.at(w.node(u, v, 1))] - x[blk.at(w.node(u, v, 0))]).length();
                        let sign = if side == 0 { 1.0 } else { -1.0 };
                        let s = w.slot(u, v);
                        w.target[s] = nrm * (h * sign);
                        w.spacing[s] = h;
                        if let Some(p) = base.get(&blk.at(w.node(u, v, 1))) {
                            w.control[s] = *p;
                            w.base_off[s] = p[axis];
                        }
                    }
                }
                walls.push(w);
            }
        }
    }
    let walls_of: Vec<Vec<usize>> = (0..blocks.len())
        .map(|b| (0..walls.len()).filter(|&w| walls[w].block == b).collect())
        .collect();
    // offsets of the edge and vertex nodes from their neighbours' mean in
    // the TFI grid, frozen like the distribution controls, so the TFI grid
    // is a fixed point there too and stretching survives
    let edge_offset: Vec<Vec3> = edges
        .iter()
        .map(|&(id, prev, next)| x[id] - off_edge_mean(x, &neighbours[&id], prev, next))
        .collect();
    // a vertex through which a block edge runs straight on (two neighbours
    // in opposite directions, as where the segments of a swept O-grid
    // meet) slides only normal to that line, like an edge node
    let vertex_line: Vec<(usize, usize)> = vertices
        .iter()
        .map(|&id| {
            let nb = &neighbours[&id];
            let dir = |q: usize| {
                let d = x[q] - x[id];
                d / d.length()
            };
            let mut best = (usize::MAX, usize::MAX, -0.9);
            for (a, &p) in nb.iter().enumerate() {
                for &q in &nb[a + 1..] {
                    let c = dir(p).dot(dir(q));
                    if c < best.2 {
                        best = (p, q, c);
                    }
                }
            }
            (best.0, best.1)
        })
        .collect();
    let vertex_offset: Vec<Vec3> = vertices
        .iter()
        .zip(&vertex_line)
        .map(|(&id, &(p, q))| x[id] - off_edge_mean(x, &neighbours[&id], p, q))
        .collect();
    // what the system with the distribution controls alone would still
    // move on the TFI grid (the curvature of grid lines across their own
    // direction, as round a bend), frozen so the TFI grid is a fixed point
    let mut residual: HashMap<usize, Vec3> = HashMap::new();
    {
        for &(id, b, c) in &structured {
            let Some((r, rr)) = stencil.derivatives(x, b, c) else {
                continue;
            };
            let p = base[&id];
            let g = metric(&r);
            let centre = x[id];
            let mut num = Vec3::ZERO;
            let mut den = 0.0;
            for a in 0..3 {
                num = num + (rr[a][a] + centre * 2.0) * g[a][a] + r[a] * (g[a][a] * p[a]);
                den += 2.0 * g[a][a];
                for e in a + 1..3 {
                    num = num + rr[a][e] * (2.0 * g[a][e]);
                }
            }
            residual.insert(id, centre - num / den);
        }
    }
    let (angle_max_before, angle_mean_before) = wall_angles(x, &blocks, &walls);

    // a sliding node lands back on its plane(s) after every update
    let land = |id: usize, p: Vec3| match slides_on.get(&id) {
        Some(pls) => pls.iter().fold(p, |q, &pl| planes[pl].project(q)),
        None => p,
    };
    let mut report = Smoothing {
        iterations: 0,
        converged: structured.is_empty() && edges.is_empty() && vertices.is_empty(),
        last_change: 0.0,
        block_nodes: structured
            .iter()
            .filter(|(_, b, c)| (0..3).all(|a| c[a] > 0 && c[a] < blocks[*b].n[a]))
            .count(),
        face_nodes: 0,
        edge_nodes: edges.len(),
        vertex_nodes: vertices.len(),
        sliding_patches: planes.len(),
        sliding_nodes: slides_on.len(),
        wall_faces: walls.len(),
        wall_angle_max_before: angle_max_before,
        wall_angle_max_after: angle_max_before,
        wall_angle_mean_before: angle_mean_before,
        wall_angle_mean_after: angle_mean_before,
        first_spacing_change: 0.0,
    };
    report.face_nodes = structured.len() - report.block_nodes;
    // the control functions at a structured node: frozen distribution
    // controls plus the wall controls decaying into the block
    let controls = |walls: &[Wall], id: usize, b: usize, c: [usize; 3]| -> [f64; 3] {
        let mut p = base[&id];
        for &wi in &walls_of[b] {
            let w = &walls[wi];
            let depth = if w.side == 0 {
                c[w.axis]
            } else {
                w.n[w.axis] - c[w.axis]
            };
            if depth == 0 {
                continue;
            }
            let weight = (-o.decay * (depth as f64 - 1.0)).exp();
            let clamp = |v: usize, r: std::ops::RangeInclusive<usize>| {
                v.clamp(*r.start(), (*r.end()).max(*r.start()))
            };
            let (u, v) = (clamp(c[w.t[0]], w.range(0)), clamp(c[w.t[1]], w.range(1)));
            let s = w.slot(u, v);
            let above = base
                .get(&blocks[b].at(w.node(u, v, 1)))
                .copied()
                .unwrap_or([0.0; 3]);
            for a in 0..3 {
                p[a] += weight * (w.control[s][a] - above[a]);
            }
        }
        p
    };
    // the grid lines of structured nodes along each axis, by home block
    let lines: [Vec<Line>; 3] = std::array::from_fn(|axis| {
        let mut by: HashMap<(usize, Index), Vec<(usize, Index)>> = HashMap::new();
        for &(id, b, c) in &structured {
            let mut key = c;
            key[axis] = 0;
            by.entry((b, key)).or_default().push((id, c));
        }
        let mut out = Vec::new();
        let mut keys: Vec<_> = by.keys().copied().collect();
        keys.sort_unstable();
        for key in keys {
            let mut nodes = by.remove(&key).unwrap_or_default();
            nodes.sort_unstable_by_key(|n| n.1[axis]);
            let mut run: Vec<(usize, [usize; 3])> = Vec::new();
            for n in nodes {
                if run.last().is_some_and(|l| l.1[axis] + 1 != n.1[axis]) {
                    out.push(Line {
                        block: key.0,
                        nodes: std::mem::take(&mut run),
                    });
                }
                run.push(n);
            }
            if !run.is_empty() {
                out.push(Line {
                    block: key.0,
                    nodes: run,
                });
            }
        }
        out
    });
    while !report.converged && report.iterations < o.max_iterations {
        report.iterations += 1;
        let mut spacing_error = 0.0f64;
        for w in &mut walls {
            spacing_error = spacing_error.max(update_wall(x, &blocks[w.block], w, o));
        }
        let mut change = 0.0f64;
        // line Gauss–Seidel: the nodes of each grid line along `axis` are
        // solved together (tridiagonal in the line direction, the rest of
        // the stencil from the current grid); the axis turns every sweep,
        // so information crosses strongly stretched cells in one sweep
        let axis = report.iterations % 3;
        for line in &lines[axis] {
            let b = line.block;
            let m = line.nodes.len();
            let (mut lo, mut di, mut hi, mut rhs) = (
                vec![0.0; m],
                vec![0.0; m],
                vec![0.0; m],
                vec![Vec3::ZERO; m],
            );
            let mut ok = true;
            for (k, &(id, c)) in line.nodes.iter().enumerate() {
                let p = controls(&walls, id, b, c);
                let Some((r, rr)) = stencil.derivatives(x, b, c) else {
                    ok = false;
                    break;
                };
                let g = metric(&r);
                let centre = x[id];
                let mut num = Vec3::ZERO;
                let mut den = 0.0;
                for a in 0..3 {
                    num = num + (rr[a][a] + centre * 2.0) * g[a][a] + r[a] * (g[a][a] * p[a]);
                    den += 2.0 * g[a][a];
                    for e in a + 1..3 {
                        num = num + rr[a][e] * (2.0 * g[a][e]);
                    }
                }
                if den.is_nan() || den <= 0.0 {
                    return Err(format!("the metric degenerates at node {id} of block {b}"));
                }
                // the line terms g (1 ∓ P/2) x∓ move to the left
                let (cm, cp) = (
                    g[axis][axis] * (1.0 - 0.5 * p[axis]),
                    g[axis][axis] * (1.0 + 0.5 * p[axis]),
                );
                let mut d = [0isize; 3];
                d[axis] = 1;
                let plus = stencil.point(x, b, c, d);
                d[axis] = -1;
                let minus = stencil.point(x, b, c, d);
                let (Some(plus), Some(minus)) = (plus, minus) else {
                    ok = false;
                    break;
                };
                let mut f = num - plus * cp - minus * cm
                    + residual.get(&id).copied().unwrap_or(Vec3::ZERO) * den;
                di[k] = den;
                if k > 0 {
                    lo[k] = -cm;
                } else {
                    f = f + minus * cm;
                }
                if k + 1 < m {
                    hi[k] = -cp;
                } else {
                    f = f + plus * cp;
                }
                rhs[k] = f;
            }
            if !ok {
                continue;
            }
            // Thomas algorithm
            for k in 1..m {
                let w = lo[k] / di[k - 1];
                di[k] -= w * hi[k - 1];
                rhs[k] = rhs[k] - rhs[k - 1] * w;
            }
            let mut sol = vec![Vec3::ZERO; m];
            sol[m - 1] = rhs[m - 1] / di[m - 1];
            for k in (0..m - 1).rev() {
                sol[k] = (rhs[k] - sol[k + 1] * hi[k]) / di[k];
            }
            for (k, &(id, _)) in line.nodes.iter().enumerate() {
                let centre = x[id];
                let new = land(id, sol[k]);
                change = change.max((new - centre).length() / h_node[id]);
                x[id] = new;
            }
        }
        for (k, &(id, prev, next)) in edges.iter().enumerate() {
            let target = off_edge_mean(x, &neighbours[&id], prev, next) + edge_offset[k];
            let t = x[next] - x[prev];
            let t = t / t.length();
            let d = target - x[id];
            let new = land(id, x[id] + (d - t * d.dot(t)));
            change = change.max((new - x[id]).length() / h_node[id]);
            x[id] = new;
        }
        for (k, &id) in vertices.iter().enumerate() {
            let (prev, next) = vertex_line[k];
            let target = off_edge_mean(x, &neighbours[&id], prev, next) + vertex_offset[k];
            let new = if prev == usize::MAX {
                target
            } else {
                let t = x[next] - x[prev];
                let t = t / t.length();
                let d = target - x[id];
                x[id] + (d - t * d.dot(t))
            };
            let new = land(id, new);
            change = change.max((new - x[id]).length() / h_node[id]);
            x[id] = new;
        }
        report.last_change = change;
        if !report.last_change.is_finite() {
            return Err("elliptic smoothing diverged".into());
        }
        // converged when nodes stop moving and the first spacing is held
        report.converged = report.last_change < o.tolerance && spacing_error < SPACING_TOLERANCE;
    }
    let (max_after, mean_after) = wall_angles(x, &blocks, &walls);
    report.wall_angle_max_after = max_after;
    report.wall_angle_mean_after = mean_after;
    for w in &walls {
        let b = &blocks[w.block];
        for v in w.range(1) {
            for u in w.range(0) {
                let h = (x[b.at(w.node(u, v, 1))] - x[b.at(w.node(u, v, 0))]).length();
                report.first_spacing_change = report
                    .first_spacing_change
                    .max((h / w.spacing[w.slot(u, v)] - 1.0).abs());
            }
        }
    }
    Ok(report)
}

/// The mean of the neighbours other than `prev` and `next` (the
/// neighbours along an edge would pull a curved edge towards its centre of
/// curvature).
fn off_edge_mean(x: &[Vec3], neighbours: &[usize], prev: usize, next: usize) -> Vec3 {
    let mut sum = Vec3::ZERO;
    let mut count = 0;
    for &q in neighbours {
        if q != prev && q != next {
            sum = sum + x[q];
            count += 1;
        }
    }
    sum / count as f64
}

/// The wall controls: the system evaluated at the wall with the target
/// off-wall derivative, solved for (P₁, P₂, P₃) and relaxed toward it.
/// Returns the largest |ln(h / h_target)| of the first spacing.
fn update_wall(x: &[Vec3], b: &Block, w: &mut Wall, o: &SmoothOptions) -> f64 {
    let mut spacing_error = 0.0f64;
    let (a, ta, tb) = (w.axis, w.t[0], w.t[1]);
    let sign = if w.side == 0 { 1.0 } else { -1.0 };
    for v in w.range(1) {
        for u in w.range(0) {
            let s = w.slot(u, v);
            let node = |d: usize| x[b.at(w.node(u, v, d))];
            let (ui, vi) = (u as isize, v as isize);
            let tan = |du: isize, dv: isize| w.point(x, b, ui + du, vi + dv, 0);
            let target = |du: isize, dv: isize| w.target_at(ui + du, vi + dv);
            let rn = w.target[s];
            let mut r = [Vec3::ZERO; 3];
            let mut rr = [[Vec3::ZERO; 3]; 3];
            r[ta] = (tan(1, 0) - tan(-1, 0)) * 0.5;
            r[tb] = (tan(0, 1) - tan(0, -1)) * 0.5;
            r[a] = rn;
            rr[ta][ta] = tan(1, 0) + tan(-1, 0) - node(0) * 2.0;
            rr[tb][tb] = tan(0, 1) + tan(0, -1) - node(0) * 2.0;
            let cross = (tan(1, 1) - tan(1, -1) - tan(-1, 1) + tan(-1, -1)) * 0.25;
            rr[ta][tb] = cross;
            rr[tb][ta] = cross;
            // off-wall second derivative of a cubic with the target slope
            rr[a][a] = (node(1) * 8.0 - node(0) * 7.0 - node(2)) * 0.5 - rn * (3.0 * sign);
            let ru = (target(1, 0) - target(-1, 0)) * 0.5;
            let rv = (target(0, 1) - target(0, -1)) * 0.5;
            rr[a][ta] = ru;
            rr[ta][a] = ru;
            rr[a][tb] = rv;
            rr[tb][a] = rv;
            let g = metric(&r);
            let mut rhs = Vec3::ZERO;
            for i in 0..3 {
                for j in 0..3 {
                    rhs = rhs - rr[i][j] * g[i][j];
                }
            }
            let cols = [r[0] * g[0][0], r[1] * g[1][1], r[2] * g[2][2]];
            let m = [
                [cols[0].x, cols[1].x, cols[2].x],
                [cols[0].y, cols[1].y, cols[2].y],
                [cols[0].z, cols[1].z, cols[2].z],
            ];
            // first spacing too large: cluster harder (P toward the wall)
            let h = (node(1) - node(0)).length();
            spacing_error = spacing_error.max((h / w.spacing[s]).ln().abs());
            w.spacing_feedback[s] = (w.spacing_feedback[s]
                - o.spacing_gain * sign * (h / w.spacing[s]).ln())
            .clamp(-o.max_control, o.max_control);
            if let Some(mut p) = solve3(m, [rhs.x, rhs.y, rhs.z]) {
                // The system with the target slope gives the tangential
                // controls (orthogonality). Its off-wall control depends on
                // how the one-sided r_ζζ is estimated and holds any spacing
                // with a matching ratio, so the off-wall control is the
                // TFI distribution control plus the spacing feedback.
                p[a] = w.base_off[s] + w.spacing_feedback[s];
                for k in 0..3 {
                    let want = p[k].clamp(-o.max_control, o.max_control);
                    w.control[s][k] += o.control_relaxation * (want - w.control[s][k]);
                }
            }
        }
    }
    spacing_error
}
