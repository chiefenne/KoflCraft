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

use mesh_core::block::{BLOCK_FACES, Surface};
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
    /// Let the nodes of curved wall lines (and the block vertices on them)
    /// slide along the wall contour toward the foot of their first
    /// off-wall node; sharp corners and chain ends stay.
    pub slide_walls: bool,
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
            slide_walls: false,
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
    /// Wall nodes that slide along a contour, and the largest distance one
    /// moved from its TFI position.
    pub contour_nodes: usize,
    pub contour_slide_max: f64,
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

/// The largest angle in degrees between two wall lines that still continue
/// one contour; a sharper corner (a trailing edge) stays fixed.
const CONTOUR_CORNER_DEG: f64 = 30.0;
/// Samples of the contour polyline per grid cell.
const CONTOUR_SAMPLES: usize = 32;
/// Fraction of the way to the foot point a contour node moves per sweep.
const CONTOUR_RELAXATION: f64 = 0.5;
/// A contour node slides at most this fraction of its TFI cell toward
/// either neighbour; the Neumann foot alone drifts the whole distribution
/// toward a convex leading edge.
const CONTOUR_WINDOW: f64 = 0.5;
/// Weight of the TFI position against the foot point for a sliding node on
/// a projection surface. There the face nodes follow the rim and carry
/// their off-wall nodes along, so the foot point follows the node it is
/// for; without the anchor the rim drifts until it folds (elbow).
const SURFACE_ANCHOR: f64 = 1.0;

/// A chain of wall nodes along a curved wall line, across blocks; the end
/// nodes stay, the others slide on a dense polyline of the spline through
/// the TFI nodes.
struct Contour {
    nodes: Vec<usize>,
    /// The first node off the wall for each node (`usize::MAX`: none).
    off: Vec<usize>,
    pts: Vec<Vec3>,
    /// Arc length at each polyline point.
    arc: Vec<f64>,
    /// Arc length of each node.
    at: Vec<f64>,
    /// Arc length of each node in the TFI grid.
    at_tfi: Vec<f64>,
}

impl Contour {
    fn free(&self) -> usize {
        (1..self.nodes.len().saturating_sub(1))
            .filter(|&k| self.off[k] != usize::MAX)
            .count()
    }

    fn eval(&self, s: f64) -> Vec3 {
        let i = match self.arc.binary_search_by(|a| a.total_cmp(&s)) {
            Ok(i) => return self.pts[i],
            Err(i) => i.clamp(1, self.arc.len() - 1),
        };
        let t = (s - self.arc[i - 1]) / (self.arc[i] - self.arc[i - 1]);
        self.pts[i - 1] + (self.pts[i] - self.pts[i - 1]) * t
    }

    /// The arc length of the polyline point nearest to `p` within [lo, hi].
    fn project(&self, p: Vec3, lo: f64, hi: f64) -> f64 {
        let start = self.arc.partition_point(|&a| a < lo).saturating_sub(1);
        let end = self
            .arc
            .partition_point(|&a| a <= hi)
            .min(self.arc.len() - 1);
        let mut best = (f64::INFINITY, lo);
        for i in start..end {
            let d = self.pts[i + 1] - self.pts[i];
            let t = ((p - self.pts[i]).dot(d) / d.dot(d)).clamp(0.0, 1.0);
            let dist = (p - (self.pts[i] + d * t)).length();
            if dist < best.0 {
                best = (dist, self.arc[i] + t * (self.arc[i + 1] - self.arc[i]));
            }
        }
        best.1.clamp(lo, hi)
    }

    /// Move each free node part of the way to the foot of its off-wall
    /// node, keeping the order along the contour; the largest move over
    /// the node's shortest edge.
    fn slide(&mut self, x: &mut [Vec3], h_node: &[f64]) -> f64 {
        let mut change = 0.0f64;
        for k in 1..self.nodes.len() - 1 {
            if self.off[k] == usize::MAX {
                continue;
            }
            let (lo, hi) = (self.at[k - 1], self.at[k + 1]);
            let margin = 0.1 * (hi - lo);
            let t = &self.at_tfi;
            let lo = (lo + margin).max(t[k] - CONTOUR_WINDOW * (t[k] - t[k - 1]));
            let hi = (hi - margin).min(t[k] + CONTOUR_WINDOW * (t[k + 1] - t[k]));
            let foot = self.project(x[self.off[k]], lo, hi);
            let s = self.at[k] + CONTOUR_RELAXATION * (foot - self.at[k]);
            let s = s.clamp(lo, hi);
            let p = self.eval(s);
            let id = self.nodes[k];
            change = change.max((p - x[id]).length() / h_node[id]);
            x[id] = p;
            self.at[k] = s;
        }
        change
    }
}

/// Chains of wall nodes along the curved direction of wall faces whose
/// other direction is straight (extruded sections), joined across blocks
/// where the lines continue smoothly.
fn wall_contours(
    x: &[Vec3],
    hexes: &[[usize; 8]],
    boundary: &[([usize; 4], usize, usize)],
    blocks: &[Block],
    walls: &[Wall],
    wall_patches: &[usize],
    on_surface: &HashMap<usize, usize>,
) -> Result<Vec<Contour>, String> {
    let bend = |ids: &[usize]| -> f64 {
        let (a, b) = (x[ids[0]], x[ids[ids.len() - 1]]);
        let len = (b - a).length();
        ids[1..ids.len() - 1]
            .iter()
            .map(|&i| {
                let t = (x[i] - a).dot(b - a) / (len * len);
                (x[i] - (a + (b - a) * t)).length() / len
            })
            .fold(0.0, f64::max)
    };
    let mut lines: Vec<Vec<usize>> = Vec::new();
    for w in walls {
        let b = &blocks[w.block];
        let mid = b.at(w.node(w.n[w.t[0]] / 2, w.n[w.t[1]] / 2, 0));
        if on_surface.contains_key(&mid) {
            continue;
        }
        let line = |along: usize, other: usize| -> Vec<usize> {
            (0..=w.n[w.t[along]])
                .map(|i| {
                    let (u, v) = if along == 0 { (i, other) } else { (other, i) };
                    b.at(w.node(u, v, 0))
                })
                .collect()
        };
        let curved = |along: usize| -> f64 {
            let other = 1 - along;
            (0..=w.n[w.t[other]])
                .map(|j| bend(&line(along, j)))
                .fold(0.0, f64::max)
        };
        let (c0, c1) = (curved(0), curved(1));
        let along = if c0 >= c1 { 0 } else { 1 };
        let (big, small) = if along == 0 { (c0, c1) } else { (c1, c0) };
        // only ruled faces: one curved direction, the other straight
        if big < 1e-9 || small > 1e-9 {
            continue;
        }
        for j in 0..=w.n[w.t[1 - along]] {
            lines.push(line(along, j));
        }
    }
    // join lines that meet end to end and continue smoothly
    let dir = |a: usize, b: usize| {
        let d = x[b] - x[a];
        d / d.length()
    };
    let cos_limit = CONTOUR_CORNER_DEG.to_radians().cos();
    let mut ends: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, l) in lines.iter().enumerate() {
        ends.entry(l[0]).or_default().push(i);
        ends.entry(l[l.len() - 1]).or_default().push(i);
    }
    let mut used = vec![false; lines.len()];
    let mut chains: Vec<Vec<usize>> = Vec::new();
    for start in 0..lines.len() {
        if used[start] {
            continue;
        }
        used[start] = true;
        let mut chain = lines[start].clone();
        for _ in 0..2 {
            loop {
                let e = chain[chain.len() - 1];
                let prev = chain[chain.len() - 2];
                let Some(next) = ends.get(&e).and_then(|c| {
                    let open: Vec<usize> = c.iter().copied().filter(|&i| !used[i]).collect();
                    (open.len() == 1 && c.len() == 2).then(|| open[0])
                }) else {
                    break;
                };
                let mut l = lines[next].clone();
                if l[0] != e {
                    l.reverse();
                }
                if dir(prev, e).dot(dir(e, l[1])) < cos_limit || l[l.len() - 1] == chain[0] {
                    break;
                }
                used[next] = true;
                chain.extend_from_slice(&l[1..]);
            }
            chain.reverse();
        }
        chains.push(chain);
    }
    // the node one step off the wall: the neighbour not on a wall whose
    // edge is most normal to the contour
    let wall_set: std::collections::HashSet<usize> = boundary
        .iter()
        .filter(|(_, p, _)| wall_patches.contains(p))
        .flat_map(|(q, _, _)| *q)
        .collect();
    let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
    for h in hexes {
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
                if wall_set.contains(&p) && !wall_set.contains(&q) {
                    let e = adj.entry(p).or_default();
                    if !e.contains(&q) {
                        e.push(q);
                    }
                }
            }
        }
    }
    let mut out = Vec::new();
    for nodes in chains {
        if nodes.len() < 3 {
            continue;
        }
        let off: Vec<usize> = (0..nodes.len())
            .map(|k| {
                let id = nodes[k];
                let t = dir(
                    nodes[k.saturating_sub(1)],
                    nodes[(k + 1).min(nodes.len() - 1)],
                );
                adj.get(&id)
                    .and_then(|c| {
                        c.iter().copied().min_by(|&p, &q| {
                            dir(id, p).dot(t).abs().total_cmp(&dir(id, q).dot(t).abs())
                        })
                    })
                    .unwrap_or(usize::MAX)
            })
            .collect();
        let (a, b) = (x[nodes[0]], x[nodes[nodes.len() - 1]]);
        let through: Vec<Vec3> = nodes[1..nodes.len() - 1].iter().map(|&i| x[i]).collect();
        let curve = mesh_core::block::Curve::Spline { through };
        let ev = super::distribution::Evaluator::new(a, b, &curve)?;
        let m = CONTOUR_SAMPLES * (nodes.len() - 1);
        let pts: Vec<Vec3> = (0..=m)
            .map(|i| match i {
                0 => a,
                _ if i == m => b,
                _ => ev.at(i as f64 / m as f64),
            })
            .collect();
        let mut arc = vec![0.0];
        for w in pts.windows(2) {
            arc.push(arc[arc.len() - 1] + (w[1] - w[0]).length());
        }
        let mut c = Contour {
            nodes,
            off,
            pts,
            arc,
            at: Vec::new(),
            at_tfi: Vec::new(),
        };
        let total = c.arc[c.arc.len() - 1];
        let mut lo = 0.0;
        c.at = c
            .nodes
            .iter()
            .enumerate()
            .map(|(k, &id)| {
                let s = if k == 0 {
                    0.0
                } else if k == c.nodes.len() - 1 {
                    total
                } else {
                    c.project(x[id], lo, total)
                };
                lo = s;
                s
            })
            .collect();
        c.at_tfi = c.at.clone();
        out.push(c);
    }
    Ok(out)
}

/// The wall nodes whose boundary quads are all wall quads of one
/// projection surface.
fn surface_nodes(
    boundary: &[([usize; 4], usize, usize)],
    boundary_surfaces: &[Option<usize>],
    wall_patches: &[usize],
) -> HashMap<usize, usize> {
    let mut of: HashMap<usize, Option<usize>> = HashMap::new();
    for ((q, p, _), s) in boundary.iter().zip(boundary_surfaces) {
        let s = if wall_patches.contains(p) { *s } else { None };
        for &v in q {
            let e = of.entry(v).or_insert(s);
            if *e != s {
                *e = None;
            }
        }
    }
    of.into_iter()
        .filter_map(|(v, s)| s.map(|s| (v, s)))
        .collect()
}

/// A wall node on a block edge of a projection surface (cylinder, torus):
/// it slides on the surface toward the foot of its off-wall node, inside a
/// window of `CONTOUR_WINDOW` of each TFI wall edge at it.
struct RimSlide {
    surface: Surface,
    id: usize,
    off: usize,
    tfi: Vec3,
    /// TFI edges from the node to its wall neighbours.
    edges: Vec<Vec3>,
}

impl RimSlide {
    fn slide(&self, x: &mut [Vec3], h_node: &[f64]) -> f64 {
        let here = x[self.id];
        let foot = self.surface.project(x[self.off]);
        let target = (foot + self.tfi * SURFACE_ANCHOR) / (1.0 + SURFACE_ANCHOR);
        let d = here + (target - here) * CONTOUR_RELAXATION - self.tfi;
        let mut t = 1.0f64;
        for e in &self.edges {
            let (a, cap) = (d.dot(*e), CONTOUR_WINDOW * e.dot(*e));
            if a > cap {
                t = t.min(cap / a);
            }
        }
        let p = self.surface.project(self.tfi + d * t);
        let change = (p - here).length() / h_node[self.id];
        x[self.id] = p;
        change
    }
}

/// The wall nodes of one block face on a projection surface: the nodes
/// inside follow the displacement of its rim by transfinite interpolation
/// at their TFI arc-length fractions, and stay on the surface.
struct SurfaceFace {
    surface: Surface,
    nu: usize,
    nv: usize,
    ids: Vec<usize>,
    tfi: Vec<Vec3>,
    /// Arc-length fractions along the u and the v lines.
    xi: Vec<f64>,
    eta: Vec<f64>,
}

impl SurfaceFace {
    fn follow(&self, x: &mut [Vec3], h_node: &[f64]) -> f64 {
        let at = |u: usize, v: usize| u + (self.nu + 1) * v;
        let d = |x: &[Vec3], u: usize, v: usize| x[self.ids[at(u, v)]] - self.tfi[at(u, v)];
        let (nu, nv) = (self.nu, self.nv);
        let mut change = 0.0f64;
        for v in 1..nv {
            for u in 1..nu {
                let k = at(u, v);
                let (s, t) = (self.xi[k], self.eta[k]);
                let e = d(x, 0, v) * (1.0 - s)
                    + d(x, nu, v) * s
                    + d(x, u, 0) * (1.0 - t)
                    + d(x, u, nv) * t
                    - (d(x, 0, 0) * ((1.0 - s) * (1.0 - t))
                        + d(x, nu, 0) * (s * (1.0 - t))
                        + d(x, 0, nv) * ((1.0 - s) * t)
                        + d(x, nu, nv) * (s * t));
                let id = self.ids[k];
                let p = self.surface.project(self.tfi[k] + e);
                change = change.max((p - x[id]).length() / h_node[id]);
                x[id] = p;
            }
        }
        change
    }
}

/// The sliding of the walls on projection surfaces: rim nodes and faces.
fn surface_sliding(
    x: &[Vec3],
    hexes: &[[usize; 8]],
    boundary: &[([usize; 4], usize, usize)],
    surfaces: &[Surface],
    on_surface: &HashMap<usize, usize>,
    blocks: &[Block],
    walls: &[Wall],
) -> (Vec<RimSlide>, Vec<SurfaceFace>) {
    let mut faces = Vec::new();
    let mut rim: std::collections::BTreeSet<usize> = Default::default();
    for w in walls {
        let b = &blocks[w.block];
        let (nu, nv) = (w.n[w.t[0]], w.n[w.t[1]]);
        let mid = b.at(w.node(nu / 2, nv / 2, 0));
        let Some(&sn) = on_surface.get(&mid) else {
            continue;
        };
        let mut ids = Vec::with_capacity((nu + 1) * (nv + 1));
        for v in 0..=nv {
            for u in 0..=nu {
                let id = b.at(w.node(u, v, 0));
                ids.push(id);
                if (u == 0 || u == nu || v == 0 || v == nv) && on_surface.contains_key(&id) {
                    rim.insert(id);
                }
            }
        }
        let at = |u: usize, v: usize| u + (nu + 1) * v;
        let mut xi = vec![0.0; ids.len()];
        let mut eta = vec![0.0; ids.len()];
        for v in 0..=nv {
            let mut acc = vec![0.0];
            for u in 1..=nu {
                acc.push(acc[u - 1] + (x[ids[at(u, v)]] - x[ids[at(u - 1, v)]]).length());
            }
            for u in 0..=nu {
                xi[at(u, v)] = acc[u] / acc[nu];
            }
        }
        for u in 0..=nu {
            let mut acc = vec![0.0];
            for v in 1..=nv {
                acc.push(acc[v - 1] + (x[ids[at(u, v)]] - x[ids[at(u, v - 1)]]).length());
            }
            for v in 0..=nv {
                eta[at(u, v)] = acc[v] / acc[nv];
            }
        }
        let tfi = ids.iter().map(|&i| x[i]).collect();
        faces.push(SurfaceFace {
            surface: surfaces[sn],
            nu,
            nv,
            ids,
            tfi,
            xi,
            eta,
        });
    }
    let mut on_boundary: std::collections::HashSet<usize> = Default::default();
    let mut wall_nbrs: HashMap<usize, Vec<usize>> = HashMap::new();
    for (q, _, _) in boundary {
        on_boundary.extend(q.iter().copied());
        for k in 0..4 {
            let (a, b) = (q[k], q[(k + 1) % 4]);
            for (p, r) in [(a, b), (b, a)] {
                if rim.contains(&p) {
                    let e = wall_nbrs.entry(p).or_default();
                    if !e.contains(&r) {
                        e.push(r);
                    }
                }
            }
        }
    }
    let mut off_nbrs: HashMap<usize, Vec<usize>> = HashMap::new();
    for h in hexes {
        for (a, b) in HEX_EDGES {
            for (p, q) in [(h[a], h[b]), (h[b], h[a])] {
                if rim.contains(&p) && !on_boundary.contains(&q) {
                    let e = off_nbrs.entry(p).or_default();
                    if !e.contains(&q) {
                        e.push(q);
                    }
                }
            }
        }
    }
    let mut slides = Vec::new();
    for id in rim {
        let (Some(nb), Some(offs)) = (wall_nbrs.get(&id), off_nbrs.get(&id)) else {
            continue;
        };
        let edges: Vec<Vec3> = nb.iter().map(|&r| x[r] - x[id]).collect();
        // the off-wall neighbour most normal to the wall edges
        let score = |c: usize| {
            let d = x[c] - x[id];
            let d = d / d.length();
            edges
                .iter()
                .map(|e| (d.dot(*e) / e.length()).abs())
                .fold(0.0, f64::max)
        };
        let Some(off) = offs
            .iter()
            .copied()
            .min_by(|&p, &q| score(p).total_cmp(&score(q)))
        else {
            continue;
        };
        slides.push(RimSlide {
            surface: surfaces[on_surface[&id]],
            id,
            off,
            tfi: x[id],
            edges,
        });
    }
    (slides, faces)
}

const HEX_EDGES: [(usize, usize); 12] = [
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
];

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
    // walls on a projection surface (pipe, elbow): the block-edge nodes
    // slide on it, anchored to their TFI position, and the face nodes
    // follow; the other curved walls (extruded sections) slide node by
    // node along their contour lines
    let surface_nodes = if o.slide_walls {
        surface_nodes(&grid.boundary, &grid.boundary_surfaces, &o.wall_patches)
    } else {
        Default::default()
    };
    let (rim_slides, surface_faces) = if o.slide_walls {
        surface_sliding(
            x,
            &grid.hexes,
            &grid.boundary,
            &grid.surfaces,
            &surface_nodes,
            &blocks,
            &walls,
        )
    } else {
        Default::default()
    };
    let mut contours = if o.slide_walls {
        wall_contours(
            x,
            &grid.hexes,
            &grid.boundary,
            &blocks,
            &walls,
            &o.wall_patches,
            &surface_nodes,
        )?
    } else {
        Vec::new()
    };
    let tfi_points: Vec<Vec3> = if o.slide_walls {
        x.to_vec()
    } else {
        Vec::new()
    };
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
        contour_nodes: contours.iter().map(|c| c.free()).sum::<usize>()
            + rim_slides.len()
            + surface_faces
                .iter()
                .map(|f| (f.nu - 1) * (f.nv - 1))
                .sum::<usize>(),
        contour_slide_max: 0.0,
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
        if !contours.is_empty() || !rim_slides.is_empty() {
            for c in &mut contours {
                change = change.max(c.slide(x, &h_node));
            }
            for r in &rim_slides {
                change = change.max(r.slide(x, &h_node));
            }
            for f in &surface_faces {
                change = change.max(f.follow(x, &h_node));
            }
            // the wall moved: the target slopes follow its new normals
            for w in &mut walls {
                let b = &blocks[w.block];
                let sign = if w.side == 0 { 1.0 } else { -1.0 };
                for v in 0..=w.n[w.t[1]] {
                    for u in 0..=w.n[w.t[0]] {
                        if let Some(nrm) = inward_normal(x, b, w, u, v) {
                            let s = w.slot(u, v);
                            w.target[s] = nrm * (w.spacing[s] * sign);
                        }
                    }
                }
            }
        }
        report.last_change = change;
        if !report.last_change.is_finite() {
            return Err("elliptic smoothing diverged".into());
        }
        // converged when nodes stop moving and the first spacing is held
        report.converged = report.last_change < o.tolerance && spacing_error < SPACING_TOLERANCE;
    }
    for f in &surface_faces {
        for (k, &id) in f.ids.iter().enumerate() {
            report.contour_slide_max = report.contour_slide_max.max((x[id] - f.tfi[k]).length());
        }
    }
    for c in &contours {
        for &id in &c.nodes {
            report.contour_slide_max = report
                .contour_slide_max
                .max((x[id] - tfi_points[id]).length());
        }
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
            let push = w.spacing_feedback[s] - o.spacing_gain * sign * (h / w.spacing[s]).ln();
            // A node whose feedback sits at the control bound and still
            // pushes outward holds the closest spacing the bounded controls
            // reach (short blunt faces such as a trailing edge); it does not
            // keep the grid from converging.
            if push.abs() <= o.max_control {
                spacing_error = spacing_error.max((h / w.spacing[s]).ln().abs());
            }
            w.spacing_feedback[s] = push.clamp(-o.max_control, o.max_control);
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
