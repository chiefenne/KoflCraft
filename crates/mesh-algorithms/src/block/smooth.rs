//! Elliptic smoothing of a multi-block grid with wall orthogonality
//! control (BS-09, MB1b).
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
//! edge is kept. Boundary nodes never move (BS-08). The maintainer chose
//! smoothing across block faces on a card (2026-11-24).
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

/// Options of the elliptic smoother.
#[derive(Clone, Debug, PartialEq)]
pub struct SmoothOptions {
    /// Patches whose faces get orthogonality control.
    pub wall_patches: Vec<usize>,
    pub max_iterations: usize,
    /// Stop when the largest node move of a sweep is below this fraction of
    /// the shortest TFI edge.
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
            max_iterations: 3000,
            tolerance: 1e-7,
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
    /// Largest move of the last sweep over the shortest TFI edge.
    pub last_change: f64,
    /// Nodes that moved: inside blocks, on shared faces, on inner block
    /// edges.
    pub block_nodes: usize,
    pub face_nodes: usize,
    pub edge_nodes: usize,
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
}

impl Wall {
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
}

struct Stencil<'a> {
    blocks: &'a [Block],
    /// (block, face) → (face node → node one step into that block).
    off: &'a HashMap<(usize, usize), HashMap<usize, usize>>,
}

impl Stencil<'_> {
    /// The node at offset `d` from index `c` of block `b`, crossing at most
    /// one shared face.
    fn node(&self, b: usize, c: [usize; 3], d: [isize; 3]) -> Option<usize> {
        let blk = &self.blocks[b];
        let mut q = [0usize; 3];
        let mut out = None;
        for a in 0..3 {
            let v = c[a] as isize + d[a];
            if v < 0 || v > blk.n[a] as isize {
                if out.is_some() {
                    return None;
                }
                out = Some((a, usize::from(v > 0)));
                q[a] = c[a];
            } else {
                q[a] = v as usize;
            }
        }
        match out {
            None => Some(blk.at(q)),
            Some((a, side)) => {
                let f = 2 * a + side;
                let (other, g) = blk.across[f]?;
                self.off.get(&(other, g))?.get(&blk.at(q)).copied()
            }
        }
    }

    fn derivatives(
        &self,
        x: &[Vec3],
        b: usize,
        c: [usize; 3],
    ) -> Option<([Vec3; 3], [[Vec3; 3]; 3])> {
        let p = |d: [isize; 3]| self.node(b, c, d).map(|v| x[v]);
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
        for v in 1..w.n[w.t[1]] {
            for u in 1..w.n[w.t[0]] {
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
    let mut role = vec![None::<Role>; grid.points.len()];
    for (q, _, _) in &grid.boundary {
        for &v in q {
            role[v] = Some(Role::Fixed);
        }
    }
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
                    let at_limit: Vec<usize> =
                        (0..3).filter(|&a| c[a] == 0 || c[a] == n[a]).collect();
                    role[id] = Some(match at_limit.len() {
                        0 | 1 => Role::Structured { block: b, c },
                        2 => {
                            let a = (0..3).find(|a| !at_limit.contains(a)).unwrap_or(0);
                            let (mut p, mut q) = (c, c);
                            p[a] -= 1;
                            q[a] += 1;
                            Role::Edge {
                                prev: blk.at(p),
                                next: blk.at(q),
                            }
                        }
                        _ => Role::Fixed,
                    });
                }
            }
        }
    }
    let stencil = Stencil {
        blocks: &blocks,
        off: &off,
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
                if matches!(role[p], Some(Role::Edge { .. })) {
                    let e = neighbours.entry(p).or_default();
                    if !e.contains(&q) {
                        e.push(q);
                    }
                }
            }
        }
    }
    let mut h_min = f64::INFINITY;
    for h in &grid.hexes {
        for (a, b) in [(0, 1), (1, 2), (0, 4)] {
            h_min = h_min.min((x[h[b]] - x[h[a]]).length());
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
    let (angle_max_before, angle_mean_before) = wall_angles(x, &blocks, &walls);

    let mut report = Smoothing {
        iterations: 0,
        converged: structured.is_empty() && edges.is_empty(),
        last_change: 0.0,
        block_nodes: structured
            .iter()
            .filter(|(_, b, c)| (0..3).all(|a| c[a] > 0 && c[a] < blocks[*b].n[a]))
            .count(),
        face_nodes: 0,
        edge_nodes: edges.len(),
        wall_faces: walls.len(),
        wall_angle_max_before: angle_max_before,
        wall_angle_max_after: angle_max_before,
        wall_angle_mean_before: angle_mean_before,
        wall_angle_mean_after: angle_mean_before,
        first_spacing_change: 0.0,
    };
    report.face_nodes = structured.len() - report.block_nodes;
    while !report.converged && report.iterations < o.max_iterations {
        report.iterations += 1;
        for w in &mut walls {
            update_wall(x, &blocks[w.block], w, o);
        }
        let mut change = 0.0f64;
        for &(id, b, c) in &structured {
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
                let clamp = |v: usize, m: usize| v.clamp(1, m.saturating_sub(1).max(1));
                let (u, v) = (clamp(c[w.t[0]], w.n[w.t[0]]), clamp(c[w.t[1]], w.n[w.t[1]]));
                let s = w.slot(u, v);
                let above = base
                    .get(&blocks[b].at(w.node(u, v, 1)))
                    .copied()
                    .unwrap_or([0.0; 3]);
                for a in 0..3 {
                    p[a] += weight * (w.control[s][a] - above[a]);
                }
            }
            let Some((r, rr)) = stencil.derivatives(x, b, c) else {
                continue;
            };
            let g = metric(&r);
            let centre = x[id];
            let mut num = Vec3::ZERO;
            let mut den = 0.0;
            for a in 0..3 {
                // r_aa = (r₊ + r₋) − 2 r: the −2 r g_aa term goes left
                num = num + (rr[a][a] + centre * 2.0) * g[a][a] + r[a] * (g[a][a] * p[a]);
                den += 2.0 * g[a][a];
                for e in a + 1..3 {
                    num = num + rr[a][e] * (2.0 * g[a][e]);
                }
            }
            if den.is_nan() || den <= 0.0 {
                return Err(format!("the metric degenerates at node {id} of block {b}"));
            }
            let new = num / den;
            change = change.max((new - centre).length());
            x[id] = new;
        }
        for &(id, prev, next) in &edges {
            let nb = &neighbours[&id];
            let mut avg = Vec3::ZERO;
            for &q in nb {
                avg = avg + x[q];
            }
            let avg = avg / nb.len() as f64;
            let t = x[next] - x[prev];
            let t = t / t.length();
            let d = avg - x[id];
            let new = x[id] + (d - t * d.dot(t));
            change = change.max((new - x[id]).length());
            x[id] = new;
        }
        report.last_change = change / h_min;
        if !report.last_change.is_finite() {
            return Err("elliptic smoothing diverged".into());
        }
        report.converged = report.last_change < o.tolerance;
    }
    let (max_after, mean_after) = wall_angles(x, &blocks, &walls);
    report.wall_angle_max_after = max_after;
    report.wall_angle_mean_after = mean_after;
    for w in &walls {
        let b = &blocks[w.block];
        for v in 1..w.n[w.t[1]] {
            for u in 1..w.n[w.t[0]] {
                let h = (x[b.at(w.node(u, v, 1))] - x[b.at(w.node(u, v, 0))]).length();
                report.first_spacing_change = report
                    .first_spacing_change
                    .max((h / w.spacing[w.slot(u, v)] - 1.0).abs());
            }
        }
    }
    Ok(report)
}

/// The wall controls: the system evaluated at the wall with the target
/// off-wall derivative, solved for (P₁, P₂, P₃) and relaxed toward it.
fn update_wall(x: &[Vec3], b: &Block, w: &mut Wall, o: &SmoothOptions) {
    let (a, ta, tb) = (w.axis, w.t[0], w.t[1]);
    let sign = if w.side == 0 { 1.0 } else { -1.0 };
    for v in 1..w.n[tb] {
        for u in 1..w.n[ta] {
            let s = w.slot(u, v);
            let node = |d: usize| x[b.at(w.node(u, v, d))];
            let shift =
                |du: isize, dv: isize| ((u as isize + du) as usize, (v as isize + dv) as usize);
            let tan = |du: isize, dv: isize| {
                let (uu, vv) = shift(du, dv);
                x[b.at(w.node(uu, vv, 0))]
            };
            let target = |du: isize, dv: isize| {
                let (uu, vv) = shift(du, dv);
                w.target[w.slot(uu, vv)]
            };
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
}
