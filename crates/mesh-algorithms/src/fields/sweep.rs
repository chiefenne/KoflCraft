//! Fill a sweep region of a [region proposal](super::regions) with a
//! swept O-grid (BS-18, block-dominant meshing).
//!
//! The stations of the sweep are level sets of the flow potential p. Level
//! sets of a harmonic function never cross, they meet every wall at right
//! angles (no flux), and in a torus bend p is the bend angle, so its level
//! sets are exactly the planes of the elbow template. The method follows
//! the harmonic volume parameterisation of generalised cylinders (Martin,
//! Cohen & Kirby 2008, *Volumetric parameterization and trivariate B-spline
//! fitting using harmonic functions*, SPM; Li, Guo & Qin 2007, *Harmonic
//! volumetric mapping for solid modeling applications*, SPM) and
//! skeleton-driven sweeping of tubular shapes (Livesu, Muntoni, Puppo &
//! Scateni 2016, *Skeleton-driven adaptive hexahedral meshing of tubular
//! shapes*, CGF 35(7)), with the potential's level sets in place of planes
//! normal to the skeleton (planes of a tight bend cross inside the duct;
//! level sets do not).
//!
//! 1. Each station is the closed loop where the level set meets the wall
//!    (marching triangles on the scaffold wall), projected onto the STL.
//! 2. The section is the O-grid of the pipe template ("butterfly"): an
//!    H-core and four ring blocks. In a section with four sharp corners (a
//!    rectangular duct) the ring corners sit on them; in a smooth section
//!    they sit at 45° + 90° m of a rotation-minimising frame (double
//!    reflection, Wang, Jüttler, Zheng & Liu 2008, *Computation of rotation
//!    minimizing frames*, ACM TOG 27(1)) carried along the stations.
//! 3. The first rows of the ring are offsets of the wall at the layer
//!    heights (first height, growth), the rest grows to the core.
//! 4. Interior nodes are moved along the station normal onto the station's
//!    level set (a bracketed root search, at most part of the way to the
//!    next station); at ports and cuts the section is planar and every node
//!    is put into its plane, so a cut can join the polyhedral mesher
//!    exactly. A station is pinned to every sharp wall edge across the flow
//!    (a mitre corner); its loop is snapped onto the edges, and when it is
//!    flat the station is the mitre plane.
//! 5. Neighbouring sections are joined node by node into hexahedra.

#![allow(clippy::needless_range_loop)] // parallel per-node arrays

use super::{Role, TetMesh, submesh, tet_gradients};
use crate::remesh::ClosestSet;
use crate::sizing::surface::Primitive;
use mesh_core::geometry::Vec3;
use std::collections::HashMap;

/// The settings of a sweep fill.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SweepOptions {
    /// Cells around a section (a multiple of 4, at least 8).
    pub cells_around: usize,
    /// Height of the first cell at the wall; NaN: 0.5 % of the first
    /// section's diameter.
    pub first_height: f64,
    /// Growth of the layer rows.
    pub growth: f64,
    /// Rows at the layer heights (wall offsets).
    pub layers: usize,
    /// Rows from the wall to the core; `None`: the fewest rows with a growth
    /// of at most [`MAX_CORE_GROWTH`] after the layers.
    pub radial_cells: Option<usize>,
    /// Corners of the core at this share of the way from the section's
    /// centre to the ring corners.
    pub core_fraction: f64,
    /// Axial cell length; `None`: the wall spacing around the first section.
    pub axial_size: Option<f64>,
    /// Largest distance of a cut's level set from its plane, relative to the
    /// section's diameter, for the cut to count as calm.
    pub calm_flatness: f64,
    /// Largest angle between a cut's plane and the plane normal to the
    /// centreline, for the cut to count as calm.
    pub calm_tilt_deg: f64,
}

impl Default for SweepOptions {
    fn default() -> Self {
        Self {
            cells_around: 32,
            first_height: 0.01,
            growth: 1.2,
            layers: 5,
            radial_cells: None,
            core_fraction: 0.5,
            axial_size: None,
            calm_flatness: 0.001,
            calm_tilt_deg: 0.5,
        }
    }
}

/// Largest growth between neighbouring rows from the layers to the core
/// when the number of rows is chosen automatically.
pub const MAX_CORE_GROWTH: f64 = 1.25;

/// Level sets sampled along a sweep to find its cuts and its length.
const SAMPLES: usize = 96;

/// How a sweep ends.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum End {
    /// At an inlet or outlet port (p = 0 or 1).
    Port,
    /// At a calm level set before a region of another kind.
    Cut,
}

/// A station: one level set of the potential and its section.
#[derive(Clone, Debug)]
pub struct Station {
    pub level: f64,
    pub centroid: Vec3,
    /// Unit normal along increasing p.
    pub normal: Vec3,
    pub area: f64,
    /// Largest distance of the wall loop from its best plane over the
    /// section's diameter.
    pub flatness: f64,
    /// Sharp corners of the wall loop.
    pub corners: usize,
    /// The nodes were put into the best plane (port or cut).
    pub planar: bool,
}

/// A boundary quadrilateral of a sweep mesh.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    Wall,
    /// The first section (lowest p).
    Start,
    /// The last section.
    Finish,
}

/// A swept O-grid.
#[derive(Clone, Debug)]
pub struct SweepMesh {
    pub points: Vec<Vec3>,
    /// VTK order: a section quad counter-clockwise about the sweep
    /// direction, then the same quad on the next section.
    pub hexes: Vec<[usize; 8]>,
    /// Outward boundary quadrilaterals with their side and owner.
    pub boundary: Vec<([usize; 4], Side, usize)>,
    pub stations: Vec<Station>,
    pub ends: [End; 2],
    /// Whether each cut met the calm limits (ports always do).
    pub calm_ends: [bool; 2],
    /// Stations pinned to a sharp wall edge across the flow (mitre corners).
    pub feature_stations: usize,
    /// Sharp cross-flow wall edges no station could follow (the edge is not
    /// close to one level set).
    pub missed_features: usize,
    /// Cells around, rows from wall to core, cells along a core side.
    pub section: [usize; 3],
    /// Layer heights from the wall (first is 0).
    pub layer_heights: Vec<f64>,
    pub min_scaled_jacobian: f64,
    pub volume: f64,
    /// Volume of the region's tetrahedra between the two end levels.
    pub region_volume: f64,
    /// Interior section nodes no root search within reach could put on
    /// their level set (left where the section template put them).
    pub unlifted: usize,
}

/// A station candidate: wall loop, normal, centroid, area, flatness.
type Probe = (Vec<Vec3>, Vec3, Vec3, f64, f64);

/// A level-set loop on the wall.
struct Loop {
    points: Vec<Vec3>,
    closed: bool,
}

/// The loops where the level set `level` of `f` crosses the triangles.
/// A node exactly at the level counts as above it.
fn level_loops(points: &[Vec3], tris: &[[usize; 3]], f: &[f64], level: f64) -> Vec<Loop> {
    // segment: two crossed edges (sorted node pairs) and the triangle
    let mut segs: Vec<([[usize; 2]; 2], usize)> = Vec::new();
    for (t, tri) in tris.iter().enumerate() {
        let above = tri.map(|v| f[v] >= level);
        if above[0] == above[1] && above[1] == above[2] {
            continue;
        }
        let mut cut = Vec::with_capacity(2);
        for i in 0..3 {
            let (a, b) = (tri[i], tri[(i + 1) % 3]);
            if above[i] != above[(i + 1) % 3] {
                cut.push(if a < b { [a, b] } else { [b, a] });
            }
        }
        if let [e0, e1] = cut[..] {
            segs.push(([e0, e1], t));
        }
    }
    let at = |e: [usize; 2]| {
        let (a, b) = (e[0], e[1]);
        let (sa, sb) = (f[a] - level, f[b] - level);
        let t = sa / (sa - sb);
        points[a] + (points[b] - points[a]) * t
    };
    let mut by_edge: HashMap<[usize; 2], Vec<usize>> = HashMap::new();
    for (s, (e, _)) in segs.iter().enumerate() {
        for k in e {
            by_edge.entry(*k).or_default().push(s);
        }
    }
    let mut used = vec![false; segs.len()];
    let mut loops = Vec::new();
    for start in 0..segs.len() {
        if used[start] {
            continue;
        }
        // walk forward from edge 1 of the start segment, then backward
        let mut chain: Vec<([usize; 2], usize)> = vec![(segs[start].0[0], segs[start].1)];
        used[start] = true;
        let mut closed = false;
        let mut edge = segs[start].0[1];
        let mut tri = segs[start].1;
        loop {
            if edge == segs[start].0[0] {
                closed = true;
                break;
            }
            chain.push((edge, tri));
            let next = by_edge
                .get(&edge)
                .and_then(|v| v.iter().copied().find(|&s| !used[s]));
            let Some(s) = next else { break };
            used[s] = true;
            tri = segs[s].1;
            edge = if segs[s].0[0] == edge {
                segs[s].0[1]
            } else {
                segs[s].0[0]
            };
        }
        if !closed {
            // extend backward from the start edge
            let mut back: Vec<([usize; 2], usize)> = Vec::new();
            let mut edge = segs[start].0[0];
            while let Some(s) = by_edge
                .get(&edge)
                .and_then(|v| v.iter().copied().find(|&s| !used[s]))
            {
                used[s] = true;
                edge = if segs[s].0[0] == edge {
                    segs[s].0[1]
                } else {
                    segs[s].0[0]
                };
                back.push((edge, segs[s].1));
            }
            back.reverse();
            back.extend(chain);
            chain = back;
        }
        let mut pts: Vec<Vec3> = Vec::with_capacity(chain.len());
        for (e, _) in chain {
            let x = at(e);
            if pts.last().is_some_and(|&q| (q - x).length() == 0.0) {
                continue;
            }
            pts.push(x);
        }
        if closed && pts.len() > 1 && (pts[0] - pts[pts.len() - 1]).length() == 0.0 {
            pts.pop();
        }
        loops.push(Loop {
            points: pts,
            closed,
        });
    }
    loops
}

/// Newell normal (twice the area vector), area-weighted centroid.
fn polygon(points: &[Vec3]) -> (Vec3, Vec3, f64) {
    let n = points.len();
    let mean = points.iter().fold(Vec3::ZERO, |s, &v| s + v) / n as f64;
    let mut area_vec = Vec3::ZERO;
    let mut moment = Vec3::ZERO;
    let mut weight = 0.0;
    for i in 0..n {
        let (a, b) = (points[i], points[(i + 1) % n]);
        let tri = (a - mean).cross(b - mean);
        area_vec = area_vec + tri;
        weight += tri.length();
        moment = moment + (mean + a + b) / 3.0 * tri.length();
    }
    let len = area_vec.length();
    let normal = if len > 0.0 { area_vec / len } else { area_vec };
    let centroid = if weight > 0.0 { moment / weight } else { mean };
    (normal, centroid, len / 2.0)
}

fn unit(v: Vec3) -> Vec3 {
    let l = v.length();
    if l > 0.0 { v / l } else { v }
}

/// Section metrics of a closed loop: normal (oriented along `along`),
/// centroid, area, flatness.
fn metrics(points: &[Vec3], along: Vec3) -> (Vec3, Vec3, f64, f64) {
    let (mut n, c, area) = polygon(points);
    if n.dot(along) < 0.0 {
        n = n * -1.0;
    }
    let diameter = 2.0 * (area / std::f64::consts::PI).sqrt();
    let off = points
        .iter()
        .map(|&p| (p - c).dot(n).abs())
        .fold(0.0, f64::max);
    let flatness = if diameter > 0.0 {
        off / diameter
    } else {
        f64::INFINITY
    };
    (n, c, area, flatness)
}

/// Point location in a set of tetrahedra on a uniform grid of their
/// bounding boxes.
struct Locator<'a> {
    mesh: &'a TetMesh,
    origin: Vec3,
    cell: f64,
    dims: [usize; 3],
    bins: Vec<Vec<usize>>,
}

impl<'a> Locator<'a> {
    fn new(mesh: &'a TetMesh, h: f64) -> Self {
        let mut lo = Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut hi = lo * -1.0;
        for p in &mesh.points {
            lo = Vec3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Vec3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        let ext = hi - lo;
        let cell = (2.0 * h).max(ext.max_abs_component() / 200.0).max(1e-12);
        let dims = [
            ((ext.x / cell) as usize + 1).min(400),
            ((ext.y / cell) as usize + 1).min(400),
            ((ext.z / cell) as usize + 1).min(400),
        ];
        let mut bins = vec![Vec::new(); dims[0] * dims[1] * dims[2]];
        let index =
            |v: f64, o: f64, d: usize| (((v - o) / cell).floor().max(0.0) as usize).min(d - 1);
        for (t, tet) in mesh.tets.iter().enumerate() {
            let ps = tet.map(|v| mesh.points[v]);
            let mut a = ps[0];
            let mut b = ps[0];
            for p in &ps[1..] {
                a = Vec3::new(a.x.min(p.x), a.y.min(p.y), a.z.min(p.z));
                b = Vec3::new(b.x.max(p.x), b.y.max(p.y), b.z.max(p.z));
            }
            for i in index(a.x, lo.x, dims[0])..=index(b.x, lo.x, dims[0]) {
                for j in index(a.y, lo.y, dims[1])..=index(b.y, lo.y, dims[1]) {
                    for k in index(a.z, lo.z, dims[2])..=index(b.z, lo.z, dims[2]) {
                        bins[(i * dims[1] + j) * dims[2] + k].push(t);
                    }
                }
            }
        }
        Self {
            mesh,
            origin: lo,
            cell,
            dims,
            bins,
        }
    }

    /// The tetrahedron holding `x` and its barycentric coordinates.
    fn locate(&self, x: Vec3) -> Option<(usize, [f64; 4])> {
        let rel = (x - self.origin) / self.cell;
        let ijk = [rel.x, rel.y, rel.z];
        let mut idx = [0usize; 3];
        for a in 0..3 {
            if !(ijk[a] >= -1e-9 && ijk[a] < self.dims[a] as f64 + 1e-9) {
                return None;
            }
            idx[a] = (ijk[a].floor().max(0.0) as usize).min(self.dims[a] - 1);
        }
        let mut best: Option<(usize, [f64; 4], f64)> = None;
        for &t in &self.bins[(idx[0] * self.dims[1] + idx[1]) * self.dims[2] + idx[2]] {
            let ps = self.mesh.tets[t].map(|v| self.mesh.points[v]);
            let Some((g, _)) = tet_gradients(ps) else {
                continue;
            };
            let mut b = [0.0; 4];
            for a in 0..4 {
                b[a] = if a == 0 { 0.0 } else { g[a].dot(x - ps[0]) };
            }
            b[0] = 1.0 - b[1] - b[2] - b[3];
            let worst = b.iter().copied().fold(f64::INFINITY, f64::min);
            if best.as_ref().is_none_or(|q| worst > q.2) {
                best = Some((t, b, worst));
            }
        }
        best.filter(|q| q.2 >= -1e-9).map(|q| (q.0, q.1))
    }
}

/// Geometric fractions of `n` cells from a first step `first` over `length`:
/// cumulative positions 0..=length.
fn geometric(first: f64, length: f64, n: usize) -> Vec<f64> {
    let total = |r: f64| {
        if (r - 1.0).abs() < 1e-12 {
            first * n as f64
        } else {
            first * (r.powi(n as i32) - 1.0) / (r - 1.0)
        }
    };
    let (mut lo, mut hi): (f64, f64) = (1e-3, 1e3);
    for _ in 0..200 {
        let mid = (lo * hi).sqrt();
        if total(mid) < length {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let r = (lo * hi).sqrt();
    let mut out = Vec::with_capacity(n + 1);
    let mut s = 0.0;
    let mut step = first;
    out.push(0.0);
    for _ in 0..n {
        s += step;
        step *= r;
        out.push(s);
    }
    let scale = length / s;
    out.iter_mut().for_each(|v| *v *= scale);
    out
}

/// The growth of a geometric law of `n` cells from `first` over `length`
/// (one cell: its length over `first`).
fn growth_of(first: f64, length: f64, n: usize) -> f64 {
    if n < 2 {
        return length / first;
    }
    let pos = geometric(first, length, n);
    (pos[2] - pos[1]) / (pos[1] - pos[0])
}

/// The section nodes of one station (see [`section_nodes`] for the order).
struct SectionGrid {
    nodes: Vec<Vec3>,
    /// Node is on the wall.
    wall: Vec<bool>,
}

/// Index of a node in a section: core (n+1)^2, then four ring blocks of
/// n columns × rows (the last row is the core boundary, the last column
/// the next block's first).
#[derive(Clone, Copy)]
struct Layout {
    n: usize,
    rows: usize,
}

impl Layout {
    fn count(&self) -> usize {
        (self.n + 1) * (self.n + 1) + 4 * self.n * self.rows
    }
    fn core(&self, a: usize, b: usize) -> usize {
        a * (self.n + 1) + b
    }
    /// Core boundary node under column i of ring block m.
    fn core_side(&self, m: usize, i: usize) -> usize {
        let n = self.n;
        match m % 4 {
            0 => self.core(i, 0),
            1 => self.core(n, i),
            2 => self.core(n - i, n),
            _ => self.core(0, n - i),
        }
    }
    /// Ring node (block m, column i 0..=n, row k 0..=rows; row 0 the wall).
    fn ring(&self, m: usize, i: usize, k: usize) -> usize {
        if i == self.n {
            return self.ring(m + 1, 0, k);
        }
        if k == self.rows {
            return self.core_side(m, i);
        }
        (self.n + 1) * (self.n + 1) + ((m % 4) * self.n + i) * self.rows + k
    }
    /// Section quads counter-clockwise about the sweep direction.
    fn quads(&self) -> Vec<[usize; 4]> {
        let mut q = Vec::with_capacity(self.n * self.n + 4 * self.n * self.rows);
        for m in 0..4 {
            for i in 0..self.n {
                for k in 0..self.rows {
                    q.push([
                        self.ring(m, i, k),
                        self.ring(m, i + 1, k),
                        self.ring(m, i + 1, k + 1),
                        self.ring(m, i, k + 1),
                    ]);
                }
            }
        }
        for a in 0..self.n {
            for b in 0..self.n {
                q.push([
                    self.core(a, b),
                    self.core(a + 1, b),
                    self.core(a + 1, b + 1),
                    self.core(a, b + 1),
                ]);
            }
        }
        q
    }
}

/// The wall loop of one station, ready for the section grid.
struct Ring {
    level: f64,
    points: Vec<Vec3>,
    normal: Vec3,
    centroid: Vec3,
    area: f64,
    flatness: f64,
    /// Indices of sharp corners in `points`.
    corners: Vec<usize>,
    planar: bool,
}

/// Indices of sharp corners of a closed polygon: local maxima of the
/// turning angle over a window of `window` of the perimeter, above 45°.
fn sharp_corners(points: &[Vec3], normal: Vec3) -> Vec<usize> {
    let n = points.len();
    if n < 4 {
        return Vec::new();
    }
    let mut s = vec![0.0; n + 1];
    for i in 0..n {
        s[i + 1] = s[i] + (points[(i + 1) % n] - points[i]).length();
    }
    let perimeter = s[n];
    let window = 0.03 * perimeter;
    let at = |t: f64| -> Vec3 {
        let t = t.rem_euclid(perimeter);
        let i = s.partition_point(|&v| v <= t).clamp(1, n) - 1;
        let len = s[i + 1] - s[i];
        let u = if len > 0.0 { (t - s[i]) / len } else { 0.0 };
        points[i] + (points[(i + 1) % n] - points[i]) * u
    };
    let turn: Vec<f64> = (0..n)
        .map(|i| {
            let p = points[i];
            let a = unit(p - at(s[i] - window));
            let b = unit(at(s[i] + window) - p);
            let c = a.cross(b).dot(normal);
            c.atan2(a.dot(b)).abs()
        })
        .collect();
    let mut out = Vec::new();
    for i in 0..n {
        if turn[i] < 45f64.to_radians() {
            continue;
        }
        // the largest turn within the window on both sides
        let mut is_max = true;
        for j in 0..n {
            if j == i {
                continue;
            }
            let d = (s[j] - s[i]).abs();
            let d = d.min(perimeter - d);
            if d <= window && (turn[j] > turn[i] || (turn[j] == turn[i] && j < i)) {
                is_max = false;
                break;
            }
        }
        if is_max {
            out.push(i);
        }
    }
    out
}

/// Project `p` onto the STL, and for a planar station into the plane; a
/// few alternations converge onto the line where both meet.
fn on_wall(stl: &ClosestSet, p: Vec3, plane: Option<(Vec3, Vec3)>) -> Vec3 {
    let mut x = p;
    let rounds = if plane.is_some() { 4 } else { 1 };
    for _ in 0..rounds {
        if let Some((q, _)) = stl.closest(x) {
            x = q;
        }
        if let Some((c, n)) = plane {
            x = x - n * (x - c).dot(n);
        }
    }
    x
}

/// The point at arc length `t` along a closed polyline with cumulative
/// lengths `s`.
fn along(points: &[Vec3], s: &[f64], t: f64) -> Vec3 {
    let n = points.len();
    let perimeter = s[n];
    let t = t.rem_euclid(perimeter);
    let i = s.partition_point(|&v| v <= t).clamp(1, n) - 1;
    let len = s[i + 1] - s[i];
    let u = if len > 0.0 { (t - s[i]) / len } else { 0.0 };
    points[i] + (points[(i + 1) % n] - points[i]) * u
}

/// Where the ray from `c` along `dir` (in the plane) leaves the polygon,
/// as an arc length.
fn ray_hit(points: &[Vec3], s: &[f64], c: Vec3, dir: Vec3, normal: Vec3) -> Option<f64> {
    let n = points.len();
    let side = normal.cross(dir);
    let mut best: Option<(f64, f64)> = None;
    for i in 0..n {
        let (a, b) = (points[i], points[(i + 1) % n]);
        let (da, db) = ((a - c).dot(side), (b - c).dot(side));
        if (da > 0.0) == (db > 0.0) {
            continue;
        }
        let u = da / (da - db);
        let x = a + (b - a) * u;
        let dist = (x - c).dot(dir);
        if dist <= 0.0 {
            continue;
        }
        if best.is_none_or(|q| dist < q.0) {
            best = Some((dist, s[i] + u * (s[i + 1] - s[i])));
        }
    }
    best.map(|q| q.1)
}

/// One step of the double-reflection rotation-minimising frame from the
/// frame (x0, t0, r0) to the point x1 with tangent t1: the new reference.
fn transport(x0: Vec3, t0: Vec3, r0: Vec3, x1: Vec3, t1: Vec3) -> Vec3 {
    let v1 = x1 - x0;
    let c1 = v1.dot(v1);
    if c1 <= 0.0 {
        return r0;
    }
    let rl = r0 - v1 * (2.0 / c1 * v1.dot(r0));
    let tl = t0 - v1 * (2.0 / c1 * v1.dot(t0));
    let v2 = t1 - tl;
    let c2 = v2.dot(v2);
    let r1 = if c2 > 0.0 {
        rl - v2 * (2.0 / c2 * v2.dot(rl))
    } else {
        rl
    };
    unit(r1 - t1 * r1.dot(t1))
}

/// The section nodes of a station: wall loop sides between the four ring
/// corners, layer rows as wall offsets, the rest grown to the core.
#[allow(clippy::too_many_arguments)]
fn section_nodes(
    ring: &Ring,
    stl: &ClosestSet,
    corner_at: [f64; 4],
    sharp: bool,
    layout: Layout,
    heights: &[f64],
    opts: &SweepOptions,
) -> Result<SectionGrid, String> {
    let n = layout.n;
    let pts = &ring.points;
    let m_pts = pts.len();
    let mut s = vec![0.0; m_pts + 1];
    for i in 0..m_pts {
        s[i + 1] = s[i] + (pts[(i + 1) % m_pts] - pts[i]).length();
    }
    let perimeter = s[m_pts];
    let plane = ring.planar.then_some((ring.centroid, ring.normal));
    // wall nodes: 4 n around, side m from corner m to corner m + 1
    let mut wall = Vec::with_capacity(4 * n);
    for m in 0..4 {
        let a = corner_at[m];
        let mut b = corner_at[(m + 1) % 4];
        if b <= a {
            b += perimeter;
        }
        for i in 0..n {
            let t = a + (b - a) * i as f64 / n as f64;
            wall.push(on_wall(stl, along(pts, &s, t), plane));
        }
    }
    let total = 4 * n;
    let normal = ring.normal;
    // inward offset direction at each wall node (unit distance from the wall)
    let mut offset = vec![Vec3::ZERO; total];
    for i in 0..total {
        let prev = wall[(i + total - 1) % total];
        let next = wall[(i + 1) % total];
        let here = wall[i];
        let corner = i % n == 0;
        if corner && sharp {
            let na = unit(normal.cross(here - prev));
            let nb = unit(normal.cross(next - here));
            let d = 1.0 + na.dot(nb);
            if d <= 1e-6 {
                return Err("a section corner folds back on itself".into());
            }
            offset[i] = (na + nb) / d;
        } else {
            offset[i] = unit(normal.cross(next - prev));
        }
    }
    let c = ring.centroid;
    let corner_pts: Vec<Vec3> = (0..4).map(|m| wall[m * n]).collect();
    let core_pts: Vec<Vec3> = corner_pts
        .iter()
        .map(|&o| c + (o - c) * opts.core_fraction)
        .collect();
    let mut nodes = vec![Vec3::ZERO; layout.count()];
    let mut on_wall_flag = vec![false; layout.count()];
    // core: bilinear in the core quadrilateral q0 q1 q2 q3
    for a in 0..=n {
        for b in 0..=n {
            let (u, v) = (a as f64 / n as f64, b as f64 / n as f64);
            let x = core_pts[0] * ((1.0 - u) * (1.0 - v))
                + core_pts[1] * (u * (1.0 - v))
                + core_pts[2] * (u * v)
                + core_pts[3] * ((1.0 - u) * v);
            nodes[layout.core(a, b)] = x;
        }
    }
    let rows = layout.rows;
    let nl = heights.len() - 1;
    for m in 0..4 {
        for i in 0..n {
            let w = wall[m * n + i];
            let d = offset[m * n + i];
            let inner = nodes[layout.core_side(m, i)];
            nodes[layout.ring(m, i, 0)] = w;
            on_wall_flag[layout.ring(m, i, 0)] = true;
            let top = w + d * heights[nl];
            for k in 1..=nl.min(rows - 1) {
                nodes[layout.ring(m, i, k)] = w + d * heights[k];
            }
            let rest = rows - nl;
            let length = (inner - top).length();
            if (inner - top).dot(d) <= 0.0 || length <= 0.0 {
                return Err(format!(
                    "the wall layers ({:.3e} thick) reach the core of the section at p = {:.4}; use fewer or thinner layers or a smaller core_fraction",
                    heights[nl], ring.level
                ));
            }
            let first = if nl == 0 {
                opts.first_height
            } else {
                (heights[nl] - heights[nl - 1]) * opts.growth
            };
            let pos = geometric(first, length, rest);
            for k in 1..rest {
                nodes[layout.ring(m, i, nl + k)] = top + (inner - top) * (pos[k] / length);
            }
        }
    }
    Ok(SectionGrid {
        nodes,
        wall: on_wall_flag,
    })
}

/// Scaled Jacobian of a hexahedron: the smallest corner determinant of
/// unit edge vectors.
pub fn scaled_jacobian(p: [Vec3; 8]) -> f64 {
    const NB: [[usize; 4]; 8] = [
        [0, 1, 3, 4],
        [1, 2, 0, 5],
        [2, 3, 1, 6],
        [3, 0, 2, 7],
        [4, 7, 5, 0],
        [5, 4, 6, 1],
        [6, 5, 7, 2],
        [7, 6, 4, 3],
    ];
    NB.iter()
        .map(|&[c, a, b, d]| {
            let e1 = unit(p[a] - p[c]);
            let e2 = unit(p[b] - p[c]);
            let e3 = unit(p[d] - p[c]);
            e1.dot(e2.cross(e3))
        })
        .fold(f64::INFINITY, f64::min)
}

/// Volume of a hexahedron (12 tetrahedra about its centre, exact for
/// planar faces, the usual fan for bilinear ones).
pub fn hex_volume(p: [Vec3; 8]) -> f64 {
    const FACES: [[usize; 4]; 6] = [
        [0, 4, 7, 3],
        [1, 2, 6, 5],
        [0, 1, 5, 4],
        [3, 7, 6, 2],
        [0, 3, 2, 1],
        [4, 5, 6, 7],
    ];
    let c = p.iter().fold(Vec3::ZERO, |s, &v| s + v) / 8.0;
    let mut v = 0.0;
    for f in FACES {
        let q = f.map(|i| p[i]);
        let fc = (q[0] + q[1] + q[2] + q[3]) / 4.0;
        for k in 0..4 {
            let (a, b) = (q[k], q[(k + 1) % 4]);
            v += (a - c).dot((b - c).cross(fc - c)) / 6.0;
        }
    }
    v
}

/// Fill one sweep region: `region` marks its tetrahedra in `mesh`, `p` is
/// the flow potential on the mesh nodes, `stl` the source wall for
/// projection.
pub fn fill(
    mesh: &TetMesh,
    p: &[f64],
    region: &[bool],
    adjacent: &[[Option<usize>; 4]],
    stl: &ClosestSet,
    opts: SweepOptions,
) -> Result<SweepMesh, String> {
    if opts.cells_around < 8 || !opts.cells_around.is_multiple_of(4) {
        return Err("cells_around must be a multiple of 4 and at least 8".into());
    }
    if !opts.first_height.is_nan() && !(opts.first_height > 0.0 && opts.first_height.is_finite()) {
        return Err("first_height must be positive".into());
    }
    if !(opts.growth >= 1.0 && opts.growth.is_finite()) {
        return Err("growth must be at least 1".into());
    }
    if !(opts.core_fraction > 0.1 && opts.core_fraction < 0.9) {
        return Err("core_fraction must lie between 0.1 and 0.9".into());
    }
    let (sub, old) = submesh(mesh, region, adjacent);
    let f: Vec<f64> = old.iter().map(|&v| p[v]).collect();
    let walls: Vec<[usize; 3]> = sub
        .faces
        .iter()
        .filter(|(_, r)| *r == Role::Wall)
        .map(|(t, _)| *t)
        .collect();
    if walls.is_empty() {
        return Err("the region has no wall".into());
    }
    let has = |role: Role| sub.faces.iter().any(|(_, r)| *r == role);
    let (pmin, pmax) = f
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| {
            (a.min(v), b.max(v))
        });
    let span = pmax - pmin;
    if span.is_nan() || span <= 0.0 {
        return Err("the potential is constant in the region".into());
    }
    // one closed loop at a level, with its metrics
    let probe = |level: f64| -> Option<Probe> {
        let loops = level_loops(&sub.points, &walls, &f, level);
        if loops.len() != 1 || !loops[0].closed || loops[0].points.len() < 3 {
            return None;
        }
        let pts = loops[0].points.clone();
        // orientation: grad p in the tets is not at hand; use the region's
        // mean flow direction from the probes (fixed below)
        let (n, c, a, fl) = metrics(&pts, Vec3::ZERO);
        Some((pts, n, c, a, fl))
    };
    let start_port = has(Role::Inlet) && pmin.abs() < 1e-9;
    let end_port = has(Role::Outlet) && (pmax - 1.0).abs() < 1e-9;
    let eps = 1e-7 * span;
    let samples: Vec<f64> = (0..SAMPLES)
        .map(|j| pmin + span * (j as f64 + 0.5) / SAMPLES as f64)
        .collect();
    let probes: Vec<Option<Probe>> = samples.iter().map(|&l| probe(l)).collect();
    // the flow direction: from the first to the last valid probe centroid
    let valid: Vec<usize> = (0..SAMPLES).filter(|&j| probes[j].is_some()).collect();
    let (Some(&first_ok), Some(&last_ok)) = (valid.first(), valid.last()) else {
        return Err("no level set of the potential crosses the region in one closed loop".into());
    };
    // tilt of each probe: its normal against the centreline tangent
    let tilt: Vec<f64> = (0..SAMPLES)
        .map(|j| {
            let Some(q) = &probes[j] else {
                return f64::INFINITY;
            };
            let c = |k: usize| probes.get(k).and_then(|p| p.as_ref()).map(|p| p.2);
            let (a, b) = (
                if j > 0 { c(j - 1) } else { None }.unwrap_or(q.2),
                c(j + 1).unwrap_or(q.2),
            );
            let t = unit(b - a);
            if t.length() == 0.0 {
                return f64::INFINITY;
            }
            q.1.dot(t).abs().clamp(0.0, 1.0).acos().to_degrees()
        })
        .collect();
    let score = |j: usize| {
        probes[j].as_ref().map_or(f64::INFINITY, |q| {
            (q.4 / opts.calm_flatness).max(tilt[j] / opts.calm_tilt_deg)
        })
    };
    // a cut: the first calm probe from the end, else the calmest probe of
    // that half (reported as not calm)
    let cut_at = |range: Vec<usize>| -> Option<(usize, bool)> {
        if let Some(&j) = range.iter().find(|&&j| score(j) <= 1.0) {
            return Some((j, true));
        }
        let half = &range[..range.len().div_ceil(2)];
        half.iter()
            .copied()
            .filter(|&j| probes[j].is_some())
            .min_by(|&a, &b| score(a).total_cmp(&score(b)))
            .map(|j| (j, false))
    };
    let mut calm_ends = [true, true];
    let lo = if start_port {
        pmin + eps
    } else {
        let (j, calm) =
            cut_at((first_ok..=last_ok).collect()).ok_or("no level set to cut the sweep at")?;
        calm_ends[0] = calm;
        samples[j]
    };
    let hi = if end_port {
        pmax - eps
    } else {
        let (j, calm) = cut_at((first_ok..=last_ok).rev().collect())
            .ok_or("no level set to cut the sweep at")?;
        calm_ends[1] = calm;
        samples[j]
    };
    if hi <= lo {
        return Err("the calm cuts leave nothing to sweep".into());
    }
    // every level between the cuts must give one closed loop
    for (j, &l) in samples.iter().enumerate() {
        if l > lo && l < hi && probes[j].is_none() {
            return Err(format!(
                "the level set p = {l:.4} meets the wall in more than one loop or an open one; the region is not one duct"
            ));
        }
    }
    // centreline length s(p) through the probe centroids (spacing by the
    // largest wall distance between level sets was tried: around a mitre it
    // crowds the stations at the inner corner into slivers)
    let mut curve: Vec<(f64, Vec3)> = Vec::new();
    let mut push = |l: f64| {
        if let Some(q) = probe(l) {
            curve.push((l, q.2));
        }
    };
    push(lo);
    for &l in &samples {
        if l > lo && l < hi {
            push(l);
        }
    }
    push(hi);
    if curve.len() < 2 {
        return Err("the sweep is too short to sample".into());
    }
    let mut s_of = vec![0.0; curve.len()];
    for k in 1..curve.len() {
        s_of[k] = s_of[k - 1] + (curve[k].1 - curve[k - 1].1).length();
    }
    let length = s_of[curve.len() - 1];
    let flow = unit(curve[curve.len() - 1].1 - curve[0].1);
    // first section: its wall spacing sets the default axial size
    let first = probe(lo).ok_or("no section at the start of the sweep")?;
    let mut perimeter = 0.0;
    for i in 0..first.0.len() {
        perimeter += (first.0[(i + 1) % first.0.len()] - first.0[i]).length();
    }
    // first height: given, or 0.5 % of the first section's diameter
    let opts = SweepOptions {
        first_height: if opts.first_height.is_nan() {
            0.01 * (first.3 / std::f64::consts::PI).sqrt()
        } else {
            opts.first_height
        },
        ..opts
    };
    let axial = opts
        .axial_size
        .unwrap_or(perimeter / opts.cells_around as f64);
    if !(axial > 0.0 && axial.is_finite()) {
        return Err("axial_size must be positive".into());
    }
    let cells_axial = ((length / axial).round() as usize).max(1);
    // station levels at equal centreline length
    let level_at = |t: f64| -> f64 {
        let k = s_of.partition_point(|&v| v <= t).clamp(1, s_of.len() - 1);
        let (s0, s1) = (s_of[k - 1], s_of[k]);
        let u = if s1 > s0 { (t - s0) / (s1 - s0) } else { 0.0 };
        curve[k - 1].0 + (curve[k].0 - curve[k - 1].0) * u
    };
    let mut levels: Vec<f64> = (0..=cells_axial)
        .map(|k| level_at(length * k as f64 / cells_axial as f64))
        .collect();
    levels[0] = lo;
    levels[cells_axial] = hi;
    // sharp wall edges across the flow (mitre corners): pin a station to
    // each, so no wall face between two stations cuts the corner
    let h_sub = mesh_h(&sub);
    let tangent_at = |l: f64| -> Vec3 {
        let k = curve
            .partition_point(|q| q.0 <= l)
            .clamp(1, curve.len() - 1);
        unit(curve[k].1 - curve[k - 1].1)
    };
    let mut pinned: Vec<Option<Vec<Primitive>>> = vec![None; levels.len()];
    let mut feature_stations = 0;
    let mut missed_features = 0;
    for ch in cross_flow_chains(&sub.points, &walls, &f, &tangent_at) {
        if ch.mean <= lo || ch.mean >= hi {
            continue;
        }
        let Some(k) = (1..cells_axial).min_by(|&a, &b| {
            (levels[a] - ch.mean)
                .abs()
                .total_cmp(&(levels[b] - ch.mean).abs())
        }) else {
            missed_features += 1;
            continue;
        };
        let gap = 0.5 * (levels[k + 1] - levels[k - 1]);
        if ch.max - ch.min > 0.5 * gap {
            missed_features += 1;
            continue;
        }
        match &mut pinned[k] {
            None => {
                levels[k] = ch.mean;
                pinned[k] = Some(ch.segments);
                feature_stations += 1;
            }
            Some(segs) if (levels[k] - ch.mean).abs() <= 0.25 * gap => segs.extend(ch.segments),
            Some(_) => missed_features += 1,
        }
    }
    if levels.windows(2).any(|w| w[1] <= w[0]) {
        return Err("the stations pinned to sharp wall edges are out of order".into());
    }
    // the wall loops of the stations
    let mut rings: Vec<Ring> = Vec::with_capacity(levels.len());
    for (k, &l) in levels.iter().enumerate() {
        let (pts, _, _, _, _) =
            probe(l).ok_or_else(|| format!("no closed wall loop at p = {l:.4}"))?;
        let along_dir = if k + 1 < levels.len() {
            probe(levels[k + 1]).map(|q| q.2).unwrap_or(Vec3::ZERO)
                - probe(l).map(|q| q.2).unwrap_or(Vec3::ZERO)
        } else {
            flow
        };
        let along_dir = if along_dir.length() > 0.0 {
            along_dir
        } else {
            flow
        };
        // a station pinned to mitre edges: snap the loop onto them first,
        // so its plane runs through the edges
        let mut pts = pts;
        if let Some(segs) = &pinned[k] {
            let edges = ClosestSet::new(segs.clone());
            for x in pts.iter_mut() {
                *x = on_wall(stl, *x, None);
                if let Some((q, d)) = edges.closest(*x)
                    && d <= 0.75 * h_sub
                {
                    *x = q;
                }
            }
            let tiny = 1e-9 * h_sub;
            pts.dedup_by(|a, b| (*a - *b).length() <= tiny);
            while pts.len() > 3 && (pts[0] - pts[pts.len() - 1]).length() <= tiny {
                pts.pop();
            }
        }
        let (n, c, area, flatness) = metrics(&pts, along_dir);
        // ends are planar; so is a station pinned to a mitre edge when its
        // loop is (nearly) the mitre plane
        let planar = k == 0 || k == levels.len() - 1 || (pinned[k].is_some() && flatness <= 0.02);
        // project the loop onto the STL (and the plane of a planar station)
        let plane = planar.then_some((c, n));
        let mut pts: Vec<Vec3> = pts.iter().map(|&x| on_wall(stl, x, plane)).collect();
        // counter-clockwise about n
        let (pn, _, _) = polygon(&pts);
        if pn.dot(n) < 0.0 {
            pts.reverse();
        }
        let corners = sharp_corners(&pts, n);
        rings.push(Ring {
            level: l,
            points: pts,
            normal: n,
            centroid: c,
            area,
            flatness,
            corners,
            planar,
        });
    }
    let sharp = rings[0].corners.len() == 4;
    for r in &rings {
        let k = r.corners.len();
        if !(k == 0 || k == 4) || (k == 4) != sharp {
            return Err(format!(
                "the section at p = {:.4} has {k} sharp corners; a sweep needs a smooth section or one with four corners throughout",
                r.level
            ));
        }
    }
    // ring corners as arc lengths
    let arc = |r: &Ring| -> Vec<f64> {
        let n = r.points.len();
        let mut s = vec![0.0; n + 1];
        for i in 0..n {
            s[i + 1] = s[i] + (r.points[(i + 1) % n] - r.points[i]).length();
        }
        s
    };
    let mut corner_at: Vec<[f64; 4]> = Vec::with_capacity(rings.len());
    if sharp {
        let mut prev: Option<[Vec3; 4]> = None;
        for r in &rings {
            let s = arc(r);
            let c: Vec<usize> = r.corners.clone();
            let mut best_shift = 0;
            if let Some(pv) = prev {
                let mut best = f64::INFINITY;
                for shift in 0..4 {
                    let d: f64 = (0..4)
                        .map(|m| (r.points[c[(m + shift) % 4]] - pv[m]).length())
                        .sum();
                    if d < best {
                        best = d;
                        best_shift = shift;
                    }
                }
            }
            let at: [f64; 4] = std::array::from_fn(|m| s[c[(m + best_shift) % 4]]);
            prev = Some(std::array::from_fn(|m| r.points[c[(m + best_shift) % 4]]));
            corner_at.push(at);
        }
    } else {
        // reference direction: towards the far end of the sweep, or the
        // coordinate axis least along the normal
        let r0 = &rings[0];
        let mut u = (rings[rings.len() - 1].centroid - r0.centroid)
            - r0.normal * (rings[rings.len() - 1].centroid - r0.centroid).dot(r0.normal);
        if u.length() < 0.05 * length.max(1e-300) {
            let axes = [
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
            ];
            let a = axes
                .iter()
                .copied()
                .min_by(|a, b| a.dot(r0.normal).abs().total_cmp(&b.dot(r0.normal).abs()))
                .unwrap_or(axes[0]);
            u = a - r0.normal * a.dot(r0.normal);
        }
        let mut u = unit(u);
        for (k, r) in rings.iter().enumerate() {
            if k > 0 {
                let q = &rings[k - 1];
                u = transport(q.centroid, q.normal, u, r.centroid, r.normal);
            }
            let s = arc(r);
            let v = r.normal.cross(u);
            let mut at = [0.0; 4];
            for m in 0..4 {
                let ang = (45.0 + 90.0 * m as f64).to_radians();
                let dir = u * ang.cos() + v * ang.sin();
                at[m] = ray_hit(&r.points, &s, r.centroid, dir, r.normal).ok_or_else(|| {
                    format!(
                        "the section at p = {:.4} is not star-shaped about its centre",
                        r.level
                    )
                })?;
            }
            // keep the corners in increasing arc length from corner 0
            let perimeter = s[s.len() - 1];
            for m in 1..4 {
                while at[m] <= at[m - 1] {
                    at[m] += perimeter;
                }
            }
            if at[3] - at[0] >= perimeter {
                return Err(format!(
                    "the O-grid corners of the section at p = {:.4} are out of order",
                    r.level
                ));
            }
            corner_at.push(at);
        }
    }
    // layer heights and rows
    let n = opts.cells_around / 4;
    let mut heights = vec![0.0];
    for k in 0..opts.layers {
        heights.push(heights[k] + opts.first_height * opts.growth.powi(k as i32));
    }
    let rows = match opts.radial_cells {
        Some(r) => {
            if r <= opts.layers {
                return Err("radial_cells must exceed layers".into());
            }
            r
        }
        None => {
            // shortest radial run of the first section: wall corner to core
            let r0 = &rings[0];
            let s = arc(r0);
            let shortest = (0..4)
                .map(|m| {
                    let o = along(&r0.points, &s, corner_at[0][m]);
                    let mid = along(
                        &r0.points,
                        &s,
                        0.5 * (corner_at[0][m]
                            + corner_at[0][(m + 1) % 4]
                            + if m == 3 { s[s.len() - 1] } else { 0.0 }),
                    );
                    let core_mid = r0.centroid
                        + ((o - r0.centroid) * opts.core_fraction
                            + (along(&r0.points, &s, corner_at[0][(m + 1) % 4]) - r0.centroid)
                                * opts.core_fraction)
                            / 2.0;
                    (mid - core_mid).length()
                })
                .fold(f64::INFINITY, f64::min);
            let rest = shortest - heights[opts.layers];
            if rest.is_nan() || rest <= 0.0 {
                return Err(
                    "the wall layers fill the whole section; use fewer or thinner layers".into(),
                );
            }
            let first = if opts.layers == 0 {
                opts.first_height
            } else {
                (heights[opts.layers] - heights[opts.layers - 1]) * opts.growth
            };
            let mut extra = 1;
            while extra < 200 && growth_of(first, rest, extra) > MAX_CORE_GROWTH {
                extra += 1;
            }
            opts.layers + extra
        }
    };
    let layout = Layout { n, rows };
    let per = layout.count();
    // section grids; interior nodes of a curved station are moved along
    // the station normal onto its level set, at most 0.45 of the way to the
    // same node of the neighbouring stations (a bracketed root search, so a
    // node never jumps where grad p vanishes, as at the outer corner of a
    // mitre)
    let grids = rings
        .iter()
        .enumerate()
        .map(|(k, r)| section_nodes(r, stl, corner_at[k], sharp, layout, &heights, &opts))
        .collect::<Result<Vec<_>, _>>()?;
    let locator = Locator::new(&sub, mesh_h(&sub));
    let value = |x: Vec3| -> Option<f64> {
        locator
            .locate(x)
            .map(|(t, b)| (0..4).map(|a| b[a] * f[sub.tets[t][a]]).sum::<f64>())
    };
    let mut points = Vec::with_capacity(per * rings.len());
    let mut unlifted = 0;
    for (k, r) in rings.iter().enumerate() {
        let grid = &grids[k];
        for (i, &x0) in grid.nodes.iter().enumerate() {
            if r.planar {
                points.push(x0 - r.normal * (x0 - r.centroid).dot(r.normal));
                continue;
            }
            if grid.wall[i] {
                points.push(x0);
                continue;
            }
            let reach = [k.checked_sub(1), (k + 1 < grids.len()).then_some(k + 1)]
                .into_iter()
                .flatten()
                .map(|j| (grids[j].nodes[i] - x0).length())
                .fold(f64::INFINITY, f64::min)
                * 0.45;
            match lift(&value, x0, r.normal, r.level, reach, 1e-12 * span) {
                Some(x) => points.push(x),
                None => {
                    points.push(x0);
                    unlifted += 1;
                }
            }
        }
    }
    // hexahedra and boundary
    let quads = layout.quads();
    let mut hexes = Vec::with_capacity(quads.len() * (rings.len() - 1));
    let mut boundary = Vec::new();
    let wall_quads: Vec<bool> = (0..quads.len())
        .map(|q| q < 4 * n * rows && (q % rows) == 0)
        .collect();
    for k in 0..rings.len() - 1 {
        let (a, b) = (k * per, (k + 1) * per);
        for (qi, q) in quads.iter().enumerate() {
            let h = [
                a + q[0],
                a + q[1],
                a + q[2],
                a + q[3],
                b + q[0],
                b + q[1],
                b + q[2],
                b + q[3],
            ];
            let cell = hexes.len();
            if wall_quads[qi] {
                boundary.push(([h[0], h[1], h[5], h[4]], Side::Wall, cell));
            }
            if k == 0 {
                boundary.push(([h[0], h[3], h[2], h[1]], Side::Start, cell));
            }
            if k == rings.len() - 2 {
                boundary.push(([h[4], h[5], h[6], h[7]], Side::Finish, cell));
            }
            hexes.push(h);
        }
    }
    let mut min_sj = f64::INFINITY;
    let mut volume = 0.0;
    for h in &hexes {
        let ps = h.map(|v| points[v]);
        min_sj = min_sj.min(scaled_jacobian(ps));
        volume += hex_volume(ps);
    }
    // the region's volume between the end levels (tets split by level)
    let region_volume = slab_volume(&sub, &f, lo, hi);
    let stations = rings
        .iter()
        .map(|r| Station {
            level: r.level,
            centroid: r.centroid,
            normal: r.normal,
            area: r.area,
            flatness: r.flatness,
            corners: r.corners.len(),
            planar: r.planar,
        })
        .collect();
    Ok(SweepMesh {
        points,
        hexes,
        boundary,
        stations,
        ends: [
            if start_port { End::Port } else { End::Cut },
            if end_port { End::Port } else { End::Cut },
        ],
        calm_ends,
        feature_stations,
        missed_features,
        section: [opts.cells_around, rows, n],
        layer_heights: heights,
        min_scaled_jacobian: min_sj,
        volume,
        region_volume,
        unlifted,
    })
}

/// The cap a cut gives the polyhedral mesher of the region behind it
/// (MD-25, BS-17): the cut section from the layer top inwards, each quad
/// split by its Delaunay diagonal (Lawson's angle test), as the planar port
/// cap of that region. A split about the quad centre would give the thin
/// wall rows flat triangles whose circles hold the next row, which no
/// tetrahedralization has as faces: the boundary recovery of the scaffold
/// then splits them without end.
///
/// The polyhedral mesher shrinks a port cap by one in-plane similarity
/// about the mean c of its rim to make room for the wall layers: every cap
/// point x goes to c + (x - c) s with s = 1 - H / r (H the layer stack, r
/// the smallest distance of c from a rim edge), and the layer nodes of a
/// rim point w at height h lie at c + (w - c) (1 - h / r). The cap built
/// here is the preimage of the cut section under that similarity, and the
/// layer rows of the cut section are moved to the mesher's layer nodes, so
/// both sides meet in the same points.
#[derive(Clone, Debug)]
pub struct CutCap {
    /// 0: the start of the sweep, 1: its end.
    pub end: usize,
    /// Unit normal out of the region behind the cut (into the sweep).
    pub normal: Vec3,
    /// Mean of the rim points.
    pub center: Vec3,
    /// Smallest distance of the centre from a rim edge.
    pub radius: f64,
    /// Rim points first (the wall nodes of the cut, in order), then the
    /// inner points, before the similarity.
    pub points: Vec<Vec3>,
    /// The rim points.
    pub rim: usize,
    /// Triangles outward of the region behind the cut.
    pub triangles: Vec<[usize; 3]>,
    /// The sweep node of every cap point.
    pub node: Vec<usize>,
}

impl SweepMesh {
    fn layout(&self) -> Layout {
        Layout {
            n: self.section[0] / 4,
            rows: self.section[1],
        }
    }

    /// The nodes of station `k`: `points[k * per..(k + 1) * per]`.
    pub fn per_station(&self) -> usize {
        self.layout().count()
    }

    /// Prepare the cut at `end` (0 start, 1 end) for the polyhedral region
    /// behind it: move its layer rows onto the polyhedral mesher's layer
    /// nodes and return the cap (see [`CutCap`]).
    pub fn cut_cap(&mut self, end: usize) -> Result<CutCap, String> {
        if end > 1 || self.ends[end] != End::Cut {
            return Err("only a cut gets a cap".into());
        }
        let layout = self.layout();
        let (n, rows) = (layout.n, layout.rows);
        let layers = self.layer_heights.len() - 1;
        if layers >= rows {
            return Err("the layer rows fill the whole section".into());
        }
        let station = if end == 0 { 0 } else { self.stations.len() - 1 };
        let base = station * layout.count();
        let sweep_normal = self.stations[station].normal;
        let normal = if end == 0 {
            sweep_normal
        } else {
            sweep_normal * -1.0
        };
        let wall: Vec<usize> = (0..4)
            .flat_map(|m| (0..n).map(move |i| layout.ring(m, i, 0)))
            .collect();
        let rim: Vec<Vec3> = wall.iter().map(|&v| self.points[base + v]).collect();
        let center = rim.iter().fold(Vec3::ZERO, |s, &p| s + p) / rim.len() as f64;
        let mut radius = f64::INFINITY;
        for i in 0..rim.len() {
            let (a, b) = (rim[i], rim[(i + 1) % rim.len()]);
            let e = b - a;
            radius = radius.min((center - a).cross(e).length() / e.length());
        }
        let stack = self.layer_heights[layers];
        let shrink = 1.0 - stack / radius;
        if !(shrink > 0.1 && radius.is_finite()) {
            return Err(format!(
                "the layer stack ({stack:.3e}) closes the cut (inner radius {radius:.3e})"
            ));
        }
        // the layer rows of the cut: the mesher's layer nodes
        for m in 0..4 {
            for i in 0..n {
                let w = self.points[base + layout.ring(m, i, 0)];
                for k in 1..=layers {
                    let t = 1.0 - self.layer_heights[k] / radius;
                    self.points[base + layout.ring(m, i, k)] = center + (w - center) * t;
                }
            }
        }
        // cap points: the rim, then the section inside the layer top (before
        // the similarity); each quad split by its Delaunay diagonal
        let section: Vec<Vec3> = self.points[base..base + layout.count()].to_vec();
        let mut quads: Vec<([usize; 4], [usize; 4])> = Vec::new();
        for m in 0..4 {
            for i in 0..n {
                for k in layers..rows {
                    let post = [
                        layout.ring(m, i, k),
                        layout.ring(m, i + 1, k),
                        layout.ring(m, i + 1, k + 1),
                        layout.ring(m, i, k + 1),
                    ];
                    // the layer top is the rim before the similarity
                    let pre = if k == layers {
                        [
                            layout.ring(m, i, 0),
                            layout.ring(m, i + 1, 0),
                            post[2],
                            post[3],
                        ]
                    } else {
                        post
                    };
                    quads.push((pre, post));
                }
            }
        }
        for a in 0..n {
            for b in 0..n {
                let q = [
                    layout.core(a, b),
                    layout.core(a + 1, b),
                    layout.core(a + 1, b + 1),
                    layout.core(a, b + 1),
                ];
                quads.push((q, q));
            }
        }
        let before = |x: Vec3| center + (x - center) / shrink;
        let mut points: Vec<Vec3> = rim.clone();
        let mut node: Vec<usize> = wall.iter().map(|&v| base + v).collect();
        let mut index: HashMap<usize, usize> =
            wall.iter().enumerate().map(|(j, &v)| (v, j)).collect();
        let mut triangles = Vec::with_capacity(2 * quads.len());
        for (pre, _) in quads {
            let mut ids = [0; 4];
            for (c, &v) in pre.iter().enumerate() {
                ids[c] = match index.get(&v) {
                    Some(&j) => j,
                    None => {
                        points.push(before(section[v]));
                        node.push(base + v);
                        index.insert(v, points.len() - 1);
                        points.len() - 1
                    }
                };
            }
            // diagonal 0-2 unless the angles at 1 and 3 sum to more than pi
            let angle = |c: usize| {
                let (p, a, b) = (
                    points[ids[c]],
                    points[ids[(c + 3) % 4]],
                    points[ids[(c + 1) % 4]],
                );
                let (u, w) = (a - p, b - p);
                u.cross(w).length().atan2(u.dot(w))
            };
            if angle(1) + angle(3) <= angle(0) + angle(2) {
                triangles.push([ids[0], ids[1], ids[2]]);
                triangles.push([ids[0], ids[2], ids[3]]);
            } else {
                triangles.push([ids[0], ids[1], ids[3]]);
                triangles.push([ids[1], ids[2], ids[3]]);
            }
        }
        if end == 1 {
            for t in triangles.iter_mut() {
                t.swap(1, 2);
            }
        }
        Ok(CutCap {
            end,
            normal,
            center,
            radius,
            points,
            rim: wall.len(),
            triangles,
            node,
        })
    }
}

/// Connected sharp wall edges that run across the flow, with the range and
/// mean of the potential on their nodes.
struct Chain {
    min: f64,
    max: f64,
    mean: f64,
    segments: Vec<Primitive>,
}

/// Sharp wall edges (normals more than 30° apart) more than 60° from the
/// flow direction, grouped into connected chains.
fn cross_flow_chains(
    points: &[Vec3],
    walls: &[[usize; 3]],
    f: &[f64],
    tangent_at: &dyn Fn(f64) -> Vec3,
) -> Vec<Chain> {
    let mut by_edge: HashMap<[usize; 2], Vec<usize>> = HashMap::new();
    for (t, tri) in walls.iter().enumerate() {
        for i in 0..3 {
            let (a, b) = (tri[i], tri[(i + 1) % 3]);
            by_edge
                .entry(if a < b { [a, b] } else { [b, a] })
                .or_default()
                .push(t);
        }
    }
    let normal = |t: usize| {
        let [a, b, c] = walls[t].map(|v| points[v]);
        unit((b - a).cross(c - a))
    };
    let cos_sharp = 30f64.to_radians().cos();
    let mut edges: Vec<[usize; 2]> = by_edge
        .iter()
        .filter(|(e, ts)| {
            if ts.len() != 2 || normal(ts[0]).dot(normal(ts[1])).abs() > cos_sharp {
                return false;
            }
            let dir = unit(points[e[1]] - points[e[0]]);
            dir.dot(tangent_at(0.5 * (f[e[0]] + f[e[1]]))).abs() <= 0.5
        })
        .map(|(e, _)| *e)
        .collect();
    edges.sort_unstable();
    // components by shared nodes
    let mut parent: HashMap<usize, usize> = HashMap::new();
    fn root(parent: &mut HashMap<usize, usize>, mut v: usize) -> usize {
        while let Some(&p) = parent.get(&v) {
            if p == v {
                break;
            }
            v = p;
        }
        v
    }
    for e in &edges {
        for &v in e {
            parent.entry(v).or_insert(v);
        }
        let (a, b) = (root(&mut parent, e[0]), root(&mut parent, e[1]));
        if a != b {
            parent.insert(a.max(b), a.min(b));
        }
    }
    let mut groups: std::collections::BTreeMap<usize, (Vec<usize>, Vec<Primitive>)> =
        std::collections::BTreeMap::new();
    for e in &edges {
        let r = root(&mut parent, e[0]);
        let g = groups.entry(r).or_default();
        g.0.extend_from_slice(e);
        g.1.push(Primitive::Segment([points[e[0]], points[e[1]]]));
    }
    groups
        .into_values()
        .map(|(mut nodes, segments)| {
            nodes.sort_unstable();
            nodes.dedup();
            let vals: Vec<f64> = nodes.iter().map(|&v| f[v]).collect();
            Chain {
                min: vals.iter().copied().fold(f64::INFINITY, f64::min),
                max: vals.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                mean: vals.iter().sum::<f64>() / vals.len() as f64,
                segments,
            }
        })
        .collect()
}

/// Mean edge length of a tet mesh.
/// The point `x0 + dir t`, |t| <= `reach`, where `value` equals `level`:
/// samples outwards on both sides for a sign change, then regula falsi
/// (Illinois). `None` when no crossing lies within reach.
fn lift(
    value: &dyn Fn(Vec3) -> Option<f64>,
    x0: Vec3,
    dir: Vec3,
    level: f64,
    reach: f64,
    tol: f64,
) -> Option<Vec3> {
    let g = |t: f64| value(x0 + dir * t).map(|v| v - level);
    let g0 = g(0.0)?;
    if g0.abs() <= tol {
        return Some(x0);
    }
    if !(reach > 0.0 && reach.is_finite()) {
        return None;
    }
    // p grows along the station normal: search the expected side first
    let first = if g0 > 0.0 { -1.0 } else { 1.0 };
    let mut bracket = None;
    'sides: for side in [first, -first] {
        let (mut ta, mut ga) = (0.0, g0);
        for j in 1..=8 {
            let tb = side * reach * j as f64 / 8.0;
            let Some(gb) = g(tb) else { continue 'sides };
            if gb.abs() <= tol {
                return Some(x0 + dir * tb);
            }
            if (ga < 0.0) != (gb < 0.0) {
                bracket = Some((ta, ga, tb, gb));
                break 'sides;
            }
            (ta, ga) = (tb, gb);
        }
    }
    let (mut a, mut fa, mut b, mut fb) = bracket?;
    let mut last = 0;
    for _ in 0..60 {
        let t = (a * fb - b * fa) / (fb - fa);
        let ft = g(t)?;
        if ft.abs() <= tol || (b - a).abs() <= 1e-14 * reach {
            return Some(x0 + dir * t);
        }
        if (ft < 0.0) == (fa < 0.0) {
            (a, fa) = (t, ft);
            if last == -1 {
                fb *= 0.5;
            }
            last = -1;
        } else {
            (b, fb) = (t, ft);
            if last == 1 {
                fa *= 0.5;
            }
            last = 1;
        }
    }
    Some(x0 + dir * (a * fb - b * fa) / (fb - fa))
}

fn mesh_h(mesh: &TetMesh) -> f64 {
    let mut sum = 0.0;
    let mut n = 0usize;
    for t in &mesh.tets {
        for a in 0..4 {
            for b in a + 1..4 {
                sum += (mesh.points[t[a]] - mesh.points[t[b]]).length();
                n += 1;
            }
        }
    }
    if n > 0 { sum / n as f64 } else { 1.0 }
}

/// Volume of the part of the tetrahedra where lo <= f <= hi (f linear),
/// by the clipped volume of each tetrahedron.
fn slab_volume(mesh: &TetMesh, f: &[f64], lo: f64, hi: f64) -> f64 {
    let mut total = 0.0;
    for t in &mesh.tets {
        let ps = t.map(|v| mesh.points[v]);
        let vs = t.map(|v| f[v]);
        total += below(ps, vs, hi) - below(ps, vs, lo);
    }
    total
}

/// Volume of the part of a linear tetrahedron where f <= level, clipped
/// exactly: one vertex below gives a corner tetrahedron, three below the
/// tetrahedron less a corner, two below a prism of three tetrahedra (the
/// clipped part is convex). Every edge ratio has a positive denominator,
/// so equal vertex values need no special case.
fn below(p: [Vec3; 4], f: [f64; 4], level: f64) -> f64 {
    let tet = |a: Vec3, b: Vec3, c: Vec3, d: Vec3| (b - a).dot((c - a).cross(d - a)).abs() / 6.0;
    let vol = tet(p[0], p[1], p[2], p[3]);
    let mut o = [0, 1, 2, 3];
    o.sort_by(|&a, &b| f[a].total_cmp(&f[b]));
    let s = o.map(|i| f[i]);
    let v = o.map(|i| p[i]);
    let t = level;
    if t <= s[0] {
        return 0.0;
    }
    if t >= s[3] {
        return vol;
    }
    // the point on edge a-b where f = t (f(a) < t <= f(b) or the reverse)
    let cut = |a: usize, b: usize| v[a] + (v[b] - v[a]) * ((t - s[a]) / (s[b] - s[a]));
    if t <= s[1] {
        let r = (1..4).map(|j| (t - s[0]) / (s[j] - s[0])).product::<f64>();
        vol * r
    } else if t >= s[2] {
        let r = (0..3).map(|j| (s[3] - t) / (s[3] - s[j])).product::<f64>();
        vol * (1.0 - r)
    } else {
        // prism: v0, p02, p03 below v1, p12, p13
        let (a0, a1, a2) = (v[0], cut(0, 2), cut(0, 3));
        let (b0, b1, b2) = (v[1], cut(1, 2), cut(1, 3));
        tet(a0, a1, a2, b0) + tet(a1, a2, b0, b1) + tet(a2, b0, b1, b2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_tet() -> [Vec3; 4] {
        [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ]
    }

    /// Keep: the slab volume behind the sweep's volume check is exact for
    /// a linear field, with distinct and with equal vertex values. A wrong
    /// case split or prism shifts every reported volume difference; the
    /// fill tests only bound that difference loosely.
    #[test]
    fn volume_below_a_level_is_exact() {
        let p = unit_tet();
        // f = x + 2y + 3z: {f <= t} is the tet with intercepts t, t/2, t/3
        for t in [0.1, 0.5, 0.9] {
            let v = below(p, [0.0, 1.0, 2.0, 3.0], t);
            assert!((v - t * t * t / 36.0).abs() < 1e-15, "t = {t}: {v}");
        }
        // two vertices below (the prism) against the complement
        for t in [1.2, 1.5, 1.9] {
            let v = below(p, [0.0, 1.0, 2.0, 3.0], t);
            let w = below(p, [0.0, -1.0, -2.0, -3.0], -t);
            assert!((v + w - 1.0 / 6.0).abs() < 1e-15, "t = {t}: {v} + {w}");
        }
        // f = y + z (two pairs of equal values): the integral of
        // (1 - u) u du over u in [0, t] = (3t^2 - 2t^3) / 6
        for t in [0.3, 0.6] {
            let v = below(p, [0.0, 0.0, 1.0, 1.0], t);
            let exact = (3.0 * t * t - 2.0 * t * t * t) / 6.0;
            assert!((v - exact).abs() < 1e-15, "t = {t}: {v} vs {exact}");
        }
        // f = x (three equal values): 1/6 (1 - (1 - t)^3)
        for t in [0.2, 0.7] {
            let v = below(p, [0.0, 1.0, 0.0, 0.0], t);
            let exact = (1.0 - (1.0 - t) * (1.0 - t) * (1.0 - t)) / 6.0;
            assert!((v - exact).abs() < 1e-15, "t = {t}: {v} vs {exact}");
        }
        assert_eq!(below(p, [0.0, 1.0, 2.0, 3.0], -1.0), 0.0);
        assert!((below(p, [0.0, 1.0, 2.0, 3.0], 4.0) - 1.0 / 6.0).abs() < 1e-15);
    }

    /// Keep: a node is moved onto its level set only within reach. A lift
    /// without the bound let nodes jump where grad p vanishes (the outer
    /// corner of a mitre) and inverted cells; the elbow test sees the
    /// result, this one the bound itself.
    #[test]
    fn lift_stays_within_reach() {
        let field = |x: Vec3| -> Option<f64> { (x.x.abs() <= 1.0).then_some(x.x) };
        let x0 = Vec3::new(0.2, 0.3, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let x = lift(&field, x0, dir, 0.5, 0.4, 1e-14).unwrap();
        assert!((x - Vec3::new(0.5, 0.3, 0.0)).length() < 1e-12, "{x:?}");
        // behind the node, on the other side of the search
        let x = lift(&field, x0, dir, -0.1, 0.4, 1e-14).unwrap();
        assert!((x.x + 0.1).abs() < 1e-12, "{x:?}");
        assert!(lift(&field, x0, dir, 0.5, 0.2, 1e-14).is_none());
        // no gradient along the search line: no crossing
        let across = Vec3::new(0.0, 1.0, 0.0);
        assert!(lift(&field, x0, across, 0.5, 0.4, 1e-14).is_none());
    }

    /// Keep: the wall rows start at the first height and end on the core.
    /// A wrong ratio search moves the first row off the layer height the
    /// polyhedral join will rely on; the fill tests do not measure rows.
    #[test]
    fn geometric_rows_keep_first_height_and_length() {
        let pos = geometric(0.1, 2.0, 5);
        assert_eq!(pos.len(), 6);
        assert!((pos[1] - 0.1).abs() < 1e-9, "{pos:?}");
        assert!((pos[5] - 2.0).abs() < 1e-12);
        let r = (pos[2] - pos[1]) / (pos[1] - pos[0]);
        for k in 2..5 {
            let q = (pos[k + 1] - pos[k]) / (pos[k] - pos[k - 1]);
            assert!((q - r).abs() < 1e-9);
        }
        assert!((growth_of(0.1, 2.0, 5) - r).abs() < 1e-12);
    }

    /// Keep: a unit cube has volume 1 and scaled Jacobian 1, and a mirrored
    /// one a negative volume; the fill reports and validation build on
    /// both measures.
    #[test]
    fn hex_measures_of_a_cube() {
        let c = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 1.0),
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(0.0, 1.0, 1.0),
        ];
        assert!((hex_volume(c) - 1.0).abs() < 1e-14);
        assert!((scaled_jacobian(c) - 1.0).abs() < 1e-14);
        let m = c.map(|v| Vec3::new(-v.x, v.y, v.z));
        assert!((hex_volume(m) + 1.0).abs() < 1e-14);
        assert!(scaled_jacobian(m) < 0.0);
    }
}
