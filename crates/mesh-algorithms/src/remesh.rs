//! Adaptive isotropic remeshing of a triangulated wall onto its source
//! surface, with a bound on the chord deviation.
//!
//! The method is the remeshing loop of Botsch & Kobbelt, *A Remeshing
//! Approach to Multiresolution Modeling* (SGP 2004): split edges longer than
//! 4/3 of the target length, collapse edges shorter than 4/5 of it, flip
//! edges towards valence six, move every vertex tangentially to the area
//! centroid of its ring and project it back onto the source surface. The
//! target length follows the curvature as in Dunyach, Vanderhaeghe, Barthe &
//! Botsch, *Adaptive Remeshing for Real-Time Mesh Deformation* (Eurographics
//! 2013): a triangle of edge length L on a sphere of radius r = 1/κ deviates
//! from it by at most ε when L = sqrt(6ε/κ − 3ε²). The curvature is the
//! normal-cycle tensor of the source (`sizing::surface`), and the lengths
//! are gradient limited.
//!
//! Constraints: facets marked `keep` (the port caps) stay exactly as they
//! are, with their rim vertices and edges. Source edges sharper than the
//! feature angle are kept as feature lines: feature edges are split along
//! the line and never flipped, a feature vertex moves only along its line
//! and is projected onto the source feature edges, and feature corners
//! (one, or three or more, feature edges, or a sharp turn of the line) stay
//! fixed. A final pass splits every facet whose sampled distance to the
//! source still exceeds the tolerance, so the reported deviation is the
//! measured one, not only the estimate.

#![allow(clippy::needless_range_loop)] // parallel per-vertex and per-face arrays

use crate::sizing::surface::{
    Primitive, SurfaceGeometry, Weighted, WeightedSet, normal_cycle_curvature,
};
use mesh_core::geometry::Vec3;
use mesh_geometry::patches::closest_point;
use std::collections::{HashMap, HashSet};

/// The settings of a wall remesh.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RemeshOptions {
    /// Largest distance of a facet from the source surface (chord error).
    pub tolerance: f64,
    pub min_size: f64,
    pub max_size: f64,
    /// Gradient limit of the target length (length per unit distance).
    pub rate: f64,
    /// Edges whose facet normals turn more than this (radians) are features.
    pub feature_angle: f64,
    pub iterations: usize,
    /// Target length at the rims of the kept facets (`None`: no rim source).
    pub rim_size: Option<f64>,
}

impl Default for RemeshOptions {
    fn default() -> Self {
        Self {
            tolerance: 0.05,
            min_size: 0.65,
            max_size: 6.0,
            rate: 0.4,
            feature_angle: 30f64.to_radians(),
            iterations: 10,
            rim_size: None,
        }
    }
}

/// A remeshed closed surface.
#[derive(Clone, Debug, PartialEq)]
pub struct Remeshed {
    pub vertices: Vec<Vec3>,
    pub triangles: Vec<[usize; 3]>,
    /// The source facet each kept triangle came from; `None` for new wall.
    pub source: Vec<Option<usize>>,
    pub stats: RemeshStats,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RemeshStats {
    pub feature_edges: usize,
    pub corners: usize,
    pub splits: usize,
    pub collapses: usize,
    pub flips: usize,
    /// Facets split in the final deviation pass.
    pub deviation_splits: usize,
    /// Largest sampled distance of a new facet from the source surface.
    pub max_deviation: f64,
    /// Largest sampled distance of the source wall from the new wall.
    pub max_reverse_deviation: f64,
    pub min_edge: f64,
    pub max_edge: f64,
}

/// Closest-point queries on triangles or segments.
pub struct ClosestSet {
    items: Vec<Primitive>,
    nodes: Vec<Node>,
}

#[derive(Clone, Copy)]
struct Node {
    min: Vec3,
    max: Vec3,
    start: usize,
    end: usize,
    children: Option<(usize, usize)>,
}

fn corners(p: &Primitive) -> &[Vec3] {
    match p {
        Primitive::Triangle(v) => v,
        Primitive::Segment(v) => v,
    }
}

fn closest_on(p: Vec3, primitive: &Primitive) -> Vec3 {
    match *primitive {
        Primitive::Triangle(t) => closest_point(p, t),
        Primitive::Segment([a, b]) => {
            let ab = b - a;
            let t = ((p - a).dot(ab) / ab.dot(ab)).clamp(0.0, 1.0);
            a + ab * t
        }
    }
}

impl ClosestSet {
    pub fn new(mut items: Vec<Primitive>) -> Self {
        let mut nodes = Vec::new();
        if !items.is_empty() {
            let n = items.len();
            build_closest(&mut items, &mut nodes, 0, n);
        }
        Self { items, nodes }
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// The closest point of the set to `p` and its distance.
    pub fn closest(&self, p: Vec3) -> Option<(Vec3, f64)> {
        let mut best: Option<(Vec3, f64)> = None;
        let mut best_sq = f64::INFINITY;
        if self.nodes.is_empty() {
            return None;
        }
        let mut stack = vec![0usize];
        while let Some(i) = stack.pop() {
            let node = self.nodes[i];
            if box_sq(p, node.min, node.max) >= best_sq {
                continue;
            }
            match node.children {
                None => {
                    for item in &self.items[node.start..node.end] {
                        let q = closest_on(p, item);
                        let d = (p - q).dot(p - q);
                        if d < best_sq {
                            best_sq = d;
                            best = Some((q, d.sqrt()));
                        }
                    }
                }
                Some((l, r)) => {
                    let (dl, dr) = (
                        box_sq(p, self.nodes[l].min, self.nodes[l].max),
                        box_sq(p, self.nodes[r].min, self.nodes[r].max),
                    );
                    if dl <= dr {
                        stack.push(r);
                        stack.push(l);
                    } else {
                        stack.push(l);
                        stack.push(r);
                    }
                }
            }
        }
        best
    }
}

fn build_closest(
    items: &mut [Primitive],
    nodes: &mut Vec<Node>,
    start: usize,
    end: usize,
) -> usize {
    let mut min = corners(&items[start])[0];
    let mut max = min;
    for item in &items[start..end] {
        for &v in corners(item) {
            min = Vec3::new(min.x.min(v.x), min.y.min(v.y), min.z.min(v.z));
            max = Vec3::new(max.x.max(v.x), max.y.max(v.y), max.z.max(v.z));
        }
    }
    let id = nodes.len();
    nodes.push(Node {
        min,
        max,
        start,
        end,
        children: None,
    });
    if end - start > 4 {
        let e = max - min;
        let axis = if e.x >= e.y && e.x >= e.z {
            0
        } else if e.y >= e.z {
            1
        } else {
            2
        };
        let key = |p: &Primitive| {
            let s = corners(p).iter().fold(Vec3::ZERO, |s, &v| s + v);
            [s.x, s.y, s.z][axis]
        };
        let mid = start + (end - start) / 2;
        items[start..end].select_nth_unstable_by(mid - start, |a, b| key(a).total_cmp(&key(b)));
        let l = build_closest(items, nodes, start, mid);
        let r = build_closest(items, nodes, mid, end);
        nodes[id].children = Some((l, r));
    }
    id
}

fn box_sq(q: Vec3, min: Vec3, max: Vec3) -> f64 {
    let d = |v: f64, lo: f64, hi: f64| (lo - v).max(v - hi).max(0.0);
    let (x, y, z) = (
        d(q.x, min.x, max.x),
        d(q.y, min.y, max.y),
        d(q.z, min.z, max.z),
    );
    x * x + y * y + z * z
}

/// Edge length at which an equilateral triangle on a sphere of curvature
/// `k` deviates from it by `eps` (Dunyach et al. 2013).
pub fn chord_length(k: f64, eps: f64) -> f64 {
    if k <= 0.0 {
        return f64::INFINITY;
    }
    let s = 6.0 * eps / k - 3.0 * eps * eps;
    if s <= 0.0 { 0.0 } else { s.sqrt() }
}

/// Points at which a facet's distance from the source is sampled: the
/// corners are on the source, so the centroid, edge midpoints and the
/// midpoints between them and the centroid.
pub fn samples(t: [Vec3; 3]) -> [Vec3; 7] {
    let [a, b, c] = t;
    let g = (a + b + c) / 3.0;
    let (ab, bc, ca) = ((a + b) / 2.0, (b + c) / 2.0, (c + a) / 2.0);
    [g, ab, bc, ca, (g + a) / 2.0, (g + b) / 2.0, (g + c) / 2.0]
}

struct Work<'a> {
    p: Vec<Vec3>,
    t: Vec<[usize; 3]>,
    alive: Vec<bool>,
    keep: Vec<bool>,
    vf: Vec<Vec<usize>>,
    locked: Vec<bool>,
    corner: Vec<bool>,
    /// Target edge length at each vertex.
    h: Vec<f64>,
    feature: HashSet<(usize, usize)>,
    surface: &'a ClosestSet,
    lines: &'a ClosestSet,
    size: &'a WeightedSet,
    opts: RemeshOptions,
    stats: RemeshStats,
}

fn key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

fn normal(t: [Vec3; 3]) -> Vec3 {
    (t[1] - t[0]).cross(t[2] - t[0])
}

/// Triangle shape 4√3 · area / Σ l²: 1 for an equilateral triangle, 0 for
/// a degenerate one.
fn shape(t: [Vec3; 3]) -> f64 {
    let l2 = (t[1] - t[0]).dot(t[1] - t[0])
        + (t[2] - t[1]).dot(t[2] - t[1])
        + (t[0] - t[2]).dot(t[0] - t[2]);
    if l2 <= 0.0 {
        return 0.0;
    }
    2.0 * 3f64.sqrt() * normal(t).length() / l2
}

/// New facets must not be worse than this, or than the worst they replace.
const MIN_SHAPE: f64 = 0.1;

fn unit(v: Vec3) -> Vec3 {
    let l = v.length();
    if l > 0.0 { v / l } else { v }
}

impl Work<'_> {
    fn target(&self, x: Vec3) -> f64 {
        self.size
            .envelope(x, self.opts.rate, 0.0, self.opts.max_size)
            .max(self.opts.min_size)
    }

    /// Target length of edge (a, b): the mean of the vertex targets.
    fn edge_target(&self, a: usize, b: usize) -> f64 {
        0.5 * (self.h[a] + self.h[b])
    }

    fn refresh(&mut self, v: usize) {
        let h = self.target(self.p[v]);
        if v < self.h.len() {
            self.h[v] = h;
        } else {
            self.h.push(h);
        }
    }

    fn corners_of(&self, f: usize) -> [Vec3; 3] {
        self.t[f].map(|v| self.p[v])
    }

    fn edge_faces(&self, a: usize, b: usize) -> Vec<usize> {
        self.vf[a]
            .iter()
            .copied()
            .filter(|&f| self.t[f].contains(&b))
            .collect()
    }

    fn neighbours(&self, a: usize) -> Vec<usize> {
        let mut n: Vec<usize> = self.vf[a]
            .iter()
            .flat_map(|&f| self.t[f])
            .filter(|&v| v != a)
            .collect();
        n.sort_unstable();
        n.dedup();
        n
    }

    fn feature_count(&self, a: usize) -> usize {
        self.neighbours(a)
            .into_iter()
            .filter(|&b| self.feature.contains(&key(a, b)))
            .count()
    }

    fn is_feature_vertex(&self, a: usize) -> bool {
        self.feature_count(a) > 0
    }

    fn fixed(&self, a: usize) -> bool {
        self.locked[a] || self.corner[a]
    }

    /// The edges of the wall facets, each once.
    fn edges(&self) -> Vec<(usize, usize)> {
        let mut set = HashSet::new();
        let mut out = Vec::new();
        for (f, t) in self.t.iter().enumerate() {
            if !self.alive[f] || self.keep[f] {
                continue;
            }
            for k in 0..3 {
                let e = key(t[k], t[(k + 1) % 3]);
                if set.insert(e) {
                    out.push(e);
                }
            }
        }
        out
    }

    /// Whether the edge borders a kept facet (a rim edge).
    fn touches_keep(&self, faces: &[usize]) -> bool {
        faces.iter().any(|&f| self.keep[f])
    }

    fn project(&self, a: usize, x: Vec3) -> Vec3 {
        let set = if self.is_feature_vertex(a) {
            self.lines
        } else {
            self.surface
        };
        set.closest(x).map_or(x, |(q, _)| q)
    }

    fn add_face(&mut self, t: [usize; 3]) -> usize {
        let f = self.t.len();
        self.t.push(t);
        self.alive.push(true);
        self.keep.push(false);
        for v in t {
            self.vf[v].push(f);
        }
        f
    }

    fn remove_face(&mut self, f: usize) {
        self.alive[f] = false;
        for v in self.t[f] {
            self.vf[v].retain(|&g| g != f);
        }
    }

    /// Split edge (a, b) at its midpoint projected onto the source.
    fn split(&mut self, a: usize, b: usize) -> bool {
        let faces = self.edge_faces(a, b);
        if faces.len() != 2 || self.touches_keep(&faces) {
            return false;
        }
        let is_feature = self.feature.contains(&key(a, b));
        let mid = (self.p[a] + self.p[b]) / 2.0;
        let m = self.p.len();
        self.p.push(mid);
        self.vf.push(Vec::new());
        self.locked.push(false);
        self.corner.push(false);
        if is_feature {
            self.feature.remove(&key(a, b));
            self.feature.insert(key(a, m));
            self.feature.insert(key(m, b));
        }
        for f in faces {
            let t = self.t[f];
            let i = (0..3).find(|&i| t[i] != a && t[i] != b).unwrap_or(0);
            let (u, v, w) = (t[(i + 1) % 3], t[(i + 2) % 3], t[i]);
            self.remove_face(f);
            self.add_face([u, m, w]);
            self.add_face([m, v, w]);
        }
        let q = self.project(m, mid);
        self.p[m] = q;
        self.refresh(m);
        self.stats.splits += 1;
        true
    }

    /// Collapse vertex a into b, if topology, features and geometry allow.
    fn collapse(&mut self, a: usize, b: usize, hi: f64) -> bool {
        if self.fixed(a) {
            return false;
        }
        let faces = self.edge_faces(a, b);
        if faces.len() != 2 || self.touches_keep(&faces) {
            return false;
        }
        let edge_feature = self.feature.contains(&key(a, b));
        if self.is_feature_vertex(a) && !edge_feature {
            return false;
        }
        if edge_feature && self.feature_count(a) != 2 {
            return false;
        }
        // link condition
        let na = self.neighbours(a);
        let nb = self.neighbours(b);
        let common = na.iter().filter(|v| nb.contains(v)).count();
        if common != 2 {
            return false;
        }
        let target = self.p[b];
        let mut new_faces = Vec::new();
        for &f in &self.vf[a] {
            if faces.contains(&f) {
                continue;
            }
            let old = self.corners_of(f);
            let t = self.t[f].map(|v| if v == a { b } else { v });
            let new = t.map(|v| self.p[v]);
            let (n0, n1) = (normal(old), normal(new));
            if n1.length() <= 1e-12 * old[0].length().max(1.0)
                || unit(n0).dot(unit(n1)) < 0.5
                || shape(new) < MIN_SHAPE.min(shape(old))
            {
                return false;
            }
            for v in t {
                if v != b && (self.p[v] - target).length() > hi * self.edge_target(v, b) {
                    return false;
                }
            }
            new_faces.push((f, t));
        }
        // features: the other feature edge of a moves onto b
        let moved: Vec<usize> = na
            .iter()
            .copied()
            .filter(|&v| v != b && self.feature.contains(&key(a, v)))
            .collect();
        for v in &moved {
            self.feature.remove(&key(a, *v));
            self.feature.insert(key(b, *v));
        }
        self.feature.remove(&key(a, b));
        for f in faces {
            self.remove_face(f);
        }
        for (f, t) in new_faces {
            self.remove_face(f);
            self.add_face(t);
        }
        self.stats.collapses += 1;
        true
    }

    fn valence_target(&self, v: usize) -> i64 {
        if self.vf[v].iter().any(|&f| self.keep[f]) {
            4
        } else {
            6
        }
    }

    fn flip(&mut self, a: usize, b: usize) -> bool {
        if self.feature.contains(&key(a, b)) {
            return false;
        }
        let faces = self.edge_faces(a, b);
        if faces.len() != 2 || self.touches_keep(&faces) {
            return false;
        }
        // orient: f0 = [a, b, c], f1 = [b, a, d]
        let (mut f0, mut f1) = (faces[0], faces[1]);
        let has =
            |t: [usize; 3], x: usize, y: usize| (0..3).any(|i| t[i] == x && t[(i + 1) % 3] == y);
        if !has(self.t[f0], a, b) {
            std::mem::swap(&mut f0, &mut f1);
        }
        if !has(self.t[f0], a, b) || !has(self.t[f1], b, a) {
            return false;
        }
        let other = |t: [usize; 3]| t.into_iter().find(|&v| v != a && v != b);
        let (Some(c), Some(d)) = (other(self.t[f0]), other(self.t[f1])) else {
            return false;
        };
        if c == d || self.neighbours(c).contains(&d) {
            return false;
        }
        let val = |v: usize| self.neighbours(v).len() as i64;
        let dev = |v: usize, delta: i64| (val(v) + delta - self.valence_target(v)).abs();
        let before = dev(a, 0) + dev(b, 0) + dev(c, 0) + dev(d, 0);
        let after = dev(a, -1) + dev(b, -1) + dev(c, 1) + dev(d, 1);
        if after >= before {
            return false;
        }
        let n0 = unit(normal(self.corners_of(f0)));
        let n1 = unit(normal(self.corners_of(f1)));
        if n0.dot(n1) < (20f64).to_radians().cos() {
            return false;
        }
        let (t0, t1) = ([c, a, d], [d, b, c]);
        let worst_old = shape(self.corners_of(f0)).min(shape(self.corners_of(f1)));
        let worst_new = shape(t0.map(|v| self.p[v])).min(shape(t1.map(|v| self.p[v])));
        if worst_new < MIN_SHAPE.min(worst_old) {
            return false;
        }
        let m = unit(n0 + n1);
        for t in [t0, t1] {
            let n = normal(t.map(|v| self.p[v]));
            if n.length() <= 0.0 || unit(n).dot(m) < 0.8 {
                return false;
            }
        }
        self.remove_face(f0);
        self.remove_face(f1);
        self.add_face(t0);
        self.add_face(t1);
        self.stats.flips += 1;
        true
    }

    fn smooth(&mut self) {
        let n = self.p.len();
        let mut next = self.p.clone();
        for v in 0..n {
            if self.vf[v].is_empty() || self.fixed(v) || self.vf[v].iter().any(|&f| self.keep[f]) {
                continue;
            }
            let x = self.p[v];
            let lines: Vec<usize> = self
                .neighbours(v)
                .into_iter()
                .filter(|&u| self.feature.contains(&key(u, v)))
                .collect();
            let moved = if lines.len() == 2 {
                let c = (self.p[lines[0]] + self.p[lines[1]]) / 2.0;
                let dir = unit(self.p[lines[1]] - self.p[lines[0]]);
                x + dir * (c - x).dot(dir) * 0.5
            } else if lines.is_empty() {
                let (mut c, mut w, mut nn) = (Vec3::ZERO, 0.0, Vec3::ZERO);
                for &f in &self.vf[v] {
                    let t = self.corners_of(f);
                    let cross = normal(t);
                    let area = cross.length();
                    c = c + (t[0] + t[1] + t[2]) / 3.0 * area;
                    w += area;
                    nn = nn + cross;
                }
                if w <= 0.0 {
                    continue;
                }
                let nn = unit(nn);
                let d = c / w - x;
                x + (d - nn * d.dot(nn)) * 0.5
            } else {
                continue;
            };
            next[v] = self.project(v, moved);
        }
        for v in 0..n {
            if next[v] == self.p[v] {
                continue;
            }
            // keep every incident facet's orientation
            let ok = self.vf[v].iter().all(|&f| {
                let old = self.corners_of(f);
                let new = self.t[f].map(|u| if u == v { next[v] } else { self.p[u] });
                let (a, b) = (normal(old), normal(new));
                b.length() > 0.0
                    && unit(a).dot(unit(b)) > 0.7
                    && shape(new) >= MIN_SHAPE.min(shape(old))
            });
            if ok {
                self.p[v] = next[v];
                self.refresh(v);
            }
        }
    }

    fn facet_deviation(&self, f: usize) -> f64 {
        samples(self.corners_of(f))
            .into_iter()
            .map(|x| self.surface.closest(x).map_or(0.0, |(_, d)| d))
            .fold(0.0, f64::max)
    }
}

/// Remesh the facets of a closed surface that are not marked `keep`.
pub fn remesh(
    vertices: &[Vec3],
    triangles: &[[usize; 3]],
    keep: &[bool],
    opts: RemeshOptions,
) -> Result<Remeshed, String> {
    if keep.len() != triangles.len() || triangles.iter().flatten().any(|&v| v >= vertices.len()) {
        return Err("invalid remesh input".into());
    }
    if !(opts.tolerance > 0.0
        && opts.min_size > 0.0
        && opts.max_size >= opts.min_size
        && opts.rate > 0.0
        && opts.tolerance.is_finite()
        && opts.max_size.is_finite()
        && opts.rate.is_finite()
        && opts.feature_angle > 0.0
        && opts.feature_angle < std::f64::consts::PI)
    {
        return Err("invalid remesh options".into());
    }
    let geometry = SurfaceGeometry::from_triangles(vertices.to_vec(), triangles.to_vec());
    // curvature sizes on wall facets
    let kappa = normal_cycle_curvature(&geometry, opts.feature_angle);
    let mut sizes: Vec<Weighted> = Vec::new();
    for (f, &k) in kappa.iter().enumerate() {
        if keep[f] || !k.is_finite() {
            continue;
        }
        let h = chord_length(k, opts.tolerance).clamp(opts.min_size, opts.max_size);
        if h < opts.max_size {
            sizes.push(Weighted {
                primitive: Primitive::Triangle(geometry.corners(f)),
                size: h,
                reach: 0.0,
            });
        }
    }
    // feature edges of the source and the rims of the kept facets
    let mut feature = HashSet::new();
    let mut line_items = Vec::new();
    for e in &geometry.edges {
        let [f, g] = e.facets;
        let [a, b] = e.vertices;
        if keep[f] != keep[g] {
            if let Some(h) = opts.rim_size {
                sizes.push(Weighted {
                    primitive: Primitive::Segment([vertices[a], vertices[b]]),
                    size: h.clamp(opts.min_size, opts.max_size),
                    reach: 0.0,
                });
            }
        } else if !keep[f] && e.angle > opts.feature_angle {
            feature.insert(key(a, b));
            line_items.push(Primitive::Segment([vertices[a], vertices[b]]));
        }
    }
    let wall: Vec<Primitive> = triangles
        .iter()
        .zip(keep)
        .filter(|(_, k)| !**k)
        .map(|(t, _)| Primitive::Triangle(t.map(|v| vertices[v])))
        .collect();
    if wall.is_empty() {
        return Err("no wall facets to remesh".into());
    }
    let surface = ClosestSet::new(wall);
    let lines = ClosestSet::new(line_items);
    let size = WeightedSet::new(sizes);
    let mut vf = vec![Vec::new(); vertices.len()];
    for (f, t) in triangles.iter().enumerate() {
        for &v in t {
            vf[v].push(f);
        }
    }
    let mut locked = vec![false; vertices.len()];
    for (f, t) in triangles.iter().enumerate() {
        if keep[f] {
            for &v in t {
                locked[v] = true;
            }
        }
    }
    let mut w = Work {
        p: vertices.to_vec(),
        t: triangles.to_vec(),
        alive: vec![true; triangles.len()],
        keep: keep.to_vec(),
        vf,
        corner: vec![false; vertices.len()],
        h: Vec::new(),
        locked,
        feature,
        surface: &surface,
        lines: &lines,
        size: &size,
        opts,
        stats: RemeshStats::default(),
    };
    // corners: feature valence other than two, or a sharp turn of the line
    for v in 0..vertices.len() {
        let lines: Vec<usize> = w
            .neighbours(v)
            .into_iter()
            .filter(|&u| w.feature.contains(&key(u, v)))
            .collect();
        let turn = lines.len() == 2 && {
            let (d0, d1) = (unit(w.p[lines[0]] - w.p[v]), unit(w.p[v] - w.p[lines[1]]));
            d0.dot(d1) < (45f64).to_radians().cos()
        };
        if (!lines.is_empty() && lines.len() != 2) || turn {
            w.corner[v] = true;
        }
    }
    for v in 0..vertices.len() {
        w.refresh(v);
    }
    w.stats.feature_edges = w.feature.len();
    w.stats.corners = w.corner.iter().filter(|&&c| c).count();
    let (hi, lo) = (4.0 / 3.0, 4.0 / 5.0);
    for _ in 0..opts.iterations {
        for (a, b) in w.edges() {
            if w.vf[a].is_empty() || w.vf[b].is_empty() || w.edge_faces(a, b).len() != 2 {
                continue;
            }
            let l = (w.p[a] - w.p[b]).length();
            if l > hi * w.edge_target(a, b) {
                w.split(a, b);
            }
        }
        for (a, b) in w.edges() {
            if w.vf[a].is_empty() || w.vf[b].is_empty() || w.edge_faces(a, b).len() != 2 {
                continue;
            }
            let l = (w.p[a] - w.p[b]).length();
            if l < lo * w.edge_target(a, b) && !w.collapse(a, b, hi) {
                w.collapse(b, a, hi);
            }
        }
        for (a, b) in w.edges() {
            if w.edge_faces(a, b).len() == 2 {
                w.flip(a, b);
            }
        }
        w.smooth();
    }
    // flips can leave long diagonals after the last split pass
    for (a, b) in w.edges() {
        if w.vf[a].is_empty() || w.vf[b].is_empty() || w.edge_faces(a, b).len() != 2 {
            continue;
        }
        let l = (w.p[a] - w.p[b]).length();
        if l > hi * w.edge_target(a, b) {
            w.split(a, b);
        }
    }
    // sliver cleanup: collapse edges far below the target that survived
    // (needles left by the source triangulation) and flip towards valence
    for (a, b) in w.edges() {
        if w.vf[a].is_empty() || w.vf[b].is_empty() || w.edge_faces(a, b).len() != 2 {
            continue;
        }
        let l = (w.p[a] - w.p[b]).length();
        if l < 0.25 * w.edge_target(a, b) && !w.collapse(a, b, hi) {
            w.collapse(b, a, hi);
        }
    }
    // deviation pass: split the longest edge of each facet still too far
    for _ in 0..8 {
        let bad: Vec<(usize, usize)> = (0..w.t.len())
            .filter(|&f| w.alive[f] && !w.keep[f] && w.facet_deviation(f) > opts.tolerance)
            .map(|f| {
                let t = w.t[f];
                let k = (0..3)
                    .max_by(|&i, &j| {
                        let li = (w.p[t[i]] - w.p[t[(i + 1) % 3]]).length();
                        let lj = (w.p[t[j]] - w.p[t[(j + 1) % 3]]).length();
                        li.total_cmp(&lj)
                    })
                    .unwrap_or(0);
                key(t[k], t[(k + 1) % 3])
            })
            .collect::<HashSet<_>>()
            .into_iter()
            .collect();
        if bad.is_empty() {
            break;
        }
        for (a, b) in bad {
            if w.split(a, b) {
                w.stats.deviation_splits += 1;
            }
        }
    }
    // compact
    let mut map = vec![usize::MAX; w.p.len()];
    let mut out_v = Vec::new();
    let mut out_t = Vec::new();
    let mut source = Vec::new();
    let mut origin: HashMap<[usize; 3], usize> = HashMap::new();
    for (f, t) in triangles.iter().enumerate() {
        if keep[f] {
            origin.insert(*t, f);
        }
    }
    for (f, t) in w.t.iter().enumerate() {
        if !w.alive[f] {
            continue;
        }
        let tt = t.map(|v| {
            if map[v] == usize::MAX {
                map[v] = out_v.len();
                out_v.push(w.p[v]);
            }
            map[v]
        });
        out_t.push(tt);
        source.push(if w.keep[f] {
            origin.get(t).copied()
        } else {
            None
        });
    }
    let mut stats = w.stats;
    stats.max_deviation = (0..w.t.len())
        .filter(|&f| w.alive[f] && !w.keep[f])
        .map(|f| w.facet_deviation(f))
        .fold(0.0, f64::max);
    let new_wall = ClosestSet::new(
        out_t
            .iter()
            .zip(&source)
            .filter(|(_, s)| s.is_none())
            .map(|(t, _)| Primitive::Triangle(t.map(|v| out_v[v])))
            .collect(),
    );
    stats.max_reverse_deviation = triangles
        .iter()
        .zip(keep)
        .filter(|(_, k)| !**k)
        .flat_map(|(t, _)| samples(t.map(|v| vertices[v])))
        .map(|x| new_wall.closest(x).map_or(0.0, |(_, d)| d))
        .fold(0.0, f64::max);
    let (mut min_e, mut max_e) = (f64::INFINITY, 0.0f64);
    for t in out_t
        .iter()
        .zip(&source)
        .filter(|(_, s)| s.is_none())
        .map(|(t, _)| t)
    {
        for k in 0..3 {
            let l = (out_v[t[k]] - out_v[t[(k + 1) % 3]]).length();
            min_e = min_e.min(l);
            max_e = max_e.max(l);
        }
    }
    stats.min_edge = min_e;
    stats.max_edge = max_e;
    Ok(Remeshed {
        vertices: out_v,
        triangles: out_t,
        source,
        stats,
    })
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::sizing::shapes;

    /// Each edge used once in each direction: a closed, consistently
    /// oriented manifold.
    pub(crate) fn assert_closed(t: &[[usize; 3]]) {
        let mut half = HashSet::new();
        for f in t {
            for k in 0..3 {
                assert!(
                    half.insert((f[k], f[(k + 1) % 3])),
                    "edge used twice in one direction"
                );
            }
        }
        for &(a, b) in &half {
            assert!(half.contains(&(b, a)), "open edge");
        }
    }

    pub(crate) fn volume(v: &[Vec3], t: &[[usize; 3]]) -> f64 {
        t.iter()
            .map(|f| v[f[0]].dot(v[f[1]].cross(v[f[2]])) / 6.0)
            .sum()
    }

    pub(crate) fn sphere(tolerance: f64) -> Remeshed {
        let (v, t) = shapes::icosphere(Vec3::ZERO, 10.0, 5);
        let keep = vec![false; t.len()];
        let opts = RemeshOptions {
            tolerance,
            min_size: 0.2,
            max_size: 4.0,
            iterations: 8,
            ..Default::default()
        };
        remesh(&v, &t, &keep, opts).unwrap()
    }

    // Keep rule (SW-18): the remeshed sphere is closed and every sampled
    // point of every facet is within the tolerance of the source, both ways,
    // with edges near the chord length sqrt(6 eps r) of the target. A
    // remesher that drops the projection, the deviation pass or the
    // curvature size breaks the bound or the edge range; no other test
    // measures the wall against its source.
    #[test]
    fn sphere_meets_the_chord_tolerance() {
        let eps = 0.02;
        let r = sphere(eps);
        assert_closed(&r.triangles);
        assert!(r.stats.max_deviation <= eps, "{:?}", r.stats);
        assert!(r.stats.max_reverse_deviation <= 1.5 * eps, "{:?}", r.stats);
        let chord = chord_length(0.1, eps);
        let mut lengths: Vec<f64> = r
            .triangles
            .iter()
            .flat_map(|t| (0..3).map(|k| (r.vertices[t[k]] - r.vertices[t[(k + 1) % 3]]).length()))
            .collect();
        lengths.sort_by(f64::total_cmp);
        let median = lengths[lengths.len() / 2];
        assert!(
            (median / chord - 1.0).abs() < 0.3,
            "median {median}, chord {chord}"
        );
        // fewer facets than the 20 480 of the source
        assert!(r.triangles.len() < 10_000, "{}", r.triangles.len());
    }

    // Keep rule (SW-18): on a box with one kept side, the twelve edges
    // survive as feature lines (the volume stays exact and every vertex
    // stays on a side), the eight corners stay, and the kept side keeps its
    // two triangles. Moving a feature vertex off its line, flipping a
    // feature edge or touching a kept facet changes the volume or the kept
    // triangles.
    #[test]
    fn box_keeps_features_and_kept_facets() {
        let (v, t) = shapes::cuboid(Vec3::new(3.0, 2.0, 1.0));
        let mut keep = vec![false; t.len()];
        keep[0] = true;
        keep[1] = true;
        let opts = RemeshOptions {
            tolerance: 0.01,
            min_size: 0.2,
            max_size: 0.4,
            rim_size: Some(0.4),
            ..Default::default()
        };
        let r = remesh(&v, &t, &keep, opts).unwrap();
        assert_closed(&r.triangles);
        assert!((volume(&r.vertices, &r.triangles) - 6.0).abs() < 1e-9);
        let kept: Vec<usize> = r.source.iter().filter_map(|s| *s).collect();
        assert_eq!(kept.len(), 2);
        for corner in &v {
            assert!(r.vertices.iter().any(|p| (*p - *corner).length() < 1e-12));
        }
        for p in &r.vertices {
            let on = [p.x, p.x - 3.0, p.y, p.y - 2.0, p.z, p.z - 1.0]
                .iter()
                .any(|d| d.abs() < 1e-9);
            assert!(on, "{p:?} off the box");
        }
        assert_eq!(r.stats.corners, 8);
        assert!(r.triangles.len() > 100);
    }
}
