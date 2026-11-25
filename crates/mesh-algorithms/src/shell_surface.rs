//! Zero-thickness shells inside the fluid as surfaces of the tetrahedral
//! scaffold (SHL-01 to SHL-03). The maintainer's choice of 2026-11-25: the
//! shell goes into the scaffold before the dual is built, as baffles are
//! meshed in STAR-CCM+ and Fluent, so no finished cell is ever cut.
//!
//! 1. [`imprint`]: the wall STL is cut along the shell plane inside the
//!    shell polygon (a planar section of the triangles, Kobbelt's group's
//!    plane-based view of mesh Booleans, Trettner, Nehring-Wirxel & Kobbelt
//!    2022). The cut is a labelled feature line of the wall remesh
//!    (`remesh::remesh_lines`), so the remeshed wall has edges exactly where
//!    the shell meets it, and a baffle closes the duct without a gap.
//! 2. [`mesh_shell`]: the part of the polygon inside the fluid is meshed in
//!    its plane: a constrained Delaunay triangulation (`planar`) of the
//!    wall line and of the free edges resampled to the local spacing, the
//!    regions inside the fluid kept by their winding number (Jacobson,
//!    Kavan & Sorkine-Hornung 2013), then Delaunay refinement by
//!    circumcentres clear of the boundary (Chew 1989; Ruppert 1995;
//!    Shewchuk 2002). A boundary edge a refinement point would crowd is
//!    split instead (Ruppert's segment split): free-edge pieces as the size
//!    asks, wall-line edges once (the caller bisects the wall triangles
//!    there, see `ShellSurface::wall_splits`). The size is the shell's own
//!    (`size_items`, in the build's size field) and is graded from the
//!    wall-line edges, which keep every kink of a polygonal wall. No interior edge
//!    joins two free-edge vertices, so the two sides of the shell share
//!    nothing but the free edges once the scaffold is cracked.

#![allow(clippy::needless_range_loop)] // parallel per-vertex arrays

use crate::delaunay::{HashMap, HashSet};
use crate::duct::shell::Frame;
use crate::planar::{self, Cdt, P2};
use crate::remesh::ClosestSet;
use crate::sizing::SizeField;
use crate::sizing::surface::{Primitive, Weighted, WeightedSet};
use mesh_core::geometry::Vec3;
use std::collections::VecDeque;

/// A wall cut along the shell planes (see [`imprint`]).
#[derive(Debug)]
pub struct Imprint {
    pub vertices: Vec<Vec3>,
    pub triangles: Vec<[usize; 3]>,
    pub keep: Vec<bool>,
    pub group: Vec<usize>,
    /// Edges on a shell's line on the wall, with the shell's index.
    pub lines: Vec<((usize, usize), usize)>,
}

fn key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

/// A triangle soup with directed-edge lookup, for local splits.
struct Soup {
    v: Vec<Vec3>,
    t: Vec<[usize; 3]>,
    parent: Vec<usize>,
    alive: Vec<bool>,
    edge: HashMap<(usize, usize), usize>,
}

impl Soup {
    fn add(&mut self, t: [usize; 3], parent: usize) {
        let f = self.t.len();
        for k in 0..3 {
            self.edge.insert((t[k], t[(k + 1) % 3]), f);
        }
        self.t.push(t);
        self.parent.push(parent);
        self.alive.push(true);
    }

    fn remove(&mut self, f: usize) {
        self.alive[f] = false;
        let t = self.t[f];
        for k in 0..3 {
            self.edge.remove(&(t[k], t[(k + 1) % 3]));
        }
    }

    /// Split both triangles at edge (a, b) at the new vertex m.
    fn split_edge(&mut self, a: usize, b: usize, m: usize) {
        for (x, y) in [(a, b), (b, a)] {
            if let Some(&f) = self.edge.get(&(x, y)) {
                let t = self.t[f];
                let c = t.into_iter().find(|&c| c != x && c != y).unwrap_or(x);
                let parent = self.parent[f];
                self.remove(f);
                self.add([x, m, c], parent);
                self.add([m, y, c], parent);
            }
        }
    }

    fn split_face(&mut self, f: usize, m: usize) {
        let [a, b, c] = self.t[f];
        let parent = self.parent[f];
        self.remove(f);
        self.add([a, b, m], parent);
        self.add([b, c, m], parent);
        self.add([c, a, m], parent);
    }
}

/// Barycentric coordinates of x in the triangle (projected on its plane).
fn barycentric(p: [Vec3; 3], x: Vec3) -> [f64; 3] {
    let n = (p[1] - p[0]).cross(p[2] - p[0]);
    let nn = n.dot(n);
    if nn <= 0.0 {
        return [-1.0; 3];
    }
    let l = |a: Vec3, b: Vec3| (b - a).cross(x - a).dot(n) / nn;
    [l(p[1], p[2]), l(p[2], p[0]), l(p[0], p[1])]
}

/// Cut the wall along the plane of every shell polygon, inside the polygon.
/// Facets marked `keep` (port caps) may not be cut. Returns the refined
/// wall with the cut edges as lines labelled by shell.
pub fn imprint(
    vertices: &[Vec3],
    triangles: &[[usize; 3]],
    keep: &[bool],
    group: &[usize],
    shells: &[(&str, &[Vec3])],
) -> Result<Imprint, String> {
    if keep.len() != triangles.len() || group.len() != triangles.len() || vertices.is_empty() {
        return Err("invalid imprint input".into());
    }
    let frames = shells
        .iter()
        .map(|(name, p)| Frame::new(p).map_err(|e| format!("shell {name:?}: {e}")))
        .collect::<Result<Vec<_>, _>>()?;
    for i in 0..frames.len() {
        for j in i + 1..frames.len() {
            if polygons_meet(&frames[i], shells[i].1, &frames[j], shells[j].1) {
                return Err(format!(
                    "shells {:?} and {:?} meet; shells that touch each other are not supported yet",
                    shells[i].0, shells[j].0
                ));
            }
        }
    }
    let mut lo = vertices[0];
    let mut hi = lo;
    for &v in vertices {
        lo = Vec3::new(lo.x.min(v.x), lo.y.min(v.y), lo.z.min(v.z));
        hi = Vec3::new(hi.x.max(v.x), hi.y.max(v.y), hi.z.max(v.z));
    }
    let extent = (hi - lo).max_abs_component().max(f64::MIN_POSITIVE);
    let eps = 1e-10 * extent;
    let mut soup = Soup {
        v: vertices.to_vec(),
        t: Vec::new(),
        parent: Vec::new(),
        alive: Vec::new(),
        edge: HashMap::default(),
    };
    for (f, &t) in triangles.iter().enumerate() {
        soup.add(t, f);
    }
    let mut lines: HashMap<(usize, usize), usize> = HashMap::default();
    for (k, frame) in frames.iter().enumerate() {
        let name = shells[k].0;
        let mut s: Vec<f64> = soup
            .v
            .iter()
            .map(|&p| {
                let d = frame.distance(p);
                if d.abs() <= eps { 0.0 } else { d }
            })
            .collect();
        let inside = |p: Vec3| frame.inside(frame.uv(p)) >= -eps;
        // 1. points where the polygon outline leaves the wall section
        let mut ends: Vec<(usize, Vec3)> = Vec::new();
        for f in 0..soup.t.len() {
            if !soup.alive[f] {
                continue;
            }
            let t = soup.t[f];
            if t.iter().all(|&v| s[v] == 0.0) {
                return Err(format!("shell {name:?} lies in the plane of a wall facet"));
            }
            let Some((p, q)) = section(t.map(|v| (soup.v[v], s[v]))) else {
                continue;
            };
            let (a, b) = (frame.uv(p), frame.uv(q));
            let (mut t0, mut t1) = (0.0_f64, 1.0_f64);
            let mut empty = false;
            for i in 0..frame.poly.len() {
                let (sa, sb) = (frame.edge_side(i, a), frame.edge_side(i, b));
                if sa < 0.0 && sb < 0.0 {
                    empty = true;
                    break;
                }
                if sa < 0.0 {
                    t0 = t0.max(sa / (sa - sb));
                } else if sb < 0.0 {
                    t1 = t1.min(sa / (sa - sb));
                }
            }
            let length = (q - p).length();
            if empty || (t1 - t0) * length <= eps {
                continue;
            }
            if t0 * length > eps {
                ends.push((f, p + (q - p) * t0));
            }
            if (1.0 - t1) * length > eps {
                ends.push((f, p + (q - p) * t1));
            }
        }
        for (f0, x) in ends {
            // the live triangle that holds x (f0 or a piece of it)
            let found = (0..soup.t.len())
                .filter(|&f| soup.alive[f] && (f == f0 || soup.parent[f] == soup.parent[f0]))
                .map(|f| (f, barycentric(soup.t[f].map(|v| soup.v[v]), x)))
                .filter(|(_, l)| l.iter().all(|&c| c >= -1e-9))
                .max_by(|a, b| {
                    let m = |l: &[f64; 3]| l.iter().copied().fold(f64::INFINITY, f64::min);
                    m(&a.1).total_cmp(&m(&b.1))
                });
            let Some((f, l)) = found else {
                return Err(format!(
                    "shell {name:?}: free edge end not found on the wall"
                ));
            };
            if keep[soup.parent[f]] {
                return Err(cap_error(name));
            }
            let t = soup.t[f];
            let small: Vec<usize> = (0..3).filter(|&i| l[i] < 1e-9).collect();
            match small.len() {
                0 => {
                    let m = soup.v.len();
                    soup.v.push(x);
                    s.push(0.0);
                    soup.split_face(f, m);
                }
                1 => {
                    let i = small[0];
                    let (a, b) = (t[(i + 1) % 3], t[(i + 2) % 3]);
                    if lines.contains_key(&key(a, b)) {
                        return Err(format!("shell {name:?} meets the line of another shell"));
                    }
                    let m = soup.v.len();
                    soup.v.push(x);
                    s.push(0.0);
                    soup.split_edge(a, b, m);
                }
                _ => {
                    // at a vertex: that vertex is the end
                    let i = (0..3).find(|i| !small.contains(i)).unwrap_or(0);
                    s[t[i]] = 0.0;
                }
            }
        }
        // 2. every edge across the plane inside the polygon gets its point
        let mut crossing: Vec<(usize, usize)> = Vec::new();
        for f in 0..soup.t.len() {
            if !soup.alive[f] {
                continue;
            }
            let t = soup.t[f];
            for i in 0..3 {
                let (a, b) = (t[i], t[(i + 1) % 3]);
                if a < b && s[a] * s[b] < 0.0 {
                    crossing.push((a, b));
                }
            }
        }
        crossing.sort_unstable();
        crossing.dedup();
        for (a, b) in crossing {
            let x = soup.v[a] + (soup.v[b] - soup.v[a]) * (s[a] / (s[a] - s[b]));
            if !inside(x) {
                continue;
            }
            for e in [(a, b), (b, a)] {
                if let Some(&f) = soup.edge.get(&e)
                    && keep[soup.parent[f]]
                {
                    return Err(cap_error(name));
                }
            }
            if lines.contains_key(&key(a, b)) {
                return Err(format!("shell {name:?} meets the line of another shell"));
            }
            let m = soup.v.len();
            soup.v.push(x);
            s.push(0.0);
            soup.split_edge(a, b, m);
        }
        // 3. the cut: edges on the plane inside the polygon
        for f in 0..soup.t.len() {
            if !soup.alive[f] {
                continue;
            }
            let t = soup.t[f];
            for i in 0..3 {
                let (a, b) = (t[i], t[(i + 1) % 3]);
                if s[a] == 0.0 && s[b] == 0.0 && inside((soup.v[a] + soup.v[b]) / 2.0) {
                    if keep[soup.parent[f]] {
                        return Err(cap_error(name));
                    }
                    if let Some(&other) = lines.get(&key(a, b))
                        && other != k
                    {
                        return Err(format!("shell {name:?} meets the line of another shell"));
                    }
                    lines.insert(key(a, b), k);
                }
            }
        }
    }
    let mut out = Imprint {
        vertices: soup.v,
        triangles: Vec::new(),
        keep: Vec::new(),
        group: Vec::new(),
        lines: lines.into_iter().collect(),
    };
    out.lines.sort_unstable();
    for f in 0..soup.t.len() {
        if soup.alive[f] {
            out.triangles.push(soup.t[f]);
            out.keep.push(keep[soup.parent[f]]);
            out.group.push(group[soup.parent[f]]);
        }
    }
    Ok(out)
}

fn cap_error(name: &str) -> String {
    format!(
        "shell {name:?} meets an inlet or outlet cap; shells may meet walls only for now (SHL-03 follow-up)"
    )
}

/// The segment where a triangle (corners with their signed distances)
/// meets the plane, when it crosses it.
fn section(c: [(Vec3, f64); 3]) -> Option<(Vec3, Vec3)> {
    let mut points = Vec::new();
    for i in 0..3 {
        let (p, sp) = c[i];
        let (q, sq) = c[(i + 1) % 3];
        if sp == 0.0 {
            points.push(p);
        }
        if sp * sq < 0.0 {
            points.push(p + (q - p) * (sp / (sp - sq)));
        }
    }
    match points[..] {
        [p, q] => Some((p, q)),
        _ => None,
    }
}

/// Whether two convex planar polygons touch or cross.
fn polygons_meet(fa: &Frame, a: &[Vec3], fb: &Frame, b: &[Vec3]) -> bool {
    let tol = 1e-9 * fa.scale.max(fb.scale);
    let pierces = |f: &Frame, poly: &[Vec3]| {
        (0..poly.len()).any(|i| {
            let (p, q) = (poly[i], poly[(i + 1) % poly.len()]);
            let (sp, sq) = (f.distance(p), f.distance(q));
            if sp.abs() <= tol && sq.abs() <= tol {
                // in the plane: either end inside
                return f.inside(f.uv(p)) >= -tol || f.inside(f.uv(q)) >= -tol;
            }
            if sp * sq > 0.0 && sp.abs() > tol && sq.abs() > tol {
                return false;
            }
            let x = if (sp - sq).abs() > 0.0 {
                p + (q - p) * (sp / (sp - sq))
            } else {
                p
            };
            f.inside(f.uv(x)) >= -tol
        })
    };
    pierces(fb, a) || pierces(fa, b)
}

/// The spacing a shell is meshed to: the scaffold's interior spacing, the
/// size of the wall triangles and of the edges of the shells' lines on the
/// wall grown by `rate` with the distance from them, bounded by `max_size`
/// and the size field (which holds each shell's own sizes, see
/// [`size_items`]). The line edges count on their own because the wall
/// line keeps every kink of a polygonal wall: a line edge far shorter
/// than the wall triangles beside it (where the line ends next to a kink)
/// asks for short free-edge pieces there, or its triangle on the shell is
/// a needle.
pub struct Spacing<'a> {
    set: WeightedSet,
    rate: f64,
    max_size: f64,
    field: Option<&'a SizeField>,
}

impl<'a> Spacing<'a> {
    pub fn new(
        vertices: &[Vec3],
        triangles: &[[usize; 3]],
        lines: &[(usize, usize)],
        rate: f64,
        max_size: f64,
        field: Option<&'a SizeField>,
    ) -> Self {
        let line_items = lines.iter().filter_map(|&(a, b)| {
            let l = (vertices[a] - vertices[b]).length();
            (l > 0.0).then(|| Weighted {
                primitive: Primitive::Segment([vertices[a], vertices[b]]),
                size: l.min(max_size),
                reach: 0.0,
            })
        });
        let items = triangles
            .iter()
            .map(|t| {
                let p = t.map(|v| vertices[v]);
                let l = ((p[0] - p[1]).length() + (p[1] - p[2]).length() + (p[2] - p[0]).length())
                    / 3.0;
                Weighted {
                    primitive: Primitive::Triangle(p),
                    size: l.min(max_size),
                    reach: 0.0,
                }
            })
            .chain(line_items)
            .collect();
        Self {
            set: WeightedSet::new(items),
            rate,
            max_size,
            field,
        }
    }

    pub fn at(&self, x: Vec3) -> f64 {
        let wall = self.set.envelope(x, self.rate, 0.0, self.max_size);
        self.field.map_or(wall, |f| wall.min(f.size(x)))
    }
}

/// Check a shell polygon: at least three finite corners, planar, convex,
/// counter-clockwise about its normal.
pub fn check_polygon(polygon: &[Vec3]) -> Result<(), String> {
    Frame::new(polygon).map(|_| ())
}

/// Whether two shell polygons touch or cross (shells that meet are not
/// supported yet).
pub fn polygons_touch(a: &[Vec3], b: &[Vec3]) -> Result<bool, String> {
    let (fa, fb) = (Frame::new(a)?, Frame::new(b)?);
    Ok(polygons_meet(&fa, a, &fb, b))
}

/// The narrow width of the part of a shell polygon inside the box
/// `lo`..`hi` (the smallest distance between two parallel lines that hold
/// it), or `None` when no part lies inside. The default shell size is an
/// eighth of it (the maintainer's choice of 2026-11-25: about eight cells
/// across a shell).
pub fn narrow_width(polygon: &[Vec3], lo: Vec3, hi: Vec3) -> Result<Option<f64>, String> {
    let frame = Frame::new(polygon)?;
    // Sutherland-Hodgman against the six faces of the box
    let c = |p: Vec3, axis: usize| [p.x, p.y, p.z][axis];
    let mut poly: Vec<Vec3> = polygon.to_vec();
    for axis in 0..3 {
        for (bound, sign) in [(c(lo, axis), 1.0), (c(hi, axis), -1.0)] {
            let side = |p: Vec3| sign * (c(p, axis) - bound);
            let mut out = Vec::with_capacity(poly.len() + 1);
            for i in 0..poly.len() {
                let (a, b) = (poly[i], poly[(i + 1) % poly.len()]);
                let (da, db) = (side(a), side(b));
                if da >= 0.0 {
                    out.push(a);
                }
                if (da >= 0.0) != (db >= 0.0) {
                    out.push(a + (b - a) * (da / (da - db)));
                }
            }
            poly = out;
            if poly.len() < 3 {
                return Ok(None);
            }
        }
    }
    let q: Vec<[f64; 2]> = poly.iter().map(|&p| frame.uv(p)).collect();
    let mut width = f64::INFINITY;
    for i in 0..q.len() {
        let (a, b) = (q[i], q[(i + 1) % q.len()]);
        let e = [b[0] - a[0], b[1] - a[1]];
        let l = e[0].hypot(e[1]);
        if l <= 1e-12 * frame.scale {
            continue;
        }
        let far = q
            .iter()
            .map(|p| ((p[0] - a[0]) * e[1] - (p[1] - a[1]) * e[0]).abs() / l)
            .fold(0.0, f64::max);
        width = width.min(far);
    }
    Ok((width.is_finite() && width > 0.0).then_some(width))
}

/// The size primitives a shell adds to the build's size field: its
/// polygon at `size` (a fan of triangles) and its free edges, the parts of
/// its outline inside the closed, outward wall, at `edge_size`. In the
/// field they refine the wall along the shell's line (so the line's edges
/// are split to the shell size), the scaffold around the shell, and the
/// shell itself.
pub fn size_items(
    polygon: &[Vec3],
    size: f64,
    edge_size: f64,
    wall_vertices: &[Vec3],
    wall_triangles: &[[usize; 3]],
) -> Result<Vec<Weighted>, String> {
    Frame::new(polygon)?;
    if !(size.is_finite() && size > 0.0 && edge_size.is_finite() && edge_size > 0.0) {
        return Err("shell sizes must be positive and finite".into());
    }
    let mut items: Vec<Weighted> = (1..polygon.len() - 1)
        .map(|i| Weighted {
            primitive: Primitive::Triangle([polygon[0], polygon[i], polygon[i + 1]]),
            size,
            reach: 0.0,
        })
        .collect();
    // the outline in pieces of at most four edge sizes (64 per side at
    // most); runs of pieces inside the fluid are the free edges
    for i in 0..polygon.len() {
        let (a, b) = (polygon[i], polygon[(i + 1) % polygon.len()]);
        let n = (((b - a).length() / (4.0 * edge_size)).ceil() as usize).clamp(1, 64);
        let point = |j: usize| a + (b - a) * (j as f64 / n as f64);
        let mut run: Option<usize> = None;
        for j in 0..=n {
            let inside = j < n
                && winding(
                    (point(j) + point(j + 1)) * 0.5,
                    wall_vertices,
                    wall_triangles,
                )
                .abs()
                    > 0.5;
            match (inside, run) {
                (true, None) => run = Some(j),
                (false, Some(start)) => {
                    items.push(Weighted {
                        primitive: Primitive::Segment([point(start), point(j)]),
                        size: edge_size,
                        reach: 0.0,
                    });
                    run = None;
                }
                _ => {}
            }
        }
    }
    Ok(items)
}

/// The generalized winding number of a closed, outward triangle mesh at p.
pub fn winding(p: Vec3, vertices: &[Vec3], triangles: &[[usize; 3]]) -> f64 {
    let mut sum = 0.0;
    for t in triangles {
        let [a, b, c] = t.map(|v| vertices[v] - p);
        let (la, lb, lc) = (a.length(), b.length(), c.length());
        let num = a.dot(b.cross(c));
        let den = la * lb * lc + a.dot(b) * lc + b.dot(c) * la + c.dot(a) * lb;
        sum += 2.0 * num.atan2(den);
    }
    sum / (4.0 * std::f64::consts::PI)
}

/// A meshed shell.
pub struct ShellSurface {
    /// New vertices; their ids start at the `base` given to [`mesh_shell`].
    pub vertices: Vec<Vec3>,
    /// Counter-clockwise about the shell normal; ids below `base` are wall
    /// vertices (on the line where the shell meets the wall).
    pub triangles: Vec<[usize; 3]>,
    /// Edges of the wall line split for the shell, as given in `line`,
    /// with the new vertices on each from its first vertex to its second.
    /// The wall triangles at these edges must be split at them too.
    pub wall_splits: Vec<((usize, usize), Vec<usize>)>,
    pub stats: ShellSurfaceStats,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ShellSurfaceStats {
    pub triangles: usize,
    pub area: f64,
    /// Length of the shell's edges inside the fluid.
    pub free_edge_length: f64,
    /// Length of the line where the shell meets the wall.
    pub wall_line_length: f64,
    /// Ends of the wall line (where a free edge meets the wall).
    pub free_edge_ends: usize,
    pub min_angle_degrees: f64,
    /// Pieces of the free edges inside the fluid.
    pub free_edge_pieces: usize,
    /// Shortest, median and longest edge of the shell triangles.
    pub edge_lengths: [f64; 3],
    /// Wall-line edges split (each once at most) and free-edge pieces
    /// split by the refinement.
    pub wall_line_splits: usize,
    pub free_edge_splits: usize,
}

/// A boundary segment of the shell's triangulation.
#[derive(Clone, Copy)]
enum Segment {
    /// A piece of a free edge.
    Free,
    /// A piece of the wall-line edge `edge`, between its parameters `t`
    /// (0 at the edge's first vertex).
    Wall { edge: (usize, usize), t: [f64; 2] },
}

#[derive(Clone, Copy, PartialEq)]
enum Node {
    Wall(usize),
    New,
}

/// The points of a shell's triangulation before refinement.
#[derive(Default)]
struct Points {
    nodes: Vec<Node>,
    uv: Vec<P2>,
    position: Vec<Vec3>,
    /// On the polygon outline (a free edge or one of its ends).
    free: Vec<bool>,
    of_wall: HashMap<usize, usize>,
}

impl Points {
    fn add(&mut self, frame: &Frame, node: Node, p: Vec3, on_outline: bool) -> usize {
        if let Node::Wall(w) = node
            && let Some(&i) = self.of_wall.get(&w)
        {
            self.free[i] |= on_outline;
            return i;
        }
        let i = self.nodes.len();
        self.nodes.push(node);
        self.uv.push(frame.uv(p));
        self.position.push(p);
        self.free.push(on_outline);
        if let Node::Wall(w) = node {
            self.of_wall.insert(w, i);
        }
        i
    }
}

/// Mesh the part of the shell `polygon` inside the fluid bounded by the
/// closed, outward `wall` (vertices and triangles, ports included). `line`
/// is the shell's line on the wall (wall vertex ids), `size` the spacing,
/// `base` the id of the first new vertex.
pub fn mesh_shell(
    polygon: &[Vec3],
    wall_vertices: &[Vec3],
    wall_triangles: &[[usize; 3]],
    line: &[(usize, usize)],
    size: &dyn Fn(Vec3) -> f64,
    base: usize,
) -> Result<ShellSurface, String> {
    let frame = Frame::new(polygon)?;
    let tol = 1e-7 * frame.scale;
    let mut stats = ShellSurfaceStats {
        min_angle_degrees: 180.0,
        ..Default::default()
    };
    // the line: degree at most two; its ends lie on the outline
    let mut degree: HashMap<usize, usize> = HashMap::default();
    for &(a, b) in line {
        if a.max(b) >= wall_vertices.len() {
            return Err("invalid shell line".into());
        }
        *degree.entry(a).or_default() += 1;
        *degree.entry(b).or_default() += 1;
        stats.wall_line_length += (wall_vertices[a] - wall_vertices[b]).length();
    }
    if degree.values().any(|&d| d > 2) {
        return Err("the shell's line on the wall branches".into());
    }
    let mut ends: Vec<(f64, usize)> = Vec::new();
    let mut line_vertices: Vec<usize> = degree.keys().copied().collect();
    line_vertices.sort_unstable();
    for &v in &line_vertices {
        if degree[&v] == 1 {
            let q = frame.uv(wall_vertices[v]);
            let pos = frame.perimeter(q, tol).ok_or(
                "the shell's line on the wall ends inside the polygon (the wall remesh lost a corner)",
            )?;
            ends.push((pos, v));
        }
    }
    stats.free_edge_ends = ends.len();
    // the outline: corners and line ends, in order
    let k = frame.poly.len();
    let mut outline: Vec<(f64, Option<usize>)> = (0..k).map(|i| (i as f64, None)).collect();
    for &(pos, v) in &ends {
        let at = outline
            .iter()
            .position(|&(p, w)| w.is_none() && corner_gap(&frame, p, pos) <= tol);
        match at {
            Some(i) => outline[i] = (pos, Some(v)),
            None => outline.push((pos, Some(v))),
        }
    }
    outline.sort_by(|a, b| a.0.total_cmp(&b.0));
    // CDT points: line vertices, outline breakpoints, resampled arcs
    let mut pts = Points::default();
    for &v in &line_vertices {
        pts.add(&frame, Node::Wall(v), wall_vertices[v], false);
    }
    let breaks: Vec<usize> = outline
        .iter()
        .map(|&(pos, w)| match w {
            Some(v) => pts.add(&frame, Node::Wall(v), wall_vertices[v], true),
            None => pts.add(
                &frame,
                Node::New,
                frame.point(frame.poly[pos as usize % k]),
                true,
            ),
        })
        .collect();
    let mut arcs: Vec<Vec<usize>> = Vec::new();
    for i in 0..breaks.len() {
        let (a, b) = (breaks[i], breaks[(i + 1) % breaks.len()]);
        let (pa, pb) = (pts.position[a], pts.position[b]);
        let length = (pb - pa).length();
        // ∫ ds / h by the trapezoid rule, on pieces a quarter of the
        // smallest spacing met on a first pass (the spacing can drop
        // steeply towards an end, next to a short wall-line edge)
        let at = |t: f64| size(pa + (pb - pa) * t).max(1e-12);
        let smallest = (0..=16)
            .map(|j| at(j as f64 / 16.0))
            .fold(f64::INFINITY, f64::min);
        let pieces = ((4.0 * length / smallest).ceil() as usize).clamp(16, 4096);
        let h: Vec<f64> = (0..=pieces).map(|j| at(j as f64 / pieces as f64)).collect();
        let mut cumulative = vec![0.0];
        for j in 0..pieces {
            let ds = length / pieces as f64;
            cumulative.push(cumulative[j] + 0.5 * ds * (1.0 / h[j] + 1.0 / h[j + 1]));
        }
        let total = cumulative[pieces];
        let n = total.round().max(1.0) as usize;
        let mut arc = vec![a];
        for m in 1..n {
            let target = total * m as f64 / n as f64;
            let j = cumulative[1..]
                .partition_point(|&c| c < target)
                .min(pieces - 1);
            let share = (target - cumulative[j]) / (cumulative[j + 1] - cumulative[j]).max(1e-300);
            let t = (j as f64 + share) / pieces as f64;
            arc.push(pts.add(&frame, Node::New, pa + (pb - pa) * t, true));
        }
        arc.push(b);
        arcs.push(arc);
    }
    let Points {
        mut nodes,
        uv,
        mut position,
        mut free,
        of_wall,
    } = pts;
    // triangulate
    let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    for q in &uv {
        for d in 0..2 {
            lo[d] = lo[d].min(q[d]);
            hi[d] = hi[d].max(q[d]);
        }
    }
    let mut cdt = Cdt::new(lo, hi);
    let mut id = Vec::with_capacity(uv.len());
    for (i, &q) in uv.iter().enumerate() {
        let c = cdt
            .insert(q)
            .map_err(|e| format!("shell point {i} not inserted: {e:?}"))?;
        id.push(c);
    }
    // back from CDT ids to node ids
    let mut node_of: HashMap<usize, usize> = (0..id.len()).map(|i| (id[i], i)).collect();
    for &(a, b) in line {
        cdt.constrain(id[of_wall[&a]], id[of_wall[&b]])?;
    }
    for arc in &arcs {
        for w in arc.windows(2) {
            cdt.constrain(id[w[0]], id[w[1]])?;
        }
    }
    cdt.make_delaunay();
    let regions = cdt.label_regions();
    // a region is fluid when the winding number of the wall at its largest
    // triangle's centroid says so
    let mut fluid = vec![false; regions];
    for r in 1..regions {
        let best = (0..cdt.tris.len())
            .filter(|&t| cdt.alive[t] && cdt.region[t] == r)
            .max_by(|&a, &b| area2(&cdt, a).total_cmp(&area2(&cdt, b)));
        if let Some(t) = best {
            let [a, b, c] = cdt.tris[t].map(|v| cdt.points[v]);
            let g = [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0];
            fluid[r] = winding(frame.point(g), wall_vertices, wall_triangles).abs() > 0.5;
        }
    }
    if !fluid.iter().any(|&f| f) {
        return Err("the shell polygon lies outside the fluid".into());
    }
    for arc in &arcs {
        for w in arc.windows(2) {
            let (a, b) = (id[w[0]], id[w[1]]);
            let touches = (0..cdt.tris.len()).any(|t| {
                cdt.alive[t]
                    && fluid[cdt.region[t]]
                    && cdt.tris[t].contains(&a)
                    && cdt.tris[t].contains(&b)
            });
            if touches {
                stats.free_edge_pieces += 1;
                stats.free_edge_length += (position[w[0]] - position[w[1]]).length();
            }
        }
    }
    // refine: circumcentres of large or poor triangles, clear of the boundary
    let mut segments = Vec::new();
    for &(a, b) in line {
        segments.push(Primitive::Segment([wall_vertices[a], wall_vertices[b]]));
    }
    for arc in &arcs {
        for w in arc.windows(2) {
            segments.push(Primitive::Segment([position[w[0]], position[w[1]]]));
        }
    }
    let boundary = ClosestSet::new(segments);
    // the same segments as CDT edges, for splits (a split leaves the
    // boundary's geometry, so `boundary` stays valid)
    let mut segs: Vec<(usize, usize, Segment)> = line
        .iter()
        .map(|&(a, b)| {
            (
                id[of_wall[&a]],
                id[of_wall[&b]],
                Segment::Wall {
                    edge: (a, b),
                    t: [0.0, 1.0],
                },
            )
        })
        .collect();
    for arc in &arcs {
        for w in arc.windows(2) {
            segs.push((id[w[0]], id[w[1]], Segment::Free));
        }
    }
    let mut wall_points: Vec<((usize, usize), f64, usize)> = Vec::new();
    const RADIUS_SHARE: f64 = 0.7;
    const RADIUS_EDGE: f64 = 1.5;
    const CLEARANCE: f64 = 0.45;
    const MAX_POINTS: usize = 2_000_000;
    const SPLIT_SHARE: f64 = 0.6;
    let mut queue: VecDeque<usize> = (0..cdt.tris.len())
        .filter(|&t| cdt.alive[t] && fluid[cdt.region[t]])
        .collect();
    let mut tried: HashSet<[usize; 3]> = HashSet::default();
    while let Some(t) = queue.pop_front() {
        if !cdt.alive[t] || !fluid[cdt.region[t]] || cdt.points.len() > MAX_POINTS {
            continue;
        }
        let tri = cdt.tris[t];
        let p = tri.map(|v| cdt.points[v]);
        let Some((c, r)) = planar::circumcircle(p[0], p[1], p[2]) else {
            continue;
        };
        let shortest = (0..3)
            .map(|i| planar::distance(p[i], p[(i + 1) % 3]))
            .fold(f64::INFINITY, f64::min);
        let g = [
            (p[0][0] + p[1][0] + p[2][0]) / 3.0,
            (p[0][1] + p[1][1] + p[2][1]) / 3.0,
        ];
        let h = size(frame.point(g));
        if r <= RADIUS_SHARE * h && r <= RADIUS_EDGE * shortest {
            continue;
        }
        let mut sorted_tri = tri;
        sorted_tri.sort_unstable();
        if !tried.insert(sorted_tri) {
            continue;
        }
        let Some(home) = cdt.locate(c) else { continue };
        if !fluid[cdt.region[home]] {
            continue;
        }
        let x = frame.point(c);
        let clear = CLEARANCE * size(x).min(h);
        if boundary.closest(x).is_none_or(|(_, d)| d < clear) {
            // Ruppert (1995): split the segment the centre crowds instead,
            // when it is long for the spacing at its middle; a wall-line
            // edge once at most (its wall triangles are bisected with it)
            let Some(i) = nearest_segment(&cdt, &segs, c) else {
                continue;
            };
            let (sa, sb, kind) = segs[i];
            let (pa, pb) = (cdt.points[sa], cdt.points[sb]);
            let mid = [(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0];
            let long = planar::distance(pa, pb) > SPLIT_SHARE * size(frame.point(mid));
            let allowed = match kind {
                Segment::Free => true,
                Segment::Wall { t, .. } => t[1] - t[0] > 0.75,
            };
            if !(long && allowed) {
                continue;
            }
            let Ok(v) = cdt.split_constrained(sa, sb) else {
                continue;
            };
            let (p3, first, second, on_free) = match kind {
                Segment::Free => {
                    stats.free_edge_splits += 1;
                    (frame.point(mid), Segment::Free, Segment::Free, true)
                }
                Segment::Wall { edge, t } => {
                    stats.wall_line_splits += 1;
                    let tm = 0.5 * (t[0] + t[1]);
                    let (wa, wb) = (wall_vertices[edge.0], wall_vertices[edge.1]);
                    wall_points.push((edge, tm, v));
                    (
                        wa + (wb - wa) * tm,
                        Segment::Wall {
                            edge,
                            t: [t[0], tm],
                        },
                        Segment::Wall {
                            edge,
                            t: [tm, t[1]],
                        },
                        false,
                    )
                }
            };
            node_of.insert(v, nodes.len());
            nodes.push(Node::New);
            position.push(p3);
            free.push(on_free);
            segs[i] = (sa, v, first);
            segs.push((v, sb, second));
            queue.extend(cdt.created.iter().copied());
            tried.remove(&sorted_tri);
            queue.push_back(t);
            continue;
        }
        if cdt.insert(c).is_ok() {
            queue.extend(cdt.created.iter().copied());
            queue.push_back(t);
        }
    }
    // no interior edge between two outline vertices: split it at its middle
    let is_free = |cdt: &Cdt, node_of: &HashMap<usize, usize>, v: usize| {
        node_of.get(&v).is_some_and(|&n| free[n]) && v < cdt.points.len()
    };
    for _ in 0..32 {
        let mut chords: Vec<(usize, usize)> = Vec::new();
        for t in 0..cdt.tris.len() {
            if !cdt.alive[t] || !fluid[cdt.region[t]] {
                continue;
            }
            for i in 0..3 {
                let [a, b] = cdt.edge(t, i);
                if a < b
                    && !cdt.is_constrained(a, b)
                    && is_free(&cdt, &node_of, a)
                    && is_free(&cdt, &node_of, b)
                {
                    chords.push((a, b));
                }
            }
        }
        if chords.is_empty() {
            break;
        }
        chords.sort_unstable();
        chords.dedup();
        for (a, b) in chords {
            if !cdt.has_edge(a, b) {
                continue;
            }
            let (pa, pb) = (cdt.points[a], cdt.points[b]);
            let m = [(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0];
            cdt.insert(m)
                .map_err(|e| format!("shell chord not split: {e:?}"))?;
        }
    }
    // a few smoothing sweeps of the interior points, kept only where every
    // triangle around stays positive, then Delaunay flips
    let fixed: HashSet<usize> = node_of.keys().copied().collect();
    for _ in 0..3 {
        for v in planar::CORNERS..cdt.points.len() {
            if fixed.contains(&v) {
                continue;
            }
            let around = cdt.around(v);
            if around.is_empty() || around.iter().any(|&t| !fluid[cdt.region[t]]) {
                continue;
            }
            let mut sum = [0.0, 0.0];
            let mut count = 0.0;
            for &t in &around {
                for w in cdt.tris[t] {
                    if w != v {
                        sum[0] += cdt.points[w][0];
                        sum[1] += cdt.points[w][1];
                        count += 1.0;
                    }
                }
            }
            let old = cdt.points[v];
            cdt.points[v] = [sum[0] / count, sum[1] / count];
            let folded = around.iter().any(|&t| {
                let [a, b, c] = cdt.tris[t].map(|w| cdt.points[w]);
                planar::orient(a, b, c) != std::cmp::Ordering::Greater
            });
            let near_boundary = boundary
                .closest(frame.point(cdt.points[v]))
                .is_none_or(|(_, d)| d < 0.3 * size(frame.point(cdt.points[v])));
            if folded || near_boundary {
                cdt.points[v] = old;
            }
        }
        cdt.make_delaunay();
    }
    // the triangles in the fluid regions
    let mut out_id: HashMap<usize, usize> = HashMap::default();
    let mut vertices = Vec::new();
    let mut triangles = Vec::new();
    for t in 0..cdt.tris.len() {
        if !cdt.alive[t] || !fluid[cdt.region[t]] {
            continue;
        }
        let tri = cdt.tris[t].map(|v| {
            if let Some(&n) = node_of.get(&v)
                && let Node::Wall(w) = nodes[n]
            {
                return w;
            }
            *out_id.entry(v).or_insert_with(|| {
                let p = match node_of.get(&v) {
                    Some(&n) => position[n],
                    None => frame.point(cdt.points[v]),
                };
                vertices.push(p);
                base + vertices.len() - 1
            })
        });
        let q = cdt.tris[t].map(|v| cdt.points[v]);
        stats.area += 0.5 * area2(&cdt, t);
        for i in 0..3 {
            let (a, b, c) = (q[i], q[(i + 1) % 3], q[(i + 2) % 3]);
            let (u, w) = ([b[0] - a[0], b[1] - a[1]], [c[0] - a[0], c[1] - a[1]]);
            let angle = (u[0] * w[1] - u[1] * w[0]).atan2(u[0] * w[0] + u[1] * w[1]);
            stats.min_angle_degrees = stats.min_angle_degrees.min(angle.to_degrees());
        }
        triangles.push(tri);
    }
    stats.triangles = triangles.len();
    let mut wall_splits: Vec<((usize, usize), Vec<usize>)> = Vec::new();
    wall_points.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
    for (edge, _, v) in wall_points {
        let Some(&out) = out_id.get(&v) else {
            return Err("a split wall-line point lies on no shell triangle".into());
        };
        match wall_splits.last_mut() {
            Some((e, list)) if *e == edge => list.push(out),
            _ => wall_splits.push((edge, vec![out])),
        }
    }
    let mut edges: Vec<f64> = Vec::new();
    let mut seen: HashSet<(usize, usize)> = HashSet::default();
    let at = |v: usize| {
        if v < base {
            wall_vertices[v]
        } else {
            vertices[v - base]
        }
    };
    for t in &triangles {
        for i in 0..3 {
            let (a, b) = (t[i], t[(i + 1) % 3]);
            if seen.insert(key(a, b)) {
                edges.push((at(a) - at(b)).length());
            }
        }
    }
    edges.sort_by(f64::total_cmp);
    if let (Some(&lo), Some(&hi)) = (edges.first(), edges.last()) {
        stats.edge_lengths = [lo, edges[edges.len() / 2], hi];
    }
    Ok(ShellSurface {
        vertices,
        triangles,
        wall_splits,
        stats,
    })
}

/// The segment nearest to `c` (in the plane).
fn nearest_segment(cdt: &Cdt, segs: &[(usize, usize, Segment)], c: P2) -> Option<usize> {
    let d = |&(a, b, _): &(usize, usize, Segment)| {
        let (p, q) = (cdt.points[a], cdt.points[b]);
        let e = [q[0] - p[0], q[1] - p[1]];
        let l2 = e[0] * e[0] + e[1] * e[1];
        let t = if l2 > 0.0 {
            (((c[0] - p[0]) * e[0] + (c[1] - p[1]) * e[1]) / l2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        planar::distance(c, [p[0] + e[0] * t, p[1] + e[1] * t])
    };
    (0..segs.len()).min_by(|&i, &j| d(&segs[i]).total_cmp(&d(&segs[j])))
}

fn area2(cdt: &Cdt, t: usize) -> f64 {
    let [a, b, c] = cdt.tris[t].map(|v| cdt.points[v]);
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

/// Distance along the outline between a corner place and a point place.
fn corner_gap(frame: &Frame, corner: f64, pos: f64) -> f64 {
    let k = frame.poly.len();
    let q = frame.poly[corner as usize % k];
    let i = pos.floor() as usize % k;
    let (a, b) = (frame.poly[i], frame.poly[(i + 1) % k]);
    let t = pos - pos.floor();
    let p = [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t];
    planar::distance(p, q)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sizing::shapes;

    // Keep rule (SW-18): a square plate inside a box (no wall contact) is
    // meshed over its whole area at the requested spacing, with no interior
    // edge between two free-edge vertices; a baffle wider than the box is
    // clipped to the box section, its outline running along the cut wall.
    // Breaks when the region test, the outline or the chord split goes
    // wrong; the duct tests only see the end result.
    #[test]
    fn plate_and_spanning_baffle_in_a_box() {
        let (v, t) = shapes::cuboid(Vec3::new(4.0, 2.0, 2.0));
        let keep = vec![false; t.len()];
        let group = vec![3; t.len()];
        let plate = vec![
            Vec3::new(1.0, 0.5, 0.7),
            Vec3::new(1.8, 0.5, 0.7),
            Vec3::new(1.8, 1.5, 0.7),
            Vec3::new(1.0, 1.5, 0.7),
        ];
        let baffle = vec![
            Vec3::new(3.0, -1.0, -1.0),
            Vec3::new(3.0, 3.0, -1.0),
            Vec3::new(3.0, 3.0, 1.2),
            Vec3::new(3.0, -1.0, 1.2),
        ];
        let im = imprint(
            &v,
            &t,
            &keep,
            &group,
            &[("plate", &plate), ("baffle", &baffle)],
        )
        .unwrap();
        assert!(im.lines.iter().all(|&(_, k)| k == 1), "{:?}", im.lines);
        // the baffle's line: the bottom and both sides of the box section
        let length: f64 = im
            .lines
            .iter()
            .map(|&((a, b), _)| (im.vertices[a] - im.vertices[b]).length())
            .sum();
        assert!((length - (2.0 + 2.0 * 1.2)).abs() < 1e-9, "{length}");
        let size = |_: Vec3| 0.2;
        let p = mesh_shell(
            &plate,
            &im.vertices,
            &im.triangles,
            &[],
            &size,
            im.vertices.len(),
        )
        .unwrap();
        assert!((p.stats.area - 0.8).abs() < 1e-9, "{:?}", p.stats);
        assert!(
            (p.stats.free_edge_length - 3.6).abs() < 1e-9,
            "{:?}",
            p.stats
        );
        assert!(p.stats.min_angle_degrees > 15.0, "{:?}", p.stats);
        let line: Vec<(usize, usize)> = im.lines.iter().map(|&(e, _)| e).collect();
        let b = mesh_shell(
            &baffle,
            &im.vertices,
            &im.triangles,
            &line,
            &size,
            im.vertices.len(),
        )
        .unwrap();
        assert!((b.stats.area - 2.0 * 1.2).abs() < 1e-9, "{:?}", b.stats);
        assert!(
            (b.stats.free_edge_length - 2.0).abs() < 1e-9,
            "{:?}",
            b.stats
        );
        assert_eq!(b.stats.free_edge_ends, 2);
        // no interior edge with both ends on a free edge
        let on_line: HashSet<(usize, usize)> = line.iter().copied().collect();
        for s in [&p, &b] {
            let mut uses: HashMap<(usize, usize), usize> = HashMap::default();
            for t in &s.triangles {
                for k in 0..3 {
                    *uses.entry(key(t[k], t[(k + 1) % 3])).or_default() += 1;
                }
            }
            let free_edge: HashSet<usize> = uses
                .iter()
                .filter(|&(e, &n)| n == 1 && !on_line.contains(e))
                .flat_map(|(&(a, b), _)| [a, b])
                .collect();
            for (&(a, b), &n) in &uses {
                assert!(
                    n == 1 || !(free_edge.contains(&a) && free_edge.contains(&b)),
                    "chord {a}-{b}"
                );
            }
        }
    }

    // Keep rule (SW-18): what the first shell slice refuses, with a reason
    // the engineer can act on: a shell crossing a port cap, two shells that
    // meet, a shell lying in a wall facet. Breaks when one of them slips
    // through to the scaffold, where it fails without a name.
    #[test]
    fn caps_crossing_and_coplanar_shells_are_refused() {
        let (v, t) = shapes::cuboid(Vec3::new(4.0, 2.0, 2.0));
        let group = vec![3; t.len()];
        let square = |x: f64| {
            vec![
                Vec3::new(x, 0.5, 0.5),
                Vec3::new(x, 1.5, 0.5),
                Vec3::new(x, 1.5, 1.5),
                Vec3::new(x, 0.5, 1.5),
            ]
        };
        // every facet a cap: a baffle across the box meets one
        let caps = vec![true; t.len()];
        let across = vec![
            Vec3::new(2.0, -1.0, -1.0),
            Vec3::new(2.0, 3.0, -1.0),
            Vec3::new(2.0, 3.0, 3.0),
            Vec3::new(2.0, -1.0, 3.0),
        ];
        let e = imprint(&v, &t, &caps, &group, &[("baffle", &across)]).unwrap_err();
        assert!(
            e.contains("inlet or outlet cap") && e.contains("baffle"),
            "{e}"
        );
        // two crossing plates
        let keep = vec![false; t.len()];
        let flat = vec![
            Vec3::new(1.5, 0.2, 1.0),
            Vec3::new(2.5, 0.2, 1.0),
            Vec3::new(2.5, 1.8, 1.0),
            Vec3::new(1.5, 1.8, 1.0),
        ];
        let e = imprint(&v, &t, &keep, &group, &[("a", &square(2.0)), ("b", &flat)]).unwrap_err();
        assert!(e.contains("\"a\"") && e.contains("\"b\""), "{e}");
        // a plate in the wall facet x = 0
        let e = imprint(&v, &t, &keep, &group, &[("c", &square(0.0))]).unwrap_err();
        assert!(e.contains("\"c\""), "{e}");
        // a plate apart from the others is fine
        imprint(
            &v,
            &t,
            &keep,
            &group,
            &[("a", &square(1.0)), ("d", &square(3.0))],
        )
        .unwrap();
    }

    // Keep rule (SW-18): the default shell size comes from the narrow
    // width of the polygon inside the STL's box, and a shell's size
    // primitives are its polygon plus only the parts of its outline inside
    // the fluid. Breaks when the clip, the width or the winding test of
    // the free edges goes wrong; the duct test sees only the end result.
    #[test]
    fn narrow_width_and_size_items_of_a_shell() {
        // a 4 x 1 plate half outside the box [0, 2]^3: the clipped part is
        // 1 x 1 (from x = 1 to 2)
        let plate = vec![
            Vec3::new(1.0, 0.5, 1.0),
            Vec3::new(5.0, 0.5, 1.0),
            Vec3::new(5.0, 1.5, 1.0),
            Vec3::new(1.0, 1.5, 1.0),
        ];
        let (lo, hi) = (Vec3::ZERO, Vec3::new(2.0, 2.0, 2.0));
        let w = narrow_width(&plate, lo, hi).unwrap().unwrap();
        assert!((w - 1.0).abs() < 1e-12, "{w}");
        let far: Vec<Vec3> = plate
            .iter()
            .map(|&p| p + Vec3::new(10.0, 0.0, 0.0))
            .collect();
        assert_eq!(narrow_width(&far, lo, hi).unwrap(), None);
        let (v, t) = shapes::cuboid(hi);
        let items = size_items(&plate, 0.2, 0.1, &v, &t).unwrap();
        let triangles = items
            .iter()
            .filter(|w| matches!(w.primitive, Primitive::Triangle(_)))
            .count();
        assert_eq!(triangles, 2);
        // free edges inside the box: y = 0.5 and y = 1.5 from x = 1 to 2,
        // and the side at x = 1 (in pieces of four edge sizes, so a piece
        // may reach 0.4 past the wall)
        let length: f64 = items
            .iter()
            .filter_map(|w| match w.primitive {
                Primitive::Segment([a, b]) => {
                    assert_eq!(w.size, 0.1);
                    assert!(a.x <= 2.4 + 1e-12 && b.x <= 2.4 + 1e-12, "{a:?} {b:?}");
                    Some((b - a).length())
                }
                Primitive::Triangle(_) => None,
            })
            .sum();
        assert!((3.0 - 1e-9..=3.8 + 1e-9).contains(&length), "{length}");
        assert!(size_items(&plate, 0.0, 0.1, &v, &t).is_err());
    }
}
