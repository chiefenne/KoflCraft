//! Overlay of two polygon tilings of one planar region (BS-16, step 3 of the
//! 3D tetrahedron and prism insertion): the common refinement of the grids
//! that the two sides of a non-conformal block face put on that face.
//!
//! Both sides of a block face are grids on the same surface, the face's
//! transfinite patch, so they can be overlaid in the patch's own
//! parameters (s, t) without any projection. Every face of the overlay lies
//! in one cell `a` of the first tiling and one cell `b` of the second, and
//! becomes the internal face between those two cells. Because each overlap
//! polygon is one face with one area seen from both sides, the flux through
//! the interface is conserved to round-off. This is the supermesh of
//! Farrell, Piggott, Pain, Gorman and Wilson ("Conservative interpolation
//! between unstructured meshes via supermesh construction", Computer
//! Methods in Applied Mechanics and Engineering 198, 2009), built here as a
//! planar arrangement of the two edge sets in the manner of Cherchi,
//! Livesu, Scateni and Attene ("Fast and robust mesh arrangements using
//! floating-point arithmetic", ACM Transactions on Graphics 39(6), 2020):
//! every topological decision (does a point lie on a segment, do two
//! segments cross) is taken by the exact orientation predicate, and only
//! the coordinates of new crossing points are rounded. No code was read.
//!
//! **Steps.** (1) Points of the two tilings closer than the snap distance
//! are merged, (2) a point on or within the snap distance of the other
//! tiling's edge is put into that edge, (3) edges that cross properly are
//! split at their crossing, (4) the faces of the planar graph are traced
//! and (5) each face takes the cell of each tiling on its left. A face
//! that is bounded only by edges of the other tiling inherits the cell
//! from its neighbours across those edges.
//!
//! **Snapping.** With `snap = 0` the overlay is exact: points merge only
//! when they are equal and lie on an edge only when the predicate says
//! collinear. With `snap > 0` a point merges with or joins a feature closer
//! than `snap` times the shortest edge at the point, so two grid lines that
//! nearly coincide give one line instead of a strip of sliver faces. The
//! largest distance snapped over is reported. Two points on the region's
//! border are judged by the border edges at them alone, because the border
//! is shared with the faces around the region, whose own overlays must
//! decide the same.
//!
//! The tilings need not be convex cells; they must cover the same region
//! and each must be conforming (every interior edge shared by two cells
//! with opposite orientation).

#![allow(clippy::needless_range_loop)] // parallel per-point arrays

use std::cmp::Ordering;
use std::collections::{HashMap, VecDeque};

use crate::predicates::orient2d;

/// A point of the plane.
pub type P2 = [f64; 2];

/// A tiling: points and counterclockwise cell rings.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Tiling {
    pub points: Vec<P2>,
    pub cells: Vec<Vec<usize>>,
}

/// Where a point of the overlay comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Origin {
    /// A point of the first tiling only.
    A(usize),
    /// A point of the second tiling only.
    B(usize),
    /// A point of each tiling, merged.
    Both(usize, usize),
    /// The crossing of an edge of each tiling.
    Crossing,
}

/// One face of the overlay: the part of cell `a` of the first tiling that
/// lies in cell `b` of the second, as a counterclockwise ring.
#[derive(Clone, Debug, PartialEq)]
pub struct OverlapFace {
    pub a: usize,
    pub b: usize,
    pub ring: Vec<usize>,
}

/// The overlay of two tilings.
#[derive(Clone, Debug, PartialEq)]
pub struct Overlay {
    pub points: Vec<P2>,
    pub origin: Vec<Origin>,
    /// The overlay point of every point of the first tiling.
    pub a_point: Vec<usize>,
    /// The overlay point of every point of the second tiling.
    pub b_point: Vec<usize>,
    pub faces: Vec<OverlapFace>,
    /// For every edge (u, v), u < v, of the first tiling: the overlay
    /// points strictly inside it, in order from u to v. These points must
    /// be put into every other face that has the edge, so the cells close.
    pub a_edges: HashMap<(usize, usize), Vec<usize>>,
    /// The same for the second tiling.
    pub b_edges: HashMap<(usize, usize), Vec<usize>>,
    /// The largest distance a point was snapped over (0 for an exact
    /// overlay).
    pub max_snap: f64,
    /// The largest |sum of a cell's overlay faces / cell area − 1| over the
    /// cells of both tilings.
    pub max_area_error: f64,
}

/// Why two tilings could not be overlaid.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OverlayError(pub String);

impl std::fmt::Display for OverlayError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for OverlayError {}

fn fail<T>(message: impl Into<String>) -> Result<T, OverlayError> {
    Err(OverlayError(message.into()))
}

/// Twice the signed area of a ring.
pub fn ring_area2(points: &[P2], ring: &[usize]) -> f64 {
    let n = ring.len();
    (0..n)
        .map(|i| {
            let (p, q) = (points[ring[i]], points[ring[(i + 1) % n]]);
            p[0] * q[1] - p[1] * q[0]
        })
        .sum()
}

fn sub(a: P2, b: P2) -> P2 {
    [a[0] - b[0], a[1] - b[1]]
}

fn dot(a: P2, b: P2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

fn dist(a: P2, b: P2) -> f64 {
    let d = sub(a, b);
    dot(d, d).sqrt()
}

/// The directed edges of a tiling with the cell on their left, checked for
/// conformity.
fn left_cells(t: &Tiling, name: &str) -> Result<HashMap<(usize, usize), usize>, OverlayError> {
    let mut left = HashMap::new();
    for (c, ring) in t.cells.iter().enumerate() {
        if ring.len() < 3 || ring.iter().any(|&v| v >= t.points.len()) {
            return fail(format!("{name} cell {c} is not a polygon of its points"));
        }
        if ring_area2(&t.points, ring) <= 0.0 {
            return fail(format!("{name} cell {c} is not counterclockwise"));
        }
        for i in 0..ring.len() {
            let e = (ring[i], ring[(i + 1) % ring.len()]);
            if e.0 == e.1 {
                return fail(format!("{name} cell {c} repeats a point"));
            }
            if left.insert(e, c).is_some() {
                return fail(format!(
                    "{name} edge {}–{} is used twice in one direction",
                    e.0, e.1
                ));
            }
        }
    }
    Ok(left)
}

/// Union–find over the merged points.
struct Merge(Vec<usize>);

impl Merge {
    fn find(&mut self, mut i: usize) -> usize {
        while self.0[i] != i {
            self.0[i] = self.0[self.0[i]];
            i = self.0[i];
        }
        i
    }
}

/// A uniform grid of bins over the bounding box, for candidate pairs.
struct Bins {
    lo: P2,
    size: f64,
    n: [usize; 2],
    items: Vec<Vec<usize>>,
}

impl Bins {
    fn new(points: &[P2], size: f64) -> Self {
        let mut lo = [f64::INFINITY; 2];
        let mut hi = [f64::NEG_INFINITY; 2];
        for p in points {
            for k in 0..2 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        let size = if size > 0.0 && size.is_finite() {
            size
        } else {
            1.0
        };
        let n = [0, 1].map(|k| (((hi[k] - lo[k]) / size).ceil() as usize).clamp(1, 4096));
        Self {
            lo,
            size,
            n,
            items: vec![Vec::new(); n[0] * n[1]],
        }
    }

    fn range(&self, lo: P2, hi: P2) -> [std::ops::RangeInclusive<usize>; 2] {
        [0, 1].map(|k| {
            let a = ((lo[k] - self.lo[k]) / self.size).floor().max(0.0) as usize;
            let b = ((hi[k] - self.lo[k]) / self.size).floor().max(0.0) as usize;
            a.min(self.n[k] - 1)..=b.min(self.n[k] - 1)
        })
    }

    fn insert(&mut self, item: usize, lo: P2, hi: P2) {
        let [ri, rj] = self.range(lo, hi);
        for j in rj {
            for i in ri.clone() {
                self.items[j * self.n[0] + i].push(item);
            }
        }
    }

    fn query(&self, lo: P2, hi: P2, out: &mut Vec<usize>) {
        out.clear();
        let [ri, rj] = self.range(lo, hi);
        for j in rj {
            for i in ri.clone() {
                out.extend_from_slice(&self.items[j * self.n[0] + i]);
            }
        }
        out.sort_unstable();
        out.dedup();
    }
}

fn bbox(a: P2, b: P2, pad: f64) -> (P2, P2) {
    (
        [a[0].min(b[0]) - pad, a[1].min(b[1]) - pad],
        [a[0].max(b[0]) + pad, a[1].max(b[1]) + pad],
    )
}

/// An edge of the arrangement from `p` to `q`, with the directed edge of
/// each tiling it runs along (in the direction p → q), if any.
#[derive(Clone, Debug)]
struct Edge {
    p: usize,
    q: usize,
    a: Option<(usize, usize)>,
    b: Option<(usize, usize)>,
}

/// Overlay two tilings of one region (see the module notes). `snap` is
/// relative to the shortest edge at a point; 0 gives the exact overlay.
pub fn overlay(a: &Tiling, b: &Tiling, snap: f64) -> Result<Overlay, OverlayError> {
    if !(0.0..0.25).contains(&snap) {
        return fail(format!("snap {snap} must be at least 0 and below 0.25"));
    }
    let a_left = left_cells(a, "first tiling")?;
    let b_left = left_cells(b, "second tiling")?;
    let na = a.points.len();
    let mut pts: Vec<P2> = a.points.iter().chain(&b.points).copied().collect();
    let side = |i: usize| usize::from(i >= na); // 0 = A, 1 = B
    // the original segments, one per undirected edge, in global point ids
    let mut segs: Vec<(usize, usize, usize)> = Vec::new(); // (p, q, side)
    for (s, left) in [(0, &a_left), (1, &b_left)] {
        let off = if s == 0 { 0 } else { na };
        let mut keys: Vec<(usize, usize)> =
            left.keys().map(|&(u, v)| (u.min(v), u.max(v))).collect();
        keys.sort_unstable();
        keys.dedup();
        segs.extend(keys.into_iter().map(|(u, v)| (u + off, v + off, s)));
    }
    // the shortest edge at every point, and the border level (1 on the
    // region's border, 2 at a corner of it)
    let mut local = vec![f64::INFINITY; pts.len()];
    for &(p, q, _) in &segs {
        let l = dist(pts[p], pts[q]);
        if l <= 0.0 {
            return fail("a tiling has an edge of zero length");
        }
        local[p] = local[p].min(l);
        local[q] = local[q].min(l);
    }
    let mut border: Vec<Vec<usize>> = vec![Vec::new(); pts.len()];
    for (s, left) in [(0, &a_left), (1, &b_left)] {
        let off = if s == 0 { 0 } else { na };
        for &(u, v) in left.keys() {
            if !left.contains_key(&(v, u)) {
                border[u + off].push(v + off);
                border[v + off].push(u + off);
            }
        }
    }
    let level: Vec<u8> = (0..pts.len())
        .map(|i| match border[i].as_slice() {
            [] => 0,
            [x, y] if orient2d(pts[*x], pts[i], pts[*y]) == Ordering::Equal => 1,
            _ => 2,
        })
        .collect();
    let mean = segs
        .iter()
        .map(|&(p, q, _)| dist(pts[p], pts[q]))
        .sum::<f64>()
        / segs.len().max(1) as f64;
    let eps = |i: usize| snap * local[i];
    // the shortest border edge at a border point
    let along: Vec<f64> = (0..pts.len())
        .map(|i| {
            border[i]
                .iter()
                .map(|&j| dist(pts[i], pts[j]))
                .fold(f64::INFINITY, f64::min)
        })
        .collect();
    let mut max_snap: f64 = 0.0;

    // (1) merge points of the two tilings
    let mut merge = Merge((0..pts.len()).collect());
    {
        let mut bins = Bins::new(&pts, mean);
        for j in na..pts.len() {
            bins.insert(j, pts[j], pts[j]);
        }
        let mut cand = Vec::new();
        for i in 0..na {
            let r = eps(i) + snap * mean;
            let (lo, hi) = bbox(pts[i], pts[i], r);
            bins.query(lo, hi, &mut cand);
            let best = cand
                .iter()
                .map(|&j| (j, dist(pts[i], pts[j])))
                // two border points are judged by the border edges alone:
                // a border is shared with other faces, whose overlays must
                // decide alike
                .filter(|&(j, d)| {
                    let tol = if level[i] > 0 && level[j] > 0 {
                        snap * along[i].min(along[j])
                    } else {
                        eps(i).min(eps(j))
                    };
                    d <= tol
                })
                .min_by(|x, y| x.1.total_cmp(&y.1));
            if let Some((j, _)) = best {
                let (ri, rj) = (merge.find(i), merge.find(j));
                if ri == rj {
                    continue;
                }
                if side(ri) == side(rj) {
                    // j already merged with another point of A
                    continue;
                }
                // the merged point keeps the more constrained position
                let keep = if level[j] > level[i] { j } else { i };
                let other = if keep == i { j } else { i };
                max_snap = max_snap.max(dist(pts[i], pts[j]));
                pts[other] = pts[keep];
                merge.0[ri] = keep;
                merge.0[rj] = keep;
                merge.0[other] = keep;
                merge.0[keep] = keep;
            }
        }
    }
    let rep: Vec<usize> = (0..pts.len()).map(|i| merge.find(i)).collect();
    // which sides a representative carries
    let mut carries = vec![[false; 2]; pts.len()];
    for i in 0..pts.len() {
        carries[rep[i]][side(i)] = true;
    }

    // (2) points on edges of the other tiling
    let mut chains: Vec<Vec<usize>> = Vec::with_capacity(segs.len());
    {
        let mut bins = Bins::new(&pts, mean);
        let reps: Vec<usize> = (0..pts.len()).filter(|&i| rep[i] == i).collect();
        for &w in &reps {
            bins.insert(w, pts[w], pts[w]);
        }
        let mut cand = Vec::new();
        for &(p0, q0, s) in &segs {
            let (p, q) = (rep[p0], rep[q0]);
            if p == q {
                return fail("two points of one edge snapped together; lower snap");
            }
            let (pp, qq) = (pts[p], pts[q]);
            let d = sub(qq, pp);
            let len2 = dot(d, d);
            let pad = snap * local[p0].min(local[q0]);
            let (lo, hi) = bbox(pp, qq, pad);
            bins.query(lo, hi, &mut cand);
            let mut on: Vec<(f64, usize)> = Vec::new();
            for &w in &cand {
                if w == p || w == q || carries[w][s] {
                    continue;
                }
                let x = pts[w];
                let t = dot(sub(x, pp), d) / len2;
                let hit = if snap == 0.0 {
                    orient2d(pp, qq, x) == Ordering::Equal
                        && t > 0.0
                        && t < 1.0
                        && x[0] >= pp[0].min(qq[0])
                        && x[0] <= pp[0].max(qq[0])
                        && x[1] >= pp[1].min(qq[1])
                        && x[1] <= pp[1].max(qq[1])
                } else {
                    let foot = [pp[0] + t * d[0], pp[1] + t * d[1]];
                    let off = dist(x, foot);
                    let tol = eps(w).min(pad);
                    let inside = t > 0.0 && t < 1.0;
                    if inside && off <= tol {
                        max_snap = max_snap.max(off);
                        true
                    } else {
                        false
                    }
                };
                if hit {
                    on.push((t, w));
                }
            }
            on.sort_by(|x, y| x.0.total_cmp(&y.0));
            let mut chain = Vec::with_capacity(on.len() + 2);
            chain.push(p);
            chain.extend(on.into_iter().map(|(_, w)| w));
            chain.push(q);
            chains.push(chain);
        }
    }

    // sub-edges, one per undirected pair, with the tiling edge(s) along them
    let mut edges: Vec<Edge> = Vec::new();
    let mut edge_of: HashMap<(usize, usize), usize> = HashMap::new();
    let original = |seg: &(usize, usize, usize)| -> (usize, usize) {
        let off = if seg.2 == 0 { 0 } else { na };
        (seg.0 - off, seg.1 - off)
    };
    for (si, chain) in chains.iter().enumerate() {
        let seg = segs[si];
        let (u, v) = original(&seg);
        for w in chain.windows(2) {
            let (x, y) = (w[0], w[1]);
            let key = (x.min(y), x.max(y));
            let id = *edge_of.entry(key).or_insert_with(|| {
                edges.push(Edge {
                    p: key.0,
                    q: key.1,
                    a: None,
                    b: None,
                });
                edges.len() - 1
            });
            let dir = if x == edges[id].p { (u, v) } else { (v, u) };
            let slot = if seg.2 == 0 {
                &mut edges[id].a
            } else {
                &mut edges[id].b
            };
            if slot.is_some_and(|d| d != dir) {
                return fail("an edge of one tiling overlaps another of the same tiling");
            }
            *slot = Some(dir);
        }
    }

    // (3) proper crossings between edges of A only and edges of B only
    let mut cross: Vec<Vec<(f64, usize)>> = vec![Vec::new(); edges.len()];
    {
        let mut bins = Bins::new(&pts, mean);
        for (i, e) in edges.iter().enumerate() {
            if e.a.is_some() && e.b.is_none() {
                let (lo, hi) = bbox(pts[e.p], pts[e.q], 0.0);
                bins.insert(i, lo, hi);
            }
        }
        let mut cand = Vec::new();
        for j in 0..edges.len() {
            let ej = edges[j].clone();
            if !(ej.b.is_some() && ej.a.is_none()) {
                continue;
            }
            let (lo, hi) = bbox(pts[ej.p], pts[ej.q], 0.0);
            bins.query(lo, hi, &mut cand);
            for &i in &cand {
                let ei = &edges[i];
                if ei.p == ej.p || ei.p == ej.q || ei.q == ej.p || ei.q == ej.q {
                    continue;
                }
                let (p, q, r, s) = (pts[ei.p], pts[ei.q], pts[ej.p], pts[ej.q]);
                let o1 = orient2d(p, q, r);
                let o2 = orient2d(p, q, s);
                let o3 = orient2d(r, s, p);
                let o4 = orient2d(r, s, q);
                let proper = o1 != Ordering::Equal
                    && o2 != Ordering::Equal
                    && o1 != o2
                    && o3 != Ordering::Equal
                    && o4 != Ordering::Equal
                    && o3 != o4;
                if !proper {
                    // touching points were put into the edges in step 2
                    continue;
                }
                let d1 = sub(q, p);
                let d2 = sub(s, r);
                let den = d1[0] * d2[1] - d1[1] * d2[0];
                let w = sub(r, p);
                let t = ((w[0] * d2[1] - w[1] * d2[0]) / den).clamp(0.0, 1.0);
                let u = ((w[0] * d1[1] - w[1] * d1[0]) / den).clamp(0.0, 1.0);
                let x = [p[0] + t * d1[0], p[1] + t * d1[1]];
                pts.push(x);
                let id = pts.len() - 1;
                cross[i].push((t, id));
                cross[j].push((u, id));
            }
        }
    }
    let n_points = pts.len();

    // the final edges, split at their crossings
    let mut fin: Vec<Edge> = Vec::with_capacity(edges.len());
    let mut pieces: Vec<Vec<usize>> = Vec::with_capacity(edges.len()); // chain p..q per sub-edge
    for (i, e) in edges.iter().enumerate() {
        let mut c = std::mem::take(&mut cross[i]);
        c.sort_by(|x, y| x.0.total_cmp(&y.0));
        let mut chain = vec![e.p];
        chain.extend(c.into_iter().map(|(_, id)| id));
        chain.push(e.q);
        for w in chain.windows(2) {
            fin.push(Edge {
                p: w[0],
                q: w[1],
                a: e.a,
                b: e.b,
            });
        }
        pieces.push(chain);
    }

    // (4) half-edges and faces: half-edge 2k runs p → q of fin[k], 2k + 1
    // back; around every point the outgoing half-edges by angle
    let mut out: Vec<Vec<usize>> = vec![Vec::new(); n_points];
    for (k, e) in fin.iter().enumerate() {
        out[e.p].push(2 * k);
        out[e.q].push(2 * k + 1);
    }
    let from = |h: usize| {
        let e = &fin[h / 2];
        if h.is_multiple_of(2) { e.p } else { e.q }
    };
    let to = |h: usize| {
        let e = &fin[h / 2];
        if h.is_multiple_of(2) { e.q } else { e.p }
    };
    let mut pos = vec![0usize; 2 * fin.len()];
    for (v, list) in out.iter_mut().enumerate() {
        let o = pts[v];
        list.sort_by(|&x, &y| {
            let (dx, dy) = (sub(pts[to(x)], o), sub(pts[to(y)], o));
            dx[1].atan2(dx[0]).total_cmp(&dy[1].atan2(dy[0]))
        });
        for (i, &h) in list.iter().enumerate() {
            pos[h] = i;
        }
    }
    let next = |h: usize| {
        let v = to(h);
        let list = &out[v];
        let twin = h ^ 1;
        list[(pos[twin] + list.len() - 1) % list.len()]
    };
    let mut face_of = vec![usize::MAX; 2 * fin.len()];
    let mut loops: Vec<Vec<usize>> = Vec::new();
    for h0 in 0..2 * fin.len() {
        if face_of[h0] != usize::MAX {
            continue;
        }
        let f = loops.len();
        let mut lp = Vec::new();
        let mut h = h0;
        while face_of[h] == usize::MAX {
            face_of[h] = f;
            lp.push(h);
            h = next(h);
            if lp.len() > 2 * fin.len() {
                return fail("face tracing did not close");
            }
        }
        if h != h0 {
            return fail("face tracing met a half-edge of another face");
        }
        loops.push(lp);
    }
    let rings: Vec<Vec<usize>> = loops
        .iter()
        .map(|lp| lp.iter().map(|&h| from(h)).collect())
        .collect();
    let area2: Vec<f64> = rings.iter().map(|r| ring_area2(&pts, r)).collect();
    let outer: Vec<bool> = area2.iter().map(|&x| x <= 0.0).collect();
    let n_outer = outer.iter().filter(|&&o| o).count();
    if n_outer != 1 {
        return fail(format!(
            "the overlay has {n_outer} unbounded or empty faces; the tilings must cover one region without holes"
        ));
    }

    // (5) labels
    let label = |which: usize| -> Result<Vec<Option<usize>>, OverlayError> {
        let left = if which == 0 { &a_left } else { &b_left };
        let mut lab: Vec<Option<usize>> = vec![None; loops.len()];
        for (f, lp) in loops.iter().enumerate() {
            for &h in lp {
                let e = &fin[h / 2];
                let dir = if which == 0 { e.a } else { e.b };
                let Some((u, v)) = dir else { continue };
                let d = if h.is_multiple_of(2) { (u, v) } else { (v, u) };
                let c = left.get(&d).copied();
                match (c, outer[f]) {
                    (Some(c), false) => {
                        if lab[f].is_some_and(|x| x != c) {
                            return fail("an overlay face lies in two cells of one tiling");
                        }
                        lab[f] = Some(c);
                    }
                    (None, false) => {
                        return fail(
                            "an overlay face lies outside one tiling; the tilings must cover the same region",
                        );
                    }
                    (Some(_), true) => {
                        return fail("a cell lies outside the region of the other tiling");
                    }
                    (None, true) => {}
                }
            }
        }
        let mut queue: VecDeque<usize> = (0..loops.len()).filter(|&f| lab[f].is_some()).collect();
        while let Some(f) = queue.pop_front() {
            for &h in &loops[f] {
                let e = &fin[h / 2];
                let along = if which == 0 { e.a } else { e.b };
                if along.is_some() {
                    continue;
                }
                let g = face_of[h ^ 1];
                if outer[g] {
                    return fail("a cell of one tiling reaches outside the other");
                }
                if lab[g].is_none() {
                    lab[g] = lab[f];
                    queue.push_back(g);
                }
            }
        }
        Ok(lab)
    };
    let la = label(0)?;
    let lb = label(1)?;

    // compact the points that are used
    let mut used = vec![usize::MAX; n_points];
    let mut points = Vec::new();
    let mut origin = Vec::new();
    let mut take = |v: usize, points: &mut Vec<P2>, origin: &mut Vec<Origin>| {
        if used[v] == usize::MAX {
            used[v] = points.len();
            points.push(pts[v]);
            origin.push(Origin::Crossing);
        }
        used[v]
    };
    let mut faces = Vec::new();
    for (f, ring) in rings.iter().enumerate() {
        if outer[f] {
            continue;
        }
        let (Some(ca), Some(cb)) = (la[f], lb[f]) else {
            return fail("an overlay face has no cell in one tiling");
        };
        let ring = ring
            .iter()
            .map(|&v| take(v, &mut points, &mut origin))
            .collect();
        faces.push(OverlapFace { a: ca, b: cb, ring });
    }
    let mut a_point = vec![usize::MAX; na];
    let mut b_point = vec![usize::MAX; b.points.len()];
    for i in 0..rep.len() {
        let v = take(rep[i], &mut points, &mut origin);
        if i < na {
            a_point[i] = v;
        } else {
            b_point[i - na] = v;
        }
    }
    for (i, &v) in a_point.iter().enumerate() {
        origin[v] = match origin[v] {
            Origin::B(j) => Origin::Both(i, j),
            _ => Origin::A(i),
        };
    }
    for (j, &v) in b_point.iter().enumerate() {
        origin[v] = match origin[v] {
            Origin::A(i) => Origin::Both(i, j),
            _ => Origin::B(j),
        };
    }

    // the points inside every tiling edge, from its first point
    let mut a_edges = HashMap::new();
    let mut b_edges = HashMap::new();
    for (si, chain) in chains.iter().enumerate() {
        let seg = segs[si];
        let mut inner: Vec<usize> = Vec::new();
        for (k, w) in chain.windows(2).enumerate() {
            let (x, y) = (w[0], w[1]);
            let id = edge_of[&(x.min(y), x.max(y))];
            let piece = &pieces[id];
            let mid: Vec<usize> = if piece[0] == x {
                piece[1..piece.len() - 1].to_vec()
            } else {
                piece[1..piece.len() - 1].iter().rev().copied().collect()
            };
            if k > 0 {
                inner.push(x);
            }
            inner.extend(mid);
        }
        let inner: Vec<usize> = inner
            .into_iter()
            .map(|v| take(v, &mut points, &mut origin))
            .collect();
        let key = original(&seg);
        if seg.2 == 0 {
            a_edges.insert(key, inner);
        } else {
            b_edges.insert(key, inner);
        }
    }

    // area check per cell
    let mut sum_a = vec![0.0; a.cells.len()];
    let mut sum_b = vec![0.0; b.cells.len()];
    for f in &faces {
        let ar = ring_area2(&points, &f.ring);
        sum_a[f.a] += ar;
        sum_b[f.b] += ar;
    }
    let mut max_area_error: f64 = 0.0;
    for (t, sums) in [(a, &sum_a), (b, &sum_b)] {
        for (c, ring) in t.cells.iter().enumerate() {
            let ar = ring_area2(&t.points, ring);
            max_area_error = max_area_error.max((sums[c] / ar - 1.0).abs());
        }
    }
    Ok(Overlay {
        points,
        origin,
        a_point,
        b_point,
        faces,
        a_edges,
        b_edges,
        max_snap,
        max_area_error,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A tensor grid of the unit square with the given lines.
    fn grid(xs: &[f64], ys: &[f64]) -> Tiling {
        let mut points = Vec::new();
        for &y in ys {
            for &x in xs {
                points.push([x, y]);
            }
        }
        let n = xs.len();
        let mut cells = Vec::new();
        for j in 0..ys.len() - 1 {
            for i in 0..n - 1 {
                let v = j * n + i;
                cells.push(vec![v, v + 1, v + n + 1, v + n]);
            }
        }
        Tiling { points, cells }
    }

    fn lines(n: usize) -> Vec<f64> {
        (0..=n).map(|i| i as f64 / n as f64).collect()
    }

    fn total(o: &Overlay) -> f64 {
        o.faces
            .iter()
            .map(|f| ring_area2(&o.points, &f.ring))
            .sum::<f64>()
            / 2.0
    }

    /// Invariant: two tensor grids with distinct lines give one face per
    /// pair of overlapping cells, every face convex and inside both cells,
    /// and the areas add up per cell. A wrong `next` rule or label would
    /// change the count or the sums; nothing else tests the kernel.
    #[test]
    fn distinct_grids_give_the_product_of_columns_and_rows() {
        let a = grid(&lines(3), &lines(3));
        let b = grid(&lines(4), &lines(5));
        let o = overlay(&a, &b, 0.0).unwrap();
        // x lines: 1/3, 2/3 and 1/4, 1/2, 3/4 → 6 columns; y: 2 + 4 → 7 rows
        assert_eq!(o.faces.len(), 42);
        assert!((total(&o) - 1.0).abs() < 1e-14);
        assert!(o.max_area_error < 1e-12, "{}", o.max_area_error);
        assert_eq!(o.max_snap, 0.0);
        for f in &o.faces {
            assert!(ring_area2(&o.points, &f.ring) > 0.0);
        }
    }

    /// Invariant: lines both grids share give no extra faces and no
    /// crossing points on them (exact mode); merged points are `Both`.
    #[test]
    fn shared_lines_merge_exactly() {
        let a = grid(&lines(2), &lines(2));
        let b = grid(&lines(4), &lines(4));
        let o = overlay(&a, &b, 0.0).unwrap();
        assert_eq!(o.faces.len(), 16);
        assert_eq!(o.points.len(), 25);
        assert!(o.origin.iter().all(|x| !matches!(x, Origin::Crossing)));
        assert_eq!(
            o.origin
                .iter()
                .filter(|x| matches!(x, Origin::Both(..)))
                .count(),
            9
        );
        // the middle edge of A along y = 0 holds B's point at x = 1/4
        let e = &o.a_edges[&(0, 1)];
        assert_eq!(e.len(), 1);
        assert_eq!(o.points[e[0]], [0.25, 0.0]);
    }

    /// Invariant: sheared, non-axis-aligned cells (the shape of
    /// tri-center sub-faces in a face's parameters) overlay with
    /// conserved areas, and every edge's inner points lie on it in order.
    #[test]
    fn skewed_cells_conserve_area_and_order_edge_points() {
        let mut a = grid(&lines(5), &lines(4));
        for p in &mut a.points {
            if p[0] > 0.0 && p[0] < 1.0 && p[1] > 0.0 && p[1] < 1.0 {
                p[0] += 0.07 * (7.0 * p[1]).sin();
                p[1] += 0.05 * (5.0 * p[0]).cos();
            }
        }
        let b = grid(&lines(7), &lines(3));
        let o = overlay(&a, &b, 0.0).unwrap();
        assert!((total(&o) - 1.0).abs() < 1e-13);
        assert!(o.max_area_error < 1e-12, "{}", o.max_area_error);
        for (edges, t) in [(&o.a_edges, &a), (&o.b_edges, &b)] {
            let map = if std::ptr::eq(t, &a) {
                &o.a_point
            } else {
                &o.b_point
            };
            for (&(u, v), inner) in edges {
                let (p, q) = (o.points[map[u]], o.points[map[v]]);
                let d = sub(q, p);
                let mut last = 0.0;
                for &w in inner {
                    let s = dot(sub(o.points[w], p), d) / dot(d, d);
                    assert!(s > last && s < 1.0, "edge {u}-{v}");
                    last = s;
                }
            }
        }
    }

    /// Invariant: the three quadrilaterals of a tri-center split of a
    /// triangle-shaped face (its flat corner at (1, 1) of the face
    /// parameters) overlay a neighbour's full grid with each overlay face in
    /// one quadrilateral, areas conserved: the interface of a tet piece.
    #[test]
    fn tri_center_sub_faces_overlay_a_full_grid() {
        let a = Tiling {
            points: vec![
                [0.0, 0.0],
                [0.5, 0.0],
                [1.0, 0.0],
                [0.0, 0.5],
                [0.45, 0.45],
                [1.0, 1.0],
                [0.0, 1.0],
            ],
            cells: vec![vec![0, 1, 4, 3], vec![1, 2, 5, 4], vec![3, 4, 5, 6]],
        };
        let b = grid(&lines(7), &lines(6));
        let o = overlay(&a, &b, 0.02).unwrap();
        assert!((total(&o) - 1.0).abs() < 1e-12);
        assert!(o.max_area_error < 1e-9, "{}", o.max_area_error);
        for c in 0..3 {
            assert!(o.faces.iter().any(|f| f.a == c));
        }
        assert!(o.faces.iter().all(|f| ring_area2(&o.points, &f.ring) > 0.0));
    }

    /// Invariant: each cell, with the inner points of its edges put into
    /// its ring, has the area of its overlay faces: the closure the 3D
    /// cells need. Breaking `a_edges` (a point missing) breaks it.
    #[test]
    fn cells_with_edge_points_close_on_their_faces() {
        let a = grid(&lines(3), &lines(2));
        let b = grid(&lines(5), &lines(4));
        let o = overlay(&a, &b, 0.0).unwrap();
        for (c, ring) in a.cells.iter().enumerate() {
            let mut full = Vec::new();
            for i in 0..ring.len() {
                let (u, v) = (ring[i], ring[(i + 1) % ring.len()]);
                full.push(o.a_point[u]);
                let inner = &o.a_edges[&(u.min(v), u.max(v))];
                if u < v {
                    full.extend(inner);
                } else {
                    full.extend(inner.iter().rev());
                }
            }
            // every point of the cell's faces on its border is in `full`
            for f in o.faces.iter().filter(|f| f.a == c) {
                for &w in &f.ring {
                    let p = o.points[w];
                    let on_border = ring.iter().enumerate().any(|(i, &u)| {
                        let v = ring[(i + 1) % ring.len()];
                        orient2d(a.points[u], a.points[v], p) == Ordering::Equal
                    });
                    if on_border {
                        assert!(full.contains(&w), "cell {c} misses point {w}");
                    }
                }
            }
        }
    }

    /// Invariant: with snapping, two lines 1e-5 apart give one line, not a
    /// strip of slivers; without, the slivers are there.
    #[test]
    fn snap_removes_near_coincident_slivers() {
        let a = grid(&lines(2), &lines(2));
        let b = grid(&[0.0, 0.5 + 1e-5, 1.0], &[0.0, 0.3, 1.0]);
        let exact = overlay(&a, &b, 0.0).unwrap();
        let snapped = overlay(&a, &b, 0.02).unwrap();
        assert_eq!(exact.faces.len(), 9);
        assert_eq!(snapped.faces.len(), 6);
        assert!(snapped.max_snap > 0.0 && snapped.max_snap < 2e-5);
        let min = |o: &Overlay| {
            o.faces
                .iter()
                .map(|f| ring_area2(&o.points, &f.ring))
                .fold(f64::INFINITY, f64::min)
        };
        assert!(min(&exact) < 1e-4);
        assert!(min(&snapped) > 0.1);
        assert!((total(&snapped) - 1.0).abs() < 1e-12);
    }

    /// Invariant: tilings of different regions are refused, not overlaid
    /// into a partial result (AG-08).
    #[test]
    fn different_regions_are_refused() {
        let a = grid(&lines(2), &lines(2));
        let b = grid(&[0.0, 0.5], &lines(2));
        assert!(overlay(&a, &b, 0.0).is_err());
        let mut cw = grid(&lines(1), &lines(1));
        cw.cells[0].reverse();
        assert!(overlay(&a, &cw, 0.0).is_err());
    }
}
