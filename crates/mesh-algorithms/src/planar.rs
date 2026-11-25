//! A constrained Delaunay triangulation in the plane, for meshing a planar
//! region such as a zero-thickness shell (SHL-01).
//!
//! Method: incremental Bowyer–Watson insertion (Bowyer 1981, Watson 1981)
//! inside a box, the cavity grown until it is star-shaped from the new point
//! under exact orientations (Shewchuk 1997), so a rounded in-circle test can
//! never leave a folded triangle; a cavity never crosses a constraint, which
//! makes each insertion a constrained Delaunay update (Chew 1989).
//! Constraint edges are recovered by edge flips (Sloan 1993; Anglada 1997)
//! and the triangulation is made constrained Delaunay again by Lawson flips.
//! Regions are the edge-connected sets of triangles between constraints.

use crate::delaunay::{HashMap, HashSet};
use crate::predicates::orient3d;
use mesh_core::geometry::Vec3;
use std::cmp::Ordering;
use std::collections::VecDeque;

pub type P2 = [f64; 2];
pub const NONE: usize = usize::MAX;

/// Why an insertion did not happen.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Refused {
    /// The point is already a vertex (its index).
    Duplicate(usize),
    Outside,
    Constrained,
}

/// Exact orientation: `Greater` when a, b, c turn counter-clockwise.
pub fn orient(a: P2, b: P2, c: P2) -> Ordering {
    orient3d(
        Vec3::new(a[0], a[1], 0.0),
        Vec3::new(b[0], b[1], 0.0),
        Vec3::new(c[0], c[1], 0.0),
        Vec3::new(a[0], a[1], 1.0),
    )
}

/// Positive when d lies inside the circumcircle of the counter-clockwise
/// triangle a, b, c (rounded; the cavity repair makes it safe).
fn incircle(a: P2, b: P2, c: P2, d: P2) -> f64 {
    let (adx, ady) = (a[0] - d[0], a[1] - d[1]);
    let (bdx, bdy) = (b[0] - d[0], b[1] - d[1]);
    let (cdx, cdy) = (c[0] - d[0], c[1] - d[1]);
    let al = adx * adx + ady * ady;
    let bl = bdx * bdx + bdy * bdy;
    let cl = cdx * cdx + cdy * cdy;
    adx * (bdy * cl - cdy * bl) - ady * (bdx * cl - cdx * bl) + al * (bdx * cdy - bdy * cdx)
}

/// Circumcentre and circumradius of a triangle.
pub fn circumcircle(a: P2, b: P2, c: P2) -> Option<(P2, f64)> {
    let (bx, by) = (b[0] - a[0], b[1] - a[1]);
    let (cx, cy) = (c[0] - a[0], c[1] - a[1]);
    let d = 2.0 * (bx * cy - by * cx);
    if d.abs() <= f64::MIN_POSITIVE {
        return None;
    }
    let (b2, c2) = (bx * bx + by * by, cx * cx + cy * cy);
    let ux = (cy * b2 - by * c2) / d;
    let uy = (bx * c2 - cx * b2) / d;
    let r = (ux * ux + uy * uy).sqrt();
    (r.is_finite()).then_some(([a[0] + ux, a[1] + uy], r))
}

pub fn distance(a: P2, b: P2) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt()
}

fn key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

pub struct Cdt {
    pub points: Vec<P2>,
    /// Counter-clockwise triangles.
    pub tris: Vec<[usize; 3]>,
    /// `nbr[t][i]`: the triangle across the edge opposite corner i.
    nbr: Vec<[usize; 3]>,
    pub alive: Vec<bool>,
    /// Region label of every triangle (see [`Cdt::label_regions`]); new
    /// triangles take the label of the cavity they fill.
    pub region: Vec<usize>,
    constrained: HashSet<(usize, usize)>,
    /// A live triangle at every vertex.
    at: Vec<usize>,
    last: usize,
    near: f64,
    /// The triangles the last insertion created.
    pub created: Vec<usize>,
}

/// The number of box corners at the start of `points`.
pub const CORNERS: usize = 4;

impl Cdt {
    /// Two triangles over the box around `min`/`max`, enlarged by its
    /// diagonal on every side.
    pub fn new(min: P2, max: P2) -> Self {
        let d = distance(min, max).max(1e-300);
        let (lo, hi) = ([min[0] - d, min[1] - d], [max[0] + d, max[1] + d]);
        Self {
            points: vec![lo, [hi[0], lo[1]], hi, [lo[0], hi[1]]],
            tris: vec![[0, 1, 2], [0, 2, 3]],
            nbr: vec![[NONE, 1, NONE], [NONE, NONE, 0]],
            alive: vec![true, true],
            region: vec![0, 0],
            constrained: HashSet::default(),
            at: vec![0, 0, 0, 1],
            last: 0,
            near: 1e-12 * d,
            created: Vec::new(),
        }
    }

    /// The edge opposite corner i of triangle t.
    pub fn edge(&self, t: usize, i: usize) -> [usize; 2] {
        let tri = self.tris[t];
        [tri[(i + 1) % 3], tri[(i + 2) % 3]]
    }

    pub fn neighbour(&self, t: usize, i: usize) -> usize {
        self.nbr[t][i]
    }

    pub fn is_constrained(&self, a: usize, b: usize) -> bool {
        self.constrained.contains(&key(a, b))
    }

    /// A triangle that contains `p` (closed), by a visibility walk.
    pub fn locate(&self, p: P2) -> Option<usize> {
        let mut t = if self.alive.get(self.last).copied().unwrap_or(false) {
            self.last
        } else {
            self.alive.iter().position(|&a| a)?
        };
        let limit = 4 * self.tris.len() + 64;
        'walk: for step in 0..limit {
            for k in 0..3 {
                let i = (step + k) % 3;
                let [a, b] = self.edge(t, i);
                if orient(self.points[a], self.points[b], p) == Ordering::Less {
                    let n = self.nbr[t][i];
                    if n == NONE {
                        return None;
                    }
                    t = n;
                    continue 'walk;
                }
            }
            return Some(t);
        }
        (0..self.tris.len()).find(|&t| {
            self.alive[t]
                && (0..3).all(|i| {
                    let [a, b] = self.edge(t, i);
                    orient(self.points[a], self.points[b], p) != Ordering::Less
                })
        })
    }

    fn near_vertex(&self, t: usize, p: P2) -> Option<usize> {
        self.tris[t]
            .iter()
            .copied()
            .find(|&v| distance(self.points[v], p) <= self.near)
    }

    /// Insert `p`; its index.
    pub fn insert(&mut self, p: P2) -> Result<usize, Refused> {
        if !(p[0].is_finite() && p[1].is_finite()) {
            return Err(Refused::Outside);
        }
        let start = self.locate(p).ok_or(Refused::Outside)?;
        if let Some(v) = self.near_vertex(start, p) {
            return Err(Refused::Duplicate(v));
        }
        let mut cavity: HashSet<usize> = HashSet::default();
        cavity.insert(start);
        let mut stack = vec![start];
        while let Some(t) = stack.pop() {
            for i in 0..3 {
                let n = self.nbr[t][i];
                if n == NONE || cavity.contains(&n) {
                    continue;
                }
                let [a, b] = self.edge(t, i);
                if self.is_constrained(a, b) {
                    continue;
                }
                let [x, y, z] = self.tris[n].map(|v| self.points[v]);
                if incircle(x, y, z, p) > 0.0 {
                    cavity.insert(n);
                    stack.push(n);
                }
            }
        }
        // grow until every boundary edge sees p strictly on its left
        loop {
            let mut grow = Vec::new();
            for &t in &cavity {
                for i in 0..3 {
                    let n = self.nbr[t][i];
                    if n != NONE && cavity.contains(&n) {
                        continue;
                    }
                    let [a, b] = self.edge(t, i);
                    if orient(self.points[a], self.points[b], p) != Ordering::Greater {
                        if n == NONE {
                            return Err(Refused::Outside);
                        }
                        if self.is_constrained(a, b) {
                            return Err(Refused::Constrained);
                        }
                        grow.push(n);
                    }
                }
            }
            if grow.is_empty() {
                break;
            }
            cavity.extend(grow);
        }
        let mut cells: Vec<usize> = cavity.iter().copied().collect();
        cells.sort_unstable();
        for &t in &cells {
            if let Some(v) = self.near_vertex(t, p) {
                return Err(Refused::Duplicate(v));
            }
        }
        let region = self.region[start];
        let index = self.points.len();
        self.points.push(p);
        self.at.push(NONE);
        let mut boundary = Vec::new();
        for &t in &cells {
            for i in 0..3 {
                let n = self.nbr[t][i];
                if n == NONE || !cavity.contains(&n) {
                    let [a, b] = self.edge(t, i);
                    boundary.push((a, b, n, t));
                }
            }
        }
        for &t in &cells {
            self.alive[t] = false;
        }
        self.created.clear();
        let mut starting: HashMap<usize, usize> = HashMap::default();
        let mut ending: HashMap<usize, usize> = HashMap::default();
        for &(a, b, outer, old) in &boundary {
            let t = self.tris.len();
            self.tris.push([a, b, index]);
            self.nbr.push([NONE, NONE, outer]);
            self.alive.push(true);
            self.region.push(region);
            if outer != NONE {
                for j in 0..3 {
                    if self.nbr[outer][j] == old {
                        let [x, y] = self.edge(outer, j);
                        if x == b && y == a {
                            self.nbr[outer][j] = t;
                        }
                    }
                }
            }
            starting.insert(a, t);
            ending.insert(b, t);
            self.created.push(t);
        }
        for &t in &self.created {
            let [a, b, _] = self.tris[t];
            // edge (b, p) is edge (p, b) of the triangle that starts at b
            self.nbr[t][0] = starting.get(&b).copied().unwrap_or(NONE);
            // edge (p, a) is edge (a, p) of the triangle that ends at a
            self.nbr[t][1] = ending.get(&a).copied().unwrap_or(NONE);
        }
        for &t in &self.created {
            for v in self.tris[t] {
                self.at[v] = t;
            }
        }
        self.last = self.created[0];
        Ok(index)
    }

    /// The live triangles around vertex v.
    pub fn around(&self, v: usize) -> Vec<usize> {
        let start = self.at[v];
        if start == NONE || !self.alive[start] {
            return (0..self.tris.len())
                .filter(|&t| self.alive[t] && self.tris[t].contains(&v))
                .collect();
        }
        let mut seen = vec![start];
        let mut i = 0;
        while i < seen.len() {
            let t = seen[i];
            for k in 0..3 {
                let n = self.nbr[t][k];
                if n != NONE && self.tris[n].contains(&v) && !seen.contains(&n) {
                    seen.push(n);
                }
            }
            i += 1;
        }
        seen
    }

    /// The triangle and corner opposite the directed edge a → b.
    fn find_edge(&self, a: usize, b: usize) -> Option<(usize, usize)> {
        for t in self.around(a) {
            for i in 0..3 {
                if self.edge(t, i) == [a, b] {
                    return Some((t, i));
                }
            }
        }
        None
    }

    pub fn has_edge(&self, a: usize, b: usize) -> bool {
        self.find_edge(a, b).is_some() || self.find_edge(b, a).is_some()
    }

    /// Flip the edge opposite corner i of t; false when the quadrilateral
    /// is not strictly convex.
    fn flip(&mut self, t: usize, i: usize) -> bool {
        let n = self.nbr[t][i];
        if n == NONE {
            return false;
        }
        let [a, b] = self.edge(t, i);
        let c = self.tris[t][i];
        let Some(j) = (0..3).find(|&k| self.tris[n][k] != a && self.tris[n][k] != b) else {
            return false;
        };
        let d = self.tris[n][j];
        let p = |v: usize| self.points[v];
        if orient(p(c), p(a), p(d)) != Ordering::Greater
            || orient(p(c), p(d), p(b)) != Ordering::Greater
        {
            return false;
        }
        let opposite = |s: &Self, tri: usize, v: usize| {
            let k = (0..3).find(|&k| s.tris[tri][k] == v).unwrap_or(0);
            s.nbr[tri][k]
        };
        let t_ob = opposite(self, t, b); // across (c, a)
        let t_oa = opposite(self, t, a); // across (b, c)
        let n_ob = opposite(self, n, b); // across (a, d)
        let n_oa = opposite(self, n, a); // across (d, b)
        self.tris[t] = [c, a, d];
        self.nbr[t] = [n_ob, n, t_ob];
        self.tris[n] = [c, d, b];
        self.nbr[n] = [n_oa, t_oa, t];
        let relink = |s: &mut Self, outer: usize, from: usize, to: usize| {
            if outer != NONE {
                for k in 0..3 {
                    if s.nbr[outer][k] == from {
                        s.nbr[outer][k] = to;
                    }
                }
            }
        };
        relink(self, n_ob, n, t);
        relink(self, t_oa, t, n);
        for v in [a, c, d] {
            self.at[v] = t;
        }
        self.at[b] = n;
        self.last = t;
        true
    }

    /// Make a → b an edge of the triangulation and mark it as a
    /// constraint. A vertex on the open segment splits it in two.
    pub fn constrain(&mut self, a: usize, b: usize) -> Result<(), String> {
        if a == b {
            return Err("degenerate constraint".into());
        }
        if self.has_edge(a, b) {
            self.constrained.insert(key(a, b));
            return Ok(());
        }
        let (pa, pb) = (self.points[a], self.points[b]);
        let side = |s: &Self, v: usize| orient(pa, pb, s.points[v]);
        // the first crossed edge: the corner at a whose wedge holds b
        let mut first = None;
        for t in self.around(a) {
            let k = (0..3).find(|&k| self.tris[t][k] == a).unwrap_or(0);
            let (u, v) = (self.tris[t][(k + 1) % 3], self.tris[t][(k + 2) % 3]);
            for w in [u, v] {
                if side(self, w) == Ordering::Equal && between(pa, pb, self.points[w]) {
                    self.constrain(a, w)?;
                    return self.constrain(w, b);
                }
            }
            if orient(pa, self.points[u], pb) == Ordering::Greater
                && orient(pa, pb, self.points[v]) == Ordering::Greater
            {
                first = Some((t, k, u, v));
                break;
            }
        }
        let Some((mut t, k, mut right, mut left)) = first else {
            return Err("constraint start not found".into());
        };
        let mut crossed = vec![(right, left)];
        let mut across = self.nbr[t][k];
        let mut guard = 0;
        loop {
            guard += 1;
            if across == NONE || guard > self.tris.len() + 8 {
                return Err("constraint leaves the triangulation".into());
            }
            t = across;
            let Some(w) = self.tris[t].into_iter().find(|&w| w != right && w != left) else {
                return Err("broken triangle".into());
            };
            if w == b {
                break;
            }
            match side(self, w) {
                Ordering::Equal => {
                    self.constrain(a, w)?;
                    return self.constrain(w, b);
                }
                Ordering::Greater => left = w,
                Ordering::Less => right = w,
            }
            crossed.push((right, left));
            // the next triangle: across the edge (right, left) of t
            across = (0..3)
                .find(|&i| {
                    let [x, y] = self.edge(t, i);
                    (x == right && y == left) || (x == left && y == right)
                })
                .map_or(NONE, |i| self.nbr[t][i]);
        }
        let mut queue: VecDeque<(usize, usize)> = crossed.into();
        let crosses = |s: &Self, u: usize, v: usize| {
            u != a
                && u != b
                && v != a
                && v != b
                && side(s, u) != side(s, v)
                && side(s, u) != Ordering::Equal
                && side(s, v) != Ordering::Equal
                && orient(s.points[u], s.points[v], pa) != orient(s.points[u], s.points[v], pb)
        };
        let mut guard = 0usize;
        while let Some((u, v)) = queue.pop_front() {
            guard += 1;
            if guard > 1_000_000 {
                return Err("constraint recovery did not converge".into());
            }
            let Some((t, i)) = self.find_edge(u, v).or_else(|| self.find_edge(v, u)) else {
                continue;
            };
            if !crosses(self, u, v) {
                continue;
            }
            if self.flip(t, i) {
                // the new diagonal
                let [c, _, d] = self.tris[t];
                if crosses(self, c, d) {
                    queue.push_back((c, d));
                }
            } else {
                queue.push_back((u, v));
            }
        }
        if !self.has_edge(a, b) {
            return Err("constraint not recovered".into());
        }
        self.constrained.insert(key(a, b));
        Ok(())
    }

    /// Lawson flips of every unconstrained edge that is not locally
    /// Delaunay.
    pub fn make_delaunay(&mut self) {
        let mut stack: Vec<(usize, usize)> = (0..self.tris.len())
            .filter(|&t| self.alive[t])
            .flat_map(|t| (0..3).map(move |i| (t, i)))
            .collect();
        let mut guard = 0usize;
        while let Some((t, i)) = stack.pop() {
            guard += 1;
            if guard > 10_000_000 {
                break;
            }
            if !self.alive[t] {
                continue;
            }
            let n = self.nbr[t][i];
            if n == NONE {
                continue;
            }
            let [a, b] = self.edge(t, i);
            if self.is_constrained(a, b) {
                continue;
            }
            let Some(d) = self.tris[n].into_iter().find(|&w| w != a && w != b) else {
                continue;
            };
            let [x, y, z] = self.tris[t].map(|v| self.points[v]);
            let w = self.points[d];
            // cocircular corners (a regular grid) flip only on a clear
            // in-circle test, or rounding would flip them back and forth
            let m = [x, y, z]
                .iter()
                .map(|q| (q[0] - w[0]).abs().max((q[1] - w[1]).abs()))
                .fold(0.0, f64::max);
            if incircle(x, y, z, w) > 1e-10 * m.powi(4) && self.flip(t, i) {
                for k in 0..3 {
                    stack.push((t, k));
                    stack.push((n, k));
                }
            }
        }
    }

    /// Label the regions between constraints: 0 for the one at the box
    /// corners, then 1, 2, ...; the number of labels.
    pub fn label_regions(&mut self) -> usize {
        let mut label = vec![NONE; self.tris.len()];
        let mut next = 0;
        let mut seeds: Vec<usize> = (0..self.tris.len())
            .filter(|&t| self.alive[t] && self.tris[t].iter().any(|&v| v < CORNERS))
            .collect();
        seeds.extend((0..self.tris.len()).filter(|&t| self.alive[t]));
        for seed in seeds {
            if label[seed] != NONE {
                continue;
            }
            label[seed] = next;
            let mut stack = vec![seed];
            while let Some(t) = stack.pop() {
                for i in 0..3 {
                    let n = self.nbr[t][i];
                    let [a, b] = self.edge(t, i);
                    if n == NONE || label[n] != NONE || self.is_constrained(a, b) {
                        continue;
                    }
                    label[n] = next;
                    stack.push(n);
                }
            }
            // the box corners form one region whatever the constraints
            if next > 0 || !self.tris[seed].iter().any(|&v| v < CORNERS) {
                next += 1;
            } else {
                next = 1;
            }
        }
        for (t, l) in label.into_iter().enumerate() {
            if self.alive[t] {
                self.region[t] = l;
            }
        }
        next
    }
}

/// Whether c lies strictly between a and b on their line (collinear c).
fn between(a: P2, b: P2, c: P2) -> bool {
    let d = [b[0] - a[0], b[1] - a[1]];
    let s = (c[0] - a[0]) * d[0] + (c[1] - a[1]) * d[1];
    s > 0.0 && s < d[0] * d[0] + d[1] * d[1]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check(cdt: &Cdt) {
        let mut edges: HashMap<(usize, usize), usize> = HashMap::default();
        for t in 0..cdt.tris.len() {
            if !cdt.alive[t] {
                continue;
            }
            let [a, b, c] = cdt.tris[t].map(|v| cdt.points[v]);
            assert_eq!(orient(a, b, c), Ordering::Greater, "triangle {t}");
            for i in 0..3 {
                let [u, v] = cdt.edge(t, i);
                assert!(edges.insert((u, v), t).is_none(), "edge twice");
                let n = cdt.nbr[t][i];
                if n != NONE {
                    assert!(cdt.alive[n]);
                    assert!(
                        (0..3).any(|j| cdt.edge(n, j) == [v, u] && cdt.nbr[n][j] == t),
                        "neighbour link"
                    );
                }
            }
        }
    }

    // Keep rule (SW-18): random points and a zigzag of constraints give a
    // valid, linked, counter-clockwise triangulation that contains every
    // constraint, and a square loop of constraints splits it into two
    // regions with the right areas. Breaks when an insertion folds a
    // triangle, a flip loses a link or recovery misses an edge; the shell
    // mesher depends on all of these and no other test isolates them.
    #[test]
    fn constraints_are_recovered_and_regions_split() {
        let mut cdt = Cdt::new([0.0, 0.0], [1.0, 1.0]);
        let mut rng = 0x1234_5678_9abc_def0u64;
        let mut random = || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            (rng >> 11) as f64 / (1u64 << 53) as f64
        };
        let square = [[0.2, 0.2], [0.8, 0.2], [0.8, 0.8], [0.2, 0.8]];
        let ids: Vec<usize> = square.iter().map(|&p| cdt.insert(p).unwrap()).collect();
        for _ in 0..400 {
            let p = [random(), random()];
            let _ = cdt.insert(p);
        }
        check(&cdt);
        for k in 0..4 {
            cdt.constrain(ids[k], ids[(k + 1) % 4]).unwrap();
        }
        // a long diagonal constraint across many triangles, collinear points
        let (a, b) = (
            cdt.insert([0.05, 0.5]).unwrap(),
            cdt.insert([0.15, 0.95]).unwrap(),
        );
        cdt.constrain(a, b).unwrap();
        cdt.make_delaunay();
        check(&cdt);
        for k in 0..4 {
            assert!(cdt.has_edge(ids[k], ids[(k + 1) % 4]));
        }
        assert!(cdt.has_edge(a, b));
        let n = cdt.label_regions();
        assert!(n >= 2, "{n}");
        let area = |r: usize| -> f64 {
            (0..cdt.tris.len())
                .filter(|&t| cdt.alive[t] && cdt.region[t] == r)
                .map(|t| {
                    let [p, q, s] = cdt.tris[t].map(|v| cdt.points[v]);
                    0.5 * ((q[0] - p[0]) * (s[1] - p[1]) - (q[1] - p[1]) * (s[0] - p[0]))
                })
                .sum()
        };
        let inner: f64 = (1..n).map(area).filter(|&x| (x - 0.36).abs() < 1e-12).sum();
        assert!((inner - 0.36).abs() < 1e-12, "{inner}");
        // a refinement point inside the square keeps the region label
        let t = cdt.locate([0.5, 0.5]).unwrap();
        let r = cdt.region[t];
        cdt.insert([0.5, 0.5]).unwrap();
        assert!(cdt.created.iter().all(|&t| cdt.region[t] == r));
        check(&cdt);
    }
}
