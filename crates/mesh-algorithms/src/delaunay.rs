//! Incremental 3D Delaunay tetrahedralization with constrained faces.
//!
//! Points are inserted by the Bowyer–Watson cavity method (Bowyer 1981;
//! Watson 1981): the tetrahedra whose circumspheres contain the new point
//! are removed and the cavity is filled by joining its boundary faces to
//! the point. The cavity grows until every boundary face sees the point
//! strictly from the front, so new tetrahedra are positively oriented even
//! for cospherical input. Points are located by a visibility walk. A set of
//! constrained faces is never crossed by a cavity; a point whose cavity
//! cannot be made star-shaped without crossing one is refused, which keeps
//! a recovered boundary in place during interior refinement. The signs come
//! from the exact predicates of `predicates`.
//!
//! Eight corners of an enlarged bounding box start the triangulation, so
//! every point inserted inside the box has a containing tetrahedron.

#![allow(clippy::needless_range_loop)] // parallel per-vertex and per-face arrays

use crate::predicates::{insphere, orient3d};
use mesh_core::geometry::Vec3;
use std::cmp::Ordering;
use std::collections::hash_map::DefaultHasher;
use std::hash::BuildHasherDefault;

/// Hash containers with a fixed hasher: insertion and walk order depend on
/// their iteration order, so a per-process random seed would make the
/// scaffold differ from run to run.
pub(crate) type HashMap<K, V> =
    std::collections::HashMap<K, V, BuildHasherDefault<DefaultHasher>>;
/// See [`HashMap`].
pub(crate) type HashSet<K> = std::collections::HashSet<K, BuildHasherDefault<DefaultHasher>>;

pub const NONE: usize = usize::MAX;

/// Faces opposite each local vertex, ordered so that the face with the
/// opposite vertex forms a positive tetrahedron: (f0, f1, f2, opposite).
pub const FACES: [[usize; 3]; 4] = [[1, 3, 2], [0, 2, 3], [0, 3, 1], [0, 1, 2]];

pub struct Delaunay {
    pub points: Vec<Vec3>,
    pub tets: Vec<[usize; 4]>,
    pub neighbours: Vec<[usize; 4]>,
    pub alive: Vec<bool>,
    free: Vec<usize>,
    last: usize,
    pub constrained: HashSet<[usize; 3]>,
    /// The number of bounding-box corner points at the start of `points`.
    pub corners: usize,
    /// The tetrahedra the last successful insertion created.
    pub created: Vec<usize>,
    rng: u64,
}

/// Why an insertion did not happen.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Refused {
    /// The point is already a vertex (its index).
    Duplicate(usize),
    Outside,
    Constrained,
}

fn sorted(mut f: [usize; 3]) -> [usize; 3] {
    f.sort_unstable();
    f
}

impl Delaunay {
    /// A triangulation of the eight corners of the box around `min`/`max`,
    /// enlarged by its diagonal on every side.
    pub fn new(min: Vec3, max: Vec3) -> Self {
        let d = (max - min).length().max(1.0);
        let lo = min - Vec3::new(d, d, d);
        let hi = max + Vec3::new(d, d, d);
        let mut points = Vec::new();
        for i in 0..8 {
            points.push(Vec3::new(
                if i & 1 == 0 { lo.x } else { hi.x },
                if i & 2 == 0 { lo.y } else { hi.y },
                if i & 4 == 0 { lo.z } else { hi.z },
            ));
        }
        // Kuhn subdivision of the cube into six tetrahedra along 0–7.
        let paths = [[1, 3], [1, 5], [2, 3], [2, 6], [4, 5], [4, 6]];
        let mut tets = Vec::new();
        for [a, b] in paths {
            let mut t = [0, a, b, 7];
            if orient3d(points[t[0]], points[t[1]], points[t[2]], points[t[3]]) == Ordering::Less {
                t.swap(1, 2);
            }
            tets.push(t);
        }
        let mut dt = Self {
            points,
            neighbours: vec![[NONE; 4]; tets.len()],
            alive: vec![true; tets.len()],
            tets,
            free: Vec::new(),
            last: 0,
            constrained: HashSet::default(),
            corners: 8,
            created: Vec::new(),
            rng: 0x9E37_79B9_7F4A_7C15,
        };
        dt.link_all();
        dt
    }

    fn link_all(&mut self) {
        let mut faces: HashMap<[usize; 3], (usize, usize)> = HashMap::default();
        for t in 0..self.tets.len() {
            if !self.alive[t] {
                continue;
            }
            for i in 0..4 {
                let key = sorted(FACES[i].map(|k| self.tets[t][k]));
                if let Some((u, j)) = faces.remove(&key) {
                    self.neighbours[t][i] = u;
                    self.neighbours[u][j] = t;
                } else {
                    faces.insert(key, (t, i));
                }
            }
        }
    }

    fn next_random(&mut self) -> u64 {
        self.rng ^= self.rng << 13;
        self.rng ^= self.rng >> 7;
        self.rng ^= self.rng << 17;
        self.rng
    }

    fn pos(&self, t: usize) -> [Vec3; 4] {
        self.tets[t].map(|v| self.points[v])
    }

    /// A tetrahedron containing `p` (or `None` outside the box).
    pub fn locate(&mut self, p: Vec3) -> Option<usize> {
        let mut t = if self.last < self.tets.len() && self.alive[self.last] {
            self.last
        } else {
            self.alive.iter().position(|&a| a)?
        };
        let limit = 4 * self.tets.len() + 64;
        'walk: for _ in 0..limit {
            let start = (self.next_random() % 4) as usize;
            let q = self.pos(t);
            for k in 0..4 {
                let i = (start + k) % 4;
                let f = FACES[i];
                if orient3d(q[f[0]], q[f[1]], q[f[2]], p) == Ordering::Less {
                    // p is beyond face i (the positive side faces inward)
                    let n = self.neighbours[t][i];
                    if n == NONE {
                        return None;
                    }
                    t = n;
                    continue 'walk;
                }
            }
            self.last = t;
            return Some(t);
        }
        None
    }

    /// Insert `p` and return its index.
    pub fn insert(&mut self, p: Vec3) -> Result<usize, Refused> {
        let start = self.locate(p).ok_or(Refused::Outside)?;
        if let Some(&v) = self.tets[start].iter().find(|&&v| self.points[v] == p) {
            return Err(Refused::Duplicate(v));
        }
        let mut in_cavity: HashSet<usize> = HashSet::default();
        let mut stack = vec![start];
        in_cavity.insert(start);
        while let Some(t) = stack.pop() {
            for i in 0..4 {
                let n = self.neighbours[t][i];
                if n == NONE || in_cavity.contains(&n) {
                    continue;
                }
                let face = sorted(FACES[i].map(|k| self.tets[t][k]));
                if self.constrained.contains(&face) {
                    continue;
                }
                let q = self.pos(n);
                if insphere(q[0], q[1], q[2], q[3], p) == Ordering::Greater {
                    in_cavity.insert(n);
                    stack.push(n);
                }
            }
        }
        // make the cavity star-shaped from p
        loop {
            let mut grow = Vec::new();
            for &t in &in_cavity {
                let q = self.pos(t);
                for i in 0..4 {
                    let n = self.neighbours[t][i];
                    if n != NONE && in_cavity.contains(&n) {
                        continue;
                    }
                    let f = FACES[i];
                    // p must be strictly on the inner side of the face
                    if orient3d(q[f[0]], q[f[1]], q[f[2]], p) != Ordering::Greater {
                        let face = sorted(f.map(|k| self.tets[t][k]));
                        if n == NONE {
                            return Err(Refused::Outside);
                        }
                        if self.constrained.contains(&face) {
                            return Err(Refused::Constrained);
                        }
                        grow.push(n);
                    }
                }
            }
            if grow.is_empty() {
                break;
            }
            in_cavity.extend(grow);
        }
        let index = self.points.len();
        self.points.push(p);
        // boundary faces: (old tet, local face)
        let mut boundary = Vec::new();
        for &t in &in_cavity {
            for i in 0..4 {
                let n = self.neighbours[t][i];
                if n == NONE || !in_cavity.contains(&n) {
                    boundary.push((t, i, n));
                }
            }
        }
        let mut created = Vec::with_capacity(boundary.len());
        let mut edges: HashMap<(usize, usize), (usize, usize)> = HashMap::default();
        for &(t, i, n) in &boundary {
            let f = FACES[i].map(|k| self.tets[t][k]);
            let tet = [f[0], f[1], f[2], index];
            let id = if let Some(slot) = self.free.pop() {
                self.tets[slot] = tet;
                self.neighbours[slot] = [NONE; 4];
                self.alive[slot] = true;
                slot
            } else {
                self.tets.push(tet);
                self.neighbours.push([NONE; 4]);
                self.alive.push(true);
                self.tets.len() - 1
            };
            created.push(id);
            // the face opposite p (local 3) is the old boundary face
            self.neighbours[id][3] = n;
            if n != NONE {
                let back = (0..4).find(|&j| self.neighbours[n][j] == t);
                if let Some(j) = back {
                    self.neighbours[n][j] = id;
                }
            }
            // faces opposite local 0, 1, 2 contain p and an edge of f
            for local in 0..3 {
                let others: Vec<usize> = (0..3).filter(|&k| k != local).map(|k| f[k]).collect();
                let e = (others[0].min(others[1]), others[0].max(others[1]));
                if let Some((u, j)) = edges.remove(&e) {
                    self.neighbours[id][local] = u;
                    self.neighbours[u][j] = id;
                } else {
                    edges.insert(e, (id, local));
                }
            }
        }
        for t in in_cavity {
            if !created.contains(&t) {
                self.alive[t] = false;
                self.free.push(t);
            }
        }
        if let Some(&t) = created.first() {
            self.last = t;
        }
        self.created = created;
        Ok(index)
    }

    /// Replace the tetrahedra around the interior edge (c, d), found from
    /// `start`, by a fan triangulation of their ring of vertices joined to c
    /// and d (edge removal, Shewchuk 2002, *Two Discrete Optimization
    /// Algorithms for the Topological Improvement of Tetrahedral Meshes*),
    /// choosing the fan apex whose worst `quality` is best. It happens only
    /// when every new tetrahedron is positive and the worst quality rises;
    /// no constrained face is removed (the ring faces through c–d must be
    /// unconstrained; the outer faces stay). Returns the new tetrahedra.
    pub fn remove_edge(
        &mut self,
        start: usize,
        c: usize,
        d: usize,
        quality: impl Fn([Vec3; 4]) -> f64,
    ) -> Option<Vec<usize>> {
        // walk the ring: in tet t, the other vertices u, v; go across the
        // face opposite u (it holds c, d, v)
        let mut ring = Vec::new();
        let mut tets = Vec::new();
        let mut t = start;
        let first = self.tets[t];
        if !first.contains(&c) || !first.contains(&d) {
            return None;
        }
        let others = |tet: [usize; 4]| {
            let o: Vec<usize> = tet.into_iter().filter(|&v| v != c && v != d).collect();
            (o[0], o[1])
        };
        let (u0, _) = others(first);
        ring.push(u0);
        let mut cur = u0;
        loop {
            tets.push(t);
            if tets.len() > 32 {
                return None;
            }
            let tet = self.tets[t];
            let (a, b) = others(tet);
            let next = if a == cur { b } else { a };
            let iu = tet.iter().position(|&x| x == cur)?;
            if self
                .constrained
                .contains(&sorted(FACES[iu].map(|k| tet[k])))
            {
                return None;
            }
            let n = self.neighbours[t][iu];
            if n == NONE || !self.alive[n] {
                return None;
            }
            if next == u0 {
                if n != start {
                    return None;
                }
                break;
            }
            ring.push(next);
            cur = next;
            t = n;
        }
        if ring.len() < 3 || ring.len() != tets.len() {
            return None;
        }
        let old_worst = tets
            .iter()
            .map(|&t| quality(self.pos(t)))
            .fold(f64::INFINITY, f64::min);
        let k = ring.len();
        let vol = |q: [Vec3; 4]| (q[1] - q[0]).dot((q[2] - q[0]).cross(q[3] - q[0]));
        let before: f64 = tets.iter().map(|&t| vol(self.pos(t))).sum();
        let mut best: Option<(f64, Vec<[usize; 4]>)> = None;
        // Optimal triangulation of the ring polygon by dynamic programming
        // over sub-polygons (Klincsek 1980), maximising the worst quality,
        // for each orientation sign of the ring.
        for sigma in [Ordering::Greater, Ordering::Less] {
            let tri = |i: usize, m: usize, j: usize| -> Option<(f64, [[usize; 4]; 2])> {
                let (a, b, e) = (ring[i], ring[m], ring[j]);
                let up = [a, b, e, c];
                let down = [b, a, e, d];
                let qu = up.map(|x| self.points[x]);
                let qd = down.map(|x| self.points[x]);
                if orient3d(qu[0], qu[1], qu[2], qu[3]) != sigma
                    || orient3d(qd[0], qd[1], qd[2], qd[3]) != sigma
                {
                    return None;
                }
                let fix = |mut t: [usize; 4]| {
                    if sigma == Ordering::Less {
                        t.swap(0, 1);
                    }
                    t
                };
                let (fu, fd) = (fix(up), fix(down));
                let q =
                    quality(fu.map(|x| self.points[x])).min(quality(fd.map(|x| self.points[x])));
                Some((q, [fu, fd]))
            };
            let mut score = vec![vec![f64::INFINITY; k]; k];
            let mut split = vec![vec![usize::MAX; k]; k];
            for len in 2..k {
                for i in 0..k - len {
                    let j = i + len;
                    let mut best_v = f64::NEG_INFINITY;
                    for m in i + 1..j {
                        let Some((q, _)) = tri(i, m, j) else { continue };
                        let v = q.min(score[i][m]).min(score[m][j]);
                        if v > best_v {
                            best_v = v;
                            split[i][j] = m;
                        }
                    }
                    score[i][j] = best_v;
                }
            }
            if score[0][k - 1] <= old_worst || !score[0][k - 1].is_finite() {
                continue;
            }
            let mut new = Vec::new();
            let mut stack = vec![(0, k - 1)];
            let mut ok = true;
            while let Some((i, j)) = stack.pop() {
                if j - i < 2 {
                    continue;
                }
                let m = split[i][j];
                match (m != usize::MAX).then(|| tri(i, m, j)).flatten() {
                    Some((_, pair)) => new.extend(pair),
                    None => ok = false,
                }
                stack.push((i, m));
                stack.push((m, j));
            }
            if ok && best.as_ref().is_none_or(|(w, _)| score[0][k - 1] > *w) {
                best = Some((score[0][k - 1], new));
            }
        }
        let (_, new) = best?;
        // the new tetrahedra must fill the same volume
        let after: f64 = new.iter().map(|t| vol(t.map(|x| self.points[x]))).sum();
        if (after - before).abs() > 1e-9 * before.abs() {
            return None;
        }
        Some(self.replace(&tets, &new))
    }

    /// Replace live tetrahedra `old` by `new` filling the same region.
    fn replace(&mut self, old: &[usize], new: &[[usize; 4]]) -> Vec<usize> {
        let old_set: HashSet<usize> = old.iter().copied().collect();
        let mut outer: HashMap<[usize; 3], usize> = HashMap::default();
        for &t in old {
            for i in 0..4 {
                let n = self.neighbours[t][i];
                if n == NONE || !old_set.contains(&n) {
                    outer.insert(sorted(FACES[i].map(|k| self.tets[t][k])), n);
                }
            }
        }
        for &t in old {
            self.alive[t] = false;
            self.free.push(t);
        }
        let mut ids = Vec::new();
        let mut inner: HashMap<[usize; 3], (usize, usize)> = HashMap::default();
        for &tet in new {
            let id = if let Some(slot) = self.free.pop() {
                self.tets[slot] = tet;
                self.neighbours[slot] = [NONE; 4];
                self.alive[slot] = true;
                slot
            } else {
                self.tets.push(tet);
                self.neighbours.push([NONE; 4]);
                self.alive.push(true);
                self.tets.len() - 1
            };
            ids.push(id);
            for i in 0..4 {
                let key = sorted(FACES[i].map(|k| tet[k]));
                if let Some(&n) = outer.get(&key) {
                    self.neighbours[id][i] = n;
                    if n != NONE {
                        for j in 0..4 {
                            let back = sorted(FACES[j].map(|k| self.tets[n][k]));
                            if back == key {
                                self.neighbours[n][j] = id;
                            }
                        }
                    }
                } else if let Some((u, j)) = inner.remove(&key) {
                    self.neighbours[id][i] = u;
                    self.neighbours[u][j] = id;
                } else {
                    inner.insert(key, (id, i));
                }
            }
        }
        self.last = ids[0];
        ids
    }

    /// The faces of the live tetrahedra, sorted, each once.
    pub fn face_set(&self) -> HashSet<[usize; 3]> {
        let mut out = HashSet::default();
        for (t, tet) in self.tets.iter().enumerate() {
            if self.alive[t] {
                for f in FACES {
                    out.insert(sorted(f.map(|k| tet[k])));
                }
            }
        }
        out
    }

    pub fn edge_set(&self) -> HashSet<(usize, usize)> {
        let mut out = HashSet::default();
        for (t, tet) in self.tets.iter().enumerate() {
            if self.alive[t] {
                for i in 0..4 {
                    for j in i + 1..4 {
                        out.insert((tet[i].min(tet[j]), tet[i].max(tet[j])));
                    }
                }
            }
        }
        out
    }
}

/// Circumcentre and circumradius of a tetrahedron.
pub fn circumsphere(p: [Vec3; 4]) -> Option<(Vec3, f64)> {
    let (a, b, c) = (p[1] - p[0], p[2] - p[0], p[3] - p[0]);
    let det = 2.0 * a.dot(b.cross(c));
    if det == 0.0 || !det.is_finite() {
        return None;
    }
    let num = b.cross(c) * a.dot(a) + c.cross(a) * b.dot(b) + a.cross(b) * c.dot(c);
    let o = num / det;
    Some((p[0] + o, o.length()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check(dt: &Delaunay) {
        for (t, tet) in dt.tets.iter().enumerate() {
            if !dt.alive[t] {
                continue;
            }
            let q = tet.map(|v| dt.points[v]);
            assert_eq!(orient3d(q[0], q[1], q[2], q[3]), Ordering::Greater);
            for i in 0..4 {
                let n = dt.neighbours[t][i];
                if n != NONE {
                    assert!(dt.alive[n]);
                    assert!(dt.neighbours[n].contains(&t));
                }
            }
        }
    }

    /// Invariant: every tetrahedron is positive, adjacency is symmetric and
    /// no input point lies strictly inside any circumsphere (empty sphere).
    /// Grid points are massively cospherical, the case where a cavity
    /// without the star-shape repair creates flat tetrahedra.
    #[test]
    fn grid_is_delaunay() {
        let mut dt = Delaunay::new(Vec3::ZERO, Vec3::new(3., 3., 3.));
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    dt.insert(Vec3::new(i as f64, j as f64, k as f64)).unwrap();
                }
            }
        }
        check(&dt);
        for (t, tet) in dt.tets.iter().enumerate() {
            if !dt.alive[t] {
                continue;
            }
            let q = tet.map(|v| dt.points[v]);
            for p in &dt.points[dt.corners..] {
                assert_ne!(insphere(q[0], q[1], q[2], q[3], *p), Ordering::Greater);
            }
        }
        let volume: f64 = (0..dt.tets.len())
            .filter(|&t| dt.alive[t])
            .map(|t| {
                let q = dt.pos(t);
                (q[1] - q[0]).dot((q[2] - q[0]).cross(q[3] - q[0])) / 6.0
            })
            .sum();
        let d = (Vec3::new(3., 3., 3.)).length();
        let side = 3.0 + 2.0 * d;
        assert!((volume / side.powi(3) - 1.0).abs() < 1e-12);
    }
}
