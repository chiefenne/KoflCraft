//! Nearest-site queries for Voronoi clipping (ALG-04): a k-d tree over the
//! sites (Bentley 1975; Friedman, Bentley & Finkel 1977), queried for the k
//! nearest sites of a site in order of distance, ties by site index, so the
//! result is a pure function of the sites.

use mesh_core::geometry::Vec3;

/// Leaves hold at most this many sites.
const LEAF: usize = 8;

/// A static k-d tree over points.
#[derive(Clone, Debug)]
pub struct KdTree {
    points: Vec<Vec3>,
    /// Point indices, permuted so every node owns a contiguous range.
    order: Vec<usize>,
    nodes: Vec<Node>,
}

#[derive(Clone, Copy, Debug)]
struct Node {
    start: usize,
    end: usize,
    /// Split axis and value for inner nodes; children follow in `nodes`.
    split: Option<(usize, f64, usize, usize)>,
    min: Vec3,
    max: Vec3,
}

fn axis(p: Vec3, a: usize) -> f64 {
    match a {
        0 => p.x,
        1 => p.y,
        _ => p.z,
    }
}

impl KdTree {
    pub fn new(points: &[Vec3]) -> Self {
        let mut tree = Self {
            points: points.to_vec(),
            order: (0..points.len()).collect(),
            nodes: Vec::with_capacity(2 * points.len() / LEAF + 1),
        };
        if !points.is_empty() {
            tree.build(0, points.len());
        }
        tree
    }

    fn build(&mut self, start: usize, end: usize) -> usize {
        let (mut min, mut max) = (
            self.points[self.order[start]],
            self.points[self.order[start]],
        );
        for &i in &self.order[start..end] {
            let p = self.points[i];
            min = Vec3::new(min.x.min(p.x), min.y.min(p.y), min.z.min(p.z));
            max = Vec3::new(max.x.max(p.x), max.y.max(p.y), max.z.max(p.z));
        }
        let id = self.nodes.len();
        self.nodes.push(Node {
            start,
            end,
            split: None,
            min,
            max,
        });
        if end - start > LEAF {
            let extent = max - min;
            let a = if extent.x >= extent.y && extent.x >= extent.z {
                0
            } else if extent.y >= extent.z {
                1
            } else {
                2
            };
            let mid = start + (end - start) / 2;
            let points = &self.points;
            self.order[start..end].select_nth_unstable_by(mid - start, |&i, &j| {
                axis(points[i], a)
                    .total_cmp(&axis(points[j], a))
                    .then(i.cmp(&j))
            });
            let value = axis(self.points[self.order[mid]], a);
            let left = self.build(start, mid);
            let right = self.build(mid, end);
            self.nodes[id].split = Some((a, value, left, right));
        }
        id
    }

    /// The `k` points nearest to point `of` (itself excluded), as
    /// `(squared distance, index)` sorted by distance and then index.
    pub fn nearest(&self, of: usize, k: usize) -> Vec<(f64, usize)> {
        match self.points.get(of) {
            Some(&q) => self.search(q, Some(of), k),
            None => Vec::new(),
        }
    }

    /// The `k` points nearest to `q`, sorted like [`KdTree::nearest`].
    pub fn nearest_to(&self, q: Vec3, k: usize) -> Vec<(f64, usize)> {
        self.search(q, None, k)
    }

    fn search(&self, q: Vec3, skip: Option<usize>, k: usize) -> Vec<(f64, usize)> {
        let mut best: Vec<(f64, usize)> = Vec::with_capacity(k + 1);
        if k == 0 || self.nodes.is_empty() {
            return best;
        }
        let of = skip.unwrap_or(usize::MAX);
        let mut stack = vec![0usize];
        while let Some(n) = stack.pop() {
            let node = self.nodes[n];
            if best.len() == k && box_distance(q, node.min, node.max) > best[k - 1].0 {
                continue;
            }
            match node.split {
                None => {
                    for &i in &self.order[node.start..node.end] {
                        if i == of {
                            continue;
                        }
                        let d = (self.points[i] - q).dot(self.points[i] - q);
                        let key = (d, i);
                        if best.len() < k || less(key, best[k - 1]) {
                            let at = best.partition_point(|&b| less(b, key));
                            best.insert(at, key);
                            best.truncate(k);
                        }
                    }
                }
                Some((a, value, left, right)) => {
                    // visit the near side first (pushed last)
                    if axis(q, a) < value {
                        stack.push(right);
                        stack.push(left);
                    } else {
                        stack.push(left);
                        stack.push(right);
                    }
                }
            }
        }
        best
    }
}

impl KdTree {
    /// The point nearest to `q`, ties by index; `None` for an empty tree.
    pub fn closest(&self, q: Vec3) -> Option<usize> {
        let mut best: Option<(f64, usize)> = None;
        if self.nodes.is_empty() {
            return None;
        }
        let mut stack = vec![0usize];
        while let Some(n) = stack.pop() {
            let node = self.nodes[n];
            if best.is_some_and(|b| box_distance(q, node.min, node.max) > b.0) {
                continue;
            }
            match node.split {
                None => {
                    for &i in &self.order[node.start..node.end] {
                        let d = (self.points[i] - q).dot(self.points[i] - q);
                        if best.is_none_or(|b| less((d, i), b)) {
                            best = Some((d, i));
                        }
                    }
                }
                Some((a, value, left, right)) => {
                    if axis(q, a) < value {
                        stack.push(right);
                        stack.push(left);
                    } else {
                        stack.push(left);
                        stack.push(right);
                    }
                }
            }
        }
        best.map(|b| b.1)
    }
}

fn less(a: (f64, usize), b: (f64, usize)) -> bool {
    a.0 < b.0 || (a.0 == b.0 && a.1 < b.1)
}

/// Squared distance from `q` to the box `[min, max]`.
fn box_distance(q: Vec3, min: Vec3, max: Vec3) -> f64 {
    let d = |v: f64, lo: f64, hi: f64| {
        if v < lo {
            lo - v
        } else if v > hi {
            v - hi
        } else {
            0.0
        }
    };
    let (x, y, z) = (
        d(q.x, min.x, max.x),
        d(q.y, min.y, max.y),
        d(q.z, min.z, max.z),
    );
    x * x + y * y + z * z
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Invariant: the tree returns exactly the k nearest of a brute-force
    /// sort, ties by index, on a lattice full of equal distances and on
    /// scattered points. A pruning bound that is too tight drops a
    /// neighbour and fails here before any cell is clipped wrongly.
    #[test]
    fn nearest_equals_brute_force() {
        let mut points = Vec::new();
        for i in 0..6 {
            for j in 0..5 {
                for k in 0..4 {
                    points.push(Vec3::new(i as f64, j as f64, k as f64));
                }
            }
        }
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        for _ in 0..300 {
            let mut r = || {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                (state >> 11) as f64 / (1u64 << 53) as f64
            };
            points.push(Vec3::new(r() * 7.0, r() * 3.0, r() * 0.5));
        }
        let tree = KdTree::new(&points);
        for of in (0..points.len()).step_by(7) {
            let mut all: Vec<(f64, usize)> = (0..points.len())
                .filter(|&i| i != of)
                .map(|i| ((points[i] - points[of]).dot(points[i] - points[of]), i))
                .collect();
            all.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            for k in [1, 5, 26, 64, points.len()] {
                let want: Vec<_> = all.iter().take(k).copied().collect();
                assert_eq!(tree.nearest(of, k), want, "point {of}, k {k}");
            }
            // A query point equal to a lattice point finds it (distance 0).
            assert_eq!(
                tree.closest(points[of]),
                Some(
                    of.min(
                        (0..points.len())
                            .find(|&i| points[i] == points[of])
                            .unwrap_or(of)
                    )
                )
            );
            let q = points[of] + Vec3::new(0.31, -0.17, 0.05);
            let want = (0..points.len())
                .map(|i| ((points[i] - q).dot(points[i] - q), i))
                .min_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)))
                .map(|b| b.1);
            assert_eq!(tree.closest(q), want);
        }
    }
}
