//! Conforming wall seeds for the Voronoi core under wall layers (M4).
//!
//! Every vertex of the innermost layer surface carries a ball. The seed of a
//! wall triangle is the point on the fluid side that lies on the three balls
//! of its corners. Two seeds over triangles that share an edge are both at
//! the ball radius from each end of the edge, so their bisector holds the
//! edge; if no other site lies inside a ball, the Voronoi cell of a seed
//! cuts the wall in exactly its triangle (for a point of the triangle the
//! seed is at most a ball radius away at the corners, any other site
//! farther, and a bisector half-space that holds no corner holds no point
//! of the triangle). The core cell over a wall triangle then meets the
//! innermost layer prism in the whole triangle: no cut wall polygons, no
//! sliver columns.
//!
//! This is the seed placement of VoroCrust (Abdelkader, Bajaj, Ebeida,
//! Mahmoud, Mitchell, Owens and Rushdi, "VoroCrust: Voronoi meshing without
//! clipping", ACM Trans. Graph. 39(3), 2020, and "Sampling conditions for
//! conforming Voronoi meshing by the VoroCrust algorithm", SoCG 2018), used
//! on one side of a given triangulation instead of on a fresh sampling. A
//! seed must lie outside every other ball: the triangles for which that
//! holds are those of the weighted Delaunay (regular) triangulation of the
//! balls (Edelsbrunner and Shah, "Incremental topological flipping works for
//! regular triangulations", Algorithmica 15, 1996), reached from the wall
//! triangulation by edge flips with the empty-ball test in space. With
//! balls that cover the wall every regular triangle has a seed. Built from
//! the papers alone (LIC-13); research/m4-wall-columns.md has the benchmark.

use crate::neighbours::KdTree;
use mesh_core::geometry::Vec3;
use std::collections::HashMap;

/// Ball radius of a vertex over the largest circumradius of its triangles.
/// Just above one the balls cover the wall and the seeds stay low, which
/// keeps them on their own side of a convex wall edge.
const COVER: f64 = 1.05;
/// Largest drop of the radius along a wall edge, per edge length.
const RADIUS_RATE: f64 = 0.3;
/// Relative spread put on the squared radii, so four corners of two
/// triangles are never exactly on one sphere (two seeds in one point).
const JITTER: f64 = 0.01;
/// A seed counts as inside a ball below this share of the squared radius.
const INSIDE: f64 = 1e-9;
/// Lowest seed over its triangle, in triangle sizes.
const LOWEST: f64 = 0.03;
/// Seeds of two triangles on one edge must be this share of the edge apart;
/// below it the direction of their bisector is lost to rounding.
const APART: f64 = 2e-3;
/// Share by which the ball opposite a pair of close seeds shrinks.
const SHRINK: f64 = 0.03;
/// Rounds of shrinking and flipping for close seeds.
const ROUNDS: usize = 12;
/// Free sites keep this share of the radius away from a ball centre. Well
/// above one, so a free site next to a wall vertex whose triangles are not
/// all seeded leaves the seeds there a fair part of the unseeded triangles,
/// not a slice, and no wall cell is flattened by a site just over it.
pub(super) const CLEARANCE: f64 = 1.5;

const NONE: usize = usize::MAX;

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct SeedStats {
    /// Seeds kept.
    pub seeds: usize,
    /// Edge flips to the regular triangulation of the balls.
    pub flips: usize,
    /// Triangles without a seed (three balls without a common point, a seed
    /// inside a fourth ball, or one of two seeds in nearly one point).
    pub dropped: usize,
    /// Balls shrunk to part two close seeds.
    pub shrunk: usize,
}

pub(super) struct WallSeeds {
    /// Seed and size of its triangle.
    pub seeds: Vec<(Vec3, f64)>,
    /// The wall vertices whose balls each seed lies on (its triangle after
    /// the flips), parallel to `seeds`.
    pub corners: Vec<[usize; 3]>,
    /// Ball radius by wall vertex; zero where no kept seed uses the vertex.
    pub radius: Vec<f64>,
    pub stats: SeedStats,
}

fn key(u: usize, v: usize) -> (usize, usize) {
    (u.min(v), u.max(v))
}

fn opposite(t: [usize; 3], e: (usize, usize)) -> usize {
    t.into_iter()
        .find(|&v| v != e.0 && v != e.1)
        .unwrap_or(t[0])
}

struct Balls<'a> {
    points: &'a [Vec3],
    /// Squared radii.
    weight: Vec<f64>,
}

impl Balls<'_> {
    fn normal(&self, t: [usize; 3]) -> Vec3 {
        let [a, b, c] = t.map(|v| self.points[v]);
        (b - a).cross(c - a)
    }

    fn size(&self, t: [usize; 3]) -> f64 {
        (0.5 * self.normal(t).length() * 4.0 / 3f64.sqrt()).sqrt()
    }

    /// The point on the three balls of `t` on the side its normal points
    /// away from, if it is high enough over the triangle.
    fn seed(&self, t: [usize; 3]) -> Option<Vec3> {
        let [a, b, c] = t.map(|v| self.points[v]);
        let [wa, wb, wc] = t.map(|v| self.weight[v]);
        let (u, v) = (b - a, c - a);
        let n = u.cross(v);
        let area2 = n.length();
        let (uu, uv, vv) = (u.dot(u), u.dot(v), v.dot(v));
        let det = uu * vv - uv * uv;
        if area2.is_nan() || area2 <= 0.0 || det <= 0.0 {
            return None;
        }
        // foot: equal power to the three balls, in the triangle plane
        let r1 = (uu - (wb - wa)) / 2.0;
        let r2 = (vv - (wc - wa)) / 2.0;
        let d = u * ((r1 * vv - r2 * uv) / det) + v * ((r2 * uu - r1 * uv) / det);
        let h2 = wa - d.dot(d);
        let low = LOWEST * self.size(t);
        (h2 >= low * low).then(|| a + d - n * (h2.sqrt() / area2))
    }

    /// Whether the seed of `t` is missing or inside the ball of `x`.
    fn covered(&self, t: [usize; 3], x: usize) -> bool {
        match self.seed(t) {
            None => true,
            Some(s) => {
                let d = s - self.points[x];
                d.dot(d) < self.weight[x] * (1.0 - INSIDE)
            }
        }
    }
}

struct Flips {
    tris: Vec<[usize; 3]>,
    sides: HashMap<(usize, usize), [usize; 2]>,
}

impl Flips {
    fn link(&mut self, i: usize) {
        let t = self.tris[i];
        for k in 0..3 {
            let slot = self
                .sides
                .entry(key(t[k], t[(k + 1) % 3]))
                .or_insert([NONE; 2]);
            if slot[0] == NONE {
                slot[0] = i;
            } else {
                slot[1] = i;
            }
        }
    }

    fn unlink(&mut self, i: usize) {
        let t = self.tris[i];
        for k in 0..3 {
            let e = key(t[k], t[(k + 1) % 3]);
            if let Some(slot) = self.sides.get_mut(&e) {
                if slot[0] == i {
                    *slot = [slot[1], NONE];
                } else if slot[1] == i {
                    slot[1] = NONE;
                }
                if slot[0] == NONE {
                    self.sides.remove(&e);
                }
            }
        }
    }

    /// Flip edges whose triangles have a seed inside the ball across the
    /// edge until the triangulation is regular for the balls; returns the
    /// number of flips.
    fn regular(&mut self, balls: &Balls, mut stack: Vec<(usize, usize)>) -> usize {
        let mut flips = 0;
        let mut budget = 30 * self.sides.len().max(stack.len());
        while let Some(e) = stack.pop() {
            if budget == 0 {
                break;
            }
            budget -= 1;
            let Some(&[a, b]) = self.sides.get(&e) else {
                continue;
            };
            if a == NONE || b == NONE {
                continue;
            }
            let (ta, tb) = (self.tris[a], self.tris[b]);
            let (k, l) = (opposite(ta, e), opposite(tb, e));
            if !(balls.covered(ta, l) || balls.covered(tb, k)) {
                continue;
            }
            let Some(q) = ta.iter().position(|&v| v == k) else {
                continue;
            };
            let (i, j) = (ta[(q + 1) % 3], ta[(q + 2) % 3]);
            let (n1, n2) = ([i, l, k], [l, j, k]);
            if self.sides.contains_key(&key(k, l)) {
                continue;
            }
            // the four corners must be a convex quadrilateral seen from
            // the wall side, and the new seeds outside the corners they lose
            let outward = balls.normal(ta) + balls.normal(tb);
            if balls.normal(n1).dot(outward) <= 0.0
                || balls.normal(n2).dot(outward) <= 0.0
                || balls.covered(n1, j)
                || balls.covered(n2, i)
            {
                continue;
            }
            self.unlink(a);
            self.unlink(b);
            self.tris[a] = n1;
            self.tris[b] = n2;
            self.link(a);
            self.link(b);
            flips += 1;
            stack.extend([key(i, l), key(l, j), key(j, k), key(k, i)]);
        }
        flips
    }
}

/// Seeds for the wall triangles `triangles` over `points` (normals out of
/// the fluid), and the balls free sites must stay out of.
pub(super) fn wall_seeds(points: &[Vec3], triangles: &[[usize; 3]]) -> WallSeeds {
    let count = points.len();
    let mut stats = SeedStats::default();
    // Radii: the balls of a triangle's corners cover it when each is at
    // least its circumradius; the radius falls off along the wall at a
    // bounded rate, so neighbouring balls stay comparable.
    let mut radius = vec![0.0_f64; count];
    let mut adjacent: Vec<Vec<usize>> = vec![Vec::new(); count];
    for t in triangles {
        let [a, b, c] = t.map(|v| points[v]);
        let area2 = (b - a).cross(c - a).length();
        let longest = (b - a).length().max((c - b).length()).max((a - c).length());
        let circum = (b - a).length() * (c - b).length() * (a - c).length() / (2.0 * area2);
        let r = if circum.is_finite() {
            circum.min(4.0 * longest)
        } else {
            4.0 * longest
        };
        for k in 0..3 {
            radius[t[k]] = radius[t[k]].max(COVER * r);
            adjacent[t[k]].push(t[(k + 1) % 3]);
            adjacent[t[(k + 1) % 3]].push(t[k]);
        }
    }
    let mut order: Vec<usize> = (0..count).collect();
    order.sort_by(|&a, &b| radius[b].total_cmp(&radius[a]).then(a.cmp(&b)));
    let mut changed = true;
    while changed {
        changed = false;
        for &v in &order {
            for &u in &adjacent[v] {
                let least = radius[v] - RADIUS_RATE * (points[u] - points[v]).length();
                if least > radius[u] * (1.0 + 1e-12) {
                    radius[u] = least;
                    changed = true;
                }
            }
        }
    }
    let mut balls = Balls {
        points,
        weight: radius
            .iter()
            .enumerate()
            .map(|(v, r)| {
                let spread = (v as f64 * 0.618_033_988_749_894_9).fract() * 2.0 - 1.0;
                r * r * (1.0 + JITTER * spread)
            })
            .collect(),
    };

    let mut mesh = Flips {
        tris: triangles.to_vec(),
        sides: HashMap::new(),
    };
    for i in 0..mesh.tris.len() {
        mesh.link(i);
    }
    let all_edges = |mesh: &Flips| {
        let mut edges: Vec<(usize, usize)> = mesh.sides.keys().copied().collect();
        edges.sort_unstable();
        edges
    };
    stats.flips = mesh.regular(&balls, all_edges(&mesh));
    // Two seeds in nearly one point: the four corners are nearly on one
    // sphere. A smaller ball at a corner off the common edge moves its
    // seed away from the edge.
    for _ in 0..ROUNDS {
        let seeds: Vec<Option<Vec3>> = mesh.tris.iter().map(|&t| balls.seed(t)).collect();
        let mut shrink: Vec<usize> = Vec::new();
        for e in all_edges(&mesh) {
            let [a, b] = mesh.sides[&e];
            if a == NONE || b == NONE {
                continue;
            }
            if let (Some(sa), Some(sb)) = (seeds[a], seeds[b])
                && (sa - sb).length() < APART * (points[e.0] - points[e.1]).length()
            {
                shrink.push(opposite(mesh.tris[a], e));
            }
        }
        shrink.sort_unstable();
        shrink.dedup();
        if shrink.is_empty() {
            break;
        }
        for &v in &shrink {
            balls.weight[v] *= 1.0 - SHRINK;
        }
        stats.shrunk += shrink.len();
        stats.flips += mesh.regular(&balls, all_edges(&mesh));
    }

    // The seeds that are outside every other ball and apart from their
    // neighbours.
    let tree = KdTree::new(points);
    let mut seeds: Vec<Option<Vec3>> = mesh
        .tris
        .iter()
        .map(|&t| {
            let s = balls.seed(t)?;
            let inside = tree
                .nearest_to(s, 24)
                .into_iter()
                .any(|(d2, x)| !t.contains(&x) && d2 < balls.weight[x] * (1.0 - INSIDE));
            (!inside).then_some(s)
        })
        .collect();
    for e in all_edges(&mesh) {
        let [a, b] = mesh.sides[&e];
        if a == NONE || b == NONE {
            continue;
        }
        if let (Some(sa), Some(sb)) = (seeds[a], seeds[b])
            && (sa - sb).length() < APART * (points[e.0] - points[e.1]).length()
        {
            seeds[b] = None;
        }
    }
    let mut out = Vec::new();
    let mut corners = Vec::new();
    let mut active = vec![false; count];
    for (i, s) in seeds.iter().enumerate() {
        if let Some(s) = *s {
            out.push((s, balls.size(mesh.tris[i])));
            corners.push(mesh.tris[i]);
            for v in mesh.tris[i] {
                active[v] = true;
            }
        }
    }
    stats.seeds = out.len();
    stats.dropped = mesh.tris.len() - out.len();
    WallSeeds {
        seeds: out,
        corners,
        radius: (0..count)
            .map(|v| {
                if active[v] {
                    balls.weight[v].sqrt()
                } else {
                    0.0
                }
            })
            .collect(),
        stats,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A strip of a cylinder wall of radius `r` about the z axis, seen from
    /// inside (`inside`: the fluid is inside, as in a pipe) or from outside
    /// (a post in the flow). Rows are offset so the triangles are irregular.
    fn cylinder(r: f64, inside: bool) -> (Vec<Vec3>, Vec<[usize; 3]>) {
        let (rows, columns) = (5, 9);
        let mut points = Vec::new();
        for i in 0..rows {
            for j in 0..columns {
                let wobble = ((i * 7 + j * 3) % 5) as f64 * 0.02;
                let angle = 0.11 * j as f64 + 0.04 * (i % 2) as f64 + wobble;
                points.push(Vec3::new(
                    r * angle.cos(),
                    r * angle.sin(),
                    0.1 * i as f64 + wobble,
                ));
            }
        }
        let mut triangles = Vec::new();
        for i in 0..rows - 1 {
            for j in 0..columns - 1 {
                let (a, b) = (i * columns + j, i * columns + j + 1);
                let (c, d) = (a + columns, b + columns);
                for t in [[a, b, d], [a, d, c]] {
                    // normals out of the fluid
                    let n = (points[t[1]] - points[t[0]]).cross(points[t[2]] - points[t[0]]);
                    let outward = n.dot(Vec3::new(points[t[0]].x, points[t[0]].y, 0.0)) > 0.0;
                    triangles.push(if outward == inside {
                        t
                    } else {
                        [t[0], t[2], t[1]]
                    });
                }
            }
        }
        (points, triangles)
    }

    /// Keep rule: the seed of a triangle is on the balls of its three
    /// corners, on the fluid side. Everything else rests on this: two seeds
    /// over one edge are then both at the ball radius from its ends, and
    /// their bisector holds the edge. A sign error in the foot (weights
    /// swapped) or the height taken along the normal instead of against it
    /// breaks it; the mesh tests would only show more sliver cells.
    #[test]
    fn seed_is_on_the_three_balls_on_the_fluid_side() {
        let points = [
            Vec3::new(10.0, 20.0, 5.0),
            Vec3::new(11.0, 20.1, 5.2),
            Vec3::new(10.3, 21.2, 4.9),
        ];
        let balls = Balls {
            points: &points,
            weight: vec![1.3, 1.1, 1.6],
        };
        let seed = balls
            .seed([0, 1, 2])
            .expect("the balls meet over the triangle");
        for (v, &corner) in points.iter().enumerate() {
            let d = seed - corner;
            assert!((d.dot(d) - balls.weight[v]).abs() < 1e-12, "corner {v}");
        }
        assert!((seed - points[0]).dot(balls.normal([0, 1, 2])) < 0.0);
        // balls too small to meet have no seed
        let small = Balls {
            points: &points,
            weight: vec![0.1, 0.1, 0.1],
        };
        assert!(small.seed([0, 1, 2]).is_none());
    }

    /// Keep rule: on a curved wall, concave or convex, every point of a
    /// wall triangle is nearest to the seed of its own triangle, so the
    /// Voronoi cell of a seed cuts the wall in exactly that triangle (no
    /// cut wall polygons, no sliver columns). A seed inside a foreign ball,
    /// a flip that is not carried to the seeds, or a radius below the
    /// circumradius breaks it. The mesh tests count seeds but cannot tell
    /// whether the cells conform.
    #[test]
    fn seed_cells_cut_the_wall_in_their_triangles() {
        for inside in [true, false] {
            let (points, triangles) = cylinder(1.0, inside);
            let found = wall_seeds(&points, &triangles);
            assert_eq!(found.stats.dropped, 0, "inside {inside}");
            assert_eq!(found.seeds.len(), triangles.len());
            // the triangles after the flips, in seed order
            let mut mesh = Flips {
                tris: triangles.clone(),
                sides: HashMap::new(),
            };
            for i in 0..mesh.tris.len() {
                mesh.link(i);
            }
            let weight: Vec<f64> = found.radius.iter().map(|r| r * r).collect();
            let balls = Balls {
                points: &points,
                weight,
            };
            let mut edges: Vec<(usize, usize)> = mesh.sides.keys().copied().collect();
            edges.sort_unstable();
            assert_eq!(mesh.regular(&balls, edges), found.stats.flips);
            for (i, t) in mesh.tris.iter().enumerate() {
                let own = found.seeds[i].0;
                for v in 0..3 {
                    let d = own - points[t[v]];
                    assert!((d.dot(d).sqrt() / found.radius[t[v]] - 1.0).abs() < 1e-9);
                }
                let [a, b, c] = t.map(|v| points[v]);
                for [la, lb] in [
                    [1.0 / 3.0, 1.0 / 3.0],
                    [0.9, 0.05],
                    [0.05, 0.9],
                    [0.05, 0.05],
                ] {
                    let x = a * la + b * lb + c * (1.0 - la - lb);
                    let reach = (own - x).length();
                    for (j, other) in found.seeds.iter().enumerate() {
                        assert!(
                            j == i || (other.0 - x).length() > reach,
                            "inside {inside}: seed {j} is nearer to triangle {i}"
                        );
                    }
                }
            }
        }
    }

    /// Keep rule: a pair of triangles that is not regular for the balls is
    /// flipped, and the seeds are those of the flipped pair. Four corners
    /// of a flat kite cut along its long diagonal: each seed of the given
    /// pair lies inside the ball of the corner across the diagonal. Without
    /// the flip both triangles lose their seed (two cut wall triangles);
    /// the cylinder of the test above needs no flip.
    #[test]
    fn pair_that_is_not_regular_is_flipped() {
        let points = [
            Vec3::new(-2.0, 0.0, 0.0),
            Vec3::new(0.0, -0.5, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(0.0, 0.5, 0.0),
        ];
        // fluid over the plane (z > 0): normals point down
        let triangles = [[0, 2, 1], [0, 3, 2]];
        let found = wall_seeds(&points, &triangles);
        assert_eq!(found.stats.flips, 1);
        assert_eq!((found.seeds.len(), found.stats.dropped), (2, 0));
        for (seed, _) in &found.seeds {
            assert!(seed.z > 0.0);
            // on the balls of the short diagonal's ends, which both
            // flipped triangles have
            for v in [1, 3] {
                let d = *seed - points[v];
                assert!((d.dot(d).sqrt() / found.radius[v] - 1.0).abs() < 1e-9);
            }
        }
    }

    /// Keep rule: a seed inside the ball of a wall vertex that is not a
    /// corner of its triangle is dropped, and a vertex no kept seed uses
    /// has no ball. Two wall sheets closer than the seed height (a gap
    /// thinner than the wall mesh): each seed is inside the balls of the
    /// sheet across. Keeping such a seed puts a site where another seed's
    /// cell must reach its corner, and the cells no longer cut the wall in
    /// their triangles. The fixtures have no gap this thin.
    #[test]
    fn seed_inside_a_ball_across_a_thin_gap_is_dropped() {
        let sheet = |z: f64| {
            [
                Vec3::new(0.0, 0.0, z),
                Vec3::new(1.0, 0.0, z),
                Vec3::new(0.5, 0.9, z),
            ]
        };
        let points: Vec<Vec3> = sheet(0.0).into_iter().chain(sheet(0.05)).collect();
        // fluid between the sheets: the lower normal points down, the
        // upper one up
        let triangles = [[0, 2, 1], [3, 4, 5]];
        let found = wall_seeds(&points, &triangles);
        assert_eq!((found.seeds.len(), found.stats.dropped), (0, 2));
        assert!(found.radius.iter().all(|&r| r == 0.0));
        // far apart, both keep their seeds
        let points: Vec<Vec3> = sheet(0.0).into_iter().chain(sheet(5.0)).collect();
        let found = wall_seeds(&points, &triangles);
        assert_eq!((found.seeds.len(), found.stats.dropped), (2, 0));
        assert!(found.radius.iter().all(|&r| r > 0.0));
    }
}
