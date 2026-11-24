//! Local thickness of a closed triangulated surface of any shape, on both
//! sides of the wall: the gap between facing walls (inside) and the solid
//! between back-to-back walls (outside).
//!
//! The thickness at a wall point p with unit normal n is the diameter of
//! the medial ball there: the largest ball tangent to the wall at p that
//! the opposite walls leave free. It is found by the shrinking-ball method
//! of Ma, Bae & Choi, *3D medial axis point approximation using
//! nearest neighbors and the normal field* (The Visual Computer 2012):
//! start with a large ball centred at c = p + r·u (u = −n inside, +n
//! outside), find the wall point q closest to c, and while q lies inside
//! the ball replace r by the radius of the ball tangent at p through q,
//! r = |q − p|² / (2 (q − p)·u). Here the closest point is exact on the
//! triangles (a bounding-volume hierarchy), not on point samples.
//!
//! Only opposite walls count, as for the convex gap source of
//! `sizing::surface`: a facet stops the ball when its normal makes an angle
//! above 180° − `wedge_angle` with n. Without this rule the thickness would
//! collapse to zero at every sharp edge (the local feature size does), and
//! the facets of a curved wall next to p would stop the ball early.

#![allow(clippy::needless_range_loop)] // parallel per-facet arrays

use super::surface::{Primitive, SurfaceGeometry, Weighted};
use mesh_core::geometry::Vec3;
use mesh_geometry::patches::closest_point;

#[derive(Clone, Copy)]
struct Item {
    corners: [Vec3; 3],
    normal: Vec3,
}

#[derive(Clone, Copy)]
struct Node {
    min: Vec3,
    max: Vec3,
    start: usize,
    end: usize,
    children: Option<(usize, usize)>,
}

/// Oriented triangles for medial-ball queries.
pub struct OrientedSet {
    items: Vec<Item>,
    nodes: Vec<Node>,
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

fn build(items: &mut [Item], nodes: &mut Vec<Node>, start: usize, end: usize) -> usize {
    let mut min = items[start].corners[0];
    let mut max = min;
    for item in &items[start..end] {
        for v in item.corners {
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
        let key = |p: &Item| {
            let s = p.corners.iter().fold(Vec3::ZERO, |s, &v| s + v);
            [s.x, s.y, s.z][axis]
        };
        let mid = start + (end - start) / 2;
        items[start..end].select_nth_unstable_by(mid - start, |a, b| key(a).total_cmp(&key(b)));
        let l = build(items, nodes, start, mid);
        let r = build(items, nodes, mid, end);
        nodes[id].children = Some((l, r));
    }
    id
}

/// Shrinking stops after this many balls (each strictly smaller).
const MAX_BALLS: usize = 64;

impl OrientedSet {
    /// All facets of the surface with their outward normals; degenerate
    /// facets are left out.
    pub fn new(geometry: &SurfaceGeometry) -> Self {
        let mut items: Vec<Item> = (0..geometry.triangles.len())
            .filter(|&f| geometry.normals[f].is_finite() && geometry.areas[f] > 0.0)
            .map(|f| Item {
                corners: geometry.corners(f),
                normal: geometry.normals[f],
            })
            .collect();
        let mut nodes = Vec::new();
        if !items.is_empty() {
            let n = items.len();
            build(&mut items, &mut nodes, 0, n);
        }
        Self { items, nodes }
    }

    /// The point closest to `c`, nearer than `sqrt(limit_sq)`, on a facet
    /// whose normal m has n·m < `cos_limit`.
    fn closest_opposite(&self, c: Vec3, limit_sq: f64, n: Vec3, cos_limit: f64) -> Option<Vec3> {
        let mut best = None;
        let mut best_sq = limit_sq;
        if self.nodes.is_empty() {
            return None;
        }
        let mut stack = vec![0usize];
        while let Some(i) = stack.pop() {
            let node = self.nodes[i];
            if box_sq(c, node.min, node.max) >= best_sq {
                continue;
            }
            match node.children {
                None => {
                    for item in &self.items[node.start..node.end] {
                        if n.dot(item.normal) >= cos_limit {
                            continue;
                        }
                        let q = closest_point(c, item.corners);
                        let d = (c - q).dot(c - q);
                        if d < best_sq {
                            best_sq = d;
                            best = Some(q);
                        }
                    }
                }
                Some((l, r)) => {
                    let (dl, dr) = (
                        box_sq(c, self.nodes[l].min, self.nodes[l].max),
                        box_sq(c, self.nodes[r].min, self.nodes[r].max),
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

    /// The diameter of the medial ball at wall point `p` on the side `u`
    /// (a unit vector, −n for the inside and +n for the outside of a wall
    /// with outward normal `n`), or `None` when the opposite walls leave a
    /// ball of diameter `cap` free.
    pub fn thickness(&self, p: Vec3, n: Vec3, u: Vec3, wedge_angle: f64, cap: f64) -> Option<f64> {
        let cos_limit = -wedge_angle.cos();
        let mut r = cap / 2.0;
        let mut hit = false;
        for _ in 0..MAX_BALLS {
            let c = p + u * r;
            // a point counts as inside when it is clearly nearer than r
            let inside = r * (1.0 - 1e-9);
            let Some(q) = self.closest_opposite(c, inside * inside, n, cos_limit) else {
                break;
            };
            let d = q - p;
            let along = d.dot(u);
            if along <= 0.0 {
                break;
            }
            let next = d.dot(d) / (2.0 * along);
            if next.is_nan() || next >= r {
                break;
            }
            hit = true;
            r = next;
        }
        hit.then_some(2.0 * r)
    }
}

/// The settings of a thickness size source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThicknessOptions {
    /// Edge lengths asked across the local thickness: size = t / across.
    pub across: f64,
    pub wedge_angle: f64,
    pub min_size: f64,
    pub max_size: f64,
}

/// What a thickness source found.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ThicknessStats {
    /// Smallest gap between facing walls (`INFINITY`: none below the cap).
    pub min_gap: f64,
    /// Smallest solid thickness between back-to-back walls.
    pub min_wall: f64,
    /// Pieces whose size was held at `min_size`.
    pub clamped: usize,
    pub pieces: usize,
}

/// Pieces are split in four while the thickness varies by more than this
/// ratio over a piece or the piece is longer than this many thicknesses.
const RATIO: f64 = 1.25;
const SPAN: f64 = 2.0;
const DEPTH: usize = 6;

/// Primitives of a two-sided thickness source on the facets not marked
/// `skip`: each wall piece asks for its smallest thickness / `across`, held
/// over the thickness so the opposite wall gets the same size.
pub fn thickness_primitives(
    geometry: &SurfaceGeometry,
    skip: &[bool],
    options: ThicknessOptions,
) -> (Vec<Weighted>, ThicknessStats) {
    let set = OrientedSet::new(geometry);
    let cap = options.max_size * options.across;
    let mut stats = ThicknessStats {
        min_gap: f64::INFINITY,
        min_wall: f64::INFINITY,
        ..Default::default()
    };
    let mut out = Vec::new();
    for f in 0..geometry.triangles.len() {
        let n = geometry.normals[f];
        if skip[f] || !n.is_finite() || geometry.areas[f] <= 0.0 {
            continue;
        }
        let mut stack = vec![(geometry.corners(f), 0usize)];
        while let Some((tri, depth)) = stack.pop() {
            let [a, b, c] = tri;
            let mut lo = f64::INFINITY;
            let mut hi: f64 = 0.0;
            let mut missing = false;
            for p in [a, b, c, (a + b + c) / 3.0] {
                let gap = set.thickness(p, n, n * -1.0, options.wedge_angle, cap);
                let wall = set.thickness(p, n, n, options.wedge_angle, cap);
                stats.min_gap = stats.min_gap.min(gap.unwrap_or(f64::INFINITY));
                stats.min_wall = stats.min_wall.min(wall.unwrap_or(f64::INFINITY));
                match (gap, wall) {
                    (None, None) => missing = true,
                    (g, w) => {
                        let t = g.unwrap_or(f64::INFINITY).min(w.unwrap_or(f64::INFINITY));
                        lo = lo.min(t);
                        hi = hi.max(t);
                    }
                }
            }
            if !lo.is_finite() {
                continue;
            }
            let longest = (b - a).length().max((c - b).length()).max((a - c).length());
            let floor = options.min_size * options.across;
            let uneven = missing || hi > RATIO * lo || longest > SPAN * lo.max(floor);
            if uneven && depth < DEPTH && longest > floor {
                let (ab, bc, ca) = ((a + b) / 2.0, (b + c) / 2.0, (c + a) / 2.0);
                stack.push(([a, ab, ca], depth + 1));
                stack.push(([ab, b, bc], depth + 1));
                stack.push(([ca, bc, c], depth + 1));
                stack.push(([ab, bc, ca], depth + 1));
                continue;
            }
            let ask = lo / options.across;
            let size = ask.clamp(options.min_size, options.max_size);
            if size < options.max_size {
                stats.pieces += 1;
                stats.clamped += usize::from(ask < options.min_size);
                out.push(Weighted {
                    primitive: Primitive::Triangle(tri),
                    size,
                    reach: lo,
                });
            }
        }
    }
    (out, stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A closed box [0,x]×[0,y]×[0,z] with outward normals, 12 facets.
    fn cuboid(min: Vec3, max: Vec3) -> (Vec<Vec3>, Vec<[usize; 3]>) {
        let v: Vec<Vec3> = (0..8)
            .map(|i| {
                Vec3::new(
                    if i & 1 == 0 { min.x } else { max.x },
                    if i & 2 == 0 { min.y } else { max.y },
                    if i & 4 == 0 { min.z } else { max.z },
                )
            })
            .collect();
        let quads = [
            [0, 2, 3, 1],
            [4, 5, 7, 6],
            [0, 1, 5, 4],
            [2, 6, 7, 3],
            [0, 4, 6, 2],
            [1, 3, 7, 5],
        ];
        let t = quads
            .iter()
            .flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]])
            .collect();
        (v, t)
    }

    #[test]
    fn slab_gap_is_its_height_and_edges_do_not_collapse_it() {
        let (v, t) = cuboid(Vec3::ZERO, Vec3::new(10.0, 8.0, 0.5));
        let g = SurfaceGeometry::from_triangles(v, t);
        let set = OrientedSet::new(&g);
        let wedge = 45f64.to_radians();
        let up = Vec3::new(0.0, 0.0, 1.0);
        for p in [
            Vec3::new(5.0, 4.0, 0.5),
            Vec3::new(0.0, 0.0, 0.5),
            Vec3::new(9.99, 0.01, 0.5),
        ] {
            let gap = set.thickness(p, up, up * -1.0, wedge, 100.0).unwrap();
            assert!((gap - 0.5).abs() < 1e-9, "{gap}");
            assert_eq!(set.thickness(p, up, up, wedge, 100.0), None);
        }
        // the side wall sees the wall 10 away, not the slab faces at its edge
        let side = Vec3::new(-1.0, 0.0, 0.0);
        let gap = set
            .thickness(Vec3::new(0.0, 4.0, 0.25), side, side * -1.0, wedge, 100.0)
            .unwrap();
        assert!((gap - 10.0).abs() < 1e-9, "{gap}");
        // beyond the cap nothing is reported
        assert_eq!(
            set.thickness(Vec3::new(0.0, 4.0, 0.25), side, side * -1.0, wedge, 5.0),
            None
        );
    }

    /// An L-shaped (nonconvex) fluid: two boxes sharing a face region are
    /// not a closed manifold, so use a fluid box with a thin solid plate
    /// standing in it: two separate closed shells, box outward and plate
    /// inward (its normals point into the plate).
    #[test]
    fn thin_plate_inside_a_box_is_seen_from_the_outside_of_the_wall() {
        let (mut v, mut t) = cuboid(Vec3::ZERO, Vec3::new(10.0, 10.0, 10.0));
        let (pv, pt) = cuboid(Vec3::new(4.0, 2.0, 2.0), Vec3::new(4.3, 8.0, 8.0));
        let base = v.len();
        v.extend(pv);
        // flipped winding: the fluid's outward normal points into the plate
        t.extend(pt.iter().map(|q| [q[0] + base, q[2] + base, q[1] + base]));
        let g = SurfaceGeometry::from_triangles(v, t);
        let set = OrientedSet::new(&g);
        let wedge = 45f64.to_radians();
        // on the plate face x = 4, the fluid's outward normal is +x
        let n = Vec3::new(1.0, 0.0, 0.0);
        let p = Vec3::new(4.0, 5.0, 5.0);
        let wall = set.thickness(p, n, n, wedge, 100.0).unwrap();
        assert!((wall - 0.3).abs() < 1e-9, "{wall}");
        let gap = set.thickness(p, n, n * -1.0, wedge, 100.0).unwrap();
        assert!((gap - 4.0).abs() < 1e-9, "{gap}");

        let skip = vec![false; g.triangles.len()];
        let (prims, stats) = thickness_primitives(
            &g,
            &skip,
            ThicknessOptions {
                across: 1.0,
                wedge_angle: wedge,
                min_size: 0.05,
                max_size: 3.0,
            },
        );
        assert!((stats.min_wall - 0.3).abs() < 1e-9);
        // the plate's 0.3-wide rim faces see the box wall 2 away, and the
        // plate 6 across behind them, not the 0.3 between their own edges
        assert!((stats.min_gap - 2.0).abs() < 1e-9, "{}", stats.min_gap);
        let smallest = prims.iter().map(|w| w.size).fold(f64::INFINITY, f64::min);
        assert!((smallest - 0.3).abs() < 1e-9);
    }
}
