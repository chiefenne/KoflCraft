//! Facet selection on a closed triangulated shell, the geometry behind the
//! boundary patch tools (BC-01): plane, box and normal-direction
//! predicates, feature-angle region growing and splitting into regions.
//!
//! Every operation is an exact predicate on the facets in input coordinates
//! or a flood fill over the shell's edge adjacency; nothing is sampled or
//! approximated. The facet normal is the unit normal of the wound triangle
//! (outward on a positively oriented shell) and the feature angle between
//! two facets is the angle between their normals.

use std::collections::HashMap;

use mesh_core::geometry::Vec3;

use crate::ClosedSurface;

/// Per-facet geometry and edge adjacency of a closed shell.
#[derive(Clone, Debug)]
pub struct FacetGeometry {
    pub normals: Vec<Vec3>,
    pub centroids: Vec<Vec3>,
    pub areas: Vec<f64>,
    /// The facet across each edge (a→b, b→c, c→a) of a facet.
    pub neighbours: Vec<[usize; 3]>,
    corners: Vec<[Vec3; 3]>,
}

impl FacetGeometry {
    /// The geometry of a validated shell, whose edges are each shared by
    /// exactly two oppositely wound facets.
    pub fn new(surface: &ClosedSurface) -> Self {
        let vertices = surface.vertices();
        let facets = surface.facets();
        let mut directed = HashMap::with_capacity(3 * facets.len());
        for (index, facet) in facets.iter().enumerate() {
            let [a, b, c] = facet.triangle.vertices;
            for edge in [(a, b), (b, c), (c, a)] {
                directed.insert(edge, index);
            }
        }
        let mut geometry = Self {
            normals: Vec::with_capacity(facets.len()),
            centroids: Vec::with_capacity(facets.len()),
            areas: Vec::with_capacity(facets.len()),
            neighbours: Vec::with_capacity(facets.len()),
            corners: Vec::with_capacity(facets.len()),
        };
        for (index, facet) in facets.iter().enumerate() {
            let [a, b, c] = facet.triangle.vertices;
            let [p, q, r] = [a, b, c].map(|id| vertices[id]);
            let cross = (q - p).cross(r - p);
            let length = cross.length();
            geometry.normals.push(cross / length);
            geometry.areas.push(length / 2.0);
            geometry.centroids.push((p + q + r) / 3.0);
            geometry.corners.push([p, q, r]);
            // the closed shell pairs every edge, so the reverse edge exists;
            // fall back to the facet itself rather than panic
            geometry.neighbours.push(
                [(b, a), (c, b), (a, c)].map(|edge| directed.get(&edge).copied().unwrap_or(index)),
            );
        }
        geometry
    }

    pub fn len(&self) -> usize {
        self.normals.len()
    }

    pub fn is_empty(&self) -> bool {
        self.normals.is_empty()
    }

    /// The three corners of a facet in input coordinates.
    pub fn corners(&self, facet: usize) -> [Vec3; 3] {
        self.corners[facet]
    }

    /// Facets whose three corners lie within `tolerance` of the plane
    /// through `point` with unit normal `normal`.
    pub fn on_plane(&self, point: Vec3, normal: Vec3, tolerance: f64) -> Vec<bool> {
        self.corners
            .iter()
            .map(|c| {
                c.iter()
                    .all(|&v| (v - point).dot(normal).abs() <= tolerance)
            })
            .collect()
    }

    /// Facets in the closed box `[min, max]`: by centroid, or with all
    /// three corners inside when `all_corners`.
    pub fn in_box(&self, min: Vec3, max: Vec3, all_corners: bool) -> Vec<bool> {
        let inside = |v: Vec3| {
            v.x >= min.x
                && v.x <= max.x
                && v.y >= min.y
                && v.y <= max.y
                && v.z >= min.z
                && v.z <= max.z
        };
        (0..self.len())
            .map(|f| {
                if all_corners {
                    self.corners[f].iter().all(|&v| inside(v))
                } else {
                    inside(self.centroids[f])
                }
            })
            .collect()
    }

    /// Facets whose normal is within `max_angle` (radians) of the unit
    /// vector `direction`.
    pub fn facing(&self, direction: Vec3, max_angle: f64) -> Vec<bool> {
        let limit = cos_limit(max_angle);
        self.normals
            .iter()
            .map(|n| n.dot(direction) >= limit)
            .collect()
    }

    /// Whether the edge between two neighbouring facets is smooth: the
    /// angle between their normals is at most `feature_angle` (radians).
    pub fn smooth(&self, a: usize, b: usize, feature_angle: f64) -> bool {
        self.normals[a].dot(self.normals[b]) >= cos_limit(feature_angle)
    }

    /// Region growing: the facets reachable from `seed` across smooth
    /// edges, through facets with `allowed` set. The seed must be allowed.
    pub fn grow(&self, seed: usize, feature_angle: f64, allowed: &[bool]) -> Vec<bool> {
        let mut region = vec![false; self.len()];
        if !allowed.get(seed).copied().unwrap_or(false) {
            return region;
        }
        region[seed] = true;
        let mut stack = vec![seed];
        while let Some(f) = stack.pop() {
            for &n in &self.neighbours[f] {
                if allowed[n] && !region[n] && self.smooth(f, n, feature_angle) {
                    region[n] = true;
                    stack.push(n);
                }
            }
        }
        region
    }

    /// The regions of the facets with `members` set: connected across
    /// smooth edges (every edge when `feature_angle` is `None`), each sorted,
    /// ordered by their lowest facet.
    pub fn regions(&self, members: &[bool], feature_angle: Option<f64>) -> Vec<Vec<usize>> {
        let angle = feature_angle.unwrap_or(std::f64::consts::PI);
        let mut seen = vec![false; self.len()];
        let mut regions = Vec::new();
        for start in 0..self.len() {
            if !members[start] || seen[start] {
                continue;
            }
            seen[start] = true;
            let mut region = vec![start];
            let mut stack = vec![start];
            while let Some(f) = stack.pop() {
                for &n in &self.neighbours[f] {
                    if members[n] && !seen[n] && self.smooth(f, n, angle) {
                        seen[n] = true;
                        region.push(n);
                        stack.push(n);
                    }
                }
            }
            region.sort_unstable();
            regions.push(region);
        }
        regions
    }

    /// The facet nearest to `point`, by exact point–triangle distance; ties
    /// go to the lower index.
    pub fn nearest(&self, point: Vec3) -> usize {
        let mut best = (f64::INFINITY, 0);
        for (f, c) in self.corners.iter().enumerate() {
            let d = (closest_point(point, *c) - point).length();
            if d < best.0 {
                best = (d, f);
            }
        }
        best.1
    }
}

/// The cosine below which an angle exceeds `angle`, with a few ulps of
/// slack so that coplanar facets compare as parallel.
fn cos_limit(angle: f64) -> f64 {
    angle.clamp(0.0, std::f64::consts::PI).cos() - 64.0 * f64::EPSILON
}

/// The closest point to `p` on the triangle `[a, b, c]`, by the Voronoi
/// regions of its vertices, edges and face (Ericson, *Real-Time Collision
/// Detection*, 2005, section 5.1.5).
pub fn closest_point(p: Vec3, [a, b, c]: [Vec3; 3]) -> Vec3 {
    let ab = b - a;
    let ac = c - a;
    let ap = p - a;
    let d1 = ab.dot(ap);
    let d2 = ac.dot(ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = p - b;
    let d3 = ab.dot(bp);
    let d4 = ac.dot(bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return a + ab * (d1 / (d1 - d3));
    }
    let cp = p - c;
    let d5 = ab.dot(cp);
    let d6 = ac.dot(cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return a + ac * (d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    let denom = 1.0 / (va + vb + vc);
    a + ab * (vb * denom) + ac * (vc * denom)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Invariant: the closest point is a vertex, an edge point or the
    /// projection as the point's region says. A swapped region test returns
    /// a point off the triangle and fails here.
    #[test]
    fn closest_point_by_region() {
        let t = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let near = |p: Vec3, q: Vec3| (p - q).length() < 1e-15;
        assert!(near(closest_point(Vec3::new(-1.0, -1.0, 2.0), t), t[0]));
        assert!(near(closest_point(Vec3::new(3.0, -0.5, 0.0), t), t[1]));
        assert!(near(closest_point(Vec3::new(-0.2, 4.0, 1.0), t), t[2]));
        assert!(near(
            closest_point(Vec3::new(0.5, -2.0, 1.0), t),
            Vec3::new(0.5, 0.0, 0.0)
        ));
        assert!(near(
            closest_point(Vec3::new(1.0, 1.0, -3.0), t),
            Vec3::new(0.5, 0.5, 0.0)
        ));
        assert!(near(
            closest_point(Vec3::new(0.2, 0.3, 5.0), t),
            Vec3::new(0.2, 0.3, 0.0)
        ));
    }
}
