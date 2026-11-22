//! Geometry-derived size sources (GRD-04): wall curvature, gaps between
//! opposite walls, and feature edges of a closed triangulated surface.
//!
//! Each source becomes a set of weighted primitives (triangles or segments,
//! each with its own size and reach). The source asks for
//!
//! h(x) = min_k [ h_k + r · max(0, dist(x, P_k) − reach_k − d) ]
//!
//! which is the closed-form gradient-limited field of the datum sources
//! (`field.rs`) applied per primitive: every term is r-Lipschitz, so the
//! minimum is too. A bounding-volume hierarchy with the smallest size and the
//! largest reach of each node prunes the minimum exactly (branch and bound).
//!
//! The estimators and their benchmark are in `research/grd-04-geometry-sizing.md`:
//!
//! - curvature: the edge-based curvature tensor of the normal cycle
//!   (Cohen-Steiner & Morvan 2003), averaged over the smooth vertex ring of
//!   each facet; size = angle / κ_max;
//! - gaps: the local thickness at a wall point is the diameter of the largest
//!   ball tangent to the wall there that the opposite walls leave free (the
//!   medial ball of the shrinking-ball method, Ma, Bae & Choi 2012), exact in
//!   closed form for a convex domain; size = thickness / cells across;
//! - features: edges whose facet normals turn more than the feature angle
//!   (the dihedral-angle rule, Botsch et al. 2010).

use mesh_core::geometry::Vec3;
use mesh_geometry::ClosedSurface;

/// A primitive of a geometry-derived source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Primitive {
    Triangle([Vec3; 3]),
    Segment([Vec3; 2]),
}

impl Primitive {
    /// The Euclidean distance from `p` to the primitive.
    pub fn distance(&self, p: Vec3) -> f64 {
        match *self {
            Self::Triangle([a, b, c]) => (p - closest_on_triangle(p, a, b, c)).length(),
            Self::Segment([a, b]) => {
                let ab = b - a;
                let t = ((p - a).dot(ab) / ab.dot(ab)).clamp(0.0, 1.0);
                (p - (a + ab * t)).length()
            }
        }
    }

    fn corners(&self) -> &[Vec3] {
        match self {
            Self::Triangle(v) => v,
            Self::Segment(v) => v,
        }
    }

    fn is_valid(&self) -> bool {
        let corners = self.corners();
        corners.iter().all(|v| v.is_finite())
            && match *self {
                Self::Triangle([a, b, c]) => (b - a).cross(c - a).length() > 0.0,
                Self::Segment([a, b]) => (b - a).length() > 0.0,
            }
    }
}

/// The closest point of triangle `abc` to `p`, by Voronoi regions of the
/// vertices, edges and face (Ericson 2005, section 5.1.5).
fn closest_on_triangle(p: Vec3, a: Vec3, b: Vec3, c: Vec3) -> Vec3 {
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let (d1, d2) = (ab.dot(ap), ac.dot(ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = p - b;
    let (d3, d4) = (ab.dot(bp), ac.dot(bp));
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return a + ab * (d1 / (d1 - d3));
    }
    let cp = p - c;
    let (d5, d6) = (ab.dot(cp), ac.dot(cp));
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return a + ac * (d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    let denominator = 1.0 / (va + vb + vc);
    a + ab * (vb * denominator) + ac * (vc * denominator)
}

/// A primitive with the size it asks for and the distance over which that
/// size holds before it grows.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Weighted {
    pub primitive: Primitive,
    pub size: f64,
    pub reach: f64,
}

/// Weighted primitives in a bounding-volume hierarchy.
#[derive(Clone, Debug, PartialEq)]
pub struct WeightedSet {
    items: Vec<Weighted>,
    nodes: Vec<Node>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Node {
    min: Vec3,
    max: Vec3,
    smallest: f64,
    reach: f64,
    start: usize,
    end: usize,
    children: Option<(usize, usize)>,
}

const LEAF: usize = 4;

impl WeightedSet {
    pub fn new(mut items: Vec<Weighted>) -> Self {
        let mut nodes = Vec::with_capacity(2 * items.len() / LEAF + 1);
        if !items.is_empty() {
            let len = items.len();
            build(&mut items, &mut nodes, 0, len);
        }
        Self { items, nodes }
    }

    pub fn items(&self) -> &[Weighted] {
        &self.items
    }

    /// Finite corners, nondegenerate primitives, positive sizes and
    /// nonnegative reaches.
    pub fn is_valid(&self) -> bool {
        self.items.iter().all(|w| {
            w.primitive.is_valid()
                && w.size.is_finite()
                && w.size > 0.0
                && w.reach.is_finite()
                && w.reach >= 0.0
        })
    }

    /// The smallest size of any primitive.
    pub fn smallest(&self) -> Option<f64> {
        self.nodes.first().map(|n| n.smallest)
    }

    /// `min(cap, min_k [h_k + rate · max(0, dist_k − reach_k − extra)])`.
    pub fn envelope(&self, p: Vec3, rate: f64, extra: f64, cap: f64) -> f64 {
        let mut best = cap;
        if self.nodes.is_empty() {
            return best;
        }
        let bound = |n: &Node| {
            n.smallest + rate * (box_distance(p, n.min, n.max).sqrt() - n.reach - extra).max(0.0)
        };
        let mut stack = vec![0usize];
        while let Some(index) = stack.pop() {
            let node = self.nodes[index];
            if bound(&node) >= best {
                continue;
            }
            match node.children {
                None => {
                    for w in &self.items[node.start..node.end] {
                        let beyond = (w.primitive.distance(p) - w.reach - extra).max(0.0);
                        best = best.min(w.size + rate * beyond);
                    }
                }
                Some((left, right)) => {
                    let (l, r) = (bound(&self.nodes[left]), bound(&self.nodes[right]));
                    // Visit the more promising child first (pushed last).
                    if l <= r {
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

    /// `min_k max(0, dist_k − reach_k)`: zero within the reach of a
    /// primitive, infinite for an empty set.
    pub fn reach_distance(&self, p: Vec3) -> f64 {
        let mut best = f64::INFINITY;
        if self.nodes.is_empty() {
            return best;
        }
        let mut stack = vec![0usize];
        while let Some(index) = stack.pop() {
            let node = self.nodes[index];
            if box_distance(p, node.min, node.max).sqrt() - node.reach >= best {
                continue;
            }
            match node.children {
                None => {
                    for w in &self.items[node.start..node.end] {
                        best = best.min((w.primitive.distance(p) - w.reach).max(0.0));
                    }
                }
                Some((left, right)) => {
                    stack.push(right);
                    stack.push(left);
                }
            }
        }
        best
    }
}

fn build(items: &mut [Weighted], nodes: &mut Vec<Node>, start: usize, end: usize) -> usize {
    let slice = &items[start..end];
    let mut min = slice[0].primitive.corners()[0];
    let mut max = min;
    let (mut smallest, mut reach) = (f64::INFINITY, 0.0f64);
    for w in slice {
        for &v in w.primitive.corners() {
            min = Vec3::new(min.x.min(v.x), min.y.min(v.y), min.z.min(v.z));
            max = Vec3::new(max.x.max(v.x), max.y.max(v.y), max.z.max(v.z));
        }
        smallest = smallest.min(w.size);
        reach = reach.max(w.reach);
    }
    let id = nodes.len();
    nodes.push(Node {
        min,
        max,
        smallest,
        reach,
        start,
        end,
        children: None,
    });
    if end - start > LEAF {
        let extent = max - min;
        let axis = if extent.x >= extent.y && extent.x >= extent.z {
            0
        } else if extent.y >= extent.z {
            1
        } else {
            2
        };
        let key = |w: &Weighted| {
            let c = w.primitive.corners();
            let s = c.iter().fold(Vec3::ZERO, |sum, &v| sum + v);
            [s.x, s.y, s.z][axis]
        };
        let mid = start + (end - start) / 2;
        items[start..end].select_nth_unstable_by(mid - start, |a, b| key(a).total_cmp(&key(b)));
        let left = build(items, nodes, start, mid);
        let right = build(items, nodes, mid, end);
        nodes[id].children = Some((left, right));
    }
    id
}

fn box_distance(q: Vec3, min: Vec3, max: Vec3) -> f64 {
    let d = |v: f64, lo: f64, hi: f64| (lo - v).max(v - hi).max(0.0);
    let (x, y, z) = (
        d(q.x, min.x, max.x),
        d(q.y, min.y, max.y),
        d(q.z, min.z, max.z),
    );
    x * x + y * y + z * z
}

/// One edge of the surface with its two facets and the angle between their
/// normals (0 for coplanar facets, π/2 for a cube edge).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Edge {
    pub vertices: [usize; 2],
    pub facets: [usize; 2],
    pub angle: f64,
}

/// The facet geometry of a closed surface in world coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceGeometry {
    pub vertices: Vec<Vec3>,
    pub triangles: Vec<[usize; 3]>,
    /// Outward unit normals.
    pub normals: Vec<Vec3>,
    pub areas: Vec<f64>,
    pub edges: Vec<Edge>,
    /// The edges of each facet.
    pub facet_edges: Vec<[usize; 3]>,
    /// The facets around each vertex.
    pub vertex_facets: Vec<Vec<usize>>,
}

impl SurfaceGeometry {
    pub fn new(surface: &ClosedSurface) -> Self {
        let vertices = surface.vertices().to_vec();
        let triangles: Vec<[usize; 3]> = surface
            .facets()
            .iter()
            .map(|f| f.triangle.vertices)
            .collect();
        let normals: Vec<Vec3> = surface.facets().iter().map(|f| f.plane.normal()).collect();
        let areas = triangles
            .iter()
            .map(|&[a, b, c]| {
                (vertices[b] - vertices[a])
                    .cross(vertices[c] - vertices[a])
                    .length()
                    / 2.0
            })
            .collect();
        Self::assemble(vertices, triangles, normals, areas)
    }

    /// The geometry of oriented triangles that form a closed manifold
    /// (checked by `ClosedSurface`; here only for the benchmark and tests).
    pub fn from_triangles(vertices: Vec<Vec3>, triangles: Vec<[usize; 3]>) -> Self {
        let mut normals = Vec::with_capacity(triangles.len());
        let mut areas = Vec::with_capacity(triangles.len());
        for &[a, b, c] in &triangles {
            let cross = (vertices[b] - vertices[a]).cross(vertices[c] - vertices[a]);
            normals.push(cross / cross.length());
            areas.push(cross.length() / 2.0);
        }
        Self::assemble(vertices, triangles, normals, areas)
    }

    fn assemble(
        vertices: Vec<Vec3>,
        triangles: Vec<[usize; 3]>,
        normals: Vec<Vec3>,
        areas: Vec<f64>,
    ) -> Self {
        let mut slots = std::collections::BTreeMap::<(usize, usize), Vec<usize>>::new();
        for (f, t) in triangles.iter().enumerate() {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                slots.entry((a.min(b), a.max(b))).or_default().push(f);
            }
        }
        let mut edges = Vec::with_capacity(slots.len());
        let mut facet_edges = vec![[usize::MAX; 3]; triangles.len()];
        let mut filled = vec![0usize; triangles.len()];
        for ((a, b), facets) in slots {
            let [f, g] = [facets[0], *facets.get(1).unwrap_or(&facets[0])];
            let (n, m) = (normals[f], normals[g]);
            let angle = n.cross(m).length().atan2(n.dot(m));
            let id = edges.len();
            edges.push(Edge {
                vertices: [a, b],
                facets: [f, g],
                angle,
            });
            for facet in [f, g] {
                if filled[facet] < 3 && !facet_edges[facet][..filled[facet]].contains(&id) {
                    facet_edges[facet][filled[facet]] = id;
                    filled[facet] += 1;
                }
            }
        }
        let mut vertex_facets = vec![Vec::new(); vertices.len()];
        for (f, t) in triangles.iter().enumerate() {
            for &v in t {
                vertex_facets[v].push(f);
            }
        }
        Self {
            vertices,
            triangles,
            normals,
            areas,
            edges,
            facet_edges,
            vertex_facets,
        }
    }

    pub fn corners(&self, facet: usize) -> [Vec3; 3] {
        self.triangles[facet].map(|v| self.vertices[v])
    }

    pub fn centroid(&self, facet: usize) -> Vec3 {
        let [a, b, c] = self.corners(facet);
        (a + b + c) / 3.0
    }

    fn edge_vector(&self, edge: &Edge) -> Vec3 {
        self.vertices[edge.vertices[1]] - self.vertices[edge.vertices[0]]
    }

    /// The facets that share a vertex with `facet` and whose normals turn
    /// at most `feature_angle` from its normal (the facet included).
    pub fn smooth_ring(&self, facet: usize, feature_angle: f64) -> Vec<usize> {
        let n = self.normals[facet];
        let mut ring: Vec<usize> = self.triangles[facet]
            .iter()
            .flat_map(|&v| self.vertex_facets[v].iter().copied())
            .filter(|&g| {
                let m = self.normals[g];
                n.cross(m).length().atan2(n.dot(m)) <= feature_angle
            })
            .collect();
        ring.sort_unstable();
        ring.dedup();
        ring
    }
}

/// The largest principal curvature of every facet by the normal-cycle
/// tensor (Cohen-Steiner & Morvan 2003): over a region B,
/// T = (1 / area(B)) Σ_e β_e |e ∩ B| ê êᵀ, with the edges' normal turn β_e.
/// Each edge counts half for each of its facets, B is the facet's smooth
/// vertex ring, and edges sharper than `feature_angle` are features, not
/// curvature. For a convex surface β_e ≥ 0 and κ_max is the largest
/// eigenvalue of T.
pub fn normal_cycle_curvature(geometry: &SurfaceGeometry, feature_angle: f64) -> Vec<f64> {
    // Per-facet sums S_f = ½ Σ_{e ∈ f} β_e |e| ê êᵀ (upper triangle).
    let sums: Vec<[f64; 6]> = (0..geometry.triangles.len())
        .map(|f| {
            let mut s = [0.0; 6];
            for &e in &geometry.facet_edges[f] {
                let edge = &geometry.edges[e];
                if edge.angle > feature_angle {
                    continue;
                }
                let v = geometry.edge_vector(edge);
                let length = v.length();
                let u = v / length;
                let w = 0.5 * edge.angle * length;
                add_outer(&mut s, u, w);
            }
            s
        })
        .collect();
    (0..geometry.triangles.len())
        .map(|f| {
            let ring = geometry.smooth_ring(f, feature_angle);
            let mut t = [0.0; 6];
            let mut area = 0.0;
            for &g in &ring {
                for k in 0..6 {
                    t[k] += sums[g][k];
                }
                area += geometry.areas[g];
            }
            for value in &mut t {
                *value /= area;
            }
            largest_eigenvalue(t).max(0.0)
        })
        .collect()
}

fn add_outer(s: &mut [f64; 6], u: Vec3, w: f64) {
    s[0] += w * u.x * u.x;
    s[1] += w * u.x * u.y;
    s[2] += w * u.x * u.z;
    s[3] += w * u.y * u.y;
    s[4] += w * u.y * u.z;
    s[5] += w * u.z * u.z;
}

/// The largest eigenvalue of the symmetric matrix [xx, xy, xz, yy, yz, zz]
/// by the trigonometric solution of the characteristic cubic (Smith 1961).
pub fn largest_eigenvalue(m: [f64; 6]) -> f64 {
    let [a, b, c, d, e, f] = m;
    let p1 = b * b + c * c + e * e;
    let q = (a + d + f) / 3.0;
    if p1 <= 1e-30 * (a * a + d * d + f * f).max(f64::MIN_POSITIVE) {
        return a.max(d).max(f);
    }
    let p2 = (a - q).powi(2) + (d - q).powi(2) + (f - q).powi(2) + 2.0 * p1;
    let p = (p2 / 6.0).sqrt();
    // B = (A − qI) / p; r = det(B) / 2.
    let (ba, bd, bf) = ((a - q) / p, (d - q) / p, (f - q) / p);
    let (bb, bc, be) = (b / p, c / p, e / p);
    let det = ba * (bd * bf - be * be) - bb * (bb * bf - be * bc) + bc * (bb * be - bd * bc);
    let phi = (det / 2.0).clamp(-1.0, 1.0).acos() / 3.0;
    q + 2.0 * p * phi.cos()
}

/// How the local thickness of a gap is measured at a wall point.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Thickness {
    /// The diameter of the largest ball tangent to the wall at the point
    /// that the opposite walls leave free (medial ball).
    MedialBall,
    /// The distance along the inward normal to the opposite wall (ray cast).
    NormalRay,
}

/// Planes of a convex surface for exact thickness queries.
#[derive(Clone, Debug, PartialEq)]
pub struct ConvexPlanes {
    normals: Vec<Vec3>,
    offsets: Vec<f64>,
}

impl ConvexPlanes {
    pub fn new(geometry: &SurfaceGeometry) -> Self {
        let offsets = (0..geometry.triangles.len())
            .map(|f| geometry.normals[f].dot(geometry.vertices[geometry.triangles[f][0]]))
            .collect();
        Self {
            normals: geometry.normals.clone(),
            offsets,
        }
    }

    /// The local thickness at `p` on a wall with outward normal `n`.
    /// Only walls facing it count: their normals make an angle above
    /// 180° − `wedge_angle` with `n`, so walls that meet at an interior angle
    /// below `wedge_angle` form a gap and a cube edge does not. For a convex
    /// domain the answer is exact: the ball centred at p − r·n lies inside
    /// every opposite half-space when r ≤ D_g / (1 − n·n_g), with D_g the
    /// distance of p from plane g; the ray leaves through plane g at
    /// D_g / (−n·n_g). `None` when no wall faces it.
    pub fn thickness(&self, p: Vec3, n: Vec3, wedge_angle: f64, method: Thickness) -> Option<f64> {
        let limit = -wedge_angle.cos();
        let mut best = f64::INFINITY;
        for (g, &m) in self.normals.iter().enumerate() {
            let c = n.dot(m);
            if c >= limit {
                continue;
            }
            let d = (self.offsets[g] - m.dot(p)).max(0.0);
            let t = match method {
                Thickness::MedialBall => 2.0 * d / (1.0 - c),
                Thickness::NormalRay => d / -c,
            };
            best = best.min(t);
        }
        best.is_finite().then_some(best)
    }
}

/// Primitives of a curvature source: each facet with κ > 0 asks for
/// `angle / κ`, clamped to `[min_size, max_size]`; facets at `max_size`
/// are left out because they ask for nothing.
pub fn curvature_primitives(
    geometry: &SurfaceGeometry,
    angle: f64,
    feature_angle: f64,
    min_size: f64,
    max_size: f64,
) -> Vec<Weighted> {
    normal_cycle_curvature(geometry, feature_angle)
        .into_iter()
        .enumerate()
        .filter(|&(_, k)| k > 0.0)
        .map(|(f, k)| Weighted {
            primitive: Primitive::Triangle(geometry.corners(f)),
            size: (angle / k).clamp(min_size, max_size),
            reach: 0.0,
        })
        .filter(|w| w.size < max_size)
        .collect()
}

/// Primitives of a feature-edge source: each edge sharper than
/// `feature_angle` asks for `size`.
pub fn feature_primitives(
    geometry: &SurfaceGeometry,
    feature_angle: f64,
    size: f64,
) -> Vec<Weighted> {
    geometry
        .edges
        .iter()
        .filter(|e| e.angle > feature_angle)
        .map(|e| Weighted {
            primitive: Primitive::Segment(e.vertices.map(|v| geometry.vertices[v])),
            size,
            reach: 0.0,
        })
        .collect()
}

/// The gap source's settings.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GapOptions {
    pub cells_across: f64,
    pub wedge_angle: f64,
    pub method: Thickness,
    pub min_size: f64,
    pub max_size: f64,
}

/// Facets are split in four until the thickness varies by at most this
/// ratio over a piece, or the depth limit is reached.
const GAP_RATIO: f64 = 1.1;
const GAP_DEPTH: usize = 5;

/// Primitives of a gap source: every wall piece asks for its smallest
/// thickness / `cells_across`, held over half the thickness so the size
/// reaches across the gap. The thickness is concave along a facet of a
/// convex domain (a minimum of affine functions), so its minimum over a
/// piece is at a corner: the size is exact on the corners and conservative
/// inside.
pub fn gap_primitives(
    geometry: &SurfaceGeometry,
    planes: &ConvexPlanes,
    options: GapOptions,
) -> Vec<Weighted> {
    let mut out = Vec::new();
    for f in 0..geometry.triangles.len() {
        let n = geometry.normals[f];
        let thickness = |p: Vec3| planes.thickness(p, n, options.wedge_angle, options.method);
        let mut stack = vec![(geometry.corners(f), 0usize)];
        while let Some((tri, depth)) = stack.pop() {
            let values: Vec<Option<f64>> = tri.iter().map(|&p| thickness(p)).collect();
            let known: Vec<f64> = values.iter().flatten().copied().collect();
            if known.is_empty() {
                continue;
            }
            let lo = known.iter().copied().fold(f64::INFINITY, f64::min);
            let hi = known.iter().copied().fold(0.0, f64::max);
            let uneven = known.len() < 3 || hi > GAP_RATIO * lo;
            if uneven && depth < GAP_DEPTH {
                let [a, b, c] = tri;
                let (ab, bc, ca) = ((a + b) / 2.0, (b + c) / 2.0, (c + a) / 2.0);
                stack.push(([a, ab, ca], depth + 1));
                stack.push(([ab, b, bc], depth + 1));
                stack.push(([ca, bc, c], depth + 1));
                stack.push(([ab, bc, ca], depth + 1));
                continue;
            }
            let size = (lo / options.cells_across).clamp(options.min_size, options.max_size);
            if size < options.max_size {
                out.push(Weighted {
                    primitive: Primitive::Triangle(tri),
                    size,
                    reach: lo / 2.0,
                });
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::sizing::shapes;

    fn geometry((v, f): shapes::Shape) -> SurfaceGeometry {
        SurfaceGeometry::from_triangles(v, f)
    }

    pub fn icosphere(r: f64, levels: usize) -> SurfaceGeometry {
        geometry(shapes::icosphere(Vec3::ZERO, r, levels))
    }

    pub fn strip_cylinder(r: f64, l: f64, n: usize) -> SurfaceGeometry {
        geometry(shapes::strip_cylinder(0.0, 0.0, r, l, n))
    }

    pub fn cuboid(max: Vec3) -> SurfaceGeometry {
        geometry(shapes::cuboid(max))
    }

    // Keep rule (SW-18): the closest-point routine and the pruned envelope
    // equal brute force. A wrong Voronoi-region branch or a node bound that
    // is too high (it would skip the primitive that sets h) breaks this; the
    // mesh tests only see sizes after relaxation.
    #[test]
    fn distances_and_envelope_equal_brute_force() {
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut r = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let tri = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.2, 0.0),
            Vec3::new(0.3, 0.9, 0.4),
        ];
        for _ in 0..500 {
            let p = Vec3::new(r() * 3.0 - 1.0, r() * 3.0 - 1.0, r() * 2.0 - 1.0);
            let mut best = f64::INFINITY;
            for i in 0..=200 {
                for j in 0..=(200 - i) {
                    let (u, v) = (i as f64 / 200.0, j as f64 / 200.0);
                    let q = tri[0] + (tri[1] - tri[0]) * u + (tri[2] - tri[0]) * v;
                    best = best.min((p - q).length());
                }
            }
            let exact = Primitive::Triangle(tri).distance(p);
            assert!(exact <= best + 1e-12 && best - exact < 0.01, "{p:?}");
        }
        let items: Vec<Weighted> = (0..300)
            .map(|i| {
                let c = Vec3::new(r() * 4.0, r() * 4.0, r());
                let primitive = if i % 3 == 0 {
                    Primitive::Segment([c, c + Vec3::new(r() * 0.3, r() * 0.3, 0.1)])
                } else {
                    Primitive::Triangle([
                        c,
                        c + Vec3::new(0.2 + r() * 0.2, 0.0, r() * 0.1),
                        c + Vec3::new(0.0, 0.2 + r() * 0.2, r() * 0.1),
                    ])
                };
                Weighted {
                    primitive,
                    size: 0.02 + r() * 0.2,
                    reach: r() * 0.1,
                }
            })
            .collect();
        let set = WeightedSet::new(items.clone());
        for _ in 0..400 {
            let p = Vec3::new(r() * 5.0 - 0.5, r() * 5.0 - 0.5, r() * 2.0 - 0.5);
            let (rate, extra, cap) = (0.3, 0.05, 0.8);
            let brute = items
                .iter()
                .map(|w| w.size + rate * (w.primitive.distance(p) - w.reach - extra).max(0.0))
                .fold(cap, f64::min);
            assert_eq!(set.envelope(p, rate, extra, cap), brute);
            let reach = items
                .iter()
                .map(|w| (w.primitive.distance(p) - w.reach).max(0.0))
                .fold(f64::INFINITY, f64::min);
            assert_eq!(set.reach_distance(p), reach);
        }
        assert_eq!(
            WeightedSet::new(vec![]).envelope(Vec3::ZERO, 0.2, 0.0, 0.7),
            0.7
        );
    }

    #[test]
    fn eigenvalues_of_known_matrices() {
        assert!((largest_eigenvalue([2.0, 0.0, 0.0, 5.0, 0.0, 1.0]) - 5.0).abs() < 1e-12);
        // [[2,1,0],[1,2,0],[0,0,0]] has eigenvalues 3, 1, 0.
        assert!((largest_eigenvalue([2.0, 1.0, 0.0, 2.0, 0.0, 0.0]) - 3.0).abs() < 1e-12);
        assert!((largest_eigenvalue([1.0, 1.0, 1.0, 1.0, 1.0, 1.0]) - 3.0).abs() < 1e-12);
    }

    // Keep rule (SW-18): the normal-cycle curvature recovers 1/R on an
    // icosphere and on a CAD-style strip cylinder, and zero on flat caps and
    // a box, with cube edges treated as features. A wrong edge weight, a
    // missing ½ per facet or feature edges counted as curvature breaks this;
    // the mesh test only checks sizes near one wall.
    #[test]
    fn normal_cycle_curvature_recovers_known_radii() {
        let feature = 30f64.to_radians();
        let sphere = icosphere(0.5, 3);
        let k = normal_cycle_curvature(&sphere, feature);
        for value in &k {
            assert!((value * 0.5 - 1.0).abs() < 0.08, "{value}");
        }
        let cylinder = strip_cylinder(0.25, 1.0, 48);
        let k = normal_cycle_curvature(&cylinder, feature);
        for (f, value) in k.iter().enumerate() {
            if cylinder.normals[f].z.abs() > 0.5 {
                assert_eq!(*value, 0.0, "cap facet {f}");
            } else {
                // Chord facets: θ/w = (π/n) / (R sin(π/n)), 1.0007/R here.
                assert!((value * 0.25 - 1.0).abs() < 0.01, "side facet {f}: {value}");
            }
        }
        let cube = cuboid(Vec3::new(1.0, 1.0, 1.0));
        assert!(
            normal_cycle_curvature(&cube, feature)
                .iter()
                .all(|&k| k == 0.0)
        );
        assert_eq!(feature_primitives(&cube, feature, 0.1).len(), 12);
        assert_eq!(feature_primitives(&cylinder, feature, 0.1).len(), 96);
    }

    // Keep rule (SW-18): the medial-ball and ray thickness are exact on a
    // slab and a sphere, cube edges are not gaps, and a gap source on a thin
    // slab asks for thickness / cells across over the whole gap. A sign error
    // in the plane distance, the wrong opposite-wall test or a missing reach
    // breaks this.
    #[test]
    fn thickness_and_gap_sizes() {
        let wedge = 45f64.to_radians();
        let slab = cuboid(Vec3::new(1.0, 1.0, 0.1));
        let planes = ConvexPlanes::new(&slab);
        let top = Vec3::new(0.0, 0.0, 1.0);
        for method in [Thickness::MedialBall, Thickness::NormalRay] {
            for p in [Vec3::new(0.5, 0.5, 0.1), Vec3::new(0.999, 0.001, 0.1)] {
                let t = planes.thickness(p, top, wedge, method).unwrap();
                assert!((t - 0.1).abs() < 1e-12, "{method:?} {t}");
            }
        }
        // A side wall at x = 1: the opposite wall is x = 0, not the top.
        let t = planes
            .thickness(
                Vec3::new(1.0, 0.5, 0.099),
                Vec3::new(1.0, 0.0, 0.0),
                wedge,
                Thickness::MedialBall,
            )
            .unwrap();
        assert!((t - 1.0).abs() < 1e-12);
        let sphere = icosphere(0.5, 2);
        let planes_s = ConvexPlanes::new(&sphere);
        let c = sphere.centroid(0);
        let t = planes_s
            .thickness(c, sphere.normals[0], wedge, Thickness::MedialBall)
            .unwrap();
        assert!((t - 1.0).abs() < 0.03, "{t}");
        let options = GapOptions {
            cells_across: 4.0,
            wedge_angle: wedge,
            method: Thickness::MedialBall,
            min_size: 0.001,
            max_size: 0.2,
        };
        let pieces = gap_primitives(&slab, &planes, options);
        // Top and bottom ask for 0.025 over the gap; the sides are 1 apart.
        let thin: Vec<_> = pieces.iter().filter(|w| w.size < 0.1).collect();
        assert_eq!(thin.len(), 4);
        for w in thin {
            assert!((w.size - 0.025).abs() < 1e-12 && (w.reach - 0.05).abs() < 1e-12);
        }
        let set = WeightedSet::new(pieces);
        let middle = set.envelope(Vec3::new(0.5, 0.5, 0.05), 0.3, 0.0, 0.2);
        assert!((middle - 0.025).abs() < 1e-12);
    }
}
