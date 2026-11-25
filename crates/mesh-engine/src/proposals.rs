//! Patch proposals on a bare closed STL (BC-03): KoflCraft splits the
//! surface into the faces an engineer would name, says what each one is
//! (plane, cylinder, sphere, freeform), which ones look like openings and
//! which lie on a side of the bounding box, and writes a patch map an agent
//! or a person only has to rename and give roles, instead of picking
//! triangles.
//!
//! Method (clean room, from the papers only):
//! 1. Regions: facets connected across smooth edges, the edges whose
//!    dihedral turn is at most the feature angle (sharp-edge segmentation,
//!    as in the feature-edge step of hierarchical face clustering, Garland,
//!    Willmott and Heckbert 2001).
//! 2. Primitive per region from the area-weighted Gauss map, the tensor
//!    N = Σ a n nᵀ of the facet normals (Attene, Falcidieno and Spagnuolo
//!    2006, "Hierarchical mesh segmentation based on fitting primitives";
//!    the normal-cone tests of Lukács, Martin and Marshall 1998):
//!    - plane: every normal within `plane_angle_deg` of the mean;
//!    - cylinder: the normals are perpendicular to one axis (the eigenvector
//!      of the smallest eigenvalue of N), the normal lines meet the axis and
//!      the corners keep one distance from it;
//!    - sphere: the normal lines meet in one point and the corners keep one
//!      distance from it;
//!    - otherwise freeform.
//!
//!    The axis and centre are the least-squares points of the normal lines,
//!    the radius the mean corner distance; the relative spread of the corner
//!    distances is the evidence of the fit.
//! 3. Openings are the planar regions the port detection (ports.rs) ranks as
//!    round pipe end caps (score at least 0.9), or planar caps whose whole
//!    rim turns away at 60° to 120° into one surrounding wall region (the
//!    end of a pipe of any section). Box sides are the other planes on a
//!    side of the bounding box. Both get no role: which opening is the
//!    inlet, and whether a box side is a wall, a symmetry plane or a far
//!    field, is the engineer's call.
//!
//! Every result is deterministic: regions are ordered by their lowest facet,
//! names by kind and then area.

use std::collections::BTreeSet;

use mesh_core::geometry::Vec3;
use serde::Serialize;

use crate::artifact::point;
use crate::patches::{
    PATCH_MAP_FORMAT, PATCH_MAP_FORMAT_VERSION, PatchDef, PatchMap, PatchRole, PatchSurface,
};
use crate::ports::PortOptions;
use crate::{EngineError, ErrorCode};

pub const PROPOSAL_FORMAT: &str = "koflcraft.patch-proposal";
pub const DEFAULT_FEATURE_ANGLE_DEG: f64 = 30.0;
pub const DEFAULT_PLANE_ANGLE_DEG: f64 = 1.0;
/// Largest relative spread of the corner distances of a cylinder or sphere.
pub const DEFAULT_FIT_TOLERANCE: f64 = 0.02;
/// Smallest port-detection score (circularity × pipe-end rim share) of an
/// opening that is not capped by a single surrounding wall: 0.9 keeps round
/// pipe ends (≈ 1) and leaves out squares (π/4).
const OPENING_SCORE: f64 = 0.9;

/// Options of a proposal.
#[derive(Clone, Copy, Debug)]
pub struct ProposalOptions {
    pub feature_angle_deg: f64,
    pub plane_angle_deg: f64,
    pub fit_tolerance: f64,
    /// Join every region that is neither an opening nor a box side into
    /// one patch `walls`.
    pub merge_walls: bool,
}

impl Default for ProposalOptions {
    fn default() -> Self {
        Self {
            feature_angle_deg: DEFAULT_FEATURE_ANGLE_DEG,
            plane_angle_deg: DEFAULT_PLANE_ANGLE_DEG,
            fit_tolerance: DEFAULT_FIT_TOLERANCE,
            merge_walls: true,
        }
    }
}

/// What a region is.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Kind {
    Plane,
    Cylinder,
    Sphere,
    Freeform,
}

/// What the proposal thinks a region is for.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Suggestion {
    /// A planar disc at the end of a pipe: an inlet or an outlet.
    Opening,
    /// A plane on a side of the bounding box: a symmetry plane, a far
    /// field, an opening of a channel, or a wall.
    BoxSide,
    /// Everything else.
    Wall,
}

/// One region and its evidence.
#[derive(Clone, Debug, Serialize)]
pub struct Region {
    /// Index in file order of the lowest facet.
    pub id: usize,
    /// The patch the region was put in.
    pub patch: String,
    pub kind: Kind,
    pub suggestion: Suggestion,
    pub facets: usize,
    pub area: f64,
    /// Area-weighted centroid.
    pub centre: [f64; 3],
    /// Outward normal of a plane.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub normal: Option<[f64; 3]>,
    /// Unit axis of a cylinder.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub axis: Option<[f64; 3]>,
    /// A point on the axis of a cylinder, the centre of a sphere.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub axis_point: Option<[f64; 3]>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub radius: Option<f64>,
    /// Cylinder: length along the axis.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub length: Option<f64>,
    /// Largest normal deviation (plane) or relative radius spread
    /// (cylinder, sphere): how well the primitive fits.
    pub fit_error: f64,
    /// `xmin` … `zmax` for a box side.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub box_side: Option<&'static str>,
    /// Port-detection score of an opening (1 for a round pipe end).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub opening_score: Option<f64>,
    pub bounds: [[f64; 3]; 2],
    /// Ids of the regions across its rim.
    pub neighbours: Vec<usize>,
    #[serde(skip)]
    pub members: Vec<usize>,
}

/// The proposal: the regions and the map built from them.
#[derive(Clone, Debug, Serialize)]
pub struct Proposal {
    pub format: &'static str,
    pub feature_angle_deg: f64,
    pub regions: Vec<Region>,
    /// Patches whose role the proposal could not guess.
    pub needs_role: Vec<String>,
    #[serde(skip)]
    pub map: PatchMap,
}

fn invalid(m: impl std::fmt::Display) -> EngineError {
    EngineError::new(ErrorCode::InvalidArguments, m)
}

/// Eigenvalues (ascending) and unit eigenvectors of a symmetric 3×3
/// matrix by cyclic Jacobi rotations.
#[allow(clippy::needless_range_loop)] // rows and columns p, q of one matrix
pub(crate) fn symmetric_eigen(m: [[f64; 3]; 3]) -> ([f64; 3], [Vec3; 3]) {
    let mut a = m;
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..64 {
        let off = a[0][1].powi(2) + a[0][2].powi(2) + a[1][2].powi(2);
        let scale = a[0][0].powi(2) + a[1][1].powi(2) + a[2][2].powi(2);
        if off <= 1e-30 * scale.max(f64::MIN_POSITIVE) {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if a[p][q] == 0.0 {
                continue;
            }
            let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let t = if theta == 0.0 { 1.0 } else { t };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            for k in 0..3 {
                let (akp, akq) = (a[k][p], a[k][q]);
                a[k][p] = c * akp - s * akq;
                a[k][q] = s * akp + c * akq;
            }
            for k in 0..3 {
                let (apk, aqk) = (a[p][k], a[q][k]);
                a[p][k] = c * apk - s * aqk;
                a[q][k] = s * apk + c * aqk;
            }
            for row in &mut v {
                let (vp, vq) = (row[p], row[q]);
                row[p] = c * vp - s * vq;
                row[q] = s * vp + c * vq;
            }
        }
    }
    let mut order = [0, 1, 2];
    order.sort_by(|&i, &j| a[i][i].total_cmp(&a[j][j]));
    let values = order.map(|i| a[i][i]);
    let vectors = order.map(|i| Vec3::new(v[0][i], v[1][i], v[2][i]));
    (values, vectors)
}

fn outer(u: Vec3, w: f64, m: &mut [[f64; 3]; 3]) {
    let u = [u.x, u.y, u.z];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] += w * u[i] * u[j];
        }
    }
}

fn solve(m: [[f64; 3]; 3], b: Vec3) -> Option<Vec3> {
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(m);
    let scale = m.iter().flatten().map(|x| x.abs()).fold(0.0, f64::max);
    if d.is_nan() || d.abs() <= 1e-12 * scale.powi(3) {
        return None;
    }
    let b = [b.x, b.y, b.z];
    let mut x = [0.0; 3];
    for (k, xk) in x.iter_mut().enumerate() {
        let mut mk = m;
        for i in 0..3 {
            mk[i][k] = b[i];
        }
        *xk = det(mk) / d;
    }
    Some(Vec3::new(x[0], x[1], x[2]))
}

/// The point nearest (least squares, area-weighted) to the lines through
/// the facet centroids along the facet normals, restricted to the plane
/// perpendicular to `axis` when given (then the point is on the axis).
fn normal_line_point(
    surface: &PatchSurface,
    members: &[usize],
    axis: Option<Vec3>,
) -> Option<Vec3> {
    let g = &surface.geometry;
    let mut m = [[0.0; 3]; 3];
    let mut b = Vec3::ZERO;
    for &f in members {
        let n = g.normals[f];
        let a = g.areas[f];
        // projector I − n nᵀ (− d dᵀ)
        let mut p = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        outer(n, -1.0, &mut p);
        if let Some(d) = axis {
            outer(d, -1.0, &mut p);
        }
        let c = g.centroids[f];
        let pc = Vec3::new(
            p[0][0] * c.x + p[0][1] * c.y + p[0][2] * c.z,
            p[1][0] * c.x + p[1][1] * c.y + p[1][2] * c.z,
            p[2][0] * c.x + p[2][1] * c.y + p[2][2] * c.z,
        );
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] += a * p[i][j];
            }
        }
        b = b + pc * a;
    }
    if let Some(d) = axis {
        // pin the free direction along the axis to the area centroid
        let total: f64 = members.iter().map(|&f| g.areas[f]).sum();
        let centre = members
            .iter()
            .fold(Vec3::ZERO, |s, &f| s + g.centroids[f] * g.areas[f])
            / total;
        outer(d, total, &mut m);
        b = b + d * (total * d.dot(centre));
    }
    solve(m, b)
}

impl PatchSurface {
    /// Split the surface into named, classified regions and build a map.
    pub fn propose_patches(&self, options: ProposalOptions) -> Result<Proposal, EngineError> {
        let ProposalOptions {
            feature_angle_deg,
            plane_angle_deg,
            fit_tolerance,
            merge_walls,
        } = options;
        if !feature_angle_deg.is_finite() || !(1.0..=90.0).contains(&feature_angle_deg) {
            return Err(invalid("feature_angle_deg must be from 1 to 90"));
        }
        if !plane_angle_deg.is_finite() || !(0.0..=10.0).contains(&plane_angle_deg) {
            return Err(invalid("plane_angle_deg must be from 0 to 10"));
        }
        if !fit_tolerance.is_finite() || !(0.0..=0.5).contains(&fit_tolerance) {
            return Err(invalid("fit_tolerance must be from 0 to 0.5"));
        }
        let g = &self.geometry;
        let n = g.len();
        let all = vec![true; n];
        let groups = g.regions(&all, Some(feature_angle_deg.to_radians()));
        let mut region_of = vec![0usize; n];
        for (r, members) in groups.iter().enumerate() {
            for &f in members {
                region_of[f] = r;
            }
        }
        // openings: planar discs the port detection scores as pipe ends
        let detection = self.detect_ports(PortOptions {
            planar_angle_deg: plane_angle_deg.max(0.1),
            candidates: usize::MAX / 2,
        })?;
        let bounds = self.surface.bounds();
        let (lo, hi) = (bounds.min(), bounds.max());
        let extent = bounds.extent().max_abs_component();
        let tol = 1e-6 * extent;
        let cos_plane = plane_angle_deg.to_radians().cos() - 64.0 * f64::EPSILON;

        let mut regions = Vec::with_capacity(groups.len());
        for (r, members) in groups.into_iter().enumerate() {
            let mut area = 0.0;
            let mut weighted = Vec3::ZERO;
            let mut normal_sum = Vec3::ZERO;
            let mut tensor = [[0.0; 3]; 3];
            let (mut bmin, mut bmax) = (
                Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
                Vec3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
            );
            let mut neighbours = BTreeSet::new();
            for &f in &members {
                let a = g.areas[f];
                area += a;
                weighted = weighted + g.centroids[f] * a;
                normal_sum = normal_sum + g.normals[f] * a;
                outer(g.normals[f], a, &mut tensor);
                for v in g.corners(f) {
                    bmin = Vec3::new(bmin.x.min(v.x), bmin.y.min(v.y), bmin.z.min(v.z));
                    bmax = Vec3::new(bmax.x.max(v.x), bmax.y.max(v.y), bmax.z.max(v.z));
                }
                for &m in &g.neighbours[f] {
                    if region_of[m] != r {
                        neighbours.insert(region_of[m]);
                    }
                }
            }
            let centre = weighted / area;
            let mean = normal_sum / normal_sum.length().max(f64::MIN_POSITIVE);
            let min_cos = members
                .iter()
                .map(|&f| g.normals[f].dot(mean))
                .fold(1.0, f64::min);
            let corners = members.clone();
            let corner_spread = |distance: &dyn Fn(Vec3) -> f64| {
                let mut ds = Vec::with_capacity(3 * corners.len());
                for &f in &corners {
                    for v in g.corners(f) {
                        ds.push(distance(v));
                    }
                }
                let mean = ds.iter().sum::<f64>() / ds.len() as f64;
                let spread = ds.iter().map(|d| (d - mean).abs()).fold(0.0, f64::max);
                (mean, spread / mean.max(f64::MIN_POSITIVE))
            };
            let mut region = Region {
                id: members[0],
                patch: String::new(),
                kind: Kind::Freeform,
                suggestion: Suggestion::Wall,
                facets: members.len(),
                area,
                centre: point(centre),
                normal: None,
                axis: None,
                axis_point: None,
                radius: None,
                length: None,
                fit_error: 0.0,
                box_side: None,
                opening_score: None,
                bounds: [point(bmin), point(bmax)],
                neighbours: neighbours.into_iter().collect(),
                members,
            };
            if min_cos >= cos_plane {
                region.kind = Kind::Plane;
                region.normal = Some(point(mean));
                region.fit_error = min_cos.clamp(-1.0, 1.0).acos().to_degrees();
                for (axis, name) in [
                    (0usize, ["xmin", "xmax"]),
                    (1, ["ymin", "ymax"]),
                    (2, ["zmin", "zmax"]),
                ] {
                    let comp = |v: Vec3| [v.x, v.y, v.z][axis];
                    let flat = comp(bmax) - comp(bmin) <= tol;
                    if flat && comp(mean).abs() >= cos_plane {
                        if comp(mean) < 0.0 && (comp(bmin) - comp(lo)).abs() <= tol {
                            region.box_side = Some(name[0]);
                        } else if comp(mean) > 0.0 && (comp(bmax) - comp(hi)).abs() <= tol {
                            region.box_side = Some(name[1]);
                        }
                    }
                }
            } else {
                let (values, vectors) = symmetric_eigen(tensor);
                let axis = vectors[0];
                // normals perpendicular to one axis: mean sin² of the angle
                let perpendicular = values[0] / area;
                let (sin_tol, rel) = (plane_angle_deg.to_radians().sin().max(0.01), fit_tolerance);
                let mut fitted = false;
                if perpendicular <= sin_tol * sin_tol
                    && let Some(p) = normal_line_point(self, &region.members, Some(axis))
                {
                    let radial = |v: Vec3| {
                        let w = v - p;
                        (w - axis * w.dot(axis)).length()
                    };
                    let (radius, spread) = corner_spread(&radial);
                    if spread <= rel {
                        let along: Vec<f64> = region
                            .members
                            .iter()
                            .flat_map(|&f| g.corners(f))
                            .map(|v| (v - p).dot(axis))
                            .collect();
                        let lo = along.iter().copied().fold(f64::INFINITY, f64::min);
                        let hi = along.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                        region.kind = Kind::Cylinder;
                        region.axis = Some(point(axis));
                        region.axis_point = Some(point(p + axis * ((lo + hi) / 2.0)));
                        region.radius = Some(radius);
                        region.length = Some(hi - lo);
                        region.fit_error = spread;
                        fitted = true;
                    }
                }
                if !fitted && let Some(p) = normal_line_point(self, &region.members, None) {
                    let (radius, spread) = corner_spread(&|v: Vec3| (v - p).length());
                    if spread <= rel {
                        region.kind = Kind::Sphere;
                        region.axis_point = Some(point(p));
                        region.radius = Some(radius);
                        region.fit_error = spread;
                    }
                }
            }
            regions.push(region);
        }
        for c in &detection.candidates {
            let r = region_of[c.members[0]];
            // a round pipe end, or a cap whose whole rim turns into one
            // surrounding wall (a pipe of any section)
            let capped = c.pipe_end_rim >= 1.0 && regions[r].neighbours.len() == 1;
            if c.score < OPENING_SCORE && !capped {
                continue;
            }
            if regions[r].kind == Kind::Plane && regions[r].members == c.members {
                regions[r].suggestion = Suggestion::Opening;
                regions[r].opening_score = Some(c.score);
            }
        }
        for r in &mut regions {
            if r.suggestion == Suggestion::Wall && r.box_side.is_some() {
                r.suggestion = Suggestion::BoxSide;
            }
        }

        // names: openings by area, box sides by side, the rest by kind
        let mut order: Vec<usize> = (0..regions.len()).collect();
        let rank = |r: &Region| r.suggestion as usize;
        order.sort_by(|&a, &b| {
            rank(&regions[a])
                .cmp(&rank(&regions[b]))
                .then(regions[b].area.total_cmp(&regions[a].area))
                .then(regions[a].id.cmp(&regions[b].id))
        });
        let mut patches: Vec<PatchDef> = Vec::new();
        fn index_of(name: String, role: PatchRole, patches: &mut Vec<PatchDef>) -> usize {
            match patches.iter().position(|p| p.name == name) {
                Some(i) => i,
                None => {
                    patches.push(PatchDef {
                        name,
                        role,
                        source_name: None,
                    });
                    patches.len() - 1
                }
            }
        }
        let mut patch_of = vec![0usize; regions.len()];
        let (mut openings, mut counts) = (0usize, [0usize; 4]);
        let mut sides_used: Vec<&str> = Vec::new();
        for &r in &order {
            let region = &regions[r];
            let (name, role) = match region.suggestion {
                Suggestion::Opening => {
                    openings += 1;
                    (format!("opening_{openings}"), PatchRole::Unspecified)
                }
                Suggestion::BoxSide => {
                    let side = region.box_side.unwrap_or("side");
                    let k = sides_used.iter().filter(|&&s| s == side).count();
                    sides_used.push(side);
                    let name = if k == 0 {
                        side.to_string()
                    } else {
                        format!("{side}_{}", k + 1)
                    };
                    (name, PatchRole::Unspecified)
                }
                Suggestion::Wall if merge_walls => ("walls".to_string(), PatchRole::Wall),
                Suggestion::Wall => {
                    let k = region.kind as usize;
                    counts[k] += 1;
                    let kind = ["plane", "cylinder", "sphere", "wall"][k];
                    (format!("{kind}_{}", counts[k]), PatchRole::Wall)
                }
            };
            patch_of[r] = index_of(name, role, &mut patches);
        }
        // patch order: openings, box sides, walls; each by area
        let mut facet_patch = vec![0usize; n];
        for (r, region) in regions.iter_mut().enumerate() {
            region.patch = patches[patch_of[r]].name.clone();
            for &f in &region.members {
                facet_patch[f] = patch_of[r];
            }
        }
        let needs_role = patches
            .iter()
            .filter(|p| p.role == PatchRole::Unspecified)
            .map(|p| p.name.clone())
            .collect();
        let map = PatchMap {
            format: PATCH_MAP_FORMAT.into(),
            format_version: PATCH_MAP_FORMAT_VERSION.into(),
            source_sha256: self.sha256.clone(),
            facet_count: n,
            patches,
            facet_patch,
        };
        self.check(&map)?;
        Ok(Proposal {
            format: PROPOSAL_FORMAT,
            feature_angle_deg,
            regions,
            needs_role,
            map,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eigen_of_a_diagonal_and_a_rotated_matrix() {
        let (v, _) = symmetric_eigen([[3.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 2.0]]);
        assert_eq!(v, [1.0, 2.0, 3.0]);
        // [[2,1,0],[1,2,0],[0,0,5]]: 1, 3, 5
        let (v, e) = symmetric_eigen([[2.0, 1.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 5.0]]);
        for (a, b) in v.iter().zip([1.0, 3.0, 5.0]) {
            assert!((a - b).abs() < 1e-12);
        }
        let s = std::f64::consts::FRAC_1_SQRT_2;
        assert!((e[0].dot(Vec3::new(s, -s, 0.0)).abs() - 1.0).abs() < 1e-12);
    }
}
