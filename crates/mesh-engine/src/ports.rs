//! Port detection on a bare closed STL (BC-03, ports): find the planar end caps of
//! pipes, the openings where a duct's inlet and outlet sit, and write a patch
//! map with `walls`, `inlet` and `outlet`, so a nonconvex duct meshes without
//! hand-written facet ranges.
//!
//! The surface is split into coplanar regions: a flood fill from each seed
//! facet across edges to facets whose normal is within `planar_angle_deg` of
//! the seed's and whose corners lie on the seed's plane. Comparing with the
//! seed, not the neighbour, keeps a finely tessellated cylinder from
//! creeping into one region. A region is a candidate when it is a disc (its
//! boundary edges form one connected loop). Its evidence is the
//! isoperimetric quotient 4πA/P² (1 for a circle, 0.785 for a square) and
//! the rim: the share of rim edges where the wall turns away from the cap
//! (convex) at 60° to 120°, as at the flat end of a pipe. The score is the
//! product. Ties keep the larger area, then the lower facet.

use std::collections::HashMap;

use mesh_core::geometry::Vec3;
use serde::Serialize;

use crate::artifact::point;
use crate::patches::{
    PATCH_MAP_FORMAT, PATCH_MAP_FORMAT_VERSION, PatchDef, PatchMap, PatchRole, PatchSurface,
};
use crate::{EngineError, ErrorCode};

/// Default largest angle between a facet and the seed of its planar region.
pub const DEFAULT_PLANAR_ANGLE_DEG: f64 = 1.0;
/// Default number of candidates reported.
pub const DEFAULT_CANDIDATES: usize = 8;
/// Corner distance from the seed plane, relative to the longest extent.
const RELATIVE_PLANE_TOLERANCE: f64 = 1e-5;

/// Options of a detection.
#[derive(Clone, Copy, Debug)]
pub struct PortOptions {
    pub planar_angle_deg: f64,
    pub candidates: usize,
}

impl Default for PortOptions {
    fn default() -> Self {
        Self {
            planar_angle_deg: DEFAULT_PLANAR_ANGLE_DEG,
            candidates: DEFAULT_CANDIDATES,
        }
    }
}

/// One planar disc that may be a port.
#[derive(Clone, Debug, Serialize)]
pub struct PortCandidate {
    /// Rank, 0 the most port-like.
    pub id: usize,
    pub facets: usize,
    /// Lowest and highest facet index (the facets need not be contiguous).
    pub facet_span: [usize; 2],
    pub area: f64,
    /// Area-weighted centroid.
    pub centre: [f64; 3],
    /// Outward normal.
    pub normal: [f64; 3],
    /// Radius of the circle of equal area.
    pub equivalent_radius: f64,
    /// 4πA/P²: 1 for a circle.
    pub circularity: f64,
    pub rim_edges: usize,
    /// Share of rim edges where the wall turns away at 60° to 120°.
    pub pipe_end_rim: f64,
    /// `[min, max]` angle between the cap and the wall across the rim.
    pub rim_angle_deg: [f64; 2],
    pub score: f64,
    #[serde(skip)]
    pub members: Vec<usize>,
}

/// The candidates, best first, and the number of planar regions examined.
#[derive(Clone, Debug, Serialize)]
pub struct Detection {
    pub planar_regions: usize,
    pub disc_regions: usize,
    pub candidates: Vec<PortCandidate>,
}

fn invalid(m: impl std::fmt::Display) -> EngineError {
    EngineError::new(ErrorCode::InvalidArguments, m)
}

impl PatchSurface {
    /// The planar disc regions of the surface ranked as pipe end caps.
    pub fn detect_ports(&self, options: PortOptions) -> Result<Detection, EngineError> {
        let angle = options.planar_angle_deg;
        if !angle.is_finite() || !(0.0..=30.0).contains(&angle) {
            return Err(invalid("planar_angle_deg must be from 0 to 30"));
        }
        if options.candidates == 0 {
            return Err(invalid("candidates must be at least 1"));
        }
        let cos = angle.to_radians().cos() - 64.0 * f64::EPSILON;
        let g = &self.geometry;
        let facets = self.surface.facets();
        let tol = RELATIVE_PLANE_TOLERANCE * self.surface.bounds().extent().max_abs_component();
        let n = g.len();
        let mut region = vec![usize::MAX; n];
        let mut regions: Vec<Vec<usize>> = Vec::new();
        for seed in 0..n {
            if region[seed] != usize::MAX {
                continue;
            }
            let id = regions.len();
            let normal = g.normals[seed];
            let origin = g.corners(seed)[0];
            let on_plane = |f: usize| {
                g.normals[f].dot(normal) >= cos
                    && g.corners(f)
                        .iter()
                        .all(|&v| (v - origin).dot(normal).abs() <= tol)
            };
            region[seed] = id;
            let mut members = vec![seed];
            let mut stack = vec![seed];
            while let Some(f) = stack.pop() {
                for &m in &g.neighbours[f] {
                    if region[m] == usize::MAX && on_plane(m) {
                        region[m] = id;
                        members.push(m);
                        stack.push(m);
                    }
                }
            }
            members.sort_unstable();
            regions.push(members);
        }
        let planar_regions = regions.len();
        let mut candidates = Vec::new();
        for (id, members) in regions.into_iter().enumerate() {
            // boundary edges and their connected components (loops)
            let mut parent: HashMap<usize, usize> = HashMap::new();
            fn find(p: &mut HashMap<usize, usize>, x: usize) -> usize {
                let mut r = x;
                while let Some(&q) = p.get(&r) {
                    if q == r {
                        break;
                    }
                    r = q;
                }
                p.insert(x, r);
                r
            }
            let mut area = 0.0;
            let mut weighted = Vec3::ZERO;
            let mut normal = Vec3::ZERO;
            let mut perimeter = 0.0;
            let mut rim = 0usize;
            let mut pipe_end = 0usize;
            let mut rim_angle = [f64::INFINITY, 0.0f64];
            for &f in &members {
                area += g.areas[f];
                weighted = weighted + g.centroids[f] * g.areas[f];
                normal = normal + g.normals[f] * g.areas[f];
            }
            let normal = normal / normal.length().max(f64::MIN_POSITIVE);
            for &f in &members {
                let v = facets[f].triangle.vertices;
                let c = g.corners(f);
                for k in 0..3 {
                    let m = g.neighbours[f][k];
                    if region[m] == id {
                        continue;
                    }
                    rim += 1;
                    let (a, b) = (v[k], v[(k + 1) % 3]);
                    for x in [a, b] {
                        parent.entry(x).or_insert(x);
                    }
                    let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                    if ra != rb {
                        parent.insert(ra, rb);
                    }
                    let (pa, pb) = (c[k], c[(k + 1) % 3]);
                    perimeter += (pb - pa).length();
                    let turn = g.normals[m]
                        .dot(normal)
                        .clamp(-1.0, 1.0)
                        .acos()
                        .to_degrees();
                    rim_angle = [rim_angle[0].min(turn), rim_angle[1].max(turn)];
                    let behind = (g.centroids[m] - (pa + pb) * 0.5).dot(normal) < 0.0;
                    if behind && (60.0..=120.0).contains(&turn) {
                        pipe_end += 1;
                    }
                }
            }
            if rim == 0 {
                continue;
            }
            let keys: Vec<usize> = parent.keys().copied().collect();
            let mut roots: Vec<usize> = keys.into_iter().map(|x| find(&mut parent, x)).collect();
            roots.sort_unstable();
            roots.dedup();
            if roots.len() != 1 {
                continue;
            }
            let circularity = 4.0 * std::f64::consts::PI * area / (perimeter * perimeter);
            let share = pipe_end as f64 / rim as f64;
            candidates.push(PortCandidate {
                id: 0,
                facets: members.len(),
                facet_span: [members[0], *members.last().unwrap()],
                area,
                centre: point(weighted / area),
                normal: point(normal),
                equivalent_radius: (area / std::f64::consts::PI).sqrt(),
                circularity,
                rim_edges: rim,
                pipe_end_rim: share,
                rim_angle_deg: rim_angle,
                score: circularity * share,
                members,
            });
        }
        let disc_regions = candidates.len();
        candidates.sort_by(|a, b| {
            b.score
                .total_cmp(&a.score)
                .then(b.area.total_cmp(&a.area))
                .then(a.facet_span[0].cmp(&b.facet_span[0]))
        });
        candidates.truncate(options.candidates.max(2));
        for (i, c) in candidates.iter_mut().enumerate() {
            c.id = i;
        }
        Ok(Detection {
            planar_regions,
            disc_regions,
            candidates,
        })
    }

    /// A map with `walls`, `inlet` (candidate `inlet`) and `outlet`
    /// (candidate `outlet`).
    pub fn port_map(
        &self,
        detection: &Detection,
        inlet: usize,
        outlet: usize,
    ) -> Result<PatchMap, EngineError> {
        let pick = |i: usize, what: &str| {
            detection.candidates.get(i).ok_or_else(|| {
                EngineError::new(
                    ErrorCode::EmptySelection,
                    format!(
                        "{what} candidate {i} does not exist; {} candidates were found",
                        detection.candidates.len()
                    ),
                )
            })
        };
        if inlet == outlet {
            return Err(invalid("inlet and outlet must be different candidates"));
        }
        let (a, b) = (pick(inlet, "inlet")?, pick(outlet, "outlet")?);
        let mut facet_patch = vec![0; self.facet_count()];
        for &f in &a.members {
            facet_patch[f] = 1;
        }
        for &f in &b.members {
            facet_patch[f] = 2;
        }
        let def = |name: &str, role| PatchDef {
            name: name.into(),
            role,
            source_name: None,
        };
        let map = PatchMap {
            format: PATCH_MAP_FORMAT.into(),
            format_version: PATCH_MAP_FORMAT_VERSION.into(),
            source_sha256: self.sha256.clone(),
            facet_count: self.facet_count(),
            patches: vec![
                def("walls", PatchRole::Wall),
                def("inlet", PatchRole::Inlet),
                def("outlet", PatchRole::Outlet),
            ],
            facet_patch,
        };
        self.check(&map)?;
        Ok(map)
    }
}
