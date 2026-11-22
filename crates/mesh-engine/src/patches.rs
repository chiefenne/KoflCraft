//! Boundary patches defined after import (BC-01, BC-02): a patch map gives
//! every facet of a closed STL shell exactly one named patch with a role
//! (BC-04), and each operation (`patch.*` tools) reads one map and returns
//! a new one, so a sequence of calls is a journal that replays (AG-13).
//!
//! The map is checked against its surface on every read (BC-06): the STL
//! digest, the facet count, one valid patch index per facet, unique names
//! that are FoehnCraft boundary names (FC-04), and no empty patch. Because
//! the names already are hand-off names, MSH 4.1, the viewing VTU, the
//! FoehnCraft VTU and OpenFOAM all carry the same names.

use std::collections::{BTreeMap, BTreeSet};

use mesh_core::geometry::Vec3;
use mesh_geometry::{ClosedSurface, FacetGeometry, StlOptions, import_closed_stl};
pub use mesh_io::handoff::PatchRole;
use mesh_io::handoff::{handoff_name, patches_from_names};
use serde::{Deserialize, Serialize};

use crate::artifact::{point, sha256, vector};
use crate::{EngineError, ErrorCode};

pub const PATCH_MAP_FORMAT: &str = "koflcraft.patch-map";
pub const PATCH_MAP_FORMAT_VERSION: &str = "1.0";
pub const PATCH_SUMMARY_FORMAT: &str = "koflcraft.patch-summary";
/// Longest patch name.
pub const MAX_NAME_LEN: usize = 64;
/// Default plane tolerance, relative to the longest extent of the surface.
pub const RELATIVE_PLANE_TOLERANCE: f64 = 1e-6;

/// The `koflcraft.patch-map` file.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PatchMap {
    pub format: String,
    pub format_version: String,
    /// SHA-256 of the STL file the map belongs to.
    pub source_sha256: String,
    pub facet_count: usize,
    pub patches: Vec<PatchDef>,
    /// The index into `patches` of every facet, in file order.
    pub facet_patch: Vec<usize>,
}

/// One patch of a map.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PatchDef {
    pub name: String,
    pub role: PatchRole,
    /// The STL solid name the patch came from, when it differs.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_name: Option<String>,
}

/// A surface the patch tools work on: the closed shell (convex or not),
/// its facet geometry and its digest.
pub struct PatchSurface {
    surface: ClosedSurface,
    geometry: FacetGeometry,
    sha256: String,
}

fn err(code: ErrorCode, message: impl std::fmt::Display) -> EngineError {
    EngineError::new(code, message)
}

impl PatchSurface {
    /// Import and check a closed STL shell (any shape) for patch work.
    pub fn import(data: &[u8], facet_limit: Option<usize>) -> Result<Self, EngineError> {
        let surface = import_closed_stl(
            data,
            StlOptions {
                facet_patches: None,
                facet_limit,
            },
        )
        .map_err(|e| err(ErrorCode::InvalidDomain, e))?;
        let geometry = FacetGeometry::new(&surface);
        Ok(Self {
            surface,
            geometry,
            sha256: sha256(data),
        })
    }

    pub fn facet_count(&self) -> usize {
        self.geometry.len()
    }

    /// `[min, max]` of the surface.
    pub fn surface_bounds(&self) -> [[f64; 3]; 2] {
        let b = self.surface.bounds();
        [point(b.min()), point(b.max())]
    }

    /// A map from the STL solid names: one patch per solid (binary STL and
    /// unnamed solids: `surface`), names made hand-off names, roles from
    /// the names.
    pub fn initial_map(&self) -> PatchMap {
        let patches = patches_from_names(self.surface.patch_names())
            .into_iter()
            .map(|p| PatchDef {
                name: p.name,
                role: p.role,
                source_name: p.source_name,
            })
            .collect();
        PatchMap {
            format: PATCH_MAP_FORMAT.into(),
            format_version: PATCH_MAP_FORMAT_VERSION.into(),
            source_sha256: self.sha256.clone(),
            facet_count: self.facet_count(),
            patches,
            facet_patch: self
                .surface
                .facets()
                .iter()
                .map(|f| f.triangle.patch.0)
                .collect(),
        }
    }

    /// Check that `map` belongs to this surface and assigns every facet
    /// to exactly one valid patch (BC-06).
    pub fn check(&self, map: &PatchMap) -> Result<(), EngineError> {
        let bad = |m: String| Err(err(ErrorCode::InvalidPatchMap, m));
        if map.format != PATCH_MAP_FORMAT || map.format_version != PATCH_MAP_FORMAT_VERSION {
            return bad(format!(
                "expected format {PATCH_MAP_FORMAT} {PATCH_MAP_FORMAT_VERSION}, got {} {}",
                map.format, map.format_version
            ));
        }
        if map.source_sha256 != self.sha256 {
            return bad("the patch map belongs to another STL file (source_sha256 differs)".into());
        }
        if map.facet_count != self.facet_count() || map.facet_patch.len() != self.facet_count() {
            return bad(format!(
                "the STL has {} facets; the map declares {} and assigns {}",
                self.facet_count(),
                map.facet_count,
                map.facet_patch.len()
            ));
        }
        let mut names = BTreeSet::new();
        for p in &map.patches {
            check_name(&p.name)?;
            if !names.insert(p.name.as_str()) {
                return bad(format!("patch {} is defined twice", p.name));
            }
        }
        let mut used = vec![0usize; map.patches.len()];
        for (facet, &p) in map.facet_patch.iter().enumerate() {
            match used.get_mut(p) {
                Some(n) => *n += 1,
                None => {
                    return bad(format!(
                        "facet {facet} has patch index {p}; the map has {} patches",
                        map.patches.len()
                    ));
                }
            }
        }
        if let Some(p) = used.iter().position(|&n| n == 0) {
            return bad(format!("patch {} has no facets", map.patches[p].name));
        }
        Ok(())
    }

    /// The patch name of every facet, for meshing with this map.
    pub fn facet_names(&self, map: &PatchMap) -> Result<Vec<String>, EngineError> {
        self.check(map)?;
        Ok(map
            .facet_patch
            .iter()
            .map(|&p| map.patches[p].name.clone())
            .collect())
    }

    /// Geometric evidence of every patch and the BC-06 check.
    pub fn summary(&self, map: &PatchMap) -> Result<PatchSummary, EngineError> {
        self.check(map)?;
        let g = &self.geometry;
        let n = map.patches.len();
        let mut stats: Vec<PatchStats> = map
            .patches
            .iter()
            .map(|p| PatchStats {
                name: p.name.clone(),
                role: p.role,
                foehncraft_type: p.role.foehncraft_type(),
                source_name: p.source_name.clone(),
                facets: 0,
                area: 0.0,
                normal: None,
                max_normal_deviation_deg: 0.0,
                bounds: None,
                regions: 0,
                neighbours: Vec::new(),
            })
            .collect();
        let mut normal_sum = vec![Vec3::ZERO; n];
        let mut neighbours = vec![BTreeSet::new(); n];
        let mut bounds: Vec<Option<(Vec3, Vec3)>> = vec![None; n];
        for (f, &p) in map.facet_patch.iter().enumerate() {
            stats[p].facets += 1;
            stats[p].area += g.areas[f];
            normal_sum[p] = normal_sum[p] + g.normals[f] * g.areas[f];
            for v in g.corners(f) {
                let b = bounds[p].get_or_insert((v, v));
                b.0 = Vec3::new(b.0.x.min(v.x), b.0.y.min(v.y), b.0.z.min(v.z));
                b.1 = Vec3::new(b.1.x.max(v.x), b.1.y.max(v.y), b.1.z.max(v.z));
            }
            for &m in &g.neighbours[f] {
                let q = map.facet_patch[m];
                if q != p {
                    neighbours[p].insert(q);
                }
            }
        }
        for (p, s) in stats.iter_mut().enumerate() {
            let length = normal_sum[p].length();
            // a closed or nearly closed patch has no meaningful mean normal
            if length > 1e-9 * s.area {
                let mean = normal_sum[p] / length;
                s.normal = Some(point(mean));
                s.max_normal_deviation_deg = map
                    .facet_patch
                    .iter()
                    .enumerate()
                    .filter(|&(_, &q)| q == p)
                    .map(|(f, _)| g.normals[f].dot(mean).clamp(-1.0, 1.0).acos().to_degrees())
                    .fold(0.0, f64::max);
            }
            s.bounds = bounds[p].map(|(min, max)| [point(min), point(max)]);
            let members: Vec<bool> = map.facet_patch.iter().map(|&q| q == p).collect();
            s.regions = g.regions(&members, None).len();
            s.neighbours = neighbours[p]
                .iter()
                .map(|&q| map.patches[q].name.clone())
                .collect();
        }
        Ok(PatchSummary {
            format: PATCH_SUMMARY_FORMAT,
            format_version: PATCH_MAP_FORMAT_VERSION,
            source_sha256: self.sha256.clone(),
            facet_count: self.facet_count(),
            area: g.areas.iter().sum(),
            every_facet_in_one_patch: true,
            patches: stats,
        })
    }

    fn allowed(&self, map: &PatchMap, within: Option<&[String]>) -> Result<Vec<bool>, EngineError> {
        match within {
            None => Ok(vec![true; self.facet_count()]),
            Some(names) => {
                let ids: BTreeSet<usize> = names
                    .iter()
                    .map(|n| patch_index(map, n))
                    .collect::<Result<_, _>>()?;
                Ok(map.facet_patch.iter().map(|p| ids.contains(p)).collect())
            }
        }
    }

    /// Move the facets a selection picks (within `target.within`, if given)
    /// into the patch `target.name`, creating it when it is new.
    pub fn select(
        &self,
        map: &PatchMap,
        selection: &Selection,
        target: &Target,
    ) -> Result<(PatchMap, Change), EngineError> {
        self.check(map)?;
        let allowed = self.allowed(map, target.within.as_deref())?;
        let g = &self.geometry;
        let picked = match selection {
            Selection::Plane {
                point,
                normal,
                tolerance,
            } => {
                let normal = unit(*normal, "normal")?;
                let tol = match tolerance {
                    Some(t) if t.is_finite() && *t >= 0.0 => *t,
                    Some(_) => {
                        return Err(err(
                            ErrorCode::InvalidArguments,
                            "tolerance must be finite and not negative",
                        ));
                    }
                    None => {
                        RELATIVE_PLANE_TOLERANCE
                            * self.surface.bounds().extent().max_abs_component()
                    }
                };
                g.on_plane(finite(*point, "point")?, normal, tol)
            }
            Selection::Box { min, max, contain } => {
                let (min, max) = (finite(*min, "min")?, finite(*max, "max")?);
                if min.x > max.x || min.y > max.y || min.z > max.z {
                    return Err(err(
                        ErrorCode::InvalidArguments,
                        "box min must not exceed max",
                    ));
                }
                g.in_box(min, max, matches!(contain, Contain::AllCorners))
            }
            Selection::Normal {
                direction,
                max_angle_deg,
            } => g.facing(unit(*direction, "direction")?, angle(*max_angle_deg)?),
            Selection::Grow {
                seed_facet,
                seed_point,
                feature_angle_deg,
            } => {
                let seed = match (seed_facet, seed_point) {
                    (Some(f), None) if *f < self.facet_count() => *f,
                    (Some(f), None) => {
                        return Err(err(
                            ErrorCode::InvalidArguments,
                            format!(
                                "seed_facet {f} is not below the facet count {}",
                                self.facet_count()
                            ),
                        ));
                    }
                    (None, Some(p)) => g.nearest(finite(*p, "seed_point")?),
                    _ => {
                        return Err(err(
                            ErrorCode::InvalidArguments,
                            "give exactly one of seed_facet and seed_point",
                        ));
                    }
                };
                if !allowed[seed] {
                    return Err(err(
                        ErrorCode::EmptySelection,
                        format!("the seed facet {seed} is not in the patches of `within`"),
                    ));
                }
                g.grow(seed, angle(*feature_angle_deg)?, &allowed)
            }
        };
        let chosen: Vec<usize> = (0..self.facet_count())
            .filter(|&f| picked[f] && allowed[f])
            .collect();
        if chosen.is_empty() {
            return Err(err(
                ErrorCode::EmptySelection,
                "the selection contains no facets",
            ));
        }
        assign(map, &chosen, &target.name, target.role)
    }

    /// Split a patch into its regions: connected across edges whose
    /// feature angle is at most `feature_angle_deg` (all edges when
    /// `None`). Parts are named `names`, or `<patch>_1`, `<patch>_2`, …
    /// in the order of their lowest facet, and keep the patch's role.
    pub fn split(
        &self,
        map: &PatchMap,
        patch: &str,
        feature_angle_deg: Option<f64>,
        names: Option<&[String]>,
    ) -> Result<(PatchMap, Change), EngineError> {
        self.check(map)?;
        let id = patch_index(map, patch)?;
        let members: Vec<bool> = map.facet_patch.iter().map(|&q| q == id).collect();
        let regions = self
            .geometry
            .regions(&members, feature_angle_deg.map(angle).transpose()?);
        if regions.len() < 2 {
            return Err(err(
                ErrorCode::EmptySelection,
                format!("patch {patch} is one region; there is nothing to split"),
            ));
        }
        let names: Vec<String> = match names {
            Some(n) if n.len() == regions.len() => n.to_vec(),
            Some(n) => {
                return Err(err(
                    ErrorCode::InvalidArguments,
                    format!(
                        "patch {patch} splits into {} regions, {} names given",
                        regions.len(),
                        n.len()
                    ),
                ));
            }
            None => (1..=regions.len())
                .map(|k| format!("{patch}_{k}"))
                .collect(),
        };
        let mut unique = BTreeSet::new();
        for name in &names {
            check_name(name)?;
            if !unique.insert(name.as_str())
                || map
                    .patches
                    .iter()
                    .any(|p| &p.name == name && p.name != patch)
            {
                return Err(err(
                    ErrorCode::InvalidPatchName,
                    format!("patch name {name} is used already"),
                ));
            }
        }
        let old = &map.patches[id];
        let mut patches: Vec<PatchDef> = Vec::new();
        let mut index = Vec::with_capacity(map.patches.len());
        for (q, p) in map.patches.iter().enumerate() {
            index.push(patches.len());
            if q == id {
                for name in &names {
                    patches.push(PatchDef {
                        name: name.clone(),
                        role: old.role,
                        source_name: old.source_name.clone(),
                    });
                }
            } else {
                patches.push(p.clone());
            }
        }
        let mut facet_patch: Vec<usize> = map.facet_patch.iter().map(|&q| index[q]).collect();
        for (k, region) in regions.iter().enumerate() {
            for &f in region {
                facet_patch[f] = index[id] + k;
            }
        }
        let changed = regions[1..].iter().map(Vec::len).sum();
        let new = PatchMap {
            patches,
            facet_patch,
            ..map.clone()
        };
        Ok((
            new,
            Change {
                facets_changed: changed,
                created: names,
                removed: vec![patch.to_string()],
            },
        ))
    }
}

/// Merge patches into one named `name` (one of them or a new name), at
/// the place of the first; its role is `role`, else the first's.
pub fn merge(
    map: &PatchMap,
    patches: &[String],
    name: &str,
    role: Option<PatchRole>,
) -> Result<(PatchMap, Change), EngineError> {
    check_name(name)?;
    if patches.is_empty() {
        return Err(err(
            ErrorCode::InvalidArguments,
            "name at least one patch to merge",
        ));
    }
    let ids: BTreeSet<usize> = patches
        .iter()
        .map(|n| patch_index(map, n))
        .collect::<Result<_, _>>()?;
    if let Some(other) = map
        .patches
        .iter()
        .enumerate()
        .find(|(q, p)| p.name == name && !ids.contains(q))
    {
        return Err(err(
            ErrorCode::InvalidPatchName,
            format!(
                "patch {} exists and is not merged; add it to `patches`",
                other.1.name
            ),
        ));
    }
    let first = patch_index(map, &patches[0])?;
    let chosen: Vec<usize> = (0..map.facet_patch.len())
        .filter(|&f| ids.contains(&map.facet_patch[f]))
        .collect();
    let role = role.unwrap_or(map.patches[first].role);
    // rebuild: the merged patch sits where the first one was
    let mut defs = Vec::new();
    let mut index = vec![usize::MAX; map.patches.len()];
    let mut target = usize::MAX;
    for (q, p) in map.patches.iter().enumerate() {
        if ids.contains(&q) {
            if q == first {
                target = defs.len();
                defs.push(PatchDef {
                    name: name.to_string(),
                    role,
                    source_name: if ids.len() == 1 {
                        p.source_name.clone()
                    } else {
                        None
                    },
                });
            }
        } else {
            index[q] = defs.len();
            defs.push(p.clone());
        }
    }
    let facet_patch = map
        .facet_patch
        .iter()
        .map(|&q| if ids.contains(&q) { target } else { index[q] })
        .collect();
    let removed: Vec<String> = map
        .patches
        .iter()
        .enumerate()
        .filter(|(q, p)| ids.contains(q) && p.name != name)
        .map(|(_, p)| p.name.clone())
        .collect();
    let created = if map.patches.iter().any(|p| p.name == name) {
        Vec::new()
    } else {
        vec![name.to_string()]
    };
    let facets_changed = chosen
        .iter()
        .filter(|&&f| map.patches[map.facet_patch[f]].name != name)
        .count();
    if facets_changed == 0 && map.patches[first].role == role {
        return Err(err(ErrorCode::EmptySelection, "the merge changes nothing"));
    }
    Ok((
        PatchMap {
            patches: defs,
            facet_patch,
            ..map.clone()
        },
        Change {
            facets_changed,
            created,
            removed,
        },
    ))
}

/// Rename a patch; the new name must be free.
pub fn rename(map: &PatchMap, from: &str, to: &str) -> Result<(PatchMap, Change), EngineError> {
    check_name(to)?;
    let id = patch_index(map, from)?;
    if map.patches.iter().any(|p| p.name == to) {
        return Err(err(
            ErrorCode::InvalidPatchName,
            format!("patch {to} exists; use patch.merge to join patches"),
        ));
    }
    let mut new = map.clone();
    new.patches[id].name = to.to_string();
    Ok((
        new,
        Change {
            facets_changed: 0,
            created: vec![to.to_string()],
            removed: vec![from.to_string()],
        },
    ))
}

/// Set the role of a patch (BC-04).
pub fn set_role(
    map: &PatchMap,
    name: &str,
    role: PatchRole,
) -> Result<(PatchMap, Change), EngineError> {
    let id = patch_index(map, name)?;
    if map.patches[id].role == role {
        return Err(err(
            ErrorCode::EmptySelection,
            format!("patch {name} has this role already"),
        ));
    }
    let mut new = map.clone();
    new.patches[id].role = role;
    Ok((new, Change::default()))
}

/// What a selection picks.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Selection {
    /// Facets with all corners within `tolerance` of a plane.
    Plane {
        point: [f64; 3],
        normal: [f64; 3],
        tolerance: Option<f64>,
    },
    /// Facets in a closed axis-aligned box.
    Box {
        min: [f64; 3],
        max: [f64; 3],
        contain: Contain,
    },
    /// Facets whose outward normal is within an angle of a direction.
    Normal {
        direction: [f64; 3],
        max_angle_deg: f64,
    },
    /// Feature-angle region growing from a seed facet or the facet nearest
    /// to a point.
    Grow {
        seed_facet: Option<usize>,
        seed_point: Option<[f64; 3]>,
        feature_angle_deg: f64,
    },
}

/// Which facets a box contains.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Contain {
    /// The centroid is in the box.
    #[default]
    Centroid,
    /// All three corners are in the box.
    AllCorners,
}

/// Where selected facets go.
#[derive(Clone, Debug, PartialEq)]
pub struct Target {
    pub name: String,
    /// The role; a new patch without one gets the role its name suggests.
    pub role: Option<PatchRole>,
    /// Only facets of these patches are moved.
    pub within: Option<Vec<String>>,
}

/// What an operation changed.
#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct Change {
    pub facets_changed: usize,
    /// Patches the operation created.
    pub created: Vec<String>,
    /// Patches that no longer exist (emptied, merged, split or renamed).
    pub removed: Vec<String>,
}

/// The evidence of a patch map (`koflcraft.patch-summary`, BC-03 evidence
/// for agents, BC-06 check).
#[derive(Clone, Debug, Serialize)]
pub struct PatchSummary {
    pub format: &'static str,
    pub format_version: &'static str,
    pub source_sha256: String,
    pub facet_count: usize,
    pub area: f64,
    /// BC-06: every facet belongs to exactly one patch (the map was checked).
    pub every_facet_in_one_patch: bool,
    pub patches: Vec<PatchStats>,
}

/// Geometric evidence of one patch.
#[derive(Clone, Debug, Serialize)]
pub struct PatchStats {
    pub name: String,
    pub role: PatchRole,
    /// The FoehnCraft boundary-condition type the role maps to.
    pub foehncraft_type: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_name: Option<String>,
    pub facets: usize,
    pub area: f64,
    /// Area-weighted mean outward normal; `null` when the normals cancel.
    pub normal: Option<[f64; 3]>,
    /// Largest angle between a facet normal and the mean (0 for a plane).
    pub max_normal_deviation_deg: f64,
    /// `[min, max]` of the patch's corners.
    pub bounds: Option<[[f64; 3]; 2]>,
    /// Edge-connected regions of the patch.
    pub regions: usize,
    /// Patches that share an edge with this one.
    pub neighbours: Vec<String>,
}

/// Check that a patch name is a FoehnCraft boundary name: 1 to 64 ASCII
/// letters, digits, `_` or `-` (FC-04).
pub fn check_name(name: &str) -> Result<(), EngineError> {
    if name.is_empty() || name.len() > MAX_NAME_LEN || handoff_name(name) != name {
        return Err(err(
            ErrorCode::InvalidPatchName,
            format!(
                "patch name {name:?} must be 1 to {MAX_NAME_LEN} ASCII letters, digits, `_` or `-`"
            ),
        ));
    }
    Ok(())
}

fn patch_index(map: &PatchMap, name: &str) -> Result<usize, EngineError> {
    map.patches
        .iter()
        .position(|p| p.name == name)
        .ok_or_else(|| {
            let known: Vec<&str> = map.patches.iter().map(|p| p.name.as_str()).collect();
            err(
                ErrorCode::UnknownPatch,
                format!("no patch {name}; the map has {}", known.join(", ")),
            )
        })
}

/// Move `chosen` facets into the patch `name`; drop patches left empty.
fn assign(
    map: &PatchMap,
    chosen: &[usize],
    name: &str,
    role: Option<PatchRole>,
) -> Result<(PatchMap, Change), EngineError> {
    check_name(name)?;
    let mut patches = map.patches.clone();
    let mut created = Vec::new();
    let target = match patches.iter().position(|p| p.name == name) {
        Some(t) => t,
        None => {
            patches.push(PatchDef {
                name: name.to_string(),
                role: role.unwrap_or_else(|| PatchRole::from_name(name)),
                source_name: None,
            });
            created.push(name.to_string());
            patches.len() - 1
        }
    };
    let role_changed = role.is_some_and(|r| patches[target].role != r);
    if let Some(r) = role {
        patches[target].role = r;
    }
    let mut facet_patch = map.facet_patch.clone();
    let mut changed = 0;
    for &f in chosen {
        if facet_patch[f] != target {
            facet_patch[f] = target;
            changed += 1;
        }
    }
    if changed == 0 && !role_changed {
        return Err(err(
            ErrorCode::EmptySelection,
            format!("every selected facet is in patch {name} already"),
        ));
    }
    // drop the patches the move emptied, keeping the order of the rest
    let mut count = vec![0usize; patches.len()];
    for &p in &facet_patch {
        count[p] += 1;
    }
    let mut index = vec![usize::MAX; patches.len()];
    let mut kept = Vec::new();
    let mut removed = Vec::new();
    for (q, p) in patches.into_iter().enumerate() {
        if count[q] == 0 {
            removed.push(p.name);
        } else {
            index[q] = kept.len();
            kept.push(p);
        }
    }
    let facet_patch = facet_patch.into_iter().map(|q| index[q]).collect();
    Ok((
        PatchMap {
            patches: kept,
            facet_patch,
            ..map.clone()
        },
        Change {
            facets_changed: changed,
            created,
            removed,
        },
    ))
}

/// The roles of a surface's patches (in the surface's order) from a map.
pub(crate) fn roles_by_name(map: &PatchMap, names: &[String]) -> Vec<PatchRole> {
    let roles: BTreeMap<&str, PatchRole> = map
        .patches
        .iter()
        .map(|p| (p.name.as_str(), p.role))
        .collect();
    names
        .iter()
        .map(|n| {
            roles
                .get(n.as_str())
                .copied()
                .unwrap_or(PatchRole::Unspecified)
        })
        .collect()
}

fn finite(p: [f64; 3], what: &str) -> Result<Vec3, EngineError> {
    let v = vector(p);
    if v.is_finite() {
        Ok(v)
    } else {
        Err(err(
            ErrorCode::InvalidArguments,
            format!("{what} must be finite"),
        ))
    }
}

fn unit(p: [f64; 3], what: &str) -> Result<Vec3, EngineError> {
    let v = finite(p, what)?;
    let length = v.length();
    if length > 0.0 && length.is_normal() {
        Ok(v / length)
    } else {
        Err(err(
            ErrorCode::InvalidArguments,
            format!("{what} must not be zero"),
        ))
    }
}

fn angle(degrees: f64) -> Result<f64, EngineError> {
    if degrees.is_finite() && (0.0..=180.0).contains(&degrees) {
        Ok(degrees.to_radians())
    } else {
        Err(err(
            ErrorCode::InvalidArguments,
            "angles are degrees from 0 to 180",
        ))
    }
}
