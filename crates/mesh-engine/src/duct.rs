//! Mesh an annotated STL using a validated, explicitly supplied partition.
mod scaffold;
use crate::{
    EngineError, ErrorCode, MeshDirectory,
    artifact::{digest, require_schema, sha256, vector},
};
use mesh_algorithms::duct::{self, Region};
use mesh_geometry::{StlOptions, import_closed_stl};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::ffi::OsStr;
use std::time::Instant;

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Partition {
    pub schema_version: u32,
    pub source_sha256: String,
    pub regions: Vec<RegionConfig>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub volume_scaffold: Option<ScaffoldConfig>,
    pub first_thickness: f64,
    pub growth: f64,
    pub layers: usize,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ScaffoldConfig {
    pub file: String,
    sha256: String,
    max_relative_geometry_error: f64,
    /// Polyhedral core built on the scaffold.
    #[serde(default, skip_serializing_if = "CoreKind::is_default")]
    pub core: CoreKind,
    /// Lloyd iterations of the Voronoi core.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lloyd_iterations: Option<usize>,
}

/// The polyhedral core of a scaffold mesh.
#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum CoreKind {
    /// Clipped Voronoi diagram restricted to the scaffold, Lloyd-relaxed (M4).
    #[default]
    Voronoi,
    /// Vertex-centred median dual of the scaffold (M1b).
    MedianDual,
}

impl CoreKind {
    fn is_default(&self) -> bool {
        *self == Self::default()
    }
}

/// Most Lloyd iterations a partition may request.
pub const MAX_LLOYD_ITERATIONS: usize = 200;

impl Partition {
    pub fn scaffold_file(&self, mode: Mode) -> Option<&str> {
        if matches!(mode, Mode::Polyhedral) {
            self.volume_scaffold.as_ref().map(|s| s.file.as_str())
        } else {
            None
        }
    }
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct RegionConfig {
    kernel: [f64; 3],
    facets: Vec<[usize; 2]>,
}

#[derive(Deserialize, Serialize)]
pub struct PatchMap {
    schema_version: u32,
    source_sha256: String,
    pub facet_count: usize,
    default_patch: Patch,
    ports: [Port; 2],
}
#[derive(Deserialize, Serialize)]
struct Patch {
    name: String,
    role: String,
}
#[derive(Deserialize, Serialize)]
struct Port {
    name: String,
    role: String,
    facet_range: [usize; 2],
}
#[derive(Clone, Copy, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Mode {
    Bulk,
    Layered,
    Polyhedral,
}

pub fn generate(
    data: &[u8],
    patch_map: &PatchMap,
    plan: &Partition,
    mode: Mode,
    scaffold_data: Option<&[u8]>,
    dir: &MeshDirectory,
    name: &OsStr,
) -> Result<Value, EngineError> {
    let invalid = |s: &str| EngineError::new(ErrorCode::InvalidDomain, s);
    let start = Instant::now();
    let (surface, patches) = prepare(data, patch_map, plan)?;
    let source_digest = sha256(data);
    let mut steps = Vec::new();
    if matches!(mode, Mode::Layered | Mode::Polyhedral) {
        if plan.layers == 0
            || !plan.growth.is_finite()
            || plan.growth <= 0.
            || !plan.first_thickness.is_finite()
            || plan.first_thickness <= 0.
        {
            return Err(invalid(
                "layer count, growth and thickness must be positive",
            ));
        }
        plan.layers
            .checked_mul(surface.facets().len())
            .and_then(|count| count.checked_mul(std::mem::size_of::<duct::Cell>()))
            .ok_or_else(|| invalid("layer mesh size overflow"))?;
        steps
            .try_reserve_exact(plan.layers)
            .map_err(|_| invalid("layer count exceeds addressable memory"))?;
        let mut total = 0.0;
        let mut step = plan.first_thickness;
        for _ in 0..plan.layers {
            if !step.is_finite()
                || step <= 0.
                || total + step <= total
                || total + step >= surface.bounds().extent().max_abs_component()
            {
                return Err(invalid("unrepresentable or oversized layer thickness"));
            }
            total += step;
            steps.push(step);
            step *= plan.growth;
        }
    }
    let regions: Vec<_> = plan
        .regions
        .iter()
        .map(|r| Region {
            kernel: vector(r.kernel),
            facets: r.facets.clone(),
        })
        .collect();
    let validated = |e: String| EngineError::new(ErrorCode::ValidationFailed, e);
    let mut scaffold_volume = None;
    let mut voronoi: Option<duct::VoronoiStats> = None;
    let mut source_areas = [0.; 3];
    let scale = surface.bounds().extent().max_abs_component();
    for (f, &tag) in surface.facets().iter().zip(&patches) {
        source_areas[tag - 1] += f.normalized_area * scale * scale;
    }
    let mesh = if let (Mode::Polyhedral, Some(spec)) = (mode, &plan.volume_scaffold) {
        let bytes = scaffold_data.ok_or_else(|| invalid("missing configured volume scaffold"))?;
        if sha256(bytes) != spec.sha256 {
            return Err(invalid("volume scaffold digest mismatch"));
        }
        let limit = spec.max_relative_geometry_error;
        if !limit.is_finite() || !(0.0..1.0).contains(&limit) {
            return Err(invalid("invalid scaffold geometry error limit"));
        }
        let input = scaffold::read(bytes).map_err(validated)?;
        let (mesh, shell) = match spec.core {
            CoreKind::MedianDual => duct::build_volume_dual(input, &steps).map_err(validated)?,
            CoreKind::Voronoi => {
                let lloyd_iterations = spec
                    .lloyd_iterations
                    .unwrap_or(duct::DEFAULT_LLOYD_ITERATIONS);
                if lloyd_iterations > MAX_LLOYD_ITERATIONS {
                    return Err(invalid("too many Lloyd iterations"));
                }
                let (mesh, shell, stats) = duct::build_voronoi_core(
                    input,
                    &steps,
                    duct::VoronoiOptions { lloyd_iterations },
                )
                .map_err(validated)?;
                voronoi = Some(stats);
                (mesh, shell)
            }
        };
        scaffold_volume = Some(shell.volume());
        if (shell.volume() / surface.volume() - 1.).abs() > limit
            || source_areas
                .iter()
                .zip(mesh.measurements.areas)
                .any(|(&a, b)| (b / a - 1.).abs() > limit)
        {
            return Err(invalid(
                "remeshed volume or patch areas exceed configured STL error limit",
            ));
        }
        for tag in 1..=2 {
            let f = surface
                .facets()
                .iter()
                .zip(&patches)
                .find(|(_, p)| **p == tag)
                .unwrap()
                .0;
            let origin = surface.vertices()[f.triangle.vertices[0]];
            let normal = f.plane.normal();
            for face in shell
                .facets()
                .iter()
                .filter(|f| f.triangle.patch.0 == tag - 1)
            {
                if face.plane.normal().dot(normal) < 1. - 1e-9
                    || face
                        .triangle
                        .vertices
                        .iter()
                        .any(|&v| (shell.vertices()[v] - origin).dot(normal).abs() > scale * 1e-9)
                {
                    return Err(invalid(
                        "scaffold opening differs from STL port plane or winding",
                    ));
                }
            }
        }
        mesh
    } else if matches!(mode, Mode::Polyhedral) {
        duct::build_polyhedral(&surface, &patches, &regions, &steps).map_err(validated)?
    } else {
        duct::build(&surface, &patches, &regions, &steps).map_err(validated)?
    };
    let m = &mesh.measurements;
    let mut tetrahedra = 0;
    let mut prisms = 0;
    let mut polyhedra = 0;
    let mut face_range = [usize::MAX, 0];
    for cell in &mesh.cells {
        match &cell.element {
            duct::Element::Tet(_) => tetrahedra += 1,
            duct::Element::Wedge(_) => prisms += 1,
            duct::Element::Polyhedron(p) => {
                polyhedra += 1;
                face_range[0] = face_range[0].min(p.faces.len());
                face_range[1] = face_range[1].max(p.faces.len());
            }
        }
    }
    let report = json!({
        "schema_version": 1,
        "algorithm": if voronoi.is_some() {"volume_clipped_voronoi_cvt_v1"} else if scaffold_volume.is_some() {"volume_median_dual_v2"} else if polyhedra > 0 {"stl_partitioned_dual_v1"} else {"stl_partitioned_star_v1"},
        "mode": mode,
        "source_sha256": source_digest,
        "config_sha256": digest(&(patch_map, plan, mode))?,
        "regions": if scaffold_volume.is_some(){1}else{regions.len()},
        "scaffold_volume": scaffold_volume,
        "stl_patch_areas": {"inlet":source_areas[0],"outlet":source_areas[1],"walls":source_areas[2]},
        "source_facets": surface.facets().len(),
        "tetrahedra": tetrahedra,
        "polyhedra": polyhedra,
        "core_sphericity": m.core_sphericity.map(|[min,p05,median]| json!({"minimum":min,"p05":p05,"median":median})),
        "polyhedron_face_range": (polyhedra > 0).then_some(face_range),
        "prism_core_faces":m.prism_core_faces,
        "prisms": prisms,
        "source_volume": surface.volume(),
        "volume": m.volume,
        "minimum_cell_volume": m.minimum_volume,
        "minimum_tet_mean_ratio": (tetrahedra > 0).then_some(m.minimum_tet_mean_ratio),
        "tetrahedra_below_mean_ratio_0_1": m.tetrahedra_below_mean_ratio_0_1,
        "minimum_wedge_jacobian": m.minimum_wedge_jacobian,
        "minimum_kernel_clearance": m.minimum_kernel_clearance,
        "patch_areas": {"inlet": m.areas[0], "outlet": m.areas[1], "walls": m.areas[2]},
        "port_areas_by_layer": {"inlet": m.port_layer_areas[0], "outlet": m.port_layer_areas[1]},
        "thicknesses": steps,
        "normal_thickness_ranges": m.normal_thickness_ranges,
        "voronoi_core": voronoi.as_ref().map(voronoi_json),
        "generation_seconds": start.elapsed().as_secs_f64()
    });
    let mut volume = Vec::new();
    let mut gmsh = Vec::new();
    let mut faces: [Vec<u8>; 3] = Default::default();
    let io = |e| EngineError::new(ErrorCode::ExportFailed, e);
    mesh_io::duct::volume(&mesh, &mut volume).map_err(io)?;
    mesh_io::duct::gmsh(&mesh, mesh_io::GmshVersion::V4_1, &mut gmsh).map_err(io)?;
    for (i, out) in faces.iter_mut().enumerate() {
        mesh_io::duct::patch(&mesh, i + 1, out).map_err(io)?;
    }
    let report_json = serde_json::to_vec_pretty(&report)
        .map_err(|e| EngineError::new(ErrorCode::SerializationFailed, e))?;
    let handoff = mesh_io::duct::handoff(&mesh).map_err(crate::handoff::export_error)?;
    let foehn = crate::handoff::foehncraft_vtu(&handoff, &[])?;
    let mut files: Vec<(&str, &[u8])> = vec![
        ("volume.vtu", &volume),
        ("mesh.msh", &gmsh),
        ("mesh.vtm", mesh_io::channel::MULTIBLOCK),
        ("inlet.vtu", &faces[0]),
        ("outlet.vtu", &faces[1]),
        ("walls.vtu", &faces[2]),
        ("report.json", &report_json),
        (crate::handoff::FOEHNCRAFT_VTU, &foehn),
    ];
    let generator = report["algorithm"].as_str().unwrap_or("duct").to_string();
    let evidence = crate::handoff::evidence(&handoff, &generator, &files)?;
    files.push((crate::handoff::MESH_EVIDENCE, &evidence));
    dir.write_bundle(name, &files)
        .map_err(|e| EngineError::new(ErrorCode::IoError, e))?;
    Ok(report)
}

/// The surface and per-facet patch ids (1 inlet, 2 outlet, 3 wall) of a
/// checked STL, patch map and partition.
fn prepare(
    data: &[u8],
    patch_map: &PatchMap,
    plan: &Partition,
) -> Result<(mesh_geometry::ClosedSurface, Vec<usize>), EngineError> {
    let invalid = |s: &str| EngineError::new(ErrorCode::InvalidDomain, s);
    require_schema(patch_map.schema_version)?;
    require_schema(plan.schema_version)?;
    let source_digest = sha256(data);
    if source_digest != patch_map.source_sha256 || source_digest != plan.source_sha256 {
        return Err(invalid(
            "patch map or partition digest differs from input STL",
        ));
    }
    let surface = import_closed_stl(
        data,
        StlOptions {
            facet_limit: Some(patch_map.facet_count),
            facet_patches: None,
        },
    )
    .map_err(|e| EngineError::new(ErrorCode::InvalidDomain, e))?;
    if surface.facets().len() != patch_map.facet_count
        || patch_map.default_patch.role != "wall"
        || patch_map.default_patch.name != "walls"
    {
        return Err(invalid(
            "expected a matching facet count and a walls default patch",
        ));
    }
    let mut patches = vec![3; surface.facets().len()];
    for port in &patch_map.ports {
        let id = match (port.name.as_str(), port.role.as_str()) {
            ("inlet", "inlet") => 1,
            ("outlet", "outlet") => 2,
            _ => return Err(invalid("expected inlet and outlet roles")),
        };
        let [a, b] = port.facet_range;
        let range = patches
            .get_mut(a..b)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| invalid("invalid port range"))?;
        for tag in range {
            if *tag != 3 {
                return Err(invalid("overlapping ports"));
            }
            *tag = id;
        }
    }
    if !patches.contains(&1) || !patches.contains(&2) || !patches.contains(&3) {
        return Err(invalid("expected distinct inlet, outlet and wall patches"));
    }
    Ok((surface, patches))
}

/// The number of cells `generate` builds (AG-15), counted from the inputs:
/// bulk: one tetrahedron per facet and two per edge of a region interface;
/// layered: plus one prism per wall facet and layer; polyhedral: one
/// polyhedron per scaffold node plus one (Voronoi core) or six (median
/// dual) prisms per wall triangle of the scaffold and layer. Inputs
/// `generate` refuses are refused here too.
pub fn cell_count(
    data: &[u8],
    patch_map: &PatchMap,
    plan: &Partition,
    mode: Mode,
    scaffold_data: Option<&[u8]>,
) -> Result<usize, EngineError> {
    let invalid = |s: &str| EngineError::new(ErrorCode::InvalidDomain, s);
    let (surface, patches) = prepare(data, patch_map, plan)?;
    let layers = if matches!(mode, Mode::Bulk) {
        0
    } else {
        plan.layers
    };
    if layers == 0 && !matches!(mode, Mode::Bulk) {
        return Err(invalid(
            "layer count, growth and thickness must be positive",
        ));
    }
    if let (Mode::Polyhedral, Some(config)) = (mode, &plan.volume_scaffold) {
        let bytes = scaffold_data.ok_or_else(|| invalid("missing configured volume scaffold"))?;
        let scaffold =
            scaffold::read(bytes).map_err(|e| EngineError::new(ErrorCode::ValidationFailed, e))?;
        let walls = scaffold.boundary.iter().filter(|(_, p)| *p == 3).count();
        // median dual: six subprisms per wall triangle and layer; Voronoi
        // core: one wedge (or split-cap polyhedron) per wall triangle and
        // layer, one cell per site (the scaffold nodes the tetrahedra use)
        let (per_wall, core) = match config.core {
            CoreKind::MedianDual => (6, scaffold.points.len()),
            CoreKind::Voronoi => {
                let mut used = vec![false; scaffold.points.len()];
                for t in &scaffold.tetrahedra {
                    for &n in t {
                        used[n] = true;
                    }
                }
                (1, used.iter().filter(|&&u| u).count())
            }
        };
        return walls
            .checked_mul(per_wall * layers)
            .and_then(|v| v.checked_add(core))
            .ok_or_else(|| invalid("cell count overflow"));
    }
    // region of each facet, from the partition's facet ranges
    let mut region = vec![usize::MAX; surface.facets().len()];
    for (r, config) in plan.regions.iter().enumerate() {
        for &[a, b] in &config.facets {
            for slot in region
                .get_mut(a..b)
                .ok_or_else(|| invalid("invalid region range"))?
            {
                *slot = r;
            }
        }
    }
    let mut edges = std::collections::HashMap::new();
    for (i, f) in surface.facets().iter().enumerate() {
        let v = f.triangle.vertices;
        for (a, b) in [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
            edges
                .entry((a.min(b), a.max(b)))
                .or_insert_with(Vec::new)
                .push(region[i]);
        }
    }
    let seam = edges
        .values()
        .filter(|r| r.len() == 2 && r[0] != r[1])
        .count();
    let walls = patches.iter().filter(|&&p| p == 3).count();
    surface
        .facets()
        .len()
        .checked_add(2 * seam)
        .and_then(|v| v.checked_add(walls.checked_mul(layers)?))
        .ok_or_else(|| invalid("cell count overflow"))
}

fn voronoi_json(v: &duct::VoronoiStats) -> Value {
    let mut out = json!({
        "sites": v.sites,
        "core_cells": v.core_cells,
        "split_sites": v.split_sites,
        "empty_sites": v.empty_sites,
        "pieces": v.pieces,
        "lloyd_iterations": v.lloyd_iterations,
        "last_relative_shift": v.last_relative_shift,
        "merged_core_faces": v.merged_core_faces,
    });
    let extra = json!({
        "reentrant_splits": v.reentrant_splits,
        "reassigned_pieces": v.reassigned_pieces,
        "detached_pieces": v.detached_pieces,
        "added_sites": v.added_sites,
        "whole_wall_triangles": v.whole_wall_triangles,
        "split_wall_triangles": v.split_wall_triangles,
        "layer_pieces": v.layer_pieces,
    });
    if let (Some(o), Value::Object(e)) = (out.as_object_mut(), extra) {
        o.extend(e);
    }
    out
}
