//! Mesh an annotated STL using a validated, explicitly supplied partition.
mod scaffold;
use crate::handoff::HandoffFormat;
use crate::{
    EngineError, ErrorCode, MeshDirectory,
    artifact::{digest, require_schema, sha256, vector},
};
use mesh_algorithms::duct::{self, Region};
use mesh_algorithms::sizing::surface::SurfaceGeometry;
use mesh_algorithms::{remesh, scaffold as scaffold_build};
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
    /// A Gmsh MSH 2.2 scaffold next to the partition (or `build`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    sha256: Option<String>,
    /// Build the scaffold in Rust from the STL: remesh the wall onto it to a
    /// chord tolerance, then tetrahedralize and refine the interior.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub build: Option<ScaffoldBuild>,
    max_relative_geometry_error: f64,
    /// Largest allowed sampled distance between the scaffold wall and the
    /// STL wall, both ways, in model units.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_wall_deviation: Option<f64>,
    /// Polyhedral core built on the scaffold.
    #[serde(default, skip_serializing_if = "CoreKind::is_default")]
    pub core: CoreKind,
    /// Lloyd iterations of the Voronoi core.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lloyd_iterations: Option<usize>,
    /// How the wall layers meet the Voronoi core.
    #[serde(default, skip_serializing_if = "WallJoinKind::is_default")]
    pub wall_join: WallJoinKind,
}

/// How the wall layers meet the Voronoi core.
#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum WallJoinKind {
    /// A poly-prism column under each core cell and wall plane.
    #[default]
    Columns,
    /// One innermost polyhedron with a split cap per shared wall triangle.
    SplitCap,
}

impl WallJoinKind {
    fn is_default(&self) -> bool {
        *self == Self::default()
    }
}

/// The settings of a scaffold built in Rust (lengths in model units).
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ScaffoldBuild {
    /// Largest distance of a wall facet from the STL (chord tolerance).
    pub tolerance: f64,
    pub min_size: f64,
    pub max_size: f64,
    /// Growth of the spacing per unit distance (wall and interior).
    pub rate: f64,
    /// Wall spacing at the port rims.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rim_size: Option<f64>,
    /// Source edges sharper than this (degrees) are kept as feature lines.
    #[serde(default = "default_feature_angle")]
    pub feature_angle: f64,
    #[serde(default = "default_remesh_iterations")]
    pub iterations: usize,
    /// Wall spacing from the local thickness of gaps and thin walls.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub thickness: Option<ThicknessBuild>,
}

/// Wall spacing = local thickness / `across`, not below `min_size` (which
/// may lie below the build's `min_size`). The thickness is the medial-ball
/// diameter between opposite walls, on the fluid and on the solid side.
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ThicknessBuild {
    pub across: f64,
    pub min_size: f64,
    /// Walls within this angle (degrees) of opposite count as a gap.
    #[serde(default = "default_wedge_angle")]
    pub wedge_angle: f64,
}

fn default_wedge_angle() -> f64 {
    45.0
}

fn default_feature_angle() -> f64 {
    30.0
}
fn default_remesh_iterations() -> usize {
    10
}

/// Limits of a scaffold build.
const MAX_REMESH_ITERATIONS: usize = 50;
const MAX_SCAFFOLD_POINTS: usize = 2_000_000;

/// The polyhedral core of a scaffold mesh.
#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum CoreKind {
    /// Clipped Voronoi diagram restricted to the scaffold, Lloyd-relaxed
    /// (M4). On hold since the solver rules decide (2026-11-25): its
    /// wall-join cells can break them.
    Voronoi,
    /// Vertex-centred median dual of the scaffold (M1b), cells that break
    /// the solver rules split into valid parts. The default, as the dual of
    /// a tet mesh is in production polyhedral meshers.
    #[default]
    MedianDual,
}

impl CoreKind {
    fn is_default(&self) -> bool {
        *self == Self::default()
    }
}

/// Layer cells per wall triangle and layer assumed by the estimate of a
/// Voronoi core with wall columns (measured 0.70 on the elbow, 0.98 on
/// Manifold). The margin also covers the core cells of the conforming wall
/// seeds, one per seeded wall triangle in place of one per wall vertex.
const COLUMN_CELLS_PER_WALL: usize = 2;

/// Rebuilds of a Voronoi core that FoehnCraft would reject.
const VORONOI_RETRIES: usize = 2;

/// Most Lloyd iterations a partition may request.
pub const MAX_LLOYD_ITERATIONS: usize = 200;

impl Partition {
    pub fn scaffold_file(&self, mode: Mode) -> Option<&str> {
        if matches!(mode, Mode::Polyhedral) {
            self.volume_scaffold
                .as_ref()
                .and_then(|s| s.file.as_deref())
        } else {
            None
        }
    }

    /// Whether a build has wall columns under a Voronoi core, whose layer
    /// cell count only the built core decides.
    pub fn wall_columns(&self, mode: Mode) -> bool {
        matches!(mode, Mode::Polyhedral)
            && self.volume_scaffold.as_ref().is_some_and(|s| {
                s.core == CoreKind::Voronoi && s.wall_join == WallJoinKind::Columns
            })
    }
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct RegionConfig {
    kernel: [f64; 3],
    facets: Vec<[usize; 2]>,
}

/// The boundary of a duct: a port file (`schema_version` 1: a `walls`
/// default patch and one inlet and one outlet facet range) or a
/// `koflcraft.patch-map` from the `patch.*` tools, whose patches are one
/// inlet, one outlet and walls (BC-01). A patch map names any facets, so the
/// ports need not be contiguous in the STL.
#[derive(Serialize)]
#[serde(untagged)]
pub enum PatchMap {
    Ports(PortFile),
    Map(crate::patches::PatchMap),
}

impl<'de> Deserialize<'de> for PatchMap {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        use serde::de::Error;
        let value = Value::deserialize(d)?;
        if value.get("format").is_some() {
            serde_json::from_value(value)
                .map(Self::Map)
                .map_err(D::Error::custom)
        } else {
            serde_json::from_value(value)
                .map(Self::Ports)
                .map_err(D::Error::custom)
        }
    }
}

impl PatchMap {
    pub fn facet_count(&self) -> usize {
        match self {
            Self::Ports(p) => p.facet_count,
            Self::Map(m) => m.facet_count,
        }
    }
}

#[derive(Deserialize, Serialize)]
pub struct PortFile {
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

#[allow(clippy::too_many_arguments)]
pub fn generate(
    data: &[u8],
    patch_map: &PatchMap,
    plan: &Partition,
    mode: Mode,
    scaffold_data: Option<&[u8]>,
    dir: &MeshDirectory,
    name: &OsStr,
    handoff_format: HandoffFormat,
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
    let mut scaffold_report: Option<Value> = None;
    let mut voronoi: Option<duct::VoronoiStats> = None;
    let mut source_areas = [0.; 3];
    let scale = surface.bounds().extent().max_abs_component();
    for (f, &tag) in surface.facets().iter().zip(&patches) {
        source_areas[tag - 1] += f.normalized_area * scale * scale;
    }
    let mesh = if let (Mode::Polyhedral, Some(spec)) = (mode, &plan.volume_scaffold) {
        let limit = spec.max_relative_geometry_error;
        if !limit.is_finite() || !(0.0..1.0).contains(&limit) {
            return Err(invalid("invalid scaffold geometry error limit"));
        }
        let (input, built) = load_scaffold(&surface, &patches, spec, scaffold_data)?;
        let deviation = wall_deviation(&surface, &patches, &input);
        if let Some(max) = spec.max_wall_deviation {
            if !max.is_finite() || max <= 0.0 {
                return Err(invalid("invalid wall deviation limit"));
            }
            if deviation[0].max(deviation[1]) > max {
                return Err(validated(format!(
                    "scaffold wall deviates from the STL by {:.4} (limit {max})",
                    deviation[0].max(deviation[1])
                )));
            }
        }
        scaffold_report = Some(json!({
            "source": if built.is_some() {"rust_remesh_delaunay_v1"} else {"file"},
            "build": built,
            "wall_deviation": {"max_scaffold_to_stl": deviation[0], "max_stl_to_scaffold": deviation[1]},
        }));
        let retry_input = input.clone();
        let (mesh, shell) = match spec.core {
            CoreKind::MedianDual => duct::build_volume_dual(input, &steps).map_err(validated)?,
            CoreKind::Voronoi => {
                let lloyd_iterations = spec
                    .lloyd_iterations
                    .unwrap_or(duct::DEFAULT_LLOYD_ITERATIONS);
                if lloyd_iterations > MAX_LLOYD_ITERATIONS {
                    return Err(invalid("too many Lloyd iterations"));
                }
                // The repairs at reentrant wall edges are local and can
                // miss a face on a rare site layout; FoehnCraft's own test
                // (one face per cell pair, S·d > 0) decides, and a failed
                // layout is rebuilt with one more Lloyd iteration.
                let mut input = Some(input);
                let mut attempt = 0;
                loop {
                    let scaffold = match input.take() {
                        Some(s) => s,
                        None => retry_input.clone(),
                    };
                    let (mesh, shell, mut stats) = duct::build_voronoi_core(
                        scaffold,
                        &steps,
                        duct::VoronoiOptions {
                            lloyd_iterations: lloyd_iterations + attempt,
                            wall_join: match spec.wall_join {
                                WallJoinKind::Columns => duct::WallJoin::Columns,
                                WallJoinKind::SplitCap => duct::WallJoin::SplitCap,
                            },
                        },
                    )
                    .map_err(validated)?;
                    let valid =
                        mesh_io::duct::handoff(&mesh).is_ok_and(|h| h.inverted_faces() == 0);
                    if valid || attempt == VORONOI_RETRIES {
                        stats.retries = attempt;
                        voronoi = Some(stats);
                        break (mesh, shell);
                    }
                    attempt += 1;
                }
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
                // f32 STL: a tilted cap is planar to the coordinate rounding
                if face.plane.normal().dot(normal) < 1. - 1e-6
                    || face
                        .triangle
                        .vertices
                        .iter()
                        .any(|&v| (shell.vertices()[v] - origin).dot(normal).abs() > scale * 1e-6)
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
        "scaffold": scaffold_report,
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
    let report_json = serde_json::to_vec_pretty(&report)
        .map_err(|e| EngineError::new(ErrorCode::SerializationFailed, e))?;
    let handoff = mesh_io::duct::handoff(&mesh).map_err(crate::handoff::export_error)?;
    let generator = report["algorithm"].as_str().unwrap_or("duct").to_string();
    // each file goes to disk as it is written, hashed on the way
    let fail = crate::handoff::stream_error;
    let mut bundle = dir
        .bundle(name)
        .map_err(|e| EngineError::new(ErrorCode::IoError, e))?;
    bundle
        .stream("volume.vtu", |w| mesh_io::duct::volume(&mesh, w))
        .map_err(fail)?;
    bundle
        .stream("mesh.msh", |w| {
            mesh_io::duct::gmsh(&mesh, mesh_io::GmshVersion::V4_1, w)
        })
        .map_err(fail)?;
    bundle
        .stream("mesh.vtm", |w| w.write_all(mesh_io::channel::MULTIBLOCK))
        .map_err(fail)?;
    for (i, file) in ["inlet.vtu", "outlet.vtu", "walls.vtu"]
        .into_iter()
        .enumerate()
    {
        bundle
            .stream(file, |w| mesh_io::duct::patch(&mesh, i + 1, w))
            .map_err(fail)?;
    }
    bundle
        .stream("report.json", |w| w.write_all(&report_json))
        .map_err(fail)?;
    bundle
        .stream(handoff_format.file_name(), |w| {
            handoff_format.write(&handoff, &[], w)
        })
        .map_err(fail)?;
    let evidence = crate::handoff::evidence(&handoff, &generator, bundle.files())?;
    bundle
        .stream(crate::handoff::MESH_EVIDENCE, |w| w.write_all(&evidence))
        .map_err(fail)?;
    bundle.finish();
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
    require_schema(plan.schema_version)?;
    if sha256(data) != plan.source_sha256 {
        return Err(invalid(
            "patch map or partition digest differs from input STL",
        ));
    }
    surface_tags(data, patch_map)
}

/// The surface and per-facet patch ids (1 inlet, 2 outlet, 3 wall) of a
/// checked STL and patch map.
fn surface_tags(
    data: &[u8],
    patch_map: &PatchMap,
) -> Result<(mesh_geometry::ClosedSurface, Vec<usize>), EngineError> {
    let invalid = |s: &str| EngineError::new(ErrorCode::InvalidDomain, s);
    let source_digest = sha256(data);
    let (map_digest, facet_count) = match patch_map {
        PatchMap::Ports(p) => {
            require_schema(p.schema_version)?;
            (&p.source_sha256, p.facet_count)
        }
        PatchMap::Map(m) => (&m.source_sha256, m.facet_count),
    };
    if source_digest != *map_digest {
        return Err(invalid(
            "patch map or partition digest differs from input STL",
        ));
    }
    let surface = import_closed_stl(
        data,
        StlOptions {
            facet_limit: Some(facet_count),
            facet_patches: None,
        },
    )
    .map_err(|e| EngineError::new(ErrorCode::InvalidDomain, e))?;
    if surface.facets().len() != facet_count {
        return Err(invalid("expected a matching facet count"));
    }
    let patches = match patch_map {
        PatchMap::Ports(p) => port_tags(p, facet_count)?,
        PatchMap::Map(m) => map_tags(m)?,
    };
    if !patches.contains(&1) || !patches.contains(&2) || !patches.contains(&3) {
        return Err(invalid("expected distinct inlet, outlet and wall patches"));
    }
    Ok((surface, patches))
}

/// What `setup` may change; lengths in model units.
#[derive(Clone, Copy, Debug)]
pub struct SetupOptions {
    /// Multiplies every derived length; below 1 is finer.
    pub size_factor: f64,
    pub layers: Option<usize>,
    pub first_thickness: Option<f64>,
    pub growth: Option<f64>,
}

impl Default for SetupOptions {
    fn default() -> Self {
        Self {
            size_factor: 1.0,
            layers: None,
            first_thickness: None,
            growth: None,
        }
    }
}

/// The partition of a polyhedral duct for an STL and its patch map: a
/// scaffold built in Rust and wall layers, every length a fixed fraction of
/// the smaller port's equivalent radius r (the fractions of the Manifold
/// fixture, r = 25: chord tolerance r/500, wall sizes r/38 to r/4, rim size
/// r/29, first layer r/312, growth 1.4, three layers), times `size_factor`.
/// The regions are empty: the scaffold route does not use them. Returns the
/// partition and r.
pub fn setup(
    data: &[u8],
    patch_map: &PatchMap,
    options: SetupOptions,
) -> Result<(Partition, f64), EngineError> {
    let invalid = |s: &str| EngineError::new(ErrorCode::InvalidArguments, s);
    let k = options.size_factor;
    if !k.is_finite() || k <= 0.0 {
        return Err(invalid("size_factor must be positive"));
    }
    let (surface, patches) = surface_tags(data, patch_map)?;
    let scale = surface.bounds().extent().max_abs_component();
    let mut areas = [0.0; 2];
    for (f, &tag) in surface.facets().iter().zip(&patches) {
        if tag < 3 {
            areas[tag - 1] += f.normalized_area * scale * scale;
        }
    }
    let r = (areas[0].min(areas[1]) / std::f64::consts::PI).sqrt();
    let h = r * k;
    let layers = options.layers.unwrap_or(3);
    let first_thickness = options.first_thickness.unwrap_or(h / 312.5);
    let growth = options.growth.unwrap_or(1.4);
    if layers == 0
        || !(first_thickness.is_finite() && first_thickness > 0.0)
        || !(growth.is_finite() && growth > 0.0)
    {
        return Err(invalid(
            "layers, first_thickness and growth must be positive",
        ));
    }
    // The wall cannot follow the STL more closely than the STL's own facets
    // follow the shape they were made from.
    let facet_error = stl_facet_error(&surface, &patches);
    let tolerance = (h / 500.0).max(1.25 * facet_error);
    let round = |v: f64| {
        let digits = 10f64.powi(3 - v.abs().log10().ceil() as i32);
        (v * digits).round() / digits
    };
    let plan = Partition {
        schema_version: 1,
        source_sha256: sha256(data),
        regions: Vec::new(),
        volume_scaffold: Some(ScaffoldConfig {
            file: None,
            sha256: None,
            build: Some(ScaffoldBuild {
                tolerance: round(tolerance),
                min_size: round(h / 38.5),
                max_size: round(h / 4.17),
                rate: 0.4,
                rim_size: Some(round(h / 29.4)),
                feature_angle: default_feature_angle(),
                iterations: default_remesh_iterations(),
                thickness: Some(ThicknessBuild {
                    across: 1.0,
                    min_size: round(h / 154.0),
                    wedge_angle: default_wedge_angle(),
                }),
            }),
            max_relative_geometry_error: 0.002,
            max_wall_deviation: Some(round(1.5 * tolerance)),
            core: CoreKind::default(),
            lloyd_iterations: None,
            wall_join: WallJoinKind::default(),
        }),
        first_thickness: round(first_thickness),
        growth,
        layers,
    };
    Ok((plan, r))
}

/// The chord error of the STL's own facets on its smooth walls: where two
/// wall facets meet at a turn β below the feature angle, facets of width w
/// across the edge cut the smooth wall they came from by about β·w/8 (the
/// sagitta of a chord of length w that turns by β). The value is the
/// length-weighted 99th percentile over these edges, so a few odd facets do
/// not set it.
pub fn stl_facet_error(surface: &mesh_geometry::ClosedSurface, patches: &[usize]) -> f64 {
    let g = SurfaceGeometry::new(surface);
    let limit = default_feature_angle().to_radians();
    let mut items: Vec<(f64, f64)> = Vec::new();
    for e in &g.edges {
        let [f, k] = e.facets;
        if patches[f] != 3 || patches[k] != 3 || e.angle >= limit || e.angle < 1e-3 {
            continue;
        }
        let [a, b] = e.vertices;
        let length = (g.vertices[a] - g.vertices[b]).length();
        if length <= 0.0 {
            continue;
        }
        let width = (g.areas[f] + g.areas[k]) / length;
        items.push((e.angle * width / 8.0, length));
    }
    items.sort_by(|a, b| a.0.total_cmp(&b.0));
    let total: f64 = items.iter().map(|i| i.1).sum();
    let mut seen = 0.0;
    for (error, length) in items {
        seen += length;
        if seen >= 0.99 * total {
            return error;
        }
    }
    0.0
}

/// Patch ids (1 inlet, 2 outlet, 3 wall) of a port file.
fn port_tags(map: &PortFile, facet_count: usize) -> Result<Vec<usize>, EngineError> {
    let invalid = |s: &str| EngineError::new(ErrorCode::InvalidDomain, s);
    if map.default_patch.role != "wall" || map.default_patch.name != "walls" {
        return Err(invalid("expected a walls default patch"));
    }
    let mut patches = vec![3; facet_count];
    for port in &map.ports {
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
    Ok(patches)
}

/// Patch ids (1 inlet, 2 outlet, 3 wall) of a patch map: exactly one patch
/// with role inlet, one with role outlet, every other patch a wall. The
/// duct mesh names them `inlet`, `outlet` and `walls`.
fn map_tags(map: &crate::patches::PatchMap) -> Result<Vec<usize>, EngineError> {
    use crate::patches::PatchRole;
    let bad = |s: String| EngineError::new(ErrorCode::InvalidPatchMap, s);
    if map.format != crate::patches::PATCH_MAP_FORMAT
        || map.format_version != crate::patches::PATCH_MAP_FORMAT_VERSION
        || map.facet_patch.len() != map.facet_count
    {
        return Err(bad(format!(
            "expected a {} {} map with one patch per facet",
            crate::patches::PATCH_MAP_FORMAT,
            crate::patches::PATCH_MAP_FORMAT_VERSION
        )));
    }
    let mut tag = Vec::with_capacity(map.patches.len());
    let (mut inlets, mut outlets) = (0, 0);
    for p in &map.patches {
        tag.push(match p.role {
            PatchRole::Inlet => {
                inlets += 1;
                1
            }
            PatchRole::Outlet => {
                outlets += 1;
                2
            }
            PatchRole::Wall => 3,
            other => {
                return Err(bad(format!(
                    "a duct takes inlet, outlet and wall patches; {} has role {other:?}",
                    p.name
                )));
            }
        });
    }
    if inlets != 1 || outlets != 1 {
        return Err(bad(format!(
            "a duct needs exactly one inlet and one outlet patch, the map has {inlets} and {outlets}"
        )));
    }
    map.facet_patch
        .iter()
        .enumerate()
        .map(|(f, &p)| {
            tag.get(p).copied().ok_or_else(|| {
                bad(format!(
                    "facet {f} has patch index {p}; the map has {} patches",
                    tag.len()
                ))
            })
        })
        .collect()
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
        let (scaffold, _) = load_scaffold(&surface, &patches, config, scaffold_data)?;
        let walls = scaffold.boundary.iter().filter(|(_, p)| *p == 3).count();
        // median dual: six subprisms per wall triangle and layer; Voronoi
        // core: one wedge (or split-cap polyhedron) per wall triangle and
        // layer, one cell per site (the scaffold nodes the tetrahedra use).
        // Wall columns follow the built core: two per wall triangle and
        // layer bound the fixtures from above (elbow 0.8, Manifold 2.1).
        let (per_wall, core) = match config.core {
            CoreKind::MedianDual => (6, scaffold.points.len()),
            CoreKind::Voronoi => {
                let mut used = vec![false; scaffold.points.len()];
                for t in &scaffold.tetrahedra {
                    for &n in t {
                        used[n] = true;
                    }
                }
                let per_wall = match config.wall_join {
                    WallJoinKind::Columns => COLUMN_CELLS_PER_WALL,
                    WallJoinKind::SplitCap => 1,
                };
                (per_wall, used.iter().filter(|&&u| u).count())
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

/// The scaffold of a polyhedral run: the configured file (digest checked)
/// or, with `build`, the STL wall remeshed onto itself to the chord
/// tolerance and tetrahedralized in Rust; the build's statistics.
fn load_scaffold(
    surface: &mesh_geometry::ClosedSurface,
    patches: &[usize],
    spec: &ScaffoldConfig,
    data: Option<&[u8]>,
) -> Result<(duct::Scaffold, Option<Value>), EngineError> {
    let invalid = |s: &str| EngineError::new(ErrorCode::InvalidDomain, s);
    let validated = |e: String| EngineError::new(ErrorCode::ValidationFailed, e);
    match (&spec.file, &spec.build) {
        (Some(_), None) => {
            let bytes = data.ok_or_else(|| invalid("missing configured volume scaffold"))?;
            if spec.sha256.as_deref() != Some(sha256(bytes).as_str()) {
                return Err(invalid("volume scaffold digest mismatch"));
            }
            Ok((scaffold::read(bytes).map_err(validated)?, None))
        }
        (None, Some(b)) => {
            if spec.sha256.is_some() {
                return Err(invalid("a built scaffold takes no digest"));
            }
            let positive = |v: f64| v.is_finite() && v > 0.0;
            if !(positive(b.tolerance)
                && positive(b.min_size)
                && positive(b.max_size)
                && b.max_size >= b.min_size
                && positive(b.rate)
                && b.rim_size.is_none_or(positive)
                && b.feature_angle > 0.0
                && b.feature_angle < 180.0
                && b.iterations <= MAX_REMESH_ITERATIONS
                && b.thickness.is_none_or(|t| {
                    positive(t.across)
                        && positive(t.min_size)
                        && t.min_size <= b.max_size
                        && t.wedge_angle > 0.0
                        && t.wedge_angle < 90.0
                }))
            {
                return Err(invalid("invalid scaffold build settings"));
            }
            let start = Instant::now();
            let triangles: Vec<[usize; 3]> = surface
                .facets()
                .iter()
                .map(|f| f.triangle.vertices)
                .collect();
            let keep: Vec<bool> = patches.iter().map(|&p| p != 3).collect();
            let wall = remesh::remesh(
                surface.vertices(),
                &triangles,
                &keep,
                remesh::RemeshOptions {
                    tolerance: b.tolerance,
                    min_size: b.min_size,
                    max_size: b.max_size,
                    rate: b.rate,
                    feature_angle: b.feature_angle.to_radians(),
                    iterations: b.iterations,
                    rim_size: b.rim_size,
                    thickness: b.thickness.map(|t| remesh::ThicknessSizing {
                        across: t.across,
                        min_size: t.min_size,
                        wedge_angle: t.wedge_angle.to_radians(),
                    }),
                },
            )
            .map_err(validated)?;
            let remesh_seconds = start.elapsed().as_secs_f64();
            if mesh_algorithms::trace_enabled() {
                let w = wall.stats;
                eprintln!(
                    "kofl trace: wall remesh {remesh_seconds:.1} s, {} rounds: {} triangles, {} sharp-edge splits, {} merges, edges {:.3e} to {:.3e}, min gap {:.3e}, min wall {:.3e}, {} thickness pieces at the floor",
                    w.rounds,
                    wall.triangles.len(),
                    w.feature_ball_splits,
                    w.feature_ball_merges,
                    w.min_edge,
                    w.max_edge,
                    w.min_gap,
                    w.min_wall,
                    w.thickness_clamped
                );
            }
            let tags: Vec<usize> = wall
                .source
                .iter()
                .map(|s| s.map_or(3, |f| patches[f]))
                .collect();
            let (scaffold, st) = scaffold_build::build(
                &wall.vertices,
                &wall.triangles,
                &tags,
                scaffold_build::ScaffoldOptions {
                    max_size: b.max_size,
                    rate: b.rate,
                    max_points: MAX_SCAFFOLD_POINTS,
                    ..Default::default()
                },
            )
            .map_err(validated)?;
            let w = wall.stats;
            Ok((
                scaffold,
                Some(json!({
                    "settings": b,
                    "wall_triangles": tags.iter().filter(|&&t| t == 3).count(),
                    "feature_edges": w.feature_edges,
                    "feature_corners": w.corners,
                    "wall_edge_range": [w.min_edge, w.max_edge],
                    "deviation_splits": w.deviation_splits,
                    "remesh_rounds": w.rounds,
                    "short_edges": w.short_edges,
                    "sharp_edge_splits": w.feature_ball_splits,
                    "min_gap": w.min_gap.is_finite().then_some(w.min_gap),
                    "min_wall": w.min_wall.is_finite().then_some(w.min_wall),
                    "thickness_clamped": w.thickness_clamped,
                    "surface_points": st.surface_points,
                    "recovery_points": st.recovery_points,
                    "interior_points": st.interior_points,
                    "refused_points": st.refused,
                    "tetrahedra": st.tetrahedra,
                    "remesh_seconds": remesh_seconds,
                    "tetrahedralize_seconds": start.elapsed().as_secs_f64() - remesh_seconds,
                })),
            ))
        }
        _ => Err(invalid(
            "volume scaffold needs exactly one of file or build",
        )),
    }
}

/// The largest sampled distance from the scaffold wall to the STL wall and
/// from the STL wall to the scaffold wall (each facet at its centroid,
/// edge midpoints and three interior points).
fn wall_deviation(
    surface: &mesh_geometry::ClosedSurface,
    patches: &[usize],
    scaffold: &duct::Scaffold,
) -> [f64; 2] {
    use mesh_algorithms::sizing::surface::Primitive;
    let stl: Vec<[mesh_core::geometry::Vec3; 3]> = surface
        .facets()
        .iter()
        .zip(patches)
        .filter(|(_, p)| **p == 3)
        .map(|(f, _)| f.triangle.vertices.map(|v| surface.vertices()[v]))
        .collect();
    let wall: Vec<[mesh_core::geometry::Vec3; 3]> = scaffold
        .boundary
        .iter()
        .filter(|(_, p)| *p == 3)
        .map(|(t, _)| t.map(|v| scaffold.points[v]))
        .collect();
    let one_way = |from: &[[mesh_core::geometry::Vec3; 3]],
                   to: &[[mesh_core::geometry::Vec3; 3]]| {
        let set = remesh::ClosestSet::new(to.iter().map(|&t| Primitive::Triangle(t)).collect());
        from.iter()
            .flat_map(|&t| remesh::samples(t))
            .map(|x| set.closest(x).map_or(f64::INFINITY, |(_, d)| d))
            .fold(0.0, f64::max)
    };
    [one_way(&wall, &stl), one_way(&stl, &wall)]
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
        "wall_join": match v.wall_join {
            duct::WallJoin::Columns => "columns",
            duct::WallJoin::SplitCap => "split_cap",
        },
    });
    let extra = json!({
        "reentrant_splits": v.reentrant_splits,
        "torn_contacts": v.torn_contacts,
        "wall_seeds": v.wall_seeds,
        "seed_flips": v.seed_flips,
        "unseeded_wall_triangles": v.unseeded_wall_triangles,
        "reassigned_pieces": v.reassigned_pieces,
        "detached_pieces": v.detached_pieces,
        "added_sites": v.added_sites,
        "merged_slivers": v.merged_slivers,
        "retries": v.retries,
        "whole_wall_triangles": v.whole_wall_triangles,
        "split_wall_triangles": v.split_wall_triangles,
        "layer_pieces": v.layer_pieces,
        "layer_columns": v.layer_columns,
        "merged_column_units": v.merged_column_units,
        "split_columns": v.split_columns,
        "absorbed_columns": v.absorbed_columns,
        "collapsed_pieces": v.collapsed_pieces,
        "welded_vertices": v.welded_vertices,
        "merged_films": v.merged_films,
        "turned_wall_rings": v.turned_wall_rings,
        "merged_wall_rings": v.merged_wall_rings,
    });
    if let (Some(o), Value::Object(e)) = (out.as_object_mut(), extra) {
        o.extend(e);
    }
    out
}
