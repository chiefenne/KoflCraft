//! Mesh an annotated STL using a validated, explicitly supplied partition.
mod scaffold;
use crate::handoff::HandoffFormat;
use crate::{
    EngineError, ErrorCode, MeshDirectory,
    artifact::{digest, require_schema, sha256, vector},
};
use mesh_algorithms::duct::{self, Region};
use mesh_algorithms::hexcore::{self, HexCore};
use mesh_algorithms::sizing::surface::SurfaceGeometry;
use mesh_algorithms::{remesh, scaffold as scaffold_build, shell_surface};
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
    /// Names of wall patches that grow no layers (planar patches only): the
    /// layers of the walls around run straight into them.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub no_layers: Vec<String>,
    /// Names of planar wall patches whose layers yield: they keep their
    /// layers, and every wall that crosses them (at a right angle) runs its
    /// own layers straight through down to them, with no mitre.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub layers_yield: Vec<String>,
    /// Zero-thickness shells inside the fluid (baffles, plates, fins), cut
    /// into the finished mesh as walls with two sides.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub shells: Vec<ShellConfig>,
}

/// A planar zero-thickness shell: a convex polygon, its corners
/// counter-clockwise about the side named `name` (the other side is
/// `<name>_back`). A shell that should span the duct may reach past the
/// wall; the part outside the fluid is ignored.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct ShellConfig {
    pub name: String,
    pub polygon: Vec<[f64; 3]>,
    /// Cell size on the shell. Default: the build's `max_size` or an
    /// eighth of the shell's narrow width inside the STL's box, whichever
    /// is smaller (about eight cells across).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub size: Option<f64>,
    /// Cell size along the shell's free edges. Default: half the shell's
    /// size.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub edge_size: Option<f64>,
}

/// The sizes a shell is built to, and whether each was set or chosen.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellSizes {
    pub size: f64,
    pub edge_size: f64,
    pub size_set: bool,
    pub edge_size_set: bool,
}

impl ShellConfig {
    /// The shell's sizes in a build (see the fields' defaults); `lo`, `hi`
    /// bound the STL.
    pub fn sizes(
        &self,
        build: &ScaffoldBuild,
        lo: mesh_core::geometry::Vec3,
        hi: mesh_core::geometry::Vec3,
    ) -> Result<ShellSizes, EngineError> {
        let invalid = |m: String| EngineError::new(ErrorCode::InvalidArguments, m);
        let positive = |v: f64| v.is_finite() && v > 0.0;
        let polygon: Vec<mesh_core::geometry::Vec3> =
            self.polygon.iter().map(|&p| vector(p)).collect();
        let size = match self.size {
            Some(v) if positive(v) => v.min(build.max_size),
            Some(_) => {
                return Err(invalid(format!(
                    "shell {}: size must be positive",
                    self.name
                )));
            }
            None => shell_surface::narrow_width(&polygon, lo, hi)
                .map_err(|e| invalid(format!("shell {}: {e}", self.name)))?
                .map_or(build.max_size, |w| (w / 8.0).min(build.max_size)),
        };
        let edge_size = match self.edge_size {
            Some(v) if positive(v) => v.min(build.max_size),
            Some(_) => {
                return Err(invalid(format!(
                    "shell {}: edge_size must be positive",
                    self.name
                )));
            }
            None => size / 2.0,
        };
        Ok(ShellSizes {
            size,
            edge_size,
            size_set: self.size.is_some(),
            edge_size_set: self.edge_size.is_some(),
        })
    }
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
    /// Settings of the hex core.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hex_core: Option<HexCoreConfig>,
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
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
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
    /// Local sizes from datum objects and coordinate regions (GRD-02,
    /// GRD-03): wall and interior spacing follow min(own sizes, grading).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub grading: Option<Grading>,
}

/// Datum objects and coordinate regions of a duct build, the `sources` and
/// `frames` of an STL size field (kinds `datum` and `region` only). The
/// field's growth limit is the build's `rate` and its largest size the
/// build's `max_size`.
#[derive(Clone, Debug, Default, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Grading {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub frames: Vec<crate::sizing::FrameConfig>,
    pub sources: Vec<crate::sizing::SourceConfig>,
}

impl Grading {
    /// The size field of the sources within the build's limits.
    pub fn field(
        &self,
        build: &ScaffoldBuild,
    ) -> Result<mesh_algorithms::sizing::SizeField, EngineError> {
        let config = crate::sizing::SizeFieldConfig {
            min_size: self
                .sources
                .iter()
                .filter_map(|s| s.size())
                .fold(build.min_size, f64::min),
            max_size: build.max_size,
            gradation: 1.0 + build.rate,
            relaxation_iterations: None,
            frames: self.frames.clone(),
            sources: self.sources.clone(),
        };
        config.explicit_field()
    }
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
    /// Octree cubes in the bulk, the median dual between them and the wall
    /// layers ("hex core with layers", prototype; research/hex-core.md).
    /// Needs a built scaffold.
    HexCore,
}

/// Settings of the hex core (`"core": "hex_core"`).
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct HexCoreConfig {
    /// Cube edge over the local scaffold spacing (default 1).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub size_factor: Option<f64>,
    /// Gap between the layer stack and the cubes, in cube edges (default 1).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub margin: Option<f64>,
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
    // KOFL_TRACE=1: the wall-clock seconds of each stage on standard error
    let mut stage = Instant::now();
    let mut lap = |name: &str| {
        if mesh_algorithms::trace_enabled() {
            eprintln!(
                "kofl trace: stage {name}: {:.1} s",
                stage.elapsed().as_secs_f64()
            );
        }
        stage = Instant::now();
    };
    let (surface, patches, mut names) = prepare(data, patch_map, plan)?;
    // the STL's own patches; shells built into the scaffold follow them
    let stl_patches = names.len();
    let scaffold_shells = !plan.shells.is_empty()
        && matches!(mode, Mode::Polyhedral)
        && scaffold_data.is_none()
        && plan
            .volume_scaffold
            .as_ref()
            .is_some_and(|s| s.build.is_some() && matches!(s.core, CoreKind::MedianDual));
    if scaffold_shells {
        if !plan.layers_yield.is_empty() {
            return Err(invalid(
                "shells and layers_yield in one build are not supported yet",
            ));
        }
        for s in &plan.shells {
            for name in [s.name.clone(), format!("{}_back", s.name)] {
                if name.is_empty() || names.iter().any(|(n, _)| *n == name) {
                    return Err(invalid(&format!(
                        "shell patch name {name:?} is empty or already used"
                    )));
                }
                names.push((name, "wall".to_string()));
            }
        }
    }
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
    // Walls whose layers yield: build on a domain with their planes moved
    // in by the layer height, then extrude the slab back out.
    let full_source = surface.clone();
    let mut yielding_ids: Vec<usize> = Vec::new();
    let mut yield_planes = Vec::new();
    let surface = if plan.layers_yield.is_empty() {
        surface
    } else {
        let core_ok = matches!(mode, Mode::Polyhedral)
            && plan
                .volume_scaffold
                .as_ref()
                .is_some_and(|s| s.build.is_some() && matches!(s.core, CoreKind::MedianDual));
        if !core_ok || scaffold_data.is_some() {
            return Err(invalid(
                "layers_yield needs the polyhedral mode with a built scaffold and the median dual core",
            ));
        }
        yielding_ids = (1..=names.len())
            .filter(|&id| plan.layers_yield.contains(&names[id - 1].0))
            .collect();
        // flush openings inside a yielding wall yield with it
        let geometry = SurfaceGeometry::new(&surface);
        for id in duct::WALL..=names.len() {
            if names[id - 1].1 != "wall"
                && geometry.edges.iter().any(|e| {
                    let [f, g] = e.facets;
                    let (a, b) = (duct::patch_index(patches[f]), duct::patch_index(patches[g]));
                    (a == id && yielding_ids.contains(&b)) || (b == id && yielding_ids.contains(&a))
                })
            {
                yielding_ids.push(id);
            }
        }
        let yielding: Vec<bool> = patches
            .iter()
            .map(|&p| duct::is_slip(p) && yielding_ids.contains(&duct::patch_index(p)))
            .collect();
        let (short, planes) = duct::shorten_yielding(&surface, &yielding, steps.iter().sum())
            .map_err(|e| invalid(&e))?;
        yield_planes = planes;
        short
    };
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
    let mut shell_report: Vec<(
        ShellConfig,
        shell_surface::ShellSurfaceStats,
        duct::CrackStats,
        Option<ShellSizes>,
    )> = Vec::new();
    let mut grading_field = None;
    let mut voronoi: Option<duct::VoronoiStats> = None;
    let mut source_areas = vec![0.; names.len()];
    let scale = surface.bounds().extent().max_abs_component();
    for (f, &tag) in surface.facets().iter().zip(&patches) {
        source_areas[duct::patch_index(tag) - 1] += f.normalized_area * scale * scale;
    }
    let mesh = if let (Mode::Polyhedral, Some(spec)) = (mode, &plan.volume_scaffold) {
        let limit = spec.max_relative_geometry_error;
        if !limit.is_finite() || !(0.0..1.0).contains(&limit) {
            return Err(invalid("invalid scaffold geometry error limit"));
        }
        let shell_configs: &[ShellConfig] = if scaffold_shells { &plan.shells } else { &[] };
        let (input, built, hex, shell_surfaces) = load_scaffold(
            &surface,
            &patches,
            spec,
            scaffold_data,
            steps.iter().sum(),
            shell_configs,
        )?;
        lap("scaffold");
        let deviation = wall_deviation(&surface, &patches, &input);
        lap("wall deviation");
        // the shells' two sides become boundary faces (patches without
        // layers: the wall layers run into a shell as into a no_layers wall)
        let sides: Vec<(usize, usize)> = (0..shell_surfaces.len())
            .map(|k| {
                (
                    stl_patches + 1 + 2 * k + duct::SLIP,
                    stl_patches + 2 + 2 * k + duct::SLIP,
                )
            })
            .collect();
        let (input, crack) = duct::crack_shells(input, &sides).map_err(validated)?;
        shell_report = shell_surfaces
            .iter()
            .enumerate()
            .map(|(k, st)| {
                let s = &plan.shells[k];
                let sizes = spec
                    .build
                    .as_ref()
                    .map(|b| s.sizes(b, surface.bounds().min(), surface.bounds().max()))
                    .transpose()?;
                Ok((s.clone(), *st, crack, sizes))
            })
            .collect::<Result<_, EngineError>>()?;
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
        grading_field = match &spec.build {
            Some(b) => b.grading.as_ref().map(|g| g.field(b)).transpose()?,
            None => None,
        };
        scaffold_report = Some(json!({
            "source": if built.is_some() {"rust_remesh_delaunay_v1"} else {"file"},
            "build": built,
            "wall_deviation": {"max_scaffold_to_stl": deviation[0], "max_stl_to_scaffold": deviation[1]},
        }));
        let retry_input = input.clone();
        let (mesh, shell) = match spec.core {
            CoreKind::MedianDual => duct::build_volume_dual(input, &steps).map_err(validated)?,
            CoreKind::HexCore => {
                let core = hex
                    .as_ref()
                    .ok_or_else(|| invalid("the hex core needs a built scaffold"))?;
                duct::build_volume_dual_hex(input, &steps, core).map_err(validated)?
            }
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
        lap("core, layers and checks");
        let mut mesh = mesh;
        unmark_slip(&mut mesh, names.len());
        scaffold_volume = Some(shell.volume());
        let volume_error = shell.volume() / surface.volume() - 1.;
        if volume_error.abs() > limit {
            return Err(invalid(&format!(
                "remeshed volume differs from the STL by {:.3}% (limit {:.3}%); a smaller build tolerance or a larger max_relative_geometry_error helps",
                100. * volume_error,
                100. * limit
            )));
        }
        if let Some((i, (a, b))) = source_areas
            .iter()
            .zip(&mesh.measurements.areas)
            .enumerate()
            .find(|&(_, (&a, &b))| a > 0. && (b / a - 1.).abs() > limit)
        {
            return Err(invalid(&format!(
                "remeshed area of patch {:?} differs from the STL by {:.3}% (limit {:.3}%); a smaller build tolerance or a larger max_relative_geometry_error helps",
                names[i].0,
                100. * (b / a - 1.),
                100. * limit
            )));
        }
        for tag in 1..duct::WALL {
            let Some((f, _)) = surface
                .facets()
                .iter()
                .zip(&patches)
                .find(|(_, p)| **p == tag)
            else {
                continue; // no cap: the opening is flush and layered
            };
            let origin = surface.vertices()[f.triangle.vertices[0]];
            let normal = f.plane.normal();
            let id = tag.to_string();
            for face in shell
                .facets()
                .iter()
                .filter(|f| shell.patch_names()[f.triangle.patch.0] == id)
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
        if !yield_planes.is_empty() {
            let shell_patches: Vec<usize> = shell
                .facets()
                .iter()
                .map(|f| {
                    shell.patch_names()[f.triangle.patch.0]
                        .parse::<usize>()
                        .map(duct::patch_index)
                        .unwrap_or(0)
                })
                .collect();
            let full = duct::extrude_yielding(
                &mut mesh,
                &shell,
                &shell_patches,
                &yield_planes,
                &yielding_ids,
                &steps,
            )
            .map_err(validated)?;
            scaffold_volume = Some(full.volume());
        }
        mesh
    } else if matches!(mode, Mode::Polyhedral) {
        duct::build_polyhedral(&surface, &patches, &regions, &steps).map_err(validated)?
    } else {
        duct::build(&surface, &patches, &regions, &steps).map_err(validated)?
    };
    let mut mesh = mesh;
    mesh.patches = names.clone();
    unmark_slip(&mut mesh, names.len());
    let shells = if scaffold_shells {
        scaffold_shell_report(&mesh, stl_patches, &shell_report)
    } else {
        // the exact cut of the finished cells (convex cores)
        imprint_shells(&mut mesh, &plan.shells)?
    };
    let unlayered: Vec<usize> = patches
        .iter()
        .filter(|&&p| duct::is_slip(p))
        .map(|&p| duct::patch_index(p))
        .filter(|id| !yielding_ids.contains(id))
        .collect();
    // the report names the areas of the full STL, also when yielding walls
    // were built on a shortened domain
    let mut source_areas = vec![0.; names.len()];
    for (f, &tag) in full_source.facets().iter().zip(&patches) {
        source_areas[duct::patch_index(tag) - 1] += f.normalized_area * scale * scale;
    }
    let m = &mesh.measurements;
    // the walls of the STL; shells are reported on their own
    let walls_of = |a: &[f64]| a.iter().take(stl_patches).skip(duct::WALL - 1).sum::<f64>();
    let patch_list: Vec<Value> = names
        .iter()
        .enumerate()
        .filter(|&(i, _)| source_areas.get(i).is_some_and(|&a| a > 0.) || i >= stl_patches)
        .map(|(i, (name, role))| {
            json!({"id": i + 1, "name": name, "role": role,
                   "layered": duct::is_wall(i + 1) && !unlayered.contains(&(i + 1)) && i < stl_patches,
                   "stl_area": source_areas.get(i).copied().unwrap_or(0.), "area": m.areas.get(i).copied().unwrap_or(0.)})
        })
        .collect();
    let mut tetrahedra = 0;
    let mut prisms = 0;
    let mut polyhedra = 0;
    // wall layer cells written as polyhedra: polygon prisms, wall columns
    let mut layer_polyhedra = 0;
    let mut face_range = [usize::MAX, 0];
    for cell in &mesh.cells {
        match &cell.element {
            duct::Element::Tet(_) => tetrahedra += 1,
            duct::Element::Wedge(_) => prisms += 1,
            duct::Element::Polyhedron(p) => {
                polyhedra += 1;
                layer_polyhedra += usize::from(cell.layer > 0);
                face_range[0] = face_range[0].min(p.faces.len());
                face_range[1] = face_range[1].max(p.faces.len());
            }
        }
    }
    let grading = grading_field
        .as_ref()
        .map(|field| grading_evidence(field, &mesh));
    let report = json!({
        "schema_version": 1,
        "grading": grading,
        "shells": shells,
        "algorithm": if voronoi.is_some() {"volume_clipped_voronoi_cvt_v1"} else if plan.volume_scaffold.as_ref().is_some_and(|s| s.core == CoreKind::HexCore) {"volume_hex_core_v1"} else if scaffold_volume.is_some() {"volume_median_dual_v2"} else if polyhedra > 0 {"stl_partitioned_dual_v1"} else {"stl_partitioned_star_v1"},
        "mode": mode,
        "source_sha256": source_digest,
        "config_sha256": digest(&(patch_map, plan, mode))?,
        "regions": if scaffold_volume.is_some(){1}else{regions.len()},
        "scaffold_volume": scaffold_volume,
        "scaffold": scaffold_report,
        "stl_patch_areas": {"inlet":source_areas[0],"outlet":source_areas[1],"walls":walls_of(&source_areas)},
        "patches": patch_list,
        "source_facets": surface.facets().len(),
        "tetrahedra": tetrahedra,
        "polyhedra": polyhedra,
        "layer_polyhedra": layer_polyhedra,
        "core_sphericity": m.core_sphericity.map(|[min,p05,median]| json!({"minimum":min,"p05":p05,"median":median})),
        "polyhedron_face_range": (polyhedra > 0).then_some(face_range),
        "prism_core_faces":m.prism_core_faces,
        "prisms": prisms,
        "source_volume": full_source.volume(),
        "volume": m.volume,
        "minimum_cell_volume": m.minimum_volume,
        "minimum_tet_mean_ratio": (tetrahedra > 0).then_some(m.minimum_tet_mean_ratio),
        "tetrahedra_below_mean_ratio_0_1": m.tetrahedra_below_mean_ratio_0_1,
        "minimum_wedge_jacobian": m.minimum_wedge_jacobian,
        "minimum_kernel_clearance": m.minimum_kernel_clearance,
        "patch_areas": {"inlet": m.areas[0], "outlet": m.areas[1], "walls": walls_of(&m.areas)},
        "port_areas_by_layer": {"inlet": m.port_layer_areas[0], "outlet": m.port_layer_areas[1]},
        "thicknesses": steps,
        "normal_thickness_ranges": m.normal_thickness_ranges,
        "voronoi_core": voronoi.as_ref().map(voronoi_json),
        "generation_seconds": start.elapsed().as_secs_f64()
    });
    let report_json = serde_json::to_vec_pretty(&report)
        .map_err(|e| EngineError::new(ErrorCode::SerializationFailed, e))?;
    lap("report");
    let generator = report["algorithm"].as_str().unwrap_or("duct").to_string();
    // each file goes to disk as it is written, hashed on the way; the
    // volume, Gmsh and patch files are written side by side while the
    // hand-off mesh is built
    let fail = crate::handoff::stream_error;
    let mut bundle = dir
        .bundle(name)
        .map_err(|e| EngineError::new(ErrorCode::IoError, e))?;
    let mesh_ref = &mesh;
    let mut files: Vec<(&str, mesh_io::Content<'_>)> = vec![
        (
            "volume.vtu",
            Box::new(move |w: &mut dyn std::io::Write| mesh_io::duct::volume(mesh_ref, w)),
        ),
        (
            "mesh.msh",
            Box::new(move |w: &mut dyn std::io::Write| {
                mesh_io::duct::gmsh(mesh_ref, mesh_io::GmshVersion::V4_1, w)
            }),
        ),
        (
            "mesh.vtm",
            Box::new(|w: &mut dyn std::io::Write| w.write_all(mesh_io::channel::MULTIBLOCK)),
        ),
    ];
    for (i, file) in ["inlet.vtu", "outlet.vtu", "walls.vtu"]
        .into_iter()
        .enumerate()
    {
        files.push((
            file,
            Box::new(move |w: &mut dyn std::io::Write| mesh_io::duct::patch(mesh_ref, i + 1, w)),
        ));
    }
    let report_bytes = &report_json;
    files.push((
        "report.json",
        Box::new(move |w: &mut dyn std::io::Write| w.write_all(report_bytes)),
    ));
    let (handoff, written) = std::thread::scope(|scope| {
        let handoff = scope.spawn(|| mesh_io::duct::handoff(mesh_ref));
        let written = bundle.stream_all(files);
        (handoff.join(), written)
    });
    // the mesh's own error first, then the file system's
    let handoff = handoff
        .map_err(|_| crate::handoff::export_error("the hand-off worker stopped"))?
        .map_err(crate::handoff::export_error)?;
    written.map_err(fail)?;
    lap("hand-off mesh, volume, Gmsh and patch files");
    bundle
        .stream(handoff_format.file_name(), |w| {
            handoff_format.write(&handoff, &[], w)
        })
        .map_err(fail)?;
    lap("hand-off file");
    let evidence = crate::handoff::evidence(&handoff, &generator, bundle.files())?;
    lap("evidence");
    bundle
        .stream(crate::handoff::MESH_EVIDENCE, |w| w.write_all(&evidence))
        .map_err(fail)?;
    bundle.finish();
    lap("small files and finish");
    Ok(report)
}

/// The polyhedral mesh (median dual with wall layers) of a closed surface
/// built in memory: patch ids 1 and 2 are planar caps, 3 and up walls. The
/// region behind the cuts of a sweep fill comes this way (MD-25). Returns
/// the mesh, the scaffold's wall and the build report.
pub(crate) fn mesh_surface(
    surface: &mesh_geometry::ClosedSurface,
    patches: &[usize],
    build: ScaffoldBuild,
    steps: &[f64],
    max_relative_geometry_error: f64,
) -> Result<(duct::DuctMesh, mesh_geometry::ClosedSurface, Option<Value>), EngineError> {
    let validated = |e: String| EngineError::new(ErrorCode::ValidationFailed, e);
    let spec = ScaffoldConfig {
        file: None,
        sha256: None,
        build: Some(build),
        max_relative_geometry_error,
        max_wall_deviation: None,
        core: CoreKind::MedianDual,
        lloyd_iterations: None,
        wall_join: WallJoinKind::default(),
        hex_core: None,
    };
    let (input, built, _, _) =
        load_scaffold(surface, patches, &spec, None, steps.iter().sum(), &[])?;
    let (mesh, shell) = duct::build_volume_dual(input, steps).map_err(validated)?;
    let error = shell.volume() / surface.volume() - 1.;
    if error.abs() > max_relative_geometry_error {
        return Err(EngineError::new(
            ErrorCode::InvalidDomain,
            format!(
                "the remeshed region behind the cuts differs from its surface by {:.3}% in volume (limit {:.3}%)",
                100. * error,
                100. * max_relative_geometry_error
            ),
        ));
    }
    Ok((mesh, shell, built))
}

/// The surface and per-facet patch ids (1 inlet, 2 outlet, 3 wall) of a
/// checked STL, patch map and partition.
fn prepare(
    data: &[u8],
    patch_map: &PatchMap,
    plan: &Partition,
) -> Result<(mesh_geometry::ClosedSurface, Vec<usize>, Patches), EngineError> {
    let invalid = |s: &str| EngineError::new(ErrorCode::InvalidDomain, s);
    require_schema(plan.schema_version)?;
    if sha256(data) != plan.source_sha256 {
        return Err(invalid(
            "patch map or partition digest differs from input STL",
        ));
    }
    let (surface, mut patches, names) = surface_tags(data, patch_map)?;
    if let Some(name) = plan
        .layers_yield
        .iter()
        .find(|n| plan.no_layers.contains(n))
    {
        return Err(invalid(&format!(
            "patch {name:?} is listed under both no_layers and layers_yield"
        )));
    }
    for name in plan.no_layers.iter().chain(&plan.layers_yield) {
        let ids: Vec<usize> = (1..=names.len())
            .filter(|&id| names[id - 1].0 == *name && duct::is_wall(id))
            .collect();
        if ids.is_empty() && !names.iter().any(|(n, _)| n == name) {
            return Err(invalid(&format!(
                "no_layers or layers_yield names an unknown patch {name:?}"
            )));
        }
        for p in patches.iter_mut().filter(|p| ids.contains(p)) {
            *p += duct::SLIP;
        }
    }
    // A flush opening in a wall without layers has none either.
    if !plan.no_layers.is_empty() || !plan.layers_yield.is_empty() {
        let geometry = SurfaceGeometry::new(&surface);
        for id in duct::WALL..=names.len() {
            if names[id - 1].1 == "wall" || !patches.contains(&id) {
                continue;
            }
            let inside = geometry.edges.iter().any(|e| {
                let [f, g] = e.facets;
                (patches[f] == id && duct::is_slip(patches[g]))
                    || (patches[g] == id && duct::is_slip(patches[f]))
            });
            if inside {
                for p in patches.iter_mut().filter(|p| **p == id) {
                    *p += duct::SLIP;
                }
            }
        }
    }
    Ok((surface, patches, names))
}

/// The report entry of every shell built into the scaffold: its two
/// patches (front: the side its normal points into), the measured area of
/// each side, the shell surface and what the crack of the scaffold did.
fn scaffold_shell_report(
    mesh: &duct::DuctMesh,
    stl_patches: usize,
    shells: &[(
        ShellConfig,
        shell_surface::ShellSurfaceStats,
        duct::CrackStats,
        Option<ShellSizes>,
    )],
) -> Vec<Value> {
    let area = |id: usize| mesh.measurements.areas.get(id - 1).copied().unwrap_or(0.0);
    let polygon_area = |p: &[[f64; 3]]| {
        let v: Vec<mesh_core::geometry::Vec3> = p
            .iter()
            .map(|&[x, y, z]| mesh_core::geometry::Vec3::new(x, y, z))
            .collect();
        let mut n = mesh_core::geometry::Vec3::ZERO;
        for i in 0..v.len() {
            n = n + (v[i] - v[0]).cross(v[(i + 1) % v.len()] - v[0]);
        }
        n.length() / 2.0
    };
    shells
        .iter()
        .enumerate()
        .map(|(k, (s, st, crack, sizes))| {
            let (front, back) = (stl_patches + 1 + 2 * k, stl_patches + 2 + 2 * k);
            let faces = |id: usize| mesh.boundary.iter().filter(|f| f.patch == id).count();
            json!({
                "name": s.name,
                "method": "scaffold_baffle_v1",
                "patches": {"front": {"id": front, "name": s.name},
                            "back": {"id": back, "name": format!("{}_back", s.name)}},
                "polygon_area": polygon_area(&s.polygon),
                "surface_area": st.area,
                "front_area": area(front),
                "back_area": area(back),
                "front_faces": faces(front),
                "back_faces": faces(back),
                "shell_triangles": st.triangles,
                "free_edge_length": st.free_edge_length,
                "wall_line_length": st.wall_line_length,
                "free_edge_ends": st.free_edge_ends,
                "min_triangle_angle_degrees": st.min_angle_degrees,
                "sizes": sizes.map(|z| json!({
                    "size": z.size,
                    "size_from": if z.size_set {"set"} else {"default: max_size or narrow width / 8"},
                    "edge_size": z.edge_size,
                    "edge_size_from": if z.edge_size_set {"set"} else {"default: size / 2"},
                })),
                "achieved": {
                    "shell_edge_min": st.edge_lengths[0],
                    "shell_edge_median": st.edge_lengths[1],
                    "shell_edge_max": st.edge_lengths[2],
                    "free_edge_piece_mean": (st.free_edge_pieces > 0)
                        .then(|| st.free_edge_length / st.free_edge_pieces as f64),
                    "wall_line_splits": st.wall_line_splits,
                    "free_edge_splits": st.free_edge_splits,
                },
                "scaffold_node_copies": crack.copies,
                "free_edge_nodes": crack.single_nodes,
            })
        })
        .collect()
}

/// Cut the shells into the finished mesh, in their order; the report entry
/// of each.
fn imprint_shells(
    mesh: &mut duct::DuctMesh,
    shells: &[ShellConfig],
) -> Result<Vec<Value>, EngineError> {
    let mut out = Vec::with_capacity(shells.len());
    for s in shells {
        let shell = duct::PlanarShell {
            name: s.name.clone(),
            polygon: s
                .polygon
                .iter()
                .map(|&[x, y, z]| mesh_core::geometry::Vec3::new(x, y, z))
                .collect(),
        };
        let st = duct::imprint_planar_shell(mesh, &shell)
            .map_err(|e| EngineError::new(ErrorCode::InvalidDomain, e))?;
        out.push(json!({
            "name": s.name,
            "patches": {"front": {"id": st.patches[0], "name": s.name},
                        "back": {"id": st.patches[1], "name": format!("{}_back", s.name)}},
            "polygon_area": st.polygon_area,
            "front_area": st.front_area,
            "back_area": st.back_area,
            "front_faces": st.front_faces,
            "back_faces": st.back_faces,
            "cut_cells": st.cut_cells,
            "through_cuts": st.through_cuts,
            "snapped_vertices": st.snapped_vertices,
            "doubled_vertices": st.doubled_vertices,
        }));
    }
    Ok(out)
}

/// Drop the [`duct::SLIP`] mark of the patches without layers from the
/// boundary faces and the measured areas of a built mesh with `n` patches.
fn unmark_slip(mesh: &mut duct::DuctMesh, n: usize) {
    for face in &mut mesh.boundary {
        face.patch = duct::patch_index(face.patch);
    }
    let areas = &mut mesh.measurements.areas;
    for i in duct::SLIP..areas.len() {
        areas[i - duct::SLIP] += areas[i];
    }
    areas.truncate(n.max(duct::WALL));
}

/// Name and role of every patch id, indexed by id - 1.
type Patches = Vec<(String, String)>;

/// The surface, per-facet patch ids and patch names of a checked STL and
/// patch map: 1 the inlet cap, 2 the outlet cap, 3 and up the patches that
/// grow wall layers ([`duct::is_wall`]). An inlet or outlet whose rim is
/// flush with the wall around it (every rim edge turns by less than the
/// feature angle) cannot contract in its plane; it grows layers like the
/// wall and keeps its name and role.
fn surface_tags(
    data: &[u8],
    patch_map: &PatchMap,
) -> Result<(mesh_geometry::ClosedSurface, Vec<usize>, Patches), EngineError> {
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
    let (mut patches, mut names) = match patch_map {
        PatchMap::Ports(p) => (port_tags(p, facet_count)?, duct::classic_patches()),
        PatchMap::Map(m) => map_tags(m)?,
    };
    let geometry = SurfaceGeometry::new(&surface);
    let limit = default_feature_angle().to_radians();
    for cap in 1..duct::WALL {
        let mut rim = 0;
        let mut flush = true;
        for e in &geometry.edges {
            let [f, g] = e.facets;
            if (patches[f] == cap) != (patches[g] == cap) {
                rim += 1;
                flush &= e.angle < limit;
            }
        }
        if rim > 0 && flush {
            names.push(names[cap - 1].clone());
            let id = names.len();
            for p in patches.iter_mut().filter(|p| **p == cap) {
                *p = id;
            }
        }
    }
    if !patches.iter().any(|&p| duct::is_wall(p)) {
        return Err(invalid("expected a wall patch"));
    }
    Ok((surface, patches, names))
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
    let (surface, patches, names) = surface_tags(data, patch_map)?;
    let scale = surface.bounds().extent().max_abs_component();
    // the openings: the inlet and outlet caps, or flush openings
    let mut areas = vec![0.; names.len()];
    for (f, &tag) in surface.facets().iter().zip(&patches) {
        areas[tag - 1] += f.normalized_area * scale * scale;
    }
    let smallest = names
        .iter()
        .zip(&areas)
        .filter(|((_, role), a)| role != "wall" && **a > 0.)
        .map(|(_, &a)| a)
        .fold(f64::INFINITY, f64::min);
    if !smallest.is_finite() {
        return Err(invalid(
            "setup sizes the mesh from the smallest inlet or outlet; the map has none",
        ));
    }
    let r = (smallest / std::f64::consts::PI).sqrt();
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
        shells: Vec::new(),
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
                grading: None,
            }),
            max_relative_geometry_error: 0.002,
            max_wall_deviation: Some(round(1.5 * tolerance)),
            core: CoreKind::default(),
            lloyd_iterations: None,
            wall_join: WallJoinKind::default(),
            hex_core: None,
        }),
        first_thickness: round(first_thickness),
        growth,
        layers,
        no_layers: Vec::new(),
        layers_yield: Vec::new(),
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
        if duct::is_cap(patches[f])
            || duct::is_cap(patches[k])
            || e.angle >= limit
            || e.angle < 1e-3
        {
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

/// Patch ids and names of a patch map: the patch with role inlet is 1 and
/// the one with role outlet 2 (at most one of each), and the wall patches
/// get 3, 4, ... in map order. Ids without a patch keep a placeholder name.
fn map_tags(map: &crate::patches::PatchMap) -> Result<(Vec<usize>, Patches), EngineError> {
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
    let mut names = duct::classic_patches();
    names.truncate(duct::WALL - 1);
    let mut tag = Vec::with_capacity(map.patches.len());
    let (mut inlets, mut outlets) = (0, 0);
    for p in &map.patches {
        tag.push(match p.role {
            PatchRole::Inlet => {
                inlets += 1;
                names[0] = (p.name.clone(), "inlet".into());
                1
            }
            PatchRole::Outlet => {
                outlets += 1;
                names[1] = (p.name.clone(), "outlet".into());
                2
            }
            PatchRole::Wall => {
                names.push((p.name.clone(), "wall".into()));
                names.len()
            }
            other => {
                return Err(bad(format!(
                    "a duct takes inlet, outlet and wall patches; {} has role {other:?}",
                    p.name
                )));
            }
        });
    }
    if inlets > 1 || outlets > 1 {
        return Err(bad(format!(
            "a duct takes at most one inlet and one outlet patch, the map has {inlets} and {outlets}"
        )));
    }
    let ids = map
        .facet_patch
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
        .collect::<Result<_, _>>()?;
    Ok((ids, names))
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
    let (surface, patches, _) = prepare(data, patch_map, plan)?;
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
        let height = if plan.growth == 1.0 {
            plan.first_thickness * layers as f64
        } else {
            plan.first_thickness * (plan.growth.powi(layers as i32) - 1.) / (plan.growth - 1.)
        };
        let (scaffold, _, hex, _) =
            load_scaffold(&surface, &patches, config, scaffold_data, height, &[])?;
        let walls = scaffold
            .boundary
            .iter()
            .filter(|(_, p)| duct::is_wall(*p))
            .count();
        // median dual: one polygon prism per sector of a wall vertex and
        // layer, one cell per scaffold node (and per cube of a hex core);
        // Voronoi core: one wedge (or
        // split-cap polyhedron) per wall triangle and layer, one cell per
        // site (the scaffold nodes the tetrahedra use). Wall columns follow
        // the built core: two per wall triangle and layer bound the
        // fixtures from above (elbow 0.8, Manifold 2.1).
        let (per_wall, core) = match config.core {
            CoreKind::MedianDual | CoreKind::HexCore => {
                let (triangles, patches): (Vec<[usize; 3]>, Vec<Option<usize>>) = scaffold
                    .boundary
                    .iter()
                    .filter(|(_, p)| duct::is_wall(*p))
                    .map(|&(t, p)| (t, Some(p)))
                    .unzip();
                let sectors = duct::prism_sectors(&scaffold.points, &triangles, &patches);
                let prisms = (0..sectors.len()).filter(|&i| sectors[i] == i).count();
                let cubes = hex.as_ref().map_or(0, |c| c.cubes.len());
                return prisms
                    .checked_mul(layers)
                    .and_then(|v| v.checked_add(scaffold.points.len() + cubes))
                    .ok_or_else(|| invalid("cell count overflow"));
            }
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
    let walls = patches.iter().filter(|&&p| duct::is_wall(p)).count();
    surface
        .facets()
        .len()
        .checked_add(2 * seam)
        .and_then(|v| v.checked_add(walls.checked_mul(layers)?))
        .ok_or_else(|| invalid("cell count overflow"))
}

/// The scaffold of a polyhedral run: the configured file (digest checked)
/// or, with `build`, the STL wall remeshed onto itself to the chord
/// tolerance and tetrahedralized in Rust. Returns the scaffold, the
/// build's statistics, the hex core (when asked) and the surface of each
/// shell built into it (`shells` only with a built scaffold).
#[allow(clippy::type_complexity)]
fn load_scaffold(
    surface: &mesh_geometry::ClosedSurface,
    patches: &[usize],
    spec: &ScaffoldConfig,
    data: Option<&[u8]>,
    layer_height: f64,
    shells: &[ShellConfig],
) -> Result<
    (
        duct::Scaffold,
        Option<Value>,
        Option<HexCore>,
        Vec<shell_surface::ShellSurfaceStats>,
    ),
    EngineError,
> {
    let invalid = |s: &str| EngineError::new(ErrorCode::InvalidDomain, s);
    let validated = |e: String| EngineError::new(ErrorCode::ValidationFailed, e);
    match (&spec.file, &spec.build) {
        (Some(_), None) => {
            let bytes = data.ok_or_else(|| invalid("missing configured volume scaffold"))?;
            if spec.sha256.as_deref() != Some(sha256(bytes).as_str()) {
                return Err(invalid("volume scaffold digest mismatch"));
            }
            if spec.core == CoreKind::HexCore {
                return Err(invalid("the hex core needs a built scaffold"));
            }
            Ok((
                scaffold::read(bytes).map_err(validated)?,
                None,
                None,
                Vec::new(),
            ))
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
            let field = build_field(b, surface, shells)?;
            let start = Instant::now();
            let triangles: Vec<[usize; 3]> = surface
                .facets()
                .iter()
                .map(|f| f.triangle.vertices)
                .collect();
            let keep: Vec<bool> = patches.iter().map(|&p| duct::is_cap(p)).collect();
            // the line where each shell meets the wall, cut into the STL so
            // the remeshed wall has edges along it (SHL-01)
            let polygons: Vec<Vec<mesh_core::geometry::Vec3>> = shells
                .iter()
                .map(|s| s.polygon.iter().map(|&p| vector(p)).collect())
                .collect();
            let named: Vec<(&str, &[mesh_core::geometry::Vec3])> = shells
                .iter()
                .zip(&polygons)
                .map(|(s, p)| (s.name.as_str(), p.as_slice()))
                .collect();
            let imprinted = if shells.is_empty() {
                None
            } else {
                Some(
                    shell_surface::imprint(surface.vertices(), &triangles, &keep, patches, &named)
                        .map_err(|e| invalid(&e))?,
                )
            };
            let (source_v, source_t, source_keep, source_group, lines) = match &imprinted {
                Some(i) => (
                    i.vertices.as_slice(),
                    i.triangles.as_slice(),
                    i.keep.as_slice(),
                    i.group.as_slice(),
                    i.lines.as_slice(),
                ),
                None => (
                    surface.vertices(),
                    triangles.as_slice(),
                    keep.as_slice(),
                    patches,
                    &[][..],
                ),
            };
            let wall = remesh::remesh_lines(
                source_v,
                source_t,
                source_keep,
                source_group,
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
                field.as_ref(),
                lines,
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
            let mut tags: Vec<usize> = wall.group.clone();
            let mut vertices = wall.vertices.clone();
            let mut triangles = wall.triangles.clone();
            let hex_start = Instant::now();
            let core = if spec.core == CoreKind::HexCore {
                let c = spec.hex_core.unwrap_or(HexCoreConfig {
                    size_factor: None,
                    margin: None,
                });
                let defaults = hexcore::HexCoreOptions::default();
                let core = hexcore::build(
                    &wall.vertices,
                    &wall.triangles,
                    b.max_size,
                    b.rate,
                    field.as_ref(),
                    layer_height,
                    hexcore::HexCoreOptions {
                        size_factor: c.size_factor.unwrap_or(defaults.size_factor),
                        margin: c.margin.unwrap_or(defaults.margin),
                        ..defaults
                    },
                )
                .map_err(validated)?;
                let offset = vertices.len();
                vertices.extend_from_slice(&core.shell_points);
                for t in &core.shell_triangles {
                    triangles.push(t.map(|v| v + offset));
                    tags.push(duct::HEX);
                }
                Some(core)
            } else {
                None
            };
            let hex_seconds = hex_start.elapsed().as_secs_f64();
            let hole = core
                .as_ref()
                .map(|c| move |x: mesh_core::geometry::Vec3| c.contains(x));
            // each shell meshed in its plane at the scaffold's spacing
            let line_edges: Vec<(usize, usize)> = wall.lines.iter().map(|&(e, _)| e).collect();
            let spacing = shell_surface::Spacing::new(
                &wall.vertices,
                &wall.triangles,
                &line_edges,
                b.rate,
                b.max_size,
                field.as_ref(),
            );
            let mut shell_triangles = Vec::new();
            let mut shell_stats = Vec::new();
            for (k, polygon) in polygons.iter().enumerate() {
                let line: Vec<(usize, usize)> = wall
                    .lines
                    .iter()
                    .filter(|&&(_, l)| l == k)
                    .map(|&(e, _)| e)
                    .collect();
                let m = shell_surface::mesh_shell(
                    polygon,
                    &wall.vertices,
                    &wall.triangles,
                    &line,
                    &|x| spacing.at(x),
                    vertices.len(),
                )
                .map_err(|e| invalid(&format!("shell {:?}: {e}", shells[k].name)))?;
                vertices.extend(m.vertices);
                split_wall_edges(&mut triangles, &mut tags, &m.wall_splits)
                    .map_err(|e| invalid(&format!("shell {:?}: {e}", shells[k].name)))?;
                shell_triangles.extend(m.triangles.into_iter().map(|t| (t, k)));
                shell_stats.push(m.stats);
            }
            let (scaffold, st) = scaffold_build::build_with_shells(
                &vertices,
                &triangles,
                &tags,
                &shell_triangles,
                scaffold_build::ScaffoldOptions {
                    max_size: b.max_size,
                    rate: b.rate,
                    max_points: MAX_SCAFFOLD_POINTS,
                    ..Default::default()
                },
                field.as_ref(),
                hole.as_ref()
                    .map(|h| h as &dyn Fn(mesh_core::geometry::Vec3) -> bool),
            )
            .map_err(validated)?;
            let w = wall.stats;
            let hex_report = core.as_ref().map(|c| {
                json!({
                    "cubes": c.stats.cubes,
                    "cube_sizes": c.stats.sizes.iter().map(|(e, n)| json!({"edge": e, "count": n})).collect::<Vec<_>>(),
                    "removed_for_manifold_shell": c.stats.removed_for_manifold,
                    "shell_triangles": c.stats.shell_triangles,
                    "octree_seconds": hex_seconds,
                })
            });
            Ok((
                scaffold,
                Some(json!({
                    "hex_core": hex_report,
                    "settings": b,
                    "wall_triangles": tags.iter().filter(|&&t| duct::is_wall(t)).count(),
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
                core,
                shell_stats,
            ))
        }
        _ => Err(invalid(
            "volume scaffold needs exactly one of file or build",
        )),
    }
}

/// Split the wall triangles at the wall-line edges a shell split: each
/// triangle on such an edge becomes a fan from its third corner through
/// the new points (in order from the edge's first vertex). The points lie
/// on the edges, so the wall keeps its shape.
fn split_wall_edges(
    triangles: &mut Vec<[usize; 3]>,
    tags: &mut Vec<usize>,
    splits: &[((usize, usize), Vec<usize>)],
) -> Result<(), String> {
    if splits.is_empty() {
        return Ok(());
    }
    let key = |a: usize, b: usize| (a.min(b), a.max(b));
    let mut on: std::collections::HashMap<(usize, usize), Vec<usize>> =
        std::collections::HashMap::new();
    for ((a, b), _) in splits {
        on.insert(key(*a, *b), Vec::new());
    }
    for (f, t) in triangles.iter().enumerate() {
        for i in 0..3 {
            if let Some(list) = on.get_mut(&key(t[i], t[(i + 1) % 3])) {
                list.push(f);
            }
        }
    }
    for ((a, b), points) in splits {
        let faces = on[&key(*a, *b)].clone();
        if faces.len() != 2 {
            return Err(format!(
                "a split wall-line edge has {} wall triangles, not two",
                faces.len()
            ));
        }
        for f in faces {
            let t = triangles[f];
            let Some(i) = (0..3).find(|&i| key(t[i], t[(i + 1) % 3]) == key(*a, *b)) else {
                return Err("a split wall-line edge was split twice".into());
            };
            let (u, w, apex) = (t[i], t[(i + 1) % 3], t[(i + 2) % 3]);
            let mut chain = vec![u];
            if u == *a {
                chain.extend(points.iter().copied());
            } else {
                chain.extend(points.iter().rev().copied());
            }
            chain.push(w);
            triangles[f] = [chain[0], chain[1], apex];
            let mut pieces = vec![f];
            for pair in chain[1..].windows(2) {
                pieces.push(triangles.len());
                triangles.push([pair[0], pair[1], apex]);
                tags.push(tags[f]);
            }
            // another split edge of the same triangle (at a corner of the
            // line) now lies on one of the pieces
            for &g in &pieces[1..] {
                let tg = triangles[g];
                for i in 0..3 {
                    if let Some(list) = on.get_mut(&key(tg[i], tg[(i + 1) % 3])) {
                        for h in list.iter_mut() {
                            if *h == f {
                                *h = g;
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

/// The size field of a built scaffold: the grading's datum objects and
/// regions and each shell's own sizes (its polygon and free edges, see
/// [`ShellConfig::sizes`]); `None` when neither asks for anything.
fn build_field(
    b: &ScaffoldBuild,
    surface: &mesh_geometry::ClosedSurface,
    shells: &[ShellConfig],
) -> Result<Option<mesh_algorithms::sizing::SizeField>, EngineError> {
    use mesh_algorithms::sizing::{
        Frame, Shape, SizeField, SizeSource, SourceKind, surface::WeightedSet,
    };
    let grading = b.grading.as_ref().map(|g| g.field(b)).transpose()?;
    if shells.is_empty() {
        return Ok(grading);
    }
    let invalid = |m: String| EngineError::new(ErrorCode::InvalidArguments, m);
    let (lo, hi) = (surface.bounds().min(), surface.bounds().max());
    let triangles: Vec<[usize; 3]> = surface
        .facets()
        .iter()
        .map(|f| f.triangle.vertices)
        .collect();
    let mut min_size = grading.as_ref().map_or(b.min_size, |g| g.min_size());
    let mut sources: Vec<SizeSource> = grading
        .as_ref()
        .map(|g| g.sources().to_vec())
        .unwrap_or_default();
    for s in shells {
        let sizes = s.sizes(b, lo, hi)?;
        if sizes.size >= b.max_size && sizes.edge_size >= b.max_size {
            continue;
        }
        let polygon: Vec<mesh_core::geometry::Vec3> =
            s.polygon.iter().map(|&p| vector(p)).collect();
        let items = shell_surface::size_items(
            &polygon,
            sizes.size,
            sizes.edge_size,
            surface.vertices(),
            &triangles,
        )
        .map_err(|e| invalid(format!("shell {}: {e}", s.name)))?;
        let set = WeightedSet::new(items);
        let size = set.smallest().unwrap_or(b.max_size);
        min_size = min_size.min(size);
        sources.push(SizeSource {
            name: format!("shell:{}", s.name),
            kind: SourceKind::Shell,
            shape: Shape::Weighted(std::sync::Arc::new(set)),
            frame: Frame::GLOBAL,
            size,
            distance: 0.0,
            growth: None,
        });
    }
    if sources.is_empty() {
        return Ok(None);
    }
    SizeField::new(sources, min_size, b.max_size, 1.0 + b.rate)
        .map(Some)
        .map_err(|e| EngineError::new(ErrorCode::InvalidSizeField, format!("{e:?}")))
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
        .filter(|(_, p)| !duct::is_cap(**p))
        .map(|(f, _)| f.triangle.vertices.map(|v| surface.vertices()[v]))
        .collect();
    let wall: Vec<[mesh_core::geometry::Vec3; 3]> = scaffold
        .boundary
        .iter()
        .filter(|(_, p)| !duct::is_cap(*p) && !duct::is_hex(*p))
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

/// GRD-08 for the duct path: per grading source, the core cells (no layer cell)
/// whose centroid lies in its zone (within its distance; inside a region)
/// where it sets the size, and their achieved cbrt(volume); with the median
/// of the other core cells for contrast.
fn grading_evidence(field: &mesh_algorithms::sizing::SizeField, mesh: &duct::DuctMesh) -> Value {
    let Ok(h) = mesh_io::duct::handoff(mesh) else {
        return json!({"error": "the mesh has no hand-off form; sizes not measured"});
    };
    let cells = h.cell_volumes();
    let core: Vec<(f64, mesh_core::geometry::Vec3)> = if cells.len() == mesh.cells.len() {
        cells
            .into_iter()
            .zip(&mesh.cells)
            .filter(|(_, c)| c.layer == 0)
            .map(|(m, _)| m)
            .collect()
    } else {
        cells
    };
    let quantiles = |mut v: Vec<f64>| {
        if v.is_empty() {
            return Value::Null;
        }
        v.sort_by(f64::total_cmp);
        let q = |p: f64| v[((v.len() - 1) as f64 * p).round() as usize];
        json!({"min": v[0], "p10": q(0.1), "median": q(0.5), "p90": q(0.9), "max": v[v.len() - 1]})
    };
    let mut claimed = vec![false; core.len()];
    let sources: Vec<Value> = field
        .sources()
        .iter()
        .enumerate()
        .map(|(i, s)| {
            let mut sizes = Vec::new();
            for (k, &(v, c)) in core.iter().enumerate() {
                let near = field.source_distance(i, c) <= s.distance;
                let governs = field.source_size(i, c) <= field.size(c) * (1.0 + 1e-9);
                if near && governs {
                    sizes.push(v.max(0.0).cbrt());
                    claimed[k] = true;
                }
            }
            json!({
                "name": s.name, "kind": s.kind.name(), "target_size": s.size,
                "cells": sizes.len(), "achieved": quantiles(sizes),
            })
        })
        .collect();
    let others = core
        .iter()
        .zip(&claimed)
        .filter(|(_, c)| !**c)
        .map(|(&(v, _), _)| v.max(0.0).cbrt())
        .collect();
    json!({
        "method": "datum_region_size_field_on_scaffold_v1",
        "size_measure": "cbrt(cell volume), core cells only",
        "sources": sources,
        "other_core_cells": quantiles(others),
    })
}

impl Partition {
    fn build_mut(&mut self) -> Result<&mut ScaffoldBuild, EngineError> {
        self.volume_scaffold
            .as_mut()
            .and_then(|s| s.build.as_mut())
            .ok_or_else(|| {
                EngineError::new(
                    ErrorCode::InvalidArguments,
                    "grading needs a partition with a scaffold built in Rust (`volume_scaffold.build`, as mesh.setup_duct writes it)",
                )
            })
    }

    /// The grading sources of the build, by name (GRD-08).
    pub fn grading(&self) -> Option<&Grading> {
        self.volume_scaffold
            .as_ref()
            .and_then(|s| s.build.as_ref())
            .and_then(|b| b.grading.as_ref())
    }

    /// Add a datum or region source (and a frame it names, when new); the
    /// name must be free and the field must build.
    pub fn add_grading(
        &mut self,
        source: crate::sizing::SourceConfig,
        frame: Option<crate::sizing::FrameConfig>,
    ) -> Result<(), EngineError> {
        let invalid = |m: String| EngineError::new(ErrorCode::InvalidSizeField, m);
        let build = self.build_mut()?;
        let mut grading = build.grading.clone().unwrap_or_default();
        if grading.sources.iter().any(|s| s.name() == source.name()) {
            return Err(invalid(format!(
                "a grading source {} exists; remove it first or pick another name",
                source.name()
            )));
        }
        if let Some(f) = frame {
            match grading.frames.iter().find(|g| g.name == f.name) {
                Some(g) if *g != f => {
                    return Err(invalid(format!("frame {} exists with other axes", f.name)));
                }
                Some(_) => {}
                None => grading.frames.push(f),
            }
        }
        grading.sources.push(source);
        grading.field(build)?;
        build.grading = Some(grading);
        Ok(())
    }

    /// Add the zero-thickness shell `shell` (SHL-01): its name and
    /// `<name>_back` must be free, the polygon planar, convex and
    /// counter-clockwise about the `name` side, the sizes positive, and the
    /// shell may not touch another (not supported yet).
    pub fn add_shell(&mut self, shell: ShellConfig) -> Result<(), EngineError> {
        let invalid = |m: String| EngineError::new(ErrorCode::InvalidArguments, m);
        let name = shell.name.as_str();
        if name.is_empty() || name.chars().any(|c| c.is_control() || c.is_whitespace()) {
            return Err(invalid(
                "a shell name must be nonempty, without spaces or control characters".into(),
            ));
        }
        let sides = |n: &str| [n.to_string(), format!("{n}_back")];
        for other in &self.shells {
            if sides(&other.name).iter().any(|o| sides(name).contains(o)) {
                return Err(invalid(format!(
                    "shell {} exists (its sides are {} and {}_back); remove it first or pick another name",
                    other.name, other.name, other.name
                )));
            }
        }
        let polygon: Vec<mesh_core::geometry::Vec3> =
            shell.polygon.iter().map(|&p| vector(p)).collect();
        shell_surface::check_polygon(&polygon)
            .map_err(|e| invalid(format!("shell {name}: {e}")))?;
        for (v, what) in [(shell.size, "size"), (shell.edge_size, "edge_size")] {
            if v.is_some_and(|v| !(v.is_finite() && v > 0.0)) {
                return Err(invalid(format!("shell {name}: {what} must be positive")));
            }
        }
        for other in &self.shells {
            let theirs: Vec<mesh_core::geometry::Vec3> =
                other.polygon.iter().map(|&p| vector(p)).collect();
            if shell_surface::polygons_touch(&polygon, &theirs).unwrap_or(false) {
                return Err(invalid(format!(
                    "shell {name} touches shell {}; shells that meet are not supported yet",
                    other.name
                )));
            }
        }
        self.shells.push(shell);
        Ok(())
    }

    /// Remove the shell `name`.
    pub fn remove_shell(&mut self, name: &str) -> Result<(), EngineError> {
        let before = self.shells.len();
        self.shells.retain(|s| s.name != name);
        if self.shells.len() == before {
            let known: Vec<&str> = self.shells.iter().map(|s| s.name.as_str()).collect();
            return Err(EngineError::new(
                ErrorCode::InvalidArguments,
                if known.is_empty() {
                    format!("no shell {name}; the partition has no shells")
                } else {
                    format!("no shell {name}; the partition has {}", known.join(", "))
                },
            ));
        }
        Ok(())
    }

    /// Each shell with the sizes a build would use (the defaults filled in
    /// when the partition has a Rust-built scaffold; `lo`, `hi` bound the
    /// STL).
    pub fn shell_sizes(
        &self,
        lo: mesh_core::geometry::Vec3,
        hi: mesh_core::geometry::Vec3,
    ) -> Result<Vec<Option<ShellSizes>>, EngineError> {
        let build = self.volume_scaffold.as_ref().and_then(|s| s.build.as_ref());
        self.shells
            .iter()
            .map(|s| build.map(|b| s.sizes(b, lo, hi)).transpose())
            .collect()
    }

    /// The shells as `shell.list` shows them: name, sides, polygon and its
    /// area, and the sizes; with the partition's STL (`stl`) the default
    /// sizes are filled in as a build would choose them.
    pub fn shell_listing(&self, stl: Option<&[u8]>) -> Result<Value, EngineError> {
        let sizes = match stl {
            Some(data) => {
                if sha256(data) != self.source_sha256 {
                    return Err(EngineError::new(
                        ErrorCode::InvalidDomain,
                        "the STL differs from the partition's source",
                    ));
                }
                let surface = import_closed_stl(data, StlOptions::default())
                    .map_err(|e| EngineError::new(ErrorCode::InvalidDomain, e))?;
                self.shell_sizes(surface.bounds().min(), surface.bounds().max())?
            }
            None => vec![None; self.shells.len()],
        };
        let shells: Vec<Value> = self
            .shells
            .iter()
            .zip(sizes)
            .map(|(s, z)| {
                let polygon: Vec<mesh_core::geometry::Vec3> =
                    s.polygon.iter().map(|&p| vector(p)).collect();
                let mut n = mesh_core::geometry::Vec3::ZERO;
                for i in 0..polygon.len() {
                    n = n
                        + (polygon[i] - polygon[0])
                            .cross(polygon[(i + 1) % polygon.len()] - polygon[0]);
                }
                let area = n.length() / 2.0;
                let normal = if area > 0.0 { n / (2.0 * area) } else { n };
                let (size, size_from, edge, edge_from) = match z {
                    Some(z) => (
                        json!(z.size),
                        if z.size_set {
                            "set"
                        } else {
                            "default: max_size or narrow width / 8"
                        },
                        json!(z.edge_size),
                        if z.edge_size_set {
                            "set"
                        } else {
                            "default: size / 2"
                        },
                    ),
                    None => (
                        json!(s.size),
                        if s.size.is_some() {
                            "set"
                        } else {
                            "default (give `surface` to see it)"
                        },
                        json!(s.edge_size),
                        if s.edge_size.is_some() {
                            "set"
                        } else {
                            "default (give `surface` to see it)"
                        },
                    ),
                };
                json!({
                    "name": s.name,
                    "sides": {"front": s.name, "back": format!("{}_back", s.name)},
                    "polygon": s.polygon,
                    "normal": [normal.x, normal.y, normal.z],
                    "polygon_area": area,
                    "size": size, "size_from": size_from,
                    "edge_size": edge, "edge_size_from": edge_from,
                })
            })
            .collect();
        Ok(json!({ "shells": shells }))
    }

    /// Remove the source `name`; frames no source uses any more go too.
    pub fn remove_grading(&mut self, name: &str) -> Result<(), EngineError> {
        let build = self.build_mut()?;
        let Some(mut grading) = build.grading.clone() else {
            return Err(EngineError::new(
                ErrorCode::InvalidArguments,
                "the partition has no grading",
            ));
        };
        let before = grading.sources.len();
        grading.sources.retain(|s| s.name() != name);
        if grading.sources.len() == before {
            let known: Vec<&str> = grading.sources.iter().map(|s| s.name()).collect();
            return Err(EngineError::new(
                ErrorCode::InvalidArguments,
                format!(
                    "no grading source {name}; the partition has {}",
                    known.join(", ")
                ),
            ));
        }
        let used: Vec<String> = grading
            .sources
            .iter()
            .filter_map(|s| s.frame().map(str::to_string))
            .collect();
        grading.frames.retain(|f| used.contains(&f.name));
        build.grading = (!grading.sources.is_empty()).then_some(grading);
        Ok(())
    }
}
