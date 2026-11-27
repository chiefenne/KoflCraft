//! The tool registry of KoflCraft (AG-12): every engine operation a front
//! end can call, with a stable dotted id, a contract version, input and
//! output JSON Schemas and an action kind. The command line (`kofl`,
//! `mesh-cli`) and the MCP server dispatch these ids and nothing else, and
//! every call ends in the same result envelope (`koflcraft.tool-result`).
//!
//! Conventions follow FoehnCraft (FC-07): ids `area.verb`, `system.*`,
//! `diag.codes`, `format`/`format_version` on every machine format,
//! diagnostic codes as `AREA000`.

mod block;
pub mod context;
pub mod diag;
mod estimate;
mod grading;
pub mod journal;
mod markdown;
mod mesh;
mod page;
mod patch;
mod regions;
mod schemas;
mod shell;

use std::time::Instant;

use serde::Serialize;
use serde_json::{Value, json};

pub use context::{Context, Limits};
pub use diag::{CODES, Diagnostic};
pub use estimate::{ESTIMATE_FORMAT, EstimateInput};
pub use journal::{Journal, ReplayInput, call_journaled};
pub use markdown::markdown;
pub use mesh::{
    ArtifactInput, ChannelInput, ConfigInput, DuctInput, ExportInput, GenerateInput,
    GenerateStlInput, ReportInput, StlInput, SweepInput,
};

pub const RESULT_FORMAT: &str = "koflcraft.tool-result";
pub const CATALOGUE_FORMAT: &str = "koflcraft.tool-catalogue";
pub const FORMAT_VERSION: &str = "1.0";

/// What a tool does to the world (AG-12, AG-15).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Action {
    /// Reads files and computes.
    Read,
    /// Writes new files.
    WriteFiles,
    /// Starts a mesh build (and writes its files).
    StartBuild,
}

type RunFn = fn(Value, &Context) -> mesh::Run;

/// One registered tool.
pub struct Tool {
    pub id: &'static str,
    pub version: u32,
    pub action: Action,
    pub summary: &'static str,
    /// The name the tool had in `mesh-mcp` before the registry.
    pub legacy_name: Option<&'static str>,
    /// Input fields that name files the call writes (AG-15 effects).
    pub outputs: &'static [&'static str],
    input_schema: fn() -> Value,
    output_schema: fn() -> Value,
    run: RunFn,
}

impl Tool {
    pub fn input_schema(&self) -> Value {
        (self.input_schema)()
    }

    pub fn output_schema(&self) -> Value {
        (self.output_schema)()
    }

    /// The files a call with `input` would create, before it runs.
    pub fn effects(&self, input: &Value) -> Vec<String> {
        self.outputs
            .iter()
            .filter_map(|f| input[*f].as_str().map(str::to_string))
            .collect()
    }
}

fn no_input() -> Value {
    schemas::input(json!({}), &[])
}

macro_rules! out {
    ($name:literal) => {{
        fn f() -> Value {
            schemas::output_schema($name)
        }
        f
    }};
}

/// Every tool, in catalogue order.
pub static TOOLS: &[Tool] = &[
    Tool {
        id: "system.info",
        version: 1,
        action: Action::Read,
        summary: "Version, commit, platform and features of this build.",
        legacy_name: None,
        outputs: &[],
        input_schema: no_input,
        output_schema: || {
            schemas::object_schema(json!({
                "format": {"const": "koflcraft.info"}, "format_version": {"type": "string"},
                "version": {"type": "string"}, "commit": {"type": "string"},
                "platform": {"type": "string"}, "debug_assertions": {"type": "boolean"}
            }))
        },
        run: |_, _| Ok(system_info()),
    },
    Tool {
        id: "system.commands",
        version: 1,
        action: Action::Read,
        summary: "The catalogue of every tool with its input and output schemas.",
        legacy_name: None,
        outputs: &[],
        input_schema: no_input,
        output_schema: || json!({"type": "object", "description": "`koflcraft.tool-catalogue`"}),
        run: |_, _| Ok(catalogue()),
    },
    Tool {
        id: "diag.codes",
        version: 1,
        action: Action::Read,
        summary: "The diagnostic codes with their severity, meaning and fix (AG-18).",
        legacy_name: None,
        outputs: &[],
        input_schema: no_input,
        output_schema: || json!({"type": "object", "description": "`koflcraft.codes`"}),
        run: |_, _| Ok(codes()),
    },
    Tool {
        id: "domain.inspect",
        version: 1,
        action: Action::Read,
        summary: "Inspect box bounds, sites and the limitations of the reference backend.",
        legacy_name: Some("inspect_domain"),
        outputs: &[],
        input_schema: || schemas::input(json!({"config": schemas::mesh_config()}), &[]),
        output_schema: out!("inspect_domain"),
        run: mesh::inspect,
    },
    Tool {
        id: "domain.inspect_stl",
        version: 1,
        action: Action::Read,
        summary: "Validate a convex STL shell and report bounds, volume and named patches (from the STL solids, or from a patch map given as `patches`).",
        legacy_name: Some("inspect_stl_domain"),
        outputs: &[],
        input_schema: || {
            schemas::input(
                json!({"surface": schemas::path(), "config": schemas::stl_config(), "patches": schemas::path()}),
                &["patches"],
            )
        },
        output_schema: out!("inspect_stl_domain"),
        run: mesh::inspect_stl,
    },
    Tool {
        id: "patch.init",
        version: 1,
        action: Action::WriteFiles,
        summary: "Start a patch map for a closed STL (convex or not) from its solid names: one patch per solid, names made FoehnCraft boundary names, roles from the names (BC-01).",
        legacy_name: None,
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"surface": schemas::path(), "output": schemas::path()}),
                &[],
            )
        },
        output_schema: schemas::patch_result,
        run: patch::init,
    },
    Tool {
        id: "patch.list",
        version: 1,
        action: Action::Read,
        summary: "Check a patch map against its STL (every facet in exactly one patch, BC-06) and report each patch: role, facets, area, mean normal and its spread, bounds, regions and neighbouring patches.",
        legacy_name: None,
        outputs: &[],
        input_schema: || {
            schemas::input(
                json!({"surface": schemas::path(), "patches": schemas::path()}),
                &[],
            )
        },
        output_schema: schemas::patch_summary,
        run: patch::list,
    },
    Tool {
        id: "patch.select_plane",
        version: 1,
        action: Action::WriteFiles,
        summary: "Move the facets whose three corners lie on a plane (within `tolerance`) into the patch `name`, optionally only from the patches in `within`; write the new map to `output`.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::patch_plane,
        output_schema: schemas::patch_result,
        run: patch::select_plane,
    },
    Tool {
        id: "patch.select_box",
        version: 1,
        action: Action::WriteFiles,
        summary: "Move the facets inside an axis-aligned box (by centroid, or `contain: all_corners`) into the patch `name`.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::patch_box,
        output_schema: schemas::patch_result,
        run: patch::select_box,
    },
    Tool {
        id: "patch.select_normal",
        version: 1,
        action: Action::WriteFiles,
        summary: "Move the facets whose outward normal is within `max_angle_deg` of `direction` into the patch `name`.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::patch_normal,
        output_schema: schemas::patch_result,
        run: patch::select_normal,
    },
    Tool {
        id: "patch.grow",
        version: 1,
        action: Action::WriteFiles,
        summary: "Feature-angle region growing: from `seed_facet` (or the facet nearest to `seed_point`) across every edge whose dihedral angle is at most `feature_angle_deg`, and move the region into the patch `name`.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::patch_grow,
        output_schema: schemas::patch_result,
        run: patch::grow,
    },
    Tool {
        id: "patch.detect_ports",
        version: 1,
        action: Action::WriteFiles,
        summary: "Find the planar end caps of the pipes of a bare closed STL (convex or not) and write a map with `walls`, `inlet` and `outlet` (BC-03, ports): planar disc regions ranked by circularity and a pipe-end rim. Lists the candidates; `inlet` and `outlet` pick others by id.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"surface": schemas::path(), "output": schemas::path(),
                       "inlet": {"type": "integer", "minimum": 0}, "outlet": {"type": "integer", "minimum": 0},
                       "planar_angle_deg": {"type": "number", "minimum": 0, "maximum": 30},
                       "candidates": {"type": "integer", "minimum": 2}}),
                &["inlet", "outlet", "planar_angle_deg", "candidates"],
            )
        },
        output_schema: schemas::patch_detect_ports,
        run: patch::detect_ports,
    },
    Tool {
        id: "patch.propose",
        version: 1,
        action: Action::WriteFiles,
        summary: "Propose the patches of a bare closed STL (BC-03): split it at feature edges into faces, say what each is (plane, cylinder with axis and radius, sphere, freeform), mark the openings (pipe ends) and the planes on the bounding box (`xmin` … `zmax`), and write a map with `opening_1`, …, the box sides and `walls`. Openings and box sides get no role: set it with patch.set_role.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"surface": schemas::path(), "output": schemas::path(),
                       "feature_angle_deg": {"type": "number", "minimum": 1, "maximum": 90},
                       "plane_angle_deg": {"type": "number", "minimum": 0, "maximum": 10},
                       "fit_tolerance": {"type": "number", "minimum": 0, "maximum": 0.5},
                       "merge_walls": {"type": "boolean"}}),
                &[
                    "feature_angle_deg",
                    "plane_angle_deg",
                    "fit_tolerance",
                    "merge_walls",
                ],
            )
        },
        output_schema: schemas::patch_propose,
        run: patch::propose,
    },
    Tool {
        id: "bc.skeleton",
        version: 1,
        action: Action::WriteFiles,
        summary: "Write a FoehnCraft case skeleton (TOML) from a patch map (BC-05): one `[boundary.<patch>]` per patch with the FoehnCraft type of its role and `TODO` comments for the values FoehnCraft needs. KoflCraft sets no physical value; a patch without a role is left without a type.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"surface": schemas::path(), "patches": schemas::path(), "output": schemas::path(),
                       "mesh": {"type": "string", "minLength": 1}}),
                &["mesh"],
            )
        },
        output_schema: schemas::bc_skeleton,
        run: patch::bc_skeleton,
    },
    Tool {
        id: "regions.propose",
        version: 1,
        action: Action::WriteFiles,
        summary: "Propose where block templates fit (block-dominant meshing): build a tet scaffold of a closed STL, solve the wall distance (heat method) and the inlet-to-outlet potential, find the skeleton with its end and branch points, and mark each region as sweep, junction, obstacle, box or remainder (polyhedral). A second pass tries each box and remainder region again along its own long axis. Writes the report, optionally the fields (VTU) and a picture (SVG); meshes nothing. The map needs an inlet and an outlet.",
        legacy_name: None,
        outputs: &["output", "fields", "picture"],
        input_schema: || {
            schemas::input(
                json!({"surface": schemas::path(), "patches": schemas::path(), "output": schemas::path(),
                       "fields": schemas::path(), "picture": schemas::path(),
                       "size": {"type": "number", "exclusiveMinimum": 0},
                       "heat_m": {"type": "number", "exclusiveMinimum": 0, "maximum": 20},
                       "second_pass": {"type": "boolean"}}),
                &["fields", "picture", "size", "heat_m", "second_pass"],
            )
        },
        output_schema: || json!({"type": "object", "description": "`koflcraft.region-proposal`"}),
        run: regions::propose,
    },
    Tool {
        id: "patch.split",
        version: 1,
        action: Action::WriteFiles,
        summary: "Split a patch into its regions, separated by feature edges sharper than `feature_angle_deg` (or only by connectivity), named `names` or `<patch>_1`, `<patch>_2`, ….",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::patch_split,
        output_schema: schemas::patch_result,
        run: patch::split,
    },
    Tool {
        id: "patch.merge",
        version: 1,
        action: Action::WriteFiles,
        summary: "Merge the patches in `merge` into one patch `name` (one of them or a new name) with `role` or the first one's role.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::patch_merge,
        output_schema: schemas::patch_result,
        run: patch::merge,
    },
    Tool {
        id: "patch.rename",
        version: 1,
        action: Action::WriteFiles,
        summary: "Rename a patch; the new name must be free.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"surface": schemas::path(), "patches": schemas::path(), "output": schemas::path(),
                       "from": {"type": "string"}, "to": {"type": "string"}}),
                &[],
            )
        },
        output_schema: schemas::patch_result,
        run: patch::rename,
    },
    Tool {
        id: "patch.set_role",
        version: 1,
        action: Action::WriteFiles,
        summary: "Set the role of a patch (wall, inlet, outlet, symmetry, interface, unspecified), which maps to a FoehnCraft boundary-condition type and an OpenFOAM patch type (BC-04).",
        legacy_name: None,
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"surface": schemas::path(), "patches": schemas::path(), "output": schemas::path(),
                       "name": {"type": "string"}, "role": schemas::role()}),
                &[],
            )
        },
        output_schema: schemas::patch_result,
        run: patch::set_role,
    },
    Tool {
        id: "mesh.generate",
        version: 1,
        action: Action::StartBuild,
        summary: "Build and validate a clipped-Voronoi box mesh and save it as a new JSON artifact.",
        legacy_name: Some("generate_mesh"),
        outputs: &["artifact"],
        input_schema: || {
            schemas::input(
                json!({"config": schemas::mesh_config(), "artifact": schemas::path()}),
                &[],
            )
        },
        output_schema: out!("generate_mesh"),
        run: mesh::generate,
    },
    Tool {
        id: "mesh.generate_stl",
        version: 1,
        action: Action::StartBuild,
        summary: "Generate a mesh inside a convex STL surface and save it with its source geometry. With a patch map (`patches`, from the `patch.*` tools) its names and roles name the boundary faces in every export.",
        legacy_name: Some("generate_stl_mesh"),
        outputs: &["artifact"],
        input_schema: || {
            schemas::input(
                json!({"surface": schemas::path(), "config": schemas::stl_config(), "artifact": schemas::path(), "patches": schemas::path()}),
                &["patches"],
            )
        },
        output_schema: out!("generate_stl_mesh"),
        run: mesh::generate_stl,
    },
    Tool {
        id: "mesh.generate_sweep",
        version: 1,
        action: Action::StartBuild,
        summary: "Generate a graded polyhedral block with one inlet or a straight pipe with two ports into a new directory. `handoff` picks the FoehnCraft file: the blocked binary foehncraft.vtkhdf (default) or foehncraft.vtu.",
        legacy_name: Some("generate_sweep"),
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"config": schemas::sweep_schema(), "output": schemas::path(),
                       "handoff": {"enum": ["vtu", "vtkhdf"]}}),
                &["handoff"],
            )
        },
        output_schema: out!("generate_sweep"),
        run: mesh::sweep,
    },
    Tool {
        id: "mesh.generate_channel",
        version: 1,
        action: Action::StartBuild,
        summary: "Generate and check a layered rectangular channel into a new directory (VTU, Gmsh, FoehnCraft files). `handoff` picks the FoehnCraft file: the blocked binary foehncraft.vtkhdf (default) or foehncraft.vtu.",
        legacy_name: Some("generate_channel"),
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"config": schemas::channel_config(), "output": schemas::path(),
                       "handoff": {"enum": ["vtu", "vtkhdf"]}}),
                &["handoff"],
            )
        },
        output_schema: out!("generate_channel"),
        run: mesh::channel,
    },
    Tool {
        id: "grading.add",
        version: 1,
        action: Action::WriteFiles,
        summary: "Add a named size source to a duct partition (GRD-02, GRD-03): a datum object (`kind: datum`: point, segment, plane, box, slab, sphere, cylinder or cone with a `size` within `distance`, growing away at the build's rate) or a coordinate region (`kind: region`: box, slab or cylinder with a `size` inside), optionally in a new `frame`. The wall and the core follow min(own sizes, sources); the build report gives each source's achieved cell sizes.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: || {
            let field = schemas::size_field();
            let sources = &field["properties"]["sources"]["items"]["oneOf"];
            schemas::input(
                json!({"partition": schemas::path(), "output": schemas::path(),
                       "source": {"oneOf": [sources[0], sources[1]]},
                       "frame": field["properties"]["frames"]["items"]}),
                &["frame"],
            )
        },
        output_schema: schemas::grading_result,
        run: grading::add,
    },
    Tool {
        id: "grading.remove",
        version: 1,
        action: Action::WriteFiles,
        summary: "Remove the size source `name` from a duct partition; frames no source uses any more go too.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"partition": schemas::path(), "output": schemas::path(), "name": {"type": "string"}}),
                &[],
            )
        },
        output_schema: schemas::grading_result,
        run: grading::remove,
    },
    Tool {
        id: "grading.list",
        version: 1,
        action: Action::Read,
        summary: "List the size sources and frames of a duct partition.",
        legacy_name: None,
        outputs: &[],
        input_schema: || schemas::input(json!({"partition": schemas::path()}), &[]),
        output_schema: schemas::grading_listing,
        run: grading::list,
    },
    Tool {
        id: "shell.add",
        version: 1,
        action: Action::WriteFiles,
        summary: "Add a zero-thickness shell (baffle, plate, fin) to a duct partition (SHL-01): a planar convex `polygon` (corners counter-clockwise about the side `name`; it may reach past the wall), optionally its cell `size` and free-edge `edge_size`. Its two sides become the wall patches `<name>` and `<name>_back`. A repeated name, a bad polygon and a shell touching another are refused.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"partition": schemas::path(), "output": schemas::path(), "shell": schemas::shell_config()}),
                &[],
            )
        },
        output_schema: schemas::shell_result,
        run: shell::add,
    },
    Tool {
        id: "shell.remove",
        version: 1,
        action: Action::WriteFiles,
        summary: "Remove the shell `name` from a duct partition.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"partition": schemas::path(), "output": schemas::path(), "name": {"type": "string"}}),
                &[],
            )
        },
        output_schema: schemas::shell_result,
        run: shell::remove,
    },
    Tool {
        id: "shell.list",
        version: 1,
        action: Action::Read,
        summary: "List the shells of a duct partition: name, sides, polygon, normal, area and sizes. With the partition's STL as `surface`, the default sizes are filled in as a build would choose them.",
        legacy_name: None,
        outputs: &[],
        input_schema: || {
            schemas::input(
                json!({"partition": schemas::path(), "surface": schemas::path()}),
                &["surface"],
            )
        },
        output_schema: schemas::shell_listing,
        run: shell::list,
    },
    Tool {
        id: "mesh.setup_duct",
        version: 1,
        action: Action::WriteFiles,
        summary: "Write the partition of a polyhedral duct for an STL and its patch map (one inlet, one outlet, walls): a Rust-built scaffold and wall layers sized from the smaller port radius, scaled by `size_factor` (smaller is finer). For `mesh.generate_duct` in polyhedral mode.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"surface": schemas::path(), "patches": schemas::path(), "output": schemas::path(),
                       "size_factor": {"type": "number", "exclusiveMinimum": 0},
                       "layers": {"type": "integer", "minimum": 1},
                       "first_thickness": {"type": "number", "exclusiveMinimum": 0},
                       "growth": {"type": "number", "exclusiveMinimum": 0}}),
                &["size_factor", "layers", "first_thickness", "growth"],
            )
        },
        output_schema: schemas::duct_setup_result,
        run: mesh::setup_duct,
    },
    Tool {
        id: "mesh.generate_duct",
        version: 1,
        action: Action::StartBuild,
        summary: "Generate and validate a partitioned STL bulk, wall-layer or polyhedral duct volume into a new directory. `handoff` picks the FoehnCraft file: the blocked binary foehncraft.vtkhdf (default) or foehncraft.vtu.",
        legacy_name: Some("generate_duct"),
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"surface": schemas::path(), "patches": schemas::path(), "partition": schemas::path(),
                       "mode": {"enum": ["bulk", "layered", "polyhedral"]}, "output": schemas::path(),
                       "handoff": {"enum": ["vtu", "vtkhdf"]}}),
                &["handoff"],
            )
        },
        output_schema: out!("generate_duct"),
        run: mesh::duct,
    },
    Tool {
        id: "mesh.check",
        version: 1,
        action: Action::Read,
        summary: "Check a saved artifact and its mesh geometry.",
        legacy_name: Some("check_mesh"),
        outputs: &[],
        input_schema: || schemas::input(json!({"artifact": schemas::path()}), &[]),
        output_schema: out!("check_mesh"),
        run: mesh::check,
    },
    Tool {
        id: "mesh.report",
        version: 2,
        action: Action::Read,
        summary: "Report build evidence and per-cell and per-face metrics of a saved artifact. Over MCP (or with `cursor` or `page_size`) the result is a summary with one page of rows (`page`); follow `page.next_cursor` for the rest.",
        legacy_name: Some("mesh_report"),
        outputs: &[],
        input_schema: || {
            schemas::input(
                json!({"artifact": schemas::path(), "cursor": {"type": "string"},
                       "page_size": {"type": "integer", "minimum": 1}}),
                &["cursor", "page_size"],
            )
        },
        output_schema: out!("mesh_report"),
        run: mesh::report,
    },
    Tool {
        id: "mesh.export",
        version: 1,
        action: Action::WriteFiles,
        summary: "Validate a saved artifact and export it to a new file: vtu, gmsh (MSH 4.1), gmsh-2.2, foehncraft (single-file VTU), foehncraft-vtkhdf (the same as blocked binary VTKHDF) or openfoam (polyMesh directory).",
        legacy_name: Some("export_mesh"),
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"artifact": schemas::path(), "output": schemas::path(),
                       "format": {"enum": ["vtu", "gmsh", "gmsh-2.2", "foehncraft", "foehncraft-vtkhdf", "openfoam"]}}),
                &[],
            )
        },
        output_schema: out!("export_mesh"),
        run: mesh::export,
    },
    Tool {
        id: "block.template",
        version: 1,
        action: Action::WriteFiles,
        summary: "Write a block topology (`koflcraft.block-topology`, BS-02/BS-03) for a standard case: `channel` (one block, optional two-sided wall clustering) `pipe_ogrid` (H-core and four O-grid blocks, optional wall spacing) `elbow` (the pipe O-grid swept through a bend with optional straight legs, walls projected onto cylinders and the torus) `airfoil_cgrid` (an extruded C-grid of six blocks around a NACA 4-digit section or a Selig/Lednicer coordinate file, the wall a spline through the points) `multi_element_airfoil` (an extruded H-grid around a main element with an optional slat and flap, each a NACA section or coordinate file placed by chord, deflection, gap and a negative overlap; its layout is optimised and triangle blocks are split where that measurably helps), `triangle` (a triangle extruded along +z as three hexahedra around its centroid) or `tetrahedron` (four hexahedra around the centroid; both with cells per corner, MB1h). The file is checked before it is written.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"config": schemas::block_template(), "output": schemas::path()}),
                &[],
            )
        },
        output_schema: || json!({"type": "object", "description": "`koflcraft.block-result` with the written file and its check"}),
        run: block::template,
    },
    Tool {
        id: "block.check",
        version: 1,
        action: Action::Read,
        summary: "Check a block topology before any grid is built (BS-05): indices, right-handed blocks, edge curves and laws, one cell count per chain of opposite edges, faces shared by at most two blocks, every boundary face in exactly one patch. Reports cells, points and the chains, or every issue with the blocks, edge or face it concerns.",
        legacy_name: None,
        outputs: &[],
        input_schema: || schemas::input(json!({"topology": schemas::path()}), &[]),
        output_schema: || json!({"type": "object", "description": "`koflcraft.block-result` with `check`"}),
        run: block::check,
    },
    Tool {
        id: "block.grid",
        version: 1,
        action: Action::StartBuild,
        summary: "Build a block-structured hexahedral mesh from a topology (BS-09): `algorithm` `{\"kind\": \"tfi\"}` (default, transfinite interpolation) or `{\"kind\": \"elliptic\"}` (TFI, then Thompson–Thames–Mastin smoothing of every node off the boundary, across block faces, with orthogonality and first-spacing control at wall patches unless `wall_orthogonality` is false; `wall_sliding: true` lets the wall nodes slide along curved walls (section contours, cylinder and torus surfaces), corners sharper than 30° fixed). Validate it (positive corner Jacobians, conservation, conforming closed cells, BS-13) and write `volume.vtu`, `mesh.msh` (MSH 4.1), `foehncraft.vtu`, `report.json` and `mesh-evidence.json` (with each block's cell range) into a new directory (BS-14).",
        legacy_name: None,
        outputs: &["output"],
        input_schema: || {
            schemas::input(
                json!({"topology": schemas::path(), "output": schemas::path(),
                       "algorithm": schemas::block_algorithm()}),
                &["algorithm"],
            )
        },
        output_schema: || json!({"type": "object", "description": "`koflcraft.block-result` with the `koflcraft.block-report`"}),
        run: block::grid,
    },
    Tool {
        id: "block.list",
        version: 1,
        action: Action::Read,
        summary: "List a block topology for editing (BS-04): vertices with their points, blocks with vertices and cell counts, edges with curve and law, patches with role and face count, surfaces, and the check with every chain of opposite edges.",
        legacy_name: None,
        outputs: &[],
        input_schema: || schemas::input(json!({"topology": schemas::path()}), &[]),
        output_schema: || json!({"type": "object", "description": "`koflcraft.block-result` with vertices, blocks, edges, patches, surfaces, reference and check"}),
        run: block::list,
    },
    Tool {
        id: "block.set_count",
        version: 1,
        action: Action::WriteFiles,
        summary: "Set the cell count of the chain of opposite edges through one edge (given by its two vertex indices): every block edge in the chain gets the count, so the topology stays consistent (BS-05). Writes a new checked topology.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::block_set_count,
        output_schema: schemas::block_edit_result,
        run: block::set_count,
    },
    Tool {
        id: "block.set_distribution",
        version: 1,
        action: Action::WriteFiles,
        summary: "Set the law of one block edge: `uniform`, `geometric` (`first` or `last` spacing) or `tanh` (both), read from the edge's first vertex to its second as given; it is stored in the file's edge direction. Writes a new checked topology.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::block_set_distribution,
        output_schema: schemas::block_edit_result,
        run: block::set_distribution,
    },
    Tool {
        id: "block.set_edge_curve",
        version: 1,
        action: Action::WriteFiles,
        summary: "Set the curve of one block edge: `line`, `arc` through a point, or `spline` through points from the edge's first vertex to its second. The template's analytic reference volume is removed. Writes a new checked topology.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::block_set_edge_curve,
        output_schema: schemas::block_edit_result,
        run: block::set_edge_curve,
    },
    Tool {
        id: "block.move_vertex",
        version: 1,
        action: Action::WriteFiles,
        summary: "Move one vertex `to` a point or `by` an offset. A move that turns a block left-handed or takes a projected face's corner off its surface is refused and nothing is written. The analytic reference volume is removed. Writes a new checked topology.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::block_move_vertex,
        output_schema: schemas::block_edit_result,
        run: block::move_vertex,
    },
    Tool {
        id: "block.add_vertex",
        version: 1,
        action: Action::WriteFiles,
        summary: "Append a vertex to a block topology; its index is the old vertex count. Writes a new checked topology.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::block_add_vertex,
        output_schema: schemas::block_edit_result,
        run: block::add_vertex,
    },
    Tool {
        id: "block.remove_vertex",
        version: 1,
        action: Action::WriteFiles,
        summary: "Remove a vertex that no block, edge, patch or surface uses; higher vertex indices move down by one. Writes a new checked topology.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::block_remove_vertex,
        output_schema: schemas::block_edit_result,
        run: block::remove_vertex,
    },
    Tool {
        id: "block.add_block",
        version: 1,
        action: Action::WriteFiles,
        summary: "Append a hexahedral block by its eight vertices and cell counts. Its faces that meet another block's face become shared and leave their patches (a patch left empty is removed); its other faces go to `patch` (created with `role` when new). Counts must match the chains it joins (BS-05). Writes a new checked topology.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::block_add_block,
        output_schema: schemas::block_edit_result,
        run: block::add_block,
    },
    Tool {
        id: "block.remove_block",
        version: 1,
        action: Action::WriteFiles,
        summary: "Remove a block. Faces it shared with neighbours become boundary faces of `patch`; its own boundary faces leave their patches and surfaces; edge entries no block uses any more are dropped. Writes a new checked topology.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::block_remove_block,
        output_schema: schemas::block_edit_result,
        run: block::remove_block,
    },
    Tool {
        id: "block.split_block",
        version: 1,
        action: Action::WriteFiles,
        summary: "Split a block along an axis after `at` cells. Every block of the chain through that axis is split too, so faces stay point-matched (BS-05). New vertices sit on the TFI grid nodes, and the halves' curves and laws are fitted to the grid (arcs stay arcs), so the grid keeps its nodes. Writes a new checked topology.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::block_split_block,
        output_schema: schemas::block_edit_result,
        run: block::split_block,
    },
    Tool {
        id: "block.split_triangle",
        version: 1,
        action: Action::WriteFiles,
        summary: "Split a block whose cross-section is a triangle (a corner of at least 120° at `vertex`) into three blocks around a centre of valence 3, the tri-center map (MB1h): the side from corner i to corner j gets n_i + n_j cells. The new count enters as a row of `cells` cells (from the triangle's spacing when absent) that runs from two of the triangle's sides through the neighbouring blocks to the boundary, cutting each in two on its own grid. In 3D the flat corner is a block edge (from `vertex` to `along`, the flattest edge at `vertex` when absent): every block stacked along it becomes three too and the row is a sheet through the neighbours. Writes a new checked topology.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::block_split_triangle,
        output_schema: schemas::block_edit_result,
        run: block::split_triangle,
    },
    Tool {
        id: "block.optimize_layout",
        version: 1,
        action: Action::WriteFiles,
        summary: "Optimise the layout of an extruded block topology (MB1g): blocks on wall patches stay; other vertices move (boundary ones slide on straight boundaries), interior edges bend and end spacings scale so that the TFI grid comes closer to orthogonal. With `triangle_splits` (default true) every block corner of 150° or more is tried as a triangle split (`block.split_triangle`) and kept only when its blocks' worst angle drops by 2° or more and the whole grid gets no worse (worst angle, cells above 70°, folds, growth). Reports the non-orthogonality before and after and every split kept or skipped with why. Writes a new checked topology.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::block_optimize_layout,
        output_schema: || json!({"type": "object", "description": "`koflcraft.block-result` with the written file, its check and `layout` (non-orthogonality before and after, triangle splits kept and skipped)"}),
        run: block::optimize_layout,
    },
    Tool {
        id: "block.set_patch",
        version: 1,
        action: Action::WriteFiles,
        summary: "Rename a patch of a block topology and/or set its role (wall, inlet, outlet, symmetry, interface, unspecified); the faces stay. Writes a new checked topology.",
        legacy_name: None,
        outputs: &["output"],
        input_schema: schemas::block_set_patch,
        output_schema: schemas::block_edit_result,
        run: block::set_patch,
    },
    Tool {
        id: "mesh.estimate",
        version: 1,
        action: Action::Read,
        summary: "Estimate a build before it runs (AG-15): the cell count of a start-build tool's input (exact unless `cells_basis` says otherwise: polyhedral ducts), rough memory and output bytes, and the files it would write. `kofl --check COMMAND` prints it instead of running the build.",
        legacy_name: None,
        outputs: &[],
        input_schema: || {
            schemas::input(
                json!({"tool": {"type": "string"}, "input": {"type": "object"}}),
                &[],
            )
        },
        output_schema: || {
            schemas::object_schema(json!({
                "format": {"const": "koflcraft.estimate"}, "format_version": {"type": "string"},
                "tool": {"type": "string"}, "cells": {"type": "integer", "minimum": 0},
                "cells_basis": {"type": "string"}, "memory_bytes": {"type": "integer", "minimum": 0},
                "output_bytes": {"type": "integer", "minimum": 0}, "bytes_basis": {"type": "string"},
                "files": {"type": "array", "items": {"type": "string"}}
            }))
        },
        run: estimate::estimate,
    },
    Tool {
        id: "journal.replay",
        version: 1,
        action: Action::StartBuild,
        summary: "Run the calls of a journal again and report each one whose status, exit code or result digest differs from its record (AG-13). Outputs of the journal must not exist yet.",
        legacy_name: None,
        outputs: &[],
        input_schema: || {
            schemas::input(
                json!({"script": schemas::path(), "keep_going": {"type": "boolean"}}),
                &["keep_going"],
            )
        },
        output_schema: || json!({"type": "object", "description": "`koflcraft.replay`"}),
        run: journal::replay,
    },
];

/// The tool with this id (or legacy `mesh-mcp` name).
pub fn find(id: &str) -> Option<&'static Tool> {
    TOOLS
        .iter()
        .find(|t| t.id == id || t.legacy_name == Some(id))
}

fn nearest(id: &str) -> Option<&'static str> {
    // the id with the longest common prefix, as a hint
    TOOLS
        .iter()
        .map(|t| {
            let common =
                t.id.chars()
                    .zip(id.chars())
                    .take_while(|(a, b)| a == b)
                    .count();
            (common, t.id)
        })
        .max_by_key(|&(c, _)| c)
        .filter(|&(c, _)| c > 0)
        .map(|(_, id)| id)
}

/// The result envelope of a call.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct Envelope {
    pub format: String,
    pub format_version: String,
    pub tool: String,
    pub tool_version: u32,
    /// `finished` or `error`.
    pub status: String,
    pub exit_code: i32,
    pub result: Value,
    pub diagnostics: Vec<Diagnostic>,
    pub limitations: Vec<String>,
    #[serde(skip)]
    pub seconds: f64,
}

impl Envelope {
    pub fn error(tool: &str, version: u32, d: Diagnostic) -> Self {
        Self {
            format: RESULT_FORMAT.into(),
            format_version: FORMAT_VERSION.into(),
            tool: tool.into(),
            tool_version: version,
            status: "error".into(),
            exit_code: 1,
            result: Value::Null,
            diagnostics: vec![d],
            limitations: Vec::new(),
            seconds: 0.0,
        }
    }

    pub fn ok(&self) -> bool {
        self.status == "finished"
    }
}

/// Run a tool by id with a JSON input.
pub fn call(id: &str, input: Value, ctx: &Context) -> Envelope {
    let Some(tool) = find(id) else {
        let hint = nearest(id).map_or(String::new(), |n| format!(" (nearest: {n})"));
        return Envelope::error(
            id,
            0,
            Diagnostic::new("unknown_tool", format!("no tool {id}{hint}")),
        );
    };
    let start = Instant::now();
    let input = if input.is_null() { json!({}) } else { input };
    let mut envelope = match (tool.run)(input, ctx) {
        Ok(result) => Envelope {
            format: RESULT_FORMAT.into(),
            format_version: FORMAT_VERSION.into(),
            tool: tool.id.into(),
            tool_version: tool.version,
            status: "finished".into(),
            exit_code: 0,
            limitations: mesh::limitations(tool.id, &result),
            result,
            diagnostics: Vec::new(),
            seconds: 0.0,
        },
        Err(mut d) => {
            let report = d.report.take();
            let mut e = Envelope::error(tool.id, tool.version, d);
            if let Some(r) = report {
                e.result = *r;
            }
            e
        }
    };
    envelope.seconds = start.elapsed().as_secs_f64();
    envelope
}

/// `system.info`.
pub fn system_info() -> Value {
    json!({
        "format": "koflcraft.info",
        "format_version": FORMAT_VERSION,
        "version": env!("CARGO_PKG_VERSION"),
        "commit": option_env!("KOFL_COMMIT").unwrap_or("unknown"),
        "platform": format!("{}-{}", std::env::consts::OS, std::env::consts::ARCH),
        "debug_assertions": cfg!(debug_assertions),
    })
}

/// `system.commands`: the catalogue (`koflcraft.tool-catalogue`).
pub fn catalogue() -> Value {
    json!({
        "format": CATALOGUE_FORMAT,
        "format_version": FORMAT_VERSION,
        "tools": TOOLS.iter().map(|t| json!({
            "id": t.id,
            "version": t.version,
            "action": t.action,
            "summary": t.summary,
            "legacy_name": t.legacy_name,
            "input_schema": t.input_schema(),
            "output_schema": t.output_schema(),
        })).collect::<Vec<_>>()
    })
}

/// `diag.codes` (`koflcraft.codes`).
pub fn codes() -> Value {
    json!({ "format": "koflcraft.codes", "format_version": FORMAT_VERSION, "codes": CODES })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Invariant: ids are unique, dotted `area.verb`, and so are legacy
    /// names; a duplicate id would make dispatch ambiguous.
    #[test]
    fn ids_are_unique_and_dotted() {
        let mut ids: Vec<_> = TOOLS.iter().map(|t| t.id).collect();
        assert!(ids.iter().all(|id| {
            let (a, b) = id.split_once('.').unwrap();
            !a.is_empty()
                && !b.is_empty()
                && id
                    .chars()
                    .all(|c| c.is_ascii_lowercase() || c == '.' || c == '_')
        }));
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids.len(), TOOLS.len());
    }

    /// Invariant: an unknown id and a wrong input end in an error envelope
    /// with TOOL001 / TOOL002, never a panic or a success.
    #[test]
    fn unknown_tools_and_bad_inputs_are_diagnosed() {
        let ctx = Context::cli();
        let e = call("mesh.genrate", json!({}), &ctx);
        assert_eq!(
            (e.status.as_str(), e.diagnostics[0].code),
            ("error", "TOOL001")
        );
        assert!(e.diagnostics[0].message.contains("mesh.generate"));
        let e = call("mesh.check", json!({"artifact": 3}), &ctx);
        assert_eq!(e.diagnostics[0].code, "TOOL002");
        let e = call("mesh.check", json!({"artifact": "a", "extra": 1}), &ctx);
        assert_eq!(e.diagnostics[0].code, "TOOL002");
    }

    /// Invariant: every input schema is a closed object whose required
    /// fields exist, so the MCP listing and `docs/tools.md` are usable.
    #[test]
    fn schemas_are_closed_objects() {
        for t in TOOLS {
            let s = t.input_schema();
            assert_eq!(s["type"], "object", "{}", t.id);
            assert_eq!(s["additionalProperties"], false, "{}", t.id);
            for r in s["required"].as_array().unwrap() {
                assert!(s["properties"][r.as_str().unwrap()].is_object(), "{}", t.id);
            }
        }
    }
}
