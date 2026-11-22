use std::ffi::OsStr;
use std::fmt;

use mesh_engine::channel::{ChannelConfig, generate_channel};
use mesh_engine::duct::{self, Mode, Partition, PatchMap};
use mesh_engine::{
    EngineError, MeshArtifact, MeshConfig, MeshDirectory, SCHEMA_VERSION, SURFACE_ARTIFACT_VERSION,
    StlMeshConfig, build_mesh, build_stl_mesh, check_mesh, export_mesh, inspect_domain,
    inspect_stl_domain, mesh_report,
};
use serde::Deserialize;
use serde_json::{Value, json};

const MAX_CHANNEL_LAYERS: usize = 8;
const MAX_CHANNEL_AXIAL_CELLS: usize = 32;
const MAX_SITES: usize = 256;
const MAX_FILE_BYTES: usize = 8 * 1024 * 1024;
const MAX_SURFACE_FACETS: usize = 2048;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct InspectArgs {
    config: MeshConfig,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct GenerateArgs {
    config: MeshConfig,
    artifact: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct InspectStlArgs {
    surface: String,
    config: StlMeshConfig,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct GenerateStlArgs {
    surface: String,
    config: StlMeshConfig,
    artifact: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct SweepArgs {
    config: mesh_engine::sweep::SweepConfig,
    output: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct DuctArgs {
    surface: String,
    patches: String,
    partition: String,
    mode: Mode,
    output: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ChannelArgs {
    config: ChannelConfig,
    output: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ArtifactArgs {
    artifact: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ExportArgs {
    artifact: String,
    format: String,
    output: String,
}

type ToolResult = Result<Value, ToolError>;

pub fn call(directory: &MeshDirectory, name: &str, args: Value) -> Option<ToolResult> {
    Some(match name {
        "inspect_domain" => inspect(args),
        "generate_mesh" => generate(directory, args),
        "generate_channel" => channel(directory, args),
        "generate_sweep" => generate_sweep(directory, args),
        "generate_duct" => generate_duct(directory, args),
        "inspect_stl_domain" => inspect_stl(directory, args),
        "generate_stl_mesh" => generate_stl(directory, args),
        "check_mesh" => read_artifact_args(directory, args)
            .and_then(|artifact| Ok(json!(check_mesh(&artifact)?))),
        "mesh_report" => read_artifact_args(directory, args)
            .and_then(|artifact| Ok(json!(mesh_report(&artifact)?))),
        "export_mesh" => export(directory, args),
        _ => return None,
    })
}

fn inspect(args: Value) -> ToolResult {
    let args: InspectArgs = decode(args)?;
    check_site_limit(args.config.sites.len())?;
    Ok(json!(inspect_domain(&args.config)?))
}

fn generate(directory: &MeshDirectory, args: Value) -> ToolResult {
    let args: GenerateArgs = decode(args)?;
    check_site_limit(args.config.sites.len())?;
    let name = entry_name(&args.artifact)?;
    let artifact = build_mesh(args.config)?;
    artifact.write_new(directory, name)?;
    Ok(json!({
        "schema_version": SCHEMA_VERSION,
        "artifact": args.artifact,
        "build": artifact.build()
    }))
}

fn generate_sweep(directory: &MeshDirectory, args: Value) -> ToolResult {
    let args: SweepArgs = decode(args)?;
    let report =
        mesh_engine::sweep::generate(&args.config, directory, entry_name(&args.output)?, 500_000)?;
    Ok(json!({"schema_version":SCHEMA_VERSION,"report":report,"output":args.output}))
}

fn channel(directory: &MeshDirectory, args: Value) -> ToolResult {
    let args: ChannelArgs = decode(args)?;
    let name = entry_name(&args.output)?;
    if args.config.layers > MAX_CHANNEL_LAYERS || args.config.axial_cells > MAX_CHANNEL_AXIAL_CELLS
    {
        return Err(ToolError::new(
            "channel_limit",
            format_args!(
                "MCP accepts at most {MAX_CHANNEL_LAYERS} layers and {MAX_CHANNEL_AXIAL_CELLS} axial cells"
            ),
        ));
    }
    let report = generate_channel(&args.config, directory, name)?;
    Ok(json!({ "schema_version": SCHEMA_VERSION, "output": args.output, "report": report }))
}

fn generate_duct(directory: &MeshDirectory, args: Value) -> ToolResult {
    let args: DuctArgs = decode(args)?;
    let output = entry_name(&args.output)?;
    let read_json = |name: &str| -> Result<Value, ToolError> {
        serde_json::from_slice(&read_file(directory, name)?)
            .map_err(|e| ToolError::new("invalid_json", e))
    };
    let patches: PatchMap = decode(read_json(&args.patches)?)?;
    let partition: Partition = decode(read_json(&args.partition)?)?;
    if patches.facet_count > 100_000 || partition.layers > 4 || partition.regions.len() > 64 {
        return Err(ToolError::new(
            "duct_limit",
            "MCP accepts at most 100000 source facets, 4 layers and 64 regions per duct",
        ));
    }
    let scaffold = partition
        .scaffold_file(args.mode)
        .map(|file| read_file(directory, file))
        .transpose()?;
    let report = duct::generate(
        &read_file(directory, &args.surface)?,
        &patches,
        &partition,
        args.mode,
        scaffold.as_deref(),
        directory,
        output,
    )?;
    Ok(json!({"schema_version":SCHEMA_VERSION,"output":args.output,"report":report}))
}

fn inspect_stl(directory: &MeshDirectory, args: Value) -> ToolResult {
    let args: InspectStlArgs = decode(args)?;
    check_site_limit(args.config.sites.len())?;
    let bytes = read_file(directory, &args.surface)?;
    Ok(json!(inspect_stl_domain(
        &bytes,
        &args.config,
        Some(MAX_SURFACE_FACETS)
    )?))
}

fn generate_stl(directory: &MeshDirectory, args: Value) -> ToolResult {
    let args: GenerateStlArgs = decode(args)?;
    check_site_limit(args.config.sites.len())?;
    let name = entry_name(&args.artifact)?;
    let bytes = read_file(directory, &args.surface)?;
    let artifact = build_stl_mesh(&bytes, args.config, Some(MAX_SURFACE_FACETS))?;
    artifact.write_new(directory, name)?;
    Ok(json!({
        "schema_version": SCHEMA_VERSION,
        "artifact_schema_version": artifact.schema_version(),
        "artifact": args.artifact,
        "build": artifact.build()
    }))
}

fn export(directory: &MeshDirectory, args: Value) -> ToolResult {
    let args: ExportArgs = decode(args)?;
    let name = entry_name(&args.output)?;
    let format = args.format.parse()?;
    let artifact = read_artifact(directory, &args.artifact)?;
    export_mesh(&artifact, format)?.write_new(directory, name)?;
    Ok(json!({ "schema_version": SCHEMA_VERSION, "output": args.output }))
}

fn read_artifact_args(directory: &MeshDirectory, args: Value) -> Result<MeshArtifact, ToolError> {
    let args: ArtifactArgs = decode(args)?;
    read_artifact(directory, &args.artifact)
}

fn read_artifact(directory: &MeshDirectory, artifact: &str) -> Result<MeshArtifact, ToolError> {
    let bytes = read_file(directory, artifact)?;
    let artifact: MeshArtifact =
        serde_json::from_slice(&bytes).map_err(|error| ToolError::new("invalid_json", error))?;
    check_site_limit(artifact.config().sites.len())?;
    if artifact
        .surface_facet_count()
        .is_some_and(|n| n > MAX_SURFACE_FACETS)
    {
        return Err(ToolError::new(
            "facet_limit",
            format_args!("MCP accepts at most {MAX_SURFACE_FACETS} STL facets"),
        ));
    }
    Ok(artifact)
}

fn read_file(directory: &MeshDirectory, filename: &str) -> Result<Vec<u8>, ToolError> {
    directory
        .read_regular(entry_name(filename)?, MAX_FILE_BYTES)
        .map_err(|error| ToolError::new("io_error", error))
}

fn entry_name(name: &str) -> Result<&OsStr, ToolError> {
    if name.is_empty() || matches!(name, "." | "..") || name.contains(['/', '\\', ':', '\0']) {
        return Err(ToolError::new(
            "invalid_path",
            "expected a single filename within the authorized directory",
        ));
    }
    Ok(OsStr::new(name))
}

fn decode<T: serde::de::DeserializeOwned>(args: Value) -> Result<T, ToolError> {
    serde_json::from_value(args).map_err(|error| ToolError::new("invalid_arguments", error))
}

fn check_site_limit(count: usize) -> Result<(), ToolError> {
    if count > MAX_SITES {
        return Err(ToolError::new(
            "site_limit",
            format_args!("MCP accepts at most {MAX_SITES} sites"),
        ));
    }
    Ok(())
}

pub struct ToolError(Value);

impl ToolError {
    fn new(code: &str, message: impl fmt::Display) -> Self {
        Self(json!({
            "schema_version": SCHEMA_VERSION,
            "error": { "code": code, "message": message.to_string() }
        }))
    }
}

impl From<EngineError> for ToolError {
    fn from(error: EngineError) -> Self {
        Self(json!({ "schema_version": SCHEMA_VERSION, "error": error }))
    }
}

impl fmt::Display for ToolError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

fn sweep_schema() -> Value {
    let positive = json!({"type":"number","exclusiveMinimum":0});
    let growth = json!({"type":"number","minimum":1});
    let block = object_schema(json!({
        "kind":{"const":"block"},"dimensions_mm":{"type":"array","minItems":3,"maxItems":3,"items":positive},
        "inlet":{"enum":["x_min","y_min","z_min"]}
    }));
    let pipe =
        object_schema(json!({"kind":{"const":"pipe"},"diameter_mm":positive,"length_mm":positive}));
    object_schema(
        json!({"schema_version":{"const":SCHEMA_VERSION},"geometry":{"oneOf":[block,pipe]},"core_spacing_mm":positive,
        "wall_layers":{"type":"integer","minimum":1},"first_layer_mm":positive,"layer_growth":growth,"max_neighbour_size_ratio":growth}),
    )
}

pub fn listing() -> Vec<Value> {
    let point =
        json!({ "type": "array", "items": { "type": "number" }, "minItems": 3, "maxItems": 3 });
    let config = json!({
        "type": "object",
        "additionalProperties": false,
        "required": ["schema_version", "domain", "sites"],
        "properties": {
            "schema_version": { "const": SCHEMA_VERSION },
            "domain": {
                "type": "object",
                "additionalProperties": false,
                "required": ["min", "max"],
                "properties": { "min": point, "max": point }
            },
            "sites": { "type": "array", "items": point, "minItems": 1, "maxItems": MAX_SITES },
            "fixture_id": { "type": ["string", "null"] },
            "seed": { "type": ["integer", "null"], "minimum": 0 }
        }
    });
    let filename = json!({ "type": "string", "minLength": 1, "description": "Single filename in the authorized directory" });
    let stl_config = json!({
        "type": "object",
        "additionalProperties": false,
        "required": ["schema_version", "sites"],
        "properties": {
            "schema_version": config["properties"]["schema_version"],
            "sites": config["properties"]["sites"],
            "fixture_id": config["properties"]["fixture_id"],
            "seed": config["properties"]["seed"],
            "facet_patches": { "type": ["array", "null"], "items": { "type": "string" }, "maxItems": MAX_SURFACE_FACETS }
        }
    });
    let positive = json!({"type": "number", "exclusiveMinimum": 0});
    let channel_config = object_schema(json!({
        "schema_version": {"const": SCHEMA_VERSION},
        "length": positive, "width": positive, "height": positive,
        "layers": {"type":"integer", "minimum":1, "maximum":MAX_CHANNEL_LAYERS},
        "axial_cells": {"type":"integer", "minimum":1, "maximum":MAX_CHANNEL_AXIAL_CELLS},
        "first_thickness": positive, "growth": positive
    }));
    vec![
        tool(
            "generate_sweep",
            "Generate a graded polyhedral block with one inlet or a straight pipe with two ports, up to 500000 cells.",
            json!({"config":sweep_schema(),"output":filename}),
            &["config", "output"],
            false,
        ),
        tool(
            "generate_duct",
            "Generate and validate a partitioned STL bulk or wall-layer volume with exposed inlet and outlet patches.",
            json!({"surface":filename,"patches":filename,"partition":filename,"mode":{"enum":["bulk","layered","polyhedral"]},"output":filename}),
            &["surface", "patches", "partition", "mode", "output"],
            false,
        ),
        tool(
            "generate_channel",
            "Generate and check a layered rectangular channel, writing VTU and Gmsh files into a new directory.",
            json!({ "config": channel_config, "output": filename }),
            &["config", "output"],
            false,
        ),
        tool(
            "inspect_domain",
            "Inspect box bounds, sites and reference backend limitations.",
            json!({ "config": config }),
            &["config"],
            true,
        ),
        tool(
            "generate_mesh",
            "Build and validate a mesh and save a new JSON artifact.",
            json!({ "config": config, "artifact": filename }),
            &["config", "artifact"],
            false,
        ),
        tool(
            "inspect_stl_domain",
            "Validate a convex STL shell and report bounds, volume, and named patches.",
            json!({ "surface": filename, "config": stl_config }),
            &["surface", "config"],
            true,
        ),
        tool(
            "generate_stl_mesh",
            "Generate a mesh from a convex STL surface and save it with its source geometry.",
            json!({ "surface": filename, "config": stl_config, "artifact": filename }),
            &["surface", "config", "artifact"],
            false,
        ),
        tool(
            "check_mesh",
            "Check a saved artifact and its mesh geometry.",
            json!({ "artifact": filename }),
            &["artifact"],
            true,
        ),
        tool(
            "mesh_report",
            "Report build evidence and per-cell and per-face metrics.",
            json!({ "artifact": filename }),
            &["artifact"],
            true,
        ),
        tool(
            "export_mesh",
            "Validate and export to a new VTK XML or Gmsh mesh file.",
            json!({ "artifact": filename, "output": filename, "format": { "enum": ["vtu", "gmsh"] } }),
            &["artifact", "format", "output"],
            false,
        ),
    ]
}

fn tool(
    name: &str,
    description: &str,
    properties: Value,
    required: &[&str],
    read_only: bool,
) -> Value {
    json!({
        "name": name,
        "description": description,
        "inputSchema": {
            "type": "object",
            "additionalProperties": false,
            "required": required,
            "properties": properties
        },
        "outputSchema": output_schema(name),
        "annotations": { "readOnlyHint": read_only, "destructiveHint": false, "openWorldHint": false }
    })
}

fn output_schema(name: &str) -> Value {
    let version = json!({ "const": SCHEMA_VERSION });
    let string = json!({ "type": "string" });
    let count = json!({ "type": "integer", "minimum": 0 });
    let number = json!({ "type": "number" });
    let fixture = json!({ "type": ["string", "null"] });
    let seed = json!({ "type": ["integer", "null"], "minimum": 0 });
    let digest = json!({ "type": "string", "pattern": "^[a-f0-9]{64}$" });
    let point = json!({ "type": "array", "items": number, "minItems": 3, "maxItems": 3 });
    let build = object_schema(json!({
        "algorithm": string, "engine_version": string, "operating_system": string,
        "architecture": string, "debug_assertions": { "type": "boolean" }, "elapsed_micros": count,
        "input_digest": digest, "config_digest": digest, "topology_digest": digest, "mesh_digest": digest
    }));
    let check = object_schema(json!({
        "schema_version": version, "valid": { "const": true }, "vertices": count,
        "faces": count, "internal_faces": count, "boundary_faces": count, "cells": count,
        "total_volume": number, "domain_volume": number, "relative_volume_error": number,
        "conservation_tolerance": number, "geometry_tolerance": number, "topology_digest": digest
    }));
    match name {
        "inspect_domain" => object_schema(json!({
            "schema_version": version, "domain": object_schema(json!({ "min": point, "max": point })),
            "site_count": count, "fixture_id": fixture, "seed": seed, "algorithm": string,
            "limitations": { "type": "array", "items": string }
        })),
        "generate_mesh" => {
            object_schema(json!({ "schema_version": version, "artifact": string, "build": build }))
        }
        "inspect_stl_domain" => object_schema(json!({
            "schema_version": version,
            "bounds": object_schema(json!({ "min": point, "max": point })),
            "site_count": count, "facet_count": count, "vertex_count": count, "volume": number,
            "patches": { "type": "array", "items": object_schema(json!({
                "name": string, "facets": count, "area": number
            })) },
            "fixture_id": fixture, "seed": seed, "algorithm": string
        })),
        "generate_duct" => {
            let numbers = json!({"type":"array","items":number});
            let duct = object_schema(json!({
                "schema_version":version,"algorithm":string,"mode":{"enum":["bulk","layered","polyhedral"]},
                "source_sha256":digest,"config_sha256":digest,"regions":count,"source_facets":count,"scaffold_volume":{"type":["number","null"]},"stl_patch_areas":object_schema(json!({"inlet":number,"outlet":number,"walls":number})),
                "tetrahedra":count,"polyhedra":count,"core_sphericity":{"anyOf":[object_schema(json!({"minimum":number,"p05":number,"median":number})),{"type":"null"}]},"polyhedron_face_range":{"type":["array","null"], "items":count,"minItems":2,"maxItems":2},"prism_core_faces":count,"prisms":count,"source_volume":number,"volume":number,
                "minimum_cell_volume":number,"minimum_tet_mean_ratio":{"type":["number","null"]},"tetrahedra_below_mean_ratio_0_1":count,
                "minimum_wedge_jacobian":{"type":["number","null"]},"minimum_kernel_clearance":{"type":["number","null"]},
                "patch_areas":object_schema(json!({"inlet":number,"outlet":number,"walls":number})),
                "port_areas_by_layer":object_schema(json!({"inlet":numbers,"outlet":numbers})),
                "thicknesses":numbers,"normal_thickness_ranges":{"type":"array","items":{"type":"array","minItems":2,"maxItems":2,"items":number}},
                "generation_seconds":number
            }));
            object_schema(json!({"schema_version":version,"output":string,"report":duct}))
        }
        "generate_sweep" => {
            let numbers = json!({"type":"array","items":number});
            object_schema(
                json!({"schema_version":version,"output":string,"report":object_schema(json!({
                    "schema_version":version,"algorithm":string,"config":sweep_schema(),"cells":count,"core_cells":count,"wall_layer_cells":count,
                    "axial_intervals":count,"core_axial_intervals":count,"cross_section_core_cells":count,
                    "core_cell_size_mm":object_schema(json!({"min":number,"median":number,"max":number})),
                    "neighbour_size_ratio":object_schema(json!({"p95_cell_max":number,"max":number,"limit":number})),
                    "max_internal_non_orthogonality_deg":number,"max_internal_skewness":number,"wall_thicknesses_mm":numbers,
                    "inlet_area_by_layer_mm2":numbers,"outlet_area_by_layer_mm2":{"type":["array","null"],"items":number},"volume_mm3":number,
                    "analytic_volume_mm3":number,"relative_geometry_volume_error":number,"profile_sides":count
                }))}),
            )
        }
        "generate_channel" => {
            let numbers = json!({"type":"array","items":number});
            let port = object_schema(json!({
                "faces": count, "core_area": number, "layer_areas": numbers, "total_area": number
            }));
            object_schema(json!({
                "schema_version": version, "output": string,
                "report": object_schema(json!({
                    "cells": count, "prisms": count, "core_cells": count, "total_volume": number,
                    "minimum_cell_volume": number, "wall_distances": numbers, "inlet": port, "outlet": port
                }))
            }))
        }
        "generate_stl_mesh" => object_schema(json!({
            "schema_version": version, "artifact_schema_version": { "const": SURFACE_ARTIFACT_VERSION },
            "artifact": string, "build": build
        })),
        "check_mesh" => check,
        "mesh_report" => {
            let cell = object_schema(json!({
                "centroid": point, "volume": number, "face_count": count, "neighbour_count": count
            }));
            let face = object_schema(json!({
                "centroid": point, "normal": point, "area": number,
                "non_orthogonality_degrees": number, "skewness": number
            }));
            object_schema(json!({
                "schema_version": version, "fixture_id": fixture, "seed": seed,
                "build": build, "validation": check,
                "cells": { "type": "array", "items": cell },
                "faces": { "type": "array", "items": face }
            }))
        }
        "export_mesh" => object_schema(json!({ "schema_version": version, "output": string })),
        _ => unreachable!("the listing declares all tool schemas"),
    }
}

fn object_schema(properties: Value) -> Value {
    let required: Vec<_> = properties.as_object().unwrap().keys().collect();
    json!({
        "type": "object",
        "additionalProperties": false,
        "required": required,
        "properties": properties
    })
}
