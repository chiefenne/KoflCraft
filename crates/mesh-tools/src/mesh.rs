//! The engine operations as tools. Each takes its typed input (decoded from
//! JSON with unknown fields refused), runs inside the [`Context`] and
//! returns the JSON object the command line has always printed.

use std::path::Path;

use mesh_engine::channel::{ChannelConfig, generate_channel};
use mesh_engine::duct::{self, Mode, Partition, PatchMap};
use mesh_engine::sweep::SweepConfig;
use mesh_engine::{
    ExportData, MeshArtifact, MeshConfig, SCHEMA_VERSION, StlMeshConfig, build_mesh,
    build_stl_mesh, check_mesh, export_mesh, inspect_domain, inspect_stl_domain, mesh_report,
    stl_site_count,
};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use sha2::Digest;

use crate::Context;
use crate::diag::Diagnostic;

pub(crate) type Run = Result<Value, Diagnostic>;

pub(crate) fn decode<T: serde::de::DeserializeOwned>(input: Value) -> Result<T, Diagnostic> {
    serde_json::from_value(input).map_err(|e| Diagnostic::new("invalid_arguments", e))
}

pub(crate) fn json_file<T: serde::de::DeserializeOwned>(
    ctx: &Context,
    path: &str,
) -> Result<T, Diagnostic> {
    serde_json::from_slice(&ctx.read(path)?)
        .map_err(|e| Diagnostic::new("invalid_json", e).with_path(path))
}

fn sites(ctx: &Context, n: usize) -> Result<(), Diagnostic> {
    ctx.limit("sites", n, |l| l.max_sites)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ConfigInput {
    pub config: MeshConfig,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct GenerateInput {
    pub config: MeshConfig,
    pub artifact: String,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct StlInput {
    pub surface: String,
    pub config: StlMeshConfig,
    /// A patch map (`koflcraft.patch-map`) whose names replace the solids'.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub patches: Option<String>,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct GenerateStlInput {
    pub surface: String,
    pub config: StlMeshConfig,
    pub artifact: String,
    /// A patch map (`koflcraft.patch-map`) whose names and roles replace
    /// the solids' (BC-01, BC-04).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub patches: Option<String>,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SweepInput {
    pub config: SweepConfig,
    pub output: String,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ChannelInput {
    pub config: ChannelConfig,
    pub output: String,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct DuctInput {
    pub surface: String,
    pub patches: String,
    pub partition: String,
    pub mode: Mode,
    pub output: String,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ArtifactInput {
    pub artifact: String,
}

/// `mesh.report`: the full report, or (with a cursor, a page size or a
/// front end with a page limit) a summary and one page of rows (AG-16).
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ReportInput {
    pub artifact: String,
    /// `next_cursor` of the previous page.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cursor: Option<String>,
    /// Rows per page, at most the front end's limit.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub page_size: Option<usize>,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ExportInput {
    pub artifact: String,
    pub format: String,
    pub output: String,
}

fn artifact(ctx: &Context, path: &str) -> Result<MeshArtifact, Diagnostic> {
    let artifact: MeshArtifact = json_file(ctx, path)?;
    sites(ctx, artifact.config().sites.len())?;
    if let Some(n) = artifact.surface_facet_count() {
        ctx.limit("STL facets", n, |l| l.max_surface_facets)?;
    }
    Ok(artifact)
}

fn patch_map(
    ctx: &Context,
    path: Option<&str>,
) -> Result<Option<mesh_engine::patches::PatchMap>, Diagnostic> {
    path.map(|p| json_file(ctx, p)).transpose()
}

fn facet_limit(ctx: &Context) -> Option<usize> {
    ctx.limits.as_ref().map(|l| l.max_surface_facets)
}

pub(crate) fn inspect(input: Value, _: &Context) -> Run {
    let i: ConfigInput = decode(input)?;
    Ok(json!(inspect_domain(&i.config)?))
}

pub(crate) fn generate(input: Value, ctx: &Context) -> Run {
    let i: GenerateInput = decode(input)?;
    sites(ctx, i.config.sites.len())?;
    let (dir, name) = ctx.target(&i.artifact)?;
    let artifact = build_mesh(i.config)?;
    artifact.write_new(&dir, &name)?;
    Ok(json!({
        "schema_version": SCHEMA_VERSION,
        "artifact": i.artifact,
        "build": artifact.build()
    }))
}

pub(crate) fn inspect_stl(input: Value, ctx: &Context) -> Run {
    let i: StlInput = decode(input)?;
    let bytes = ctx.read(&i.surface)?;
    sites(ctx, stl_site_count(&bytes, &i.config, facet_limit(ctx))?)?;
    let patches = patch_map(ctx, i.patches.as_deref())?;
    Ok(json!(inspect_stl_domain(
        &bytes,
        &i.config,
        facet_limit(ctx),
        patches.as_ref()
    )?))
}

pub(crate) fn generate_stl(input: Value, ctx: &Context) -> Run {
    let i: GenerateStlInput = decode(input)?;
    let bytes = ctx.read(&i.surface)?;
    sites(ctx, stl_site_count(&bytes, &i.config, facet_limit(ctx))?)?;
    let (dir, name) = ctx.target(&i.artifact)?;
    let patches = patch_map(ctx, i.patches.as_deref())?;
    let artifact = build_stl_mesh(&bytes, i.config, facet_limit(ctx), patches.as_ref())?;
    artifact.write_new(&dir, &name)?;
    let mut output = json!({
        "schema_version": SCHEMA_VERSION,
        "artifact_schema_version": artifact.schema_version(),
        "artifact": i.artifact,
        "build": artifact.build()
    });
    if let Some(sizing) = artifact.sizing() {
        output["sizing"] = json!(sizing);
    }
    Ok(output)
}

pub(crate) fn sweep(input: Value, ctx: &Context) -> Run {
    let i: SweepInput = decode(input)?;
    let (dir, name) = ctx.target(&i.output)?;
    let budget = ctx
        .limits
        .as_ref()
        .map_or(usize::MAX / 1024, |l| l.max_sweep_cells);
    let report = mesh_engine::sweep::generate(&i.config, &dir, &name, budget)?;
    Ok(json!({"schema_version": SCHEMA_VERSION, "output": i.output, "report": report}))
}

pub(crate) fn channel(input: Value, ctx: &Context) -> Run {
    let i: ChannelInput = decode(input)?;
    ctx.limit("layers", i.config.layers, |l| l.max_channel_layers)?;
    ctx.limit("axial cells", i.config.axial_cells, |l| {
        l.max_channel_axial_cells
    })?;
    let (dir, name) = ctx.target(&i.output)?;
    let report = generate_channel(&i.config, &dir, &name)?;
    Ok(json!({ "schema_version": SCHEMA_VERSION, "output": i.output, "report": report }))
}

pub(crate) fn duct(input: Value, ctx: &Context) -> Run {
    let i: DuctInput = decode(input)?;
    let (dir, name) = ctx.target(&i.output)?;
    let patches: PatchMap = json_file(ctx, &i.patches)?;
    let partition: Partition = json_file(ctx, &i.partition)?;
    ctx.limit("source facets", patches.facet_count(), |l| {
        l.max_duct_facets
    })?;
    ctx.limit("layers", partition.layers, |l| l.max_duct_layers)?;
    ctx.limit("regions", partition.regions.len(), |l| l.max_duct_regions)?;
    let scaffold = partition
        .scaffold_file(i.mode)
        .map(|file| {
            let path = Path::new(&i.partition).with_file_name(file);
            ctx.read(&path.to_string_lossy())
        })
        .transpose()?;
    let report = duct::generate(
        &ctx.read(&i.surface)?,
        &patches,
        &partition,
        i.mode,
        scaffold.as_deref(),
        &dir,
        &name,
    )?;
    Ok(json!({"schema_version": SCHEMA_VERSION, "output": i.output, "report": report}))
}

pub(crate) fn check(input: Value, ctx: &Context) -> Run {
    let i: ArtifactInput = decode(input)?;
    Ok(json!(check_mesh(&artifact(ctx, &i.artifact)?)?))
}

pub(crate) fn report(input: Value, ctx: &Context) -> Run {
    let i: ReportInput = decode(input)?;
    let report = json!(mesh_report(&artifact(ctx, &i.artifact)?)?);
    let limit = ctx.limits.as_ref().map(|l| l.max_page_rows);
    if i.cursor.is_none() && i.page_size.is_none() && limit.is_none() {
        return Ok(report);
    }
    let size = match (i.page_size, limit) {
        (Some(0), _) => {
            return Err(Diagnostic::new(
                "invalid_arguments",
                "page_size must be at least 1",
            ));
        }
        (Some(n), Some(max)) if n > max => {
            return Err(Diagnostic::new(
                "limit_exceeded",
                format!("{} pages at most {max} rows, got page_size {n}", ctx.via),
            ));
        }
        (Some(n), _) => n,
        (None, Some(max)) => max,
        (None, None) => DEFAULT_PAGE_ROWS,
    };
    crate::page::report_page(report, i.cursor.as_deref(), size)
}

/// Rows per page when neither the input nor the front end sets one.
pub(crate) const DEFAULT_PAGE_ROWS: usize = 500;

pub(crate) fn export(input: Value, ctx: &Context) -> Run {
    let i: ExportInput = decode(input)?;
    let format = i.format.parse()?;
    let (dir, name) = ctx.target(&i.output)?;
    let a = artifact(ctx, &i.artifact)?;
    let data = export_mesh(&a, format)?;
    data.write_new(&dir, &name)?;
    // the digests of what was written, so a journal replay compares content
    let sha = |b: &[u8]| crate::journal::hex(&sha2::Sha256::digest(b));
    let files: Vec<Value> = match &data {
        ExportData::Vtu(b) | ExportData::Gmsh(b) | ExportData::FoehnCraft(b) => {
            vec![json!({ "path": i.output, "sha256": sha(b) })]
        }
        ExportData::Bundle(files) => files
            .iter()
            .map(|(n, b)| json!({ "path": format!("{}/{n}", i.output), "sha256": sha(b) }))
            .collect(),
    };
    Ok(
        json!({ "schema_version": SCHEMA_VERSION, "output": i.output, "format": i.format, "files": files }),
    )
}

/// Limitations of a result (AG-17): what the output does not establish.
pub(crate) fn limitations(id: &str, output: &Value) -> Vec<String> {
    let mut l: Vec<String> = output["limitations"]
        .as_array()
        .map(|a| {
            a.iter()
                .filter_map(|v| v.as_str().map(str::to_string))
                .collect()
        })
        .unwrap_or_default();
    let report = &output["report"];
    match id {
        "mesh.generate_stl" | "domain.inspect_stl" => {
            l.push("convex STL surfaces only".into());
            if let Some(sizing) = output.get("sizing") {
                l.push(
                    "size field: STL surface and closed-body datum objects (GRD-02) are open; curvature, gap and feature-edge sources use Euclidean distances, exact in the convex domains meshed so far".into(),
                );
                for source in sizing["sources"].as_array().into_iter().flatten() {
                    let derived = &source["derived"];
                    let name = source["name"].as_str().unwrap_or("?");
                    if let Some(n) = derived["at_min_size"].as_u64()
                        && n > 0
                    {
                        l.push(format!(
                            "source {name}: {n} surface pieces ask for less than min_size and are held at it, so their target is not met there"
                        ));
                    }
                    if let (Some(below), Some(samples)) = (
                        derived["samples_below_target"].as_u64(),
                        derived["samples"].as_u64(),
                    ) && below > 0
                    {
                        l.push(format!(
                            "source {name}: {below} of {samples} wall samples have fewer cells across the gap than asked"
                        ));
                    }
                }
                if let (Some(over), Some(total)) = (
                    sizing["faces_over_gradation"].as_u64(),
                    sizing["interior_faces"].as_u64(),
                ) && over > 0
                {
                    l.push(format!(
                        "{over} of {total} interior faces join cells whose size ratio exceeds the gradation; the limit bounds the size field, not each cell pair"
                    ));
                }
            }
        }
        "mesh.generate_sweep" => {
            l.push("analytic block and pipe shapes only".into());
            if let Some(e) = report["relative_geometry_volume_error"].as_f64() {
                l.push(format!(
                    "fitted volume differs from the analytic shape by {:.3} %",
                    100.0 * e
                ));
            }
        }
        "mesh.generate_channel" => l.push("planar walls only".into()),
        "mesh.generate_duct" => {
            if let (Some(v), Some(s)) =
                (report["volume"].as_f64(), report["source_volume"].as_f64())
            {
                l.push(format!(
                    "fitted volume differs from the STL by {:.3} %",
                    100.0 * (v / s - 1.0)
                ));
            }
            l.push("geometry outside the validated elbow and Manifold fixtures is untested".into());
        }
        _ => {}
    }
    if id.starts_with("mesh.generate") {
        l.push("flow accuracy on this mesh is untested".into());
    }
    l
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SetupDuctInput {
    pub surface: String,
    pub patches: String,
    pub output: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub size_factor: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub layers: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub first_thickness: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub growth: Option<f64>,
}

/// Write the partition of a polyhedral duct (`mesh.setup_duct`).
pub(crate) fn setup_duct(input: Value, ctx: &Context) -> Run {
    let i: SetupDuctInput = decode(input)?;
    let patches: PatchMap = json_file(ctx, &i.patches)?;
    let options = duct::SetupOptions {
        size_factor: i.size_factor.unwrap_or(1.0),
        layers: i.layers,
        first_thickness: i.first_thickness,
        growth: i.growth,
    };
    let (partition, radius) = duct::setup(&ctx.read(&i.surface)?, &patches, options)?;
    let mut bytes = serde_json::to_vec_pretty(&partition)
        .map_err(|e| Diagnostic::new("serialization_failed", e))?;
    bytes.push(b'\n');
    let (dir, name) = ctx.target(&i.output)?;
    dir.write_new(&name, &bytes)
        .map_err(|e| Diagnostic::new("io_error", e).with_path(&i.output))?;
    Ok(json!({
        "format": "koflcraft.duct-setup",
        "format_version": crate::FORMAT_VERSION,
        "output": i.output,
        "files": [{ "path": i.output, "sha256": crate::journal::hex(&sha2::Sha256::digest(&bytes)) }],
        "port_radius": radius,
        "size_factor": options.size_factor,
        "partition": partition,
        "note": "For mesh.generate_duct in polyhedral mode; the partition has no regions, so bulk and layered modes need a hand-made one.",
    }))
}
