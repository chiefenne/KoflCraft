//! The block-structured tools (`block.*`, BS-03, BS-04, BS-05, BS-14): a
//! template writes a topology file, `block.check` reads one and reports
//! counts and chains or every issue, `block.grid` builds and validates the
//! hexahedral mesh into a new directory. The edit tools (`block.set_*`,
//! `block.move_vertex`) read a topology, apply one edit and write a new,
//! checked topology, so every step is journalled and replayable (BS-04);
//! `block.list` only reads.

use mesh_engine::block::{
    self, CurveFile, DistributionFile, TemplateConfig, TopologyFile, VolumeAlgorithm,
};
use mesh_engine::block_edit::{self, Edit};
use mesh_engine::patches::PatchRole;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::Context;
use crate::diag::Diagnostic;
use crate::mesh::{Run, decode, json_file};

pub const BLOCK_RESULT_FORMAT: &str = "koflcraft.block-result";

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct TemplateInput {
    pub config: TemplateConfig,
    pub output: String,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct CheckInput {
    pub topology: String,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct GridInput {
    pub topology: String,
    pub output: String,
    /// The volume algorithm; TFI when absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub algorithm: Option<VolumeAlgorithm>,
}

pub(crate) fn topology(ctx: &Context, path: &str) -> Result<TopologyFile, Diagnostic> {
    json_file(ctx, path)
}

pub(crate) fn template(input: Value, ctx: &Context) -> Run {
    let i: TemplateInput = decode(input)?;
    let mut texts = std::collections::BTreeMap::new();
    for path in i.config.coordinate_files() {
        let text = String::from_utf8(ctx.read(path)?)
            .map_err(|e| Diagnostic::new("invalid_arguments", e).with_path(path))?;
        texts.insert(path.to_string(), text);
    }
    let file = block::template_with(&i.config, &texts)?;
    write_topology(ctx, &i.output, &file)
}

pub(crate) fn check(input: Value, ctx: &Context) -> Run {
    let i: CheckInput = decode(input)?;
    let file = topology(ctx, &i.topology)?;
    let check = block::check(&file).map_err(|e| Diagnostic::from(e).with_path(&i.topology))?;
    Ok(json!({
        "format": BLOCK_RESULT_FORMAT,
        "format_version": crate::FORMAT_VERSION,
        "topology": i.topology,
        "check": check,
    }))
}

pub(crate) fn grid(input: Value, ctx: &Context) -> Run {
    let i: GridInput = decode(input)?;
    let file = topology(ctx, &i.topology)?;
    let cells = block::cell_count(&file).map_err(|e| Diagnostic::from(e).with_path(&i.topology))?;
    ctx.limit("block cells", cells, |l| l.max_build_cells)?;
    let (dir, name) = ctx.target(&i.output)?;
    let report = block::generate(&file, &i.algorithm.unwrap_or_default(), &dir, &name)?;
    Ok(json!({
        "format": BLOCK_RESULT_FORMAT,
        "format_version": crate::FORMAT_VERSION,
        "output": i.output,
        "report": report,
    }))
}

/// Write a topology and report the file with its check.
fn write_topology(ctx: &Context, output: &str, file: &TopologyFile) -> Result<Value, Diagnostic> {
    let check = block::check(file)?;
    let mut bytes =
        serde_json::to_vec_pretty(file).map_err(|e| Diagnostic::new("serialization_failed", e))?;
    bytes.push(b'\n');
    let (dir, name) = ctx.target(output)?;
    dir.write_new(&name, &bytes)
        .map_err(|e| Diagnostic::new("io_error", e).with_path(output))?;
    let mut result = json!({
        "format": BLOCK_RESULT_FORMAT,
        "format_version": crate::FORMAT_VERSION,
        "output": output,
        "files": [{ "path": output, "blake3": mesh_engine::file_digest(&bytes) }],
        "check": check,
    });
    if let Some(layout) = file.reference.as_ref().and_then(|r| r.layout.as_ref()) {
        result["layout"] = layout.clone();
    }
    Ok(result)
}

pub(crate) fn list(input: Value, ctx: &Context) -> Run {
    let i: CheckInput = decode(input)?;
    let file = topology(ctx, &i.topology)?;
    let mut listed =
        block_edit::list(&file).map_err(|e| Diagnostic::from(e).with_path(&i.topology))?;
    listed["format"] = json!(BLOCK_RESULT_FORMAT);
    listed["format_version"] = json!(crate::FORMAT_VERSION);
    listed["topology"] = json!(i.topology);
    Ok(listed)
}

fn edit(ctx: &Context, topology_path: &str, output: &str, e: Edit) -> Run {
    let file = topology(ctx, topology_path)?;
    let (new, change) = block_edit::apply(&file, &e)
        .map_err(|err| Diagnostic::from(err).with_path(topology_path))?;
    let mut result = write_topology(ctx, output, &new)?;
    result["change"] = json!(change);
    Ok(result)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SetCountInput {
    pub topology: String,
    pub output: String,
    pub edge: [usize; 2],
    pub cells: usize,
}

pub(crate) fn set_count(input: Value, ctx: &Context) -> Run {
    let i: SetCountInput = decode(input)?;
    let e = Edit::SetCount {
        edge: i.edge,
        cells: i.cells,
    };
    edit(ctx, &i.topology, &i.output, e)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SetDistributionInput {
    pub topology: String,
    pub output: String,
    pub edge: [usize; 2],
    pub distribution: DistributionFile,
}

pub(crate) fn set_distribution(input: Value, ctx: &Context) -> Run {
    let i: SetDistributionInput = decode(input)?;
    let e = Edit::SetDistribution {
        edge: i.edge,
        distribution: i.distribution,
    };
    edit(ctx, &i.topology, &i.output, e)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SetCurveInput {
    pub topology: String,
    pub output: String,
    pub edge: [usize; 2],
    pub curve: CurveFile,
}

pub(crate) fn set_edge_curve(input: Value, ctx: &Context) -> Run {
    let i: SetCurveInput = decode(input)?;
    let e = Edit::SetCurve {
        edge: i.edge,
        curve: i.curve,
    };
    edit(ctx, &i.topology, &i.output, e)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct MoveVertexInput {
    pub topology: String,
    pub output: String,
    pub vertex: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub to: Option<[f64; 3]>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub by: Option<[f64; 3]>,
}

pub(crate) fn move_vertex(input: Value, ctx: &Context) -> Run {
    let i: MoveVertexInput = decode(input)?;
    let e = Edit::MoveVertex {
        vertex: i.vertex,
        to: i.to,
        by: i.by,
    };
    edit(ctx, &i.topology, &i.output, e)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SetPatchInput {
    pub topology: String,
    pub output: String,
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rename: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub role: Option<PatchRole>,
}

pub(crate) fn set_patch(input: Value, ctx: &Context) -> Run {
    let i: SetPatchInput = decode(input)?;
    let e = Edit::SetPatch {
        name: i.name,
        rename: i.rename,
        role: i.role,
    };
    edit(ctx, &i.topology, &i.output, e)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct AddVertexInput {
    pub topology: String,
    pub output: String,
    pub point: [f64; 3],
}

pub(crate) fn add_vertex(input: Value, ctx: &Context) -> Run {
    let i: AddVertexInput = decode(input)?;
    edit(
        ctx,
        &i.topology,
        &i.output,
        Edit::AddVertex { point: i.point },
    )
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct RemoveVertexInput {
    pub topology: String,
    pub output: String,
    pub vertex: usize,
}

pub(crate) fn remove_vertex(input: Value, ctx: &Context) -> Run {
    let i: RemoveVertexInput = decode(input)?;
    edit(
        ctx,
        &i.topology,
        &i.output,
        Edit::RemoveVertex { vertex: i.vertex },
    )
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct AddBlockInput {
    pub topology: String,
    pub output: String,
    pub vertices: [usize; 8],
    pub cells: [usize; 3],
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub domain: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub patch: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub role: Option<PatchRole>,
}

pub(crate) fn add_block(input: Value, ctx: &Context) -> Run {
    let i: AddBlockInput = decode(input)?;
    let e = Edit::AddBlock {
        vertices: i.vertices,
        cells: i.cells,
        domain: i.domain,
        patch: i.patch,
        role: i.role,
    };
    edit(ctx, &i.topology, &i.output, e)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct RemoveBlockInput {
    pub topology: String,
    pub output: String,
    pub block: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub patch: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub role: Option<PatchRole>,
}

pub(crate) fn remove_block(input: Value, ctx: &Context) -> Run {
    let i: RemoveBlockInput = decode(input)?;
    let e = Edit::RemoveBlock {
        block: i.block,
        patch: i.patch,
        role: i.role,
    };
    edit(ctx, &i.topology, &i.output, e)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SplitBlockInput {
    pub topology: String,
    pub output: String,
    pub block: usize,
    pub axis: usize,
    pub at: usize,
}

pub(crate) fn split_block(input: Value, ctx: &Context) -> Run {
    let i: SplitBlockInput = decode(input)?;
    let e = Edit::SplitBlock {
        block: i.block,
        axis: i.axis,
        at: i.at,
    };
    edit(ctx, &i.topology, &i.output, e)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SplitTriangleInput {
    pub topology: String,
    pub output: String,
    pub block: usize,
    pub vertex: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cells: Option<usize>,
}

pub(crate) fn split_triangle(input: Value, ctx: &Context) -> Run {
    let i: SplitTriangleInput = decode(input)?;
    let e = Edit::SplitTriangle {
        block: i.block,
        vertex: i.vertex,
        cells: i.cells,
    };
    edit(ctx, &i.topology, &i.output, e)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct OptimizeLayoutInput {
    pub topology: String,
    pub output: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub triangle_splits: Option<bool>,
}

pub(crate) fn optimize_layout(input: Value, ctx: &Context) -> Run {
    let i: OptimizeLayoutInput = decode(input)?;
    let file = topology(ctx, &i.topology)?;
    let (new, layout) =
        mesh_engine::block_triangle::optimize_layout(&file, i.triangle_splits.unwrap_or(true))
            .map_err(|e| Diagnostic::from(e).with_path(&i.topology))?;
    let mut result = write_topology(ctx, &i.output, &new)?;
    result["layout"] = layout;
    Ok(result)
}
