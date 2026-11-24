//! The block-structured tools (`block.*`, BS-03, BS-04, BS-05, BS-14): a
//! template writes a topology file, `block.check` reads one and reports
//! counts and chains or every issue, `block.grid` builds and validates the
//! hexahedral mesh into a new directory.

use mesh_engine::block::{self, TemplateConfig, TopologyFile};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use sha2::Digest;

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
}

pub(crate) fn topology(ctx: &Context, path: &str) -> Result<TopologyFile, Diagnostic> {
    json_file(ctx, path)
}

pub(crate) fn template(input: Value, ctx: &Context) -> Run {
    let i: TemplateInput = decode(input)?;
    let file = block::template(&i.config)?;
    let check = block::check(&file)?;
    let mut bytes =
        serde_json::to_vec_pretty(&file).map_err(|e| Diagnostic::new("serialization_failed", e))?;
    bytes.push(b'\n');
    let (dir, name) = ctx.target(&i.output)?;
    dir.write_new(&name, &bytes)
        .map_err(|e| Diagnostic::new("io_error", e).with_path(&i.output))?;
    Ok(json!({
        "format": BLOCK_RESULT_FORMAT,
        "format_version": crate::FORMAT_VERSION,
        "output": i.output,
        "files": [{ "path": i.output, "sha256": crate::journal::hex(&sha2::Sha256::digest(&bytes)) }],
        "check": check,
    }))
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
    let report = block::generate(&file, &dir, &name)?;
    Ok(json!({
        "format": BLOCK_RESULT_FORMAT,
        "format_version": crate::FORMAT_VERSION,
        "output": i.output,
        "report": report,
    }))
}
