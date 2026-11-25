//! Grading tools for the duct path (GRD-02, GRD-03, GRD-08): add, list and
//! remove datum objects and coordinate regions of a partition by name, so
//! an agent refines a region of a nonconvex mesh without editing JSON.
//! Each call reads a partition and writes a new one (the journal is the
//! history, AG-13).

use mesh_engine::duct::Partition;
use mesh_engine::sizing::{FrameConfig, SourceConfig};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::Context;
use crate::diag::Diagnostic;
use crate::mesh::{Run, decode, json_file};

fn write(ctx: &Context, output: &str, partition: &Partition) -> Run {
    let mut bytes = serde_json::to_vec_pretty(partition)
        .map_err(|e| Diagnostic::new("serialization_failed", e))?;
    bytes.push(b'\n');
    let (dir, name) = ctx.target(output)?;
    dir.write_new(&name, &bytes)
        .map_err(|e| Diagnostic::new("io_error", e).with_path(output))?;
    Ok(json!({
        "format": "koflcraft.grading-result",
        "format_version": crate::FORMAT_VERSION,
        "output": output,
        "files": [{ "path": output, "blake3": mesh_engine::file_digest(&bytes) }],
        "grading": listing(partition),
    }))
}

fn listing(p: &Partition) -> Value {
    match p.grading() {
        Some(g) => json!({
            "sources": g.sources,
            "frames": g.frames,
        }),
        None => json!({"sources": [], "frames": []}),
    }
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct AddInput {
    pub partition: String,
    pub output: String,
    pub source: SourceConfig,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub frame: Option<FrameConfig>,
}

pub(crate) fn add(input: Value, ctx: &Context) -> Run {
    let i: AddInput = decode(input)?;
    let mut p: Partition = json_file(ctx, &i.partition)?;
    p.add_grading(i.source, i.frame)?;
    write(ctx, &i.output, &p)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct RemoveInput {
    pub partition: String,
    pub output: String,
    pub name: String,
}

pub(crate) fn remove(input: Value, ctx: &Context) -> Run {
    let i: RemoveInput = decode(input)?;
    let mut p: Partition = json_file(ctx, &i.partition)?;
    p.remove_grading(&i.name)?;
    write(ctx, &i.output, &p)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ListInput {
    pub partition: String,
}

pub(crate) fn list(input: Value, ctx: &Context) -> Run {
    let i: ListInput = decode(input)?;
    let p: Partition = json_file(ctx, &i.partition)?;
    Ok(listing(&p))
}
