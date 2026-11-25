//! Shell tools for the duct path (SHL-01): add, list and remove the
//! zero-thickness shells (baffles, plates, fins) of a partition by name, so
//! an agent places a shell without editing JSON. Each call reads a
//! partition and writes a new one (the journal is the history, AG-13).

use mesh_engine::duct::{Partition, ShellConfig};
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
    let listing = partition.shell_listing(None)?;
    Ok(json!({
        "format": "koflcraft.shell-result",
        "format_version": crate::FORMAT_VERSION,
        "output": output,
        "files": [{ "path": output, "blake3": mesh_engine::file_digest(&bytes) }],
        "shells": listing["shells"],
    }))
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct AddInput {
    pub partition: String,
    pub output: String,
    pub shell: ShellConfig,
}

pub(crate) fn add(input: Value, ctx: &Context) -> Run {
    let i: AddInput = decode(input)?;
    let mut p: Partition = json_file(ctx, &i.partition)?;
    p.add_shell(i.shell)?;
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
    p.remove_shell(&i.name)?;
    write(ctx, &i.output, &p)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ListInput {
    pub partition: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub surface: Option<String>,
}

pub(crate) fn list(input: Value, ctx: &Context) -> Run {
    let i: ListInput = decode(input)?;
    let p: Partition = json_file(ctx, &i.partition)?;
    let stl = i.surface.as_deref().map(|s| ctx.read(s)).transpose()?;
    Ok(p.shell_listing(stl.as_deref())?)
}
