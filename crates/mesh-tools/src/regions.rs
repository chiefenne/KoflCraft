//! `regions.propose`: the field-guided region proposal for block-dominant
//! meshing (which parts of a domain a block template can fill). Reads an
//! STL and its patch map, writes the report and optionally the field file
//! (VTU) and a picture (SVG). Nothing is meshed.

use mesh_engine::regions::{self, RegionOptions};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::Context;
use crate::diag::Diagnostic;
use crate::mesh::{Run, decode};

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ProposeInput {
    pub surface: String,
    pub patches: String,
    /// The report (JSON).
    pub output: String,
    /// The fields and regions on the scaffold (VTU).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fields: Option<String>,
    /// The picture (SVG).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub picture: Option<String>,
    /// Scaffold spacing (default: about 150 000 tetrahedra).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub size: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub heat_m: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub second_pass: Option<bool>,
}

fn write(ctx: &Context, path: &str, bytes: &[u8]) -> Result<Value, Diagnostic> {
    let (dir, name) = ctx.target(path)?;
    dir.write_new(&name, bytes)
        .map_err(|e| Diagnostic::new("io_error", e).with_path(path))?;
    Ok(json!({ "path": path, "blake3": mesh_engine::file_digest(bytes) }))
}

pub(crate) fn propose(input: Value, ctx: &Context) -> Run {
    let i: ProposeInput = decode(input)?;
    let stl = ctx.read(&i.surface)?;
    let map: Value = serde_json::from_slice(&ctx.read(&i.patches)?)
        .map_err(|e| Diagnostic::new("invalid_patch_map", e).with_path(&i.patches))?;
    let d = RegionOptions::default();
    let p = regions::propose(
        &stl,
        &map,
        RegionOptions {
            size: i.size,
            heat_m: i.heat_m.unwrap_or(d.heat_m),
            second_pass: i.second_pass.unwrap_or(d.second_pass),
            ..d
        },
    )?;
    let mut files = Vec::new();
    if let Some(f) = &i.fields {
        files.push(write(ctx, f, p.vtu.as_bytes())?);
    }
    if let Some(f) = &i.picture {
        files.push(write(ctx, f, p.svg.as_bytes())?);
    }
    let mut report = p.report;
    let mut bytes = serde_json::to_vec_pretty(&report)
        .map_err(|e| Diagnostic::new("serialization_failed", e))?;
    bytes.push(b'\n');
    files.insert(0, write(ctx, &i.output, &bytes)?);
    report["output"] = json!(i.output);
    report["files"] = json!(files);
    report["next"] = json!(
        "Open the picture (or the VTU: cell field kind 0 sweep, 1 junction, 2 obstacle, 3 box, 4 remainder) and check the proposal. Nothing is meshed yet."
    );
    Ok(report)
}
