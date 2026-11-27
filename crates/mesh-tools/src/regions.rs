//! `regions.propose`: the field-guided region proposal for block-dominant
//! meshing (which parts of a domain a block template can fill). Reads an
//! STL and its patch map, writes the report and optionally the field file
//! (VTU) and a picture (SVG). Nothing is meshed.

use mesh_engine::fill::{self, FillOptions};
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

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct FillInput {
    pub surface: String,
    pub patches: String,
    /// The new bundle directory.
    pub output: String,
    /// Scaffold spacing of the fields (default: about 150 000 tetrahedra).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub size: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cells_around: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub layers: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub growth: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub first_height: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub core_fraction: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub axial_size: Option<f64>,
    /// Mesh the region behind the cuts as polyhedra and join it to the
    /// sweeps at the cuts (default true).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub join: Option<bool>,
}

pub(crate) fn fill(input: Value, ctx: &Context) -> Run {
    let i: FillInput = decode(input)?;
    let stl = ctx.read(&i.surface)?;
    let map: Value = serde_json::from_slice(&ctx.read(&i.patches)?)
        .map_err(|e| Diagnostic::new("invalid_patch_map", e).with_path(&i.patches))?;
    let d = FillOptions::default();
    let opts = FillOptions {
        region: RegionOptions {
            size: i.size,
            ..d.region
        },
        cells_around: i.cells_around.unwrap_or(d.cells_around),
        layers: i.layers.unwrap_or(d.layers),
        growth: i.growth.unwrap_or(d.growth),
        first_height: i.first_height.or(d.first_height),
        core_fraction: i.core_fraction.unwrap_or(d.core_fraction),
        axial_size: i.axial_size.or(d.axial_size),
        join: i.join.unwrap_or(d.join),
    };
    let (dir, name) = ctx.target(&i.output)?;
    let mut report = fill::generate(&stl, &map, opts, &dir, &name)?;
    report["output"] = json!(i.output);
    let whole = report["shares"]["not_meshed"].as_f64().unwrap_or(1.0) == 0.0;
    report["next"] = json!(if whole {
        "Open foehncraft.vtu: the swept hexahedra and the polyhedra behind the cuts, one conformal mesh (cell field zone and the patches)."
    } else {
        "Open foehncraft.vtu: the swept hexahedra. The rest of the domain (shares.not_meshed) is not meshed; remainder says why."
    });
    Ok(report)
}
