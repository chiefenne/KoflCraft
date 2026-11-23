//! The boundary patch tools (`patch.*`, BC-01, BC-02): each reads an STL
//! and a patch map, applies one operation and writes a new patch map, so an
//! agent defines the boundary of a FoehnCraft case step by step and the
//! journal replays every step (AG-13). `patch.list` only reads.

use mesh_engine::patches::{self, Contain, PatchMap, PatchRole, PatchSurface, Selection, Target};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use sha2::Digest;

use crate::Context;
use crate::diag::Diagnostic;
use crate::mesh::{Run, decode, json_file};

pub const PATCH_RESULT_FORMAT: &str = "koflcraft.patch-result";

fn surface(ctx: &Context, path: &str) -> Result<PatchSurface, Diagnostic> {
    let bytes = ctx.read(path)?;
    // the file size limit bounds the import; the facet limit follows
    let surface = PatchSurface::import(&bytes, None)?;
    ctx.limit("patch facets", surface.facet_count(), |l| {
        l.max_patch_facets
    })?;
    Ok(surface)
}

/// Write the new map and report it with its evidence.
fn write(
    ctx: &Context,
    surface: &PatchSurface,
    output: &str,
    (map, change): (PatchMap, patches::Change),
) -> Run {
    let summary = surface.summary(&map)?;
    let mut bytes =
        serde_json::to_vec_pretty(&map).map_err(|e| Diagnostic::new("serialization_failed", e))?;
    bytes.push(b'\n');
    let (dir, name) = ctx.target(output)?;
    dir.write_new(&name, &bytes)
        .map_err(|e| Diagnostic::new("io_error", e).with_path(output))?;
    Ok(json!({
        "format": PATCH_RESULT_FORMAT,
        "format_version": crate::FORMAT_VERSION,
        "output": output,
        "files": [{ "path": output, "sha256": crate::journal::hex(&sha2::Sha256::digest(&bytes)) }],
        "change": change,
        "summary": summary,
    }))
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct InitInput {
    pub surface: String,
    pub output: String,
}

pub(crate) fn init(input: Value, ctx: &Context) -> Run {
    let i: InitInput = decode(input)?;
    let s = surface(ctx, &i.surface)?;
    let map = s.initial_map();
    let created = map.patches.iter().map(|p| p.name.clone()).collect();
    write(
        ctx,
        &s,
        &i.output,
        (
            map,
            patches::Change {
                facets_changed: s.facet_count(),
                created,
                removed: Vec::new(),
            },
        ),
    )
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ListInput {
    pub surface: String,
    pub patches: String,
}

pub(crate) fn list(input: Value, ctx: &Context) -> Run {
    let i: ListInput = decode(input)?;
    let s = surface(ctx, &i.surface)?;
    let map: PatchMap = json_file(ctx, &i.patches)?;
    Ok(json!(s.summary(&map)?))
}

/// The fields every selection tool shares.
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SelectInput<S> {
    pub surface: String,
    pub patches: String,
    pub output: String,
    /// The patch the selected facets move to; created when new.
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub role: Option<PatchRole>,
    /// Move only facets of these patches.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub within: Option<Vec<String>>,
    #[serde(flatten)]
    pub selection: S,
}

#[derive(Deserialize, Serialize)]
pub struct PlaneArgs {
    pub point: [f64; 3],
    pub normal: [f64; 3],
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tolerance: Option<f64>,
}

#[derive(Deserialize, Serialize)]
pub struct BoxArgs {
    pub min: [f64; 3],
    pub max: [f64; 3],
    #[serde(default)]
    pub contain: Contain,
}

#[derive(Deserialize, Serialize)]
pub struct NormalArgs {
    pub direction: [f64; 3],
    pub max_angle_deg: f64,
}

#[derive(Deserialize, Serialize)]
pub struct GrowArgs {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seed_facet: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seed_point: Option<[f64; 3]>,
    pub feature_angle_deg: f64,
}

/// `serde(flatten)` does not refuse unknown fields, so compare the keys
/// with the schema's.
fn closed(input: &Value, schema: Value) -> Result<(), Diagnostic> {
    let known = &schema["properties"];
    if let Some(object) = input.as_object()
        && let Some(extra) = object.keys().find(|k| known.get(k.as_str()).is_none())
    {
        return Err(Diagnostic::new(
            "invalid_arguments",
            format!("unknown field `{extra}`"),
        ));
    }
    Ok(())
}

fn select<S: serde::de::DeserializeOwned>(
    input: Value,
    ctx: &Context,
    schema: fn() -> Value,
    selection: impl FnOnce(S) -> Selection,
) -> Run {
    closed(&input, schema())?;
    let i: SelectInput<S> = decode(input)?;
    let s = surface(ctx, &i.surface)?;
    let map: PatchMap = json_file(ctx, &i.patches)?;
    let target = Target {
        name: i.name,
        role: i.role,
        within: i.within,
    };
    let result = s.select(&map, &selection(i.selection), &target)?;
    write(ctx, &s, &i.output, result)
}

pub(crate) fn select_plane(input: Value, ctx: &Context) -> Run {
    select(input, ctx, crate::schemas::patch_plane, |a: PlaneArgs| {
        Selection::Plane {
            point: a.point,
            normal: a.normal,
            tolerance: a.tolerance,
        }
    })
}

pub(crate) fn select_box(input: Value, ctx: &Context) -> Run {
    select(input, ctx, crate::schemas::patch_box, |a: BoxArgs| {
        Selection::Box {
            min: a.min,
            max: a.max,
            contain: a.contain,
        }
    })
}

pub(crate) fn select_normal(input: Value, ctx: &Context) -> Run {
    select(input, ctx, crate::schemas::patch_normal, |a: NormalArgs| {
        Selection::Normal {
            direction: a.direction,
            max_angle_deg: a.max_angle_deg,
        }
    })
}

pub(crate) fn grow(input: Value, ctx: &Context) -> Run {
    select(input, ctx, crate::schemas::patch_grow, |a: GrowArgs| {
        Selection::Grow {
            seed_facet: a.seed_facet,
            seed_point: a.seed_point,
            feature_angle_deg: a.feature_angle_deg,
        }
    })
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SplitInput {
    pub surface: String,
    pub patches: String,
    pub output: String,
    pub patch: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub feature_angle_deg: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub names: Option<Vec<String>>,
}

pub(crate) fn split(input: Value, ctx: &Context) -> Run {
    let i: SplitInput = decode(input)?;
    let s = surface(ctx, &i.surface)?;
    let map: PatchMap = json_file(ctx, &i.patches)?;
    let result = s.split(&map, &i.patch, i.feature_angle_deg, i.names.as_deref())?;
    write(ctx, &s, &i.output, result)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct MergeInput {
    pub surface: String,
    pub patches: String,
    pub output: String,
    /// The patches to join.
    pub merge: Vec<String>,
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub role: Option<PatchRole>,
}

pub(crate) fn merge(input: Value, ctx: &Context) -> Run {
    let i: MergeInput = decode(input)?;
    let s = surface(ctx, &i.surface)?;
    let map: PatchMap = json_file(ctx, &i.patches)?;
    s.check(&map)?;
    let result = patches::merge(&map, &i.merge, &i.name, i.role)?;
    write(ctx, &s, &i.output, result)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct RenameInput {
    pub surface: String,
    pub patches: String,
    pub output: String,
    pub from: String,
    pub to: String,
}

pub(crate) fn rename(input: Value, ctx: &Context) -> Run {
    let i: RenameInput = decode(input)?;
    let s = surface(ctx, &i.surface)?;
    let map: PatchMap = json_file(ctx, &i.patches)?;
    s.check(&map)?;
    let result = patches::rename(&map, &i.from, &i.to)?;
    write(ctx, &s, &i.output, result)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SetRoleInput {
    pub surface: String,
    pub patches: String,
    pub output: String,
    pub name: String,
    pub role: PatchRole,
}

pub(crate) fn set_role(input: Value, ctx: &Context) -> Run {
    let i: SetRoleInput = decode(input)?;
    let s = surface(ctx, &i.surface)?;
    let map: PatchMap = json_file(ctx, &i.patches)?;
    s.check(&map)?;
    let result = patches::set_role(&map, &i.name, i.role)?;
    write(ctx, &s, &i.output, result)
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct DetectPortsInput {
    pub surface: String,
    pub output: String,
    /// Candidate ids for the inlet and the outlet (default 0 and 1).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inlet: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub outlet: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub planar_angle_deg: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub candidates: Option<usize>,
}

/// Find the planar end caps of a bare closed STL and write a map with
/// `walls`, `inlet` and `outlet` (BC-03, ports).
pub(crate) fn detect_ports(input: Value, ctx: &Context) -> Run {
    let i: DetectPortsInput = decode(input)?;
    let s = surface(ctx, &i.surface)?;
    let defaults = mesh_engine::ports::PortOptions::default();
    let detection = s.detect_ports(mesh_engine::ports::PortOptions {
        planar_angle_deg: i.planar_angle_deg.unwrap_or(defaults.planar_angle_deg),
        candidates: i.candidates.unwrap_or(defaults.candidates),
    })?;
    let (inlet, outlet) = (i.inlet.unwrap_or(0), i.outlet.unwrap_or(1));
    let map = s.port_map(&detection, inlet, outlet)?;
    let ports = map.facet_patch.iter().filter(|&&p| p != 0).count();
    let change = patches::Change {
        facets_changed: s.facet_count(),
        created: vec!["walls".into(), "inlet".into(), "outlet".into()],
        removed: Vec::new(),
    };
    let mut result = write(ctx, &s, &i.output, (map, change))?;
    let strength = |c: &mesh_engine::ports::PortCandidate| c.score;
    let next = detection.candidates.get(2).map(strength).unwrap_or(0.0);
    result["detection"] = json!({
        "inlet": inlet,
        "outlet": outlet,
        "port_facets": ports,
        "planar_regions": detection.planar_regions,
        "disc_regions": detection.disc_regions,
        // how clearly the two picked caps stand out from the next candidate
        "margin": strength(&detection.candidates[inlet]).min(strength(&detection.candidates[outlet])) - next,
        "candidates": detection.candidates,
        "note": "Inlet and outlet are the two most pipe-end-like caps; which is which is a guess. Swap or pick others with `inlet` and `outlet` (candidate ids).",
    });
    Ok(result)
}
