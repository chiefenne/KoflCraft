//! `mesh.estimate` and `--check` (AG-15): what a build would make, before
//! it runs. The cell count is exact: each generator counts its cells from
//! the inputs (sites, the channel grid, the swept 2D section, the duct
//! surface or scaffold) without building the volume mesh. The one exception
//! is labelled in `cells_basis`: a polyhedral duct splits or adds a few
//! cells at reentrant wall edges, which only the built wall layers decide,
//! and the wall columns under a Voronoi core follow the built core, so their
//! count is approximate. Memory and
//! output bytes are rough per-cell figures, labelled as such.

use std::path::Path;

use mesh_engine::duct::{self, Partition, PatchMap};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::diag::Diagnostic;
use crate::mesh::{
    ChannelInput, DuctInput, GenerateInput, GenerateStlInput, SweepInput, decode, json_file,
};
use crate::{Action, Context, FORMAT_VERSION, find};

pub const ESTIMATE_FORMAT: &str = "koflcraft.estimate";

/// Bytes per cell, measured on the fixtures on 2026-11-22 and rounded up:
/// written files 0.3 to 1.9 kB per cell (box artifacts are the largest),
/// working memory a few kB per cell.
const OUTPUT_BYTES_PER_CELL: u64 = 2_048;
const MEMORY_BYTES_PER_CELL: u64 = 4_096;

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct EstimateInput {
    /// The id of a build tool.
    pub tool: String,
    /// That tool's input.
    pub input: Value,
}

/// `cells_basis` of a count that is exact.
const EXACT: &str = "exact";
/// `cells_basis` of a polyhedral duct with a volume scaffold.
const SPLIT_BASIS: &str = "exact before cells are split at reentrant wall edges or where they \
    break the solver rules, which needs the built layers (median dual: 16 on the elbow fixture)";

/// `cells_basis` of a Voronoi core with wall columns.
const COLUMN_BASIS: &str = "approximate: the wall-layer columns follow the built Voronoi core, \
    counted as two cells per wall triangle and layer (built: 0.76 on the elbow fixture, 2.08 on \
    Manifold); cells are also added at reentrant wall edges";

/// The cell count of one build call and its basis.
pub(crate) fn cells(
    id: &str,
    input: Value,
    ctx: &Context,
) -> Result<(usize, &'static str), Diagnostic> {
    let tool = find(id).ok_or_else(|| Diagnostic::new("unknown_tool", format!("no tool {id}")))?;
    if tool.action != Action::StartBuild {
        return Err(Diagnostic::new(
            "invalid_arguments",
            format!(
                "{} is not a build tool; only start-build tools have an estimate",
                tool.id
            ),
        ));
    }
    let mut basis = EXACT;
    let n = match tool.id {
        "mesh.generate" => decode::<GenerateInput>(input)?.config.sites.len(),
        "mesh.generate_stl" => {
            let i: GenerateStlInput = decode(input)?;
            let facets = ctx.limits.as_ref().map(|l| l.max_surface_facets);
            mesh_engine::stl_site_count(&ctx.read(&i.surface)?, &i.config, facets)?
        }
        "mesh.generate_channel" => {
            mesh_engine::channel::cell_count(&decode::<ChannelInput>(input)?.config)?
        }
        "mesh.generate_sweep" => {
            let budget = ctx
                .limits
                .as_ref()
                .map_or(usize::MAX / 1024, |l| l.max_sweep_cells);
            mesh_engine::sweep::cell_count(&decode::<SweepInput>(input)?.config, budget)?
        }
        "mesh.generate_duct" => {
            let i: DuctInput = decode(input)?;
            let patches: PatchMap = json_file(ctx, &i.patches)?;
            let partition: Partition = json_file(ctx, &i.partition)?;
            let scaffold = partition
                .scaffold_file(i.mode)
                .map(|file| {
                    let path = Path::new(&i.partition).with_file_name(file);
                    ctx.read(&path.to_string_lossy())
                })
                .transpose()?;
            if partition.wall_columns(i.mode) {
                basis = COLUMN_BASIS;
            } else if scaffold.is_some() {
                basis = SPLIT_BASIS;
            }
            duct::cell_count(
                &ctx.read(&i.surface)?,
                &patches,
                &partition,
                i.mode,
                scaffold.as_deref(),
            )?
        }
        "block.grid" => {
            let i: crate::block::GridInput = decode(input)?;
            mesh_engine::block::cell_count(&crate::block::topology(ctx, &i.topology)?)?
        }
        "journal.replay" => {
            let (n, b) = replay_cells(input, ctx)?;
            basis = b;
            n
        }
        other => {
            return Err(Diagnostic::new(
                "invalid_arguments",
                format!("no estimate for {other}"),
            ));
        }
    };
    Ok((n, basis))
}

/// The cells of every build call a journal would replay.
fn replay_cells(input: Value, ctx: &Context) -> Result<(usize, &'static str), Diagnostic> {
    let i: crate::ReplayInput = decode(input)?;
    let rows = crate::journal::calls(ctx, &i.script)?;
    let mut total = 0usize;
    let mut basis = EXACT;
    for row in rows {
        if crate::journal::skipped(&row) {
            continue;
        }
        if find(&row.tool).is_some_and(|t| t.action == Action::StartBuild) {
            let (n, b) = cells(&row.tool, row.input, ctx)?;
            total = total.saturating_add(n);
            if b != EXACT {
                basis = b;
            }
        }
    }
    Ok((total, basis))
}

/// `mesh.estimate`.
pub(crate) fn estimate(input: Value, ctx: &Context) -> Result<Value, Diagnostic> {
    let i: EstimateInput = decode(input)?;
    let tool = find(&i.tool)
        .ok_or_else(|| Diagnostic::new("unknown_tool", format!("no tool {}", i.tool)))?;
    let files = tool.effects(&i.input);
    let (n, basis) = cells(tool.id, i.input, ctx)?;
    let n64 = n as u64;
    Ok(json!({
        "format": ESTIMATE_FORMAT,
        "format_version": FORMAT_VERSION,
        "tool": tool.id,
        "cells": n,
        "cells_basis": basis,
        "memory_bytes": n64.saturating_mul(MEMORY_BYTES_PER_CELL),
        "output_bytes": n64.saturating_mul(OUTPUT_BYTES_PER_CELL),
        "bytes_basis": "rough: per-cell figures measured on the fixtures",
        "files": files,
    }))
}
