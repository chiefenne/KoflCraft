//! The 3D trigger of the prism split (BS-20, step 2 part 2 of the 3D
//! insertion plan): find prism-shaped blocks in any block topology by
//! geometry, try the tri-center split there, and keep it only when a
//! measured trial is better (maintainer card 2026-11-27: geometry plus a
//! measured trial).
//!
//! The candidates are the flat block edges ([`flat_edges`], 150° or more
//! at both ends). A trial splits one with each row end asked for (`sheet`,
//! `t_junction`), builds the TFI grid and measures it in 3D
//! ([`algo::quality3d`]): finite-volume non-orthogonality and growth across
//! every inner face and the scaled Jacobian of every corner. It is kept when
//! the blocks it touched get a worst angle lower by [`MIN_GAIN_DEG`], their
//! worst scaled Jacobian reaches [`MIN_SCALED_JACOBIAN`] and is no worse
//! than before, and the whole grid gets no worse (worst angle, cells over
//! 70°, folds, growth) and grows by at most a share of its cells. There is
//! no 3D layout optimiser: the trial judges the TFI grid as it comes.

use std::collections::BTreeSet;

use mesh_algorithms::block as algo;
use mesh_algorithms::block::quality3d::{BlockStats, GridQuality};
use serde::Serialize;

use crate::EngineError;
use crate::block::TopologyFile;
use crate::block_edit::{Change, refused};
use crate::block_tet::{FLAT_TET_CORNER, TetSplit, split_tet, tet_blocks};
use crate::block_triangle::{FLAT_CORNER, RowStop, TriangleSplit, flat_edges, split_triangle_with};

/// The touched blocks' worst angle must drop by this much (degrees).
pub const MIN_GAIN_DEG: f64 = 2.0;
/// The touched blocks' worst scaled Jacobian after a split (HexIt's
/// threshold, from the 3D blocking survey of 2026-11-27).
pub const MIN_SCALED_JACOBIAN: f64 = 0.4;
/// The whole grid's worst angle may rise by this much at most.
const MAX_RISE_DEG: f64 = 0.05;
const GROWTH_LIMIT: f64 = 2.0;
const GROWTH_SLACK: f64 = 1.05;
/// Candidates tried per round, flattest first.
const TRIALS_PER_ROUND: usize = 4;

/// The settings of the 3D trigger.
#[derive(Clone, Debug, PartialEq)]
pub struct InsertOptions {
    /// Row ends to try for each candidate; the better trial is kept.
    pub stops: Vec<RowStop>,
    pub max_splits: usize,
    /// The cells a kept split may add, as a share of the grid's cells.
    pub max_added_share: f64,
}

impl Default for InsertOptions {
    fn default() -> Self {
        Self {
            stops: vec![RowStop::TJunction, RowStop::Sheet],
            max_splits: 4,
            max_added_share: 0.10,
        }
    }
}

/// The 3D quality in the report's words.
#[derive(Clone, Copy, Debug, Serialize, PartialEq)]
pub struct Quality {
    pub max_nonorthogonality_deg: f64,
    pub mean_nonorthogonality_deg: f64,
    pub cells_over_70_deg: usize,
    pub max_growth: f64,
    pub folded_corners: usize,
    pub min_scaled_jacobian: f64,
    pub cells: usize,
}

/// The numbers of a kept split: the touched blocks before and after, the
/// whole grid's worst angle and cells.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct Measured {
    pub touched_before_deg: f64,
    pub touched_after_deg: f64,
    pub touched_scaled_jacobian_before: f64,
    pub touched_scaled_jacobian_after: f64,
    pub max_before_deg: f64,
    pub max_after_deg: f64,
    pub cells_before: usize,
    pub cells_after: usize,
}

/// A prism split the trigger kept.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct KeptPiece {
    #[serde(flatten)]
    pub split: TriangleSplit,
    #[serde(flatten)]
    pub measured: Measured,
}

/// A tet split the trigger kept (BS-20).
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct KeptTet {
    #[serde(flatten)]
    pub split: TetSplit,
    #[serde(flatten)]
    pub measured: Measured,
}

/// A tet-shaped block the trigger tried and left, with why.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct SkippedTet {
    pub block: usize,
    pub corner: usize,
    pub corner_deg: f64,
    pub reason: String,
}

enum Kept {
    Prism(KeptPiece),
    Tet(KeptTet),
}

impl Kept {
    fn cells_after(&self) -> usize {
        match self {
            Kept::Prism(k) => k.measured.cells_after,
            Kept::Tet(k) => k.measured.cells_after,
        }
    }
}

/// A candidate the trigger tried and left, with why (each row end tried).
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct SkippedPiece {
    pub block: usize,
    pub vertex: usize,
    pub along: usize,
    pub corner_deg: f64,
    pub reason: String,
}

/// What [`insert_prism_splits`] did.
#[derive(Clone, Debug, Serialize)]
pub struct Insertion {
    #[serde(skip)]
    pub file: TopologyFile,
    pub kept: Vec<KeptPiece>,
    pub skipped: Vec<SkippedPiece>,
    pub kept_tets: Vec<KeptTet>,
    pub skipped_tets: Vec<SkippedTet>,
    pub before: Quality,
    pub after: Quality,
}

fn measure(file: &TopologyFile) -> Result<(Quality, Vec<BlockStats>), EngineError> {
    let grid = algo::build(&file.topology()?).map_err(|i| {
        refused(
            i.iter()
                .map(|x| x.message.clone())
                .collect::<Vec<_>>()
                .join("; "),
        )
    })?;
    let (q, per) = algo::quality3d::measure(&grid);
    Ok((quality(&q, grid.hexes.len()), per))
}

fn quality(q: &GridQuality, cells: usize) -> Quality {
    Quality {
        max_nonorthogonality_deg: q.max_deg,
        mean_nonorthogonality_deg: q.mean_deg,
        cells_over_70_deg: q.cells_over_70,
        max_growth: q.max_growth,
        folded_corners: q.folded_corners,
        min_scaled_jacobian: q.min_scaled_jacobian,
        cells,
    }
}

/// Worst angle and worst scaled Jacobian over some blocks.
fn worst(per: &[BlockStats], blocks: &[usize]) -> (f64, f64) {
    blocks.iter().fold((0.0, f64::INFINITY), |(d, s), &b| {
        (d.max(per[b].max_deg), s.min(per[b].min_scaled_jacobian))
    })
}

/// Add prism splits and tet splits where they measurably help, in any 3D
/// topology. Every round tries the flattest fresh candidates of both
/// kinds and keeps the best trial.
pub fn insert_prism_splits(
    file: &TopologyFile,
    options: &InsertOptions,
) -> Result<Insertion, EngineError> {
    if options.stops.is_empty() {
        return Err(refused("give at least one row end to try"));
    }
    let mut cur = file.clone();
    let (mut q, mut per) = measure(&cur)?;
    let before = q;
    let mut kept = Vec::new();
    let mut skipped = Vec::new();
    let mut kept_tets = Vec::new();
    let mut skipped_tets = Vec::new();
    let mut tried: BTreeSet<[usize; 8]> = BTreeSet::new();
    let mut tried_tets: BTreeSet<[usize; 8]> = BTreeSet::new();
    while kept.len() + kept_tets.len() < options.max_splits {
        let mut seen: BTreeSet<(usize, usize)> = BTreeSet::new();
        let fresh: Vec<_> = flat_edges(&cur, FLAT_CORNER)?
            .into_iter()
            .filter(|e| !tried.contains(&cur.blocks[e.block].vertices))
            .filter(|e| seen.insert((e.vertex.min(e.along), e.vertex.max(e.along))))
            .take(TRIALS_PER_ROUND)
            .collect();
        let tets: Vec<_> = tet_blocks(&cur, FLAT_TET_CORNER)?
            .into_iter()
            .filter(|t| !tried_tets.contains(&cur.blocks[t.block].vertices))
            .take(TRIALS_PER_ROUND)
            .collect();
        if fresh.is_empty() && tets.is_empty() {
            break;
        }
        let mut best: Option<(f64, TopologyFile, Kept)> = None;
        let mut offer = |gain: f64, f: TopologyFile, k: Kept| {
            let better = best.as_ref().is_none_or(|b| {
                gain > b.0 + 1e-9 || (gain > b.0 - 1e-9 && k.cells_after() < b.2.cells_after())
            });
            if better {
                best = Some((gain, f, k));
            }
        };
        for c in fresh {
            tried.insert(cur.blocks[c.block].vertices);
            let mut reasons = Vec::new();
            for &stop in &options.stops {
                let mut f = cur.clone();
                let tried = split_triangle_with(
                    &mut f,
                    c.block,
                    c.vertex,
                    Some(c.along),
                    None,
                    stop,
                    &mut Change::default(),
                )
                .map_err(|e| e.message)
                .and_then(|split| {
                    judge(&cur, &f, &q, &per, options)
                        .map(|(gain, measured)| (gain, Kept::Prism(KeptPiece { split, measured })))
                });
                match tried {
                    Ok((gain, k)) => offer(gain, f, k),
                    Err(r) => reasons.push(format!("{}: {r}", stop_name(stop))),
                }
            }
            if reasons.len() == options.stops.len() {
                skipped.push(SkippedPiece {
                    block: c.block,
                    vertex: c.vertex,
                    along: c.along,
                    corner_deg: c.corner_deg,
                    reason: reasons.join("; "),
                });
            }
        }
        for t in tets {
            tried_tets.insert(cur.blocks[t.block].vertices);
            let mut reasons = Vec::new();
            for &stop in &options.stops {
                let mut f = cur.clone();
                let tried = split_tet(
                    &mut f,
                    t.block,
                    Some(t.corner),
                    None,
                    stop,
                    crate::block_tet::Join::Auto,
                    &mut Change::default(),
                )
                .map_err(|e| e.message)
                .and_then(|split| {
                    judge(&cur, &f, &q, &per, options)
                        .map(|(gain, measured)| (gain, Kept::Tet(KeptTet { split, measured })))
                });
                match tried {
                    Ok((gain, k)) => offer(gain, f, k),
                    Err(r) => reasons.push(format!("{}: {r}", stop_name(stop))),
                }
            }
            if reasons.len() == options.stops.len() {
                skipped_tets.push(SkippedTet {
                    block: t.block,
                    corner: t.corner,
                    corner_deg: t.corner_deg,
                    reason: reasons.join("; "),
                });
            }
        }
        let Some((_, f, k)) = best else { break };
        (q, per) = measure(&f)?;
        cur = f;
        match k {
            Kept::Prism(k) => kept.push(k),
            Kept::Tet(k) => kept_tets.push(k),
        }
    }
    Ok(Insertion {
        file: cur,
        kept,
        skipped,
        kept_tets,
        skipped_tets,
        before,
        after: q,
    })
}

fn stop_name(s: RowStop) -> &'static str {
    match s {
        RowStop::Sheet => "sheet",
        RowStop::TJunction => "t_junction",
    }
}

/// Judge one trial `f` against `cur`: the gain in the touched blocks'
/// worst angle and the numbers, or why it is left.
fn judge(
    cur: &TopologyFile,
    f: &TopologyFile,
    q: &Quality,
    per: &[BlockStats],
    o: &InsertOptions,
) -> Result<(f64, Measured), String> {
    let (tq, tper) = measure(f).map_err(|e| e.message)?;
    let n = cur.blocks.len();
    let old: Vec<usize> = (0..n).filter(|&b| f.blocks[b] != cur.blocks[b]).collect();
    let new: Vec<usize> = (0..f.blocks.len())
        .filter(|&b| b >= n || f.blocks[b] != cur.blocks[b])
        .collect();
    let (d0, s0) = worst(per, &old);
    let (d1, s1) = worst(&tper, &new);
    let gain = d0 - d1;
    let added = tq.cells.saturating_sub(q.cells);
    let growth_cap = q.max_growth.max(GROWTH_LIMIT) * GROWTH_SLACK;
    if tq.folded_corners > 0 {
        return Err("the split folds cells".into());
    }
    if gain < MIN_GAIN_DEG {
        return Err(format!(
            "its blocks' worst angle goes from {d0:.1}° to {d1:.1}°, less than {MIN_GAIN_DEG}° better"
        ));
    }
    if s1 < MIN_SCALED_JACOBIAN || s1 < s0 {
        return Err(format!(
            "its blocks' worst scaled Jacobian goes from {s0:.2} to {s1:.2} (needs {MIN_SCALED_JACOBIAN} or more, and no worse)"
        ));
    }
    if tq.max_nonorthogonality_deg > q.max_nonorthogonality_deg + MAX_RISE_DEG {
        return Err(format!(
            "the grid's worst angle rises from {:.1}° to {:.1}°",
            q.max_nonorthogonality_deg, tq.max_nonorthogonality_deg
        ));
    }
    if tq.cells_over_70_deg > q.cells_over_70_deg {
        return Err(format!(
            "cells above 70° go from {} to {}",
            q.cells_over_70_deg, tq.cells_over_70_deg
        ));
    }
    if tq.max_growth > growth_cap {
        return Err(format!(
            "the growth reaches {:.2}, above {growth_cap:.2}",
            tq.max_growth
        ));
    }
    if added as f64 > o.max_added_share * q.cells as f64 {
        return Err(format!(
            "it adds {added} cells, more than {:.0} % of {}",
            100.0 * o.max_added_share,
            q.cells
        ));
    }
    Ok((
        gain,
        Measured {
            touched_before_deg: d0,
            touched_after_deg: d1,
            touched_scaled_jacobian_before: s0,
            touched_scaled_jacobian_after: s1,
            max_before_deg: q.max_nonorthogonality_deg,
            max_after_deg: tq.max_nonorthogonality_deg,
            cells_before: q.cells,
            cells_after: tq.cells,
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::block_triangle::tests::{turned, two_layers};

    /// Invariant (BS-20): on a topology that is not extruded along z the
    /// trigger finds the prism-shaped column by geometry, keeps the split
    /// that lowers its blocks' worst angle (48° to 30°) and lifts their
    /// worst scaled Jacobian above 0.4 (0.16 to 0.52), prefers the
    /// T-junction end with fewer cells, and the result checks and
    /// conforms. A trigger that keeps any split, ignores the Jacobian or
    /// only runs on extruded layouts breaks it; the 2D trigger tests never
    /// leave the plane.
    #[test]
    fn the_3d_trigger_keeps_a_prism_split_that_helps() {
        let f = turned(&two_layers());
        let o = InsertOptions {
            max_added_share: 2.0,
            ..Default::default()
        };
        let i = insert_prism_splits(&f, &o).unwrap();
        assert_eq!(i.kept.len(), 1, "{:?}", i.skipped);
        let k = &i.kept[0];
        assert_eq!(k.split.stop, RowStop::TJunction);
        let m = &k.measured;
        assert!(m.touched_before_deg - m.touched_after_deg > 15.0, "{k:?}");
        assert!(m.touched_scaled_jacobian_before < 0.2);
        assert!(m.touched_scaled_jacobian_after >= MIN_SCALED_JACOBIAN);
        assert!(i.after.max_nonorthogonality_deg < i.before.max_nonorthogonality_deg);
        assert_eq!(i.after.folded_corners, 0);
        crate::block::check(&i.file).unwrap();
        let (grid, _) = crate::block::grid(&i.file, &crate::block::VolumeAlgorithm::Tfi).unwrap();
        let names: Vec<&str> = i.file.patches.iter().map(|p| p.name.as_str()).collect();
        mesh_io::block::handoff(&grid, mesh_io::handoff::patches_from_names(&names)).unwrap();
    }

    /// Invariant: a split that adds more cells than the cap allows (10 %
    /// of the grid by default) is left, and the reason names the cap for
    /// each row end tried; `block.optimize_layout` off the z plane reports
    /// the 3D mode and leaves the layout as it was.
    #[test]
    fn the_3d_trigger_respects_the_cell_cap() {
        let f = turned(&two_layers());
        let (g, layout) = crate::block_triangle::optimize_layout(&f, true).unwrap();
        assert_eq!(layout["mode"], "3d");
        assert_eq!(layout["layout_optimised"], false);
        let skipped = layout["prism_splits"]["skipped"].as_array().unwrap();
        assert!(!skipped.is_empty());
        let reason = skipped[0]["reason"].as_str().unwrap();
        assert!(
            reason.contains("t_junction:") && reason.contains("sheet:"),
            "{reason}"
        );
        assert!(reason.contains("more than 10 %"), "{reason}");
        assert_eq!(g.blocks, f.blocks);
    }

    /// Invariant (BS-20): the trigger finds a tet-shaped block by geometry
    /// and keeps its split into four hexahedra when that lifts the worst
    /// scaled Jacobian (the flat corners' near 0) to 0.4 or more; with the
    /// default cap the same split is left because it adds more than 10 %
    /// of the cells.
    #[test]
    fn the_3d_trigger_keeps_a_tet_split_that_helps() {
        let f = crate::block_tet::tests::tet_on_slab();
        let o = InsertOptions {
            max_added_share: 10.0,
            ..Default::default()
        };
        let i = insert_prism_splits(&f, &o).unwrap();
        assert_eq!(i.kept_tets.len(), 1, "{:?} {:?}", i.skipped_tets, i.skipped);
        let k = &i.kept_tets[0];
        assert_eq!(k.split.neighbours, vec![1]);
        let m = &k.measured;
        assert!(m.touched_scaled_jacobian_before < 0.1, "{k:?}");
        assert!(
            m.touched_scaled_jacobian_after >= MIN_SCALED_JACOBIAN,
            "{k:?}"
        );
        assert!(
            m.touched_before_deg - m.touched_after_deg >= MIN_GAIN_DEG,
            "{k:?}"
        );
        assert_eq!(i.after.folded_corners, 0);
        crate::block::check(&i.file).unwrap();
        let i = insert_prism_splits(&f, &InsertOptions::default()).unwrap();
        assert!(i.kept_tets.is_empty());
        assert!(
            i.skipped_tets[0].reason.contains("more than 10 %"),
            "{:?}",
            i.skipped_tets
        );
    }
    /// Invariant (card A of step 3, automatic fallback): a tet whose
    /// neighbour cannot take a prism column is still split by the trigger,
    /// with a non-conformal interface, when the measured trial is better.
    #[test]
    fn the_3d_trigger_falls_back_to_a_non_conformal_tet_split() {
        let f = crate::block_tet::tests::tet_on_box();
        let o = InsertOptions {
            max_added_share: 10.0,
            ..Default::default()
        };
        let i = insert_prism_splits(&f, &o).unwrap();
        assert_eq!(i.kept_tets.len(), 1, "{:?} {:?}", i.skipped_tets, i.skipped);
        let k = &i.kept_tets[0];
        assert_eq!(k.split.join, "non_conformal");
        assert!(k.split.neighbours.is_empty());
        assert_eq!(i.file.interfaces.len(), 1);
        assert!(k.measured.touched_scaled_jacobian_after >= MIN_SCALED_JACOBIAN);
        crate::block::check(&i.file).unwrap();
    }
}
