//! The block-structured variant behind the engine (OBJ-21, BS-01..BS-14):
//! the topology file `koflcraft.block-topology` (BS-02), the templates
//! (BS-03), the check (BS-05) and the build of a validated hexahedral mesh
//! into a new bundle (BS-13, BS-14).

use std::ffi::OsStr;

use mesh_algorithms::block::{self as algo, BlockGrid, Issue, templates};
use mesh_core::block::{
    BlockSpec, BlockTopology, Curve, Distribution, EdgeEnd, EdgeSpec, PatchSpec,
};
use mesh_core::geometry::Vec3;
use mesh_io::handoff::{PatchRole, patches_from_names};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::{EngineError, ErrorCode, MeshDirectory};

pub const TOPOLOGY_FORMAT: &str = "koflcraft.block-topology";
pub const TOPOLOGY_FORMAT_VERSION: &str = "1.0";
pub const ALGORITHM: &str = "block_tfi_v1";
/// Cells and points of a sum of cell volumes and of the volume the boundary
/// encloses may differ by this fraction at most (rounding only).
pub const CONSERVATION_TOLERANCE: f64 = 1e-10;

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum CurveFile {
    Line,
    Arc { through: [f64; 3] },
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(tag = "law", rename_all = "snake_case", deny_unknown_fields)]
pub enum DistributionFile {
    Uniform,
    /// Exactly one of `first` (at the edge's first vertex) and `last`.
    Geometric {
        #[serde(default, skip_serializing_if = "Option::is_none")]
        first: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        last: Option<f64>,
    },
    Tanh {
        first: f64,
        last: f64,
    },
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct EdgeFile {
    pub vertices: [usize; 2],
    #[serde(default = "line")]
    pub curve: CurveFile,
    #[serde(default = "uniform")]
    pub distribution: DistributionFile,
}

fn line() -> CurveFile {
    CurveFile::Line
}
fn uniform() -> DistributionFile {
    DistributionFile::Uniform
}
fn fluid() -> String {
    "fluid".into()
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct BlockFile {
    pub vertices: [usize; 8],
    pub cells: [usize; 3],
    #[serde(default = "fluid")]
    pub domain: String,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct PatchFile {
    pub name: String,
    /// The role; from the name when absent (BC-04).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub role: Option<PatchRole>,
    pub faces: Vec<[usize; 4]>,
}

/// Where a topology came from and what its shape's analytic volume is.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Reference {
    pub template: String,
    pub parameters: Value,
    pub analytic_volume: f64,
}

/// The file `koflcraft.block-topology` (BS-02).
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct TopologyFile {
    pub format: String,
    pub format_version: String,
    pub vertices: Vec<[f64; 3]>,
    #[serde(default)]
    pub edges: Vec<EdgeFile>,
    pub blocks: Vec<BlockFile>,
    pub patches: Vec<PatchFile>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference: Option<Reference>,
}

/// A template and its parameters (BS-03).
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(tag = "template", rename_all = "snake_case", deny_unknown_fields)]
pub enum TemplateConfig {
    /// A rectangular channel along +x, one block.
    Channel {
        length: f64,
        width: f64,
        height: f64,
        /// Cells along x, y, z.
        cells: [usize; 3],
        /// First cell height at the four walls (two-sided tanh).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        wall_spacing: Option<f64>,
    },
    /// A circular pipe along +x: H-core and four O-grid blocks.
    PipeOgrid {
        radius: f64,
        length: f64,
        /// Radius of the core corners over the pipe radius (default 0.55).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        core_fraction: Option<f64>,
        /// Cells along a core side (a quarter of the circumference).
        core_cells: usize,
        radial_cells: usize,
        axial_cells: usize,
        /// First cell height at the wall (geometric growth inward).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        wall_spacing: Option<f64>,
    },
}

fn invalid(message: impl std::fmt::Display) -> EngineError {
    EngineError::new(ErrorCode::InvalidBlockTopology, message)
}

fn issues_error(issues: Vec<Issue>) -> EngineError {
    let text: Vec<String> = issues
        .iter()
        .map(|i| format!("[{}] {}", i.kind, i.message))
        .collect();
    invalid(text.join("; "))
}

fn v3(p: [f64; 3]) -> Vec3 {
    Vec3::new(p[0], p[1], p[2])
}

fn a3(p: Vec3) -> [f64; 3] {
    [p.x, p.y, p.z]
}

impl TopologyFile {
    /// The in-memory topology; format, version and edge laws are checked.
    pub fn topology(&self) -> Result<BlockTopology, EngineError> {
        if self.format != TOPOLOGY_FORMAT {
            return Err(invalid(format!("the file is not a {TOPOLOGY_FORMAT}")));
        }
        if self.format_version != TOPOLOGY_FORMAT_VERSION {
            return Err(EngineError::new(
                ErrorCode::UnsupportedSchema,
                format!(
                    "{TOPOLOGY_FORMAT} version {} is not supported (expected {TOPOLOGY_FORMAT_VERSION})",
                    self.format_version
                ),
            ));
        }
        let edges = self
            .edges
            .iter()
            .map(|e| {
                let distribution = match e.distribution {
                    DistributionFile::Uniform => Distribution::Uniform,
                    DistributionFile::Geometric {
                        first: Some(s),
                        last: None,
                    } => Distribution::Geometric {
                        spacing: s,
                        at: EdgeEnd::Start,
                    },
                    DistributionFile::Geometric {
                        first: None,
                        last: Some(s),
                    } => Distribution::Geometric {
                        spacing: s,
                        at: EdgeEnd::End,
                    },
                    DistributionFile::Geometric { .. } => {
                        return Err(invalid(format!(
                            "edge {:?}: a geometric law needs exactly one of first and last",
                            e.vertices
                        )));
                    }
                    DistributionFile::Tanh { first, last } => Distribution::Tanh { first, last },
                };
                Ok(EdgeSpec {
                    vertices: e.vertices,
                    curve: match e.curve {
                        CurveFile::Line => Curve::Line,
                        CurveFile::Arc { through } => Curve::Arc {
                            through: v3(through),
                        },
                    },
                    distribution,
                })
            })
            .collect::<Result<_, _>>()?;
        Ok(BlockTopology {
            vertices: self.vertices.iter().copied().map(v3).collect(),
            edges,
            blocks: self
                .blocks
                .iter()
                .map(|b| BlockSpec {
                    vertices: b.vertices,
                    cells: b.cells,
                    domain: b.domain.clone(),
                })
                .collect(),
            patches: self
                .patches
                .iter()
                .map(|p| PatchSpec {
                    name: p.name.clone(),
                    faces: p.faces.clone(),
                })
                .collect(),
        })
    }

    fn from_topology(t: &BlockTopology, reference: Option<Reference>) -> Self {
        Self {
            format: TOPOLOGY_FORMAT.into(),
            format_version: TOPOLOGY_FORMAT_VERSION.into(),
            vertices: t.vertices.iter().copied().map(a3).collect(),
            edges: t
                .edges
                .iter()
                .map(|e| EdgeFile {
                    vertices: e.vertices,
                    curve: match e.curve {
                        Curve::Line => CurveFile::Line,
                        Curve::Arc { through } => CurveFile::Arc {
                            through: a3(through),
                        },
                    },
                    distribution: match e.distribution {
                        Distribution::Uniform => DistributionFile::Uniform,
                        Distribution::Geometric { spacing, at } => DistributionFile::Geometric {
                            first: (at == EdgeEnd::Start).then_some(spacing),
                            last: (at == EdgeEnd::End).then_some(spacing),
                        },
                        Distribution::Tanh { first, last } => {
                            DistributionFile::Tanh { first, last }
                        }
                    },
                })
                .collect(),
            blocks: t
                .blocks
                .iter()
                .map(|b| BlockFile {
                    vertices: b.vertices,
                    cells: b.cells,
                    domain: b.domain.clone(),
                })
                .collect(),
            patches: t
                .patches
                .iter()
                .map(|p| PatchFile {
                    name: p.name.clone(),
                    role: Some(PatchRole::from_name(&p.name)),
                    faces: p.faces.clone(),
                })
                .collect(),
            reference,
        }
    }
}

/// Write the topology of a template (`block.template`).
pub fn template(config: &TemplateConfig) -> Result<TopologyFile, EngineError> {
    let (t, name) = match *config {
        TemplateConfig::Channel {
            length,
            width,
            height,
            cells,
            wall_spacing,
        } => (
            templates::channel([length, width, height], cells, wall_spacing),
            "channel",
        ),
        TemplateConfig::PipeOgrid {
            radius,
            length,
            core_fraction,
            core_cells,
            radial_cells,
            axial_cells,
            wall_spacing,
        } => (
            templates::pipe_ogrid(
                radius,
                length,
                core_fraction.unwrap_or(0.55),
                [core_cells, radial_cells, axial_cells],
                wall_spacing,
            ),
            "pipe_ogrid",
        ),
    };
    let t = t.map_err(invalid)?;
    let parameters = serde_json::to_value(config)
        .map_err(|e| EngineError::new(ErrorCode::SerializationFailed, e))?;
    let file = TopologyFile::from_topology(
        &t.topology,
        Some(Reference {
            template: name.into(),
            parameters,
            analytic_volume: t.analytic_volume,
        }),
    );
    // a template's topology always passes its own check
    algo::check(&t.topology).map_err(issues_error)?;
    Ok(file)
}

/// The check of a topology (`block.check`): counts, chains, faces, or every
/// issue found.
pub fn check(file: &TopologyFile) -> Result<Value, EngineError> {
    let t = file.topology()?;
    let c = algo::check(&t).map_err(issues_error)?;
    Ok(json!({
        "blocks": t.blocks.len(),
        "cells": c.cells,
        "points": c.points,
        "shared_block_faces": c.block_faces_shared,
        "boundary_block_faces": c.boundary_faces,
        "chains": c.chains.iter().map(|ch| json!({
            "cells": ch.cells, "edges": ch.edges, "blocks": ch.blocks
        })).collect::<Vec<_>>(),
    }))
}

/// The cell count of a topology, for the estimate (AG-15).
pub fn cell_count(file: &TopologyFile) -> Result<usize, EngineError> {
    let t = file.topology()?;
    algo::check(&t).map(|c| c.cells).map_err(issues_error)
}

fn validation(message: impl std::fmt::Display) -> EngineError {
    EngineError::new(ErrorCode::ValidationFailed, message)
}

/// Build, validate and write a block-structured mesh into the new
/// directory `name` (`block.grid`): `volume.vtu`, `mesh.msh` (MSH 4.1),
/// `foehncraft.vtu`, `report.json` and `mesh-evidence.json`.
pub fn generate(
    file: &TopologyFile,
    dir: &MeshDirectory,
    name: &OsStr,
) -> Result<Value, EngineError> {
    let t = file.topology()?;
    let grid = algo::build(&t).map_err(issues_error)?;
    let (report, handoff, structure) = validate(file, &grid)?;
    let json_report = serde_json::to_vec_pretty(&report)
        .map_err(|e| EngineError::new(ErrorCode::SerializationFailed, e))?;
    let fail = crate::handoff::stream_error;
    let mut bundle = dir
        .bundle(name)
        .map_err(|e| EngineError::new(ErrorCode::IoError, e))?;
    bundle
        .stream("volume.vtu", |w| mesh_io::block::write_volume(&grid, w))
        .map_err(fail)?;
    bundle
        .stream("mesh.msh", |w| {
            mesh_io::block::write_gmsh(&grid, &handoff.patches, mesh_io::GmshVersion::V4_1, w)
        })
        .map_err(fail)?;
    bundle
        .stream(crate::handoff::FOEHNCRAFT_VTU, |w| {
            crate::handoff::HandoffFormat::Vtu.write(
                &handoff,
                &[("koflcraft_algorithm", ALGORITHM.to_string())],
                w,
            )
        })
        .map_err(fail)?;
    bundle
        .stream("report.json", |w| w.write_all(&json_report))
        .map_err(fail)?;
    let evidence =
        crate::handoff::evidence_with(&handoff, ALGORITHM, bundle.files(), Some(structure))?;
    bundle
        .stream(crate::handoff::MESH_EVIDENCE, |w| w.write_all(&evidence))
        .map_err(fail)?;
    bundle.finish();
    Ok(report)
}

/// The checks of BS-13 and the report; the hand-off mesh is the conforming
/// face-based mesh every writer of FoehnCraft files uses.
fn validate(
    file: &TopologyFile,
    grid: &BlockGrid,
) -> Result<(Value, mesh_io::handoff::HandoffMesh, Value), EngineError> {
    let quality = grid.quality();
    if let Some(q) = quality
        .iter()
        .find(|q| q.min_jacobian.is_nan() || q.min_jacobian <= 0.0)
    {
        return Err(validation(format!(
            "block {} has a hexahedron with the corner Jacobian {:e}; TFI folds this block (elliptic smoothing comes in MB1b)",
            q.block, q.min_jacobian
        )));
    }
    let volumes = grid.volumes();
    let total: f64 = volumes.iter().sum();
    let enclosed = grid.boundary_volume();
    let conservation = (total / enclosed - 1.0).abs();
    if conservation.is_nan() || conservation > CONSERVATION_TOLERANCE {
        return Err(validation(format!(
            "cell volumes sum to {total:e} but the boundary encloses {enclosed:e}"
        )));
    }
    let names: Vec<&str> = file.patches.iter().map(|p| p.name.as_str()).collect();
    let mut patches = patches_from_names(&names);
    for (patch, f) in patches.iter_mut().zip(&file.patches) {
        if let Some(role) = f.role {
            patch.role = role;
        }
    }
    let handoff = mesh_io::block::handoff(grid, patches)
        .map_err(|e| EngineError::new(ErrorCode::InvalidConnectivity, e))?;
    let summary = handoff.summary();
    if summary.inverted_faces > 0 {
        return Err(validation(format!(
            "{} interior faces have S·d ≤ 0",
            summary.inverted_faces
        )));
    }
    let blocks: Vec<Value> = grid
        .blocks
        .iter()
        .zip(&quality)
        .map(|(r, q)| {
            json!({
                "block": r.block, "domain": r.domain, "cells": r.cells,
                "first_cell": r.first_cell, "cell_count": r.cells.iter().product::<usize>(),
                "min_jacobian": q.min_jacobian, "min_scaled_jacobian": q.min_scaled_jacobian,
                "max_aspect_ratio": q.max_aspect_ratio,
            })
        })
        .collect();
    let structure = json!({
        "cell_order": "block by block; inside a block i fastest, then j, then k",
        "blocks": grid.blocks.iter().map(|r| json!({
            "block": r.block, "domain": r.domain, "cells": r.cells, "first_cell": r.first_cell
        })).collect::<Vec<_>>(),
    });
    let min_volume = volumes.iter().copied().fold(f64::INFINITY, f64::min);
    let reference = file.reference.as_ref().map(|r| {
        json!({
            "template": r.template,
            "analytic_volume": r.analytic_volume,
            "relative_deviation": total / r.analytic_volume - 1.0,
        })
    });
    let report = json!({
        "format": "koflcraft.block-report",
        "format_version": "1.0",
        "algorithm": ALGORITHM,
        "cells": grid.hexes.len(),
        "points": grid.points.len(),
        "volume": total,
        "boundary_volume": enclosed,
        "relative_conservation_error": conservation,
        "min_cell_volume": min_volume,
        "min_scaled_jacobian": quality.iter().map(|q| q.min_scaled_jacobian).fold(f64::INFINITY, f64::min),
        "max_aspect_ratio": quality.iter().map(|q| q.max_aspect_ratio).fold(1.0, f64::max),
        "reference": reference,
        "blocks": blocks,
        "edges": grid.edges.iter().map(|e| json!({
            "vertices": e.vertices, "cells": e.cells, "length": e.length,
            "first_spacing": e.first_spacing, "last_spacing": e.last_spacing,
            "max_growth": e.max_growth,
        })).collect::<Vec<_>>(),
        "summary": summary,
        "checks": [
            "topology: indices, right-handed blocks, count chains, faces, patches",
            "every hexahedron corner Jacobian positive",
            "sum of cell volumes equals the volume the boundary encloses",
            "conforming faces, closed cells, positive volumes (hand-off mesh)",
            "no interior face with S·d ≤ 0"
        ],
    });
    Ok((report, handoff, structure))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Invariant: a template's file reads back to the same topology, and
    /// a geometric law with both ends set is refused.
    #[test]
    fn topology_file_round_trips() {
        let config = TemplateConfig::PipeOgrid {
            radius: 1.0,
            length: 2.0,
            core_fraction: None,
            core_cells: 4,
            radial_cells: 3,
            axial_cells: 2,
            wall_spacing: Some(0.02),
        };
        let file = template(&config).unwrap();
        let text = serde_json::to_string(&file).unwrap();
        let back: TopologyFile = serde_json::from_str(&text).unwrap();
        assert_eq!(back, file);
        assert_eq!(cell_count(&back).unwrap(), 2 * (16 + 4 * 12));
        let mut bad = back.clone();
        bad.edges[1].distribution = DistributionFile::Geometric {
            first: Some(0.1),
            last: Some(0.1),
        };
        assert!(bad.topology().is_err());
    }

    /// Invariant: the build validates and the hand-off reproduces the
    /// volume, patches and zone; the evidence names the block ranges.
    #[test]
    fn pipe_validates_for_foehncraft() {
        let file = template(&TemplateConfig::PipeOgrid {
            radius: 0.5,
            length: 2.0,
            core_fraction: None,
            core_cells: 6,
            radial_cells: 4,
            axial_cells: 5,
            wall_spacing: Some(0.01),
        })
        .unwrap();
        let grid = algo::build(&file.topology().unwrap()).unwrap();
        let (report, handoff, structure) = validate(&file, &grid).unwrap();
        let s = handoff.summary();
        assert_eq!(s.cells, 5 * (36 + 4 * 24));
        assert!((s.total_volume / report["volume"].as_f64().unwrap() - 1.0).abs() < 1e-12);
        assert_eq!(s.zones.len(), 1);
        let names: Vec<_> = s.patches.iter().map(|p| p.name.as_str()).collect();
        assert_eq!(names, ["inlet", "outlet", "walls"]);
        assert!((s.patches[0].area - s.patches[1].area).abs() < 1e-12);
        assert_eq!(structure["blocks"].as_array().unwrap().len(), 5);
        // 24 chords on the circle: the inscribed polygon is 1.14 % short
        let d = report["reference"]["relative_deviation"].as_f64().unwrap();
        assert!(d < 0.0 && d > -0.0115, "{d}");
    }
}
