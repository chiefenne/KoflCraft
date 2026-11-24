//! The block-structured variant behind the engine (OBJ-21, BS-01..BS-14):
//! the topology file `koflcraft.block-topology` (BS-02), the templates
//! (BS-03), the check (BS-05) and the build of a validated hexahedral mesh
//! into a new bundle (BS-13, BS-14).

use std::ffi::OsStr;

use mesh_algorithms::block::smooth::{SmoothOptions, smooth};
use mesh_algorithms::block::{self as algo, BlockGrid, Issue, templates};
use mesh_core::block::{
    BlockSpec, BlockTopology, Curve, Distribution, EdgeEnd, EdgeSpec, PatchSpec, Surface,
    SurfaceSpec,
};
use mesh_core::geometry::Vec3;
use mesh_io::handoff::{PatchRole, patches_from_names};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::{EngineError, ErrorCode, MeshDirectory};

pub const TOPOLOGY_FORMAT: &str = "koflcraft.block-topology";
pub const TOPOLOGY_FORMAT_VERSION: &str = "1.0";
pub const ALGORITHM: &str = "block_tfi_v1";
pub const ELLIPTIC_ALGORITHM: &str = "block_elliptic_v1";

/// The volume algorithm of `block.grid` (BS-09), a choice separate from
/// the topology.
#[derive(Clone, Debug, Default, Deserialize, Serialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum VolumeAlgorithm {
    /// Transfinite interpolation only.
    #[default]
    Tfi,
    /// TFI, then Thompson–Thames–Mastin elliptic smoothing of every node
    /// off the boundary, across block faces.
    Elliptic {
        /// Orthogonality and first-spacing control at the faces of the
        /// patches with the role `wall` (default true).
        #[serde(default = "yes")]
        wall_orthogonality: bool,
        /// Let the wall nodes slide along the curved wall contours toward
        /// orthogonal, smooth grid lines; corners sharper than 30° stay
        /// fixed (default false). Needs `wall_orthogonality`.
        #[serde(default, skip_serializing_if = "std::ops::Not::not")]
        wall_sliding: bool,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        max_iterations: Option<usize>,
    },
}

fn yes() -> bool {
    true
}
/// Cells and points of a sum of cell volumes and of the volume the boundary
/// encloses may differ by this fraction at most (rounding only).
pub const CONSERVATION_TOLERANCE: f64 = 1e-10;

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum CurveFile {
    Line,
    Arc {
        through: [f64; 3],
    },
    /// Natural cubic spline through these points between the two vertices.
    Spline {
        through: Vec<[f64; 3]>,
    },
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

/// An analytic surface in the file.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum SurfaceShape {
    Cylinder {
        point: [f64; 3],
        axis: [f64; 3],
        radius: f64,
    },
    Torus {
        centre: [f64; 3],
        axis: [f64; 3],
        major: f64,
        minor: f64,
    },
}

/// A surface and the boundary block faces projected onto it (BS-07).
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct SurfaceFile {
    pub name: String,
    pub shape: SurfaceShape,
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
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub surfaces: Vec<SurfaceFile>,
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
    /// A circular pipe elbow: inlet leg along +x ending at the origin, a
    /// bend about +z through (0, bend_radius, 0), outlet leg; the pipe
    /// O-grid swept along it, walls projected onto cylinders and the torus.
    Elbow {
        radius: f64,
        /// Radius of the bend's centreline.
        bend_radius: f64,
        /// Bend angle in degrees, at most 180.
        angle: f64,
        #[serde(default)]
        inlet_length: f64,
        #[serde(default)]
        outlet_length: f64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        core_fraction: Option<f64>,
        core_cells: usize,
        radial_cells: usize,
        /// Cells along the centreline in total, shared out by length.
        axial_cells: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        wall_spacing: Option<f64>,
    },
    /// An airfoil C-grid extruded along +z: leading edge at the origin,
    /// chord along +x, a half circle of radius `farfield` about the
    /// trailing edge and wake blocks behind it.
    AirfoilCgrid {
        /// NACA 4-digit designation, e.g. "0012" or "2412" (closed
        /// trailing edge); exactly one of `naca` and `coordinates`.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        naca: Option<String>,
        /// A coordinate file (Selig or Lednicer); a blunt trailing edge is
        /// closed by a linear thickness correction.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        coordinates: Option<String>,
        #[serde(default = "one")]
        chord: f64,
        span: f64,
        farfield: f64,
        wake_length: f64,
        /// Cells along each side of the airfoil.
        surface_cells: usize,
        normal_cells: usize,
        wake_cells: usize,
        #[serde(default = "one_cell")]
        span_cells: usize,
        /// First cell height at the airfoil and the wake cut.
        wall_spacing: f64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        leading_edge_spacing: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        trailing_edge_spacing: Option<f64>,
    },
}

fn one() -> f64 {
    1.0
}
fn one_cell() -> usize {
    1
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
                    curve: match &e.curve {
                        CurveFile::Line => Curve::Line,
                        CurveFile::Arc { through } => Curve::Arc {
                            through: v3(*through),
                        },
                        CurveFile::Spline { through } => Curve::Spline {
                            through: through.iter().copied().map(v3).collect(),
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
            surfaces: self
                .surfaces
                .iter()
                .map(|s| SurfaceSpec {
                    name: s.name.clone(),
                    surface: match s.shape {
                        SurfaceShape::Cylinder {
                            point,
                            axis,
                            radius,
                        } => Surface::Cylinder {
                            point: v3(point),
                            axis: v3(axis),
                            radius,
                        },
                        SurfaceShape::Torus {
                            centre,
                            axis,
                            major,
                            minor,
                        } => Surface::Torus {
                            centre: v3(centre),
                            axis: v3(axis),
                            major,
                            minor,
                        },
                    },
                    faces: s.faces.clone(),
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
                    curve: match &e.curve {
                        Curve::Line => CurveFile::Line,
                        Curve::Arc { through } => CurveFile::Arc {
                            through: a3(*through),
                        },
                        Curve::Spline { through } => CurveFile::Spline {
                            through: through.iter().copied().map(a3).collect(),
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
            surfaces: t
                .surfaces
                .iter()
                .map(|s| SurfaceFile {
                    name: s.name.clone(),
                    shape: match s.surface {
                        Surface::Cylinder {
                            point,
                            axis,
                            radius,
                        } => SurfaceShape::Cylinder {
                            point: a3(point),
                            axis: a3(axis),
                            radius,
                        },
                        Surface::Torus {
                            centre,
                            axis,
                            major,
                            minor,
                        } => SurfaceShape::Torus {
                            centre: a3(centre),
                            axis: a3(axis),
                            major,
                            minor,
                        },
                    },
                    faces: s.faces.clone(),
                })
                .collect(),
            reference,
        }
    }
}

/// Write the topology of a template (`block.template`).
pub fn template(config: &TemplateConfig) -> Result<TopologyFile, EngineError> {
    template_with(config, None)
}

/// As [`template`], with the text of the airfoil coordinate file that an
/// `airfoil_cgrid` config names.
pub fn template_with(
    config: &TemplateConfig,
    coordinates_text: Option<&str>,
) -> Result<TopologyFile, EngineError> {
    let (t, name) = match *config {
        TemplateConfig::AirfoilCgrid {
            ref naca,
            ref coordinates,
            chord,
            span,
            farfield,
            wake_length,
            surface_cells,
            normal_cells,
            wake_cells,
            span_cells,
            wall_spacing,
            leading_edge_spacing,
            trailing_edge_spacing,
        } => (
            match (naca, coordinates, coordinates_text) {
                (Some(n), None, _) => templates::Naca4::parse(n).map(|s| s.section()),
                (None, Some(_), Some(text)) => templates::Section::from_coordinates(text),
                (None, Some(path), None) => Err(format!("the coordinate file {path} was not read")),
                _ => Err("an airfoil needs exactly one of naca and coordinates".into()),
            }
            .and_then(|section| {
                templates::airfoil_cgrid(
                    &section,
                    &templates::AirfoilCgrid {
                        chord,
                        span,
                        farfield,
                        wake: wake_length,
                        cells: [surface_cells, normal_cells, wake_cells, span_cells],
                        wall_spacing,
                        le_spacing: leading_edge_spacing,
                        te_spacing: trailing_edge_spacing,
                    },
                )
            }),
            "airfoil_cgrid",
        ),
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
        TemplateConfig::Elbow {
            radius,
            bend_radius,
            angle,
            inlet_length,
            outlet_length,
            core_fraction,
            core_cells,
            radial_cells,
            axial_cells,
            wall_spacing,
        } => (
            templates::elbow(
                radius,
                bend_radius,
                angle,
                inlet_length,
                outlet_length,
                core_fraction.unwrap_or(0.55),
                [core_cells, radial_cells, axial_cells],
                wall_spacing,
            ),
            "elbow",
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
    algorithm: &VolumeAlgorithm,
    dir: &MeshDirectory,
    name: &OsStr,
) -> Result<Value, EngineError> {
    let (grid, smoothing) = grid(file, algorithm)?;
    let (mut report, handoff, structure) = validate(file, &grid)?;
    report["volume_algorithm"] = json!(algorithm);
    report["algorithm"] = json!(algorithm_name(algorithm));
    report["smoothing"] = smoothing;
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
                &[("koflcraft_algorithm", algorithm_name(algorithm).to_string())],
                w,
            )
        })
        .map_err(fail)?;
    bundle
        .stream("report.json", |w| w.write_all(&json_report))
        .map_err(fail)?;
    let evidence = crate::handoff::evidence_with(
        &handoff,
        algorithm_name(algorithm),
        bundle.files(),
        Some(structure),
    )?;
    bundle
        .stream(crate::handoff::MESH_EVIDENCE, |w| w.write_all(&evidence))
        .map_err(fail)?;
    bundle.finish();
    Ok(report)
}

fn algorithm_name(algorithm: &VolumeAlgorithm) -> &'static str {
    match algorithm {
        VolumeAlgorithm::Tfi => ALGORITHM,
        VolumeAlgorithm::Elliptic { .. } => ELLIPTIC_ALGORITHM,
    }
}

/// The patches with the role `wall` (given, or from the name).
fn wall_patches(file: &TopologyFile) -> Vec<usize> {
    role_patches(file, PatchRole::Wall)
}

/// The patches with a role (given, or from the name).
fn role_patches(file: &TopologyFile, role: PatchRole) -> Vec<usize> {
    file.patches
        .iter()
        .enumerate()
        .filter(|(_, p)| p.role.unwrap_or_else(|| PatchRole::from_name(&p.name)) == role)
        .map(|(i, _)| i)
        .collect()
}

/// The grid of a topology with the chosen volume algorithm, and what the
/// smoother did (for TFI: the wall angles of the TFI grid).
pub fn grid(
    file: &TopologyFile,
    algorithm: &VolumeAlgorithm,
) -> Result<(BlockGrid, Value), EngineError> {
    let t = file.topology()?;
    let mut grid = algo::build(&t).map_err(issues_error)?;
    let walls = wall_patches(file);
    let mut options = SmoothOptions {
        wall_patches: walls.clone(),
        sliding_patches: role_patches(file, PatchRole::Symmetry),
        ..SmoothOptions::default()
    };
    match *algorithm {
        VolumeAlgorithm::Tfi => options.max_iterations = 0,
        VolumeAlgorithm::Elliptic {
            wall_orthogonality,
            wall_sliding,
            max_iterations,
        } => {
            if wall_sliding && !wall_orthogonality {
                return Err(EngineError::new(
                    ErrorCode::InvalidArguments,
                    "wall_sliding needs wall_orthogonality; set wall_orthogonality to true or wall_sliding to false",
                ));
            }
            options.slide_walls = wall_sliding;
            if wall_orthogonality && walls.is_empty() {
                // BS-12: a request the algorithm cannot honour is refused
                return Err(EngineError::new(
                    ErrorCode::InvalidArguments,
                    "wall_orthogonality needs a patch with the role wall; set wall_orthogonality to false or give a patch that role",
                ));
            }
            if !wall_orthogonality {
                options.wall_patches.clear();
            }
            if let Some(n) = max_iterations {
                options.max_iterations = n;
            }
        }
    }
    let s = smooth(&mut grid, &options)
        .map_err(|e| EngineError::new(ErrorCode::GenerationFailed, e))?;
    if options.max_iterations > 0 && !s.converged {
        return Err(EngineError::new(
            ErrorCode::GenerationFailed,
            format!(
                "elliptic smoothing did not converge in {} iterations (last change {:e} of the shortest edge); raise max_iterations or use kind tfi",
                s.iterations, s.last_change
            ),
        ));
    }
    let report = json!({
        "iterations": s.iterations, "converged": s.converged, "last_change": s.last_change,
        "moved_block_nodes": if options.max_iterations > 0 { s.block_nodes } else { 0 },
        "moved_face_nodes": if options.max_iterations > 0 { s.face_nodes } else { 0 },
        "moved_edge_nodes": if options.max_iterations > 0 { s.edge_nodes } else { 0 },
        "moved_vertex_nodes": if options.max_iterations > 0 { s.vertex_nodes } else { 0 },
        "sliding_patches": s.sliding_patches,
        "moved_sliding_nodes": if options.max_iterations > 0 { s.sliding_nodes } else { 0 },
        "wall_faces": s.wall_faces,
        "wall_angle_deg": {"max_tfi": s.wall_angle_max_before, "max": s.wall_angle_max_after,
                           "mean_tfi": s.wall_angle_mean_before, "mean": s.wall_angle_mean_after},
        "first_spacing_change": s.first_spacing_change,
        "contour_nodes": s.contour_nodes,
        "contour_slide_max": s.contour_slide_max,
    });
    Ok((grid, report))
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
            "block {} has a hexahedron with the corner Jacobian {:e}; the volume algorithm folds this block",
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
        let (grid, _) = super::grid(&file, &VolumeAlgorithm::Tfi).unwrap();
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

    /// Invariant (BS-07): the elbow's surfaces survive the file, it
    /// validates for FoehnCraft, and its volume is within the chord error
    /// of the Pappus value.
    #[test]
    fn elbow_validates_with_its_surfaces() {
        let file = template(&TemplateConfig::Elbow {
            radius: 0.5,
            bend_radius: 1.5,
            angle: 90.0,
            inlet_length: 1.0,
            outlet_length: 0.5,
            core_fraction: None,
            core_cells: 6,
            radial_cells: 4,
            axial_cells: 24,
            wall_spacing: Some(0.01),
        })
        .unwrap();
        let text = serde_json::to_string(&file).unwrap();
        let back: TopologyFile = serde_json::from_str(&text).unwrap();
        assert_eq!(back, file);
        let names: Vec<_> = file.surfaces.iter().map(|s| s.name.as_str()).collect();
        assert_eq!(names, ["cylinder_0", "torus", "cylinder_4"]);
        let (grid, _) = super::grid(&file, &VolumeAlgorithm::Tfi).unwrap();
        let (report, handoff, _) = validate(&file, &grid).unwrap();
        assert_eq!(handoff.summary().cells, 24 * (36 + 4 * 24));
        // 24 chords round the circle (1.14 %) plus the chords along the
        // bend: 1.25 % short
        let d = report["reference"]["relative_deviation"].as_f64().unwrap();
        assert!(d < -0.0114 && d > -0.0135, "{d}");
    }

    /// Invariant (BS-03, BS-19): the extruded NACA 0012 C-grid writes a
    /// file with spline edges that reads back and validates for FoehnCraft
    /// with its five patches; with one span cell the elliptic smoother
    /// moves nothing and says so.
    #[test]
    fn airfoil_cgrid_validates_and_smooths() {
        let file = template(&TemplateConfig::AirfoilCgrid {
            naca: Some("0012".into()),
            coordinates: None,
            chord: 1.0,
            span: 0.1,
            farfield: 5.0,
            wake_length: 5.0,
            surface_cells: 40,
            normal_cells: 20,
            wake_cells: 20,
            span_cells: 1,
            wall_spacing: 2e-3,
            leading_edge_spacing: None,
            trailing_edge_spacing: None,
        })
        .unwrap();
        let text = serde_json::to_string(&file).unwrap();
        let back: TopologyFile = serde_json::from_str(&text).unwrap();
        assert_eq!(back, file);
        let (grid, _) = super::grid(&file, &VolumeAlgorithm::Tfi).unwrap();
        let (_, handoff, _) = validate(&file, &grid).unwrap();
        let names: Vec<_> = handoff
            .summary()
            .patches
            .iter()
            .map(|p| p.name.clone())
            .collect();
        assert_eq!(
            names,
            [
                "airfoil_wall",
                "inlet",
                "outlet",
                "symmetry_low",
                "symmetry_high"
            ]
        );
        // one cell along the span: every node is on a boundary; the nodes
        // on the two symmetry planes slide, so the smoother still works and
        // pulls the grid lines toward the wall normal
        let (smooth, report) = super::grid(
            &file,
            &VolumeAlgorithm::Elliptic {
                wall_orthogonality: true,
                wall_sliding: false,
                max_iterations: None,
            },
        )
        .unwrap();
        assert_eq!(report["converged"], true);
        assert_eq!(report["sliding_patches"], 2);
        assert!(report["moved_sliding_nodes"].as_u64().unwrap() > 0);
        let angle = &report["wall_angle_deg"];
        assert!(angle["max"].as_f64().unwrap() < 0.5 * angle["max_tfi"].as_f64().unwrap());
        assert!(angle["mean"].as_f64().unwrap() < 0.5 * angle["mean_tfi"].as_f64().unwrap());
        // the sliding nodes stay on their planes z = 0 and z = span
        for p in &smooth.points {
            let on_low = p.z.abs() < 1e-9;
            let on_high = (p.z - 0.1).abs() < 1e-9;
            assert!(on_low || on_high, "node left its symmetry plane: {p:?}");
        }
        // wall sliding: the nodes of the section slide and are reported;
        // it needs the wall controls
        let sliding = |wall_orthogonality| VolumeAlgorithm::Elliptic {
            wall_orthogonality,
            wall_sliding: true,
            max_iterations: None,
        };
        let (_, slid) = super::grid(&file, &sliding(true)).unwrap();
        assert_eq!(slid["converged"], true);
        assert!(slid["contour_nodes"].as_u64().unwrap() > 0, "{slid}");
        assert!(slid["contour_slide_max"].as_f64().unwrap() > 0.0);
        let e = super::grid(&file, &sliding(false)).unwrap_err();
        assert!(e.message.contains("wall_sliding needs wall_orthogonality"));
        let back: VolumeAlgorithm =
            serde_json::from_str(r#"{"kind": "elliptic", "wall_sliding": true}"#).unwrap();
        assert_eq!(back, sliding(true));
    }

    /// Invariant (BS-12): elliptic smoothing with wall orthogonality on a
    /// topology without a wall patch is refused, not run without it; with
    /// walls it validates and reports smaller wall angles than TFI.
    #[test]
    fn elliptic_needs_walls_and_improves_them() {
        let mut file = template(&TemplateConfig::PipeOgrid {
            radius: 0.5,
            length: 1.0,
            core_fraction: None,
            core_cells: 6,
            radial_cells: 6,
            axial_cells: 3,
            wall_spacing: Some(0.01),
        })
        .unwrap();
        let elliptic = VolumeAlgorithm::Elliptic {
            wall_orthogonality: true,
            wall_sliding: false,
            max_iterations: None,
        };
        let (grid_e, s) = super::grid(&file, &elliptic).unwrap();
        let (_, t) = super::grid(&file, &VolumeAlgorithm::Tfi).unwrap();
        validate(&file, &grid_e).unwrap();
        let angle = |v: &Value, k: &str| v["wall_angle_deg"][k].as_f64().unwrap();
        assert!(angle(&s, "max") < 0.5 * angle(&t, "max"), "{s} {t}");
        assert_eq!(t["iterations"], 0);
        file.patches[2].role = Some(PatchRole::Unspecified);
        file.patches[2].name = "rim".into();
        let e = super::grid(&file, &elliptic).unwrap_err();
        assert_eq!(e.code, ErrorCode::InvalidArguments);
        super::grid(
            &file,
            &VolumeAlgorithm::Elliptic {
                wall_orthogonality: false,
                wall_sliding: false,
                max_iterations: None,
            },
        )
        .unwrap();
    }
}
