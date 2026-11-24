//! Generate and certify the exposed port layers of an analytic straight channel.

use std::ffi::OsStr;

use mesh_algorithms::channel::{ChannelMesh, ChannelPatch, straight_channel};
use mesh_core::geometry::Vec3;
use mesh_io::channel::{MULTIBLOCK, write_gmsh, write_patch, write_volume};
use mesh_quality::{MeshMetrics, RELATIVE_CONSERVATION_TOLERANCE, validate_mesh};
use serde::{Deserialize, Serialize};

use crate::handoff::HandoffFormat;
use crate::{EngineError, ErrorCode, MeshDirectory, SCHEMA_VERSION};

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ChannelConfig {
    pub schema_version: u32,
    pub length: f64,
    pub width: f64,
    pub height: f64,
    pub axial_cells: usize,
    pub layers: usize,
    pub first_thickness: f64,
    pub growth: f64,
}

#[derive(Serialize)]
pub struct PortReport {
    pub faces: usize,
    pub core_area: f64,
    /// In wall-to-core order.
    pub layer_areas: Vec<f64>,
    pub total_area: f64,
}

#[derive(Serialize)]
pub struct ChannelReport {
    pub cells: usize,
    pub prisms: usize,
    pub core_cells: usize,
    pub total_volume: f64,
    pub minimum_cell_volume: f64,
    pub wall_distances: Vec<f64>,
    pub inlet: PortReport,
    pub outlet: PortReport,
}

fn check_port(
    channel: &ChannelMesh,
    metrics: &MeshMetrics,
    patch: ChannelPatch,
) -> Result<PortReport, EngineError> {
    let insets = channel.insets();
    let mut areas = vec![0.0; insets.len()];
    let mut faces = 0;
    for (face, owner, _) in channel.boundary().filter(|&(_, _, p)| p == patch) {
        areas[channel.cells()[owner].layer()] += metrics.faces[face].area;
        faces += 1;
    }
    let size = channel.bounds().extent();
    let core = insets.last().unwrap();
    let core_area = areas[0];
    let total_area = areas.iter().sum();
    let agrees =
        |area: f64, expected: f64| (area / expected - 1.0).abs() <= RELATIVE_CONSERVATION_TOLERANCE;
    // Difference of nested rectangle areas, factored to retain thin rings.
    let valid_layers = insets.windows(2).zip(&areas[1..]).all(|(d, &area)| {
        agrees(
            area,
            2.0 * (d[1] - d[0]) * ((size.y - 2.0 * d[0]) + (size.z - 2.0 * d[1])),
        )
    });
    if !valid_layers
        || !agrees(total_area, size.y * size.z)
        || !agrees(core_area, (size.y - 2.0 * core) * (size.z - 2.0 * core))
    {
        return Err(EngineError::new(
            ErrorCode::ValidationFailed,
            format!(
                "{} core or layer area differs from its prescribed opening",
                patch.name()
            ),
        ));
    }
    Ok(PortReport {
        faces,
        core_area,
        layer_areas: areas[1..].to_vec(),
        total_area,
    })
}

/// Write a new self-contained directory only after topology, volume and both
/// port area partitions have passed. The destination follows the CLI's create-new semantics.
/// The number of cells `generate_channel` builds (AG-15): an
/// (2 layers + 1)² grid of columns per axial cell, every column but the
/// core split into two prisms.
pub fn cell_count(config: &ChannelConfig) -> Result<usize, EngineError> {
    let overflow = || EngineError::new(ErrorCode::InvalidDomain, "grid size overflow");
    let columns = config
        .layers
        .checked_mul(2)
        .and_then(|v| v.checked_add(1))
        .and_then(|v| v.checked_pow(2))
        .ok_or_else(overflow)?;
    columns
        .checked_mul(2)
        .and_then(|v| v.checked_sub(1))
        .and_then(|v| v.checked_mul(config.axial_cells))
        .ok_or_else(overflow)
}

pub fn generate_channel(
    config: &ChannelConfig,
    dir: &MeshDirectory,
    name: &OsStr,
    handoff_format: HandoffFormat,
) -> Result<ChannelReport, EngineError> {
    if config.schema_version != SCHEMA_VERSION {
        return Err(EngineError::new(
            ErrorCode::UnsupportedSchema,
            "unsupported channel config version",
        ));
    }
    let channel = straight_channel(
        Vec3::new(config.length, config.width, config.height),
        config.axial_cells,
        config.layers,
        config.first_thickness,
        config.growth,
    )
    .map_err(|e| EngineError::new(ErrorCode::InvalidDomain, e))?;
    let metrics = validate_mesh(channel.bounds(), channel.mesh())
        .map_err(|e| EngineError::new(ErrorCode::ValidationFailed, e))?;
    let report = ChannelReport {
        cells: channel.cells().len(),
        prisms: channel.cells().len() - config.axial_cells,
        core_cells: config.axial_cells,
        total_volume: metrics.total_volume,
        minimum_cell_volume: metrics
            .cells
            .iter()
            .map(|c| c.volume)
            .fold(f64::INFINITY, f64::min),
        wall_distances: channel.insets().to_vec(),
        inlet: check_port(&channel, &metrics, ChannelPatch::Inlet)?,
        outlet: check_port(&channel, &metrics, ChannelPatch::Outlet)?,
    };
    let handoff = mesh_io::channel::handoff(&channel).map_err(crate::handoff::export_error)?;
    // each file goes to disk as it is written, hashed on the way
    let fail = crate::handoff::stream_error;
    let mut bundle = dir
        .bundle(name)
        .map_err(|e| EngineError::new(ErrorCode::IoError, e))?;
    bundle
        .stream("volume.vtu", |w| write_volume(&channel, w))
        .map_err(fail)?;
    for (patch, file) in [
        (ChannelPatch::Inlet, "inlet.vtu"),
        (ChannelPatch::Outlet, "outlet.vtu"),
        (ChannelPatch::Walls, "walls.vtu"),
    ] {
        bundle
            .stream(file, |w| write_patch(&channel, patch, w))
            .map_err(fail)?;
    }
    bundle
        .stream("channel.vtm", |w| w.write_all(MULTIBLOCK))
        .map_err(fail)?;
    bundle
        .stream("channel.msh", |w| {
            write_gmsh(&channel, mesh_io::GmshVersion::V4_1, w)
        })
        .map_err(fail)?;
    bundle
        .stream(handoff_format.file_name(), |w| {
            handoff_format.write(&handoff, &[], w)
        })
        .map_err(fail)?;
    let evidence = crate::handoff::evidence(&handoff, "straight_channel", bundle.files())?;
    bundle
        .stream(crate::handoff::MESH_EVIDENCE, |w| w.write_all(&evidence))
        .map_err(fail)?;
    bundle.finish();
    Ok(report)
}
