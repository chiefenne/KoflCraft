//! Generate the analytic block and straight-pipe development cases, in millimetres.
use crate::{EngineError, ErrorCode, MeshDirectory, SCHEMA_VERSION};
use mesh_algorithms::sweep::{Axis, Profile, Sizing, build};
use mesh_io::sweep as export;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::ffi::OsStr;

#[derive(Deserialize, Serialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum Opening {
    XMin,
    YMin,
    ZMin,
}
#[derive(Deserialize, Serialize, Clone, Copy)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Geometry {
    Block {
        dimensions_mm: [f64; 3],
        inlet: Opening,
    },
    Pipe {
        diameter_mm: f64,
        length_mm: f64,
    },
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SweepConfig {
    pub schema_version: u32,
    pub geometry: Geometry,
    pub core_spacing_mm: f64,
    pub wall_layers: usize,
    pub first_layer_mm: f64,
    pub layer_growth: f64,
    pub max_neighbour_size_ratio: f64,
}
pub fn generate(
    config: &SweepConfig,
    dir: &MeshDirectory,
    name: &OsStr,
    cell_budget: usize,
) -> Result<Value, EngineError> {
    if config.schema_version != SCHEMA_VERSION {
        return Err(EngineError::new(
            ErrorCode::UnsupportedSchema,
            "unsupported sweep configuration version",
        ));
    }
    let (profile, length, axis) = match config.geometry {
        Geometry::Block {
            dimensions_mm: [x, y, z],
            inlet,
        } => {
            let (width, height, depth, axis) = match inlet {
                Opening::XMin => (y, z, x, Axis::X),
                Opening::YMin => (z, x, y, Axis::Y),
                Opening::ZMin => (x, y, z, Axis::Z),
            };
            (Profile::BlindBlock { width, height }, depth, axis)
        }
        Geometry::Pipe {
            diameter_mm,
            length_mm,
        } => (
            Profile::Pipe {
                diameter: diameter_mm,
            },
            length_mm,
            Axis::X,
        ),
    };
    let size = Sizing {
        core: config.core_spacing_mm,
        layers: config.wall_layers,
        first: config.first_layer_mm,
        growth: config.layer_growth,
        ratio_limit: config.max_neighbour_size_ratio,
    };
    let mesh = build(profile, length, &size, cell_budget)
        .map_err(|e| EngineError::new(ErrorCode::InvalidDomain, e))?;
    let mut core_sizes = Vec::new();
    let mut ratios = Vec::with_capacity(mesh.cells());
    let mut volume = 0.;
    for k in 0..mesh.axial_layers.len() {
        for t in 0..mesh.tiles.len() {
            let v = mesh.volume(k, t);
            volume += v;
            if mesh.layer(k, t) == 0 {
                core_sizes.push(v.cbrt());
            }
            ratios.push(mesh.neighbour_ratio(k, t));
        }
    }
    if (volume / (mesh.fitted_area * length) - 1.).abs() > 1e-10 {
        return Err(EngineError::new(
            ErrorCode::ValidationFailed,
            "volume conservation failed",
        ));
    }
    core_sizes.sort_unstable_by(f64::total_cmp);
    ratios.sort_unstable_by(f64::total_cmp);
    let mut areas = vec![0.; config.wall_layers + 1];
    for tile in &mesh.tiles {
        areas[tile.layer] += tile.area;
    }
    let report = json!({
        "schema_version": SCHEMA_VERSION,
        "algorithm": "swept_voronoi_v1",
        "config": config,
        "cells": mesh.cells(),
        "core_cells": core_sizes.len(),
        "wall_layer_cells": mesh.cells() - core_sizes.len(),
        "axial_intervals": mesh.axial_layers.len(),
        "core_axial_intervals": mesh.axial_layers.iter().filter(|&&layer| layer == 0).count(),
        "cross_section_core_cells": mesh.tiles.iter().filter(|tile| tile.layer == 0).count(),
        "core_cell_size_mm": {
            "min": core_sizes[0],
            "median": core_sizes[core_sizes.len()/2],
            "max": core_sizes.last()
        },
        "neighbour_size_ratio": {
            "p95_cell_max": ratios[(ratios.len()-1)*95/100],
            "max": ratios.last(),
            "limit": size.ratio_limit
        },
        "max_internal_non_orthogonality_deg": mesh.tiles.iter().map(|t| t.non_orthogonality).fold(0., f64::max),
        "max_internal_skewness": mesh.tiles.iter().map(|t| t.skewness).fold(0., f64::max),
        "wall_thicknesses_mm": mesh.thicknesses,
        "inlet_area_by_layer_mm2": areas,
        "outlet_area_by_layer_mm2": (!profile.blind()).then_some(&areas),
        "volume_mm3": volume,
        "analytic_volume_mm3": profile.area()*length,
        "relative_geometry_volume_error": (mesh.fitted_area/profile.area()-1.).abs(),
        "profile_sides": mesh.rim_sides
    });
    let points = mesh.vertices(axis);
    let io = |e| EngineError::new(ErrorCode::ExportFailed, e);
    let mut vtu = Vec::new();
    export::volume(&mesh, &points, &mut vtu).map_err(io)?;
    let mut gmsh = Vec::new();
    export::gmsh(&mesh, &points, mesh_io::GmshVersion::V4_1, &mut gmsh).map_err(io)?;
    let [inlet, outlet, walls] = export::patches(&mesh, &points).map_err(io)?;
    let multiblock = if profile.blind() {
        export::BLIND_MULTIBLOCK
    } else {
        mesh_io::channel::MULTIBLOCK
    };
    let json = serde_json::to_vec_pretty(&report)
        .map_err(|e| EngineError::new(ErrorCode::SerializationFailed, e))?;
    let mut files: Vec<(&str, &[u8])> = vec![
        ("mesh.vtm", multiblock),
        ("volume.vtu", &vtu),
        ("mesh.msh", &gmsh),
        ("inlet.vtu", &inlet),
        ("walls.vtu", &walls),
        ("report.json", &json),
    ];
    if !profile.blind() {
        files.push(("outlet.vtu", &outlet));
    }
    dir.write_bundle(name, &files)
        .map_err(|e| EngineError::new(ErrorCode::IoError, e))?;
    Ok(report)
}
