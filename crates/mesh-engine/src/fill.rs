//! Fill the sweep regions of a region proposal with swept O-grids
//! (block-dominant meshing, BS-17, BS-18). Decided by the maintainer on a
//! card (2026-11-27): the stations of a sweep are level sets of the flow
//! potential, each meshed with the section of the pipe template (see
//! [`mesh_algorithms::fields::sweep`]). A sweep that ends at a port ends
//! there; a sweep that ends at another region is cut at the last calm
//! (flat) level set before it, so the remainder can join it on a plane.
//! The remainder (boxes, junctions, clashing parts, the stubs behind the
//! cuts) is not meshed yet: that is the polyhedral mesher with the
//! conformal join at the cuts (next step).

use std::time::Instant;

use mesh_algorithms::fields::Role;
use mesh_algorithms::fields::regions::{self, Kind};
use mesh_algorithms::fields::sweep::{self, End, Side, SweepMesh, SweepOptions};
use mesh_algorithms::remesh::ClosestSet;
use mesh_algorithms::sizing::surface::Primitive;
use mesh_core::geometry::Vec3;
use mesh_io::handoff::{CellFaces, HandoffMesh, Patch, PatchRole};
use serde_json::{Value, json};

use crate::regions::{Analysis, RegionOptions, analyse};
use crate::{EngineError, ErrorCode};

pub const FILL_FORMAT: &str = "koflcraft.sweep-fill";
pub const FILL_ALGORITHM: &str = "sweep_potential_ogrid_v1";

/// The settings of a fill.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FillOptions {
    pub region: RegionOptions,
    /// Cells around a section (multiple of 4).
    pub cells_around: usize,
    /// Rows at the layer heights.
    pub layers: usize,
    pub growth: f64,
    /// First cell height; `None`: 0.5 % of each sweep's first section
    /// diameter.
    pub first_height: Option<f64>,
    pub core_fraction: f64,
    /// Axial cell length; `None`: the wall spacing around the section.
    pub axial_size: Option<f64>,
}

impl Default for FillOptions {
    fn default() -> Self {
        let s = SweepOptions::default();
        Self {
            region: RegionOptions {
                second_pass: false,
                ..RegionOptions::default()
            },
            cells_around: s.cells_around,
            layers: s.layers,
            growth: s.growth,
            first_height: None,
            core_fraction: s.core_fraction,
            axial_size: None,
        }
    }
}

/// The swept cells and the report.
pub struct Filled {
    pub report: Value,
    /// The swept cells as one conforming mesh (walls, ports, and the cuts as
    /// interface patches); `None` when no sweep could be filled.
    pub mesh: Option<HandoffMesh>,
}

/// Outward faces of a hexahedron in VTK order (as the block mesher).
const HEX_FACES: [[usize; 4]; 6] = [
    [0, 4, 7, 3],
    [1, 2, 6, 5],
    [0, 1, 5, 4],
    [3, 7, 6, 2],
    [0, 3, 2, 1],
    [4, 5, 6, 7],
];

fn v3(p: Vec3) -> [f64; 3] {
    let r = |x: f64| (x * 1e6).round() / 1e6;
    [r(p.x), r(p.y), r(p.z)]
}

fn share(a: f64, b: f64) -> f64 {
    if b > 0.0 {
        ((a / b) * 1e4).round() / 1e4
    } else {
        0.0
    }
}

/// Largest angle between an interior face's normal and the line joining
/// its two cell centroids, in degrees, and the count above 70°.
pub fn non_orthogonality(mesh: &HandoffMesh) -> (f64, usize) {
    let cells = mesh.cell_volumes();
    let mut worst: f64 = 0.0;
    let mut over = 0;
    for f in 0..mesh.neighbour.len() {
        let ring = mesh.face(f);
        let mut n = Vec3::ZERO;
        for k in 0..ring.len() {
            n = n + mesh.points[ring[k]].cross(mesh.points[ring[(k + 1) % ring.len()]]);
        }
        let d = cells[mesh.neighbour[f]].1 - cells[mesh.owner[f]].1;
        let c = n.dot(d) / (n.length() * d.length());
        let angle = c.clamp(-1.0, 1.0).acos().to_degrees();
        worst = worst.max(angle);
        if angle > 70.0 {
            over += 1;
        }
    }
    (worst, over)
}

/// Fill the sweep regions of an STL with its patch map.
pub fn fill(stl: &[u8], map: &Value, opts: FillOptions) -> Result<Filled, EngineError> {
    let a: Analysis = analyse(stl, map, opts.region)?;
    let t1 = Instant::now();
    let c = regions::classify(
        &a.mesh,
        &a.fem,
        &a.adjacent,
        &a.d,
        &a.p,
        opts.region.classify,
    );
    let classify_seconds = t1.elapsed().as_secs_f64();
    // the STL's wall triangles for projection
    let walls: Vec<Primitive> = a
        .surface
        .facets()
        .iter()
        .zip(&a.roles)
        .filter(|(_, r)| **r == Role::Wall)
        .map(|(f, _)| Primitive::Triangle(f.triangle.vertices.map(|v| a.surface.vertices()[v])))
        .collect();
    let stl_walls = ClosestSet::new(walls);
    let t2 = Instant::now();
    let mut sweeps: Vec<(usize, SweepMesh)> = Vec::new();
    let mut skipped: Vec<Value> = Vec::new();
    for (i, r) in c.regions.iter().enumerate() {
        if r.kind != Kind::Sweep {
            continue;
        }
        let region: Vec<bool> = c.tet_region.iter().map(|&t| t == i).collect();
        let s = SweepOptions {
            cells_around: opts.cells_around,
            first_height: opts.first_height.unwrap_or(f64::NAN),
            growth: opts.growth,
            layers: opts.layers,
            radial_cells: None,
            core_fraction: opts.core_fraction,
            axial_size: opts.axial_size,
            ..SweepOptions::default()
        };
        match sweep::fill(&a.mesh, &a.p, &region, &a.adjacent, &stl_walls, s) {
            Ok(m) => sweeps.push((i, m)),
            Err(e) => skipped.push(json!({"region": i, "reason": e})),
        }
    }
    let sweep_seconds = t2.elapsed().as_secs_f64();
    // one mesh of all swept cells
    let mut names: Vec<String> = vec!["walls".into(), "inlet".into(), "outlet".into()];
    let mut roles = vec![PatchRole::Wall, PatchRole::Inlet, PatchRole::Outlet];
    let mut points: Vec<Vec3> = Vec::new();
    let mut cells: Vec<CellFaces> = Vec::new();
    let mut boundary: Vec<(Vec<usize>, usize)> = Vec::new();
    let mut sweep_json = Vec::new();
    let mut swept_volume = 0.0;
    for (i, m) in &sweeps {
        let base = points.len();
        let first_cell = cells.len();
        points.extend(m.points.iter().copied());
        for h in &m.hexes {
            cells.push(CellFaces {
                faces: HEX_FACES
                    .iter()
                    .map(|f| f.iter().map(|&s| base + h[s]).collect())
                    .collect(),
                zone: 0,
            });
        }
        let mut cut_patch = |end: usize| {
            names.push(format!(
                "cut_{i}_{}",
                if end == 0 { "start" } else { "end" }
            ));
            roles.push(PatchRole::Interface);
            names.len() - 1
        };
        let start_patch = match m.ends[0] {
            End::Port => 1,
            End::Cut => cut_patch(0),
        };
        let end_patch = match m.ends[1] {
            End::Port => 2,
            End::Cut => cut_patch(1),
        };
        for (q, side, _) in &m.boundary {
            let patch = match side {
                Side::Wall => 0,
                Side::Start => start_patch,
                Side::Finish => end_patch,
            };
            boundary.push((q.iter().map(|&v| base + v).collect(), patch));
        }
        swept_volume += m.volume;
        let st = &m.stations;
        let cut = |k: usize, end: End, calm: bool| -> Value {
            let s = &st[k];
            json!({
                "kind": match end { End::Port => "port", End::Cut => "cut" },
                "level": s.level,
                "centroid": v3(s.centroid),
                "normal": v3(s.normal),
                "area": s.area,
                "flatness": s.flatness,
                "calm": calm,
            })
        };
        sweep_json.push(json!({
            "region": i,
            "cells": m.hexes.len(),
            "first_cell": first_cell,
            "stations": st.len(),
            "section": {
                "cells_around": m.section[0],
                "rows": m.section[1],
                "core_cells": m.section[2],
                "corners": st[0].corners,
                "kind": if st[0].corners == 4 { "rectangular: ring corners on the section corners" } else { "smooth: O-grid corners at 45° + 90° m of a rotation-minimising frame" },
            },
            "layer_heights": m.layer_heights,
            "ends": [cut(0, m.ends[0], m.calm_ends[0]), cut(st.len() - 1, m.ends[1], m.calm_ends[1])],
            "feature_stations": m.feature_stations,
            "missed_feature_edges": m.missed_features,
            "max_station_flatness": st.iter().map(|s| s.flatness).fold(0.0, f64::max),
            "min_scaled_jacobian": m.min_scaled_jacobian,
            "volume": m.volume,
            "region_volume_between_ends": m.region_volume,
            "relative_volume_difference": if m.region_volume > 0.0 { m.volume / m.region_volume - 1.0 } else { 0.0 },
            "unlifted_nodes": m.unlifted,
        }));
    }
    let mut mesh = None;
    let mut validity = Value::Null;
    if !cells.is_empty() {
        let patches: Vec<Patch> = names
            .iter()
            .zip(&roles)
            .map(|(n, &role)| Patch {
                name: n.clone(),
                source_name: None,
                role,
            })
            .collect();
        let h = HandoffMesh::from_cells(points, cells, vec!["fluid".into()], &boundary, patches)
            .map_err(|e| EngineError::new(ErrorCode::InvalidConnectivity, e))?;
        let summary = h.summary();
        let (max_no, over) = non_orthogonality(&h);
        if summary.inverted_faces > 0 || summary.min_volume.is_nan() || summary.min_volume <= 0.0 {
            return Err(EngineError::new(
                ErrorCode::ValidationFailed,
                format!(
                    "the swept cells are not valid: {} inverted faces, smallest cell volume {:e}",
                    summary.inverted_faces, summary.min_volume
                ),
            ));
        }
        validity = json!({
            "cells": summary.cells,
            "faces": summary.faces,
            "points": summary.points,
            "total_volume": summary.total_volume,
            "min_volume": summary.min_volume,
            "inverted_faces": summary.inverted_faces,
            "max_relative_closure": summary.max_relative_closure,
            "max_non_orthogonality_deg": max_no,
            "faces_over_70_deg": over,
            "min_scaled_jacobian": sweeps.iter().map(|(_, m)| m.min_scaled_jacobian).fold(f64::INFINITY, f64::min),
            "patches": summary.patches,
        });
        mesh = Some(h);
    }
    let kv = c.kind_volumes();
    let report = json!({
        "format": FILL_FORMAT,
        "format_version": "1.0",
        "algorithm": FILL_ALGORITHM,
        "input": {
            "surface_sha256": a.sha,
            "facets": a.facets,
            "size": a.h,
            "cells_around": opts.cells_around,
            "layers": opts.layers,
            "growth": opts.growth,
            "first_height": opts.first_height,
            "core_fraction": opts.core_fraction,
            "axial_size": opts.axial_size,
        },
        "scaffold": {
            "tetrahedra": a.mesh.tets.len(),
            "volume": c.volume,
            "stl_volume": a.volume,
        },
        "regions": c.regions.iter().enumerate().map(|(i, r)| json!({
            "id": i,
            "kind": r.kind.name(),
            "volume_share": share(r.volume, c.volume),
            "centroid": v3(r.centroid),
        })).collect::<Vec<_>>(),
        "shares": {
            "sweep_regions": share(kv[0], c.volume),
            "swept_cells": share(swept_volume, c.volume),
            "left_for_the_polyhedral_mesher": share((c.volume - swept_volume).max(0.0), c.volume),
        },
        "sweeps": sweep_json,
        "not_swept": skipped,
        "validity": validity,
        "limitations": [
            "only the swept cells are meshed; the remainder (boxes, junctions, clashing parts and the stubs behind calm cuts) waits for the polyhedral mesher with the conformal join at the cuts",
            "interior stations follow the potential's level sets at equal centreline spacing; at a mitre corner one station is the mitre plane and the others fan around the inner corner, so cells there are short at the inner wall and long at the outer wall (see max_non_orthogonality_deg)",
            "a sweep needs one closed wall loop per level set and the same section type throughout (four sharp corners or none); other regions are listed under not_swept",
        ],
        "seconds": {
            "scaffold": a.scaffold_seconds,
            "fields": a.field_seconds,
            "classify": classify_seconds,
            "sweeps": sweep_seconds,
            "total": a.start.elapsed().as_secs_f64(),
        },
    });
    Ok(Filled { report, mesh })
}

/// Fill and write the bundle: `foehncraft.vtu` (the swept cells, walls,
/// ports and cuts), `report.json` and `mesh-evidence.json` into a new
/// directory. Fails when no sweep could be filled.
pub fn generate(
    stl: &[u8],
    map: &Value,
    opts: FillOptions,
    dir: &crate::MeshDirectory,
    name: &std::ffi::OsStr,
) -> Result<Value, EngineError> {
    let filled = fill(stl, map, opts)?;
    let Some(mesh) = filled.mesh else {
        return Err(EngineError::new(
            ErrorCode::GenerationFailed,
            format!(
                "no sweep region could be filled: {}",
                filled.report["not_swept"]
            ),
        ));
    };
    let report = filled.report;
    let json_report = serde_json::to_vec_pretty(&report)
        .map_err(|e| EngineError::new(ErrorCode::SerializationFailed, e))?;
    let fail = crate::handoff::stream_error;
    let mut bundle = dir
        .bundle(name)
        .map_err(|e| EngineError::new(ErrorCode::IoError, e))?;
    bundle
        .stream(crate::handoff::FOEHNCRAFT_VTU, |w| {
            crate::handoff::HandoffFormat::Vtu.write(
                &mesh,
                &[("koflcraft_algorithm", FILL_ALGORITHM.to_string())],
                w,
            )
        })
        .map_err(fail)?;
    bundle
        .stream("report.json", |w| w.write_all(&json_report))
        .map_err(fail)?;
    let evidence = crate::handoff::evidence(&mesh, FILL_ALGORITHM, bundle.files())?;
    bundle
        .stream(crate::handoff::MESH_EVIDENCE, |w| w.write_all(&evidence))
        .map_err(fail)?;
    bundle.finish();
    Ok(report)
}
