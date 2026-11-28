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

use mesh_algorithms::duct::{self, DuctMesh};
use mesh_algorithms::fields::Role;
use mesh_algorithms::fields::regions::{self, Kind};
use mesh_algorithms::fields::remainder::{self, Source};
use mesh_algorithms::fields::sweep::{self, CutCap, End, Side, SweepMesh, SweepOptions};
use mesh_algorithms::remesh::ClosestSet;
use mesh_algorithms::sizing::surface::Primitive;
use mesh_core::geometry::Vec3;
use mesh_io::handoff::{CellFaces, HandoffMesh, Patch, PatchRole};
use serde_json::{Value, json};

use crate::duct::{ScaffoldBuild, mesh_surface};
use crate::join;
use crate::regions::{Analysis, RegionOptions, analyse};
use crate::{EngineError, ErrorCode};
use mesh_core::geometry::PatchId;
use mesh_geometry::{ClosedSurface, SurfaceTriangle};

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
    /// Mesh the region behind the cuts with polyhedra and join it to the
    /// sweeps (MD-25).
    pub join: bool,
    /// Blended planar stations next to each mitre plane (default true);
    /// false keeps the level sets there (for the benchmark).
    pub mitre_blend: bool,
    /// Rectangular sections as an H-grid (default true); false gives them
    /// the O-grid of round sections (for the benchmark).
    pub h_grid: bool,
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
            join: true,
            mitre_blend: true,
            h_grid: true,
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
/// its two cell centroids, in degrees, and the count above 70°, for the
/// faces between swept cells (the first `swept_cells`), between the others
/// and between one of each (the cuts).
pub fn non_orthogonality(mesh: &HandoffMesh, swept_cells: usize) -> [(f64, usize); 3] {
    let cells = mesh.cell_volumes();
    let mut out = [(0.0_f64, 0usize); 3];
    for f in 0..mesh.neighbour.len() {
        let ring = mesh.face(f);
        let mut n = Vec3::ZERO;
        for k in 0..ring.len() {
            n = n + mesh.points[ring[k]].cross(mesh.points[ring[(k + 1) % ring.len()]]);
        }
        let (o, m) = (mesh.owner[f], mesh.neighbour[f]);
        let d = cells[m].1 - cells[o].1;
        let c = n.dot(d) / (n.length() * d.length());
        let angle = c.clamp(-1.0, 1.0).acos().to_degrees();
        let kind = match (o < swept_cells, m < swept_cells) {
            (true, true) => 0,
            (false, false) => 1,
            _ => 2,
        };
        out[kind].0 = out[kind].0.max(angle);
        if angle > 70.0 {
            out[kind].1 += 1;
        }
    }
    out
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
    let run = |first: Option<f64>| -> (Vec<(usize, SweepMesh)>, Vec<Value>) {
        let mut sweeps = Vec::new();
        let mut skipped = Vec::new();
        for (i, r) in c.regions.iter().enumerate() {
            if r.kind != Kind::Sweep {
                continue;
            }
            let region: Vec<bool> = c.tet_region.iter().map(|&t| t == i).collect();
            let s = SweepOptions {
                cells_around: opts.cells_around,
                first_height: first.unwrap_or(f64::NAN),
                growth: opts.growth,
                layers: opts.layers,
                radial_cells: None,
                core_fraction: opts.core_fraction,
                axial_size: opts.axial_size,
                mitre_blend: opts.mitre_blend,
                h_grid: opts.h_grid,
                ..SweepOptions::default()
            };
            match sweep::fill(&a.mesh, &a.p, &region, &a.adjacent, &stl_walls, s) {
                Ok(m) => sweeps.push((i, m)),
                Err(e) => skipped.push(json!({"region": i, "reason": e})),
            }
        }
        (sweeps, skipped)
    };
    let (mut sweeps, skipped) = run(opts.first_height);
    // the region behind the cuts grows one layer stack, so every sweep
    // takes the smallest first height
    let join =
        opts.join && opts.layers > 0 && sweeps.iter().any(|(_, m)| m.ends.contains(&End::Cut));
    let first = |m: &SweepMesh| m.layer_heights.get(1).copied().unwrap_or(0.0);
    let lowest = sweeps
        .iter()
        .map(|(_, m)| first(m))
        .fold(f64::INFINITY, f64::min);
    let (mut sweeps, skipped) =
        if join && opts.first_height.is_none() && sweeps.iter().any(|(_, m)| first(m) != lowest) {
            run(Some(lowest))
        } else {
            (std::mem::take(&mut sweeps), skipped)
        };
    // the caps of the cuts (moves the layer rows of each cut section)
    let mut caps: Vec<(usize, CutCap)> = Vec::new();
    let mut join_error: Option<String> = None;
    if join {
        'caps: for (s, (_, m)) in sweeps.iter_mut().enumerate() {
            for end in 0..2 {
                if m.ends[end] == End::Cut {
                    match m.cut_cap(end) {
                        Ok(cap) => caps.push((s, cap)),
                        Err(e) => {
                            join_error = Some(e);
                            caps.clear();
                            break 'caps;
                        }
                    }
                }
            }
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
    let mut cuts: Vec<Option<join::Cut>> = (0..caps.len()).map(|_| None).collect();
    for (s_index, (i, m)) in sweeps.iter().enumerate() {
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
        let capped = |end: usize| caps.iter().position(|(s, c)| *s == s_index && c.end == end);
        let mut patch_of_end = |end: usize| match (m.ends[end], capped(end)) {
            (End::Port, _) => Some(end + 1),
            (End::Cut, Some(_)) => None,
            (End::Cut, None) => Some(cut_patch(end)),
        };
        let start_patch = patch_of_end(0);
        let end_patch = patch_of_end(1);
        let per = m.per_station();
        for end in 0..2 {
            if let Some(k) = capped(end) {
                let station = if end == 0 { 0 } else { m.stations.len() - 1 };
                cuts[k] = Some(join::Cut {
                    patch: 0,
                    normal: caps[k].1.normal,
                    nodes: (0..per).map(|v| base + station * per + v).collect(),
                    faces: Vec::new(),
                });
            }
        }
        for (q, side, cell) in &m.boundary {
            let (patch, end) = match side {
                Side::Wall => (Some(0), None),
                Side::Start => (start_patch, Some(0)),
                Side::Finish => (end_patch, Some(1)),
            };
            match patch {
                Some(p) => boundary.push((q.iter().map(|&v| base + v).collect(), p)),
                None => {
                    if let Some(k) = end.and_then(capped)
                        && let Some(cut) = cuts[k].as_mut()
                    {
                        cut.faces
                            .push((first_cell + cell, if end == Some(0) { 4 } else { 5 }));
                    }
                }
            }
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
                "kind": if m.h_grid { "rectangular H-grid: the wall layers, then a tensor core clustered to all four walls" } else if st[0].corners == 4 { "rectangular O-grid: ring corners on the section corners" } else { "smooth: O-grid corners at 45° + 90° m of a rotation-minimising frame" },
            },
            "layer_heights": m.layer_heights,
            "ends": [cut(0, m.ends[0], m.calm_ends[0]), cut(st.len() - 1, m.ends[1], m.calm_ends[1])],
            "feature_stations": m.feature_stations,
            "missed_feature_edges": m.missed_features,
            "blended_stations": m.blended_stations,
            "unblended_mitre_sides": m.unblended_mitres,
            "max_station_flatness": st.iter().map(|s| s.flatness).fold(0.0, f64::max),
            "min_scaled_jacobian": m.min_scaled_jacobian,
            "volume": m.volume,
            "region_volume_between_ends": m.region_volume,
            "relative_volume_difference": if m.region_volume > 0.0 { m.volume / m.region_volume - 1.0 } else { 0.0 },
            "unlifted_nodes": m.unlifted,
        }));
    }
    // the region behind the cuts: polyhedra, joined at the cuts
    let sweep_cells = cells.len();
    let mut remainder_json = Value::Null;
    let mut poly_volume = 0.0;
    if !caps.is_empty() {
        let t3 = Instant::now();
        let steps: Vec<f64> = sweeps[caps[0].0]
            .1
            .layer_heights
            .windows(2)
            .map(|w| w[1] - w[0])
            .collect();
        let cap_list: Vec<CutCap> = caps.iter().map(|(_, c)| c.clone()).collect();
        match remainder_mesh(&a, &cap_list, &steps) {
            Ok((duct_mesh, cap_ids, ports, surface_json)) => {
                let cuts: Vec<join::Cut> = cuts
                    .into_iter()
                    .zip(&cap_ids)
                    .filter_map(|(c, &id)| c.map(|c| join::Cut { patch: id, ..c }))
                    .collect();
                let diameter = cap_list.iter().map(|c| 2.0 * c.radius).fold(0.0, f64::max);
                let patch_of = |id: usize| -> Option<usize> {
                    if duct::is_wall(id) {
                        Some(0)
                    } else {
                        ports.iter().find(|(p, _)| *p == id).map(|(_, m)| *m)
                    }
                };
                let joined = join::join(
                    &mut points,
                    &mut cells,
                    &mut boundary,
                    sweep_cells,
                    &cuts,
                    &duct_mesh,
                    &patch_of,
                    1e-8 * diameter,
                )
                .map_err(|e| EngineError::new(ErrorCode::InvalidConnectivity, e))?;
                poly_volume = duct_mesh.measurements.volume;
                remainder_json = json!({
                    "surface": surface_json,
                    "cells": joined.cells,
                    "cap_faces_joined": joined.pieces,
                    "edge_points_inserted": joined.edge_points,
                    "welded_points": joined.welded,
                    "max_weld_distance": joined.max_weld,
                    "max_relative_area_error": joined.area_error,
                    "volume": poly_volume,
                    "seconds": t3.elapsed().as_secs_f64(),
                });
            }
            Err(e) => {
                join_error = Some(e.message.clone());
                // the cuts stay open, as interface patches
                for (k, cut) in cuts.into_iter().enumerate() {
                    let Some(cut) = cut else { continue };
                    let (s, cap) = &caps[k];
                    names.push(format!(
                        "cut_{}_{}",
                        sweeps[*s].0,
                        if cap.end == 0 { "start" } else { "end" }
                    ));
                    roles.push(PatchRole::Interface);
                    for (c, fi) in cut.faces {
                        boundary.push((cells[c].faces[fi].clone(), names.len() - 1));
                    }
                }
            }
        }
    }
    if let Some(e) = &join_error {
        remainder_json = json!({"error": e});
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
        let no = non_orthogonality(&h, sweep_cells);
        let max_no = no.iter().map(|k| k.0).fold(0.0, f64::max);
        let over: usize = no.iter().map(|k| k.1).sum();
        if summary.inverted_faces > 0 || summary.min_volume.is_nan() || summary.min_volume <= 0.0 {
            return Err(EngineError::new(
                ErrorCode::ValidationFailed,
                format!(
                    "the cells are not valid: {} inverted faces, smallest cell volume {:e}",
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
            "non_orthogonality_by_faces": {
                "between_hexahedra": {"max_deg": no[0].0, "over_70_deg": no[0].1},
                "between_polyhedra": {"max_deg": no[1].0, "over_70_deg": no[1].1},
                "at_the_cuts": {"max_deg": no[2].0, "over_70_deg": no[2].1},
            },
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
            "polyhedral_cells": share(poly_volume, c.volume),
            "not_meshed": share((c.volume - swept_volume - poly_volume).max(0.0), c.volume),
        },
        "remainder": remainder_json,
        "sweeps": sweep_json,
        "not_swept": skipped,
        "validity": validity,
        "limitations": [
            if poly_volume > 0.0 { "the region behind the cuts is meshed as polyhedra (median dual with the sweeps' layer stack) and joined conformally at the cuts: the cut sections of the sweeps take the cap positions of the polyhedral mesh, each hex face at a cut is split into the polyhedral cap faces it holds" } else { "only the swept cells are meshed; the region behind the cuts is not (see remainder for why)" },
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

/// The polyhedral mesh behind the cuts, the patch id of every cap, the STL
/// ports inside the region as (patch id, merged patch) and a report.
type Remainder = (DuctMesh, Vec<usize>, Vec<(usize, usize)>, Value);

/// The polyhedral mesh of the region behind the cuts (MD-25): the STL
/// without the swept ducts, closed by the caps, meshed by the median dual
/// with the sweeps' layer stack. Returns the mesh, the patch id of every
/// cap, the STL ports inside the region as (patch id, merged patch) and a
/// report.
fn remainder_mesh(a: &Analysis, caps: &[CutCap], steps: &[f64]) -> Result<Remainder, EngineError> {
    let bad = |m: String| EngineError::new(ErrorCode::InvalidDomain, m);
    let triangles: Vec<[usize; 3]> = a
        .surface
        .facets()
        .iter()
        .map(|f| f.triangle.vertices)
        .collect();
    let rem = remainder::remainder_surface(a.surface.vertices(), &triangles, caps).map_err(bad)?;
    let has = |role: Role| {
        rem.source
            .iter()
            .any(|s| matches!(*s, Source::Facet(f) if a.roles[f] == role))
    };
    if has(Role::Symmetry) {
        return Err(bad(
            "symmetry planes behind the cuts are not supported by the join yet".into(),
        ));
    }
    // STL ports keep their ids (1 inlet, 2 outlet), the cuts take the free ones
    let mut free = Vec::new();
    let mut ports = Vec::new();
    for (id, role) in [(1, Role::Inlet), (2, Role::Outlet)] {
        if has(role) {
            ports.push((id, id));
        } else {
            free.push(id);
        }
    }
    // every cut cap takes a free cap id; past two openings the caps share
    // one (the mesher gives each connected cap its own similarity, and the
    // join tells caps of one id apart by their planes)
    let shared = free.first().copied().unwrap_or(1);
    let cap_ids: Vec<usize> = (0..caps.len())
        .map(|k| free.get(k).copied().unwrap_or(shared))
        .collect();
    let patches: Vec<usize> = rem
        .source
        .iter()
        .map(|s| match *s {
            Source::Facet(f) => match a.roles[f] {
                Role::Inlet => 1,
                Role::Outlet => 2,
                _ => duct::WALL,
            },
            Source::Strip(_) => duct::WALL,
            Source::Cap(k) => cap_ids[k],
        })
        .collect();
    let faces: Vec<SurfaceTriangle> = rem
        .triangles
        .iter()
        .zip(&patches)
        .map(|(t, &p)| SurfaceTriangle {
            vertices: *t,
            patch: PatchId(p - 1),
        })
        .collect();
    let surface = ClosedSurface::new(
        rem.points.clone(),
        faces,
        vec!["1".into(), "2".into(), "3".into()],
    )
    .map_err(|e| {
        bad(format!(
            "the region behind the cuts is not a closed surface: {e}"
        ))
    })?;
    // sizes: the cuts' wall spacing at the caps, the analysis spacing inside
    let spacing = caps
        .iter()
        .map(|c| {
            let rim = &c.points[..c.rim];
            (0..rim.len())
                .map(|i| (rim[(i + 1) % rim.len()] - rim[i]).length())
                .sum::<f64>()
                / rim.len() as f64
        })
        .fold(f64::INFINITY, f64::min);
    let build = ScaffoldBuild {
        tolerance: 0.01 * spacing,
        min_size: 0.125 * spacing,
        max_size: a.h.max(spacing),
        rate: 0.4,
        rim_size: Some(spacing),
        feature_angle: 30.0,
        iterations: 10,
        thickness: None,
        grading: None,
    };
    let report_build = json!({
        "tolerance": build.tolerance, "min_size": build.min_size, "max_size": build.max_size,
        "rate": build.rate, "rim_size": build.rim_size,
    });
    let (mesh, _shell, built) = mesh_surface(&surface, &patches, build, steps, 0.002)?;
    let report = json!({
        "stl_facets_dropped": rem.dropped,
        "triangles": rem.triangles.len(),
        "caps": caps.len(),
        "build": report_build,
        "scaffold": built,
    });
    Ok((mesh, cap_ids, ports, report))
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
