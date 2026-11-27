//! Field-guided region proposal for block-dominant meshing: from a closed
//! STL and its patch roles, build a tet scaffold, solve the heat-method
//! wall distance and the inlet-to-outlet potential on it, find the
//! skeleton and propose which regions a block template can fill (sweep,
//! junction, obstacle, box) and which stay for the polyhedral mesher
//! (remainder). A second pass classifies every box and remainder region
//! again along its own long axis (Fiedler vector). Nothing is meshed: the
//! result is a report, a field file and a picture.

use std::fmt::Write;
use std::time::Instant;

use mesh_algorithms::fields::regions::{self, Classification, ClassifyOptions, Kind};
use mesh_algorithms::fields::{self, Fem, Role, TetMesh};
use mesh_algorithms::{remesh, scaffold};
use mesh_core::geometry::Vec3;
use serde_json::{Value, json};

use crate::{EngineError, ErrorCode};

pub const REGION_FORMAT: &str = "koflcraft.region-proposal";
/// Target scaffold size when no spacing is given.
const DEFAULT_TETS: f64 = 150_000.0;
const MAX_POINTS: usize = 1_500_000;

/// The settings of a proposal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RegionOptions {
    /// Scaffold spacing; `None`: about 150 000 tetrahedra.
    pub size: Option<f64>,
    /// Heat-method time step factor m in t = (m h)^2.
    pub heat_m: f64,
    pub classify: ClassifyOptions,
    pub second_pass: bool,
}

impl Default for RegionOptions {
    fn default() -> Self {
        Self {
            size: None,
            heat_m: 3.0,
            classify: ClassifyOptions::default(),
            second_pass: true,
        }
    }
}

/// The proposal and its files.
pub struct Proposal {
    pub report: Value,
    pub vtu: String,
    pub svg: String,
}

fn bad(code: ErrorCode, m: impl std::fmt::Display) -> EngineError {
    EngineError::new(code, m)
}

/// The role of each facet from a patch map (`facet_patch` + `patches`)
/// or a port file (`default_patch` + `ports` with `facet_range`). Roles
/// other than inlet, outlet and symmetry count as walls.
pub fn facet_roles(map: &Value, facets: usize, sha256: &str) -> Result<Vec<Role>, EngineError> {
    let invalid = |m: String| bad(ErrorCode::InvalidPatchMap, m);
    if let Some(s) = map["source_sha256"].as_str()
        && s != sha256
    {
        return Err(invalid(
            "the patch map belongs to another STL (source_sha256 differs)".into(),
        ));
    }
    let role = |r: &str| match r {
        "inlet" => Role::Inlet,
        "outlet" => Role::Outlet,
        "symmetry" => Role::Symmetry,
        _ => Role::Wall,
    };
    if let (Some(patches), Some(fp)) = (map["patches"].as_array(), map["facet_patch"].as_array()) {
        if fp.len() != facets {
            return Err(invalid(format!(
                "the patch map has {} facets, the STL {facets}",
                fp.len()
            )));
        }
        let roles: Vec<Role> = patches
            .iter()
            .map(|p| role(p["role"].as_str().unwrap_or("wall")))
            .collect();
        return fp
            .iter()
            .enumerate()
            .map(|(f, p)| {
                p.as_u64()
                    .and_then(|p| roles.get(p as usize).copied())
                    .ok_or_else(|| invalid(format!("facet {f} has no valid patch")))
            })
            .collect();
    }
    if let Some(ports) = map["ports"].as_array() {
        let mut out = vec![role(map["default_patch"]["role"].as_str().unwrap_or("wall")); facets];
        for p in ports {
            let range = p["facet_range"].as_array();
            let (a, b) = match range.map(|r| {
                (
                    r.first().and_then(Value::as_u64),
                    r.get(1).and_then(Value::as_u64),
                )
            }) {
                Some((Some(a), Some(b))) if a < b && b as usize <= facets => {
                    (a as usize, b as usize)
                }
                _ => return Err(invalid("a port has no valid facet_range".into())),
            };
            let r = role(p["role"].as_str().unwrap_or("wall"));
            out[a..b].fill(r);
        }
        return Ok(out);
    }
    Err(invalid(
        "expected a patch map (patches, facet_patch) or a port file (ports)".into(),
    ))
}

fn tag(r: Role) -> usize {
    match r {
        Role::Wall => 0,
        Role::Inlet => 1,
        Role::Outlet => 2,
        Role::Symmetry => 3,
        Role::Cut => 4,
    }
}

fn role_of(t: usize) -> Role {
    [
        Role::Wall,
        Role::Inlet,
        Role::Outlet,
        Role::Symmetry,
        Role::Cut,
    ]
    .get(t)
    .copied()
    .unwrap_or(Role::Wall)
}

fn share(a: f64, b: f64) -> f64 {
    if b > 0.0 {
        ((a / b) * 1e4).round() / 1e4
    } else {
        0.0
    }
}

fn v3(p: Vec3) -> [f64; 3] {
    let r = |x: f64| (x * 1e6).round() / 1e6;
    [r(p.x), r(p.y), r(p.z)]
}

fn shares(c: &Classification) -> Value {
    let v = c.kind_volumes();
    let mut m = serde_json::Map::new();
    for (k, kind) in Kind::ALL.iter().enumerate() {
        m.insert(kind.name().into(), json!(share(v[k], c.volume)));
    }
    Value::Object(m)
}

fn advice(kind: Kind) -> &'static str {
    match kind {
        Kind::Sweep => {
            "swept block template along the skeleton (O-grid in a round section, H in a rectangular one)"
        }
        Kind::Junction => "junction template, or the polyhedral mesher with imprint",
        Kind::Obstacle => "O-grid about the body, joined to H blocks",
        Kind::Box => "box template (chamber), cut at its ports",
        Kind::Remainder => "polyhedral mesher with imprint",
    }
}

/// Run the proposal on an STL and the role of each of its facets.
pub fn propose(stl: &[u8], map: &Value, opts: RegionOptions) -> Result<Proposal, EngineError> {
    let start = Instant::now();
    let surface = mesh_geometry::import_closed_stl(stl, mesh_geometry::StlOptions::default())
        .map_err(|e| bad(ErrorCode::InvalidDomain, e))?;
    let sha = crate::artifact::sha256(stl);
    let facets = surface.facets().len();
    let roles = facet_roles(map, facets, &sha)?;
    if !roles.contains(&Role::Inlet) || !roles.contains(&Role::Outlet) {
        return Err(bad(
            ErrorCode::InvalidPatchMap,
            "the map needs at least one inlet and one outlet patch (set them with patch.set_role)",
        ));
    }
    let volume = surface.volume();
    let h = match opts.size {
        Some(h) if h.is_finite() && h > 0.0 => h,
        Some(_) => return Err(bad(ErrorCode::InvalidArguments, "size must be positive")),
        None => (volume / (0.118 * DEFAULT_TETS)).cbrt(),
    };
    if !(opts.heat_m > 0.0 && opts.heat_m.is_finite()) {
        return Err(bad(ErrorCode::InvalidArguments, "heat_m must be positive"));
    }
    let triangles: Vec<[usize; 3]> = surface
        .facets()
        .iter()
        .map(|f| f.triangle.vertices)
        .collect();
    let groups: Vec<usize> = roles.iter().map(|&r| tag(r)).collect();
    let wall = remesh::remesh_grouped(
        surface.vertices(),
        &triangles,
        &vec![false; facets],
        &groups,
        remesh::RemeshOptions {
            // curved walls (tubes) and gaps get finer edges down to h / 4
            tolerance: 0.02 * h,
            min_size: 0.25 * h,
            max_size: h,
            rate: 0.4,
            feature_angle: 30f64.to_radians(),
            iterations: 10,
            rim_size: None,
            thickness: Some(remesh::ThicknessSizing {
                across: 1.0,
                min_size: 0.25 * h,
                wedge_angle: 45f64.to_radians(),
            }),
        },
        None,
    )
    .map_err(|e| bad(ErrorCode::GenerationFailed, format!("wall remesh: {e}")))?;
    let (sc, st) = scaffold::build(
        &wall.vertices,
        &wall.triangles,
        &wall.group,
        scaffold::ScaffoldOptions {
            max_size: h,
            rate: 0.4,
            max_points: MAX_POINTS,
            ..Default::default()
        },
    )
    .map_err(|e| bad(ErrorCode::GenerationFailed, format!("scaffold: {e}")))?;
    let scaffold_seconds = start.elapsed().as_secs_f64();
    let mesh = TetMesh {
        points: sc.points.clone(),
        tets: sc.tetrahedra.clone(),
        faces: sc.boundary.iter().map(|&(f, t)| (f, role_of(t))).collect(),
    };
    let fem = Fem::new(&mesh).map_err(|e| bad(ErrorCode::GenerationFailed, e))?;
    let adjacent = fields::neighbours(&mesh.tets);
    let t0 = Instant::now();
    let (d, ds) = fields::heat_distance(&mesh, &fem, opts.heat_m)
        .map_err(|e| bad(ErrorCode::InvalidPatchMap, e))?;
    let (p, ps) = fields::potential(&mesh, &fem).map_err(|e| bad(ErrorCode::InvalidPatchMap, e))?;
    let field_seconds = t0.elapsed().as_secs_f64();
    if !(ds.converged && ps.converged) {
        return Err(bad(
            ErrorCode::GenerationFailed,
            format!(
                "the field solves did not converge (distance residual {:.1e}, potential {:.1e})",
                ds.residual, ps.residual
            ),
        ));
    }
    let t1 = Instant::now();
    let c = regions::classify(&mesh, &fem, &adjacent, &d, &p, opts.classify);
    let second = if opts.second_pass {
        regions::second_pass(
            &mesh,
            &adjacent,
            &d,
            &c,
            &[Kind::Box, Kind::Remainder],
            opts.classify,
        )
        .map_err(|e| bad(ErrorCode::GenerationFailed, e))?
    } else {
        Vec::new()
    };
    let classify_seconds = t1.elapsed().as_secs_f64();
    // second-pass kind of each tetrahedron of a region it treated
    let mut second_kind: Vec<i32> = vec![-1; mesh.tets.len()];
    let mut second_json = Vec::new();
    let mut gained = 0.0;
    let mut second_volume = 0.0;
    for sp in &second {
        let mut k = 0;
        for (t, &r) in c.tet_region.iter().enumerate() {
            if r == sp.region {
                let sr = sp.result.tet_region[k];
                second_kind[t] = kind_index(sp.result.regions[sr].kind);
                k += 1;
            }
        }
        let v = sp.result.kind_volumes();
        let templ = v[0] + v[2] + v[3]; // sweep, obstacle, box
        let before = c.regions[sp.region].kind;
        let won = if before == Kind::Box {
            v[0] + v[2]
        } else {
            templ
        };
        gained += won;
        second_volume += sp.result.volume;
        second_json.push(json!({
            "region": sp.region,
            "kind_before": before.name(),
            "volume_share_of_domain": share(sp.result.volume, c.volume),
            "shares": shares(&sp.result),
            "blockable_share": share(sp.result.blockable_volume, sp.result.volume),
            "template_share": share(won, sp.result.volume),
            "branch_points": sp.result.branch_points(),
            "end_points": sp.result.end_points(),
            "fiedler_iterations": sp.solve.iterations,
        }));
    }
    let kv = c.kind_volumes();
    let template_volume = kv[0] + kv[1] + kv[2] + kv[3];
    let pieces = &c.pieces;
    let mut piece_region = vec![0; pieces.len()];
    for (i, r) in c.regions.iter().enumerate() {
        for &q in &r.pieces {
            piece_region[q] = i;
        }
    }
    let report = json!({
        "format": REGION_FORMAT,
        "format_version": "1.0",
        "input": {
            "surface_sha256": sha,
            "facets": facets,
            "inlet_facets": roles.iter().filter(|&&r| r == Role::Inlet).count(),
            "outlet_facets": roles.iter().filter(|&&r| r == Role::Outlet).count(),
            "symmetry_facets": roles.iter().filter(|&&r| r == Role::Symmetry).count(),
            "size": h,
            "heat_m": opts.heat_m,
        },
        "scaffold": {
            "nodes": fem.used.iter().filter(|&&u| u).count(),
            "tetrahedra": mesh.tets.len(),
            "wall_triangles": wall.triangles.len(),
            "mean_edge": fem.h,
            "volume": c.volume,
            "stl_volume": volume,
            "recovery_points": st.recovery_points,
        },
        "fields": {
            "heat_time": (opts.heat_m * fem.h).powi(2),
            "distance": {"iterations": ds.iterations, "residual": ds.residual,
                         "max": d.iter().cloned().fold(0.0, f64::max)},
            "potential": {"iterations": ps.iterations, "residual": ps.residual},
            "streamwise_length": c.length,
        },
        "blockable_share": share(c.blockable_volume, c.volume),
        "shares": shares(&c),
        "template_share": share(template_volume, c.volume),
        "skeleton": {
            "slabs": c.slabs,
            "pieces": pieces.len(),
            "end_points": c.end_points(),
            "branch_points": c.branch_points(),
            "loops": c.loops,
            "ridge_faces": c.ridge_faces,
            "nodes": pieces.iter().zip(&piece_region).map(|(q, r)| json!({
                "at": v3(q.position()),
                "kind": if q.is_junction() { "branch" } else if q.is_end() { "end" } else { "regular" },
                "region": r,
            })).collect::<Vec<_>>(),
            "edges": c.edges,
        },
        "regions": c.regions.iter().enumerate().map(|(i, r)| json!({
            "id": i,
            "kind": r.kind.name(),
            "volume_share": share(r.volume, c.volume),
            "clash_share": share(r.clash_volume, r.volume),
            "slabs": [r.slabs.0, r.slabs.1],
            "pieces": r.pieces.len(),
            "centroid": v3(r.centroid),
            "proposal": advice(r.kind),
        })).collect::<Vec<_>>(),
        "second_pass": {
            "method": "each box and remainder region on its own: same wall distance, its Fiedler vector (long axis) in place of the flow potential. Its cut faces as new inlets and outlets would give the same potential as the first pass, because the cuts lie on level sets of the potential.",
            "regions": second_json,
            "volume_share_of_domain": share(second_volume, c.volume),
            "template_share_gained": share(gained, c.volume),
        },
        "summary": summary(&c, template_volume, gained),
        "seconds": {
            "scaffold": scaffold_seconds,
            "fields": field_seconds,
            "classify": classify_seconds,
            "total": start.elapsed().as_secs_f64(),
        },
    });
    let vtu = vtu(&mesh, &d, &p, &c, &second_kind);
    let svg = svg(&mesh, &c, &second_kind, &report);
    Ok(Proposal { report, vtu, svg })
}

fn kind_index(k: Kind) -> i32 {
    Kind::ALL
        .iter()
        .position(|&x| x == k)
        .map_or(-1, |i| i as i32)
}

fn summary(c: &Classification, template: f64, gained: f64) -> String {
    let v = c.kind_volumes();
    let pct = |x: f64| 100.0 * x / c.volume.max(f64::MIN_POSITIVE);
    let mut s = format!(
        "Blocks fit {:.0} % of the volume (sweep {:.0} %, junction {:.0} %, obstacle {:.0} %, box {:.0} %); {:.0} % stays for the polyhedral mesher.",
        pct(template),
        pct(v[0]),
        pct(v[1]),
        pct(v[2]),
        pct(v[3]),
        pct(v[4])
    );
    if gained > 0.0 {
        let _ = write!(
            s,
            " Second pass: {:.0} % of the volume inside box and remainder regions splits further into sweeps (or boxes).",
            pct(gained)
        );
    }
    s
}

const COLOURS: [&str; 5] = ["#2e7dd1", "#e0912b", "#8a5cc2", "#3a9e5c", "#c94a4a"];

fn vtu(mesh: &TetMesh, d: &[f64], p: &[f64], c: &Classification, second: &[i32]) -> String {
    let mut s = String::new();
    let n = mesh.points.len();
    let nt = mesh.tets.len();
    let _ = write!(
        s,
        "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n<Piece NumberOfPoints=\"{n}\" NumberOfCells=\"{nt}\">\n<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n"
    );
    for q in &mesh.points {
        let _ = writeln!(s, "{} {} {}", q.x, q.y, q.z);
    }
    s.push_str("</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for t in &mesh.tets {
        let _ = writeln!(s, "{} {} {} {}", t[0], t[1], t[2], t[3]);
    }
    s.push_str("</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for k in 0..nt {
        let _ = writeln!(s, "{}", 4 * (k + 1));
    }
    s.push_str("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for _ in 0..nt {
        s.push_str("10\n");
    }
    s.push_str("</DataArray>\n</Cells>\n<PointData>\n");
    for (name, f) in [("wall_distance", d), ("potential", p)] {
        let _ = writeln!(
            s,
            "<DataArray type=\"Float64\" Name=\"{name}\" format=\"ascii\">"
        );
        for v in f {
            let _ = writeln!(s, "{v}");
        }
        s.push_str("</DataArray>\n");
    }
    s.push_str(
        "</PointData>\n<CellData>\n<DataArray type=\"Float64\" Name=\"clash\" format=\"ascii\">\n",
    );
    for v in &c.clash {
        let _ = writeln!(s, "{v}");
    }
    s.push_str("</DataArray>\n");
    let ints: [(&str, Vec<i64>); 4] = [
        ("region", c.tet_region.iter().map(|&r| r as i64).collect()),
        (
            "kind",
            c.tet_region
                .iter()
                .map(|&r| kind_index(c.regions[r].kind) as i64)
                .collect(),
        ),
        ("piece", c.tet_piece.iter().map(|&r| r as i64).collect()),
        (
            "second_pass_kind",
            second.iter().map(|&k| k as i64).collect(),
        ),
    ];
    for (name, f) in ints {
        let _ = writeln!(
            s,
            "<DataArray type=\"Int64\" Name=\"{name}\" format=\"ascii\">"
        );
        for v in f {
            let _ = writeln!(s, "{v}");
        }
        s.push_str("</DataArray>\n");
    }
    s.push_str("</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    s
}

/// An orthographic picture along the shortest bounding-box axis: the wall
/// faces coloured by the region behind them (translucent, far first), the
/// skeleton on top with its end points (white) and branch points (black),
/// and a legend with the volume shares.
fn svg(mesh: &TetMesh, c: &Classification, second: &[i32], report: &Value) -> String {
    let (mut lo, mut hi) = (
        Vec3::new(f64::MAX, f64::MAX, f64::MAX),
        Vec3::new(f64::MIN, f64::MIN, f64::MIN),
    );
    for f in &mesh.faces {
        for &v in &f.0 {
            let q = mesh.points[v];
            lo = Vec3::new(lo.x.min(q.x), lo.y.min(q.y), lo.z.min(q.z));
            hi = Vec3::new(hi.x.max(q.x), hi.y.max(q.y), hi.z.max(q.z));
        }
    }
    let e = hi - lo;
    let ext = [e.x, e.y, e.z];
    let depth = (0..3)
        .min_by(|&a, &b| ext[a].total_cmp(&ext[b]))
        .unwrap_or(2);
    let (ax, ay) = match depth {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    };
    let comp = |q: Vec3, a: usize| [q.x, q.y, q.z][a];
    let (w, hgt) = (
        900.0,
        900.0 * (ext[ay] / ext[ax].max(f64::MIN_POSITIVE)).clamp(0.3, 1.6),
    );
    let (margin, legend) = (20.0, 110.0);
    let scale = (w - 2.0 * margin) / ext[ax].max(f64::MIN_POSITIVE);
    let px = |q: Vec3| {
        (
            margin + (comp(q, ax) - comp(lo, ax)) * scale,
            margin + (comp(hi, ay) - comp(q, ay)) * scale,
        )
    };
    // the tet behind each boundary face
    let mut owner = std::collections::HashMap::new();
    for (k, t) in mesh.tets.iter().enumerate() {
        for a in 0..4 {
            owner.insert(fields::face_key(t, a), k);
        }
    }
    let mut faces: Vec<(f64, usize, [usize; 3])> = mesh
        .faces
        .iter()
        .filter_map(|(f, _)| {
            let mut key = *f;
            key.sort_unstable();
            let k = *owner.get(&key)?;
            let z = f.iter().map(|&v| comp(mesh.points[v], depth)).sum::<f64>() / 3.0;
            Some((z, k, *f))
        })
        .collect();
    faces.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    let has_second = second.iter().any(|&k| k >= 0);
    let panel_h = hgt + 30.0;
    let panels = if has_second { 2.0 } else { 1.0 };
    let total_h = panels * panel_h + legend;
    let mut s = String::new();
    let _ = write!(
        s,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{total_h:.0}\" viewBox=\"0 0 {w} {total_h:.0}\" font-family=\"sans-serif\">\n<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n"
    );
    let polygons = |s: &mut String, dy: f64, colour: &dyn Fn(usize) -> (usize, f64)| {
        for (_, k, f) in &faces {
            let (kind, opacity) = colour(*k);
            let pts: Vec<String> = f
                .iter()
                .map(|&v| {
                    let (x, y) = px(mesh.points[v]);
                    format!("{x:.1},{:.1}", y + dy)
                })
                .collect();
            let _ = writeln!(
                s,
                "<polygon points=\"{}\" fill=\"{}\" fill-opacity=\"{opacity}\"/>",
                pts.join(" "),
                COLOURS[kind.min(4)]
            );
        }
    };
    // interior zones (junction, obstacle, remainder) as their projected
    // tetrahedra, so a clash zone inside the domain shows too
    let interior = |s: &mut String, dy: f64, kind_of: &dyn Fn(usize) -> Option<usize>| {
        let chosen: Vec<usize> = (0..mesh.tets.len())
            .filter(|&k| kind_of(k).is_some())
            .collect();
        let step = chosen.len().div_ceil(40_000).max(1);
        for &k in chosen.iter().step_by(step) {
            let Some(kind) = kind_of(k) else { continue };
            let q: Vec<(f64, f64)> = mesh.tets[k].iter().map(|&v| px(mesh.points[v])).collect();
            let hull = hull4(&q);
            let pts: Vec<String> = hull
                .iter()
                .map(|&(x, y)| format!("{x:.1},{:.1}", y + dy))
                .collect();
            let _ = writeln!(
                s,
                "<polygon points=\"{}\" fill=\"{}\" fill-opacity=\"0.35\"/>",
                pts.join(" "),
                COLOURS[kind.min(4)]
            );
        }
    };
    let first_kind = |k: usize| kind_index(c.regions[c.tet_region[k]].kind).max(0) as usize;
    let _ = writeln!(
        s,
        "<text x=\"{margin}\" y=\"18\" font-size=\"15\" font-weight=\"bold\">First pass: flow potential, skeleton (white: end point, black: branch point)</text>"
    );
    polygons(&mut s, 12.0, &|k| (first_kind(k), 0.18));
    interior(&mut s, 12.0, &|k| {
        let kind = first_kind(k);
        matches!(kind, 1 | 2 | 4).then_some(kind)
    });
    for &(a, b) in &c.edges {
        let (x1, y1) = px(c.pieces[a].position());
        let (x2, y2) = px(c.pieces[b].position());
        let _ = writeln!(
            s,
            "<line x1=\"{x1:.1}\" y1=\"{:.1}\" x2=\"{x2:.1}\" y2=\"{:.1}\" stroke=\"#202020\" stroke-width=\"1.6\"/>",
            y1 + 12.0,
            y2 + 12.0
        );
    }
    for q in &c.pieces {
        let (x, y) = px(q.position());
        let y = y + 12.0;
        if q.is_junction() {
            let _ = writeln!(
                s,
                "<circle cx=\"{x:.1}\" cy=\"{y:.1}\" r=\"5\" fill=\"#000000\"/>"
            );
        } else if q.is_end() {
            let _ = writeln!(
                s,
                "<circle cx=\"{x:.1}\" cy=\"{y:.1}\" r=\"5\" fill=\"#ffffff\" stroke=\"#000000\" stroke-width=\"1.5\"/>"
            );
        }
    }
    if has_second {
        let _ = writeln!(
            s,
            "<text x=\"{margin}\" y=\"{:.0}\" font-size=\"15\" font-weight=\"bold\">Second pass: each box and remainder region along its own long axis</text>",
            panel_h + 18.0
        );
        polygons(&mut s, panel_h + 12.0, &|k| {
            if second[k] >= 0 {
                (second[k] as usize, 0.3)
            } else {
                (first_kind(k), 0.06)
            }
        });
        interior(&mut s, panel_h + 12.0, &|k| {
            (second[k] >= 0).then_some(second[k] as usize)
        });
    }
    let y0 = panels * panel_h + 10.0;
    let shares = &report["shares"];
    for (i, kind) in Kind::ALL.iter().enumerate() {
        let x = margin + i as f64 * 175.0;
        let pct = 100.0 * shares[kind.name()].as_f64().unwrap_or(0.0);
        let _ = writeln!(
            s,
            "<rect x=\"{x}\" y=\"{y0}\" width=\"16\" height=\"16\" fill=\"{}\" fill-opacity=\"0.6\"/><text x=\"{}\" y=\"{}\" font-size=\"15\">{} {pct:.0} %</text>",
            COLOURS[i],
            x + 22.0,
            y0 + 13.0,
            kind.name()
        );
    }
    let sk = &report["skeleton"];
    let _ = writeln!(
        s,
        "<text x=\"{margin}\" y=\"{}\" font-size=\"14\">Skeleton: {} end points, {} branch points, {} loops. Shares are of the first pass.</text>",
        y0 + 40.0,
        sk["end_points"],
        sk["branch_points"],
        sk["loops"]
    );
    let summary = report["summary"].as_str().unwrap_or("");
    let (l1, l2) = match summary.find("; ") {
        Some(i) => (&summary[..=i], summary[i + 2..].trim()),
        None => (summary, ""),
    };
    for (i, line) in [l1, l2].iter().enumerate() {
        let _ = writeln!(
            s,
            "<text x=\"{margin}\" y=\"{}\" font-size=\"14\">{}</text>",
            y0 + 62.0 + 20.0 * i as f64,
            xml(line)
        );
    }
    s.push_str("</svg>\n");
    s
}

/// The convex hull of up to four points (gift wrapping).
fn hull4(q: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let Some(start) = (0..q.len()).min_by(|&a, &b| q[a].0.total_cmp(&q[b].0)) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    let mut p = start;
    loop {
        out.push(q[p]);
        let mut next = (p + 1) % q.len();
        for r in 0..q.len() {
            let cross =
                (q[next].0 - q[p].0) * (q[r].1 - q[p].1) - (q[next].1 - q[p].1) * (q[r].0 - q[p].0);
            if cross < 0.0 {
                next = r;
            }
        }
        p = next;
        if p == start || out.len() > q.len() {
            break;
        }
    }
    out
}

fn xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}
