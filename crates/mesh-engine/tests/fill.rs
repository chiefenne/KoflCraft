//! Sweep fill (`regions.fill`) on the elbow (a square duct with a mitre
//! bend), the Manifold (two round pipes on a box) and a box with three round
//! pipes, alone and joined to the polyhedral mesh of the box.

#[path = "../../../test-fixtures/ducts.rs"]
mod ducts;

use mesh_core::geometry::Vec3;
use mesh_engine::fill::{FillOptions, fill};
use serde_json::Value;

fn fixture(name: &str) -> (Vec<u8>, Value) {
    let dir = format!("{}/../../test-fixtures/{name}", env!("CARGO_MANIFEST_DIR"));
    let stl = std::fs::read(format!("{dir}/{name}.stl")).unwrap();
    let map: Value =
        serde_json::from_slice(&std::fs::read(format!("{dir}/patches.json")).unwrap()).unwrap();
    (stl, map)
}

fn run(name: &str, join: bool) -> (Value, mesh_io::handoff::HandoffMesh) {
    run_with(
        name,
        FillOptions {
            join,
            ..FillOptions::default()
        },
    )
}

fn run_with(name: &str, opts: FillOptions) -> (Value, mesh_io::handoff::HandoffMesh) {
    let (stl, map) = fixture(name);
    let join = opts.join;
    let f = fill(&stl, &map, opts).unwrap();
    if let Ok(dir) = std::env::var("KOFL_FILL_OUT") {
        std::fs::write(
            format!("{dir}/{name}-fill-{join}.json"),
            serde_json::to_string_pretty(&f.report).unwrap(),
        )
        .unwrap();
    }
    (f.report, f.mesh.expect("swept cells"))
}

fn num(v: &Value) -> f64 {
    v.as_f64().unwrap_or_else(|| panic!("not a number: {v}"))
}

/// Keep: the elbow is one sweep from port to port through the mitre, and
/// the station at the mitre is the mitre plane through both corner edges,
/// so the swept cells fill the duct to within 0.01 % and none is inverted.
/// Without the pinned, snapped station the wall faces cut the outer corner
/// (-0.5 % volume); with free Newton steps onto the level sets the nodes at
/// the outer corner jump and cells invert. The region tests never mesh.
#[test]
fn elbow_sweeps_through_the_mitre() {
    let (r, mesh) = run("elbow", true);
    let s = &r["sweeps"][0];
    assert_eq!(r["sweeps"].as_array().unwrap().len(), 1, "{}", r["sweeps"]);
    assert_eq!(r["not_swept"], Value::Array(vec![]));
    assert_eq!(s["ends"][0]["kind"], "port");
    assert_eq!(s["ends"][1]["kind"], "port");
    assert_eq!(s["section"]["corners"], 4, "{}", s["section"]);
    assert_eq!(s["feature_stations"], 1, "{s}");
    assert_eq!(s["missed_feature_edges"], 0);
    assert_eq!(s["unlifted_nodes"], 0);
    let dv = num(&s["relative_volume_difference"]).abs();
    assert!(dv < 1e-4, "volume off by {dv}");
    let v = &r["validity"];
    assert_eq!(v["inverted_faces"], 0);
    assert!(num(&v["min_volume"]) > 0.0);
    // the stations next to the mitre are blended planes on both sides
    assert!(num(&s["blended_stations"]) >= 4.0, "{s}");
    assert_eq!(s["unblended_mitre_sides"], 0, "{s}");
    assert!(num(&v["max_non_orthogonality_deg"]) < 52.0, "{v}");
    assert!(num(&v["min_scaled_jacobian"]) > 0.45, "{v}");
    assert!(
        s["section"]["kind"].as_str().unwrap().contains("H-grid"),
        "{}",
        s["section"]
    );
    assert_eq!(mesh.summary().inverted_faces, 0);
    assert!(num(&r["shares"]["swept_cells"]) > 0.99, "{}", r["shares"]);
}

/// Keep: the blended mitre stations (the classical structured mitre) beat
/// the level sets there (ALG-16). Without the blend the level sets crowd at
/// the inner corner: non-orthogonality 66° and min scaled Jacobian 0.26
/// against 59° and 0.37 with it. The test above only bounds the blended
/// mesh; this one fails when the blend stops helping or stops running.
#[test]
fn mitre_blend_beats_the_level_sets() {
    let opts = |mitre_blend| FillOptions {
        mitre_blend,
        ..FillOptions::default()
    };
    let (blend, _) = run_with("elbow", opts(true));
    let (level, _) = run_with("elbow", opts(false));
    assert_eq!(level["sweeps"][0]["blended_stations"], 0);
    let no = |r: &Value| num(&r["validity"]["max_non_orthogonality_deg"]);
    let sj = |r: &Value| num(&r["validity"]["min_scaled_jacobian"]);
    assert!(
        no(&blend) + 5.0 < no(&level),
        "{} vs {}",
        no(&blend),
        no(&level)
    );
    assert!(
        sj(&blend) > sj(&level) + 0.05,
        "{} vs {}",
        sj(&blend),
        sj(&level)
    );
    // the blend keeps the volume (the stations stay on the walls)
    let dv = |r: &Value| num(&r["sweeps"][0]["relative_volume_difference"]).abs();
    assert!(dv(&blend) < 1e-4, "{}", dv(&blend));
}

/// Keep: the H-grid of a rectangular section (the wall layers, then a
/// tensor core clustered to all four walls) beats its O-grid, whose columns
/// fan towards a smaller core and meet the corners at 45°: on the elbow the
/// worst non-orthogonality drops from 59.1° to 49.5° (now at the mitre;
/// 26.6° in the straight legs, at the corner cells of the layers) and the
/// smallest scaled Jacobian rises from 0.37 to 0.49. The elbow test above
/// only bounds the H-grid.
#[test]
fn h_grid_beats_the_o_grid() {
    let opts = |h_grid| FillOptions {
        h_grid,
        ..FillOptions::default()
    };
    let (h, _) = run_with("elbow", opts(true));
    let (o, _) = run_with("elbow", opts(false));
    let no = |r: &Value| num(&r["validity"]["max_non_orthogonality_deg"]);
    let sj = |r: &Value| num(&r["validity"]["min_scaled_jacobian"]);
    assert!(no(&h) + 5.0 < no(&o), "{} vs {}", no(&h), no(&o));
    assert!(sj(&h) > sj(&o) + 0.05, "{} vs {}", sj(&h), sj(&o));
    let dv = |r: &Value| num(&r["sweeps"][0]["relative_volume_difference"]).abs();
    assert!(dv(&h) < 1e-4, "{}", dv(&h));
}

/// Keep: the square T and the fork join to the box mesh at the default
/// five wall layers. Two fixes behind it each fail here first. (1) An
/// H-grid cut slides its wall nodes to where the core's lines, scaled out
/// by the polyhedral mesher's similarity, meet the wall, so the cap is
/// Delaunay; before, the scaled core leaned by up to the layer stack
/// against the rim at the corners (48 edges per cap off by up to 83°) and
/// the scaffold's boundary recovery split without end. (2) The wall layers
/// next to a cap take the cap's sideways shift, fading out over four layer
/// stacks; before, the rim of a square cap slid sideways while the wall
/// next to it only moved inward, thin tetrahedra turned over, the layers
/// yielded and a cap quad folded (a layer front collision with three to
/// five layers; two layers passed only by luck, not on CI).
#[test]
fn square_tee_and_fork_join() {
    for d in [ducts::tee(), ducts::fork()] {
        let (stl, map) = ducts::files(&d);
        let map: Value = serde_json::from_str(&map).unwrap();
        let f = fill(&stl, &map, FillOptions::default()).unwrap();
        let r = &f.report;
        let rem = &r["remainder"];
        assert!(rem.get("error").is_none(), "{rem}");
        assert_eq!(rem["surface"]["caps"], 3, "{rem}");
        let v = &r["validity"];
        assert_eq!(v["inverted_faces"], 0, "{v}");
        assert!(num(&v["max_relative_closure"]) < 1e-9, "{v}");
        let names: Vec<&str> = v["patches"]
            .as_array()
            .unwrap()
            .iter()
            .map(|p| p["name"].as_str().unwrap())
            .collect();
        assert_eq!(names, ["walls", "inlet", "outlet"], "the cuts are interior");
        assert!(
            num(&v["non_orthogonality_by_faces"]["at_the_cuts"]["max_deg"]) < 50.0,
            "{v}"
        );
    }
}

/// Keep: a pipe between two boxes is swept from cut to cut and the two
/// boxes behind the cuts are meshed as two shells and joined. Two fixes
/// each fail here first: (1) every cut's band is marked before the sweep
/// side is flooded, else the flood from one end of the pipe runs through
/// it past its other cut and drops the whole STL (0 holes left); (2) each
/// closed shell behind the cuts gets its own polyhedral mesh (the mesher
/// takes one connected shell). The other joined fixtures have one shell
/// and sweeps with one cut each. (Parallel square tubes between headers,
/// `ducts::header_tubes`, pass the same path on one machine and collide in
/// the layer front at a square cut corner on another: the square-corner
/// gap of the join, not this path.)
#[test]
fn pipe_between_two_boxes_joins_both() {
    let d = ducts::box_pipe_box();
    let (stl, map) = ducts::files(&d);
    let map: Value = serde_json::from_str(&map).unwrap();
    let f = fill(&stl, &map, FillOptions::default()).unwrap();
    let r = &f.report;
    let rem = &r["remainder"];
    assert!(rem.get("error").is_none(), "{rem}");
    let kinds: Vec<&Value> = r["regions"]
        .as_array()
        .unwrap()
        .iter()
        .map(|g| &g["kind"])
        .collect();
    assert_eq!(
        rem["shells"], 2,
        "{rem}\nnot swept: {}\nregions: {kinds:?}\nsweeps: {}",
        r["not_swept"], r["sweeps"]
    );
    let cut_to_cut = r["sweeps"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|s| s["ends"][0]["kind"] == "cut" && s["ends"][1]["kind"] == "cut")
        .count();
    assert_eq!(cut_to_cut, 1, "{}", r["sweeps"]);
    let v = &r["validity"];
    assert_eq!(v["inverted_faces"], 0, "{v}");
    assert!(num(&v["max_relative_closure"]) < 1e-9, "{v}");
    let names: Vec<&str> = v["patches"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| p["name"].as_str().unwrap())
        .collect();
    assert_eq!(names, ["walls", "inlet", "outlet"], "the cuts are interior");
}

/// Keep: in the Manifold both pipes are swept from their port to a calm
/// cut in front of the box: the cut is planar (every node in its plane, a
/// flat level set) and square to the pipe; without the join the box is
/// left out and the cuts are patches. The cuts are the interfaces of the
/// conformal join; a looser calm test once put a cut 3° askew inside the
/// bend's exit. The region tests do not cut.
#[test]
fn manifold_pipes_end_in_calm_cuts() {
    let (r, mesh) = run("manifold", false);
    let sweeps = r["sweeps"].as_array().unwrap();
    assert_eq!(sweeps.len(), 2, "{}", r["sweeps"]);
    for s in sweeps {
        assert_eq!(s["section"]["corners"], 0, "round pipe: {}", s["section"]);
        let kinds: Vec<&str> = (0..2)
            .map(|e| s["ends"][e]["kind"].as_str().unwrap())
            .collect();
        assert!(
            kinds.contains(&"port") && kinds.contains(&"cut"),
            "{kinds:?}"
        );
        for e in 0..2 {
            let end = &s["ends"][e];
            assert_eq!(end["calm"], true, "{end}");
            assert!(num(&end["flatness"]) <= 1e-3, "{end}");
        }
        let dv = num(&s["relative_volume_difference"]).abs();
        assert!(dv < 5e-3, "volume off by {dv}");
        assert!(num(&s["min_scaled_jacobian"]) > 0.5, "{s}");
    }
    assert!(num(&r["shares"]["not_meshed"]) > 0.7, "{}", r["shares"]);
    let v = &r["validity"];
    assert_eq!(v["inverted_faces"], 0);
    assert!(num(&v["max_non_orthogonality_deg"]) < 40.0, "{v}");
    // the cut patches are planar (the join's interface)
    let names: Vec<String> = v["patches"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| p["name"].as_str().unwrap().to_string())
        .collect();
    let cuts: Vec<&String> = names.iter().filter(|n| n.starts_with("cut_")).collect();
    assert_eq!(cuts.len(), 2, "{names:?}");
    for (pi, patch) in mesh.patches.iter().enumerate() {
        if !patch.name.starts_with("cut_") {
            continue;
        }
        let faces: Vec<&[usize]> = mesh.patch_faces(pi).map(|f| mesh.face(f)).collect();
        let mut n = Vec3::ZERO;
        let mut pts = Vec::new();
        for f in &faces {
            for k in 0..f.len() {
                n = n + mesh.points[f[k]].cross(mesh.points[f[(k + 1) % f.len()]]);
                pts.push(mesh.points[f[k]]);
            }
        }
        let area = 0.5 * n.length();
        let n = n / n.length();
        let c = pts.iter().fold(Vec3::ZERO, |s, &p| s + p) / pts.len() as f64;
        let off = pts
            .iter()
            .map(|&p| (p - c).dot(n).abs())
            .fold(0.0, f64::max);
        let diameter = 2.0 * (area / std::f64::consts::PI).sqrt();
        assert!(
            off <= 1e-9 * diameter,
            "{}: off its plane by {off}",
            patch.name
        );
    }
}

/// Keep: with the join the Manifold is one mesh: the box is polyhedra
/// (median dual with the sweeps' layer stack) and meets each pipe at its
/// cut in the same points, every hex face of a cut split into the
/// polyhedral cap faces it holds. The mesh is closed and valid as a whole
/// (each face shared by two cells with opposite orientation), keeps the
/// STL's volume and has no bad face at the cuts. Each of the three fixes
/// behind it fails here first: a band that takes the other pipe (more
/// holes than cuts), cap quads split about their centres (the scaffold's
/// boundary recovery splits without end) and a wall offset that only lets
/// a tetrahedron's own wall vertices yield (a tetrahedron at a box corner
/// turns over). The sweep test above runs without the join.
#[test]
fn manifold_box_joins_the_pipes() {
    let (r, mesh) = run("manifold", true);
    let rem = &r["remainder"];
    assert!(rem.get("error").is_none(), "{rem}");
    let v = &r["validity"];
    assert_eq!(v["inverted_faces"], 0, "{v}");
    assert!(num(&v["min_volume"]) > 0.0);
    assert!(num(&v["max_relative_closure"]) < 1e-9, "{v}");
    let names: Vec<&str> = v["patches"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| p["name"].as_str().unwrap())
        .collect();
    assert_eq!(names, ["walls", "inlet", "outlet"], "the cuts are interior");
    let by = &v["non_orthogonality_by_faces"];
    assert!(num(&by["at_the_cuts"]["max_deg"]) < 60.0, "{by}");
    assert_eq!(by["at_the_cuts"]["over_70_deg"], 0, "{by}");
    assert!(num(&by["between_hexahedra"]["max_deg"]) < 40.0, "{by}");
    // the cuts weld exactly and the cap faces tile the hex faces
    let d = 50.0;
    assert!(num(&rem["max_weld_distance"]) < 1e-9 * d, "{rem}");
    assert!(num(&rem["max_relative_area_error"]) < 1e-9, "{rem}");
    assert!(num(&rem["cap_faces_joined"]) > 0.0);
    let shares = &r["shares"];
    assert!(num(&shares["polyhedral_cells"]) > 0.7, "{shares}");
    assert!(num(&shares["not_meshed"]) < 0.005, "{shares}");
    // the whole mesh keeps the STL's volume (faceting of the round pipes)
    let dv = num(&v["total_volume"]) / num(&r["scaffold"]["stl_volume"]) - 1.0;
    assert!(dv.abs() < 5e-3, "volume off by {dv}");
    assert_eq!(mesh.summary().inverted_faces, 0);
}

/// Keep: behind the cuts of the three-pipe box the region has three
/// openings, all cuts (every pipe is swept from its port). Two cuts share a
/// cap id of the polyhedral mesher, which gives each connected cap its own
/// similarity, and the join tells them apart by their planes; the mesh is
/// one, closed and valid, with no cut left as a patch. Before, the fill
/// refused past two openings; with one similarity per cap id the layer
/// stack of the shared caps ran towards their mean centre; and the +x cut,
/// a level set 0.003 across a ring of pipe nodes, gave a strip of needles
/// that the scaffold's boundary recovery split without end. The Manifold
/// test above has two openings.
#[test]
fn pipe_box_joins_three_openings() {
    let (stl, map) = ducts::files(&ducts::pipe_box());
    let map: Value = serde_json::from_str(&map).unwrap();
    let f = fill(&stl, &map, FillOptions::default()).unwrap();
    let r = &f.report;
    assert_eq!(r["sweeps"].as_array().unwrap().len(), 3, "{}", r["sweeps"]);
    assert_eq!(r["not_swept"], Value::Array(vec![]));
    let rem = &r["remainder"];
    assert!(rem.get("error").is_none(), "{rem}");
    assert_eq!(rem["surface"]["caps"], 3, "{rem}");
    assert!(num(&rem["max_weld_distance"]) < 1e-12, "{rem}");
    assert!(num(&rem["max_relative_area_error"]) < 1e-9, "{rem}");
    let v = &r["validity"];
    assert_eq!(v["inverted_faces"], 0, "{v}");
    assert!(num(&v["min_volume"]) > 0.0);
    assert!(num(&v["max_relative_closure"]) < 1e-9, "{v}");
    let names: Vec<&str> = v["patches"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| p["name"].as_str().unwrap())
        .collect();
    assert_eq!(names, ["walls", "inlet", "outlet"], "the cuts are interior");
    let by = &v["non_orthogonality_by_faces"];
    assert!(num(&by["at_the_cuts"]["max_deg"]) < 60.0, "{by}");
    assert!(num(&r["shares"]["not_meshed"]) < 0.005, "{}", r["shares"]);
    assert_eq!(f.mesh.expect("cells").summary().inverted_faces, 0);
}
