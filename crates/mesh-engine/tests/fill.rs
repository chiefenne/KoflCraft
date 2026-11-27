//! Sweep fill (`regions.fill`) on the elbow (a square duct with a mitre
//! bend) and the Manifold (two round pipes on a box).

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

fn run(name: &str) -> (Value, mesh_io::handoff::HandoffMesh) {
    let (stl, map) = fixture(name);
    let f = fill(&stl, &map, FillOptions::default()).unwrap();
    if let Ok(dir) = std::env::var("KOFL_FILL_OUT") {
        std::fs::write(
            format!("{dir}/{name}-fill.json"),
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
    let (r, mesh) = run("elbow");
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
    assert!(num(&v["max_non_orthogonality_deg"]) < 70.0, "{v}");
    assert!(num(&v["min_scaled_jacobian"]) > 0.2, "{v}");
    assert_eq!(mesh.summary().inverted_faces, 0);
    assert!(num(&r["shares"]["swept_cells"]) > 0.99, "{}", r["shares"]);
}

/// Keep: in the Manifold both pipes are swept from their port to a calm
/// cut in front of the box: the cut is planar (every node in its plane, a
/// flat level set) and square to the pipe, and the box is left for the
/// polyhedral mesher. The cuts are the interfaces of the conformal join;
/// a looser calm test once put a cut 3° askew inside the bend's exit. The
/// region tests do not cut.
#[test]
fn manifold_pipes_end_in_calm_cuts() {
    let (r, mesh) = run("manifold");
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
    assert!(num(&r["shares"]["left_for_the_polyhedral_mesher"]) > 0.7);
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
