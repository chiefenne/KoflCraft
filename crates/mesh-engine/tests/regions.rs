//! The field-guided region proposal (`regions.propose`) on analytic ducts
//! and the elbow fixture.

#[path = "../../../test-fixtures/ducts.rs"]
mod ducts;

use mesh_engine::regions::{RegionOptions, propose};
use serde_json::Value;

fn run(stl: &[u8], map: &str, size: f64) -> Value {
    let map: Value = serde_json::from_str(map).unwrap();
    let p = propose(
        stl,
        &map,
        RegionOptions {
            size: Some(size),
            ..Default::default()
        },
    )
    .unwrap();
    if let Ok(dir) = std::env::var("KOFL_REGION_OUT") {
        let n = p.report["input"]["facets"].as_u64().unwrap_or(0);
        std::fs::write(format!("{dir}/report-{n}.json"), p.report.to_string()).unwrap();
        std::fs::write(format!("{dir}/regions-{n}.svg"), &p.svg).unwrap();
    }
    p.report
}

fn share(r: &Value, kind: &str) -> f64 {
    r["shares"][kind].as_f64().unwrap()
}

/// Keep: the T has one branch point at its junction and three ends, and
/// its arms are sweeps. A false merge of pieces or a lost slab link
/// changes the branch count; the straight-duct unit test cannot see it.
#[test]
fn tee_has_one_junction_and_three_ends() {
    let d = ducts::tee();
    let (stl, map) = ducts::files(&d);
    let r = run(&stl, &map, 0.15);
    let sk = &r["skeleton"];
    assert_eq!(sk["branch_points"], 1, "{sk}");
    assert_eq!(sk["end_points"], 3, "{sk}");
    assert_eq!(sk["loops"], 0);
    assert!(share(&r, "sweep") > 0.7, "{}", r["shares"]);
    assert!(share(&r, "junction") > 0.0);
}

/// Keep: the fork splits once into two arms: one branch point, three
/// ends, no loop.
#[test]
fn fork_has_one_junction() {
    let d = ducts::fork();
    let (stl, map) = ducts::files(&d);
    let r = run(&stl, &map, 0.15);
    let sk = &r["skeleton"];
    assert_eq!(sk["branch_points"], 1, "{sk}");
    assert_eq!(sk["end_points"], 3, "{sk}");
    assert!(share(&r, "sweep") > 0.6, "{}", r["shares"]);
}

/// Keep: the elbow (one bend, no branch) is one sweep with two ends; the
/// clash measure must not flag the bend.
#[test]
fn elbow_is_one_sweep() {
    let stl = std::fs::read("../../test-fixtures/elbow/elbow.stl").unwrap();
    let map = std::fs::read_to_string("../../test-fixtures/elbow/patches.json").unwrap();
    let r = run(&stl, &map, 0.12);
    assert_eq!(r["skeleton"]["branch_points"], 0);
    assert_eq!(r["skeleton"]["end_points"], 2);
    assert!(share(&r, "sweep") > 0.95, "{}", r["shares"]);
    assert!(r["blockable_share"].as_f64().unwrap() > 0.95);
}

/// Keep: a map without an outlet is refused with a hint, not run with a
/// flat potential.
#[test]
fn map_without_outlet_is_refused() {
    let d = ducts::tee();
    let (stl, map) = ducts::files(&d);
    let map = map.replace("\"outlet\"}", "\"wall\"}");
    let map: Value = serde_json::from_str(&map).unwrap();
    let e = propose(&stl, &map, RegionOptions::default()).err().unwrap();
    assert!(e.message.contains("outlet"), "{}", e.message);
}

/// Keep: the Manifold's chamber is a box between two sweeps, and the
/// second pass finds sweeps inside it along its long axis (the reason the
/// second pass exists).
#[test]
#[ignore = "full run: about 20 s in a release build"]
fn manifold_chamber_is_a_box_the_second_pass_splits() {
    let stl = std::fs::read("../../test-fixtures/manifold/manifold.stl").unwrap();
    let map = std::fs::read_to_string("../../test-fixtures/manifold/patches.json").unwrap();
    let map: Value = serde_json::from_str(&map).unwrap();
    let r = propose(&stl, &map, RegionOptions::default())
        .unwrap()
        .report;
    assert!(share(&r, "box") > 0.5, "{}", r["shares"]);
    assert!(share(&r, "sweep") > 0.2, "{}", r["shares"]);
    assert_eq!(r["skeleton"]["branch_points"], 0);
    let gained = r["second_pass"]["template_share_gained"].as_f64().unwrap();
    assert!(gained > 0.3, "{gained}");
}
