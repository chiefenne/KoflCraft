//! Port detection on bare closed STL shells (BC-03, ports) and the duct setup.

use mesh_engine::duct::{self, PatchMap as DuctPatches, SetupOptions};
use mesh_engine::patches::{PatchRole, PatchSurface};
use mesh_engine::ports::PortOptions;

const JACKET: &[u8] = include_bytes!("../../../test-fixtures/cooling-jacket/cooling-jacket.stl");
const MANIFOLD: &[u8] = include_bytes!("../../../test-fixtures/manifold/manifold.stl");
const MANIFOLD_PORTS: &str = include_str!("../../../test-fixtures/manifold/patches.json");

/// Keep rule (SW-18). Invariant: on the Manifold surface (90,894 facets,
/// rounded pipes, a box, many planar faces) the two best candidates are
/// exactly the inlet and outlet caps recorded by hand in the fixture's port
/// file, with circular rims at 90° to the pipe wall. The hand-made ranges
/// are an independent reference; a wrong planar flood, loop count or rim
/// test changes the sets or the ranking.
#[test]
fn detection_recovers_the_hand_recorded_manifold_ports() {
    let s = PatchSurface::import(MANIFOLD, None).unwrap();
    let d = s.detect_ports(PortOptions::default()).unwrap();
    let ports: serde_json::Value = serde_json::from_str(MANIFOLD_PORTS).unwrap();
    let mut want: Vec<Vec<usize>> = (0..2)
        .map(|i| {
            let r = &ports["ports"][i]["facet_range"];
            (r[0].as_u64().unwrap() as usize..r[1].as_u64().unwrap() as usize).collect()
        })
        .collect();
    let mut got: Vec<Vec<usize>> = d.candidates[..2]
        .iter()
        .map(|c| c.members.clone())
        .collect();
    want.sort();
    got.sort();
    assert_eq!(got, want);
    for c in &d.candidates[..2] {
        assert!(c.circularity > 0.999, "{c:?}");
        assert_eq!(c.pipe_end_rim, 1.0);
        assert!((c.equivalent_radius - 25.0).abs() < 0.01);
    }
    assert!(d.candidates[2].score < 0.9 * d.candidates[1].score);
}

/// Keep rule (SW-18). Invariant: on the synthetic cooling jacket (a ribbed
/// annular shell, two blocks, two bent pipes; genus 1) the two pipe end caps
/// are the two best candidates, 96-gons of radius 4 facing away from the
/// jacket, far ahead of the block faces; the written map passes BC-06 and
/// `generate-duct` reads it, and the setup sizes follow the port radius.
#[test]
fn detection_and_setup_on_the_cooling_jacket() {
    let s = PatchSurface::import(JACKET, None).unwrap();
    let d = s.detect_ports(PortOptions::default()).unwrap();
    let area = 0.5 * 96.0 * 16.0 * (2.0 * std::f64::consts::PI / 96.0).sin();
    let mut z = Vec::new();
    for c in &d.candidates[..2] {
        assert!((c.area - area).abs() < 1e-3 * area, "{c:?}");
        assert!(c.circularity > 0.999 && c.pipe_end_rim == 1.0);
        assert!((c.normal[0] - 0.5f64.sqrt()).abs() < 1e-6);
        z.push(c.normal[2].signum());
    }
    z.sort_by(f64::total_cmp);
    assert_eq!(z, [-1.0, 1.0]);
    assert!(d.candidates[2].score < 0.6);
    let map = s.port_map(&d, 0, 1).unwrap();
    let roles: Vec<_> = map.patches.iter().map(|p| p.role).collect();
    assert_eq!(
        roles,
        [PatchRole::Wall, PatchRole::Inlet, PatchRole::Outlet]
    );
    assert!(s.port_map(&d, 1, 1).is_err());
    assert!(s.port_map(&d, 0, 99).is_err());

    let ducts: DuctPatches = serde_json::from_value(serde_json::to_value(&map).unwrap()).unwrap();
    let (plan, r) = duct::setup(JACKET, &ducts, SetupOptions::default()).unwrap();
    assert!(
        (r / (area / std::f64::consts::PI).sqrt() - 1.0).abs() < 1e-5,
        "{r}"
    );

    let v = serde_json::to_value(&plan).unwrap();
    assert_eq!(v["regions"], serde_json::json!([]));
    // Keep: the tolerance is not r/500 = 0.008 here. The jacket's cylinders
    // are 96-sided, so its own facets are 0.035 off the round wall, and a
    // tolerance below that makes the remesher chase the facet lines and
    // fill the wall with needles (the recovery blowup on thin parts).
    assert_eq!(v["volume_scaffold"]["build"]["tolerance"], 0.0442);
    // Keep: thin gaps and thin walls set the wall size (one edge across).
    assert_eq!(v["volume_scaffold"]["build"]["thickness"]["across"], 1.0);
    assert_eq!(v["first_thickness"], 0.0128);
    let coarse = SetupOptions {
        size_factor: 2.0,
        ..SetupOptions::default()
    };
    let (plan2, _) = duct::setup(JACKET, &ducts, coarse).unwrap();
    let v2 = serde_json::to_value(&plan2).unwrap();
    assert_eq!(v2["volume_scaffold"]["build"]["max_size"], 1.92);
}

/// Keep: a patch map that is not one inlet, one outlet and walls is refused
/// by the duct input, with the reason.
#[test]
fn duct_input_refuses_maps_without_one_inlet_and_outlet() {
    let s = PatchSurface::import(JACKET, None).unwrap();
    let d = s.detect_ports(PortOptions::default()).unwrap();
    let mut map = s.port_map(&d, 0, 1).unwrap();
    map.patches[2].role = PatchRole::Inlet;
    let ducts: DuctPatches = serde_json::from_value(serde_json::to_value(&map).unwrap()).unwrap();
    let e = duct::setup(JACKET, &ducts, SetupOptions::default())
        .err()
        .unwrap();
    assert!(
        e.message.contains("exactly one inlet and one outlet"),
        "{}",
        e.message
    );
    map.patches[2].role = PatchRole::Symmetry;
    let ducts: DuctPatches = serde_json::from_value(serde_json::to_value(&map).unwrap()).unwrap();
    assert!(duct::setup(JACKET, &ducts, SetupOptions::default()).is_err());
}
