//! Datum objects and coordinate regions on the nonconvex duct path (GRD-02,
//! GRD-03, GRD-08).

use mesh_engine::MeshDirectory;
use mesh_engine::duct::{self, Mode, PatchMap as DuctPatches, SetupOptions};
use mesh_engine::patches::{PatchRole, PatchSurface, set_role};
use mesh_engine::proposals::ProposalOptions;

#[path = "../../../test-fixtures/shapes.rs"]
mod shapes;

fn pipe() -> (Vec<u8>, DuctPatches) {
    let stl = shapes::ascii_stl(&shapes::strip_cylinder(0.0, 0.0, 5.0, 40.0, 48));
    let s = PatchSurface::import(&stl, None).unwrap();
    let p = s.propose_patches(ProposalOptions::default()).unwrap();
    let (m, _) = set_role(&p.map, "opening_1", PatchRole::Inlet).unwrap();
    let (m, _) = set_role(&m, "opening_2", PatchRole::Outlet).unwrap();
    let ducts = serde_json::from_value(serde_json::to_value(&m).unwrap()).unwrap();
    (stl, ducts)
}

fn options() -> SetupOptions {
    SetupOptions {
        size_factor: 2.0,
        ..SetupOptions::default()
    }
}

fn run(
    stl: &[u8],
    ducts: &DuctPatches,
    grading: Option<serde_json::Value>,
    tag: &str,
) -> serde_json::Value {
    let (plan, _) = duct::setup(stl, ducts, options()).unwrap();
    let mut v = serde_json::to_value(&plan).unwrap();
    // coarse and without thickness sizing: the test is about the grading
    v["volume_scaffold"]["build"]["thickness"] = serde_json::Value::Null;
    if let Some(g) = grading {
        v["volume_scaffold"]["build"]["grading"] = g;
    }
    let plan: duct::Partition = serde_json::from_value(v.clone()).unwrap();
    let dir = std::env::temp_dir().join(format!("kofl-grading-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let d = MeshDirectory::open(&dir).unwrap();
    let r = duct::generate(
        stl,
        ducts,
        &plan,
        Mode::Polyhedral,
        None,
        &d,
        std::ffi::OsStr::new("m"),
        mesh_engine::handoff::HandoffFormat::Vtu,
    )
    .unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    r
}

/// Keep rule (SW-18). Invariant (GRD-02, GRD-03, GRD-08): on the duct
/// path a datum sphere and a coordinate box refine the core where they say:
/// the cells in their zones reach their target sizes (median cbrt(V) within
/// 25 %), well below the other core cells, and the report names each source
/// with its achieved sizes. A grading that does not reach the wall remesh
/// or the scaffold refinement leaves the zones at the background size; the
/// convex-path tests do not run the duct mesher.
#[test]
fn datum_and_region_refine_the_duct_core() {
    let (stl, ducts) = pipe();
    let g = serde_json::json!({"sources": [
        {"kind": "datum", "name": "probe", "shape": {"type": "sphere", "centre": [0, 0, 12], "radius": 1.0},
         "size": 0.25, "distance": 1.0},
        {"kind": "region", "name": "wake", "shape": {"type": "box", "min": [-2, -2, 24], "max": [2, 2, 30]},
         "size": 0.3}
    ]});
    let r = run(&stl, &ducts, Some(g), "graded");
    let grading = &r["grading"];
    let other = grading["other_core_cells"]["median"].as_f64().unwrap();
    for (k, (name, target)) in [("probe", 0.25), ("wake", 0.3)].into_iter().enumerate() {
        let s = &grading["sources"][k];
        assert_eq!(s["name"], name);
        assert!(s["cells"].as_u64().unwrap() > 50, "{s}");
        let median = s["achieved"]["median"].as_f64().unwrap();
        assert!(
            (median / target - 1.0).abs() < 0.25,
            "{name}: {median} vs {target}"
        );
        assert!(median < 0.75 * other, "{name}: {median} vs other {other}");
    }
    // without grading the report has none
    let (plan, _) = duct::setup(&stl, &ducts, options()).unwrap();
    let v = serde_json::to_value(&plan).unwrap();
    assert!(v["volume_scaffold"]["build"].get("grading").is_none());
}

/// A derived (convex-surface) source is refused on the duct path, with the reason.
#[test]
fn derived_sources_are_refused_on_the_duct_path() {
    let (stl, ducts) = pipe();
    let (plan, _) = duct::setup(&stl, &ducts, options()).unwrap();
    let mut v = serde_json::to_value(&plan).unwrap();
    v["volume_scaffold"]["build"]["grading"] =
        serde_json::json!({"sources": [{"kind": "gap", "name": "g", "cells_across": 3}]});
    let plan: duct::Partition = serde_json::from_value(v).unwrap();
    let dir = std::env::temp_dir().join(format!("kofl-grading-{}-derived", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let d = MeshDirectory::open(&dir).unwrap();
    let e = duct::generate(
        &stl,
        &ducts,
        &plan,
        Mode::Polyhedral,
        None,
        &d,
        std::ffi::OsStr::new("m"),
        mesh_engine::handoff::HandoffFormat::Vtu,
    )
    .unwrap_err();
    assert!(e.message.contains("only datum and region"), "{e}");
    let _ = std::fs::remove_dir_all(&dir);
}
