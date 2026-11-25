//! Zero-thickness shells in the duct mesh (SHL-01 to
//! SHL-03).

use mesh_engine::MeshDirectory;
use mesh_engine::duct::{self, Mode, PatchMap as DuctPatches, SetupOptions};
use mesh_engine::patches::{PatchRole, PatchSurface, set_role};
use mesh_engine::proposals::ProposalOptions;
use serde_json::json;

#[path = "../../../test-fixtures/shapes.rs"]
mod shapes;

fn pipe_with_shells(core: &str) -> serde_json::Value {
    let stl = shapes::ascii_stl(&shapes::strip_cylinder(0.0, 0.0, 5.0, 40.0, 48));
    let s = PatchSurface::import(&stl, None).unwrap();
    let p = s.propose_patches(ProposalOptions::default()).unwrap();
    let (m, _) = set_role(&p.map, "opening_1", PatchRole::Inlet).unwrap();
    let (m, _) = set_role(&m, "opening_2", PatchRole::Outlet).unwrap();
    let ducts: DuctPatches = serde_json::from_value(serde_json::to_value(&m).unwrap()).unwrap();
    let options = SetupOptions {
        size_factor: 2.0,
        ..SetupOptions::default()
    };
    let (plan, _) = duct::setup(&stl, &ducts, options).unwrap();
    let mut v = serde_json::to_value(&plan).unwrap();
    v["volume_scaffold"]["build"]["thickness"] = serde_json::Value::Null;
    v["volume_scaffold"]["core"] = json!(core);
    // the coarse 48-gon wall is 0.19 % short of the STL volume before any
    // shell; the shell line adds 0.01 %
    v["volume_scaffold"]["max_relative_geometry_error"] = json!(0.005);
    let (z0, x0) = (24.3, 0.7);
    v["shells"] = json!([
        {"name": "baffle", "polygon": [[-9.0, -9.0, z0], [x0, -9.0, z0], [x0, 9.0, z0], [-9.0, 9.0, z0]]},
        {"name": "plate", "polygon": [[-2.0, -1.5, 11.0], [2.0, -1.5, 11.8], [2.0, 1.5, 11.8], [-2.0, 1.5, 11.0]]}
    ]);
    let plan: duct::Partition = serde_json::from_value(v).unwrap();
    let dir = std::env::temp_dir().join(format!("kofl-shell-{core}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let d = MeshDirectory::open(&dir).unwrap();
    let r = duct::generate(
        &stl,
        &ducts,
        &plan,
        Mode::Polyhedral,
        None,
        &d,
        std::ffi::OsStr::new("m"),
        mesh_engine::handoff::HandoffFormat::Vtu,
    )
    .unwrap();
    let _ = std::fs::remove_dir_all(&dir);
    check_shells(&r, x0);
    r
}

fn check_shells(r: &serde_json::Value, x0: f64) {
    let shells = r["shells"].as_array().unwrap();
    assert_eq!(shells.len(), 2, "{r:#}");
    let f = |s: &serde_json::Value, k: &str| s[k].as_f64().unwrap();
    for s in shells {
        assert!(
            (f(s, "front_area") / f(s, "back_area") - 1.0).abs() < 1e-9,
            "{s}"
        );
        let cut = s["cut_cells"].as_u64().or(s["front_faces"].as_u64());
        assert!(cut.unwrap() > 10, "{s}");
    }
    // the plate lies inside the fluid: its wall is its polygon
    let plate = &shells[1];
    assert!(
        (f(plate, "front_area") / f(plate, "polygon_area") - 1.0).abs() < 1e-9,
        "{plate}"
    );
    // the baffle covers the 48-gon section on the side x < x0
    let (r0, n) = (5.0f64, 48.0);
    let ring = 0.5 * n * r0 * r0 * (std::f64::consts::TAU / n).sin();
    let segment = r0 * r0 * (x0 / r0).acos() - x0 * (r0 * r0 - x0 * x0).sqrt();
    let expected = ring - segment * ring / (std::f64::consts::PI * r0 * r0);
    let baffle = f(&shells[0], "front_area");
    assert!(
        (baffle / expected - 1.0).abs() < 0.01,
        "{baffle} vs {expected}"
    );
    assert!(
        (f(r, "volume") / f(r, "source_volume") - 1.0).abs() < 1e-2,
        "{r:#}"
    );
}

// Keep (SW-18): a half baffle (from the wall to a free edge across the
// pipe) and a tilted plate with free edges all round, built into the
// scaffold of the default median-dual core (SHL-01 to SHL-03): both sides
// of each shell are walls of the same area, the plate covers its polygon,
// the baffle the pipe section on its side, the mesh keeps its volume and
// passes the solver rules (the build refuses otherwise). Breaks when the
// wall imprint, the shell recovery, the crack or the layers at the
// junction go wrong; the unit tests see each step alone.
#[test]
fn a_half_baffle_and_a_plate_in_the_scaffold() {
    let r = pipe_with_shells("median_dual");
    for s in r["shells"].as_array().unwrap() {
        assert_eq!(s["method"], "scaffold_baffle_v1", "{s}");
    }
}

// Keep (SW-18): the same shells cut into the finished cells of the convex
// Voronoi core (the exact-cut fallback): same areas and volume. Breaks
// when the cut meets layer columns or the report and hand-off lose the
// patches; the unit tests run on a hexahedral grid only.
#[test]
fn a_half_baffle_and_a_plate_cut_into_the_voronoi_core() {
    pipe_with_shells("voronoi");
}
