//! Patch proposals on bare closed STL shells (BC-03).

use mesh_core::geometry::Vec3;
use mesh_engine::patches::{PatchRole, PatchSurface};
use mesh_engine::proposals::{Kind, ProposalOptions, Suggestion};

#[path = "../../../test-fixtures/shapes.rs"]
mod shapes;

const JACKET: &[u8] = include_bytes!("../../../test-fixtures/cooling-jacket/cooling-jacket.stl");
const MANIFOLD: &[u8] = include_bytes!("../../../test-fixtures/manifold/manifold.stl");
const MANIFOLD_PORTS: &str = include_str!("../../../test-fixtures/manifold/patches.json");
const UNIT_BOX: &[u8] = include_bytes!("../../../test-fixtures/unit-box.stl");

fn names(p: &mesh_engine::proposals::Proposal) -> Vec<(&str, PatchRole)> {
    p.map
        .patches
        .iter()
        .map(|d| (d.name.as_str(), d.role))
        .collect()
}

/// Keep rule (SW-18). Invariant: a closed round pipe (CAD strip
/// tessellation, 64 sides) is two openings and one cylindrical wall whose
/// fitted axis, radius and length are the analytic ones; a wrong Gauss-map
/// axis or normal-line point moves them.
#[test]
fn a_round_pipe_is_two_openings_and_a_cylinder() {
    let stl = shapes::ascii_stl(&shapes::strip_cylinder(1.0, -2.0, 0.5, 4.0, 64));
    let s = PatchSurface::import(&stl, None).unwrap();
    let p = s.propose_patches(ProposalOptions::default()).unwrap();
    assert_eq!(
        names(&p),
        [
            ("opening_1", PatchRole::Unspecified),
            ("opening_2", PatchRole::Unspecified),
            ("walls", PatchRole::Wall)
        ]
    );
    assert_eq!(p.needs_role, ["opening_1", "opening_2"]);
    let wall = p.regions.iter().find(|r| r.patch == "walls").unwrap();
    assert_eq!(wall.kind, Kind::Cylinder);
    assert!((wall.axis.unwrap()[2].abs() - 1.0).abs() < 1e-9);
    let c = wall.axis_point.unwrap();
    assert!((c[0] - 1.0).abs() < 1e-9 && (c[1] + 2.0).abs() < 1e-9 && (c[2] - 2.0).abs() < 1e-9);
    // corners lie on the radius; facet centroids would sit inside it
    assert!((wall.radius.unwrap() - 0.5).abs() < 1e-9, "{wall:?}");
    assert!((wall.length.unwrap() - 4.0).abs() < 1e-9);
    for r in p
        .regions
        .iter()
        .filter(|r| r.suggestion == Suggestion::Opening)
    {
        assert_eq!(r.kind, Kind::Plane);
        assert!(r.opening_score.unwrap() > 0.99);
    }
}

/// Keep rule (SW-18). Invariant: a sphere is one spherical wall with the
/// analytic centre and radius (icosphere corners lie on the sphere); the
/// box is six box sides named by side, without a role.
#[test]
fn a_sphere_and_a_box() {
    let centre = Vec3::new(0.3, -0.2, 1.5);
    let stl = shapes::ascii_stl(&shapes::icosphere(centre, 2.0, 3));
    let s = PatchSurface::import(&stl, None).unwrap();
    let p = s.propose_patches(ProposalOptions::default()).unwrap();
    assert_eq!(names(&p), [("walls", PatchRole::Wall)]);
    let r = &p.regions[0];
    assert_eq!(r.kind, Kind::Sphere);
    assert!((r.radius.unwrap() - 2.0).abs() < 1e-9);
    let c = r.axis_point.unwrap();
    assert!((Vec3::new(c[0], c[1], c[2]) - centre).length() < 1e-9);

    let s = PatchSurface::import(UNIT_BOX, None).unwrap();
    let p = s.propose_patches(ProposalOptions::default()).unwrap();
    let got: Vec<&str> = p.map.patches.iter().map(|d| d.name.as_str()).collect();
    assert_eq!(got, ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]);
    assert!(p.regions.iter().all(|r| r.kind == Kind::Plane));
    assert_eq!(p.needs_role.len(), 6);
}

/// Keep rule (SW-18). Invariant: on the Manifold surface the two openings
/// are exactly the port caps recorded by hand (independent reference), the
/// four outer box faces are box sides and the rest is one wall patch.
#[test]
fn manifold_openings_are_the_hand_recorded_ports() {
    let s = PatchSurface::import(MANIFOLD, None).unwrap();
    let p = s.propose_patches(ProposalOptions::default()).unwrap();
    let ports: serde_json::Value = serde_json::from_str(MANIFOLD_PORTS).unwrap();
    let mut want: Vec<Vec<usize>> = (0..2)
        .map(|i| {
            let r = &ports["ports"][i]["facet_range"];
            (r[0].as_u64().unwrap() as usize..r[1].as_u64().unwrap() as usize).collect()
        })
        .collect();
    let mut got: Vec<Vec<usize>> = p
        .regions
        .iter()
        .filter(|r| r.suggestion == Suggestion::Opening)
        .map(|r| r.members.clone())
        .collect();
    want.sort();
    got.sort();
    assert_eq!(got, want);
    let n = names(&p);
    assert_eq!(n.len(), 7, "{n:?}");
    assert_eq!(n.iter().filter(|(_, r)| *r == PatchRole::Wall).count(), 1);
}

/// Keep rule (SW-18). Invariant: on the cooling jacket the two openings are
/// the two 96-gon pipe ends; the block faces (squares) are not openings,
/// and the inner and outer jacket walls fit cylinders of radius 60 and 66.
#[test]
fn jacket_openings_and_cylinders() {
    let s = PatchSurface::import(JACKET, None).unwrap();
    let p = s.propose_patches(ProposalOptions::default()).unwrap();
    let openings: Vec<_> = p
        .regions
        .iter()
        .filter(|r| r.suggestion == Suggestion::Opening)
        .collect();
    assert_eq!(openings.len(), 2);
    assert!(openings.iter().all(|r| r.facets == 96));
    let mut radii: Vec<f64> = p
        .regions
        .iter()
        .filter(|r| r.kind == Kind::Cylinder)
        .map(|r| r.radius.unwrap())
        .collect();
    radii.sort_by(f64::total_cmp);
    radii.dedup_by(|a, b| (*a - *b).abs() < 0.01);
    assert_eq!(radii.len(), 3, "{radii:?}");
    for (r, want) in radii.iter().zip([60.0, 64.0, 66.0]) {
        assert!((r - want).abs() < 0.01, "{radii:?}");
    }
    // separate patches when asked
    let p = s
        .propose_patches(ProposalOptions {
            merge_walls: false,
            ..Default::default()
        })
        .unwrap();
    assert!(p.map.patches.iter().any(|d| d.name == "cylinder_1"));
    assert_eq!(p.map.patches.len(), p.regions.len());
}

/// Out-of-range options are refused.
#[test]
fn bad_options_are_refused() {
    let s = PatchSurface::import(UNIT_BOX, None).unwrap();
    for o in [
        ProposalOptions {
            feature_angle_deg: 0.0,
            ..Default::default()
        },
        ProposalOptions {
            plane_angle_deg: f64::NAN,
            ..Default::default()
        },
        ProposalOptions {
            fit_tolerance: 1.0,
            ..Default::default()
        },
    ] {
        assert!(s.propose_patches(o).is_err());
    }
}

/// Keep rule (SW-18). Invariant (BC-05): the skeleton is a TOML case that
/// names every patch with the FoehnCraft type of its role and holds no
/// physical value: inlet and outlet values are TODO comments, a patch
/// without a role is an empty table that FoehnCraft refuses, and the mesh
/// file is the one given.
#[test]
fn skeleton_names_types_and_leaves_values_to_foehncraft() {
    let stl = shapes::ascii_stl(&shapes::strip_cylinder(0.0, 0.0, 0.5, 4.0, 32));
    let s = PatchSurface::import(&stl, None).unwrap();
    let p = s.propose_patches(ProposalOptions::default()).unwrap();
    let (map, _) = mesh_engine::patches::set_role(&p.map, "opening_1", PatchRole::Inlet).unwrap();
    let k = s.bc_skeleton(&map, "foehncraft.vtu").unwrap();
    assert_eq!(k.needs_role, ["opening_2"]);
    let case: toml::Table = toml::from_str(&k.text).unwrap();
    assert_eq!(case["case_format"].as_str(), Some("0.7"));
    assert_eq!(case["mesh"]["file"].as_str(), Some("foehncraft.vtu"));
    assert!(case.get("fluid").is_none(), "no physical value");
    let b = case["boundary"].as_table().unwrap();
    assert_eq!(b["opening_1"]["type"].as_str(), Some("velocity-inlet"));
    assert_eq!(b["opening_1"].as_table().unwrap().len(), 1);
    assert!(b["opening_2"].as_table().unwrap().is_empty());
    assert_eq!(b["walls"]["type"].as_str(), Some("wall"));
    assert!(k.text.contains("# TODO normal_velocity"));
    assert!(s.bc_skeleton(&map, "").is_err());
}
