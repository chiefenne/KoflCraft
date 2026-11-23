//! Boundary patches defined after import (BC-01, BC-02, BC-04, BC-06).

use mesh_engine::patches::{
    Contain, PatchMap, PatchRole, PatchSurface, Selection, Target, merge, rename, set_role,
};
use mesh_engine::{
    ErrorCode, ExportData, ExportFormat, StlMeshConfig, build_stl_mesh, export_mesh,
};

const UNIT_BOX: &[u8] = include_bytes!("../../../test-fixtures/unit-box.stl");
const SKEWED_BOX: &[u8] = include_bytes!("../../../test-fixtures/skewed-box.stl");
const SITES: &str = include_str!("../../../test-fixtures/stl-sites.json");
const MANIFOLD: &[u8] = include_bytes!("../../../test-fixtures/manifold/manifold.stl");
const MANIFOLD_PORTS: &str = include_str!("../../../test-fixtures/manifold/patches.json");

fn target(name: &str) -> Target {
    Target {
        name: name.into(),
        role: None,
        within: None,
    }
}

fn names(map: &PatchMap) -> Vec<&str> {
    map.patches.iter().map(|p| p.name.as_str()).collect()
}

fn facets_of(map: &PatchMap, name: &str) -> Vec<usize> {
    let id = map.patches.iter().position(|p| p.name == name).unwrap();
    (0..map.facet_patch.len())
        .filter(|&f| map.facet_patch[f] == id)
        .collect()
}

/// Keep rule (SW-18). Invariant: the map starts from the STL solids with
/// hand-off names, the source names kept and roles from the names, and the
/// summary reports exact areas, normals, regions and neighbours of the
/// unit box. A change to the name conversion, the role rule, the normal or
/// the adjacency breaks it; the export tests see only the final names.
#[test]
fn solid_names_start_the_map_with_exact_evidence() {
    let s = PatchSurface::import(UNIT_BOX, None).unwrap();
    let map = s.initial_map();
    assert_eq!(
        names(&map),
        [
            "west_inlet",
            "east_outlet",
            "south_wall",
            "north_wall",
            "floor",
            "ceiling"
        ]
    );
    assert_eq!(map.patches[0].source_name.as_deref(), Some("west inlet"));
    assert_eq!(map.patches[4].source_name, None);
    let roles: Vec<_> = map.patches.iter().map(|p| p.role).collect();
    use PatchRole::*;
    assert_eq!(roles, [Inlet, Outlet, Wall, Wall, Unspecified, Unspecified]);
    let summary = s.summary(&map).unwrap();
    assert!(summary.every_facet_in_one_patch);
    assert_eq!(summary.area, 6.0);
    let west = &summary.patches[0];
    assert_eq!((west.facets, west.area, west.regions), (2, 1.0, 1));
    assert_eq!(west.normal, Some([-1.0, 0.0, 0.0]));
    assert_eq!(west.max_normal_deviation_deg, 0.0);
    assert_eq!(west.bounds, Some([[0.0, 0.0, 0.0], [0.0, 1.0, 1.0]]));
    assert_eq!(west.neighbours.len(), 4);
    assert!(!west.neighbours.contains(&"east_outlet".to_string()));
    assert_eq!(west.foehncraft_type, Some("velocity-inlet"));
}

/// Keep rule (SW-18). Invariant: region growing from a point on each port
/// of the Manifold surface (90,894 facets, sharp rims, rounded pipes)
/// returns exactly the cap facets recorded by hand in the fixture's patch
/// map, and the plane selection of the inlet plane returns the same set.
/// A wrong neighbour across an edge, a wrong angle test or a nearest-facet
/// error changes the sets; no other test has an independent reference.
#[test]
fn region_growing_and_planes_recover_the_manifold_ports() {
    let ports: serde_json::Value = serde_json::from_str(MANIFOLD_PORTS).unwrap();
    let s = PatchSurface::import(MANIFOLD, None).unwrap();
    let start = s.initial_map();
    let range = |i: usize| -> Vec<usize> {
        let r = &ports["ports"][i]["facet_range"];
        (r[0].as_u64().unwrap() as usize..r[1].as_u64().unwrap() as usize).collect()
    };
    let grow = |point: [f64; 3]| Selection::Grow {
        seed_facet: None,
        seed_point: Some(point),
        feature_angle_deg: 30.0,
    };
    let (map, change) = s
        .select(&start, &grow([0.0, 0.0, 0.0]), &target("inlet"))
        .unwrap();
    assert_eq!(facets_of(&map, "inlet"), range(0));
    assert_eq!(change.facets_changed, 510);
    let (map, _) = s
        .select(&map, &grow([-120.0, 360.0, 40.0]), &target("outlet"))
        .unwrap();
    assert_eq!(facets_of(&map, "outlet"), range(1));
    let (map, _) = rename(&map, "surface", "walls").unwrap();
    assert_eq!(names(&map), ["walls", "inlet", "outlet"]);
    // a rename keeps the role; the role is set on its own (BC-04)
    assert_eq!(map.patches[0].role, PatchRole::Unspecified);
    let (map, _) = set_role(&map, "walls", PatchRole::Wall).unwrap();
    let roles: Vec<_> = map.patches.iter().map(|p| p.role).collect();
    assert_eq!(
        roles,
        [PatchRole::Wall, PatchRole::Inlet, PatchRole::Outlet]
    );
    let summary = s.summary(&map).unwrap();
    let inlet = &summary.patches[1];
    let area = ports["ports"][0]["area"].as_f64().unwrap();
    assert!((inlet.area / area - 1.0).abs() < 1e-12, "{}", inlet.area);
    assert_eq!(inlet.normal, Some([0.0, -1.0, 0.0]));
    assert_eq!(inlet.neighbours, ["walls"]);

    let plane = Selection::Plane {
        point: [0.0, 0.0, 0.0],
        normal: [0.0, 1.0, 0.0],
        tolerance: None,
    };
    let (by_plane, _) = s.select(&start, &plane, &target("inlet")).unwrap();
    assert_eq!(facets_of(&by_plane, "inlet"), range(0));
}

/// Keep rule (SW-18). Invariant: box, normal, `within`, merge, split and
/// role edits move exactly the facets they name, drop emptied patches,
/// keep the order and give every facet one patch. A change to the box
/// containment, the `within` mask, the renumbering or the split naming
/// breaks it.
#[test]
fn selections_merge_and_split_on_the_skewed_box() {
    let s = PatchSurface::import(SKEWED_BOX, None).unwrap();
    let start = s.initial_map();
    let walls: Vec<String> = ["south_wall", "north_wall", "floor", "ceiling"]
        .map(String::from)
        .into();
    let (map, change) = merge(&start, &walls, "walls", Some(PatchRole::Wall)).unwrap();
    assert_eq!(names(&map), ["west_inlet", "east_outlet", "walls"]);
    assert_eq!(change.removed, walls);
    assert_eq!(facets_of(&map, "walls").len(), 8);

    // everything into one patch, then split back at the box edges
    let all: Vec<String> = names(&map).into_iter().map(String::from).collect();
    let (one, _) = merge(&map, &all, "all", None).unwrap();
    assert_eq!(s.summary(&one).unwrap().patches[0].regions, 1);
    let e = s.split(&one, "all", None, None).unwrap_err();
    assert_eq!(e.code, ErrorCode::EmptySelection);
    let (six, change) = s.split(&one, "all", Some(10.0), None).unwrap();
    assert_eq!(
        names(&six),
        ["all_1", "all_2", "all_3", "all_4", "all_5", "all_6"]
    );
    assert_eq!(change.facets_changed, 10);
    for p in s.summary(&six).unwrap().patches {
        assert_eq!(
            (p.facets, p.regions, p.max_normal_deviation_deg),
            (2, 1, 0.0)
        );
    }

    // the facets facing roughly +y, but only out of the walls
    let normal = Selection::Normal {
        direction: [0.0, 1.0, 0.0],
        max_angle_deg: 60.0,
    };
    let within = Target {
        name: "top".into(),
        role: Some(PatchRole::Symmetry),
        within: Some(vec!["walls".into()]),
    };
    let (top, _) = s.select(&map, &normal, &within).unwrap();
    let summary = s.summary(&top).unwrap();
    let t = summary.patches.iter().find(|p| p.name == "top").unwrap();
    assert_eq!((t.facets, t.role), (2, PatchRole::Symmetry));

    // a box around the whole surface by centroid takes everything
    let b = s.surface_bounds();
    let all_box = Selection::Box {
        min: b[0],
        max: b[1],
        contain: Contain::AllCorners,
    };
    let (boxed, change) = s.select(&map, &all_box, &target("walls")).unwrap();
    assert_eq!(names(&boxed), ["walls"]);
    assert_eq!(change.removed, ["west_inlet", "east_outlet"]);
    assert_eq!(change.facets_changed, 4);

    let (role, _) = set_role(&map, "east_outlet", PatchRole::Wall).unwrap();
    assert_eq!(role.patches[1].role, PatchRole::Wall);
}

/// Keep rule (SW-18). Invariant: a map is refused unless it belongs to the
/// STL and gives every facet exactly one patch with a valid, unique name
/// (BC-06), and an operation that selects or changes nothing writes no
/// map. Loosening any check in `PatchSurface::check` breaks it.
#[test]
fn invalid_maps_and_empty_operations_are_diagnosed() {
    let s = PatchSurface::import(UNIT_BOX, None).unwrap();
    let good = s.initial_map();
    let code = |edit: &dyn Fn(&mut PatchMap)| {
        let mut m = good.clone();
        edit(&mut m);
        s.check(&m).unwrap_err().code
    };
    assert_eq!(
        code(&|m| m.source_sha256 = "0".repeat(64)),
        ErrorCode::InvalidPatchMap
    );
    assert_eq!(code(&|m| m.facet_count = 11), ErrorCode::InvalidPatchMap);
    assert_eq!(
        code(&|m| m.facet_patch.truncate(11)),
        ErrorCode::InvalidPatchMap
    );
    assert_eq!(code(&|m| m.facet_patch[0] = 6), ErrorCode::InvalidPatchMap);
    assert_eq!(
        code(&|m| m.facet_patch[..2].fill(1)),
        ErrorCode::InvalidPatchMap
    );
    assert_eq!(
        code(&|m| m.patches[1].name = "west_inlet".into()),
        ErrorCode::InvalidPatchMap
    );
    assert_eq!(
        code(&|m| m.patches[1].name = "east outlet".into()),
        ErrorCode::InvalidPatchName
    );
    assert_eq!(
        code(&|m| m.format_version = "2.0".into()),
        ErrorCode::InvalidPatchMap
    );
    let other = PatchSurface::import(SKEWED_BOX, None).unwrap();
    assert_eq!(
        other.check(&good).unwrap_err().code,
        ErrorCode::InvalidPatchMap
    );

    let nowhere = Selection::Plane {
        point: [0.5, 0.5, 0.5],
        normal: [1.0, 0.0, 0.0],
        tolerance: None,
    };
    let e = s.select(&good, &nowhere, &target("mid")).unwrap_err();
    assert_eq!(e.code, ErrorCode::EmptySelection);
    let west = Selection::Normal {
        direction: [-1.0, 0.0, 0.0],
        max_angle_deg: 1.0,
    };
    let e = s.select(&good, &west, &target("west_inlet")).unwrap_err();
    assert_eq!(e.code, ErrorCode::EmptySelection);
    let e = s.select(&good, &west, &target("bad name")).unwrap_err();
    assert_eq!(e.code, ErrorCode::InvalidPatchName);
    let mut within = target("x");
    within.within = Some(vec!["nope".into()]);
    assert_eq!(
        s.select(&good, &west, &within).unwrap_err().code,
        ErrorCode::UnknownPatch
    );
    assert_eq!(
        rename(&good, "floor", "ceiling").unwrap_err().code,
        ErrorCode::InvalidPatchName
    );
    assert_eq!(
        rename(&good, "roof", "top").unwrap_err().code,
        ErrorCode::UnknownPatch
    );
    let zero = Selection::Normal {
        direction: [0.0, 0.0, 0.0],
        max_angle_deg: 1.0,
    };
    assert_eq!(
        s.select(&good, &zero, &target("x")).unwrap_err().code,
        ErrorCode::InvalidArguments
    );
    assert_eq!(
        set_role(&good, "floor", PatchRole::Unspecified)
            .unwrap_err()
            .code,
        ErrorCode::EmptySelection
    );
}

/// Keep rule (SW-18). Invariant: a mesh built with a patch map carries the
/// map's names in MSH 4.1 physical groups, the FoehnCraft VTU and OpenFOAM
/// `boundary`, and the map's roles in the VTU, the OpenFOAM patch types and
/// the evidence, also where a role differs from what the name suggests. A
/// change that drops the saved roles or renames patches on the way breaks
/// it; the FoehnCraft round trip only compares KoflCraft with itself.
#[test]
fn patch_map_names_and_roles_reach_every_export() {
    let s = PatchSurface::import(UNIT_BOX, None).unwrap();
    let start = s.initial_map();
    let (map, _) = s
        .select(
            &start,
            &Selection::Plane {
                point: [0.0, 0.0, 0.0],
                normal: [1.0, 0.0, 0.0],
                tolerance: None,
            },
            &target("inflow"),
        )
        .unwrap();
    let (map, _) = set_role(&map, "inflow", PatchRole::Inlet).unwrap();
    let rest: Vec<String> = ["south_wall", "north_wall", "floor"]
        .map(String::from)
        .into();
    let (map, _) = merge(&map, &rest, "walls", Some(PatchRole::Wall)).unwrap();
    let (map, _) = set_role(&map, "ceiling", PatchRole::Symmetry).unwrap();
    assert_eq!(names(&map), ["east_outlet", "walls", "ceiling", "inflow"]);

    let config: StlMeshConfig = serde_json::from_str(SITES).unwrap();
    let artifact = build_stl_mesh(UNIT_BOX, config.clone(), None, Some(&map)).unwrap();
    let text = |format| match export_mesh(&artifact, format).unwrap() {
        ExportData::Gmsh(b) | ExportData::FoehnCraft(b) | ExportData::Vtu(b) => {
            vec![("file".to_string(), String::from_utf8(b).unwrap())]
        }
        ExportData::Bundle(files) => files
            .into_iter()
            .map(|(n, b)| (n.to_string(), String::from_utf8(b).unwrap()))
            .collect(),
    };
    let mut expected: Vec<&str> = names(&map);
    expected.sort_unstable();

    let msh = &text(ExportFormat::Gmsh(mesh_engine::GmshVersion::V4_1))[0].1;
    let section = msh.split("$PhysicalNames\n").nth(1).unwrap();
    let section = section.split("$EndPhysicalNames").next().unwrap();
    let mut groups: Vec<&str> = section
        .lines()
        .skip(1)
        .filter(|l| l.starts_with("2 "))
        .map(|l| l.split('"').nth(1).unwrap())
        .collect();
    groups.sort_unstable();
    assert_eq!(groups, expected);

    let vtu = &text(ExportFormat::FoehnCraft)[0].1;
    let strings = |name: &str| -> Vec<String> {
        let body = vtu.split(&format!("Name=\"{name}\"")).nth(1).unwrap();
        let body = body.split('>').nth(1).unwrap().split('<').next().unwrap();
        let bytes: Vec<u8> = body
            .split_whitespace()
            .map(|b| b.parse().unwrap())
            .collect();
        bytes
            .split(|&b| b == 0)
            .filter(|s| !s.is_empty())
            .map(|s| String::from_utf8(s.to_vec()).unwrap())
            .map(|s| s.split_once('=').unwrap().1.to_string())
            .collect()
    };
    let (vtu_names, vtu_roles) = (strings("boundary_names"), strings("boundary_roles"));
    let mut sorted = vtu_names.clone();
    sorted.sort_unstable();
    assert_eq!(sorted, expected);
    for (name, role) in [
        ("east_outlet", "outlet"),
        ("walls", "wall"),
        ("ceiling", "symmetry"),
        ("inflow", "inlet"),
    ] {
        let i = vtu_names.iter().position(|n| n == name).unwrap();
        assert_eq!(vtu_roles[i], role, "{name}");
    }

    let foam = text(ExportFormat::OpenFoam);
    let boundary = &foam.iter().find(|(n, _)| n == "boundary").unwrap().1;
    for (name, kind) in [
        ("east_outlet", "patch"),
        ("walls", "wall"),
        ("ceiling", "symmetry"),
        ("inflow", "patch"),
    ] {
        let block = boundary.split(&format!("    {name}\n")).nth(1).unwrap();
        assert!(
            block
                .trim_start()
                .starts_with(&format!("{{\n        type        {kind};")),
            "{name}"
        );
    }
    let evidence: serde_json::Value = serde_json::from_str(
        &foam
            .iter()
            .find(|(n, _)| n == "mesh-evidence.json")
            .unwrap()
            .1,
    )
    .unwrap();
    let inflow = evidence["patches"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["name"] == "inflow")
        .unwrap();
    assert_eq!(inflow["role"], "inlet");
    assert!(
        !evidence["limitations"]
            .to_string()
            .contains("no boundary-condition role")
    );

    // the saved roles are covered by the artifact digests
    let mut edited = serde_json::to_value(&artifact).unwrap();
    edited["surface"]["patch_roles"][0] = "wall".into();
    let edited: mesh_engine::MeshArtifact = serde_json::from_value(edited).unwrap();
    let e = export_mesh(&edited, ExportFormat::FoehnCraft)
        .err()
        .unwrap();
    assert_eq!(e.code, ErrorCode::DigestMismatch);
    let _ = config;

    // a map and facet_patches together are refused
    let mut both: StlMeshConfig = serde_json::from_str(SITES).unwrap();
    both.facet_patches = Some(vec!["a".into(); 12]);
    let e = build_stl_mesh(UNIT_BOX, both, None, Some(&map)).unwrap_err();
    assert_eq!(e.code, ErrorCode::InvalidPatchMap);
}
