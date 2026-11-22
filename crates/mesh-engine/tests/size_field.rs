//! Size control of the STL mesher (ALG-05, GRD-01 to GRD-03, GRD-05, GRD-08).

use mesh_engine::{
    ErrorCode, MeshArtifact, StlMeshConfig, build_stl_mesh, check_mesh, stl_site_count,
};
use serde_json::{Value, json};

const UNIT_BOX: &[u8] = include_bytes!("../../../test-fixtures/unit-box.stl");
const SKEWED_BOX: &[u8] = include_bytes!("../../../test-fixtures/skewed-box.stl");
const CONFIG: &str = include_str!("../../../test-fixtures/stl-size-field.json");

fn config(value: Value) -> StlMeshConfig {
    serde_json::from_value(value).unwrap()
}

fn fixture() -> Value {
    serde_json::from_str(CONFIG).unwrap()
}

/// Keep rule (SW-18). Invariant: a size field places the seeds so that the
/// validated mesh follows it: cbrt(V)/h within 15 % for 80 % of the cells,
/// each source's cells at its size within 10 % (median), the refined zones
/// finer than the background, the same inputs giving the same mesh, and the
/// estimate equal to the build. A change to the density exponent, the
/// count rule, the frame transform or the evidence breaks it; the unit tests
/// check the pieces only, not the mesh they make.
#[test]
fn size_field_grades_the_stl_mesh_and_reports_each_source() {
    let artifact = build_stl_mesh(UNIT_BOX, config(fixture()), None, None).unwrap();
    let check = check_mesh(&artifact).unwrap();
    let sizing = artifact.sizing().unwrap();
    assert_eq!(check.cells, sizing.seeds);
    assert_eq!(
        stl_site_count(UNIT_BOX, &config(fixture()), None).unwrap(),
        sizing.seeds
    );
    assert_eq!(sizing.relaxation_iterations, 20);
    assert!(sizing.size_ratio.p10 > 0.85 && sizing.size_ratio.p90 < 1.15);
    assert!((sizing.size_ratio.median - 1.0).abs() < 0.05);
    assert!(sizing.neighbour_ratio.median < 1.3);
    assert_eq!(sizing.interior_faces, check.internal_faces);
    for (source, target) in sizing.sources.iter().zip([0.06, 0.1]) {
        assert_eq!(source.target_size, target);
        assert!(
            source.cells > 5 && source.governs >= source.cells,
            "{source:?}"
        );
        let achieved = source.achieved.as_ref().unwrap();
        assert!((achieved.median / target - 1.0).abs() < 0.1, "{source:?}");
        assert!(achieved.max < 0.25 * 0.8);
    }
    // The band lies along the frame's y = 0: the plane x = y of the box.
    let band = &sizing.sources[1];
    assert_eq!((band.name.as_str(), band.kind.as_str()), ("band", "region"));

    let again = build_stl_mesh(UNIT_BOX, config(fixture()), None, None).unwrap();
    assert_eq!(again.build().mesh_digest, artifact.build().mesh_digest);
    let mut other = fixture();
    other["seed"] = json!(2);
    let other = build_stl_mesh(UNIT_BOX, config(other), None, None).unwrap();
    assert_ne!(other.build().mesh_digest, artifact.build().mesh_digest);

    // The evidence is part of the config digest: editing it is detected.
    let mut text: Value = serde_json::to_value(&artifact).unwrap();
    text["sizing"]["seeds"] = json!(1);
    let edited: MeshArtifact = serde_json::from_value(text).unwrap();
    assert_eq!(
        check_mesh(&edited).unwrap_err().code,
        ErrorCode::DigestMismatch
    );
}

#[test]
fn every_shape_grades_an_angled_surface() {
    let shapes = [
        json!({"type":"segment","start":[0.2,0.2,0.2],"end":[0.8,0.2,0.2]}),
        json!({"type":"plane","point":[0.5,0.5,0.5],"normal":[0,0,1]}),
        json!({"type":"box","min":[0.1,0.1,0.1],"max":[0.3,0.3,0.3]}),
        json!({"type":"sphere","centre":[0.5,0.5,0.5],"radius":0.1}),
        json!({"type":"cylinder","start":[0.5,0.0,0.5],"end":[0.5,1.0,0.5],"radius":0.05}),
        json!({"type":"cone","start":[0.5,0.5,0.0],"end":[0.5,0.5,1.0],"start_radius":0.2,"end_radius":0.0}),
    ];
    for shape in shapes {
        let value = json!({
            "schema_version": 1,
            "size_field": {
                "min_size": 0.08, "max_size": 0.3, "gradation": 1.4, "relaxation_iterations": 8,
                "sources": [{"kind":"datum","name":"d","shape":shape,"size":0.08,"distance":0.02}]
            },
            "seed": 3
        });
        let artifact = build_stl_mesh(SKEWED_BOX, config(value), None, None).unwrap();
        check_mesh(&artifact).unwrap();
        let source = &artifact.sizing().unwrap().sources[0];
        assert!(source.governs > 0, "{shape}");
    }
}

#[test]
fn invalid_size_fields_are_diagnosed() {
    let code = |value: Value| {
        build_stl_mesh(UNIT_BOX, config(value), None, None)
            .unwrap_err()
            .code
    };
    let mut both = fixture();
    both["sites"] = json!([[0.5, 0.5, 0.5]]);
    assert_eq!(code(both), ErrorCode::InvalidSites);
    let mut small = fixture();
    small["size_field"]["sources"][0]["size"] = json!(0.01);
    assert_eq!(code(small), ErrorCode::InvalidSizeField);
    let mut sphere_region = fixture();
    sphere_region["size_field"]["sources"][1]["shape"] =
        json!({"type":"sphere","centre":[0.5,0.5,0.5],"radius":0.1});
    assert_eq!(code(sphere_region), ErrorCode::InvalidSizeField);
    let mut frame = fixture();
    frame["size_field"]["sources"][1]["frame"] = json!("nowhere");
    assert_eq!(code(frame), ErrorCode::InvalidSizeField);
    let mut parallel = fixture();
    parallel["size_field"]["frames"][0]["y_axis"] = json!([2.0, 2.0, 0.0]);
    assert_eq!(code(parallel), ErrorCode::InvalidSizeField);
    let fine = json!({
        "schema_version": 1,
        "size_field": {"min_size": 0.01, "max_size": 0.01, "gradation": 1.2}
    });
    assert_eq!(code(fine), ErrorCode::InvalidSizeField);
}
