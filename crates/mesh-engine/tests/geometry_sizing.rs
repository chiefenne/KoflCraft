//! Geometry-derived size sources of the STL mesher (GRD-04): curvature,
//! gaps and feature edges.

#[path = "../../../test-fixtures/shapes.rs"]
mod shapes;

use mesh_core::geometry::Vec3;
use mesh_engine::{ErrorCode, StlMeshConfig, build_stl_mesh, check_mesh, stl_site_count};
use serde_json::{Value, json};

fn config(sources: Value, min: f64, max: f64) -> StlMeshConfig {
    serde_json::from_value(json!({
        "schema_version": 1,
        "size_field": {"min_size": min, "max_size": max, "gradation": 1.3, "sources": sources},
        "seed": 1
    }))
    .unwrap()
}

/// Keep rule (SW-18). Invariant: a gap source puts the asked number of
/// cells across a thin slab, counted along the wall normals of the built,
/// validated mesh, where the background size alone gives one cell; the
/// estimate equals the build. Dropping the reach (the size would grow
/// towards the middle of the gap), a wrong thickness or a broken count
/// breaks this; the unit tests check thickness and sizes, not cells.
#[test]
fn gap_source_puts_cells_across_a_thin_slab() {
    let stl = shapes::ascii_stl(&shapes::cuboid(Vec3::new(1.0, 1.0, 0.1)));
    let sources = json!([{"kind": "gap", "name": "slab", "cells_across": 3}]);
    let artifact = build_stl_mesh(&stl, config(sources.clone(), 0.01, 0.1), None, None).unwrap();
    let check = check_mesh(&artifact).unwrap();
    let sizing = artifact.sizing().unwrap();
    assert_eq!(check.cells, sizing.seeds);
    assert_eq!(
        stl_site_count(&stl, &config(sources, 0.01, 0.1), None).unwrap(),
        sizing.seeds
    );
    let source = &sizing.sources[0];
    assert_eq!(source.kind, "gap");
    assert!((source.target_size - 0.1 / 3.0).abs() < 1e-12);
    let derived = source.derived.as_ref().unwrap();
    assert_eq!(derived.method, "medial_ball_thickness_v1");
    // Top and bottom, two triangles each; the sides are 1 apart (size 1/3
    // is above max_size, so they ask for nothing).
    assert_eq!((derived.primitives, derived.at_min_size), (4, 0));
    let across = derived.cells_across.as_ref().unwrap();
    assert_eq!(derived.samples, Some(4));
    assert!(across.min >= 3.0, "{across:?}");
    assert_eq!(derived.samples_below_target, Some(0));
    let achieved = source.achieved.as_ref().unwrap();
    assert!(
        (achieved.median / source.target_size - 1.0).abs() < 0.1,
        "{achieved:?}"
    );
}

/// Keep rule (SW-18). Invariant: a curvature source refines the curved
/// wall of a CAD-style strip cylinder to angle·R and leaves the flat caps
/// to the background, and the result reports the method. A wrong curvature
/// on strip triangles or curvature leaking from the rim (a feature edge)
/// breaks this; the unit test checks the estimator only.
#[test]
fn curvature_source_refines_the_curved_wall() {
    let stl = shapes::ascii_stl(&shapes::strip_cylinder(0.5, 0.5, 0.25, 0.6, 48));
    let sources = json!([{"kind": "curvature", "name": "wall", "angle": 15, "distance": 0.03}]);
    let artifact = build_stl_mesh(&stl, config(sources, 0.01, 0.15), None, None).unwrap();
    check_mesh(&artifact).unwrap();
    let source = &artifact.sizing().unwrap().sources[0];
    let target = 15f64.to_radians() * 0.25;
    // Chord facets turn by 2π/n over a width 2R·sin(π/n).
    assert!(
        (source.target_size / target - 1.0).abs() < 0.002,
        "{source:?}"
    );
    let derived = source.derived.as_ref().unwrap();
    assert_eq!(derived.method, "normal_cycle_tensor_v1");
    assert_eq!(derived.primitives, 96, "the side strips, not the caps");
    assert!(source.cells > 20, "{source:?}");
    let achieved = source.achieved.as_ref().unwrap();
    assert!((achieved.median / target - 1.0).abs() < 0.1, "{achieved:?}");
}

/// Keep rule (SW-18). Invariant: a feature-edge source finds the twelve
/// edges of a box (not the six face diagonals) and the cells along them
/// reach its size. A dihedral test with the wrong sign or threshold breaks
/// this; no other test builds a mesh from edges.
#[test]
fn feature_edges_refine_the_box_edges() {
    let stl = shapes::ascii_stl(&shapes::cuboid(Vec3::new(1.0, 1.0, 1.0)));
    let sources =
        json!([{"kind": "feature_edges", "name": "edges", "size": 0.05, "distance": 0.05}]);
    let artifact = build_stl_mesh(&stl, config(sources, 0.02, 0.2), None, None).unwrap();
    check_mesh(&artifact).unwrap();
    let source = &artifact.sizing().unwrap().sources[0];
    assert_eq!(source.target_size, 0.05);
    let derived = source.derived.as_ref().unwrap();
    assert_eq!(
        (derived.primitives, derived.method.as_str()),
        (12, "dihedral_angle_v1")
    );
    assert!(source.cells > 20, "{source:?}");
    let achieved = source.achieved.as_ref().unwrap();
    assert!((achieved.median / 0.05 - 1.0).abs() < 0.1, "{achieved:?}");
}

#[test]
fn invalid_geometry_sources_are_diagnosed() {
    let stl = shapes::ascii_stl(&shapes::cuboid(Vec3::new(1.0, 1.0, 1.0)));
    let code = |sources: Value| {
        build_stl_mesh(&stl, config(sources, 0.05, 0.2), None, None)
            .unwrap_err()
            .code
    };
    for sources in [
        json!([{"kind": "curvature", "name": "c", "angle": 0}]),
        json!([{"kind": "curvature", "name": "c", "angle": 10, "feature_angle": 200}]),
        json!([{"kind": "gap", "name": "g", "cells_across": 0.5}]),
        json!([{"kind": "gap", "name": "g", "cells_across": 3, "wedge_angle": 120}]),
        json!([{"kind": "feature_edges", "name": "e", "size": 0.01}]),
        json!([
            {"kind": "gap", "name": "x", "cells_across": 3},
            {"kind": "curvature", "name": "x", "angle": 10}
        ]),
    ] {
        assert_eq!(
            code(sources.clone()),
            ErrorCode::InvalidSizeField,
            "{sources}"
        );
    }
    // A curvature source on a box finds nothing and asks for nothing.
    let sources = json!([{"kind": "curvature", "name": "c", "angle": 10}]);
    let artifact = build_stl_mesh(&stl, config(sources, 0.05, 0.2), None, None).unwrap();
    let source = &artifact.sizing().unwrap().sources[0];
    assert_eq!((source.target_size, source.governs), (0.2, 0));
    assert_eq!(source.derived.as_ref().unwrap().primitives, 0);
}
