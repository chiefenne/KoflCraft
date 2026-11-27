//! JSON Schemas of the tool inputs and outputs (AG-6, AG-12), written by
//! hand from the engine's types. The inputs are the same for every front
//! end; configurations are given inline, files by path.

use mesh_engine::{SCHEMA_VERSION, SURFACE_ARTIFACT_VERSION};
use serde_json::{Value, json};

pub(crate) fn sweep_schema() -> Value {
    let positive = json!({"type":"number","exclusiveMinimum":0});
    let growth = json!({"type":"number","minimum":1});
    let block = object_schema(json!({
        "kind":{"const":"block"},"dimensions_mm":{"type":"array","minItems":3,"maxItems":3,"items":positive},
        "inlet":{"enum":["x_min","y_min","z_min"]}
    }));
    let pipe =
        object_schema(json!({"kind":{"const":"pipe"},"diameter_mm":positive,"length_mm":positive}));
    object_schema(
        json!({"schema_version":{"const":SCHEMA_VERSION},"geometry":{"oneOf":[block,pipe]},"core_spacing_mm":positive,
        "wall_layers":{"type":"integer","minimum":1},"first_layer_mm":positive,"layer_growth":growth,"max_neighbour_size_ratio":growth}),
    )
}

pub(crate) fn point() -> Value {
    json!({ "type": "array", "items": { "type": "number" }, "minItems": 3, "maxItems": 3 })
}

pub(crate) fn path() -> Value {
    json!({ "type": "string", "minLength": 1, "description": "A file path, relative to the working directory (or the first write root of `kofl mcp`)" })
}

pub(crate) fn mesh_config() -> Value {
    let point = point();
    json!({
        "type": "object",
        "additionalProperties": false,
        "required": ["schema_version", "domain", "sites"],
        "properties": {
            "schema_version": { "const": SCHEMA_VERSION },
            "domain": {
                "type": "object",
                "additionalProperties": false,
                "required": ["min", "max"],
                "properties": { "min": point, "max": point }
            },
            "sites": { "type": "array", "items": point, "minItems": 1 },
            "fixture_id": { "type": ["string", "null"] },
            "seed": { "type": ["integer", "null"], "minimum": 0 }
        }
    })
}

pub(crate) fn stl_config() -> Value {
    let config = mesh_config();
    json!({
        "type": "object",
        "additionalProperties": false,
        "required": ["schema_version"],
        "description": "Give either `sites` or a `size_field` that places them (ALG-05).",
        "properties": {
            "schema_version": config["properties"]["schema_version"],
            "sites": config["properties"]["sites"],
            "size_field": size_field(),
            "fixture_id": config["properties"]["fixture_id"],
            "seed": config["properties"]["seed"],
            "facet_patches": { "type": ["array", "null"], "items": { "type": "string" } }
        }
    })
}

/// The size field of `mesh_engine::sizing::SizeFieldConfig` (GRD-01 to GRD-05, GRD-08).
pub(crate) fn size_field() -> Value {
    let positive = json!({"type": "number", "exclusiveMinimum": 0});
    let nonnegative = json!({"type": "number", "minimum": 0});
    let p = point();
    let name = json!({"type": "string", "minLength": 1});
    let shape = |kind: &str, props: Value| {
        let mut props = props;
        props["type"] = json!({ "const": kind });
        object_schema(props)
    };
    let point_shape = shape("point", json!({"point": p}));
    let segment = shape("segment", json!({"start": p, "end": p}));
    let plane = shape("plane", json!({"point": p, "normal": p}));
    let cuboid = shape("box", json!({"min": p, "max": p}));
    let slab = shape(
        "slab",
        json!({"axis": {"enum": ["x", "y", "z"]}, "min": {"type": "number"}, "max": {"type": "number"}}),
    );
    let sphere = shape("sphere", json!({"centre": p, "radius": positive}));
    let cylinder = shape(
        "cylinder",
        json!({"start": p, "end": p, "radius": positive}),
    );
    let cone = shape(
        "cone",
        json!({"start": p, "end": p, "start_radius": nonnegative, "end_radius": nonnegative}),
    );
    let frame = json!({"type": "string", "description": "A frame name; `global` when absent"});
    let datum = input(
        json!({
            "kind": {"const": "datum"}, "name": name,
            "shape": {"oneOf": [point_shape, segment, plane, cuboid, slab, sphere, cylinder, cone]},
            "size": positive,
            "distance": {"type": "number", "minimum": 0, "description": "The size holds within this distance; default 0"},
            "growth": {"type": "number", "exclusiveMinimum": 1, "description": "Size growth per unit distance is growth - 1, capped at gradation - 1; default gradation"},
            "frame": frame
        }),
        &["distance", "growth", "frame"],
    );
    let region = input(
        json!({
            "kind": {"const": "region"}, "name": name,
            "shape": {"oneOf": [cuboid, slab, cylinder]},
            "size": positive, "frame": frame
        }),
        &["frame"],
    );
    let distance = json!({"type": "number", "minimum": 0, "description": "The size holds this far beyond the source; default 0"});
    let growth = json!({"type": "number", "exclusiveMinimum": 1, "description": "Size growth per unit distance is growth - 1, capped at gradation - 1; default gradation"});
    let feature_angle = json!({"type": "number", "exclusiveMinimum": 0, "maximum": 180, "description": "Degrees; edges whose facet normals turn more are feature edges, not curvature; default 30"});
    let curvature = input(
        json!({
            "kind": {"const": "curvature"}, "name": name,
            "angle": {"type": "number", "exclusiveMinimum": 0, "maximum": 90, "description": "Degrees of wall-normal turn one cell may span: size = angle / largest principal curvature (GRD-04)"},
            "distance": distance, "growth": growth, "feature_angle": feature_angle
        }),
        &["distance", "growth", "feature_angle"],
    );
    let gap = input(
        json!({
            "kind": {"const": "gap"}, "name": name,
            "cells_across": {"type": "number", "minimum": 1, "description": "Cells across the gap between opposite walls: size = local thickness / cells_across, held across the gap (GRD-04)"},
            "distance": distance, "growth": growth,
            "wedge_angle": {"type": "number", "exclusiveMinimum": 0, "maximum": 90, "description": "Degrees; walls meeting at a smaller interior angle form a gap; default 45"}
        }),
        &["distance", "growth", "wedge_angle"],
    );
    let features = input(
        json!({
            "kind": {"const": "feature_edges"}, "name": name, "size": positive,
            "distance": distance, "growth": growth, "feature_angle": feature_angle
        }),
        &["distance", "growth", "feature_angle"],
    );
    let frame_def = object_schema(json!({"name": name, "origin": p, "x_axis": p, "y_axis": p}));
    input(
        json!({
            "min_size": positive,
            "max_size": {"type": "number", "exclusiveMinimum": 0, "description": "Also the size where no source applies"},
            "gradation": {"type": "number", "exclusiveMinimum": 1, "description": "Largest neighbour size ratio; bounds the field slope by gradation - 1"},
            "relaxation_iterations": {"type": "integer", "minimum": 0, "description": "Lloyd (CVT) iterations; default 20"},
            "frames": {"type": "array", "items": frame_def},
            "sources": {"type": "array", "items": {"oneOf": [datum, region, curvature, gap, features]}}
        }),
        &["relaxation_iterations", "frames", "sources"],
    )
}

pub(crate) fn channel_config() -> Value {
    let positive = json!({"type": "number", "exclusiveMinimum": 0});
    object_schema(json!({
        "schema_version": {"const": SCHEMA_VERSION},
        "length": positive, "width": positive, "height": positive,
        "layers": {"type":"integer", "minimum":1},
        "axial_cells": {"type":"integer", "minimum":1},
        "first_thickness": positive, "growth": positive
    }))
}

/// An object schema with every property required, unless `optional`.
pub(crate) fn input(properties: Value, optional: &[&str]) -> Value {
    let required: Vec<_> = properties
        .as_object()
        .map(|o| {
            o.keys()
                .filter(|k| !optional.contains(&k.as_str()))
                .cloned()
                .collect()
        })
        .unwrap_or_default();
    json!({
        "type": "object",
        "additionalProperties": false,
        "required": required,
        "properties": properties
    })
}

/// The output schema of a mesh tool by its legacy name.
pub(crate) fn output_schema(name: &str) -> Value {
    let version = json!({ "const": SCHEMA_VERSION });
    let string = json!({ "type": "string" });
    let count = json!({ "type": "integer", "minimum": 0 });
    let number = json!({ "type": "number" });
    let fixture = json!({ "type": ["string", "null"] });
    let seed = json!({ "type": ["integer", "null"], "minimum": 0 });
    let digest = json!({ "type": "string", "pattern": "^[a-f0-9]{64}$" });
    let point = json!({ "type": "array", "items": number, "minItems": 3, "maxItems": 3 });
    let build = object_schema(json!({
        "algorithm": string, "engine_version": string, "operating_system": string,
        "architecture": string, "debug_assertions": { "type": "boolean" }, "elapsed_micros": count,
        "input_digest": digest, "config_digest": digest, "topology_digest": digest, "mesh_digest": digest
    }));
    let check = object_schema(json!({
        "schema_version": version, "valid": { "const": true }, "vertices": count,
        "faces": count, "internal_faces": count, "boundary_faces": count, "cells": count,
        "total_volume": number, "domain_volume": number, "relative_volume_error": number,
        "conservation_tolerance": number, "geometry_tolerance": number, "topology_digest": digest
    }));
    match name {
        "inspect_domain" => object_schema(json!({
            "schema_version": version, "domain": object_schema(json!({ "min": point, "max": point })),
            "site_count": count, "fixture_id": fixture, "seed": seed, "algorithm": string,
            "limitations": { "type": "array", "items": string }
        })),
        "generate_mesh" => {
            object_schema(json!({ "schema_version": version, "artifact": string, "build": build }))
        }
        "inspect_stl_domain" => object_schema(json!({
            "schema_version": version,
            "bounds": object_schema(json!({ "min": point, "max": point })),
            "site_count": count, "facet_count": count, "vertex_count": count, "volume": number,
            "patches": { "type": "array", "items": input(json!({
                "name": string, "role": role(), "facets": count, "area": number
            }), &["role"]) },
            "fixture_id": fixture, "seed": seed, "algorithm": string
        })),
        "generate_duct" => {
            let numbers = json!({"type":"array","items":number});
            let duct = object_schema(json!({
                "schema_version":version,"algorithm":string,"mode":{"enum":["bulk","layered","polyhedral"]},
                "source_sha256":digest,"config_sha256":digest,"regions":count,"source_facets":count,"scaffold_volume":{"type":["number","null"]},"stl_patch_areas":object_schema(json!({"inlet":number,"outlet":number,"walls":number})),
                "patches":{"type":"array","items":object_schema(json!({"id":count,"name":string,"role":{"enum":["inlet","outlet","wall"]},"layered":{"type":"boolean"},"stl_area":number,"area":number}))},
                "tetrahedra":count,"polyhedra":count,"layer_polyhedra":count,"core_sphericity":{"anyOf":[object_schema(json!({"minimum":number,"p05":number,"median":number})),{"type":"null"}]},"polyhedron_face_range":{"type":["array","null"], "items":count,"minItems":2,"maxItems":2},"prism_core_faces":count,"prisms":count,"source_volume":number,"volume":number,
                "minimum_cell_volume":number,"minimum_tet_mean_ratio":{"type":["number","null"]},"tetrahedra_below_mean_ratio_0_1":count,
                "minimum_wedge_jacobian":{"type":["number","null"]},"minimum_kernel_clearance":{"type":["number","null"]},
                "patch_areas":object_schema(json!({"inlet":number,"outlet":number,"walls":number})),
                "port_areas_by_layer":object_schema(json!({"inlet":numbers,"outlet":numbers})),
                "thicknesses":numbers,"normal_thickness_ranges":{"type":"array","items":{"type":"array","minItems":2,"maxItems":2,"items":number}},
                "generation_seconds":number
            }));
            object_schema(json!({"schema_version":version,"output":string,"report":duct}))
        }
        "generate_sweep" => {
            let numbers = json!({"type":"array","items":number});
            object_schema(
                json!({"schema_version":version,"output":string,"report":object_schema(json!({
                    "schema_version":version,"algorithm":string,"config":sweep_schema(),"cells":count,"core_cells":count,"wall_layer_cells":count,
                    "axial_intervals":count,"core_axial_intervals":count,"cross_section_core_cells":count,
                    "core_cell_size_mm":object_schema(json!({"min":number,"median":number,"max":number})),
                    "neighbour_size_ratio":object_schema(json!({"p95_cell_max":number,"max":number,"limit":number})),
                    "max_internal_non_orthogonality_deg":number,"max_internal_skewness":number,"wall_thicknesses_mm":numbers,
                    "inlet_area_by_layer_mm2":numbers,"outlet_area_by_layer_mm2":{"type":["array","null"],"items":number},"volume_mm3":number,
                    "analytic_volume_mm3":number,"relative_geometry_volume_error":number,"profile_sides":count
                }))}),
            )
        }
        "generate_channel" => {
            let numbers = json!({"type":"array","items":number});
            let port = object_schema(json!({
                "faces": count, "core_area": number, "layer_areas": numbers, "total_area": number
            }));
            object_schema(json!({
                "schema_version": version, "output": string,
                "report": object_schema(json!({
                    "cells": count, "prisms": count, "core_cells": count, "total_volume": number,
                    "minimum_cell_volume": number, "wall_distances": numbers, "inlet": port, "outlet": port
                }))
            }))
        }
        "generate_stl_mesh" => input(
            json!({
                "schema_version": version, "artifact_schema_version": { "const": SURFACE_ARTIFACT_VERSION },
                "artifact": string, "build": build,
                "sizing": { "type": "object", "description": "Present when a size field placed the sites: the seeding and the achieved sizes per source (GRD-08)" }
            }),
            &["sizing"],
        ),
        "check_mesh" => check,
        "mesh_report" => {
            let cell = object_schema(json!({
                "centroid": point, "volume": number, "face_count": count, "neighbour_count": count
            }));
            let face = object_schema(json!({
                "centroid": point, "normal": point, "area": number,
                "non_orthogonality_degrees": number, "skewness": number
            }));
            let full = object_schema(json!({
                "schema_version": version, "fixture_id": fixture, "seed": seed,
                "build": build, "validation": check,
                "cells": { "type": "array", "items": cell },
                "faces": { "type": "array", "items": face }
            }));
            let maybe = json!({ "type": ["number", "null"] });
            let paged = object_schema(json!({
                "schema_version": version, "fixture_id": fixture, "seed": seed,
                "build": build, "validation": check,
                "summary": object_schema(json!({
                    "cells": count, "faces": count,
                    "min_cell_volume": maybe, "max_cell_volume": maybe,
                    "max_non_orthogonality_degrees": maybe, "max_skewness": maybe
                })),
                "page": object_schema(json!({
                    "rows": { "enum": ["cells", "faces"] }, "offset": count, "count": count,
                    "items": { "type": "array", "items": { "anyOf": [cell, face] } },
                    "next_cursor": { "type": ["string", "null"] }
                }))
            }));
            json!({ "anyOf": [full, paged] })
        }
        "export_mesh" => object_schema(json!({
            "schema_version": version, "output": string, "format": string,
            "files": { "type": "array", "items": object_schema(json!({ "path": string, "blake3": digest })) }
        })),
        _ => json!({ "type": "object" }),
    }
}

/// A patch role (BC-04, FC-04).
pub(crate) fn role() -> Value {
    json!({ "enum": ["wall", "inlet", "outlet", "symmetry", "interface", "unspecified"] })
}

fn name() -> Value {
    json!({ "type": "string", "pattern": "^[A-Za-z0-9_-]{1,64}$" })
}

fn names() -> Value {
    json!({ "type": "array", "items": name() })
}

fn angle() -> Value {
    json!({ "type": "number", "minimum": 0, "maximum": 180, "description": "Degrees" })
}

/// The fields of every selection tool, with `extra` and its optional keys.
fn selection(extra: Value, optional: &[&str]) -> Value {
    let mut properties = json!({
        "surface": path(), "patches": path(), "output": path(),
        "name": name(), "role": role(), "within": names()
    });
    let mut keep = vec!["role", "within"];
    keep.extend_from_slice(optional);
    if let (Some(p), Some(e)) = (properties.as_object_mut(), extra.as_object()) {
        p.extend(e.clone());
    }
    input(properties, &keep)
}

pub(crate) fn patch_plane() -> Value {
    selection(
        json!({"point": point(), "normal": point(),
               "tolerance": {"type": "number", "minimum": 0, "description": "Largest distance of a corner from the plane; default 1e-6 of the longest extent"}}),
        &["tolerance"],
    )
}

pub(crate) fn patch_box() -> Value {
    selection(
        json!({"min": point(), "max": point(), "contain": {"enum": ["centroid", "all_corners"]}}),
        &["contain"],
    )
}

pub(crate) fn patch_normal() -> Value {
    selection(json!({"direction": point(), "max_angle_deg": angle()}), &[])
}

pub(crate) fn patch_grow() -> Value {
    selection(
        json!({"seed_facet": {"type": "integer", "minimum": 0}, "seed_point": point(),
               "feature_angle_deg": angle()}),
        &["seed_facet", "seed_point"],
    )
}

pub(crate) fn patch_split() -> Value {
    input(
        json!({"surface": path(), "patches": path(), "output": path(), "patch": name(),
               "feature_angle_deg": angle(), "names": names()}),
        &["feature_angle_deg", "names"],
    )
}

pub(crate) fn patch_merge() -> Value {
    input(
        json!({"surface": path(), "patches": path(), "output": path(),
               "merge": {"type": "array", "items": name(), "minItems": 1}, "name": name(), "role": role()}),
        &["role"],
    )
}

/// The result of a `patch.*` tool that writes a map.
pub(crate) fn patch_result() -> Value {
    let string = json!({ "type": "string" });
    let digest = json!({ "type": "string", "pattern": "^[a-f0-9]{64}$" });
    object_schema(json!({
        "format": {"const": "koflcraft.patch-result"}, "format_version": string,
        "output": string,
        "files": { "type": "array", "items": object_schema(json!({ "path": string, "blake3": digest })) },
        "change": object_schema(json!({
            "facets_changed": {"type": "integer", "minimum": 0}, "created": names(), "removed": names()
        })),
        "summary": patch_summary()
    }))
}

/// `patch.detect_ports`: a patch result plus the ranked candidates.
pub(crate) fn patch_detect_ports() -> Value {
    let mut schema = patch_result();
    let number = json!({ "type": "number" });
    let count = json!({ "type": "integer", "minimum": 0 });
    let candidate = object_schema(json!({
        "id": count, "facets": count,
        "facet_span": {"type": "array", "items": count, "minItems": 2, "maxItems": 2},
        "area": number, "centre": point(), "normal": point(), "equivalent_radius": number,
        "circularity": number, "rim_edges": count, "pipe_end_rim": number,
        "rim_angle_deg": {"type": "array", "items": number, "minItems": 2, "maxItems": 2},
        "score": number
    }));
    schema["properties"]["detection"] = object_schema(json!({
        "inlet": count, "outlet": count, "port_facets": count, "planar_regions": count,
        "disc_regions": count, "margin": number,
        "candidates": {"type": "array", "items": candidate}, "note": {"type": "string"}
    }));
    schema["required"]
        .as_array_mut()
        .unwrap()
        .push(json!("detection"));
    schema
}

/// `grading.list`: the sources and frames of a partition.
pub(crate) fn grading_listing() -> Value {
    let field = size_field();
    object_schema(json!({
        "sources": field["properties"]["sources"],
        "frames": field["properties"]["frames"],
    }))
}

/// `grading.add` and `grading.remove`: the new partition and its sources.
pub(crate) fn grading_result() -> Value {
    let string = json!({ "type": "string" });
    object_schema(json!({
        "format": {"const": "koflcraft.grading-result"}, "format_version": string,
        "output": string,
        "files": { "type": "array", "items": object_schema(json!({ "path": string, "blake3": {"type": "string", "pattern": "^[a-f0-9]{64}$"} })) },
        "grading": grading_listing()
    }))
}

/// A zero-thickness shell of a duct partition (SHL-01).
pub(crate) fn shell_config() -> Value {
    let point = json!({"type": "array", "items": {"type": "number"}, "minItems": 3, "maxItems": 3});
    let mut schema = input(
        json!({
            "name": {"type": "string", "minLength": 1, "description": "The shell's name; its sides become the wall patches `<name>` (the side its normal points into) and `<name>_back`"},
            "polygon": {"type": "array", "items": point, "minItems": 3, "description": "A planar convex polygon, corners counter-clockwise about the `name` side; it may reach past the wall (the part outside the fluid is ignored)"},
            "size": {"type": "number", "exclusiveMinimum": 0, "description": "Cell size on the shell; default: the build's max_size or an eighth of the shell's narrow width, whichever is smaller"},
            "edge_size": {"type": "number", "exclusiveMinimum": 0, "description": "Cell size along the free edges; default: half the shell's size"},
        }),
        &["size", "edge_size"],
    );
    schema["description"] = json!("A zero-thickness shell (baffle, plate, fin)");
    schema
}

/// `shell.list`: the shells of a partition with their sizes.
pub(crate) fn shell_listing() -> Value {
    let size = json!({"type": ["number", "null"]});
    let string = json!({"type": "string"});
    let point = json!({"type": "array", "items": {"type": "number"}, "minItems": 3, "maxItems": 3});
    object_schema(json!({
        "shells": {"type": "array", "items": object_schema(json!({
            "name": string,
            "sides": object_schema(json!({"front": string, "back": string})),
            "polygon": {"type": "array", "items": point},
            "normal": point,
            "polygon_area": {"type": "number"},
            "size": size, "size_from": string,
            "edge_size": size, "edge_size_from": string,
        }))}
    }))
}

/// `shell.add` and `shell.remove`: the new partition and its shells.
pub(crate) fn shell_result() -> Value {
    let string = json!({ "type": "string" });
    object_schema(json!({
        "format": {"const": "koflcraft.shell-result"}, "format_version": string,
        "output": string,
        "files": { "type": "array", "items": object_schema(json!({ "path": string, "blake3": {"type": "string", "pattern": "^[a-f0-9]{64}$"} })) },
        "shells": shell_listing()["properties"]["shells"]
    }))
}

/// `patch.propose`: a patch result plus the classified regions.
pub(crate) fn patch_propose() -> Value {
    let mut schema = patch_result();
    let number = json!({ "type": "number" });
    let count = json!({ "type": "integer", "minimum": 0 });
    let mut region = object_schema(json!({
        "id": count, "patch": name(), "kind": {"enum": ["plane", "cylinder", "sphere", "freeform"]},
        "suggestion": {"enum": ["opening", "box_side", "wall"]},
        "facets": count, "area": number, "centre": point(), "normal": point(), "axis": point(),
        "axis_point": point(), "radius": number, "length": number, "fit_error": number,
        "box_side": {"enum": ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]},
        "opening_score": number,
        "bounds": {"type": "array", "items": point(), "minItems": 2, "maxItems": 2},
        "neighbours": {"type": "array", "items": count}
    }));
    region["required"] = json!([
        "id",
        "patch",
        "kind",
        "suggestion",
        "facets",
        "area",
        "centre",
        "fit_error",
        "bounds",
        "neighbours"
    ]);
    schema["properties"]["proposal"] = object_schema(json!({
        "feature_angle_deg": number,
        "regions": {"type": "array", "items": region},
        "needs_role": names(), "next": {"type": "string"}
    }));
    schema["required"]
        .as_array_mut()
        .unwrap()
        .push(json!("proposal"));
    schema
}

/// `bc.skeleton`: the written case skeleton and what it leaves open.
pub(crate) fn bc_skeleton() -> Value {
    let string = json!({ "type": "string" });
    let patch = object_schema(json!({
        "name": name(), "role": role(), "foehncraft_type": {"type": ["string", "null"]},
        "todo": {"type": "array", "items": string}
    }));
    object_schema(json!({
        "format": {"const": "koflcraft.bc-skeleton"}, "format_version": string,
        "case_format": string, "mesh": string, "output": string,
        "files": { "type": "array", "items": object_schema(json!({ "path": string, "blake3": {"type": "string", "pattern": "^[a-f0-9]{64}$"} })) },
        "patches": {"type": "array", "items": patch},
        "needs_role": names()
    }))
}

pub(crate) fn duct_setup_result() -> Value {
    let number = json!({ "type": "number" });
    let string = json!({ "type": "string" });
    object_schema(json!({
        "format": {"const": "koflcraft.duct-setup"}, "format_version": string,
        "output": string,
        "files": { "type": "array", "items": object_schema(json!({ "path": string, "blake3": {"type": "string"} })) },
        "port_radius": number, "size_factor": number,
        "partition": {"type": "object"}, "note": string
    }))
}

/// `koflcraft.patch-summary`: evidence per patch and the BC-06 check.
pub(crate) fn patch_summary() -> Value {
    let count = json!({ "type": "integer", "minimum": 0 });
    let number = json!({ "type": "number" });
    let point = point();
    let patch = input(
        json!({
            "name": name(), "role": role(), "foehncraft_type": {"type": ["string", "null"]},
            "source_name": {"type": "string"}, "facets": count, "area": number,
            "normal": {"anyOf": [point, {"type": "null"}]}, "max_normal_deviation_deg": number,
            "bounds": {"anyOf": [{"type": "array", "items": point, "minItems": 2, "maxItems": 2}, {"type": "null"}]},
            "regions": count, "neighbours": names()
        }),
        &["source_name"],
    );
    object_schema(json!({
        "format": {"const": "koflcraft.patch-summary"}, "format_version": {"type": "string"},
        "source_sha256": {"type": "string", "pattern": "^[a-f0-9]{64}$"},
        "facet_count": count, "area": number, "every_facet_in_one_patch": {"const": true},
        "patches": {"type": "array", "items": patch}
    }))
}

pub(crate) fn object_schema(properties: Value) -> Value {
    let required: Vec<_> = properties
        .as_object()
        .map(|o| o.keys().cloned().collect())
        .unwrap_or_default();
    json!({
        "type": "object",
        "additionalProperties": false,
        "required": required,
        "properties": properties
    })
}

/// The input of `block.template` (BS-03).
pub(crate) fn block_template() -> Value {
    let positive = json!({"type": "number", "exclusiveMinimum": 0});
    let count = json!({"type": "integer", "minimum": 1});
    let channel = input(
        json!({"template": {"const": "channel"}, "length": positive, "width": positive,
               "height": positive, "cells": {"type": "array", "items": count, "minItems": 3, "maxItems": 3},
               "wall_spacing": positive}),
        &["wall_spacing"],
    );
    let pipe = input(
        json!({"template": {"const": "pipe_ogrid"}, "radius": positive, "length": positive,
               "core_fraction": {"type": "number", "exclusiveMinimum": 0.1, "exclusiveMaximum": 0.9},
               "core_cells": count, "radial_cells": count, "axial_cells": count,
               "wall_spacing": positive}),
        &["core_fraction", "wall_spacing"],
    );
    let elbow = input(
        json!({"template": {"const": "elbow"}, "radius": positive, "bend_radius": positive,
               "angle": {"type": "number", "exclusiveMinimum": 0, "maximum": 180},
               "inlet_length": {"type": "number", "minimum": 0},
               "outlet_length": {"type": "number", "minimum": 0},
               "core_fraction": {"type": "number", "exclusiveMinimum": 0.1, "exclusiveMaximum": 0.9},
               "core_cells": count, "radial_cells": count, "axial_cells": count,
               "wall_spacing": positive}),
        &[
            "inlet_length",
            "outlet_length",
            "core_fraction",
            "wall_spacing",
        ],
    );
    let airfoil = input(
        json!({"template": {"const": "airfoil_cgrid"},
               "naca": {"type": "string", "pattern": "^[0-9]{4}$"},
               "coordinates": {"type": "string"},
               "chord": positive, "span": positive, "farfield": positive, "wake_length": positive,
               "surface_cells": {"type": "integer", "minimum": 2}, "normal_cells": count,
               "wake_cells": count, "span_cells": count, "wall_spacing": positive,
               "leading_edge_spacing": positive, "trailing_edge_spacing": positive}),
        &[
            "naca",
            "coordinates",
            "chord",
            "span_cells",
            "leading_edge_spacing",
            "trailing_edge_spacing",
        ],
    );
    let naca = json!({"type": "string", "pattern": "^[0-9]{4}$"});
    let main_section = input(
        json!({"naca": naca, "coordinates": {"type": "string"}}),
        &["naca", "coordinates"],
    );
    let element = input(
        json!({"naca": naca, "coordinates": {"type": "string"},
               "chord": positive,
               "deflection": {"type": "number", "minimum": -80, "maximum": 80,
                              "description": "Degrees: slat nose down, flap trailing edge down"},
               "gap": positive,
               "overlap": {"type": "number", "exclusiveMaximum": 0,
                           "description": "Over the main chord; negative: the elements follow each other along x"}}),
        &["naca", "coordinates"],
    );
    let multi = input(
        json!({"template": {"const": "multi_element_airfoil"},
               "main": main_section, "slat": element, "flap": element,
               "chord": positive, "span": positive, "farfield": positive, "wake_length": positive,
               "surface_cells": {"type": "integer", "minimum": 8}, "normal_cells": count,
               "wake_cells": count, "span_cells": count, "thickness_cells": count,
               "gap_cells": count, "upstream_cells": count, "wall_spacing": positive,
               "trailing_edge_thickness": {"type": "number", "exclusiveMinimum": 0, "maximum": 0.05}}),
        &[
            "slat",
            "flap",
            "chord",
            "span_cells",
            "thickness_cells",
            "gap_cells",
            "upstream_cells",
            "trailing_edge_thickness",
        ],
    );
    let corner_counts = |n: usize| {
        json!({"type": "array", "items": count, "minItems": n, "maxItems": n,
               "description": "Cells per corner: the side from corner i to corner j gets cells[i] + cells[j]"})
    };
    let xy = json!({"type": "array", "items": {"type": "number"}, "minItems": 2, "maxItems": 2});
    let triangle = input(
        json!({"template": {"const": "triangle"},
               "corners": {"type": "array", "items": xy, "minItems": 3, "maxItems": 3},
               "cells": corner_counts(3), "span": positive, "span_cells": count}),
        &["corners", "span_cells"],
    );
    let tetrahedron = input(
        json!({"template": {"const": "tetrahedron"},
               "corners": {"type": "array", "items": point(), "minItems": 4, "maxItems": 4},
               "cells": corner_counts(4)}),
        &["corners"],
    );
    json!({"oneOf": [channel, pipe, elbow, airfoil, multi, triangle, tetrahedron]})
}

/// The volume algorithm of `block.grid` (BS-09).
pub(crate) fn block_algorithm() -> Value {
    let tfi = input(json!({"kind": {"const": "tfi"}}), &[]);
    let elliptic = input(
        json!({"kind": {"const": "elliptic"}, "wall_orthogonality": {"type": "boolean"},
               "wall_sliding": {"type": "boolean"},
               "max_iterations": {"type": "integer", "minimum": 1}}),
        &["wall_orthogonality", "wall_sliding", "max_iterations"],
    );
    json!({"oneOf": [tfi, elliptic]})
}

/// A block edge by its two vertex indices.
fn block_edge() -> Value {
    json!({"type": "array", "items": {"type": "integer", "minimum": 0}, "minItems": 2, "maxItems": 2,
           "description": "The edge's two vertex indices; laws and spline points run from the first to the second"})
}

/// An edge law of the topology file (BS-06).
pub(crate) fn block_distribution() -> Value {
    let positive = json!({"type": "number", "exclusiveMinimum": 0});
    json!({"oneOf": [
        input(json!({"law": {"const": "uniform"}}), &[]),
        input(json!({"law": {"const": "geometric"}, "first": positive, "last": positive}), &["first", "last"]),
        input(json!({"law": {"const": "tanh"}, "first": positive, "last": positive}), &[]),
    ]})
}

/// An edge curve of the topology file.
pub(crate) fn block_curve() -> Value {
    json!({"oneOf": [
        input(json!({"kind": {"const": "line"}}), &[]),
        input(json!({"kind": {"const": "arc"}, "through": point()}), &[]),
        input(json!({"kind": {"const": "spline"}, "through": {"type": "array", "items": point(), "minItems": 1}}), &[]),
    ]})
}

pub(crate) fn block_set_count() -> Value {
    input(
        json!({"topology": path(), "output": path(), "edge": block_edge(),
               "cells": {"type": "integer", "minimum": 1}}),
        &[],
    )
}

pub(crate) fn block_set_distribution() -> Value {
    input(
        json!({"topology": path(), "output": path(), "edge": block_edge(),
               "distribution": block_distribution()}),
        &[],
    )
}

pub(crate) fn block_set_edge_curve() -> Value {
    input(
        json!({"topology": path(), "output": path(), "edge": block_edge(), "curve": block_curve()}),
        &[],
    )
}

pub(crate) fn block_move_vertex() -> Value {
    input(
        json!({"topology": path(), "output": path(), "vertex": {"type": "integer", "minimum": 0},
               "to": point(), "by": point()}),
        &["to", "by"],
    )
}

pub(crate) fn block_set_patch() -> Value {
    input(
        json!({"topology": path(), "output": path(), "name": {"type": "string"},
               "rename": name(), "role": role()}),
        &["rename", "role"],
    )
}

pub(crate) fn block_edit_result() -> Value {
    json!({"type": "object", "description": "`koflcraft.block-result` with the written file, its check and the `change` (blocks, edges, vertices, patches, reference_removed)"})
}

fn index() -> Value {
    json!({"type": "integer", "minimum": 0})
}

pub(crate) fn block_add_vertex() -> Value {
    input(
        json!({"topology": path(), "output": path(), "point": point()}),
        &[],
    )
}

pub(crate) fn block_remove_vertex() -> Value {
    input(
        json!({"topology": path(), "output": path(), "vertex": index()}),
        &[],
    )
}

pub(crate) fn block_add_block() -> Value {
    let count = json!({"type": "integer", "minimum": 1});
    input(
        json!({"topology": path(), "output": path(),
               "vertices": {"type": "array", "items": index(), "minItems": 8, "maxItems": 8,
                            "description": "Corners in the local frame: (0,0,0) (1,0,0) (1,1,0) (0,1,0), then the same at k = 1; right-handed"},
               "cells": {"type": "array", "items": count, "minItems": 3, "maxItems": 3},
               "domain": {"type": "string"}, "patch": name(), "role": role()}),
        &["domain", "patch", "role"],
    )
}

pub(crate) fn block_remove_block() -> Value {
    input(
        json!({"topology": path(), "output": path(), "block": index(), "patch": name(), "role": role()}),
        &["patch", "role"],
    )
}

pub(crate) fn block_split_triangle() -> Value {
    input(
        json!({"topology": path(), "output": path(), "block": index(),
               "vertex": {"type": "integer", "minimum": 0, "description": "The block's flat corner (at least 120°)"},
               "along": {"type": "integer", "minimum": 0, "description": "3D: the other end of the flat block edge at `vertex`; the flattest edge there when absent"},
               "cells": {"type": "integer", "minimum": 1, "description": "Cells of the new row; from the triangle's spacing when absent"}}),
        &["along", "cells"],
    )
}

pub(crate) fn block_optimize_layout() -> Value {
    input(
        json!({"topology": path(), "output": path(),
               "triangle_splits": {"type": "boolean", "description": "Split blocks with a corner of 150° or more into three where that measurably helps (default true)"}}),
        &["triangle_splits"],
    )
}

pub(crate) fn block_split_block() -> Value {
    input(
        json!({"topology": path(), "output": path(), "block": index(),
               "axis": {"type": "integer", "minimum": 0, "maximum": 2},
               "at": {"type": "integer", "minimum": 1, "description": "Split after this many cells along the axis"}}),
        &[],
    )
}
