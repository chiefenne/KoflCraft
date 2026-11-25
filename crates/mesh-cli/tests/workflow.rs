use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use serde_json::Value;

#[test]
fn generate_check_report_and_export_round_trip() {
    let work = WorkDirectory::new();
    let inspected = work.run(&["inspect", "input.json"]);
    assert_eq!(inspected["site_count"], 4);
    assert_eq!(inspected["fixture_id"], "unit-box-4");

    let generated = work.run(&["generate", "input.json", "mesh.json"]);
    let repeat = work.run(&["generate", "input.json", "repeat.json"]);
    let topology = &generated["build"]["topology_digest"];
    assert_eq!(topology.as_str().unwrap().len(), 64);
    assert_eq!(topology, &repeat["build"]["topology_digest"]);
    assert_eq!(
        generated["build"]["mesh_digest"],
        repeat["build"]["mesh_digest"]
    );

    let checked = work.run(&["check", "mesh.json"]);
    assert_eq!(checked["valid"], true);
    assert_eq!(checked["cells"], 4);
    assert_eq!(checked["vertices"], 21);
    assert_eq!(&checked["topology_digest"], topology);
    assert!((checked["total_volume"].as_f64().unwrap() - 1.0).abs() < 1e-12);

    let report = work.run(&["report", "mesh.json"]);
    assert_eq!(report["validation"], checked);
    assert_eq!(report["build"], generated["build"]);
    assert_eq!(report["cells"].as_array().unwrap().len(), 4);
    assert_eq!(
        report["faces"].as_array().unwrap().len() as u64,
        checked["faces"].as_u64().unwrap()
    );
    assert_eq!(report["fixture_id"], "unit-box-4");

    work.run(&["export", "mesh.json", "vtu", "mesh.vtu"]);
    assert!(
        fs::read_to_string(work.path().join("mesh.vtu"))
            .unwrap()
            .contains("NumberOfCells=\"4\"")
    );
    work.run(&["export", "mesh.json", "gmsh", "mesh.msh"]);
    let msh = fs::read_to_string(work.path().join("mesh.msh")).unwrap();
    assert!(msh.starts_with("$MeshFormat\n4.1 0 8\n"));
    assert!(msh.contains("3 7 \"fluid\""));
}

#[test]
fn corrupt_artifacts_and_existing_destinations_fail() {
    let work = WorkDirectory::new();
    work.run(&["generate", "input.json", "mesh.json"]);
    let path = work.path().join("mesh.json");
    let original = fs::read(&path).unwrap();
    assert_eq!(
        work.fail(&["generate", "input.json", "mesh.json"])["code"],
        "io_error"
    );
    assert_eq!(fs::read(&path).unwrap(), original);

    let mut artifact: Value = serde_json::from_slice(&original).unwrap();
    artifact["mesh"]["vertices"][0][0] = Value::from(0.123);
    fs::write(&path, serde_json::to_vec(&artifact).unwrap()).unwrap();
    assert_eq!(
        work.fail(&["check", "mesh.json"])["code"],
        "digest_mismatch"
    );
    assert_eq!(
        work.fail(&["report", "mesh.json"])["code"],
        "digest_mismatch"
    );
    assert_eq!(
        work.fail(&["export", "mesh.json", "vtu", "bad.vtu"])["code"],
        "digest_mismatch"
    );
    assert!(!work.path().join("bad.vtu").exists());
}

#[test]
fn saved_stl_artifact_checks_and_exports_from_embedded_geometry() {
    let work = WorkDirectory::new();
    let source = work.path().join("box.stl");
    fs::write(
        &source,
        include_bytes!("../../../test-fixtures/unit-box.stl"),
    )
    .unwrap();
    fs::write(
        work.path().join("stl.json"),
        include_bytes!("../../../test-fixtures/stl-sites.json"),
    )
    .unwrap();
    let inspected = work.run(&["inspect-stl", "box.stl", "stl.json"]);
    assert_eq!(inspected["facet_count"], 12);
    assert_eq!(inspected["patches"][0]["name"], "west inlet");
    let generated = work.run(&["generate-stl", "box.stl", "stl.json", "stl.mesh.json"]);
    assert_eq!(generated["artifact_schema_version"], 2);
    fs::remove_file(source).unwrap();
    let checked = work.run(&["check", "stl.mesh.json"]);
    assert_eq!(checked["valid"], true);
    assert_eq!(checked["cells"], 4);
    work.run(&["export", "stl.mesh.json", "vtu", "stl.vtu"]);
    work.run(&["export", "stl.mesh.json", "gmsh", "stl.msh"]);
    let msh = fs::read_to_string(work.path().join("stl.msh")).unwrap();
    assert!(msh.contains("2 1 \"west inlet\""));

    let path = work.path().join("stl.mesh.json");
    let mut artifact: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    artifact["surface"]["vertices"][0][0] = Value::from(0.02);
    fs::write(path, serde_json::to_vec(&artifact).unwrap()).unwrap();
    assert_eq!(
        work.fail(&["check", "stl.mesh.json"])["code"],
        "digest_mismatch"
    );
}

struct WorkDirectory(PathBuf);

impl WorkDirectory {
    fn new() -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "koflcraft-cli-{}-{stamp}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&dir).unwrap();
        let work = Self(dir);
        fs::write(
            work.path().join("input.json"),
            include_bytes!("../../../test-fixtures/box-4.json"),
        )
        .unwrap();
        work
    }

    fn path(&self) -> &Path {
        &self.0
    }

    fn run(&self, args: &[&str]) -> Value {
        let output = Command::new(env!("CARGO_BIN_EXE_mesh-cli"))
            .args(args)
            .current_dir(self.path())
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        serde_json::from_slice(&output.stdout).unwrap()
    }

    fn fail(&self, args: &[&str]) -> Value {
        let output = Command::new(env!("CARGO_BIN_EXE_mesh-cli"))
            .args(args)
            .current_dir(self.path())
            .output()
            .unwrap();
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        serde_json::from_slice::<Value>(&output.stderr).unwrap()["error"].clone()
    }
}

impl Drop for WorkDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(self.path());
    }
}

fn vtu_array(text: &str, name: &str) -> Vec<i64> {
    let array = text.split_once(&format!("Name=\"{name}\"")).unwrap().1;
    array
        .split_once('>')
        .unwrap()
        .1
        .split_once('<')
        .unwrap()
        .0
        .split_whitespace()
        .map(|s| s.parse().unwrap())
        .collect()
}

// A type-42 count alone also accepts cells spanning the length of a pipe.
// Read the actual VTU coordinates and bound each core cell's spatial extent.
fn maximum_core_span(vtu: &str, types: &[i64]) -> f64 {
    let points = vtu
        .split_once("<Points>")
        .unwrap()
        .1
        .split_once("<DataArray")
        .unwrap()
        .1;
    let xyz: Vec<f64> = points
        .split_once('>')
        .unwrap()
        .1
        .split_once('<')
        .unwrap()
        .0
        .split_whitespace()
        .map(|v| v.parse().unwrap())
        .collect();
    let connectivity = vtu_array(vtu, "connectivity");
    let mut start = 0;
    let mut maximum = 0.0_f64;
    for (&kind, end) in types.iter().zip(vtu_array(vtu, "offsets")) {
        let end = end as usize;
        if kind == 42 {
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            for &node in &connectivity[start..end] {
                for axis in 0..3 {
                    let value = xyz[3 * node as usize + axis];
                    lo[axis] = lo[axis].min(value);
                    hi[axis] = hi[axis].max(value);
                }
            }
            let squared = (0..3).map(|i| (hi[i] - lo[i]).powi(2)).sum::<f64>();
            maximum = maximum.max(squared.sqrt());
        }
        start = end;
    }
    maximum
}

#[test]
fn channel_exports_both_full_ports_and_rejects_colliding_layers() {
    let work = WorkDirectory::new();
    let config = include_bytes!("../../../test-fixtures/straight-channel.json");
    let path = work.path().join("duct.json");
    fs::write(&path, config).unwrap();
    let result = work.run(&["generate-channel", "duct.json", "duct"]);
    let report = &result["report"];
    assert_eq!(report["prisms"], 768);
    assert_eq!(report["core_cells"], 8);
    assert!(report["minimum_cell_volume"].as_f64().unwrap() > 0.0);
    assert!((report["total_volume"].as_f64().unwrap() - 8.0).abs() < 1e-12);
    let bundle = work.path().join("duct");
    let volume = fs::read_to_string(bundle.join("volume.vtu")).unwrap();
    let types = vtu_array(&volume, "types");
    let offsets = vtu_array(&volume, "faceoffsets");
    let layers = vtu_array(&volume, "LayerIndex");
    assert_eq!(types.len(), 776);
    assert_eq!(offsets.len(), types.len());
    assert_eq!(layers.len(), types.len());
    for i in 0..types.len() {
        assert_eq!(
            (types[i], offsets[i] > 0),
            if layers[i] == 0 {
                (42, true)
            } else {
                (13, false)
            }
        );
    }

    for (patch, id) in [("inlet", 1), ("outlet", 2)] {
        let data = fs::read_to_string(bundle.join(format!("{patch}.vtu"))).unwrap();
        assert_eq!(vtu_array(&data, "PatchId"), vec![id; 85]);
        let owners = vtu_array(&data, "OwnerCellId");
        let patch_layers = vtu_array(&data, "LayerIndex");
        assert_eq!(owners.len(), 85);
        assert_eq!(patch_layers.len(), 85);
        assert!(
            owners
                .iter()
                .zip(&patch_layers)
                .all(|(&owner, &layer)| layers[owner as usize] == layer)
        );
        let mut face_counts = [0; 4];
        for layer in patch_layers {
            face_counts[layer as usize] += 1;
        }
        assert_eq!(face_counts, [1, 44, 28, 12]);
        assert_eq!(report[patch]["faces"], 85);
        assert!((report[patch]["core_area"].as_f64().unwrap() - 1.0044).abs() < 1e-12);
        for (i, expect) in [0.2336, 0.3264, 0.4356].into_iter().enumerate() {
            assert!((report[patch]["layer_areas"][i].as_f64().unwrap() - expect).abs() < 1e-12);
        }
    }
    let mut collision: Value = serde_json::from_slice(config).unwrap();
    collision["first_thickness"] = Value::from(0.5);
    fs::write(&path, serde_json::to_vec(&collision).unwrap()).unwrap();
    assert_eq!(
        work.fail(&["generate-channel", "duct.json", "collision"])["code"],
        "invalid_domain"
    );
    assert!(!work.path().join("collision").exists());
}

fn blake3(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}

/// Invariant: a bundle hands over the blocked binary `foehncraft.vtkhdf`
/// by default and `foehncraft.vtu` with `--handoff vtu`, never both, each
/// listed in `mesh-evidence.json` (format 1.1, `file_hash` blake3) with
/// its BLAKE3 hash and the same summary;
/// `export … foehncraft-vtkhdf` writes it with the artifact's digests.
/// Breaks if: the default or the option is lost between the command line
/// and the engine, the evidence misses the file, or the export drops its
/// provenance. No other test: the mesh-io tests write the file without a
/// tool.
#[test]
fn vtkhdf_handoff_is_the_default_and_vtu_an_option() {
    let work = WorkDirectory::new();
    fs::write(
        work.path().join("duct.json"),
        include_bytes!("../../../test-fixtures/straight-channel.json"),
    )
    .unwrap();
    work.run(&["generate-channel", "duct.json", "ascii", "--handoff", "vtu"]);
    work.run(&["generate-channel", "duct.json", "binary"]);
    let evidence = |bundle: &str| -> Value {
        let path = work.path().join(bundle).join("mesh-evidence.json");
        serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
    };
    let (ascii, binary) = (evidence("ascii"), evidence("binary"));
    for e in [&ascii, &binary] {
        assert_eq!(e["format_version"], "1.1");
        assert_eq!(e["file_hash"], "blake3");
    }
    assert_eq!(ascii["summary"], binary["summary"]);
    assert!(!work.path().join("binary/foehncraft.vtu").exists());
    assert!(!work.path().join("ascii/foehncraft.vtkhdf").exists());
    let vtu = fs::read(work.path().join("ascii/foehncraft.vtu")).unwrap();
    assert!(
        ascii["files"]
            .as_array()
            .unwrap()
            .iter()
            .any(|f| f[0] == "foehncraft.vtu" && f[1] == blake3(&vtu).as_str())
    );
    let bytes = fs::read(work.path().join("binary/foehncraft.vtkhdf")).unwrap();
    assert_eq!(&bytes[..8], b"\x89HDF\r\n\x1a\n");
    let listed = binary["files"].as_array().unwrap();
    assert!(
        listed
            .iter()
            .any(|f| f[0] == "foehncraft.vtkhdf" && f[1] == blake3(&bytes).as_str()),
        "{listed:?}"
    );

    let build = work.run(&["generate", "input.json", "m.json"])["build"].clone();
    let export = work.run(&["export", "m.json", "foehncraft-vtkhdf", "m.vtkhdf"]);
    let bytes = fs::read(work.path().join("m.vtkhdf")).unwrap();
    assert_eq!(export["files"][0]["blake3"], blake3(&bytes).as_str());
    let digest = build["mesh_digest"].as_str().unwrap().as_bytes();
    assert!(bytes.windows(digest.len()).any(|w| w == digest));
}

/// Keep rule (SW-18). Invariant: `"core": "voronoi"` stays selectable with
/// both wall joins, and the solver rules (verify.rs) reject what breaks them.
#[test]
fn voronoi_core_wall_joins_against_the_solver_rules() {
    // On hold since the solver rules decide (2026-11-25). Measured on the
    // mesh as exported (merged faces, exact cell centroids), the columns
    // join passes on the elbow and split caps leave 1 core cell with a face
    // pyramid behind its centroid, which the old kernel test accepted. When
    // the Voronoi core is fixed, split caps must pass too.
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test-fixtures/elbow");
    let work = WorkDirectory::new();
    for (join, failing) in [("columns", 0), ("split_cap", 1)] {
        let mut partition: Value =
            serde_json::from_slice(&fs::read(fixture.join("partition.json")).unwrap()).unwrap();
        let scaffold = &mut partition["volume_scaffold"];
        scaffold["file"] = fixture.join("volume-scaffold.msh").to_str().unwrap().into();
        scaffold["core"] = "voronoi".into();
        scaffold["wall_join"] = join.into();
        let config = work.path().join(format!("{join}.json"));
        fs::write(&config, serde_json::to_vec(&partition).unwrap()).unwrap();
        let (stl, patches) = (fixture.join("elbow.stl"), fixture.join("patches.json"));
        let args = [
            "generate-duct",
            stl.to_str().unwrap(),
            patches.to_str().unwrap(),
            config.to_str().unwrap(),
            "polyhedral",
            join,
        ];
        if failing == 0 {
            work.run(&args);
            continue;
        }
        let error = work.fail(&args);
        let message = error["message"].as_str().unwrap().to_owned();
        assert!(
            message.starts_with(&format!("{failing} of ")) && message.contains("solver rules"),
            "{join}: {message}"
        );
    }
}

/// The Voronoi core (on hold) with the columns join on Manifold passes the
/// solver rules measured on the mesh as exported (merged faces, exact cell
/// centroids; 2026-11-25), and the export finds no face turned away.
/// About 3 minutes, full run only.
#[test]
#[ignore = "full run: Manifold with the columns join (about 3 minutes)"]
fn columns_wall_join_is_valid_on_the_manifold() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test-fixtures/manifold");
    let work = WorkDirectory::new();
    let mut partition: Value =
        serde_json::from_slice(&fs::read(fixture.join("partition.json")).unwrap()).unwrap();
    partition["volume_scaffold"]
        .as_object_mut()
        .unwrap()
        .remove("wall_join");
    partition["volume_scaffold"]["core"] = "voronoi".into();
    let config = work.path().join("columns.json");
    fs::write(&config, serde_json::to_vec(&partition).unwrap()).unwrap();
    let (stl, patches) = (fixture.join("manifold.stl"), fixture.join("patches.json"));
    let args = [
        "generate-duct",
        stl.to_str().unwrap(),
        patches.to_str().unwrap(),
        config.to_str().unwrap(),
        "polyhedral",
        "columns",
    ];
    let report = work.run(&args)["report"].clone();
    assert_eq!(report["voronoi_core"]["wall_join"], "columns");
    let evidence: Value =
        serde_json::from_slice(&fs::read(work.path().join("columns/mesh-evidence.json")).unwrap())
            .unwrap();
    assert_eq!(evidence["summary"]["inverted_faces"], 0);
}

#[test]
fn median_dual_core_stays_available_and_valid_on_the_elbow() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test-fixtures/elbow");
    let work = WorkDirectory::new();
    let mut partition: Value =
        serde_json::from_slice(&fs::read(fixture.join("partition.json")).unwrap()).unwrap();
    let scaffold = &mut partition["volume_scaffold"];
    scaffold["file"] = fixture.join("volume-scaffold.msh").to_str().unwrap().into();
    scaffold["core"] = "median_dual".into();
    let config = work.path().join("median.json");
    fs::write(&config, serde_json::to_vec(&partition).unwrap()).unwrap();
    let (stl, patches) = (fixture.join("elbow.stl"), fixture.join("patches.json"));
    let args = [
        "generate-duct",
        stl.to_str().unwrap(),
        patches.to_str().unwrap(),
        config.to_str().unwrap(),
        "polyhedral",
        "median",
    ];
    let report = work.run(&args)["report"].clone();
    assert_eq!(report["algorithm"], "volume_median_dual_v2");
    // 237 before cells failing the solver rules were split into valid parts
    assert_eq!(report["polyhedra"], 238);
    assert_eq!(report["prisms"], 6768);
    let evidence: Value =
        serde_json::from_slice(&fs::read(work.path().join("median/mesh-evidence.json")).unwrap())
            .unwrap();
    assert_eq!(evidence["summary"]["inverted_faces"], 0);
}

#[test]
fn nonconvex_elbow_preserves_volume_and_exposes_each_port_layer() {
    nonconvex_duct_preserves_volume_and_exposes_each_port_layer("elbow");
}

/// About four minutes in a release build (scaffold, 108k Voronoi cells,
/// 116k prisms), so it runs only in the full CI run (`--include-ignored`).
#[test]
#[ignore = "full run: Manifold builds three duct meshes (about 4 minutes)"]
fn nonconvex_manifold_preserves_volume_and_exposes_each_port_layer() {
    nonconvex_duct_preserves_volume_and_exposes_each_port_layer("manifold");
}

fn nonconvex_duct_preserves_volume_and_exposes_each_port_layer(only: &str) {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test-fixtures");
    let work = WorkDirectory::new();
    for (dir, stl, volume, tets, wedges, port_areas) in [
        ("elbow", "elbow.stl", 5., 52, 120, [1., 1.]),
        (
            "manifold",
            "manifold.stl",
            5290738.125459458,
            93326,
            269616,
            [1963.2400669472445, 1963.2401020872057],
        ),
    ]
    .into_iter()
    .filter(|case| case.0 == only)
    {
        let fixture = fixtures.join(dir);
        let source = fixture.join(stl);
        let patches = fixture.join("patches.json");
        let config = fixture.join("partition.json");
        let cmd = [
            "generate-duct",
            source.to_str().unwrap(),
            patches.to_str().unwrap(),
            config.to_str().unwrap(),
        ];
        for mode in ["bulk", "layered", "polyhedral"] {
            let name = format!("{dir}-{mode}");
            let args = [cmd.as_slice(), &[mode, name.as_str()]].concat();
            let report = work.run(&args)["report"].clone();
            assert!((report["source_volume"].as_f64().unwrap() / volume - 1.).abs() < 1e-10);
            let fitted_volume = report["scaffold_volume"].as_f64().unwrap_or(volume);
            let mesh_volume = report["volume"].as_f64().unwrap();
            assert!((mesh_volume / fitted_volume - 1.).abs() < 1e-10);
            assert!((mesh_volume / volume - 1.).abs() < 0.002);
            // Polyhedral: a wall triangle covered by one Voronoi core cell
            // keeps its native wedges; the others are cut into poly-prism
            // columns under the core cells (none native on the elbow).
            // Manifold's Voronoi core repairs its reentrant edge from the
            // site layout, which rounding can shift between platforms, so its
            // counts are checked within a band.
            let wall_prisms = match (mode, dir) {
                ("bulk", _) => 0,
                ("polyhedral", _) => report["prisms"].as_i64().unwrap(),
                _ => wedges,
            };
            if (mode, dir) == ("polyhedral", "manifold") {
                // Manifold builds its scaffold in Rust; the median dual
                // splits every wall prism into the six of its fine wall
                // triangles (1,020,168 on 2026-11-25, 55,717 core cells).
                assert!(
                    (1_000_000..1_040_000).contains(&wall_prisms),
                    "{wall_prisms}"
                );
                // Keep: the scaffold wall follows the STL to the chord
                // tolerance both ways (the Gmsh scaffold was 0.54 off).
                // Building the wall at a fixed coarse size, dropping the
                // projection or the deviation pass breaks this; only this
                // compares the wall of a built mesh with its STL.
                let deviation = &report["scaffold"]["wall_deviation"];
                for side in ["max_scaffold_to_stl", "max_stl_to_scaffold"] {
                    assert!(deviation[side].as_f64().unwrap() < 0.06, "{side}");
                }
                assert_eq!(report["scaffold"]["source"], "rust_remesh_delaunay_v1");
            }
            assert!(report["minimum_cell_volume"].as_f64().unwrap() > 0.);
            let expected_tets = if mode == "polyhedral" { 0 } else { tets };
            assert_eq!(report["tetrahedra"], expected_tets);
            if mode == "polyhedral" {
                // The median dual is the default core since the solver
                // rules decide (2026-11-25); the Voronoi core is on hold
                // (`voronoi_core_wall_joins_against_the_solver_rules`).
                assert_eq!(report["algorithm"], "volume_median_dual_v2");
                let polyhedra = report["polyhedra"].as_i64().unwrap();
                let expected = if dir == "elbow" {
                    // Keep: 237 before the cells that break the solver rules
                    // were split into valid parts (one on the elbow).
                    assert_eq!(polyhedra, 238);
                    assert!(report["minimum_cell_volume"].as_f64().unwrap() > 1e-4);
                    238
                } else {
                    assert!((48_500..50_500).contains(&polyhedra), "{polyhedra}");
                    polyhedra
                };
                let shape = &report["core_sphericity"];
                // Median dual (2026-11-25): elbow 0.476, 0.566, 0.762.
                let thresholds = if dir == "elbow" {
                    [0.45, 0.55, 0.75]
                } else {
                    // Manifold 2026-11-25: 0.188, 0.651, 0.787
                    [0.15, 0.60, 0.75]
                };
                for (metric, minimum) in ["minimum", "p05", "median"].into_iter().zip(thresholds) {
                    let value = shape[metric].as_f64().unwrap();
                    assert!(value > minimum, "{name} {metric} {value}");
                }
                let vtu = fs::read_to_string(work.path().join(&name).join("volume.vtu")).unwrap();
                let types = vtu_array(&vtu, "types");
                if dir == "manifold" {
                    // median dual 2026-11-25: 25.9 (one cell), 99.9 % below 19.7
                    assert!(maximum_core_span(&vtu, &types) < 30.0);
                }
                let offsets = vtu_array(&vtu, "faceoffsets");
                assert_eq!(types.len() as i64, expected + wall_prisms);
                assert_eq!(types.len(), offsets.len());
                assert_eq!(types.iter().filter(|&&t| t == 42).count() as i64, expected);
                assert!(
                    types
                        .iter()
                        .zip(offsets)
                        .all(|(t, offset)| matches!((t, offset), (13, -1) | (42, 1..)))
                );
                // Verify face stream boundaries independently on the small reference mesh.
                if dir == "elbow" {
                    let faces = vtu_array(&vtu, "faces");
                    let mut cursor = 0;
                    for end in vtu_array(&vtu, "faceoffsets")
                        .into_iter()
                        .filter(|&n| n > 0)
                    {
                        let count = faces[cursor];
                        assert!(count > 4);
                        cursor += 1;
                        for _ in 0..count {
                            let n = faces[cursor];
                            assert!(n >= 3);
                            cursor += 1 + n as usize;
                        }
                        assert_eq!(cursor as i64, end);
                    }
                    assert_eq!(cursor, faces.len());
                }
            }
            assert_eq!(report["prisms"], wall_prisms);
            // Keep: FoehnCraft rejects a mesh with any interior face whose
            // area vector points against the owner-to-neighbour centroid
            // vector (S·d <= 0). Cells at reentrant wall edges had 10 such
            // faces on the elbow and 18 on Manifold until they were split
            // (median dual) or their sites split in pairs (Voronoi core);
            // undoing that brings them back. The shell checks above do not
            // look at centroids, so only this catches it.
            let evidence: Value = serde_json::from_slice(
                &fs::read(work.path().join(&name).join("mesh-evidence.json")).unwrap(),
            )
            .unwrap();
            assert_eq!(evidence["summary"]["inverted_faces"], 0, "{name}");
            // Keep: the hand-off fans every face about its node average, as
            // FoehnCraft does, so its volume equals the mesher's only if no
            // bent wall face was merged into one polygon. Wall columns merged
            // over curved Manifold walls lost 0.14% this way. The checks
            // above integrate KoflCraft's own faces and cannot see it.
            let handoff_volume = evidence["summary"]["total_volume"].as_f64().unwrap();
            assert!((handoff_volume / mesh_volume - 1.).abs() < 1e-10, "{name}");
            for (patch, area) in ["inlet", "outlet"].into_iter().zip(port_areas) {
                assert!((report["patch_areas"][patch].as_f64().unwrap() / area - 1.).abs() < 1e-10);
                let layers = report["port_areas_by_layer"][patch].as_array().unwrap();
                assert_eq!(layers.len(), if mode == "bulk" { 1 } else { 4 });
                assert!(layers.iter().all(|v| v.as_f64().unwrap() > 0.));
            }
            if mode != "bulk" {
                for (range, requested) in report["normal_thickness_ranges"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .zip(report["thicknesses"].as_array().unwrap())
                {
                    for limit in range.as_array().unwrap() {
                        assert!(
                            (limit.as_f64().unwrap() / requested.as_f64().unwrap() - 1.).abs()
                                < if mode == "polyhedral" && dir != "elbow" {
                                    0.125
                                } else {
                                    0.001
                                }
                        );
                    }
                }
                assert!(report["minimum_wedge_jacobian"].as_f64().unwrap() > 0.);
            }
        }
        if dir == "manifold" {
            let mut bad: Value = serde_json::from_slice(&fs::read(config).unwrap()).unwrap();
            bad["first_thickness"] = 1.0.into();
            fs::write(
                work.path().join("oversized.json"),
                serde_json::to_vec(&bad).unwrap(),
            )
            .unwrap();
            let failure = work.fail(&[
                "generate-duct",
                source.to_str().unwrap(),
                patches.to_str().unwrap(),
                "oversized.json",
                "layered",
                "collision",
            ]);
            assert_eq!(failure["code"], "validation_failed");
            let message = failure["message"].as_str().unwrap();
            assert!(message.contains("collision"), "{message}");
            assert!(!work.path().join("collision").exists());
        }
    }
}

#[test]
fn analytic_benchmarks_have_fine_graded_cores_and_named_openings() {
    let work = WorkDirectory::new();
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test-fixtures/analytic");
    for (name, volume, core_spacing, minimum_core) in [
        ("block", 100_000., 1., 80_000),
        (
            "pipe",
            std::f64::consts::PI * 25_f64.powi(2) * 150.,
            1.25,
            120_000,
        ),
    ] {
        let config = fixtures.join(format!("{name}.json"));
        let report =
            work.run(&["generate-sweep", config.to_str().unwrap(), name])["report"].clone();
        let dir = work.path().join(name);
        assert_eq!(
            serde_json::from_slice::<Value>(&fs::read(dir.join("report.json")).unwrap()).unwrap(),
            report
        );
        assert!(report["core_cells"].as_u64().unwrap() >= minimum_core);
        assert!(report["core_cell_size_mm"]["max"].as_f64().unwrap() <= core_spacing);
        assert!(report["neighbour_size_ratio"]["max"].as_f64().unwrap() < 1.15);
        assert!(
            (report["volume_mm3"].as_f64().unwrap() / volume - 1.).abs()
                < if name == "block" { 1e-10 } else { 0.001 }
        );
        let inlet = &report["inlet_area_by_layer_mm2"];
        assert_eq!(inlet.as_array().unwrap().len(), 6);
        assert!(
            inlet
                .as_array()
                .unwrap()
                .iter()
                .all(|v| v.as_f64().unwrap() > 0.)
        );
        let vtm = fs::read_to_string(dir.join("mesh.vtm")).unwrap();
        assert!(vtm.contains("inlet.vtu") && vtm.contains("walls.vtu"));
        assert_eq!(vtm.contains("outlet.vtu"), name == "pipe");
        assert_eq!(dir.join("outlet.vtu").exists(), name == "pipe");
        if name == "pipe" {
            assert_eq!(inlet, &report["outlet_area_by_layer_mm2"]);
        } else {
            assert_eq!(
                report["axial_intervals"].as_u64().unwrap()
                    - report["core_axial_intervals"].as_u64().unwrap(),
                5
            );
        }
    }
}

#[test]
fn sweep_refinement_reduces_cell_size_and_rejects_abrupt_grading() {
    let work = WorkDirectory::new();
    let mut c: Value =
        serde_json::from_str(include_str!("../../../test-fixtures/analytic/block.json")).unwrap();
    c["geometry"]["dimensions_mm"] = serde_json::json!([12., 12., 8.]);
    c["geometry"]["inlet"] = "x_min".into();
    c["wall_layers"] = 2.into();
    c["max_neighbour_size_ratio"] = 1.5.into();
    let mut medians = Vec::new();
    for (name, h) in [("coarse", 1.5), ("fine", 0.75)] {
        c["core_spacing_mm"] = h.into();
        c["first_layer_mm"] = (h * 0.55).into();
        fs::write(
            work.path().join("sweep.json"),
            serde_json::to_vec(&c).unwrap(),
        )
        .unwrap();
        let r = work.run(&["generate-sweep", "sweep.json", name]);
        medians.push(r["report"]["core_cell_size_mm"]["median"].as_f64().unwrap());
    }
    assert!(medians[1] < 0.6 * medians[0]);
    c["max_neighbour_size_ratio"] = 1.0.into();
    fs::write(
        work.path().join("sweep.json"),
        serde_json::to_vec(&c).unwrap(),
    )
    .unwrap();
    let error = work.fail(&["generate-sweep", "sweep.json", "steep"]);
    assert!(
        error["message"]
            .as_str()
            .unwrap()
            .contains("neighbour size ratio")
    );
    assert!(!work.path().join("steep").exists());
}

/// Keep rule (SW-18). Invariant: `generate-duct` reads a
/// `koflcraft.patch-map` (BC-01) as it reads a port file: the elbow's port
/// ranges written as a map, with other patch names, give the same bulk mesh.
/// Breaking the map branch of the duct input fails the run or changes the
/// counts. (Scattered port facets: `ports.rs` in mesh-engine.)
#[test]
fn generate_duct_reads_a_patch_map_like_a_port_file() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test-fixtures/elbow");
    let work = WorkDirectory::new();
    let ports: Value =
        serde_json::from_slice(&fs::read(fixture.join("patches.json")).unwrap()).unwrap();
    let count = ports["facet_count"].as_u64().unwrap() as usize;
    let mut facet_patch = vec![0; count];
    for (p, port) in ports["ports"].as_array().unwrap().iter().enumerate() {
        let r = &port["facet_range"];
        for f in r[0].as_u64().unwrap()..r[1].as_u64().unwrap() {
            facet_patch[f as usize] = p + 1;
        }
    }
    let map = serde_json::json!({
        "format": "koflcraft.patch-map", "format_version": "1.0",
        "source_sha256": ports["source_sha256"], "facet_count": count,
        "patches": [{"name": "walls", "role": "wall"}, {"name": "in", "role": "inlet"},
                    {"name": "out", "role": "outlet"}],
        "facet_patch": facet_patch,
    });
    let map_path = work.path().join("map.json");
    fs::write(&map_path, serde_json::to_vec(&map).unwrap()).unwrap();
    let (stl, partition) = (fixture.join("elbow.stl"), fixture.join("partition.json"));
    let run = |patches: &Path, out: &str| {
        work.run(&[
            "generate-duct",
            stl.to_str().unwrap(),
            patches.to_str().unwrap(),
            partition.to_str().unwrap(),
            "bulk",
            out,
        ])["report"]
            .clone()
    };
    let a = run(&fixture.join("patches.json"), "ports");
    let b = run(&map_path, "map");
    for key in ["tetrahedra", "volume", "patch_areas", "source_facets"] {
        assert_eq!(a[key], b[key], "{key}");
    }
    assert_eq!(b["patch_areas"]["inlet"], 1.0);
}

/// `layers_yield`: the furnace side walls keep their layers, and the tubes
/// and the box walls run their own layers straight through them down to the
/// side walls (a slab extruded off each side wall). The whole mesh passes
/// the solver rules and the strict checks on the full domain. About 10
/// seconds.
#[test]
fn yielding_side_walls_keep_layers_under_the_crossing_walls() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test-fixtures/furnace");
    let work = WorkDirectory::new();
    let path = |f: &str| fixture.join(f).to_str().unwrap().to_string();
    let partition = fs::read_to_string(fixture.join("partition.json"))
        .unwrap()
        .replace("\"no_layers\"", "\"layers_yield\"");
    let yielding = work.path().join("partition.json");
    fs::write(&yielding, partition).unwrap();
    let report = work.run(&[
        "generate-duct",
        &path("furnace.stl"),
        &path("furnace.patches.json"),
        yielding.to_str().unwrap(),
        "polyhedral",
        "furnace",
    ])["report"]
        .clone();
    let patches = report["patches"].as_array().unwrap();
    assert!(patches.iter().all(|p| p["layered"] == true));
    let area = |name: &str| {
        let p = patches.iter().find(|p| p["name"] == name).unwrap();
        p["area"].as_f64().unwrap() / p["stl_area"].as_f64().unwrap()
    };
    // the side walls are flat: the slab gives back their whole area
    assert!((area("sides") - 1.).abs() < 0.001, "{}", area("sides"));
    assert!((area("walls") - 1.).abs() < 1e-6, "{}", area("walls"));
    assert!((area("tubes") - 1.).abs() < 0.03, "{}", area("tubes"));
    let volume = report["volume"].as_f64().unwrap();
    assert!((volume / report["source_volume"].as_f64().unwrap() - 1.).abs() < 0.001);
    // more cells than with sides under no_layers: the side walls' layers
    let cells = report["polyhedra"].as_i64().unwrap() + report["prisms"].as_i64().unwrap();
    assert!(cells > 60_000, "{cells}");
}

/// The box furnace (FoehnCraft gallery 259) with one tube: its burner disc
/// and its outlet opening are flush with their walls, so they cannot
/// contract in their planes as a pipe-end cap does; they grow wall layers
/// like the walls. Every patch of the map keeps its name and role through
/// the remesh (the patch borders are feature lines), into the report, the
/// Gmsh physical groups and the FoehnCraft hand-off. About 5 seconds.
#[test]
fn flush_openings_grow_layers_and_every_patch_keeps_its_name() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test-fixtures/furnace");
    let work = WorkDirectory::new();
    let path = |f: &str| fixture.join(f).to_str().unwrap().to_string();
    let (stl, map, partition) = (
        path("furnace.stl"),
        path("furnace.patches.json"),
        path("partition.json"),
    );
    let report = work.run(&[
        "generate-duct",
        &stl,
        &map,
        &partition,
        "polyhedral",
        "furnace",
        "--handoff",
        "vtu",
    ])["report"]
        .clone();
    assert_eq!(report["algorithm"], "volume_median_dual_v2");
    let patches = report["patches"].as_array().unwrap();
    let named: Vec<(&str, &str, bool)> = patches
        .iter()
        .map(|p| {
            (
                p["name"].as_str().unwrap(),
                p["role"].as_str().unwrap(),
                p["layered"].as_bool().unwrap(),
            )
        })
        .collect();
    assert_eq!(
        named,
        [
            ("walls", "wall", true),
            ("front", "wall", true),
            ("tubes", "wall", true),
            ("sides", "wall", false),
            ("burner", "inlet", true),
            ("outlet", "outlet", true),
        ]
    );
    for p in patches {
        let (stl_area, area) = (p["stl_area"].as_f64().unwrap(), p["area"].as_f64().unwrap());
        // the burner disc loses the segments of its 0.04 rim chords (0.8%),
        // the tube those of its wall chords
        let limit = if p["name"] == "tubes" { 0.03 } else { 0.01 };
        assert!(
            (area / stl_area - 1.).abs() < limit,
            "{}: {area} {stl_area}",
            p["name"]
        );
    }
    // no cap: every boundary face is the wall face of a layer prism
    assert_eq!(report["patch_areas"]["inlet"], 0.0);
    assert_eq!(report["patch_areas"]["outlet"], 0.0);
    assert!(report["prisms"].as_i64().unwrap() > 0);
    let volume = report["volume"].as_f64().unwrap();
    assert!((volume / report["source_volume"].as_f64().unwrap() - 1.).abs() < 0.005);
    let out = work.path().join("furnace");
    let msh = fs::read_to_string(out.join("mesh.msh")).unwrap();
    for group in [
        "\"walls\"",
        "\"front\"",
        "\"tubes\"",
        "\"burner\"",
        "\"outlet\"",
        "\"sides\"",
        "\"fluid\"",
    ] {
        assert!(msh.contains(group), "{group}");
    }
    let vtu = fs::read_to_string(out.join("foehncraft.vtu")).unwrap();
    let names = vtu
        .lines()
        .find(|l| l.contains("Name=\"boundary_names\""))
        .unwrap();
    let text: String = names
        .split('>')
        .nth(1)
        .unwrap()
        .split_whitespace()
        .filter_map(|b| b.parse::<u8>().ok())
        .map(|b| if b == 0 { ' ' } else { b as char })
        .collect();
    assert_eq!(
        text.trim(),
        "0=walls 1=front 2=tubes 3=sides 4=burner 5=outlet"
    );
}

/// The hex core (`"core": "hex_core"`, research/hex-core.md) on the furnace:
/// octree cubes in the bulk, the median dual between them and the wall
/// layers; the mesh passes the solver rules and keeps the volume.
#[test]
fn hex_core_fills_the_furnace_bulk_with_cubes() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test-fixtures/furnace");
    let work = WorkDirectory::new();
    let mut partition: Value =
        serde_json::from_str(&fs::read_to_string(fixture.join("partition.json")).unwrap()).unwrap();
    partition["volume_scaffold"]["core"] = "hex_core".into();
    let hex = work.path().join("hex.json");
    fs::write(&hex, partition.to_string()).unwrap();
    let path = |f: &str| fixture.join(f).to_str().unwrap().to_string();
    let report = work.run(&[
        "generate-duct",
        &path("furnace.stl"),
        &path("furnace.patches.json"),
        hex.to_str().unwrap(),
        "polyhedral",
        "furnace",
    ])["report"]
        .clone();
    assert_eq!(report["algorithm"], "volume_hex_core_v1");
    let core = &report["scaffold"]["build"]["hex_core"];
    let cubes = core["cubes"].as_u64().unwrap();
    assert!(cubes > 100, "{core}");
    // graded: the wall spacing grows inward, so the cubes come in several sizes
    assert!(core["cube_sizes"].as_array().unwrap().len() > 1, "{core}");
    // the cubes are cells of the mesh
    assert!(report["polyhedra"].as_u64().unwrap() > cubes);
    let volume = report["volume"].as_f64().unwrap();
    assert!((volume / report["source_volume"].as_f64().unwrap() - 1.).abs() < 0.005);
    let walls = report["patch_areas"]["walls"].as_f64().unwrap();
    assert!(walls > 0.);
}
