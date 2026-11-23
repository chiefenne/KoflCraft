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
    assert_eq!(report["polyhedra"], 237);
    assert_eq!(report["prisms"], 6768);
    let evidence: Value =
        serde_json::from_slice(&fs::read(work.path().join("median/mesh-evidence.json")).unwrap())
            .unwrap();
    assert_eq!(evidence["summary"]["inverted_faces"], 0);
}

#[test]
fn nonconvex_elbow_and_solidcraft_preserve_volume_and_expose_each_port_layer() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test-fixtures");
    let work = WorkDirectory::new();
    for (dir, stl, volume, tets, wedges, port_areas) in [
        ("elbow", "elbow.stl", 5., 52, 120, [1., 1.]),
        (
            "solidcraft-vs-foehncraft",
            "SolidCraft_vs_FoehnCraft.stl",
            5290738.125459458,
            93326,
            269616,
            [1963.2400669472445, 1963.2401020872057],
        ),
    ] {
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
            // Polyhedral: native wedges only; the innermost wedges shared by
            // several Voronoi core cells are polyhedra with a split cap.
            // SolidCraft's Voronoi core repairs its reentrant edge from the
            // site layout, which rounding can shift between platforms, so its
            // counts are checked within a band (all wall triangles stay
            // layered: native wedges plus split caps are 3 per triangle).
            let wall_prisms = match (mode, dir) {
                ("bulk", _) => 0,
                ("polyhedral", "elbow") => 752,
                ("polyhedral", _) => report["prisms"].as_i64().unwrap(),
                _ => wedges,
            };
            if (mode, dir) == ("polyhedral", "solidcraft-vs-foehncraft") {
                assert!((43_300..43_500).contains(&wall_prisms), "{wall_prisms}");
            }
            assert!(report["minimum_cell_volume"].as_f64().unwrap() > 0.);
            let expected_tets = if mode == "polyhedral" { 0 } else { tets };
            assert_eq!(report["tetrahedra"], expected_tets);
            if mode == "polyhedral" {
                assert_eq!(report["algorithm"], "volume_clipped_voronoi_cvt_v1");
                let polyhedra = report["polyhedra"].as_i64().unwrap();
                let core_faces = report["prism_core_faces"].as_i64().unwrap();
                let expected = if dir == "elbow" {
                    assert_eq!((polyhedra, core_faces), (603, 1362));
                    603
                } else {
                    assert!((48_500..48_800).contains(&polyhedra), "{polyhedra}");
                    assert!((69_000..70_000).contains(&core_faces), "{core_faces}");
                    polyhedra
                };
                let shape = &report["core_sphericity"];
                // The median dual reached 0.48 (elbow) and 0.34 (SolidCraft).
                // SolidCraft's minimum is a cell detached at the box-to-pipe
                // reentrant edge.
                let thresholds = if dir == "elbow" {
                    [0.70, 0.75, 0.80]
                } else {
                    [0.30, 0.78, 0.85]
                };
                for (metric, minimum) in ["minimum", "p05", "median"].into_iter().zip(thresholds) {
                    assert!(shape[metric].as_f64().unwrap() > minimum, "{name} {metric}");
                }
                let vtu = fs::read_to_string(work.path().join(&name).join("volume.vtu")).unwrap();
                let types = vtu_array(&vtu, "types");
                if dir == "solidcraft-vs-foehncraft" {
                    assert!(maximum_core_span(&vtu, &types) < 25.0);
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
            // faces on the elbow and 18 on SolidCraft until they were split
            // (median dual) or their sites split in pairs (Voronoi core);
            // undoing that brings them back. The shell checks above do not
            // look at centroids, so only this catches it.
            let evidence: Value = serde_json::from_slice(
                &fs::read(work.path().join(&name).join("mesh-evidence.json")).unwrap(),
            )
            .unwrap();
            assert_eq!(evidence["summary"]["inverted_faces"], 0, "{name}");
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
        if dir == "solidcraft-vs-foehncraft" {
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
            assert!(failure["message"].as_str().unwrap().contains("collision"));
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
