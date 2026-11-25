//! The `kofl` front end of the registry (AG-12).

use std::process::Command;

use serde_json::Value;

fn kofl(args: &[&str]) -> (bool, Value, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_kofl"))
        .args(args)
        .output()
        .unwrap();
    let stdout = String::from_utf8(out.stdout).unwrap();
    (
        out.status.success(),
        serde_json::from_str(&stdout).unwrap_or(Value::Null),
        stdout,
    )
}

/// Invariant: `docs/tools.md` is what the registry generates, so the
/// published catalogue never drifts from the code (as FoehnCraft's
/// `tools.md`). Editing a tool without regenerating breaks it.
#[test]
fn tools_md_is_generated_from_the_registry() {
    let (ok, _, markdown) = kofl(&["commands", "--markdown"]);
    assert!(ok);
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/tools.md");
    let file = std::fs::read_to_string(path).unwrap_or_default();
    assert!(
        file == markdown,
        "docs/tools.md is out of date: run `cargo run -p mesh-cli --bin kofl -- commands --markdown > docs/tools.md`"
    );
}

/// Invariant: the catalogue lists every id `kofl call` accepts, and a call
/// returns the envelope with status, exit code and diagnostics; a failing
/// call exits non-zero with its registered code.
#[test]
fn commands_and_call_agree() {
    let (ok, catalogue, _) = kofl(&["commands"]);
    assert!(ok);
    assert_eq!(catalogue["format"], "koflcraft.tool-catalogue");
    let ids: Vec<&str> = catalogue["tools"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| t["id"].as_str().unwrap())
        .collect();
    for id in ["system.info", "diag.codes", "mesh.generate", "mesh.export"] {
        assert!(ids.contains(&id), "{id}");
    }
    let (ok, envelope, _) = kofl(&["call", "system.info"]);
    assert!(ok);
    assert_eq!(envelope["format"], "koflcraft.tool-result");
    assert_eq!(envelope["status"], "finished");
    assert_eq!(envelope["result"]["format"], "koflcraft.info");
    let (ok, envelope, _) = kofl(&["call", "mesh.check", "--input", "/nonexistent.json"]);
    assert!(!ok);
    assert!(envelope.is_null());
    let (ok, envelope, _) = kofl(&["call", "no.such"]);
    assert!(!ok);
    assert_eq!(envelope["diagnostics"][0]["code"], "TOOL001");
}

fn kofl_in(dir: &std::path::Path, args: &[&str]) -> (bool, Value) {
    let out = Command::new(env!("CARGO_BIN_EXE_kofl"))
        .current_dir(dir)
        .env_remove("KOFL_JOURNAL")
        .args(args)
        .output()
        .unwrap();
    (
        out.status.success(),
        serde_json::from_slice(&out.stdout).unwrap_or(Value::Null),
    )
}

fn fresh(name: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("kofl-journal-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(
        dir.join("box.json"),
        include_bytes!("../../../test-fixtures/box-4.json"),
    )
    .unwrap();
    dir
}

/// Invariant (AG-13): a journal recorded in one directory replays in a
/// fresh one with every call ending the same (status, exit code, result
/// digest), and a call whose input was changed is reported as TOOL005
/// with a failing exit. A digest that kept timings, or a replay that
/// compared only statuses, breaks it; no other test runs a journal.
#[test]
fn journal_replays_and_detects_a_changed_call() {
    let rec = fresh("rec");
    let journal = rec.join("j.jsonl");
    let j = journal.to_str().unwrap();
    for args in [
        &["--journal", j, "generate", "box.json", "m.json"][..],
        &["--journal", j, "export", "m.json", "foehncraft", "m.vtu"],
        &["--journal", j, "check", "missing.json"],
    ] {
        kofl_in(&rec, args);
    }
    let text = std::fs::read_to_string(&journal).unwrap();
    let rows: Vec<Value> = text
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    assert_eq!(rows[0]["format"], "koflcraft.journal");
    assert_eq!(rows.len(), 4);
    assert_eq!(rows[3]["status"], "error");

    let same = fresh("same");
    let (ok, envelope) = kofl_in(&same, &["run", "--script", j]);
    assert!(ok, "{envelope:#}");
    assert_eq!(envelope["result"]["format"], "koflcraft.replay");
    assert_eq!(envelope["result"]["replayed"], 3);
    assert_eq!(envelope["result"]["different"], 0);

    // change the export format of the second call
    let changed = text.replace("\"format\":\"foehncraft\"", "\"format\":\"vtu\"");
    assert_ne!(changed, text);
    let edited = rec.join("edited.jsonl");
    std::fs::write(&edited, changed).unwrap();
    let other = fresh("other");
    let (ok, envelope) = kofl_in(
        &other,
        &["run", "--script", edited.to_str().unwrap(), "--keep-going"],
    );
    assert!(!ok);
    assert_eq!(envelope["status"], "error");
    assert_eq!(envelope["diagnostics"][0]["code"], "TOOL005");
    assert_eq!(envelope["result"]["different"], 1);
    assert_eq!(envelope["result"]["replayed"], 3);
    for d in [rec, same, other] {
        let _ = std::fs::remove_dir_all(d);
    }
}

/// Invariant: `kofl mcp` is the registry's MCP server; it refuses to start
/// without roots, and with `--journal` records its calls with `via: mcp`.
/// Routing `kofl mcp` around `call_journaled` breaks it.
#[test]
fn kofl_mcp_serves_the_registry_and_journals() {
    use std::io::Write;
    use std::process::Stdio;
    let dir = fresh("mcp");
    let refused = Command::new(env!("CARGO_BIN_EXE_kofl"))
        .arg("mcp")
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(!refused.status.success());
    let journal = dir.join("mcp.jsonl");
    let mut server = Command::new(env!("CARGO_BIN_EXE_kofl"))
        .env_remove("KOFL_JOURNAL")
        .arg("--journal")
        .arg(&journal)
        .args(["mcp", "--read-root"])
        .arg(&dir)
        .arg("--write-root")
        .arg(&dir)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    {
        let mut input = server.stdin.take().unwrap();
        for line in [
            r#"{"jsonrpc":"2.0","id":0,"method":"initialize","params":{"protocolVersion":"2025-11-25","capabilities":{},"clientInfo":{"name":"t","version":"1"}}}"#,
            r#"{"jsonrpc":"2.0","method":"notifications/initialized"}"#,
            r#"{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"system.info","arguments":{}}}"#,
        ] {
            writeln!(input, "{line}").unwrap();
        }
    }
    let out = server.wait_with_output().unwrap();
    assert!(out.status.success());
    let last: Value = serde_json::from_str(
        String::from_utf8(out.stdout)
            .unwrap()
            .lines()
            .last()
            .unwrap(),
    )
    .unwrap();
    assert_eq!(last["result"]["structuredContent"]["tool"], "system.info");
    let text = std::fs::read_to_string(&journal).unwrap();
    let row: Value = serde_json::from_str(text.lines().last().unwrap()).unwrap();
    assert_eq!(row["via"], "mcp");
    assert_eq!(row["tool"], "system.info");
    let _ = std::fs::remove_dir_all(dir);
}

/// Invariant (AG-15): `kofl --check` predicts exactly the cells each
/// generator then builds: box sites, STL seeds placed by a size field,
/// channel grid, swept section, and the duct in all three modes. A generator change that alters its cell count
/// without its estimate breaks it; no other test compares the two.
#[test]
fn estimates_match_the_built_cell_counts() {
    let dir = fresh("estimate");
    let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../test-fixtures");
    for f in [
        "straight-channel.json",
        "unit-box.stl",
        "stl-size-field.json",
        "elbow/elbow.stl",
        "elbow/patches.json",
        "elbow/partition.json",
        "elbow/volume-scaffold.msh",
    ] {
        let name = std::path::Path::new(f).file_name().unwrap();
        std::fs::copy(format!("{root}/{f}"), dir.join(name)).unwrap();
    }
    std::fs::write(
        dir.join("sweep.json"),
        r#"{"schema_version":1,"geometry":{"kind":"pipe","diameter_mm":12,"length_mm":6},"core_spacing_mm":1.25,"wall_layers":2,"first_layer_mm":0.55,"layer_growth":1.2,"max_neighbour_size_ratio":1.4}"#,
    )
    .unwrap();
    let built = |v: &Value| -> u64 {
        if let Some(seeds) = v["sizing"]["seeds"].as_u64() {
            return seeds;
        }
        let r = &v["report"];
        if r["cells"].is_u64() {
            return r["cells"].as_u64().unwrap();
        }
        ["tetrahedra", "prisms", "polyhedra"]
            .iter()
            .map(|k| r[k].as_u64().unwrap_or(0))
            .sum()
    };
    let cases: [(&[&str], Option<u64>); 7] = [
        (&["generate", "box.json", "m.json"], Some(4)),
        (
            &[
                "generate-stl",
                "unit-box.stl",
                "stl-size-field.json",
                "sf.json",
            ],
            None,
        ),
        (&["generate-channel", "straight-channel.json", "ch"], None),
        (&["generate-sweep", "sweep.json", "sw"], None),
        (
            &[
                "generate-duct",
                "elbow.stl",
                "patches.json",
                "partition.json",
                "bulk",
                "bulk",
            ],
            None,
        ),
        (
            &[
                "generate-duct",
                "elbow.stl",
                "patches.json",
                "partition.json",
                "layered",
                "lay",
            ],
            None,
        ),
        (
            &[
                "generate-duct",
                "elbow.stl",
                "patches.json",
                "partition.json",
                "polyhedral",
                "poly",
            ],
            None,
        ),
    ];
    for (args, fixed) in cases {
        let mut check = vec!["--check"];
        check.extend_from_slice(args);
        let (ok, estimate) = kofl_in(&dir, &check);
        assert!(ok, "{args:?}: {estimate}");
        assert_eq!(estimate["format"], "koflcraft.estimate");
        assert!(
            !dir.join(args.last().unwrap()).exists(),
            "--check must not build"
        );
        let (ok, result) = kofl_in(&dir, args);
        assert!(ok, "{args:?}");
        let actual = fixed.unwrap_or_else(|| built(&result));
        let predicted = estimate["cells"].as_u64().unwrap();
        let basis = estimate["cells_basis"].as_str().unwrap();
        if basis == "exact" {
            assert_eq!(predicted, actual, "{args:?}");
        } else if basis.starts_with("approximate") {
            // Voronoi core with wall columns: the column count follows the
            // built core, so the estimate bounds it from above (AG-15 judges
            // a build by it) without being far off.
            assert_eq!(args[4], "polyhedral");
            assert!(basis.contains("reentrant"));
            assert!(
                actual <= predicted && 3 * actual >= predicted,
                "{args:?}: {actual} {predicted}"
            );
        } else {
            // Polyhedral duct: cells split at reentrant edges or where they
            // break the solver rules come after the estimate (on the elbow's
            // median dual 10 and 6); the basis must say so.
            assert_eq!(args[4], "polyhedral");
            assert!(
                estimate["cells_basis"]
                    .as_str()
                    .unwrap()
                    .contains("solver rules")
            );
            assert_eq!(actual, predicted + 11, "{args:?}");
        }
    }
    let (ok, _) = kofl_in(&dir, &["--check", "check", "m.json"]);
    assert!(!ok, "only build tools have an estimate");
    let _ = std::fs::remove_dir_all(dir);
}

/// Keep rule (SW-18). Invariant (BC-02, AG-13): the boundary of a case is
/// defined by `patch.*` calls alone, each writes a new map, the mesh built
/// with the last map names its OpenFOAM patches after it, and the journal
/// of the whole session replays in a fresh directory with equal results.
/// A patch tool missing from the registry, a nondeterministic result or a
/// `--patches` flag that is not passed on breaks it; the engine tests do
/// not go through the command line or the journal.
#[test]
fn patch_tools_define_a_boundary_that_meshes_and_replays() {
    let rec = fresh("patch");
    let stl = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../test-fixtures/unit-box.stl"
    );
    let sites = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../test-fixtures/stl-sites.json"
    );
    let journal = rec.join("j.jsonl");
    let j = journal.to_str().unwrap();
    let calls = [
        (
            "patch.init",
            serde_json::json!({"surface": stl, "output": "p0.json"}),
        ),
        (
            "patch.select_plane",
            serde_json::json!({"surface": stl, "patches": "p0.json", "output": "p1.json",
                "name": "inlet", "role": "inlet", "point": [0, 0, 0], "normal": [1, 0, 0]}),
        ),
        (
            "patch.merge",
            serde_json::json!({"surface": stl, "patches": "p1.json", "output": "p2.json",
                "merge": ["south_wall", "north_wall", "floor", "ceiling"], "name": "walls", "role": "wall"}),
        ),
        (
            "patch.rename",
            serde_json::json!({"surface": stl, "patches": "p2.json", "output": "p3.json",
                "from": "east_outlet", "to": "outlet"}),
        ),
        (
            "patch.list",
            serde_json::json!({"surface": stl, "patches": "p3.json"}),
        ),
    ];
    let mut last = Value::Null;
    for (k, (id, input)) in calls.iter().enumerate() {
        let file = rec.join(format!("in{k}.json"));
        std::fs::write(&file, input.to_string()).unwrap();
        let (ok, envelope) = kofl_in(
            &rec,
            &[
                "--journal",
                j,
                "call",
                id,
                "--input",
                file.to_str().unwrap(),
            ],
        );
        assert!(ok, "{id}: {envelope:#}");
        last = envelope;
    }
    let names: Vec<&str> = last["result"]["patches"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| p["name"].as_str().unwrap())
        .collect();
    assert_eq!(names, ["outlet", "walls", "inlet"]);
    assert_eq!(last["result"]["every_facet_in_one_patch"], true);

    let (ok, built) = kofl_in(
        &rec,
        &[
            "--journal",
            j,
            "generate-stl",
            stl,
            sites,
            "m.json",
            "--patches",
            "p3.json",
        ],
    );
    assert!(ok, "{built:#}");
    let (ok, _) = kofl_in(
        &rec,
        &["--journal", j, "export", "m.json", "openfoam", "foam"],
    );
    assert!(ok);
    let boundary = std::fs::read_to_string(rec.join("foam/boundary")).unwrap();
    for (name, kind) in [("outlet", "patch"), ("walls", "wall"), ("inlet", "patch")] {
        let block = boundary.split(&format!("    {name}\n")).nth(1).unwrap();
        assert!(block.contains(&format!("type        {kind};")), "{name}");
    }

    // an unknown field of a selection tool is refused, not ignored
    let bad = rec.join("bad.json");
    std::fs::write(
        &bad,
        serde_json::json!({"surface": stl, "patches": "p3.json", "output": "p9.json",
            "name": "x", "direction": [1, 0, 0], "max_angle_deg": 1, "angle": 3})
        .to_string(),
    )
    .unwrap();
    let (ok, envelope) = kofl_in(
        &rec,
        &[
            "call",
            "patch.select_normal",
            "--input",
            bad.to_str().unwrap(),
        ],
    );
    assert!(!ok);
    assert_eq!(envelope["diagnostics"][0]["code"], "TOOL002");

    let again = fresh("patch-again");
    let (ok, envelope) = kofl_in(&again, &["run", "--script", j]);
    assert!(ok, "{envelope:#}");
    assert_eq!(envelope["result"]["replayed"], 7);
    assert_eq!(envelope["result"]["different"], 0);
    for d in [rec, again] {
        let _ = std::fs::remove_dir_all(d);
    }
}

/// Invariant (BS-02..BS-05, BS-14): the block tools chain through files
/// from the command line: a template is written and checked, the grid is
/// built into a new directory with the FoehnCraft files and evidence that
/// names every block's cell range, and a broken chain is refused with
/// BLOCK001 before any grid is built. No other test runs the block tools
/// through the registry.
#[test]
fn block_tools_build_a_pipe_from_a_template() {
    let dir = fresh("block");
    let write = |name: &str, v: Value| std::fs::write(dir.join(name), v.to_string()).unwrap();
    write(
        "t.json",
        serde_json::json!({"config": {"template": "pipe_ogrid", "radius": 0.05, "length": 0.2,
            "core_cells": 4, "radial_cells": 3, "axial_cells": 6, "wall_spacing": 0.002},
            "output": "pipe.json"}),
    );
    let (ok, e) = kofl_in(&dir, &["call", "block.template", "--input", "t.json"]);
    assert!(ok, "{e}");
    assert_eq!(e["result"]["check"]["cells"], 6 * (16 + 4 * 12));
    write(
        "g.json",
        serde_json::json!({"topology": "pipe.json", "output": "mesh"}),
    );
    let (ok, e) = kofl_in(&dir, &["call", "block.grid", "--input", "g.json"]);
    assert!(ok, "{e}");
    let r = &e["result"]["report"];
    assert_eq!(r["summary"]["inverted_faces"], 0);
    assert!(r["relative_conservation_error"].as_f64().unwrap() < 1e-10);
    let evidence: Value =
        serde_json::from_slice(&std::fs::read(dir.join("mesh/mesh-evidence.json")).unwrap())
            .unwrap();
    assert_eq!(
        evidence["block_structure"]["blocks"]
            .as_array()
            .unwrap()
            .len(),
        5
    );
    for f in ["volume.vtu", "mesh.msh", "foehncraft.vtu", "report.json"] {
        assert!(dir.join("mesh").join(f).is_file(), "{f}");
    }
    write(
        "ge.json",
        serde_json::json!({"topology": "pipe.json", "output": "smooth", "algorithm": {"kind": "elliptic"}}),
    );
    let (ok, e) = kofl_in(&dir, &["call", "block.grid", "--input", "ge.json"]);
    assert!(ok, "{e}");
    assert_eq!(e["result"]["report"]["algorithm"], "block_elliptic_v1");
    assert_eq!(e["result"]["report"]["smoothing"]["converged"], true);
    let mut topo: Value =
        serde_json::from_slice(&std::fs::read(dir.join("pipe.json")).unwrap()).unwrap();
    topo["blocks"][2]["cells"][1] = 7.into();
    write("bad.json", topo);
    write("c.json", serde_json::json!({"topology": "bad.json"}));
    let (ok, e) = kofl_in(&dir, &["call", "block.check", "--input", "c.json"]);
    assert!(!ok);
    assert_eq!(e["diagnostics"][0]["code"], "BLOCK001");
    assert!(
        e["diagnostics"][0]["message"]
            .as_str()
            .unwrap()
            .contains("count_chain")
    );
}

/// Invariant (BS-04, AG-13): the block edit tools chain through files from
/// the command line with a journal: list, set a chain count, a law and a
/// patch role, refuse a move that turns a block inside out (no file
/// written), split a block through its chain, build the edited topology,
/// and the journal replays in a fresh directory with the same results.
#[test]
fn block_edit_tools_chain_and_replay() {
    let dir = fresh("block-edit");
    let journal = dir.join("j.jsonl");
    let j = journal.to_str().unwrap();
    let call = |k: usize, id: &str, v: Value| {
        let file = dir.join(format!("in{k}.json"));
        std::fs::write(&file, v.to_string()).unwrap();
        kofl_in(
            &dir,
            &[
                "--journal",
                j,
                "call",
                id,
                "--input",
                file.to_str().unwrap(),
            ],
        )
    };
    let (ok, e) = call(
        0,
        "block.template",
        serde_json::json!({"config": {"template": "pipe_ogrid", "radius": 0.05, "length": 0.2,
            "core_cells": 4, "radial_cells": 3, "axial_cells": 6}, "output": "t0.json"}),
    );
    assert!(ok, "{e}");
    let (ok, listed) = call(1, "block.list", serde_json::json!({"topology": "t0.json"}));
    assert!(ok, "{listed}");
    let core = listed["result"]["blocks"][0]["vertices"].clone();
    let axial = [core[0].clone(), core[4].clone()];
    let (ok, e) = call(
        2,
        "block.set_count",
        serde_json::json!({"topology": "t0.json", "output": "t1.json", "edge": axial, "cells": 10}),
    );
    assert!(ok, "{e}");
    assert_eq!(e["result"]["check"]["cells"], 10 * (16 + 4 * 12));
    assert_eq!(
        e["result"]["change"]["blocks"],
        serde_json::json!([0, 1, 2, 3, 4])
    );
    let (ok, e) = call(
        3,
        "block.set_distribution",
        serde_json::json!({"topology": "t1.json", "output": "t2.json", "edge": axial,
            "distribution": {"law": "tanh", "first": 0.005, "last": 0.005}}),
    );
    assert!(ok, "{e}");
    let (ok, e) = call(
        4,
        "block.set_patch",
        serde_json::json!({"topology": "t2.json", "output": "t3.json", "name": "walls",
            "rename": "pipe_wall", "role": "wall"}),
    );
    assert!(ok, "{e}");
    // the core corner moved across the opposite face turns the block over
    let opposite =
        listed["result"]["vertices"][core[1].as_u64().unwrap() as usize]["point"].clone();
    let mut far = opposite.clone();
    far[0] = (opposite[0].as_f64().unwrap() * 3.0).into();
    let (ok, e) = call(
        5,
        "block.move_vertex",
        serde_json::json!({"topology": "t3.json", "output": "bad.json", "vertex": core[0], "to": far}),
    );
    assert!(!ok, "{e}");
    assert!(
        e["diagnostics"][0]["message"]
            .as_str()
            .unwrap()
            .contains("right-handed"),
        "{e}"
    );
    assert!(!dir.join("bad.json").exists());
    let (ok, e) = call(
        6,
        "block.split_block",
        serde_json::json!({"topology": "t3.json", "output": "t4.json", "block": 0, "axis": 2, "at": 4}),
    );
    assert!(ok, "{e}");
    assert_eq!(e["result"]["check"]["cells"], 10 * (16 + 4 * 12));
    assert_eq!(e["result"]["check"]["blocks"], 10);
    let (ok, e) = call(
        7,
        "block.grid",
        serde_json::json!({"topology": "t4.json", "output": "mesh"}),
    );
    assert!(ok, "{e}");
    assert_eq!(e["result"]["report"]["summary"]["inverted_faces"], 0);
    let edited: Value =
        serde_json::from_slice(&std::fs::read(dir.join("t3.json")).unwrap()).unwrap();
    let wall = edited["patches"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["name"] == "pipe_wall")
        .unwrap();
    assert_eq!(wall["role"], "wall");

    let again = fresh("block-edit-replay");
    // the refused move replays as the same refusal
    let (_, envelope) = kofl_in(&again, &["run", "--script", j, "--keep-going"]);
    assert_eq!(
        envelope["result"]["format"], "koflcraft.replay",
        "{envelope:#}"
    );
    assert_eq!(envelope["result"]["replayed"], 8, "{envelope:#}");
    assert_eq!(envelope["result"]["different"], 0, "{envelope:#}");
}

/// Invariant (BS-03, BS-07): the elbow and the airfoil C-grid run through
/// the registry: the elbow file carries its projection surfaces, and an
/// airfoil read from a Selig coordinate file in the working directory
/// builds a mesh with the airfoil as wall.
#[test]
fn block_templates_build_an_elbow_and_an_airfoil_from_a_file() {
    let dir = fresh("block-templates");
    let write = |name: &str, v: Value| std::fs::write(dir.join(name), v.to_string()).unwrap();
    write(
        "e.json",
        serde_json::json!({"config": {"template": "elbow", "radius": 0.05, "bend_radius": 0.15,
            "angle": 90, "inlet_length": 0.1, "core_cells": 4, "radial_cells": 3,
            "axial_cells": 12}, "output": "elbow.json"}),
    );
    let (ok, e) = kofl_in(&dir, &["call", "block.template", "--input", "e.json"]);
    assert!(ok, "{e}");
    assert_eq!(e["result"]["check"]["blocks"], 20);
    let elbow: Value =
        serde_json::from_slice(&std::fs::read(dir.join("elbow.json")).unwrap()).unwrap();
    assert_eq!(elbow["surfaces"][1]["shape"]["kind"], "torus");
    // a NACA 0012 in Selig order, trailing edge over the upper side
    let n = 60;
    let y = |x: f64| {
        0.6 * (0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x * x + 0.2843 * x.powi(3)
            - 0.1036 * x.powi(4))
    };
    let mut dat = String::from("NACA 0012\n");
    for i in (0..=n).rev() {
        let x = 0.5 * (1.0 - (std::f64::consts::PI * i as f64 / n as f64).cos());
        dat += &format!("{x:.8} {:.8}\n", y(x));
    }
    for i in 1..=n {
        let x = 0.5 * (1.0 - (std::f64::consts::PI * i as f64 / n as f64).cos());
        dat += &format!("{x:.8} {:.8}\n", -y(x));
    }
    std::fs::write(dir.join("n0012.dat"), dat).unwrap();
    write(
        "a.json",
        serde_json::json!({"config": {"template": "airfoil_cgrid", "coordinates": "n0012.dat",
            "span": 0.1, "farfield": 5, "wake_length": 5, "surface_cells": 24,
            "normal_cells": 12, "wake_cells": 12, "wall_spacing": 0.002}, "output": "af.json"}),
    );
    let (ok, e) = kofl_in(&dir, &["call", "block.template", "--input", "a.json"]);
    assert!(ok, "{e}");
    assert_eq!(e["result"]["check"]["blocks"], 6);
    write(
        "g.json",
        serde_json::json!({"topology": "af.json", "output": "af"}),
    );
    let (ok, e) = kofl_in(&dir, &["call", "block.grid", "--input", "g.json"]);
    assert!(ok, "{e}");
    assert_eq!(e["result"]["report"]["summary"]["inverted_faces"], 0);
}
