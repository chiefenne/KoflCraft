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
        if estimate["cells_basis"] == "exact" {
            assert_eq!(predicted, actual, "{args:?}");
        } else {
            // Polyhedral duct: the elbow splits 10 dual cells at its
            // reentrant edge after the estimate (research note
            // m4-reentrant-dual-split); the basis must say so.
            assert_eq!(args[4], "polyhedral");
            assert!(
                estimate["cells_basis"]
                    .as_str()
                    .unwrap()
                    .contains("reentrant")
            );
            assert_eq!(actual, predicted + 10, "{args:?}");
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
