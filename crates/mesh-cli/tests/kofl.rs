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
