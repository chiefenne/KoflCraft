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
