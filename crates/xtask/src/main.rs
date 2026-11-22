//! Repository tasks (SW-19): `cargo xtask layering`.
//!
//! `layering` reads the `[dependencies]` of every crate under `crates/`
//! and checks the layers of section 6 of `docs/requirements.md`: a crate
//! depends only on crates of a lower layer (with the exceptions listed in
//! `ALLOWED`), and the front ends reach the engine only through
//! `mesh-tools`, in their manifests and in their source.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

/// The layer of each workspace crate, lowest first.
const LAYERS: &[(&str, u8)] = &[
    ("mesh-core", 0),
    ("mesh-geometry", 1),
    ("mesh-quality", 1),
    ("mesh-algorithms", 2),
    ("mesh-io", 2),
    ("mesh-engine", 3),
    ("mesh-tools", 4),
    ("mesh-cli", 5),
    ("mesh-mcp", 5),
    // repository tooling: may use any crate, nothing may use it
    ("mesh-acceptance", 9),
    ("xtask", 9),
];

/// Same-layer dependencies that are allowed, with the reason.
const ALLOWED: &[(&str, &str, &str)] = &[
    (
        "mesh-quality",
        "mesh-geometry",
        "validation checks meshes against the convex surface they fill",
    ),
    (
        "mesh-io",
        "mesh-algorithms",
        "the writers read the sweep, channel and duct meshes the algorithms build",
    ),
    (
        "mesh-cli",
        "mesh-mcp",
        "`kofl mcp` starts the MCP server, which itself goes through mesh-tools",
    ),
];

/// Front ends: they may not name the engine or anything below it.
const FRONT_ENDS: &[&str] = &["mesh-cli", "mesh-mcp"];
const BELOW_TOOLS: &[&str] = &[
    "mesh-core",
    "mesh-geometry",
    "mesh-quality",
    "mesh-algorithms",
    "mesh-io",
    "mesh-engine",
];

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>()
        .as_slice()
    {
        ["layering"] => {
            let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
            match layering(&root) {
                Ok(report) => {
                    println!("{report}");
                    ExitCode::SUCCESS
                }
                Err(errors) => {
                    for e in errors {
                        eprintln!("layering: {e}");
                    }
                    ExitCode::FAILURE
                }
            }
        }
        _ => {
            eprintln!("Usage: cargo xtask layering");
            ExitCode::FAILURE
        }
    }
}

/// The names in the `[dependencies]` table of a manifest (not dev or build
/// dependencies, which do not ship).
fn dependencies(manifest: &str) -> Vec<String> {
    let mut deps = Vec::new();
    let mut inside = false;
    for line in manifest.lines() {
        let line = line.trim();
        if line.starts_with('[') {
            inside = line == "[dependencies]";
            continue;
        }
        if !inside || line.is_empty() || line.starts_with('#') {
            continue;
        }
        let key = line
            .split(['=', '.'])
            .next()
            .unwrap_or("")
            .trim()
            .trim_matches('"');
        if !key.is_empty() {
            deps.push(key.to_string());
        }
    }
    deps
}

fn rust_files(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            rust_files(&path, out);
        } else if path.extension().is_some_and(|e| e == "rs") {
            out.push(path);
        }
    }
}

fn layering(root: &Path) -> Result<String, Vec<String>> {
    let layer: BTreeMap<&str, u8> = LAYERS.iter().copied().collect();
    let mut errors = Vec::new();
    let mut edges = 0;
    let mut crates: Vec<PathBuf> = std::fs::read_dir(root.join("crates"))
        .map_err(|e| vec![format!("cannot read crates/: {e}")])?
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.join("Cargo.toml").is_file())
        .collect();
    crates.sort();
    for dir in &crates {
        let name = dir
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        let Some(&own) = layer.get(name.as_str()) else {
            errors.push(format!(
                "crate {name} has no layer: add it to LAYERS in crates/xtask and to section 6 of docs/requirements.md"
            ));
            continue;
        };
        let manifest = std::fs::read_to_string(dir.join("Cargo.toml")).unwrap_or_default();
        for dep in dependencies(&manifest) {
            let Some(&theirs) = layer.get(dep.as_str()) else {
                continue; // an external crate
            };
            edges += 1;
            let allowed = ALLOWED.iter().any(|(a, b, _)| *a == name && *b == dep);
            if theirs >= own && !allowed {
                errors.push(format!(
                    "{name} (L{own}) depends on {dep} (L{theirs}); a crate may depend only on lower layers"
                ));
            }
            if FRONT_ENDS.contains(&name.as_str()) && BELOW_TOOLS.contains(&dep.as_str()) {
                errors.push(format!(
                    "front end {name} depends on {dep}; front ends reach the engine only through mesh-tools"
                ));
            }
        }
        if FRONT_ENDS.contains(&name.as_str()) {
            let mut files = Vec::new();
            rust_files(&dir.join("src"), &mut files);
            files.sort();
            for file in files {
                let text = std::fs::read_to_string(&file).unwrap_or_default();
                for below in BELOW_TOOLS {
                    let path = format!("{}::", below.replace('-', "_"));
                    if text.contains(&path) {
                        errors.push(format!(
                            "{} names {path}; front ends reach the engine only through mesh-tools",
                            file.strip_prefix(root).unwrap_or(&file).display()
                        ));
                    }
                }
            }
        }
    }
    if errors.is_empty() {
        Ok(format!(
            "layering: {} crates, {edges} workspace dependencies, all downward; front ends use mesh-tools only",
            crates.len()
        ))
    } else {
        Err(errors)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Invariant: only `[dependencies]` count, so a dev-dependency on the
    /// engine (tests compare against it) is not a layering error, and both
    /// `name.workspace = true` and `name = {…}` forms are read. A parser
    /// that read every table, or only one form, breaks it.
    #[test]
    fn reads_only_the_dependencies_table() {
        let manifest = "[package]\nname = \"x\"\n\n[dependencies]\nmesh-tools.workspace = true\nserde = { version = \"1\" }\n\n[dev-dependencies]\nmesh-engine.workspace = true\n";
        assert_eq!(dependencies(manifest), ["mesh-tools", "serde"]);
    }

    /// Invariant: the workspace as committed passes, so CI catches the
    /// first upward or front-end dependency.
    #[test]
    fn the_workspace_is_layered() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        if let Err(e) = layering(&root) {
            panic!("{e:#?}");
        }
    }
}
