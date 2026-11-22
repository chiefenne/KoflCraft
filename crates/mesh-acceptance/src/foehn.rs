//! FoehnCraft round trip (FC-06): every mesh exported for FoehnCraft is read
//! by `foehn mesh`, run as a separate process like Gmsh, and must report the
//! counts, volume and patch areas of KoflCraft's own summary
//! (`mesh-evidence.json`), every cell closed, and the same verdict on
//! inverted faces (S·d ≤ 0).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::{Value, json};

use crate::{RunResult, require};

/// Relative tolerance on volume and patch areas (FC-06).
pub const TOLERANCE: f64 = 1e-12;
/// FoehnCraft's own closure limit for a valid mesh.
const CLOSURE: f64 = 1e-10;

/// The FoehnCraft binary: `$FOEHN`, else `foehn` on PATH, with its version.
pub fn find() -> RunResult<Option<(PathBuf, String)>> {
    let path = std::env::var_os("FOEHN")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("foehn"));
    match Command::new(&path).arg("--version").output() {
        Ok(out) if out.status.success() => Ok(Some((
            path,
            String::from_utf8_lossy(&out.stdout).trim().to_string(),
        ))),
        Ok(out) => Err(format!(
            "{} --version failed: {}",
            path.display(),
            String::from_utf8_lossy(&out.stderr)
        )
        .into()),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(e) => Err(e.into()),
    }
}

/// How the file relates to KoflCraft's face-based mesh.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// `foehncraft.vtu` or `polyMesh`: the same cells and faces.
    Faces,
    /// Gmsh MSH: a tetrahedral subdivision; only volume, areas and names.
    Subdivision,
}

fn rel(a: f64, b: f64) -> f64 {
    if b == 0.0 {
        a.abs()
    } else {
        (a / b - 1.0).abs()
    }
}

/// Read `file` with `foehn mesh` and compare it with `expected` (the
/// `summary` of `mesh-evidence.json`).
pub fn check(
    foehn: &Path,
    file: &Path,
    expected: &Value,
    kind: Kind,
    out: &Path,
) -> RunResult<Value> {
    let output = Command::new(foehn)
        .arg("mesh")
        .arg(file)
        .args(["--json", "--quiet", "--out"])
        .arg(out)
        .output()?;
    let code = output.status.code().unwrap_or(-1);
    require(
        code == 0 || code == 2,
        &format!(
            "foehn mesh {}: exit {code}: {}",
            file.display(),
            String::from_utf8_lossy(&output.stderr)
        ),
    )?;
    let report: Value = serde_json::from_slice(&fs::read(out.join("mesh_report.json"))?)?;
    let q = &report["quality"];
    let num = |v: &Value| v.as_f64().ok_or("missing number in the FoehnCraft report");
    let count = |v: &Value| v.as_u64().ok_or("missing count in the FoehnCraft report");
    let name = file.display().to_string();
    let fail = |what: String| -> RunResult<()> { Err(format!("{name}: {what}").into()) };

    let mut worst: f64 = rel(num(&q["total_volume"])?, num(&expected["total_volume"])?);
    let boundaries: Vec<(String, u64, f64)> = q["boundaries"]
        .as_array()
        .ok_or("missing boundaries in the FoehnCraft report")?
        .iter()
        .map(|b| {
            Ok((
                b[0].as_str().unwrap_or_default().to_string(),
                count(&b[1])?,
                num(&b[2])?,
            ))
        })
        .collect::<RunResult<_>>()?;
    let patches = expected["patches"].as_array().ok_or("missing patches")?;
    let used: Vec<&Value> = patches.iter().filter(|p| p["faces"] != 0).collect();
    if boundaries.len() != used.len() {
        fail(format!(
            "{} boundaries, KoflCraft wrote {}",
            boundaries.len(),
            used.len()
        ))?;
    }
    for p in used {
        let pname = p["name"].as_str().unwrap_or_default();
        let Some((_, faces, area)) = boundaries.iter().find(|b| b.0 == pname) else {
            return Err(format!("{name}: FoehnCraft has no boundary {pname}").into());
        };
        if kind == Kind::Faces && *faces != count(&p["faces"])? {
            fail(format!(
                "boundary {pname} has {faces} faces, KoflCraft {}",
                p["faces"]
            ))?;
        }
        worst = worst.max(rel(*area, num(&p["area"])?));
    }
    if worst > TOLERANCE {
        fail(format!("volume or patch area differs by {worst:e}"))?;
    }
    let closure = num(&q["closure"]["max"])?;
    if closure > CLOSURE {
        fail(format!("a cell is not closed (|Σ n dS| = {closure:e})"))?;
    }
    let inverted = count(&q["inverted_faces"])?;
    let valid = code == 0;
    if kind == Kind::Faces {
        for key in ["cells", "faces", "interior_faces", "boundary_faces"] {
            if q[key] != expected[key] {
                fail(format!(
                    "{key}: FoehnCraft {}, KoflCraft {}",
                    q[key], expected[key]
                ))?;
            }
        }
        if inverted != count(&expected["inverted_faces"])? {
            fail(format!(
                "{inverted} inverted faces, KoflCraft predicted {}",
                expected["inverted_faces"]
            ))?;
        }
        if valid != (inverted == 0) {
            fail(format!(
                "FoehnCraft valid = {valid} with {inverted} inverted faces"
            ))?;
        }
    } else if !valid {
        fail("FoehnCraft calls the tetrahedral subdivision invalid".to_string())?;
    }
    Ok(json!({
        "file": name,
        "source": report["source"],
        "foehncraft_valid": valid,
        "cells": q["cells"],
        "inverted_faces": inverted,
        "max_relative_difference": worst,
        "max_closure": closure,
    }))
}
