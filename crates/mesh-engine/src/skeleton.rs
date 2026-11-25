//! The FoehnCraft boundary-condition skeleton (BC-05): a case file that
//! names every patch of a map with the FoehnCraft type its role maps to and
//! leaves every physical value to FoehnCraft. KoflCraft never writes a
//! value (no velocity, no pressure, no fluid): the lines that need one are
//! comments marked `TODO`, and a patch without a role gets an empty table,
//! so `foehn case check` refuses the file until a person or an agent has
//! decided. This keeps one source of truth for physics (FoehnCraft) and one
//! for geometry and names (KoflCraft).

use std::fmt::Write;

use serde::Serialize;

use crate::patches::{PatchMap, PatchRole, PatchSurface};
use crate::{EngineError, ErrorCode};

pub const SKELETON_FORMAT: &str = "koflcraft.bc-skeleton";
/// The FoehnCraft case format the skeleton is written in.
pub const FOEHNCRAFT_CASE_FORMAT: &str = "0.7";

/// What the skeleton says about one patch.
#[derive(Clone, Debug, Serialize)]
pub struct SkeletonPatch {
    pub name: String,
    pub role: PatchRole,
    /// The FoehnCraft type written, `null` when the role gives none.
    pub foehncraft_type: Option<&'static str>,
    /// The values FoehnCraft still needs (comments in the file).
    pub todo: Vec<&'static str>,
}

/// The skeleton text and what it leaves open.
#[derive(Clone, Debug, Serialize)]
pub struct Skeleton {
    pub format: &'static str,
    pub case_format: &'static str,
    pub mesh: String,
    pub patches: Vec<SkeletonPatch>,
    /// Patches without a FoehnCraft type: set their roles first.
    pub needs_role: Vec<String>,
    #[serde(skip)]
    pub text: String,
}

/// TOML basic string with escapes.
fn quoted(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if (c as u32) < 0x20 || c == '\u{7f}' => {
                let _ = write!(out, "\\u{:04X}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn todo(role: PatchRole) -> &'static [(&'static str, &'static str)] {
    match role {
        PatchRole::Inlet => &[(
            "normal_velocity",
            "m/s into the domain; or `velocity`, `profile`, `table` (foehn case schema)",
        )],
        PatchRole::Outlet => &[("pressure", "Pa; FoehnCraft takes 0 when it is left out")],
        _ => &[],
    }
}

fn fmt3(v: [f64; 3]) -> String {
    format!("[{:.6e}, {:.6e}, {:.6e}]", v[0], v[1], v[2])
}

impl PatchSurface {
    /// The skeleton of the boundary conditions of `map`, for a mesh file
    /// `mesh` (the path written into `[mesh] file`).
    pub fn bc_skeleton(&self, map: &PatchMap, mesh: &str) -> Result<Skeleton, EngineError> {
        if mesh.is_empty() {
            return Err(EngineError::new(
                ErrorCode::InvalidArguments,
                "mesh must name the mesh file the case reads",
            ));
        }
        let summary = self.summary(map)?;
        let mut t = String::new();
        let _ = writeln!(
            t,
            "# FoehnCraft boundary-condition skeleton, written by KoflCraft {} (BC-05).",
            env!("CARGO_PKG_VERSION")
        );
        let _ = writeln!(t, "#");
        let _ = writeln!(
            t,
            "# KoflCraft names the patches and their types; it sets no physical value."
        );
        let _ = writeln!(
            t,
            "# Fill in every TODO with FoehnCraft (`foehn case check` lists what is missing)."
        );
        let _ = writeln!(t, "# Surface sha256: {}", summary.source_sha256);
        let _ = writeln!(t);
        let _ = writeln!(t, "case_format = {}", quoted(FOEHNCRAFT_CASE_FORMAT));
        let _ = writeln!(t);
        let _ = writeln!(t, "[mesh]");
        let _ = writeln!(t, "file = {}", quoted(mesh));
        let _ = writeln!(t);
        let _ = writeln!(t, "# TODO [fluid]: density (kg/m³) and viscosity (Pa s).");
        let mut patches = Vec::with_capacity(summary.patches.len());
        let mut needs_role = Vec::new();
        for p in &summary.patches {
            let _ = writeln!(t);
            let mut evidence = format!("# {}: {} facets, area {:.6e}", p.name, p.facets, p.area);
            if let Some(n) = p.normal
                && p.max_normal_deviation_deg < 1.0
            {
                let _ = write!(evidence, ", plane with normal {}", fmt3(n));
            }
            if let Some([lo, hi]) = p.bounds {
                let mid = [0, 1, 2].map(|k| 0.5 * (lo[k] + hi[k]));
                let _ = write!(evidence, ", centre of bounds {}", fmt3(mid));
            }
            let _ = writeln!(t, "{evidence}");
            // patch names are bare TOML keys (letters, digits, `_`, `-`; FC-04)
            let _ = writeln!(t, "[boundary.{}]", p.name);
            let kind = p.role.foehncraft_type();
            let mut open = Vec::new();
            match kind {
                Some(k) => {
                    let _ = writeln!(t, "type = {}", quoted(k));
                    for (key, note) in todo(p.role) {
                        let _ = writeln!(t, "# TODO {key} = …  ({note})");
                        open.push(*key);
                    }
                }
                None => {
                    let _ = writeln!(
                        t,
                        "# TODO type: the role is {}; give the patch a role (wall, inlet, outlet, \
                         symmetry) with patch.set_role and write the skeleton again.",
                        format!("{:?}", p.role).to_lowercase()
                    );
                    open.push("type");
                    needs_role.push(p.name.clone());
                }
            }
            patches.push(SkeletonPatch {
                name: p.name.clone(),
                role: p.role,
                foehncraft_type: kind,
                todo: open,
            });
        }
        Ok(Skeleton {
            format: SKELETON_FORMAT,
            case_format: FOEHNCRAFT_CASE_FORMAT,
            mesh: mesh.to_string(),
            patches,
            needs_role,
            text: t,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::quoted;

    #[test]
    fn quoting_escapes() {
        assert_eq!(quoted("a\"b\\c\n"), "\"a\\\"b\\\\c\\u000A\"");
    }
}
