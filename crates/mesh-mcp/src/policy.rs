//! The approval policy of the MCP server (AG-15), from the TOML file given
//! with `--policy`. Each action kind is `allow`, `ask` or `refuse`
//! (`default` takes the built-in rule):
//!
//! ```toml
//! [policy]
//! start_build = "allow"     # a build within the estimate limits
//! exceed_limits = "ask"     # a build over them, or over any input limit
//! write_outside_roots = "refuse"   # fixed: AG-14
//! overwrite = "refuse"             # fixed: outputs are always new
//!
//! [limits]                  # any field of the MCP limits, for example
//! max_build_cells = 500000
//! max_output_bytes = 1073741824
//! ```
//!
//! The defaults are the maintainer's decision of 2026-11-22: builds within
//! the estimate limits run, larger ones ask.

use std::path::{Path, PathBuf};

use mesh_tools::Limits;
use serde::Deserialize;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Rule {
    Default,
    Allow,
    Ask,
    Refuse,
}

#[derive(Debug, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct Rules {
    start_build: Rule,
    exceed_limits: Rule,
    write_outside_roots: Rule,
    overwrite: Rule,
}

impl Default for Rules {
    fn default() -> Self {
        Self {
            start_build: Rule::Default,
            exceed_limits: Rule::Default,
            write_outside_roots: Rule::Default,
            overwrite: Rule::Default,
        }
    }
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct File {
    policy: Rules,
    limits: Limits,
}

/// The effective policy.
#[derive(Clone, Debug, PartialEq)]
pub struct Policy {
    pub start_build: Rule,
    pub exceed_limits: Rule,
    pub limits: Limits,
    /// The file it came from, if any.
    pub source: Option<PathBuf>,
}

impl Default for Policy {
    fn default() -> Self {
        Self {
            start_build: Rule::Allow,
            exceed_limits: Rule::Ask,
            limits: Limits::mcp(),
            source: None,
        }
    }
}

fn resolve(rule: Rule, default: Rule) -> Rule {
    if rule == Rule::Default { default } else { rule }
}

impl Policy {
    /// Read a policy file. It must not lie below a write root, where a tool
    /// call could have written it.
    pub fn load(path: &Path, write_roots: &[PathBuf]) -> Result<Self, String> {
        let canonical = path
            .canonicalize()
            .map_err(|e| format!("policy file {}: {e}", path.display()))?;
        for root in write_roots {
            if root.canonicalize().is_ok_and(|r| canonical.starts_with(r)) {
                return Err(format!(
                    "the policy file {} is inside the write root {}; keep it where tools cannot write",
                    path.display(),
                    root.display()
                ));
            }
        }
        let text = std::fs::read_to_string(&canonical)
            .map_err(|e| format!("policy file {}: {e}", path.display()))?;
        Self::parse(&text).map(|p| Self {
            source: Some(canonical),
            ..p
        })
    }

    pub fn parse(text: &str) -> Result<Self, String> {
        let file: File = toml::from_str(text).map_err(|e| format!("policy file: {e}"))?;
        for (name, rule) in [
            ("write_outside_roots", file.policy.write_outside_roots),
            ("overwrite", file.policy.overwrite),
        ] {
            if !matches!(rule, Rule::Default | Rule::Refuse) {
                return Err(format!(
                    "policy file: {name} can only be `refuse` (reads and writes stay in the roots, outputs are always new)"
                ));
            }
        }
        let d = Self::default();
        Ok(Self {
            start_build: resolve(file.policy.start_build, d.start_build),
            exceed_limits: resolve(file.policy.exceed_limits, d.exceed_limits),
            limits: file.limits,
            source: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Invariant (AG-15): an empty file gives the decided defaults (allow
    /// within the limits, ask above), limits override field by field, and
    /// the fixed kinds cannot be relaxed. Changing a default or accepting
    /// `allow` for writes outside the roots breaks it.
    #[test]
    fn defaults_overrides_and_fixed_kinds() {
        let p = Policy::parse("").unwrap();
        assert_eq!((p.start_build, p.exceed_limits), (Rule::Allow, Rule::Ask));
        assert_eq!(p.limits, Limits::mcp());
        let p = Policy::parse(
            "[policy]\nstart_build = \"ask\"\nexceed_limits = \"default\"\n[limits]\nmax_build_cells = 10\n",
        )
        .unwrap();
        assert_eq!((p.start_build, p.exceed_limits), (Rule::Ask, Rule::Ask));
        assert_eq!(p.limits.max_build_cells, 10);
        assert_eq!(p.limits.max_sites, Limits::mcp().max_sites);
        assert!(Policy::parse("[policy]\nwrite_outside_roots = \"allow\"\n").is_err());
        assert!(Policy::parse("[policy]\noverwrite = \"ask\"\n").is_err());
        assert!(Policy::parse("[limits]\nmax_cellz = 1\n").is_err());
        assert!(Policy::parse("[policy]\nstart_build = \"maybe\"\n").is_err());
    }
}
