//! Where a tool call may read and write (AG-14) and the limits of the front
//! end that made it (AG-15).

use std::ffi::OsString;
use std::path::{Path, PathBuf};

use mesh_engine::MeshDirectory;

use crate::diag::Diagnostic;

/// Limits a front end puts on tool inputs. The command line has none; the
/// MCP server uses [`Limits::mcp`] (the limits `mesh-mcp` had), which a
/// policy file can change field by field (AG-15).
#[derive(Clone, Debug, PartialEq, serde::Deserialize, serde::Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct Limits {
    pub max_sites: usize,
    pub max_surface_facets: usize,
    pub max_file_bytes: usize,
    pub max_channel_layers: usize,
    pub max_channel_axial_cells: usize,
    pub max_sweep_cells: usize,
    pub max_duct_facets: usize,
    pub max_duct_layers: usize,
    pub max_duct_regions: usize,
    /// Rows per page of per-cell and per-face detail (AG-16).
    pub max_page_rows: usize,
    /// Cells one build may make without approval (AG-15).
    pub max_build_cells: usize,
    /// Estimated output bytes one build may write without approval.
    pub max_output_bytes: u64,
}

impl Default for Limits {
    fn default() -> Self {
        Self::mcp()
    }
}

impl Limits {
    pub fn mcp() -> Self {
        Self {
            max_sites: 256,
            max_surface_facets: 2048,
            max_file_bytes: 8 * 1024 * 1024,
            max_channel_layers: 8,
            max_channel_axial_cells: 32,
            max_sweep_cells: 500_000,
            max_duct_facets: 100_000,
            max_duct_layers: 4,
            max_duct_regions: 64,
            max_page_rows: 500,
            max_build_cells: 500_000,
            max_output_bytes: 1 << 30,
        }
    }
}

/// The context of one call.
#[derive(Clone, Debug)]
pub struct Context {
    /// Relative paths are relative to this directory.
    pub base: PathBuf,
    /// Files may be read only below one of these (canonical); `None`: anywhere.
    pub read_roots: Option<Vec<PathBuf>>,
    /// Outputs may be created only below one of these; `None`: anywhere.
    pub write_roots: Option<Vec<PathBuf>>,
    pub limits: Option<Limits>,
    /// The front end: `cli`, `mcp`, `replay`.
    pub via: &'static str,
}

fn invalid(path: &Path, message: &str) -> Diagnostic {
    Diagnostic::new("invalid_path", message).with_path(path.display())
}

fn io(path: &Path, error: impl std::fmt::Display) -> Diagnostic {
    Diagnostic::new("io_error", error).with_path(path.display())
}

fn canonical_roots(roots: &[PathBuf]) -> Result<Vec<PathBuf>, Diagnostic> {
    roots
        .iter()
        .map(|r| r.canonicalize().map_err(|e| io(r, e)))
        .collect()
}

impl Context {
    /// The command line: the current directory, no roots, no limits.
    pub fn cli() -> Self {
        Self {
            base: std::env::current_dir().unwrap_or_else(|_| PathBuf::from(".")),
            read_roots: None,
            write_roots: None,
            limits: None,
            via: "cli",
        }
    }

    /// A confined context: paths relative to the first write root, reads
    /// and writes inside the roots, with limits.
    pub fn confined(
        read_roots: Vec<PathBuf>,
        write_roots: Vec<PathBuf>,
        limits: Option<Limits>,
        via: &'static str,
    ) -> Result<Self, Diagnostic> {
        let read_roots = canonical_roots(&read_roots)?;
        let write_roots = canonical_roots(&write_roots)?;
        let base = write_roots
            .first()
            .or(read_roots.first())
            .cloned()
            .ok_or_else(|| Diagnostic::new("usage", "a confined context needs a root"))?;
        Ok(Self {
            base,
            read_roots: Some(read_roots),
            write_roots: Some(write_roots),
            limits,
            via,
        })
    }

    pub fn resolve(&self, path: &str) -> PathBuf {
        let p = Path::new(path);
        if p.is_absolute() {
            p.to_path_buf()
        } else {
            self.base.join(p)
        }
    }

    /// Read a regular file, inside the read roots and the size limit.
    pub fn read(&self, path: &str) -> Result<Vec<u8>, Diagnostic> {
        let full = self.resolve(path);
        let canonical = full.canonicalize().map_err(|e| io(&full, e))?;
        if let Some(roots) = &self.read_roots
            && !roots.iter().any(|r| canonical.starts_with(r))
        {
            return Err(invalid(&full, "the file is outside the read roots"));
        }
        let (Some(parent), Some(name)) = (canonical.parent(), canonical.file_name()) else {
            return Err(invalid(&full, "expected a file"));
        };
        let limit = self
            .limits
            .as_ref()
            .map_or(usize::MAX - 1, |l| l.max_file_bytes);
        MeshDirectory::open(parent)
            .and_then(|dir| dir.read_regular(name, limit))
            .map_err(|e| io(&full, e))
    }

    /// The directory and new entry name of an output, inside the write roots.
    pub fn target(&self, path: &str) -> Result<(MeshDirectory, OsString), Diagnostic> {
        let full = self.resolve(path);
        let Some(name) = full.file_name().map(|n| n.to_os_string()) else {
            return Err(invalid(&full, "expected an output name"));
        };
        if name == ".." || name == "." {
            return Err(invalid(&full, "expected an output name"));
        }
        let parent = full
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let canonical = parent.canonicalize().map_err(|e| io(parent, e))?;
        if let Some(roots) = &self.write_roots
            && !roots.iter().any(|r| canonical.starts_with(r))
        {
            return Err(invalid(&full, "the output is outside the write roots"));
        }
        let dir = MeshDirectory::open(&canonical).map_err(|e| io(parent, e))?;
        Ok((dir, name))
    }

    pub(crate) fn limit(
        &self,
        what: &str,
        value: usize,
        max: impl Fn(&Limits) -> usize,
    ) -> Result<(), Diagnostic> {
        match &self.limits {
            Some(l) if value > max(l) => Err(Diagnostic::new(
                "limit_exceeded",
                format!(
                    "{} accepts at most {} {what}, got {value}",
                    self.via,
                    max(l)
                ),
            )),
            _ => Ok(()),
        }
    }
}
