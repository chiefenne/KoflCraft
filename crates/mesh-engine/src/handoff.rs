//! Files for FoehnCraft (requirements section 16): the single-file VTU of
//! FC-02 with the digests of FC-05 in its field data, and the side-car
//! evidence `koflcraft.mesh-evidence` that names the files of a bundle with
//! their SHA-256, the summary FoehnCraft must reproduce (FC-06) and the
//! patch roles (FC-04).

use std::io::{self, Write};

use mesh_io::StreamError;
use mesh_io::handoff::{HandoffMesh, HandoffSummary, Patch, write_foehncraft_vtu};
use serde::Serialize;

use crate::{EngineError, ErrorCode};

pub const EVIDENCE_FORMAT: &str = "koflcraft.mesh-evidence";
pub const EVIDENCE_FORMAT_VERSION: &str = "1.0";

/// The name of the FoehnCraft VTU in a bundle.
pub const FOEHNCRAFT_VTU: &str = "foehncraft.vtu";
/// The name of the evidence side-car in a bundle.
pub const MESH_EVIDENCE: &str = "mesh-evidence.json";

#[derive(Serialize)]
pub struct MeshEvidence {
    pub format: &'static str,
    pub format_version: &'static str,
    pub koflcraft_version: &'static str,
    pub generator: String,
    pub summary: HandoffSummary,
    pub patches: Vec<Patch>,
    /// Every other file of the bundle with its SHA-256.
    pub files: Vec<(String, String)>,
    pub limitations: Vec<String>,
}

pub(crate) fn export_error(error: impl std::fmt::Display) -> EngineError {
    EngineError::new(ErrorCode::ExportFailed, error)
}

/// A file that could not be streamed: the file system's errors as
/// `io_error`, the writer's as `export_failed`.
pub(crate) fn stream_error(error: StreamError) -> EngineError {
    match error {
        StreamError::Io(e) => EngineError::new(ErrorCode::IoError, e),
        StreamError::Content(e) => export_error(e),
    }
}

/// Write the FoehnCraft VTU of a mesh, with `provenance` in its field data.
pub(crate) fn write_foehncraft(
    mesh: &HandoffMesh,
    provenance: &[(&str, String)],
    out: &mut dyn Write,
) -> io::Result<()> {
    let mut metadata = vec![("koflcraft_version", env!("CARGO_PKG_VERSION").to_string())];
    metadata.extend(provenance.iter().cloned());
    write_foehncraft_vtu(mesh, &metadata, out)
}

/// The evidence of a bundle whose files (other than the evidence itself)
/// are `files`, given with their SHA-256 as the bundle wrote them.
pub(crate) fn evidence(
    mesh: &HandoffMesh,
    generator: &str,
    files: &[(String, String)],
) -> Result<Vec<u8>, EngineError> {
    let mut limitations = Vec::new();
    let empty = mesh.empty_patches();
    if !empty.is_empty() {
        limitations.push(format!(
            "patches without faces are not in the files: {}",
            empty.join(", ")
        ));
    }
    let unspecified: Vec<_> = mesh
        .patches
        .iter()
        .filter(|p| p.role.foehncraft_type().is_none())
        .map(|p| p.name.as_str())
        .collect();
    if !unspecified.is_empty() {
        limitations.push(format!(
            "no boundary-condition role could be derived for: {}",
            unspecified.join(", ")
        ));
    }
    let inverted = mesh.inverted_faces();
    if inverted > 0 {
        limitations.push(format!(
            "{inverted} interior faces have S·d ≤ 0; FoehnCraft reports such a mesh as invalid"
        ));
    }
    limitations.push("flow accuracy on this mesh is untested".to_string());
    let evidence = MeshEvidence {
        format: EVIDENCE_FORMAT,
        format_version: EVIDENCE_FORMAT_VERSION,
        koflcraft_version: env!("CARGO_PKG_VERSION"),
        generator: generator.to_string(),
        summary: mesh.summary(),
        patches: mesh.patches.clone(),
        files: files.to_vec(),
        limitations,
    };
    serde_json::to_vec_pretty(&evidence)
        .map_err(|e| EngineError::new(ErrorCode::SerializationFailed, e))
}
