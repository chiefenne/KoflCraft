//! Files for FoehnCraft (requirements section 16): the single-file VTU of
//! FC-02 with the digests of FC-05 in its field data, and the side-car
//! evidence `koflcraft.mesh-evidence` that names the files of a bundle with
//! their SHA-256, the summary FoehnCraft must reproduce (FC-06) and the
//! patch roles (FC-04).

use mesh_io::handoff::{HandoffMesh, HandoffSummary, Patch, write_foehncraft_vtu};
use serde::Serialize;

use crate::artifact::sha256;
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

/// The FoehnCraft VTU of a mesh, with `provenance` in its field data.
pub(crate) fn foehncraft_vtu(
    mesh: &HandoffMesh,
    provenance: &[(&str, String)],
) -> Result<Vec<u8>, EngineError> {
    let mut metadata = vec![("koflcraft_version", env!("CARGO_PKG_VERSION").to_string())];
    metadata.extend(provenance.iter().cloned());
    let mut vtu = Vec::new();
    write_foehncraft_vtu(mesh, &metadata, &mut vtu).map_err(export_error)?;
    Ok(vtu)
}

/// The evidence of a bundle whose files (other than the evidence itself)
/// are `files`.
pub(crate) fn evidence(
    mesh: &HandoffMesh,
    generator: &str,
    files: &[(&str, &[u8])],
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
        files: files
            .iter()
            .map(|(name, bytes)| (name.to_string(), sha256(bytes)))
            .collect(),
        limitations,
    };
    serde_json::to_vec_pretty(&evidence)
        .map_err(|e| EngineError::new(ErrorCode::SerializationFailed, e))
}
