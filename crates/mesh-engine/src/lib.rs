mod artifact;
mod assembly;
pub mod channel;
mod export;
pub mod handoff;
mod report;
mod stl;

use std::fmt;

use serde::Serialize;

pub use artifact::{
    ALGORITHM, BoxBounds, BuildEvidence, MeshArtifact, MeshConfig, SCHEMA_VERSION,
    SURFACE_ARTIFACT_VERSION, build_mesh,
};
pub use assembly::{
    GeneratedMesh, GenerationError, RELATIVE_WELD_TOLERANCE, generate_mesh, generate_surface_mesh,
};
pub use export::{ExportData, ExportFormat, export_mesh};
pub use mesh_io::{GmshVersion, MeshDirectory};
pub use report::{
    CellReport, DomainReport, FaceReport, MeshCheck, MeshReport, check_mesh, inspect_domain,
    mesh_report,
};
pub use stl::{StlDomainReport, StlMeshConfig, StlPatchReport, build_stl_mesh, inspect_stl_domain};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ErrorCode {
    UnsupportedSchema,
    InvalidDomain,
    InvalidSites,
    GenerationFailed,
    InvalidConnectivity,
    DigestMismatch,
    ValidationFailed,
    UnsupportedFormat,
    ExportFailed,
    IoError,
    SerializationFailed,
}

#[derive(Debug, Serialize)]
pub struct EngineError {
    pub code: ErrorCode,
    pub message: String,
}

impl EngineError {
    fn new(code: ErrorCode, message: impl fmt::Display) -> Self {
        Self {
            code,
            message: message.to_string(),
        }
    }
}

impl fmt::Display for EngineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for EngineError {}

pub mod duct;

pub mod sweep;
