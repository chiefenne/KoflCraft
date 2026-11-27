mod artifact;
mod assembly;
pub mod block;
pub mod block_edit;
pub mod block_triangle;
pub mod channel;
mod export;
pub mod fill;
pub mod handoff;
pub mod patches;
pub mod ports;
pub mod proposals;
pub mod regions;
mod report;
pub mod sizing;
pub mod skeleton;
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
pub use export::{Export, ExportData, ExportFormat, Written, export_mesh, prepare_export};
pub use mesh_io::{FILE_HASH, GmshVersion, MeshDirectory, file_digest};
pub use report::{
    CellReport, DomainReport, FaceReport, MeshCheck, MeshReport, check_mesh, inspect_domain,
    mesh_report,
};
pub use stl::{
    StlDomainReport, StlMeshConfig, StlPatchReport, build_stl_mesh, inspect_stl_domain,
    stl_site_count,
};

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
    InvalidSizeField,
    InvalidArguments,
    InvalidPatchMap,
    InvalidBlockTopology,
    UnknownPatch,
    InvalidPatchName,
    EmptySelection,
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
