use std::ffi::OsStr;
use std::str::FromStr;

use mesh_core::geometry::BoxSide;
use mesh_core::mesh::Mesh;
use mesh_io::handoff::{HandoffMesh, PatchRole, patches_from_names, write_polymesh};
use mesh_io::{ExportError, GmshVersion, write_gmsh_domain, write_vtu_domain};
use mesh_quality::{MeshDomain, validate_domain};

use crate::{EngineError, ErrorCode, MeshArtifact, MeshDirectory};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExportFormat {
    Vtu,
    /// Gmsh MSH in the given version: `gmsh` is 4.1 (FC-01), `gmsh-2.2` the
    /// older version kept as an option (OPEN-11).
    Gmsh(GmshVersion),
    /// The single-file VTU for FoehnCraft: cells and boundary faces with
    /// names and digests (FC-02, FC-05).
    FoehnCraft,
    /// An OpenFOAM `constant/polyMesh` directory (FC-03).
    OpenFoam,
}

impl FromStr for ExportFormat {
    type Err = EngineError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "vtu" => Ok(Self::Vtu),
            "gmsh" | "gmsh-4.1" => Ok(Self::Gmsh(GmshVersion::V4_1)),
            "gmsh-2.2" => Ok(Self::Gmsh(GmshVersion::V2_2)),
            "foehncraft" => Ok(Self::FoehnCraft),
            "openfoam" => Ok(Self::OpenFoam),
            _ => Err(EngineError::new(
                ErrorCode::UnsupportedFormat,
                "format must be vtu, gmsh (MSH 4.1), gmsh-2.2, foehncraft or openfoam",
            )),
        }
    }
}

pub enum ExportData {
    Vtu(Vec<u8>),
    Gmsh(Vec<u8>),
    FoehnCraft(Vec<u8>),
    /// A directory of named files.
    Bundle(Vec<(&'static str, Vec<u8>)>),
}

impl ExportData {
    pub fn write_new(&self, dir: &MeshDirectory, name: &OsStr) -> Result<(), EngineError> {
        match self {
            Self::Vtu(bytes) | Self::Gmsh(bytes) | Self::FoehnCraft(bytes) => {
                dir.write_new(name, bytes)
            }
            Self::Bundle(files) => {
                let files: Vec<(&str, &[u8])> =
                    files.iter().map(|(n, b)| (*n, b.as_slice())).collect();
                dir.write_bundle(name, &files)
            }
        }
        .map_err(|error| EngineError::new(ErrorCode::IoError, error))
    }
}

/// Validate a saved mesh and render its export through the shared format writers.
pub fn export_mesh(
    artifact: &MeshArtifact,
    format: ExportFormat,
) -> Result<ExportData, EngineError> {
    let (domain, mesh) = artifact.verified_connectivity()?;
    let domain = domain.as_domain();
    let data = match format {
        ExportFormat::Vtu => {
            let mut vtu = Vec::new();
            write_vtu_domain(domain, &mesh, &mut vtu).map_err(export_error)?;
            ExportData::Vtu(vtu)
        }
        ExportFormat::Gmsh(version) => {
            let mut msh = Vec::new();
            write_gmsh_domain(domain, &mesh, version, &mut msh).map_err(export_error)?;
            ExportData::Gmsh(msh)
        }
        ExportFormat::FoehnCraft => {
            let handoff = artifact_handoff(domain, &mesh, artifact.patch_roles())?;
            let build = artifact.build();
            ExportData::FoehnCraft(crate::handoff::foehncraft_vtu(
                &handoff,
                &[
                    ("koflcraft_mesh_digest", build.mesh_digest.clone()),
                    ("koflcraft_topology_digest", build.topology_digest.clone()),
                    ("koflcraft_input_digest", build.input_digest.clone()),
                    ("koflcraft_algorithm", build.algorithm.clone()),
                ],
            )?)
        }
        ExportFormat::OpenFoam => {
            let handoff = artifact_handoff(domain, &mesh, artifact.patch_roles())?;
            let mut files = write_polymesh(&handoff).map_err(crate::handoff::export_error)?;
            let all: Vec<(&str, &[u8])> = files.iter().map(|(n, b)| (*n, b.as_slice())).collect();
            let evidence = crate::handoff::evidence(&handoff, &artifact.build().algorithm, &all)?;
            files.push((crate::handoff::MESH_EVIDENCE, evidence));
            ExportData::Bundle(files)
        }
    };
    Ok(data)
}

/// The hand-off mesh; patch roles from a patch map replace the roles the
/// names suggest (BC-04).
fn artifact_handoff(
    domain: MeshDomain<'_>,
    mesh: &Mesh,
    roles: Option<&[PatchRole]>,
) -> Result<HandoffMesh, EngineError> {
    validate_domain(domain, mesh).map_err(|e| EngineError::new(ErrorCode::ValidationFailed, e))?;
    let names: Vec<String> = match domain {
        MeshDomain::Box(_) => BoxSide::ALL.iter().map(|s| s.name().to_string()).collect(),
        MeshDomain::Surface(surface) => surface.patch_names().to_vec(),
    };
    let mut patches = patches_from_names(&names);
    if let Some(roles) = roles {
        for (patch, role) in patches.iter_mut().zip(roles) {
            patch.role = *role;
        }
    }
    HandoffMesh::from_mesh(mesh, patches)
        .map_err(|e| EngineError::new(ErrorCode::InvalidConnectivity, e))
}

fn export_error(error: ExportError) -> EngineError {
    let code = match error {
        ExportError::Validation(_) => ErrorCode::ValidationFailed,
        ExportError::Io(_) | ExportError::InvalidSubdivision { .. } => ErrorCode::ExportFailed,
    };
    EngineError::new(code, error)
}
