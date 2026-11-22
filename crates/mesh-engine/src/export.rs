use std::ffi::OsStr;
use std::str::FromStr;

use mesh_io::{ExportError, PolyMeshWriters, write_gmsh_domain, write_openfoam, write_vtu_domain};
use mesh_quality::MeshDomain;

use crate::{EngineError, ErrorCode, MeshArtifact, MeshDirectory};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExportFormat {
    Vtu,
    Gmsh,
    OpenFoam,
}

impl FromStr for ExportFormat {
    type Err = EngineError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "vtu" => Ok(Self::Vtu),
            "gmsh" => Ok(Self::Gmsh),
            _ => Err(EngineError::new(
                ErrorCode::UnsupportedFormat,
                "format must be vtu or gmsh",
            )),
        }
    }
}

pub enum ExportData {
    Vtu(Vec<u8>),
    Gmsh(Vec<u8>),
    OpenFoam(OpenFoamFiles),
}

impl ExportData {
    pub fn write_new(&self, dir: &MeshDirectory, name: &OsStr) -> Result<(), EngineError> {
        match self {
            Self::Vtu(bytes) | Self::Gmsh(bytes) => dir.write_new(name, bytes),
            Self::OpenFoam(files) => dir.write_bundle(name, &files.files()),
        }
        .map_err(|error| EngineError::new(ErrorCode::IoError, error))
    }
}

#[derive(Default)]
pub struct OpenFoamFiles {
    pub points: Vec<u8>,
    pub faces: Vec<u8>,
    pub owner: Vec<u8>,
    pub neighbour: Vec<u8>,
    pub boundary: Vec<u8>,
}

impl OpenFoamFiles {
    pub fn files(&self) -> [(&'static str, &[u8]); 5] {
        [
            ("points", &self.points),
            ("faces", &self.faces),
            ("owner", &self.owner),
            ("neighbour", &self.neighbour),
            ("boundary", &self.boundary),
        ]
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
        ExportFormat::Gmsh => {
            let mut msh = Vec::new();
            write_gmsh_domain(domain, &mesh, &mut msh).map_err(export_error)?;
            ExportData::Gmsh(msh)
        }
        ExportFormat::OpenFoam => {
            let MeshDomain::Box(domain) = domain else {
                return Err(EngineError::new(
                    ErrorCode::UnsupportedFormat,
                    "surface meshes support vtu and gmsh",
                ));
            };
            let mut files = OpenFoamFiles::default();
            write_openfoam(
                domain,
                &mesh,
                PolyMeshWriters {
                    points: &mut files.points,
                    faces: &mut files.faces,
                    owner: &mut files.owner,
                    neighbour: &mut files.neighbour,
                    boundary: &mut files.boundary,
                },
            )
            .map_err(export_error)?;
            ExportData::OpenFoam(files)
        }
    };
    Ok(data)
}

fn export_error(error: ExportError) -> EngineError {
    let code = match error {
        ExportError::Validation(_) => ErrorCode::ValidationFailed,
        ExportError::Io(_) | ExportError::InvalidSubdivision { .. } => ErrorCode::ExportFailed,
    };
    EngineError::new(code, error)
}
