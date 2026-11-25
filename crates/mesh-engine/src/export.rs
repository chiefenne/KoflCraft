use std::ffi::OsStr;
use std::io::{self, Write};
use std::str::FromStr;

use mesh_core::geometry::BoxSide;
use mesh_core::mesh::Mesh;
use mesh_io::handoff::{
    HandoffMesh, POLYMESH_FILES, PatchRole, patches_from_names, write_polymesh, write_polymesh_file,
};
use mesh_io::{ExportError, GmshVersion, StreamError, write_gmsh_domain, write_vtu_domain};
use mesh_quality::{MeshDomain, validate_domain};

use crate::artifact::ArtifactDomain;
use crate::handoff::HandoffFormat;
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
    /// The same mesh, names and digests as a blocked VTKHDF file, binary
    /// and readable block by block (`docs/foehncraft.md`).
    FoehnCraftVtkhdf,
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
            "foehncraft-vtkhdf" => Ok(Self::FoehnCraftVtkhdf),
            "openfoam" => Ok(Self::OpenFoam),
            _ => Err(EngineError::new(
                ErrorCode::UnsupportedFormat,
                "format must be vtu, gmsh (MSH 4.1), gmsh-2.2, foehncraft, foehncraft-vtkhdf or openfoam",
            )),
        }
    }
}

/// An export in memory (see [`Export::to_memory`]).
pub enum ExportData {
    Vtu(Vec<u8>),
    Gmsh(Vec<u8>),
    FoehnCraft(Vec<u8>),
    /// A directory of named files.
    Bundle(Vec<(&'static str, Vec<u8>)>),
}

/// What [`Export::write_new`] wrote, with the BLAKE3 hash of every file.
pub enum Written {
    /// One file at the output path.
    File(String),
    /// A new directory: its files in order, each with its digest.
    Bundle(Vec<(String, String)>),
}

/// A validated export of a saved mesh, not yet written: streamed to disk by
/// [`Export::write_new`] (one block of each file in memory at a time) or
/// gathered by [`Export::to_memory`].
pub struct Export {
    format: ExportFormat,
    domain: ArtifactDomain,
    mesh: Mesh,
    /// The face-based mesh of the FoehnCraft and OpenFOAM formats.
    handoff: Option<HandoffMesh>,
    /// Field data of the FoehnCraft VTU (FC-05).
    provenance: Vec<(&'static str, String)>,
    algorithm: String,
}

/// Validate a saved mesh and render its export through the shared format writers.
pub fn export_mesh(
    artifact: &MeshArtifact,
    format: ExportFormat,
) -> Result<ExportData, EngineError> {
    prepare_export(artifact, format)?.to_memory()
}

/// Validate a saved mesh for an export and build what its format needs.
pub fn prepare_export(
    artifact: &MeshArtifact,
    format: ExportFormat,
) -> Result<Export, EngineError> {
    let (domain, mesh) = artifact.verified_connectivity()?;
    let handoff = match format {
        ExportFormat::FoehnCraft | ExportFormat::FoehnCraftVtkhdf | ExportFormat::OpenFoam => Some(
            artifact_handoff(domain.as_domain(), &mesh, artifact.patch_roles())?,
        ),
        ExportFormat::Vtu | ExportFormat::Gmsh(_) => None,
    };
    let build = artifact.build();
    Ok(Export {
        format,
        domain,
        mesh,
        handoff,
        provenance: vec![
            ("koflcraft_mesh_digest", build.mesh_digest.clone()),
            ("koflcraft_topology_digest", build.topology_digest.clone()),
            ("koflcraft_input_digest", build.input_digest.clone()),
            ("koflcraft_algorithm", build.algorithm.clone()),
        ],
        algorithm: build.algorithm.clone(),
    })
}

impl Export {
    /// Write the export as the new file (or, for `openfoam`, the new
    /// directory) `name` in `dir`, each file streamed and hashed as it is
    /// written. Nothing is left behind on an error.
    pub fn write_new(&self, dir: &MeshDirectory, name: &OsStr) -> Result<Written, EngineError> {
        if self.format == ExportFormat::OpenFoam {
            let handoff = self.handoff()?;
            let fail = crate::handoff::stream_error;
            let mut bundle = dir
                .bundle(name)
                .map_err(|e| EngineError::new(ErrorCode::IoError, e))?;
            for file in POLYMESH_FILES {
                bundle
                    .stream(file, |w| write_polymesh_file(handoff, file, w))
                    .map_err(fail)?;
            }
            let evidence = crate::handoff::evidence(handoff, &self.algorithm, bundle.files())?;
            bundle
                .stream(crate::handoff::MESH_EVIDENCE, |w| w.write_all(&evidence))
                .map_err(fail)?;
            return Ok(Written::Bundle(bundle.finish()));
        }
        let mut failure = None;
        let written = dir.stream_new(name, |w| self.write_single(w, &mut failure));
        match (written, failure) {
            (Ok(digest), _) => Ok(Written::File(digest)),
            (Err(StreamError::Content(_)), Some(error)) => Err(export_error(error)),
            (Err(error), _) => Err(crate::handoff::stream_error(error)),
        }
    }

    /// The export in memory, for tests and small meshes.
    pub fn to_memory(&self) -> Result<ExportData, EngineError> {
        if self.format == ExportFormat::OpenFoam {
            let handoff = self.handoff()?;
            let mut files = write_polymesh(handoff).map_err(crate::handoff::export_error)?;
            let digests: Vec<(String, String)> = files
                .iter()
                .map(|(n, b)| (n.to_string(), mesh_io::file_digest(b)))
                .collect();
            let evidence = crate::handoff::evidence(handoff, &self.algorithm, &digests)?;
            files.push((crate::handoff::MESH_EVIDENCE, evidence));
            return Ok(ExportData::Bundle(files));
        }
        let mut bytes = Vec::new();
        let mut failure = None;
        if let Err(error) = self.write_single(&mut bytes, &mut failure) {
            return Err(match failure {
                Some(error) => export_error(error),
                None => crate::handoff::export_error(error),
            });
        }
        Ok(match self.format {
            ExportFormat::Vtu => ExportData::Vtu(bytes),
            ExportFormat::Gmsh(_) => ExportData::Gmsh(bytes),
            _ => ExportData::FoehnCraft(bytes),
        })
    }

    fn handoff(&self) -> Result<&HandoffMesh, EngineError> {
        self.handoff.as_ref().ok_or_else(|| {
            EngineError::new(ErrorCode::ExportFailed, "the export has no hand-off mesh")
        })
    }

    /// The file of a single-file format into `out`. A writer error that is
    /// not plain I/O (validation, subdivision) is kept in `failure`, so it
    /// keeps its diagnostic code.
    fn write_single(
        &self,
        out: &mut dyn Write,
        failure: &mut Option<ExportError>,
    ) -> io::Result<()> {
        let mut keep = |result: Result<(), ExportError>| match result {
            Ok(()) => Ok(()),
            Err(ExportError::Io(e)) => Err(e),
            Err(error) => {
                let e = io::Error::other(error.to_string());
                *failure = Some(error);
                Err(e)
            }
        };
        let domain = self.domain.as_domain();
        match self.format {
            ExportFormat::Vtu => keep(write_vtu_domain(domain, &self.mesh, out)),
            ExportFormat::Gmsh(version) => {
                keep(write_gmsh_domain(domain, &self.mesh, version, out))
            }
            ExportFormat::FoehnCraft | ExportFormat::FoehnCraftVtkhdf => {
                let handoff_format = if self.format == ExportFormat::FoehnCraft {
                    HandoffFormat::Vtu
                } else {
                    HandoffFormat::Vtkhdf
                };
                match &self.handoff {
                    Some(handoff) => handoff_format.write(handoff, &self.provenance, out),
                    None => Err(io::Error::other("the export has no hand-off mesh")),
                }
            }
            ExportFormat::OpenFoam => Err(io::Error::other(
                "openfoam is a directory, not a single file",
            )),
        }
    }
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
