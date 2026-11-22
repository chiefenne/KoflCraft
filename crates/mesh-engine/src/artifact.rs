use std::ffi::OsStr;
use std::fmt::Write;
use std::time::Instant;

use mesh_algorithms::voronoi::{validate_sites, validate_surface_sites};
use mesh_core::geometry::{AxisAlignedBox, BoxSide, FacetId, PatchId, Vec3};
use mesh_core::mesh::{FaceAdjacency, Mesh};
use mesh_geometry::{ConvexSurface, SurfaceTriangle};
use mesh_io::handoff::PatchRole;
use mesh_quality::{MeshDomain, MeshMetrics, validate_domain};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::sizing::SizingEvidence;
use crate::{EngineError, ErrorCode, GeneratedMesh, GenerationError, MeshDirectory, generate_mesh};

pub const SCHEMA_VERSION: u32 = 1;
pub const ALGORITHM: &str = "cpu_box_clipped_voronoi_v1";
pub const SURFACE_ARTIFACT_VERSION: u32 = 2;
pub(crate) const SURFACE_ALGORITHM: &str = "cpu_convex_surface_clipped_voronoi_v1";

pub(crate) type Point = [f64; 3];

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BoxBounds {
    pub min: Point,
    pub max: Point,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MeshConfig {
    pub schema_version: u32,
    pub domain: BoxBounds,
    pub sites: Vec<Point>,
    pub fixture_id: Option<String>,
    pub seed: Option<u64>,
}

impl MeshConfig {
    fn input(&self) -> Result<(AxisAlignedBox, Vec<Vec3>), EngineError> {
        require_schema(self.schema_version)?;
        let domain = AxisAlignedBox::new(vector(self.domain.min), vector(self.domain.max))
            .ok_or_else(|| {
                EngineError::new(
                    ErrorCode::InvalidDomain,
                    "box bounds must be finite and increasing",
                )
            })?;
        let sites = self.sites.iter().copied().map(vector).collect();
        Ok((domain, sites))
    }

    pub(crate) fn checked_domain(&self) -> Result<AxisAlignedBox, EngineError> {
        let (domain, sites) = self.input()?;
        validate_sites(domain, &sites)
            .map_err(|error| EngineError::new(ErrorCode::InvalidSites, error))?;
        Ok(domain)
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BuildEvidence {
    pub algorithm: String,
    pub engine_version: String,
    pub operating_system: String,
    pub architecture: String,
    pub debug_assertions: bool,
    pub elapsed_micros: u64,
    pub input_digest: String,
    pub config_digest: String,
    pub topology_digest: String,
    pub mesh_digest: String,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MeshArtifact {
    schema_version: u32,
    config: MeshConfig,
    build: BuildEvidence,
    mesh: PackedMesh,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    surface: Option<PackedSurface>,
    /// How a size field placed the sites, and the sizes achieved (GRD-08).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    sizing: Option<SizingEvidence>,
}

pub(crate) enum ArtifactDomain {
    Box(AxisAlignedBox),
    Surface(ConvexSurface),
}

impl ArtifactDomain {
    pub(crate) fn as_domain(&self) -> MeshDomain<'_> {
        match self {
            Self::Box(bounds) => MeshDomain::Box(*bounds),
            Self::Surface(surface) => MeshDomain::Surface(surface),
        }
    }
}

impl MeshArtifact {
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    pub fn surface_facet_count(&self) -> Option<usize> {
        self.surface.as_ref().map(|surface| surface.facets.len())
    }

    pub fn config(&self) -> &MeshConfig {
        &self.config
    }

    pub fn build(&self) -> &BuildEvidence {
        &self.build
    }

    /// The patch roles saved from a patch map (BC-04), in the surface's
    /// patch order.
    pub fn patch_roles(&self) -> Option<&[PatchRole]> {
        self.surface.as_ref()?.patch_roles.as_deref()
    }

    pub fn sizing(&self) -> Option<&SizingEvidence> {
        self.sizing.as_ref()
    }

    pub fn write_new(&self, dir: &MeshDirectory, name: &OsStr) -> Result<(), EngineError> {
        let mut data = serde_json::to_vec_pretty(self)
            .map_err(|error| EngineError::new(ErrorCode::SerializationFailed, error))?;
        data.push(b'\n');
        dir.write_new(name, &data)
            .map_err(|error| EngineError::new(ErrorCode::IoError, error))
    }

    pub(crate) fn verified_mesh(&self) -> Result<(Mesh, MeshMetrics), EngineError> {
        let (domain, mesh) = self.verified_connectivity()?;
        let metrics = validate_domain(domain.as_domain(), &mesh)
            .map_err(|error| EngineError::new(ErrorCode::ValidationFailed, error))?;
        Ok((mesh, metrics))
    }

    pub(crate) fn verified_connectivity(&self) -> Result<(ArtifactDomain, Mesh), EngineError> {
        let expected_schema = if self.surface.is_some() {
            SURFACE_ARTIFACT_VERSION
        } else {
            SCHEMA_VERSION
        };
        if self.schema_version != expected_schema {
            return Err(EngineError::new(
                ErrorCode::UnsupportedSchema,
                "artifact version does not match its domain",
            ));
        }
        if self.build.input_digest != input_digest(&self.config, self.surface.as_ref())?
            || self.build.config_digest
                != config_digest(&self.config, self.surface.as_ref(), self.sizing.as_ref())?
            || self.build.topology_digest != digest(&self.mesh.topology)?
            || self.build.mesh_digest != digest(&self.mesh)?
        {
            return Err(EngineError::new(
                ErrorCode::DigestMismatch,
                "artifact digests do not match its config or mesh",
            ));
        }
        let domain = match &self.surface {
            None => ArtifactDomain::Box(self.config.checked_domain()?),
            Some(packed) => {
                let (bounds, sites) = self.config.input()?;
                let surface = packed.to_surface()?;
                if bounds != surface.bounds() {
                    return Err(EngineError::new(
                        ErrorCode::InvalidDomain,
                        "artifact bounds differ from the saved surface",
                    ));
                }
                validate_surface_sites(&surface, &sites)
                    .map_err(|error| EngineError::new(ErrorCode::InvalidSites, error))?;
                ArtifactDomain::Surface(surface)
            }
        };
        let mesh = self.mesh.to_mesh()?;
        if mesh.cell_count() != self.config.sites.len() {
            return Err(EngineError::new(
                ErrorCode::InvalidConnectivity,
                "cell count differs from the configured site count",
            ));
        }
        Ok((domain, mesh))
    }
}

/// Build a versioned artifact with a certified mesh and reproducibility evidence.
pub fn build_mesh(config: MeshConfig) -> Result<MeshArtifact, EngineError> {
    let (domain, sites) = config.input()?;
    let started = Instant::now();
    let generated = generate_mesh(domain, &sites).map_err(generation_error)?;
    pack_artifact(config, &generated, None, None, started)
}

pub(crate) fn generation_error(error: GenerationError) -> EngineError {
    let code = match error {
        GenerationError::Clip(
            mesh_algorithms::voronoi::VoronoiError::NoSites
            | mesh_algorithms::voronoi::VoronoiError::SiteOutside { .. }
            | mesh_algorithms::voronoi::VoronoiError::SitesTooClose { .. }
            | mesh_algorithms::voronoi::VoronoiError::ThinDomain,
        ) => ErrorCode::InvalidSites,
        GenerationError::Validation(_) => ErrorCode::ValidationFailed,
        _ => ErrorCode::GenerationFailed,
    };
    EngineError::new(code, error)
}

pub(crate) fn pack_artifact(
    config: MeshConfig,
    generated: &GeneratedMesh,
    surface: Option<(&ConvexSurface, Option<Vec<PatchRole>>)>,
    sizing: Option<SizingEvidence>,
    started: Instant,
) -> Result<MeshArtifact, EngineError> {
    let elapsed_micros = u64::try_from(started.elapsed().as_micros()).unwrap_or(u64::MAX);
    let (schema_version, algorithm) = if surface.is_some() {
        (SURFACE_ARTIFACT_VERSION, SURFACE_ALGORITHM)
    } else {
        (SCHEMA_VERSION, ALGORITHM)
    };
    let surface = surface.map(|(surface, roles)| PackedSurface {
        patch_roles: roles,
        ..PackedSurface::from(surface)
    });
    let mesh = PackedMesh::from(generated.mesh());
    let build = BuildEvidence {
        algorithm: algorithm.to_owned(),
        engine_version: env!("CARGO_PKG_VERSION").to_owned(),
        operating_system: std::env::consts::OS.to_owned(),
        architecture: std::env::consts::ARCH.to_owned(),
        debug_assertions: cfg!(debug_assertions),
        elapsed_micros,
        input_digest: input_digest(&config, surface.as_ref())?,
        config_digest: config_digest(&config, surface.as_ref(), sizing.as_ref())?,
        topology_digest: digest(&mesh.topology)?,
        mesh_digest: digest(&mesh)?,
    };
    Ok(MeshArtifact {
        schema_version,
        config,
        build,
        mesh,
        surface,
        sizing,
    })
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct PackedSurface {
    vertices: Vec<Point>,
    facets: Vec<PackedFacet>,
    patch_names: Vec<String>,
    /// The role of each patch from a patch map (BC-04); absent when the
    /// roles follow the names, so older artifacts keep their digests.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    patch_roles: Option<Vec<PatchRole>>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct PackedFacet {
    vertices: [usize; 3],
    patch: usize,
}

impl From<&ConvexSurface> for PackedSurface {
    fn from(surface: &ConvexSurface) -> Self {
        Self {
            vertices: surface.vertices().iter().copied().map(point).collect(),
            facets: surface
                .facets()
                .iter()
                .map(|facet| PackedFacet {
                    vertices: facet.triangle.vertices,
                    patch: facet.triangle.patch.0,
                })
                .collect(),
            patch_names: surface.patch_names().to_vec(),
            patch_roles: None,
        }
    }
}

impl PackedSurface {
    fn to_surface(&self) -> Result<ConvexSurface, EngineError> {
        if self
            .patch_roles
            .as_ref()
            .is_some_and(|r| r.len() != self.patch_names.len())
        {
            return Err(EngineError::new(
                ErrorCode::InvalidDomain,
                "the artifact has a different number of patch roles and names",
            ));
        }
        ConvexSurface::new(
            self.vertices.iter().copied().map(vector).collect(),
            self.facets
                .iter()
                .map(|facet| SurfaceTriangle {
                    vertices: facet.vertices,
                    patch: PatchId(facet.patch),
                })
                .collect(),
            self.patch_names.clone(),
        )
        .map_err(|error| EngineError::new(ErrorCode::InvalidDomain, error))
    }
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct PackedMesh {
    vertices: Vec<Point>,
    topology: PackedTopology,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct PackedTopology {
    face_offsets: Vec<usize>,
    face_vertices: Vec<usize>,
    adjacency: Vec<PackedFace>,
    cell_count: usize,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
enum PackedFace {
    Interior {
        owner: usize,
        neighbour: usize,
    },
    Boundary {
        owner: usize,
        patch: String,
    },
    Surface {
        owner: usize,
        patch: usize,
        facet: usize,
    },
}

impl From<&Mesh> for PackedMesh {
    fn from(mesh: &Mesh) -> Self {
        let mut topology = PackedTopology {
            face_offsets: vec![0],
            face_vertices: Vec::new(),
            adjacency: Vec::with_capacity(mesh.face_count()),
            cell_count: mesh.cell_count(),
        };
        for face in 0..mesh.face_count() {
            topology
                .face_vertices
                .extend(mesh.face_vertices(face).unwrap());
            topology.face_offsets.push(topology.face_vertices.len());
            topology
                .adjacency
                .push(match mesh.face_adjacency(face).unwrap() {
                    FaceAdjacency::Interior { owner, neighbour } => {
                        PackedFace::Interior { owner, neighbour }
                    }
                    FaceAdjacency::Boundary { owner, patch } => PackedFace::Boundary {
                        owner,
                        patch: patch.name().to_owned(),
                    },
                    FaceAdjacency::Surface {
                        owner,
                        patch,
                        facet,
                    } => PackedFace::Surface {
                        owner,
                        patch: patch.0,
                        facet: facet.0,
                    },
                });
        }
        Self {
            vertices: mesh.vertices().iter().copied().map(point).collect(),
            topology,
        }
    }
}

impl PackedMesh {
    fn to_mesh(&self) -> Result<Mesh, EngineError> {
        let mut adjacency = Vec::with_capacity(self.topology.adjacency.len());
        for face in &self.topology.adjacency {
            adjacency.push(match face {
                PackedFace::Interior { owner, neighbour } => FaceAdjacency::Interior {
                    owner: *owner,
                    neighbour: *neighbour,
                },
                PackedFace::Boundary { owner, patch } => {
                    let side = BoxSide::ALL
                        .into_iter()
                        .find(|side| side.name() == patch)
                        .ok_or_else(|| {
                            EngineError::new(
                                ErrorCode::InvalidConnectivity,
                                format_args!("unknown boundary patch {patch}"),
                            )
                        })?;
                    FaceAdjacency::Boundary {
                        owner: *owner,
                        patch: side,
                    }
                }
                PackedFace::Surface {
                    owner,
                    patch,
                    facet,
                } => FaceAdjacency::Surface {
                    owner: *owner,
                    patch: PatchId(*patch),
                    facet: FacetId(*facet),
                },
            });
        }
        Mesh::from_faces(
            self.vertices.iter().copied().map(vector).collect(),
            self.topology.face_offsets.clone(),
            self.topology.face_vertices.clone(),
            adjacency,
            self.topology.cell_count,
        )
        .map_err(|error| EngineError::new(ErrorCode::InvalidConnectivity, error))
    }
}

pub(crate) fn require_schema(version: u32) -> Result<(), EngineError> {
    if version != SCHEMA_VERSION {
        return Err(EngineError::new(
            ErrorCode::UnsupportedSchema,
            format_args!("schema version {version} is unsupported"),
        ));
    }
    Ok(())
}

pub(crate) fn point(value: Vec3) -> Point {
    [value.x, value.y, value.z]
}

pub(crate) fn vector([x, y, z]: Point) -> Vec3 {
    Vec3::new(x, y, z)
}

fn input_digest(
    config: &MeshConfig,
    surface: Option<&PackedSurface>,
) -> Result<String, EngineError> {
    match surface {
        Some(surface) => digest(&(surface, &config.sites)),
        None => digest(&(config.domain, &config.sites)),
    }
}

/// The sizing evidence joins the digest only when present, so artifacts
/// built from explicit sites keep their digests.
fn config_digest(
    config: &MeshConfig,
    surface: Option<&PackedSurface>,
    sizing: Option<&SizingEvidence>,
) -> Result<String, EngineError> {
    match (surface, sizing) {
        (Some(surface), Some(sizing)) => digest(&(config, surface, sizing)),
        (Some(surface), None) => digest(&(config, surface)),
        (None, _) => digest(config),
    }
}

pub(crate) fn digest(value: &impl Serialize) -> Result<String, EngineError> {
    let json = serde_json::to_vec(value)
        .map_err(|error| EngineError::new(ErrorCode::SerializationFailed, error))?;
    Ok(sha256(&json))
}

pub(crate) fn sha256(bytes: &[u8]) -> String {
    let mut hex = String::with_capacity(64);
    for byte in Sha256::digest(bytes) {
        write!(hex, "{byte:02x}").unwrap();
    }
    hex
}
