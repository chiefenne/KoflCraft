use mesh_core::mesh::Mesh;
use mesh_quality::{RELATIVE_CONSERVATION_TOLERANCE, RELATIVE_GEOMETRY_TOLERANCE};
use serde::Serialize;

use crate::artifact::{Point, point};
use crate::{
    ALGORITHM, BoxBounds, BuildEvidence, EngineError, MeshArtifact, MeshConfig, SCHEMA_VERSION,
};

#[derive(Debug, Serialize)]
pub struct DomainReport {
    pub schema_version: u32,
    pub domain: BoxBounds,
    pub site_count: usize,
    pub fixture_id: Option<String>,
    pub seed: Option<u64>,
    pub algorithm: &'static str,
    pub limitations: [&'static str; 2],
}

pub fn inspect_domain(config: &MeshConfig) -> Result<DomainReport, EngineError> {
    config.checked_domain()?;
    Ok(DomainReport {
        schema_version: SCHEMA_VERSION,
        domain: config.domain,
        site_count: config.sites.len(),
        fixture_id: config.fixture_id.clone(),
        seed: config.seed,
        algorithm: ALGORITHM,
        limitations: [
            "axis-aligned box",
            "exhaustive CPU clipping with floating-point tolerances",
        ],
    })
}

#[derive(Debug, Serialize)]
pub struct MeshCheck {
    pub schema_version: u32,
    pub valid: bool,
    pub vertices: usize,
    pub faces: usize,
    pub internal_faces: usize,
    pub boundary_faces: usize,
    pub cells: usize,
    pub total_volume: f64,
    pub domain_volume: f64,
    pub relative_volume_error: f64,
    pub conservation_tolerance: f64,
    pub geometry_tolerance: f64,
    pub topology_digest: String,
}

#[derive(Debug, Serialize)]
pub struct MeshReport {
    pub schema_version: u32,
    pub fixture_id: Option<String>,
    pub seed: Option<u64>,
    pub build: BuildEvidence,
    pub validation: MeshCheck,
    /// Ordered by cell index in the artifact.
    pub cells: Vec<CellReport>,
    /// Ordered by face index in the artifact.
    pub faces: Vec<FaceReport>,
}

#[derive(Debug, Serialize)]
pub struct CellReport {
    pub centroid: Point,
    pub volume: f64,
    pub face_count: usize,
    pub neighbour_count: usize,
}

#[derive(Debug, Serialize)]
pub struct FaceReport {
    pub centroid: Point,
    pub normal: Point,
    pub area: f64,
    pub non_orthogonality_degrees: f64,
    pub skewness: f64,
}

/// Validate a saved artifact without generating candidate cells again.
pub fn check_mesh(artifact: &MeshArtifact) -> Result<MeshCheck, EngineError> {
    let (mesh, metrics) = artifact.verified_mesh()?;
    Ok(summary(artifact, &mesh, &metrics))
}

pub fn mesh_report(artifact: &MeshArtifact) -> Result<MeshReport, EngineError> {
    let (mesh, metrics) = artifact.verified_mesh()?;
    let validation = summary(artifact, &mesh, &metrics);
    Ok(MeshReport {
        schema_version: SCHEMA_VERSION,
        fixture_id: artifact.config().fixture_id.clone(),
        seed: artifact.config().seed,
        build: artifact.build().clone(),
        validation,
        cells: metrics
            .cells
            .into_iter()
            .map(|cell| CellReport {
                centroid: point(cell.centroid),
                volume: cell.volume,
                face_count: cell.face_count,
                neighbour_count: cell.neighbour_count,
            })
            .collect(),
        faces: metrics
            .faces
            .into_iter()
            .map(|face| FaceReport {
                centroid: point(face.centroid),
                normal: point(face.normal),
                area: face.area,
                non_orthogonality_degrees: face.non_orthogonality_degrees,
                skewness: face.skewness,
            })
            .collect(),
    })
}

fn summary(artifact: &MeshArtifact, mesh: &Mesh, metrics: &mesh_quality::MeshMetrics) -> MeshCheck {
    MeshCheck {
        schema_version: SCHEMA_VERSION,
        valid: true,
        vertices: mesh.vertices().len(),
        faces: mesh.face_count(),
        internal_faces: mesh.face_count() - metrics.boundary_face_count,
        boundary_faces: metrics.boundary_face_count,
        cells: mesh.cell_count(),
        total_volume: metrics.total_volume,
        domain_volume: metrics.domain_volume,
        relative_volume_error: (metrics.total_volume / metrics.domain_volume - 1.0).abs(),
        conservation_tolerance: RELATIVE_CONSERVATION_TOLERANCE,
        geometry_tolerance: RELATIVE_GEOMETRY_TOLERANCE,
        topology_digest: artifact.build().topology_digest.clone(),
    }
}
