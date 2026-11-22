//! Reconcile candidate cells into a single mesh with globally shared faces and vertices.

use std::collections::BTreeMap;
use std::fmt;

use mesh_algorithms::clipping::{ConvexPolyhedron, FaceSource};
use mesh_algorithms::voronoi::{VoronoiError, clipped_surface_voronoi, clipped_voronoi};
use mesh_core::geometry::{AxisAlignedBox, Vec3};
use mesh_core::mesh::{ConnectivityError, FaceAdjacency, Mesh};
use mesh_geometry::ConvexSurface;
use mesh_quality::{MeshDomain, MeshMetrics, ValidationError, validate_domain};

/// Relative to the box's longest side.
pub const RELATIVE_WELD_TOLERANCE: f64 = 2048.0 * f64::EPSILON;

#[derive(Debug, Clone, PartialEq)]
pub struct GeneratedMesh {
    mesh: Mesh,
    metrics: MeshMetrics,
}

impl GeneratedMesh {
    pub fn mesh(&self) -> &Mesh {
        &self.mesh
    }

    pub fn metrics(&self) -> &MeshMetrics {
        &self.metrics
    }
}

/// Generate a certified, convex-cell Voronoi mesh in input site order.
pub fn generate_mesh(
    domain: AxisAlignedBox,
    sites: &[Vec3],
) -> Result<GeneratedMesh, GenerationError> {
    let cells = clipped_voronoi(domain, sites).map_err(GenerationError::Clip)?;
    certify(MeshDomain::Box(domain), &cells)
}

pub fn generate_surface_mesh(
    surface: &ConvexSurface,
    sites: &[Vec3],
) -> Result<GeneratedMesh, GenerationError> {
    let cells = clipped_surface_voronoi(surface, sites).map_err(GenerationError::Clip)?;
    certify(MeshDomain::Surface(surface), &cells)
}

fn certify(
    domain: MeshDomain<'_>,
    cells: &[ConvexPolyhedron],
) -> Result<GeneratedMesh, GenerationError> {
    let mesh = reconcile(domain.bounds(), cells)?;
    let metrics = validate_domain(domain, &mesh).map_err(GenerationError::Validation)?;
    Ok(GeneratedMesh { mesh, metrics })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenerationError {
    Clip(VoronoiError),
    AmbiguousVertex { cell: usize },
    UnmatchedFace { first: usize, second: usize },
    Connectivity(ConnectivityError),
    Validation(ValidationError),
}

impl fmt::Display for GenerationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Clip(error) => error.fmt(f),
            Self::AmbiguousVertex { cell } => {
                write!(f, "cell {cell} has vertices that cannot be uniquely welded")
            }
            Self::UnmatchedFace { first, second } => write!(
                f,
                "cells {first} and {second} disagree on their shared face"
            ),
            Self::Connectivity(error) => error.fmt(f),
            Self::Validation(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for GenerationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Clip(error) => Some(error),
            Self::Connectivity(error) => Some(error),
            Self::Validation(error) => Some(error),
            _ => None,
        }
    }
}

fn reconcile(domain: AxisAlignedBox, cells: &[ConvexPolyhedron]) -> Result<Mesh, GenerationError> {
    let mut vertices = Vertices {
        scale: domain.extent().max_abs_component(),
        origin: domain.min(),
        points: Vec::new(),
        buckets: BTreeMap::new(),
    };
    let mut face_vertices = Vec::new();
    let mut face_offsets = vec![0];
    let mut adjacency = Vec::new();
    let mut pending = BTreeMap::new();
    for (cell, polyhedron) in cells.iter().enumerate() {
        let mut local_ids = Vec::with_capacity(polyhedron.vertices().len());
        for &point in polyhedron.vertices() {
            let global = vertices
                .insert(point)
                .ok_or(GenerationError::AmbiguousVertex { cell })?;
            if local_ids.contains(&global) {
                return Err(GenerationError::AmbiguousVertex { cell });
            }
            local_ids.push(global);
        }
        for face in polyhedron.faces() {
            let ids = face.vertices().iter().map(|&vertex| local_ids[vertex]);
            let face_adjacency = match face.source() {
                FaceSource::Boundary(patch) => FaceAdjacency::Boundary { owner: cell, patch },
                FaceSource::Surface { facet, patch } => FaceAdjacency::Surface {
                    owner: cell,
                    patch,
                    facet,
                },
                FaceSource::Bisector { other_site } if cell < other_site => {
                    if pending
                        .insert((cell, other_site), adjacency.len())
                        .is_some()
                    {
                        return Err(GenerationError::UnmatchedFace {
                            first: cell,
                            second: other_site,
                        });
                    }
                    FaceAdjacency::Interior {
                        owner: cell,
                        neighbour: other_site,
                    }
                }
                FaceSource::Bisector { other_site } => {
                    let error = GenerationError::UnmatchedFace {
                        first: other_site,
                        second: cell,
                    };
                    let owner_face = pending.remove(&(other_site, cell)).ok_or(error)?;
                    let owner_ids =
                        &face_vertices[face_offsets[owner_face]..face_offsets[owner_face + 1]];
                    if !reverse_ring(owner_ids, ids) {
                        return Err(error);
                    }
                    continue;
                }
                FaceSource::Cut => unreachable!("Voronoi candidates label all clipping planes"),
            };
            adjacency.push(face_adjacency);
            face_vertices.extend(ids);
            face_offsets.push(face_vertices.len());
        }
    }
    if let Some(((first, second), _)) = pending.pop_first() {
        return Err(GenerationError::UnmatchedFace { first, second });
    }
    Mesh::from_faces(
        vertices.points,
        face_offsets,
        face_vertices,
        adjacency,
        cells.len(),
    )
    .map_err(GenerationError::Connectivity)
}

fn reverse_ring(owner: &[usize], mut other: impl ExactSizeIterator<Item = usize>) -> bool {
    if owner.len() != other.len() {
        return false;
    }
    let first = other.next().unwrap();
    let Some(start) = owner.iter().position(|&id| id == first) else {
        return false;
    };
    other.eq(owner[..start]
        .iter()
        .rev()
        .chain(owner[start + 1..].iter().rev())
        .copied())
}

struct Vertices {
    origin: Vec3,
    scale: f64,
    points: Vec<Vec3>,
    buckets: BTreeMap<[i64; 3], Vec<usize>>,
}

impl Vertices {
    fn insert(&mut self, point: Vec3) -> Option<usize> {
        let tol = RELATIVE_WELD_TOLERANCE;
        let relative = (point - self.origin) / self.scale;
        let [x, y, z] =
            [relative.x, relative.y, relative.z].map(|value| (value / tol).floor() as i64);
        let mut found = None;
        for i in x - 1..=x + 1 {
            for j in y - 1..=y + 1 {
                for k in z - 1..=z + 1 {
                    for &id in self.buckets.get(&[i, j, k]).into_iter().flatten() {
                        if ((point - self.points[id]) / self.scale).length() <= tol
                            && found.replace(id).is_some()
                        {
                            return None;
                        }
                    }
                }
            }
        }
        Some(found.unwrap_or_else(|| {
            let id = self.points.len();
            self.points.push(point);
            self.buckets.entry([x, y, z]).or_default().push(id);
            id
        }))
    }
}
