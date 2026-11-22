//! Independent geometry and topology checks for a convex volume mesh.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use mesh_core::geometry::{AxisAlignedBox, BoxSide, FacetId, Vec3};
use mesh_core::mesh::{FaceAdjacency, Mesh};
use mesh_geometry::ConvexSurface;

/// Relative to the box's longest side.
pub const RELATIVE_GEOMETRY_TOLERANCE: f64 = 8192.0 * f64::EPSILON;
/// Relative to the expected box volume, patch area, or cell surface area.
pub const RELATIVE_CONSERVATION_TOLERANCE: f64 = 1e-10;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FaceMetrics {
    pub centroid: Vec3,
    /// Unit normal pointing out of the owner.
    pub normal: Vec3,
    pub area: f64,
    pub non_orthogonality_degrees: f64,
    /// Normalized in-plane distance from the face centroid to the centroid-line intersection.
    /// At a boundary the intersection is the orthogonal projection of the owner centroid.
    pub skewness: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CellMetrics {
    pub centroid: Vec3,
    pub volume: f64,
    pub face_count: usize,
    pub neighbour_count: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MeshMetrics {
    pub cells: Vec<CellMetrics>,
    pub faces: Vec<FaceMetrics>,
    pub boundary_face_count: usize,
    pub total_volume: f64,
    pub domain_volume: f64,
}

#[derive(Clone, Copy)]
pub enum MeshDomain<'a> {
    Box(AxisAlignedBox),
    Surface(&'a ConvexSurface),
}

impl MeshDomain<'_> {
    pub fn bounds(self) -> AxisAlignedBox {
        match self {
            Self::Box(bounds) => bounds,
            Self::Surface(surface) => surface.bounds(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValidationError {
    DomainScale,
    Vertex { vertex: usize, reason: &'static str },
    Face { face: usize, reason: &'static str },
    Cell { cell: usize, reason: &'static str },
    BoundaryClosure,
    PatchArea { patch: BoxSide },
    FacetArea { facet: FacetId },
    VolumeConservation,
}

impl fmt::Display for ValidationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DomainScale => {
                f.write_str("domain area or volume is outside the representable range")
            }
            Self::Vertex { vertex, reason } => write!(f, "vertex {vertex}: {reason}"),
            Self::Face { face, reason } => write!(f, "face {face}: {reason}"),
            Self::Cell { cell, reason } => write!(f, "cell {cell}: {reason}"),
            Self::BoundaryClosure => f.write_str("boundary edges must have two opposite uses"),
            Self::PatchArea { patch } => write!(f, "{patch:?} area disagrees with the box"),
            Self::FacetArea { facet } => {
                write!(f, "boundary area disagrees with surface facet {}", facet.0)
            }
            Self::VolumeConservation => f.write_str("total cell volume disagrees with the domain"),
        }
    }
}

impl std::error::Error for ValidationError {}

/// Check a packed mesh against its prescribed domain and compute finite-volume metrics.
/// All geometry checks run in box-relative coordinates.
pub fn validate_mesh(domain: AxisAlignedBox, mesh: &Mesh) -> Result<MeshMetrics, ValidationError> {
    validate_domain(MeshDomain::Box(domain), mesh)
}

pub fn validate_domain(
    domain: MeshDomain<'_>,
    mesh: &Mesh,
) -> Result<MeshMetrics, ValidationError> {
    let bounds = domain.bounds();
    let origin = bounds.min();
    let scale = bounds.extent().max_abs_component();
    let extent = bounds.extent() / scale;
    let area_scale = scale * scale;
    let volume_scale = area_scale * scale;
    let (expected_volume, boundary_slots) = match domain {
        MeshDomain::Box(_) => (extent.x * extent.y * extent.z, BoxSide::ALL.len()),
        MeshDomain::Surface(surface) => (surface.normalized_volume(), surface.facets().len()),
    };
    let domain_volume = expected_volume * volume_scale;
    if !domain_volume.is_normal() || !volume_scale.is_normal() {
        return Err(ValidationError::DomainScale);
    }
    let tol = RELATIVE_GEOMETRY_TOLERANCE;
    let points: Vec<_> = mesh
        .vertices()
        .iter()
        .map(|&point| (point - origin) / scale)
        .collect();
    for (vertex, point) in points.iter().enumerate() {
        if !point.is_finite()
            || point.x < -tol
            || point.x > extent.x + tol
            || point.y < -tol
            || point.y > extent.y + tol
            || point.z < -tol
            || point.z > extent.z + tol
        {
            return Err(ValidationError::Vertex {
                vertex,
                reason: "outside the domain bounds",
            });
        }
        if let MeshDomain::Surface(surface) = domain
            && surface
                .facets()
                .iter()
                .any(|facet| facet.plane.signed_distance(*point) > tol)
        {
            return Err(ValidationError::Vertex {
                vertex,
                reason: "outside the convex surface",
            });
        }
    }

    let mut faces = Vec::with_capacity(mesh.face_count());
    let mut used = vec![false; points.len()];
    let mut face_keys = BTreeSet::new();
    let mut boundary_edges = EdgeUses::new();
    let mut patch_areas = vec![0.0; boundary_slots];
    let mut boundary_face_count = 0;
    for face in 0..mesh.face_count() {
        let ids = mesh.face_vertices(face).unwrap();
        let mut key = ids.to_vec();
        key.sort_unstable();
        if !face_keys.insert(key) {
            return Err(ValidationError::Face {
                face,
                reason: "duplicate polygon",
            });
        }
        for &id in ids {
            used[id] = true;
        }
        let metrics = face_metrics(&points, ids, face)?;
        let adjacency = mesh.face_adjacency(face).unwrap();
        if adjacency.neighbour().is_none() {
            let (slot, facing) = boundary_projection(
                domain,
                extent,
                adjacency,
                face,
                metrics.normal,
                ids,
                &points,
            )?;
            patch_areas[slot] += metrics.area * facing;
            boundary_face_count += 1;
            add_edges(&mut boundary_edges, ids, false);
        }
        faces.push(metrics);
    }
    if let Some(vertex) = used.iter().position(|&is_used| !is_used) {
        return Err(ValidationError::Vertex {
            vertex,
            reason: "unused by any face",
        });
    }
    if !closed(&boundary_edges) {
        return Err(ValidationError::BoundaryClosure);
    }
    match domain {
        MeshDomain::Box(_) => {
            let expected_areas = [
                extent.y * extent.z,
                extent.x * extent.z,
                extent.x * extent.y,
            ];
            for patch in BoxSide::ALL {
                if !conserved(patch_areas[patch as usize], expected_areas[patch.axis()]) {
                    return Err(ValidationError::PatchArea { patch });
                }
            }
        }
        MeshDomain::Surface(surface) => {
            for (index, (facet, area)) in surface.facets().iter().zip(patch_areas).enumerate() {
                if !conserved(area, facet.normalized_area) {
                    return Err(ValidationError::FacetArea {
                        facet: FacetId(index),
                    });
                }
            }
        }
    }

    let mut cells = Vec::with_capacity(mesh.cell_count());
    for cell in 0..mesh.cell_count() {
        cells.push(cell_metrics(mesh, &points, &faces, cell)?);
    }
    let total_volume = cells.iter().map(|cell| cell.volume).sum::<f64>();
    let world_volume = total_volume * volume_scale;
    if !world_volume.is_finite() || !conserved(total_volume, expected_volume) {
        return Err(ValidationError::VolumeConservation);
    }
    for (face, metrics) in faces.iter_mut().enumerate() {
        let adjacency = mesh.face_adjacency(face).unwrap();
        let owner = cells[adjacency.owner()].centroid;
        let to_face = metrics.centroid - owner;
        let delta = adjacency
            .neighbour()
            .map_or(to_face, |neighbour| cells[neighbour].centroid - owner);
        let projection = delta.dot(metrics.normal);
        let face_distance = to_face.dot(metrics.normal);
        if face_distance <= 0.0 || (adjacency.neighbour().is_some() && projection <= face_distance)
        {
            return Err(ValidationError::Face {
                face,
                reason: "inconsistent cell centroid orientation",
            });
        }
        metrics.non_orthogonality_degrees = delta
            .cross(metrics.normal)
            .length()
            .atan2(projection)
            .to_degrees();
        let offset = if adjacency.neighbour().is_some() {
            to_face - delta * (face_distance / projection)
        } else {
            to_face - metrics.normal * projection
        };
        metrics.skewness = offset.length() / delta.length();
        metrics.centroid = origin + metrics.centroid * scale;
        metrics.area *= area_scale;
        if !metrics.centroid.is_finite()
            || !metrics.area.is_finite()
            || metrics.area <= 0.0
            || !metrics.skewness.is_finite()
        {
            return Err(ValidationError::Face {
                face,
                reason: "metrics outside the representable range",
            });
        }
    }
    for (cell, metrics) in cells.iter_mut().enumerate() {
        metrics.centroid = origin + metrics.centroid * scale;
        metrics.volume *= volume_scale;
        if !metrics.centroid.is_finite() || !metrics.volume.is_finite() || metrics.volume <= 0.0 {
            return Err(ValidationError::Cell {
                cell,
                reason: "metrics outside the representable range",
            });
        }
    }
    Ok(MeshMetrics {
        cells,
        faces,
        boundary_face_count,
        total_volume: world_volume,
        domain_volume,
    })
}

fn boundary_projection(
    domain: MeshDomain<'_>,
    extent: Vec3,
    adjacency: FaceAdjacency,
    face: usize,
    normal: Vec3,
    ids: &[usize],
    points: &[Vec3],
) -> Result<(usize, f64), ValidationError> {
    let tol = RELATIVE_GEOMETRY_TOLERANCE;
    let error = ValidationError::Face {
        face,
        reason: "incorrect boundary patch, position, or winding",
    };
    match (domain, adjacency) {
        (MeshDomain::Box(_), FaceAdjacency::Boundary { patch, .. }) => {
            let axis = patch.axis();
            let component = [normal.x, normal.y, normal.z][axis];
            let (facing, offset) = if patch.is_max() {
                (component, [extent.x, extent.y, extent.z][axis])
            } else {
                (-component, 0.0)
            };
            if facing <= 0.0
                || ids.iter().any(|&id| {
                    let Vec3 { x, y, z } = points[id];
                    ([x, y, z][axis] - offset).abs() > tol
                })
            {
                return Err(error);
            }
            Ok((patch as usize, facing))
        }
        (MeshDomain::Surface(surface), FaceAdjacency::Surface { patch, facet, .. }) => {
            let boundary = surface.facets().get(facet.0).ok_or(error)?;
            let plane = boundary.plane;
            let facing = normal.dot(plane.normal());
            if patch != boundary.triangle.patch
                || facing <= 0.0
                || ids
                    .iter()
                    .any(|&id| plane.signed_distance(points[id]).abs() > tol)
            {
                return Err(error);
            }
            let [a, b, c] = boundary
                .triangle
                .vertices
                .map(|id| surface.normalized_vertices()[id]);
            for (start, end) in [(a, b), (b, c), (c, a)] {
                let edge = end - start;
                let inward = plane.normal().cross(edge / edge.length());
                if ids
                    .iter()
                    .any(|&id| (points[id] - start).dot(inward) < -tol)
                {
                    return Err(error);
                }
            }
            Ok((facet.0, facing))
        }
        _ => Err(error),
    }
}

fn face_metrics(
    points: &[Vec3],
    ids: &[usize],
    face: usize,
) -> Result<FaceMetrics, ValidationError> {
    let error = |reason| ValidationError::Face { face, reason };
    let origin = points[ids[0]];
    let triangles = || {
        ids[1..]
            .windows(2)
            .map(|pair| (points[pair[0]] - origin, points[pair[1]] - origin))
    };
    let area_vector = triangles().fold(Vec3::ZERO, |sum, (b, c)| sum + b.cross(c)) * 0.5;
    let area = area_vector.length();
    if area <= 0.0 {
        return Err(error("zero area"));
    }
    let normal = area_vector / area;
    let tol = RELATIVE_GEOMETRY_TOLERANCE;
    if ids
        .iter()
        .any(|&id| (points[id] - origin).dot(normal).abs() > tol)
    {
        return Err(error("polygon is not planar"));
    }
    for (a, b) in ring_edges(ids) {
        let start = points[a];
        let edge = points[b] - start;
        let length = edge.length();
        if length == 0.0 {
            return Err(error("zero-length edge"));
        }
        let inward = normal.cross(edge / length);
        if ids
            .iter()
            .any(|&id| (points[id] - start).dot(inward) < -tol)
        {
            return Err(error("polygon is not convex and consistently wound"));
        }
    }
    let moment = triangles().fold(Vec3::ZERO, |sum, (b, c)| {
        sum + (b + c) * (b.cross(c).dot(normal) / 6.0)
    });
    Ok(FaceMetrics {
        centroid: origin + moment / area,
        normal,
        area,
        non_orthogonality_degrees: 0.0,
        skewness: 0.0,
    })
}

fn cell_metrics(
    mesh: &Mesh,
    points: &[Vec3],
    faces: &[FaceMetrics],
    cell: usize,
) -> Result<CellMetrics, ValidationError> {
    let error = |reason| ValidationError::Cell { cell, reason };
    let cell_faces = mesh.cell_faces(cell).unwrap();
    let vertices: BTreeSet<_> = cell_faces
        .iter()
        .flat_map(|&face| mesh.face_vertices(face).unwrap())
        .copied()
        .collect();
    let origin = vertices
        .iter()
        .fold(Vec3::ZERO, |sum, &id| sum + points[id])
        / vertices.len() as f64;
    let mut edges = EdgeUses::new();
    let mut area_vector = Vec3::ZERO;
    let mut area = 0.0;
    let mut volume = 0.0;
    let mut moment = Vec3::ZERO;
    let mut neighbours = BTreeSet::new();
    for &face in cell_faces {
        let metrics = faces[face];
        let ids = mesh.face_vertices(face).unwrap();
        let adjacency = mesh.face_adjacency(face).unwrap();
        let reversed = adjacency.owner() != cell;
        let orientation = if reversed { -1.0 } else { 1.0 };
        let outward = metrics.normal * orientation;
        let on_face = points[ids[0]];
        if (origin - on_face).dot(outward) >= 0.0
            || vertices
                .iter()
                .any(|&id| (points[id] - on_face).dot(outward) > RELATIVE_GEOMETRY_TOLERANCE)
        {
            return Err(error("face points inward or cell is not convex"));
        }
        add_edges(&mut edges, ids, reversed);
        area_vector = area_vector + outward * metrics.area;
        area += metrics.area;
        if let Some(neighbour) = adjacency.neighbour() {
            neighbours.insert(if reversed {
                adjacency.owner()
            } else {
                neighbour
            });
        }
        let a = on_face - origin;
        for pair in ids[1..].windows(2) {
            let b = points[pair[0]] - origin;
            let c = points[pair[1]] - origin;
            let tetra_volume = orientation * a.dot(b.cross(c)) / 6.0;
            volume += tetra_volume;
            moment = moment + (a + b + c) * (tetra_volume / 4.0);
        }
    }
    if !closed(&edges) || area_vector.length() > area * RELATIVE_CONSERVATION_TOLERANCE {
        return Err(error("surface does not close"));
    }
    if volume <= 0.0 || !volume.is_finite() {
        return Err(error("volume must be positive and finite"));
    }
    Ok(CellMetrics {
        centroid: origin + moment / volume,
        volume,
        face_count: cell_faces.len(),
        neighbour_count: neighbours.len(),
    })
}

fn conserved(actual: f64, expected: f64) -> bool {
    expected > 0.0 && (actual / expected - 1.0).abs() <= RELATIVE_CONSERVATION_TOLERANCE
}

type EdgeUses = BTreeMap<(usize, usize), [usize; 2]>;

fn ring_edges(ids: &[usize]) -> impl Iterator<Item = (usize, usize)> + '_ {
    ids.iter().copied().zip(ids.iter().copied().cycle().skip(1))
}

fn add_edges(uses: &mut EdgeUses, ids: &[usize], reversed: bool) {
    for (a, b) in ring_edges(ids) {
        uses.entry((a.min(b), a.max(b))).or_default()[usize::from((a > b) ^ reversed)] += 1;
    }
}

fn closed(uses: &EdgeUses) -> bool {
    uses.values().all(|&count| count == [1, 1])
}
