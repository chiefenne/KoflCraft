//! Indexed convex cells and half-space intersection.

use std::collections::{BTreeMap, btree_map::Entry};
use std::fmt;

use mesh_core::geometry::{AxisAlignedBox, BoxSide, FacetId, PatchId, Plane, Vec3};
use mesh_geometry::ConvexSurface;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FaceSource {
    Boundary(BoxSide),
    Surface { facet: FacetId, patch: PatchId },
    Bisector { other_site: usize },
    Cut,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Polygon {
    vertices: Vec<usize>,
    source: FaceSource,
}

impl Polygon {
    /// Local vertex indices wound counterclockwise as seen from outside the cell.
    pub fn vertices(&self) -> &[usize] {
        &self.vertices
    }

    pub fn source(&self) -> FaceSource {
        self.source
    }
}

/// A convex cell with local vertex and face connectivity.
#[derive(Clone, Debug, PartialEq)]
pub struct ConvexPolyhedron {
    vertices: Vec<Vec3>,
    faces: Vec<Polygon>,
}

impl ConvexPolyhedron {
    pub(crate) fn from_surface(surface: &ConvexSurface) -> Self {
        Self {
            vertices: surface.normalized_vertices().to_vec(),
            faces: surface
                .facets()
                .iter()
                .enumerate()
                .map(|(index, facet)| Polygon {
                    vertices: facet.triangle.vertices.to_vec(),
                    source: FaceSource::Surface {
                        facet: FacetId(index),
                        patch: facet.triangle.patch,
                    },
                })
                .collect(),
        }
    }

    pub fn from_box(bounds: AxisAlignedBox) -> Self {
        let min = bounds.min();
        let max = bounds.max();
        let vertices = vec![
            min,
            Vec3::new(max.x, min.y, min.z),
            Vec3::new(max.x, max.y, min.z),
            Vec3::new(min.x, max.y, min.z),
            Vec3::new(min.x, min.y, max.z),
            Vec3::new(max.x, min.y, max.z),
            max,
            Vec3::new(min.x, max.y, max.z),
        ];
        let faces = [
            (BoxSide::XMin, [0, 4, 7, 3]),
            (BoxSide::XMax, [1, 2, 6, 5]),
            (BoxSide::YMin, [0, 1, 5, 4]),
            (BoxSide::YMax, [3, 7, 6, 2]),
            (BoxSide::ZMin, [0, 3, 2, 1]),
            (BoxSide::ZMax, [4, 5, 6, 7]),
        ]
        .into_iter()
        .map(|(side, vertices)| Polygon {
            vertices: vertices.to_vec(),
            source: FaceSource::Boundary(side),
        })
        .collect();
        Self { vertices, faces }
    }

    pub fn vertices(&self) -> &[Vec3] {
        &self.vertices
    }

    pub fn faces(&self) -> &[Polygon] {
        &self.faces
    }

    /// Retain the negative half-space of `plane`, using an absolute distance tolerance.
    ///
    /// Vertices within the tolerance are treated as lying on the plane. Returns
    /// `None` if there is an outside vertex and no strictly inside vertex.
    pub fn clip(
        mut self,
        plane: Plane,
        source: FaceSource,
        tolerance: f64,
    ) -> Result<Option<Self>, ClipError> {
        if !tolerance.is_finite() || tolerance < 0.0 {
            return Err(ClipError::InvalidTolerance);
        }
        let distances: Vec<_> = self
            .vertices
            .iter()
            .map(|&vertex| plane.signed_distance(vertex))
            .collect();
        if distances.iter().any(|d| !d.is_finite()) {
            return Err(ClipError::NonfiniteDistance);
        }
        if !distances.iter().any(|&d| d > tolerance) {
            return Ok(Some(self));
        }
        if !distances.iter().any(|&d| d < -tolerance) {
            return Ok(None);
        }

        let mut intersections = BTreeMap::new();
        let mut on_plane: Vec<_> = distances.iter().map(|d| d.abs() <= tolerance).collect();
        let mut retained = Vec::with_capacity(self.faces.len() + 1);
        for face in self.faces {
            if !face.vertices.iter().any(|&v| distances[v] < -tolerance) {
                continue;
            }
            let mut clipped = Vec::with_capacity(face.vertices.len() + 1);
            for (a, b) in edges(&face.vertices) {
                let da = distances[a];
                let db = distances[b];
                if da <= tolerance {
                    clipped.push(a);
                }
                if (da < -tolerance && db > tolerance) || (db < -tolerance && da > tolerance) {
                    // Both incident faces use the same intersection on their shared edge.
                    let intersection = intersections.entry(edge_key(a, b)).or_insert_with(|| {
                        let scale = da.abs().max(db.abs());
                        let fraction = (da / scale) / (da / scale - db / scale);
                        let point =
                            self.vertices[a] + (self.vertices[b] - self.vertices[a]) * fraction;
                        let index = self.vertices.len();
                        self.vertices.push(point);
                        on_plane.push(true);
                        index
                    });
                    clipped.push(*intersection);
                }
            }
            retained.push(Polygon {
                vertices: clipped,
                source: face.source,
            });
        }
        retained.push(Polygon {
            vertices: cut_boundary(&retained, &on_plane)?,
            source,
        });

        let mut indices = vec![None; self.vertices.len()];
        let mut vertices = Vec::new();
        for face in &mut retained {
            for index in &mut face.vertices {
                *index = *indices[*index].get_or_insert_with(|| {
                    let mapped = vertices.len();
                    vertices.push(self.vertices[*index]);
                    mapped
                });
            }
        }
        if vertices.iter().any(|vertex| !vertex.is_finite()) {
            return Err(ClipError::NonfiniteIntersection);
        }
        Ok(Some(Self {
            vertices,
            faces: retained,
        }))
    }

    pub(crate) fn map_vertices(&mut self, map: impl Fn(Vec3) -> Vec3) {
        for vertex in &mut self.vertices {
            *vertex = map(*vertex);
        }
    }
}

fn edges(vertices: &[usize]) -> impl Iterator<Item = (usize, usize)> + '_ {
    vertices
        .iter()
        .copied()
        .zip(vertices.iter().copied().cycle().skip(1))
}

fn edge_key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

fn cut_boundary(faces: &[Polygon], on_plane: &[bool]) -> Result<Vec<usize>, ClipError> {
    let mut exposed = BTreeMap::new();
    for face in faces {
        for (a, b) in edges(&face.vertices).filter(|&(a, b)| on_plane[a] && on_plane[b]) {
            // Opposite incidences cancel. Reverse the exposed edges for cap winding.
            match exposed.entry(edge_key(a, b)) {
                Entry::Vacant(entry) => {
                    entry.insert((b, a));
                }
                Entry::Occupied(entry) if *entry.get() == (a, b) => {
                    entry.remove();
                }
                _ => return Err(ClipError::InvalidCutBoundary),
            }
        }
    }

    let edge_count = exposed.len();
    let mut next: BTreeMap<_, _> = exposed.into_values().collect();
    if next.len() != edge_count || edge_count < 3 {
        return Err(ClipError::InvalidCutBoundary);
    }
    let (&first, _) = next
        .first_key_value()
        .ok_or(ClipError::InvalidCutBoundary)?;
    let mut boundary = Vec::with_capacity(edge_count);
    let mut vertex = first;
    for _ in 0..edge_count {
        boundary.push(vertex);
        vertex = next.remove(&vertex).ok_or(ClipError::InvalidCutBoundary)?;
    }
    if vertex != first {
        return Err(ClipError::InvalidCutBoundary);
    }
    Ok(boundary)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ClipError {
    InvalidTolerance,
    NonfiniteDistance,
    NonfiniteIntersection,
    InvalidCutBoundary,
}

impl fmt::Display for ClipError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidTolerance => "clipping tolerance must be finite and nonnegative",
            Self::NonfiniteDistance => "a vertex has a nonfinite distance to the clipping plane",
            Self::NonfiniteIntersection => "an edge intersection has nonfinite coordinates",
            Self::InvalidCutBoundary => "the cut does not form a single closed polygon",
        })
    }
}

impl std::error::Error for ClipError {}
