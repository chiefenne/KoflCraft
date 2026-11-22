//! Packed connectivity for a globally indexed finite-volume mesh.

use std::fmt;

use crate::geometry::{BoxSide, FacetId, PatchId, Vec3};

/// A face is stored once, with its vertices wound outward from its owner.
/// The neighbour sees the reverse winding.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FaceAdjacency {
    Interior {
        owner: usize,
        neighbour: usize,
    },
    Boundary {
        owner: usize,
        patch: BoxSide,
    },
    Surface {
        owner: usize,
        patch: PatchId,
        facet: FacetId,
    },
}

impl FaceAdjacency {
    pub fn owner(self) -> usize {
        match self {
            Self::Interior { owner, .. }
            | Self::Boundary { owner, .. }
            | Self::Surface { owner, .. } => owner,
        }
    }

    pub fn neighbour(self) -> Option<usize> {
        match self {
            Self::Interior { neighbour, .. } => Some(neighbour),
            Self::Boundary { .. } | Self::Surface { .. } => None,
        }
    }
}

/// Immutable connectivity with offset arrays for the variable-length face and
/// cell lists. Vertex, face, and cell IDs are zero-based.
#[derive(Clone, Debug, PartialEq)]
pub struct Mesh {
    vertices: Vec<Vec3>,
    face_offsets: Vec<usize>,
    face_vertices: Vec<usize>,
    adjacency: Vec<FaceAdjacency>,
    cell_offsets: Vec<usize>,
    cell_faces: Vec<usize>,
}

impl Mesh {
    /// Check packed indices and derive cell incidence from face ownership.
    ///
    /// Offsets include a leading zero and a terminal offset equal to the length
    /// of the vertex-index array. Winding, manifoldness, and spatial partition
    /// require geometric validation before the mesh is used for CFD.
    pub fn from_faces(
        vertices: Vec<Vec3>,
        face_offsets: Vec<usize>,
        face_vertices: Vec<usize>,
        adjacency: Vec<FaceAdjacency>,
        cell_count: usize,
    ) -> Result<Self, ConnectivityError> {
        if face_offsets.len().checked_sub(1) != Some(adjacency.len())
            || face_offsets.first() != Some(&0)
            || face_offsets.last() != Some(&face_vertices.len())
        {
            return Err(ConnectivityError::FaceOffsets);
        }
        if cell_count == 0 || cell_count > adjacency.len() / 2 {
            return Err(ConnectivityError::CellCount);
        }
        if let Some(vertex) = vertices.iter().position(|v| !v.is_finite()) {
            return Err(ConnectivityError::NonfiniteVertex { vertex });
        }

        let mut cell_sizes = vec![0; cell_count];
        for (face, (&incidence, range)) in adjacency.iter().zip(face_offsets.windows(2)).enumerate()
        {
            let indices = face_vertices
                .get(range[0]..range[1])
                .filter(|indices| indices.len() >= 3)
                .ok_or(ConnectivityError::FaceVertices { face })?;
            for (i, &vertex) in indices.iter().enumerate() {
                if vertex >= vertices.len() || indices[..i].contains(&vertex) {
                    return Err(ConnectivityError::FaceVertices { face });
                }
            }
            for cell in incident_cells(incidence) {
                let size = cell_sizes
                    .get_mut(cell)
                    .ok_or(ConnectivityError::FaceCells { face })?;
                *size += 1;
            }
            if incidence.neighbour() == Some(incidence.owner()) {
                return Err(ConnectivityError::FaceCells { face });
            }
        }

        let mut cell_offsets = vec![0; cell_count + 1];
        for (cell, size) in cell_sizes.into_iter().enumerate() {
            if size < 4 {
                return Err(ConnectivityError::CellFaces { cell });
            }
            cell_offsets[cell + 1] = cell_offsets[cell] + size;
        }
        let mut cell_faces = vec![0; cell_offsets[cell_count]];
        let mut next = cell_offsets[..cell_count].to_vec();
        for (face, &incidence) in adjacency.iter().enumerate() {
            for cell in incident_cells(incidence) {
                cell_faces[next[cell]] = face;
                next[cell] += 1;
            }
        }
        Ok(Self {
            vertices,
            face_offsets,
            face_vertices,
            adjacency,
            cell_offsets,
            cell_faces,
        })
    }

    pub fn vertices(&self) -> &[Vec3] {
        &self.vertices
    }

    pub fn face_count(&self) -> usize {
        self.adjacency.len()
    }

    pub fn cell_count(&self) -> usize {
        self.cell_offsets.len() - 1
    }

    pub fn face_vertices(&self, face: usize) -> Option<&[usize]> {
        let start = *self.face_offsets.get(face)?;
        let end = *self.face_offsets.get(face + 1)?;
        Some(&self.face_vertices[start..end])
    }

    pub fn face_adjacency(&self, face: usize) -> Option<FaceAdjacency> {
        self.adjacency.get(face).copied()
    }

    pub fn cell_faces(&self, cell: usize) -> Option<&[usize]> {
        let start = *self.cell_offsets.get(cell)?;
        let end = *self.cell_offsets.get(cell + 1)?;
        Some(&self.cell_faces[start..end])
    }
}

fn incident_cells(adjacency: FaceAdjacency) -> impl Iterator<Item = usize> {
    std::iter::once(adjacency.owner()).chain(adjacency.neighbour())
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ConnectivityError {
    FaceOffsets,
    CellCount,
    NonfiniteVertex { vertex: usize },
    FaceVertices { face: usize },
    FaceCells { face: usize },
    CellFaces { cell: usize },
}

impl fmt::Display for ConnectivityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::FaceOffsets => f.write_str("face offsets do not describe the packed face array"),
            Self::CellCount => f.write_str("cell count is inconsistent with the number of faces"),
            Self::NonfiniteVertex { vertex } => write!(f, "vertex {vertex} is not finite"),
            Self::FaceVertices { face } => write!(f, "face {face} has invalid vertex indices"),
            Self::FaceCells { face } => write!(f, "face {face} has invalid cell incidences"),
            Self::CellFaces { cell } => write!(f, "cell {cell} has fewer than four faces"),
        }
    }
}

impl std::error::Error for ConnectivityError {}
