//! Data model of the block-structured variant (BS-01, BS-02): a block
//! topology of vertices, curved edges with node distributions, hexahedral
//! blocks with cell counts, and named boundary patches made of block faces.
//!
//! The file format `koflcraft.block-topology` (SW-20) is read and written by
//! `mesh-engine`; this module holds only the in-memory model.
//!
//! Local block frame (as in the classic multi-block literature): the eight
//! vertices of a block are the corners (i, j, k) ∈ {0, 1}³ in the order
//! (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1).
//! A block must be right-handed: (e_i × e_j) · e_k > 0 at every corner.

use crate::geometry::Vec3;

/// The shape of a block edge between its two vertices.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Curve {
    Line,
    /// The circular arc from the first vertex through `through` to the second.
    Arc {
        through: Vec3,
    },
}

/// Where a one-sided spacing is prescribed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EdgeEnd {
    /// At the first vertex of the edge as listed.
    Start,
    /// At the second vertex.
    End,
}

/// The node distribution along an edge (BS-06).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Distribution {
    Uniform,
    /// Geometric growth from the first cell of length `spacing` at `at`.
    Geometric {
        spacing: f64,
        at: EdgeEnd,
    },
    /// Two-sided hyperbolic tangent (Vinokur 1983) with the first and last
    /// cell lengths given.
    Tanh {
        first: f64,
        last: f64,
    },
}

/// An edge given explicitly; every other block edge is a uniform line.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EdgeSpec {
    pub vertices: [usize; 2],
    pub curve: Curve,
    pub distribution: Distribution,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BlockSpec {
    pub vertices: [usize; 8],
    /// Cells along the local i, j and k axes.
    pub cells: [usize; 3],
    pub domain: String,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PatchSpec {
    pub name: String,
    /// Block faces by their four vertices (any order of the four).
    pub faces: Vec<[usize; 4]>,
}

#[derive(Clone, Debug, PartialEq, Default)]
pub struct BlockTopology {
    pub vertices: Vec<Vec3>,
    pub edges: Vec<EdgeSpec>,
    pub blocks: Vec<BlockSpec>,
    pub patches: Vec<PatchSpec>,
}

/// The corner (i, j, k) ∈ {0, 1}³ of the local block frame as a vertex slot.
pub const fn corner_slot(i: usize, j: usize, k: usize) -> usize {
    // (0,0)→0, (1,0)→1, (1,1)→2, (0,1)→3 in the i-j plane
    let base = match (i, j) {
        (0, 0) => 0,
        (1, 0) => 1,
        (1, 1) => 2,
        _ => 3,
    };
    base + 4 * k
}

/// The twelve block edges as (axis, the corner of the edge's start). The
/// start corner has 0 on `axis`.
pub fn block_edges() -> [(usize, [usize; 3]); 12] {
    let mut out = [(0, [0; 3]); 12];
    let mut n = 0;
    for axis in 0..3 {
        let (a, b) = ((axis + 1) % 3, (axis + 2) % 3);
        for cb in 0..2 {
            for ca in 0..2 {
                let mut c = [0; 3];
                c[a] = ca;
                c[b] = cb;
                out[n] = (axis, c);
                n += 1;
            }
        }
    }
    out
}

/// The six block faces as (fixed axis, side, first free axis, second free
/// axis); the free axes are in increasing order.
pub const BLOCK_FACES: [(usize, usize, usize, usize); 6] = [
    (0, 0, 1, 2),
    (0, 1, 1, 2),
    (1, 0, 0, 2),
    (1, 1, 0, 2),
    (2, 0, 0, 1),
    (2, 1, 0, 1),
];

/// The vertex slot of a corner given by the fixed axis and side of a face
/// and the corner's coordinates on its two free axes.
pub fn face_corner_slot(face: (usize, usize, usize, usize), ca: usize, cb: usize) -> usize {
    let (fixed, side, a, b) = face;
    let mut c = [0; 3];
    c[fixed] = side;
    c[a] = ca;
    c[b] = cb;
    corner_slot(c[0], c[1], c[2])
}

/// The key of a block face: its four vertices sorted.
pub fn face_key(mut vertices: [usize; 4]) -> [usize; 4] {
    vertices.sort_unstable();
    vertices
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn corners_follow_the_local_frame() {
        assert_eq!(corner_slot(0, 0, 0), 0);
        assert_eq!(corner_slot(1, 1, 0), 2);
        assert_eq!(corner_slot(0, 1, 1), 7);
        let edges = block_edges();
        assert_eq!(edges.len(), 12);
        assert_eq!(edges.iter().filter(|e| e.0 == 1).count(), 4);
    }
}
