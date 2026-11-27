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
#[derive(Clone, Debug, PartialEq)]
pub enum Curve {
    Line,
    /// The circular arc from the first vertex through `through` to the second.
    Arc {
        through: Vec3,
    },
    /// The natural cubic spline (chord-length parameter) from the first
    /// vertex through the points `through`, in order, to the second.
    Spline {
        through: Vec<Vec3>,
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
#[derive(Clone, Debug, PartialEq)]
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

/// An analytic surface that boundary block faces are projected onto.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Surface {
    /// The cylinder of `radius` about the line through `point` along `axis`.
    Cylinder {
        point: Vec3,
        axis: Vec3,
        radius: f64,
    },
    /// The torus about the axis through `centre` along `axis`: the circle
    /// of radius `major` in the plane normal to `axis`, swept by a circle of
    /// radius `minor`.
    Torus {
        centre: Vec3,
        axis: Vec3,
        major: f64,
        minor: f64,
    },
}

fn unit(v: Vec3) -> Option<Vec3> {
    let l = v.length();
    (l > 0.0 && l.is_finite()).then(|| v / l)
}

impl Surface {
    /// Why the surface is unusable, if it is.
    pub fn problem(&self) -> Option<&'static str> {
        let ok = |v: f64| v.is_finite() && v > 0.0;
        match *self {
            Self::Cylinder {
                point,
                axis,
                radius,
            } => {
                if !point.is_finite() || unit(axis).is_none() || !ok(radius) {
                    return Some(
                        "a cylinder needs a finite point, a non-zero axis and a positive radius",
                    );
                }
            }
            Self::Torus {
                centre,
                axis,
                major,
                minor,
            } => {
                if !centre.is_finite() || unit(axis).is_none() || !ok(major) || !ok(minor) {
                    return Some(
                        "a torus needs a finite centre, a non-zero axis and positive radii",
                    );
                }
                if minor >= major {
                    return Some("a torus needs its minor radius below its major radius");
                }
            }
        }
        None
    }

    /// The length scale of the surface (its smallest radius).
    pub fn scale(&self) -> f64 {
        match *self {
            Self::Cylinder { radius, .. } => radius,
            Self::Torus { minor, .. } => minor,
        }
    }

    /// The closest point of the surface; a point on the axis (cylinder) or
    /// on the core circle (torus) has no unique one and is returned as is.
    pub fn project(&self, p: Vec3) -> Vec3 {
        match *self {
            Self::Cylinder {
                point,
                axis,
                radius,
            } => {
                let Some(a) = unit(axis) else { return p };
                let d = p - point;
                let foot = point + a * d.dot(a);
                match unit(p - foot) {
                    Some(r) => foot + r * radius,
                    None => p,
                }
            }
            Self::Torus {
                centre,
                axis,
                major,
                minor,
            } => {
                let Some(a) = unit(axis) else { return p };
                let d = p - centre;
                let Some(radial) = unit(d - a * d.dot(a)) else {
                    return p;
                };
                let ring = centre + radial * major;
                match unit(p - ring) {
                    Some(r) => ring + r * minor,
                    None => p,
                }
            }
        }
    }
}

/// A named surface and the boundary block faces whose nodes are projected
/// onto it (BS-07).
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceSpec {
    pub name: String,
    pub surface: Surface,
    /// Block faces by their four vertices (any order of the four).
    pub faces: Vec<[usize; 4]>,
}

/// A T-junction (BS-21): the block edge `edge` passes through the
/// vertices `through`, in order from `edge[0]` to `edge[1]`. The blocks on
/// one side carry the whole edge; on the other side the pieces between the
/// vertices are edges of smaller blocks. Its cell count is the sum of the
/// pieces' counts, and the grids still match node for node.
#[derive(Clone, Debug, PartialEq)]
pub struct JunctionSpec {
    pub edge: [usize; 2],
    pub through: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Default)]
pub struct BlockTopology {
    pub vertices: Vec<Vec3>,
    pub edges: Vec<EdgeSpec>,
    pub blocks: Vec<BlockSpec>,
    pub patches: Vec<PatchSpec>,
    pub surfaces: Vec<SurfaceSpec>,
    /// T-junctions on block edges; empty for a topology whose blocks share
    /// faces whole.
    pub junctions: Vec<JunctionSpec>,
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

    #[test]
    fn projections_land_on_their_surface() {
        let torus = Surface::Torus {
            centre: Vec3::new(1.0, 2.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 2.0),
            major: 3.0,
            minor: 0.5,
        };
        let p = torus.project(Vec3::new(4.3, 2.1, 0.2));
        let rho = ((p.x - 1.0).powi(2) + (p.y - 2.0).powi(2)).sqrt();
        assert!((((rho - 3.0).powi(2) + p.z * p.z).sqrt() - 0.5).abs() < 1e-14);
        let cyl = Surface::Cylinder {
            point: Vec3::ZERO,
            axis: Vec3::new(1.0, 0.0, 0.0),
            radius: 2.0,
        };
        let q = cyl.project(Vec3::new(5.0, 0.3, -0.4));
        assert!((q.x - 5.0).abs() < 1e-15 && ((q.y * q.y + q.z * q.z).sqrt() - 2.0).abs() < 1e-15);
        assert!(torus.problem().is_none() && cyl.problem().is_none());
    }
}
