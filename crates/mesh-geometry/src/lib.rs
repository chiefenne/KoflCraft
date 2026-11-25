//! Import and validation of triangulated volume boundaries.

pub mod patches;
mod stl;

pub use patches::FacetGeometry;
pub use stl::{StlOptions, import_closed_stl, import_stl};

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::ops::Deref;

use mesh_core::geometry::{AxisAlignedBox, PatchId, Plane, Vec3};

/// Plane tolerance relative to the longest extent of the surface bounding box.
pub const RELATIVE_SURFACE_TOLERANCE: f64 = 8192.0 * f64::EPSILON;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SurfaceTriangle {
    pub vertices: [usize; 3],
    pub patch: PatchId,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Facet {
    pub triangle: SurfaceTriangle,
    /// Oriented triangle plane in bounding-box normalized coordinates.
    pub plane: Plane,
    pub normalized_area: f64,
}

/// One connected triangular shell with paired, oppositely wound edges and positive signed volume.
/// Convexity and geometric intersection checks belong to the conversion to `ConvexSurface`.
#[derive(Clone, Debug, PartialEq)]
pub struct ClosedSurface {
    vertices: Vec<Vec3>,
    normalized_vertices: Vec<Vec3>,
    facets: Vec<Facet>,
    patch_names: Vec<String>,
    bounds: AxisAlignedBox,
    normalized_volume: f64,
}

impl ClosedSurface {
    pub fn new(
        vertices: Vec<Vec3>,
        triangles: Vec<SurfaceTriangle>,
        patch_names: Vec<String>,
    ) -> Result<Self, SurfaceError> {
        Self::build(vertices, triangles, patch_names, true)
    }

    /// [`ClosedSurface::new`] for a domain boundary of several closed
    /// pieces: a meshed domain around a zero-thickness shell that touches
    /// no wall (its two sides form a closed piece of zero volume). Every
    /// other check stays.
    pub fn with_pieces(
        vertices: Vec<Vec3>,
        triangles: Vec<SurfaceTriangle>,
        patch_names: Vec<String>,
    ) -> Result<Self, SurfaceError> {
        Self::build(vertices, triangles, patch_names, false)
    }

    fn build(
        vertices: Vec<Vec3>,
        triangles: Vec<SurfaceTriangle>,
        patch_names: Vec<String>,
        connected: bool,
    ) -> Result<Self, SurfaceError> {
        let first = *vertices.first().ok_or(SurfaceError::Empty)?;
        let mut min = first;
        let mut max = first;
        for (id, &v) in vertices.iter().enumerate() {
            if !v.is_finite() {
                return Err(SurfaceError::InvalidVertex(id));
            }
            min = Vec3::new(min.x.min(v.x), min.y.min(v.y), min.z.min(v.z));
            max = Vec3::new(max.x.max(v.x), max.y.max(v.y), max.z.max(v.z));
        }
        let bounds = AxisAlignedBox::new(min, max).ok_or(SurfaceError::InvalidBounds)?;
        let scale = bounds.extent().max_abs_component();
        let normalized_vertices: Vec<_> = vertices.iter().map(|&v| (v - min) / scale).collect();
        if triangles.len() < 4 {
            return Err(SurfaceError::Empty);
        }
        let mut names = BTreeSet::new();
        for (index, name) in patch_names.iter().enumerate() {
            if name.trim().is_empty()
                || name.chars().any(|c| c.is_control() || c == '"')
                || !names.insert(name)
            {
                return Err(SurfaceError::InvalidPatch(index));
            }
        }

        let mut facets = Vec::with_capacity(triangles.len());
        let mut triangle_keys = BTreeSet::new();
        let mut used_vertices = vec![false; vertices.len()];
        let mut used_patches = vec![false; patch_names.len()];
        let mut edges = BTreeMap::<(usize, usize), [Option<usize>; 2]>::new();
        for (index, triangle) in triangles.into_iter().enumerate() {
            let [a, b, c] = triangle.vertices;
            let mut key = triangle.vertices;
            key.sort_unstable();
            if key[2] >= vertices.len() || key[0] == key[1] || key[1] == key[2] {
                return Err(SurfaceError::DegenerateFacet(index));
            }
            if !triangle_keys.insert(key) {
                return Err(SurfaceError::DuplicateFacet(index));
            }
            *used_patches
                .get_mut(triangle.patch.0)
                .ok_or(SurfaceError::InvalidPatch(triangle.patch.0))? = true;
            for id in triangle.vertices {
                used_vertices[id] = true;
            }
            let [p, q, r] = triangle.vertices.map(|id| normalized_vertices[id]);
            let cross = (q - p).cross(r - p);
            let plane = Plane::new(p, cross).ok_or(SurfaceError::DegenerateFacet(index))?;
            let length = cross.length();
            if length
                <= RELATIVE_SURFACE_TOLERANCE
                    * (q - p).length().max((r - p).length()).max((r - q).length())
            {
                return Err(SurfaceError::DegenerateFacet(index));
            }
            facets.push(Facet {
                triangle,
                plane,
                normalized_area: length / 2.0,
            });
            for (start, end) in [(a, b), (b, c), (c, a)] {
                let edge = (start.min(end), start.max(end));
                let slot = &mut edges.entry(edge).or_default()[usize::from(start > end)];
                if slot.replace(index).is_some() {
                    return Err(SurfaceError::WindingOrNonManifold {
                        first: edge.0,
                        second: edge.1,
                    });
                }
            }
        }
        if let Some(index) = used_vertices.iter().position(|&used| !used) {
            return Err(SurfaceError::InvalidVertex(index));
        }
        if let Some(index) = used_patches.iter().position(|&used| !used) {
            return Err(SurfaceError::InvalidPatch(index));
        }
        let mut neighbours = vec![Vec::new(); facets.len()];
        for ((first, second), pair) in edges {
            let [Some(left), Some(right)] = pair else {
                return Err(SurfaceError::OpenEdge { first, second });
            };
            neighbours[left].push(right);
            neighbours[right].push(left);
        }
        let mut seen = vec![false; facets.len()];
        let mut stack = vec![0];
        seen[0] = true;
        while let Some(index) = stack.pop() {
            for &next in &neighbours[index] {
                if !std::mem::replace(&mut seen[next], true) {
                    stack.push(next);
                }
            }
        }
        if connected && seen.contains(&false) {
            return Err(SurfaceError::Disconnected);
        }

        let origin = vertex_mean(&normalized_vertices);
        let normalized_volume: f64 = facets
            .iter()
            .map(|facet| {
                let [a, b, c] = facet
                    .triangle
                    .vertices
                    .map(|id| normalized_vertices[id] - origin);
                a.dot(b.cross(c)) / 6.0
            })
            .sum();
        let volume_scale = scale * scale * scale;
        if normalized_volume <= 0.0
            || !volume_scale.is_normal()
            || !(normalized_volume * volume_scale).is_normal()
        {
            return Err(SurfaceError::InvalidVolume);
        }
        Ok(Self {
            vertices,
            normalized_vertices,
            facets,
            patch_names,
            bounds,
            normalized_volume,
        })
    }

    pub fn vertices(&self) -> &[Vec3] {
        &self.vertices
    }
    pub fn normalized_vertices(&self) -> &[Vec3] {
        &self.normalized_vertices
    }
    pub fn facets(&self) -> &[Facet] {
        &self.facets
    }
    pub fn patch_names(&self) -> &[String] {
        &self.patch_names
    }
    pub fn bounds(&self) -> AxisAlignedBox {
        self.bounds
    }
    pub fn normalized_volume(&self) -> f64 {
        self.normalized_volume
    }

    pub fn volume(&self) -> f64 {
        let scale = self.bounds.extent().max_abs_component();
        self.normalized_volume * (scale * scale * scale)
    }
}

/// A closed surface whose facets support a convex volume.
#[derive(Clone, Debug, PartialEq)]
pub struct ConvexSurface(ClosedSurface);

impl ConvexSurface {
    pub fn new(
        vertices: Vec<Vec3>,
        triangles: Vec<SurfaceTriangle>,
        patch_names: Vec<String>,
    ) -> Result<Self, SurfaceError> {
        ClosedSurface::new(vertices, triangles, patch_names)?.try_into()
    }

    pub fn contains_interior(&self, point: Vec3) -> bool {
        let point = (point - self.bounds.min()) / self.bounds.extent().max_abs_component();
        self.facets
            .iter()
            .all(|facet| facet.plane.signed_distance(point) < 0.0)
    }
}

impl Deref for ConvexSurface {
    type Target = ClosedSurface;

    fn deref(&self) -> &ClosedSurface {
        &self.0
    }
}

impl TryFrom<ClosedSurface> for ConvexSurface {
    type Error = SurfaceError;

    fn try_from(surface: ClosedSurface) -> Result<Self, SurfaceError> {
        let points = surface.normalized_vertices();
        let facets = surface.facets();
        let origin = vertex_mean(points);
        for (index, facet) in facets.iter().enumerate() {
            if facet.plane.signed_distance(origin) >= 0.0
                || points
                    .iter()
                    .any(|&v| facet.plane.signed_distance(v) > RELATIVE_SURFACE_TOLERANCE)
            {
                return Err(SurfaceError::NonConvexFacet(index));
            }
        }
        for (index, facet) in facets.iter().enumerate() {
            let a = facet.triangle.vertices.map(|id| points[id]);
            for (other_index, other) in facets[..index].iter().enumerate() {
                let b = other.triangle.vertices.map(|id| points[id]);
                if facet.plane.normal().dot(other.plane.normal()) > 0.0
                    && a.iter().all(|&p| {
                        other.plane.signed_distance(p).abs() <= RELATIVE_SURFACE_TOLERANCE
                    })
                    && b.iter().all(|&p| {
                        facet.plane.signed_distance(p).abs() <= RELATIVE_SURFACE_TOLERANCE
                    })
                    && coplanar_overlap(a, b, facet.plane.normal())
                {
                    return Err(SurfaceError::OverlappingFacets {
                        first: other_index,
                        second: index,
                    });
                }
            }
        }
        Ok(Self(surface))
    }
}

fn vertex_mean(points: &[Vec3]) -> Vec3 {
    points.iter().fold(Vec3::ZERO, |sum, &v| sum + v) / points.len() as f64
}

fn coplanar_overlap(a: [Vec3; 3], b: [Vec3; 3], normal: Vec3) -> bool {
    for vertices in [a, b] {
        for i in 0..3 {
            let edge = vertices[(i + 1) % 3] - vertices[i];
            let axis = normal.cross(edge / edge.length());
            let ([amin, amax], [bmin, bmax]) = (projection(a, axis), projection(b, axis));
            if amax.min(bmax) - amin.max(bmin) <= RELATIVE_SURFACE_TOLERANCE {
                return false;
            }
        }
    }
    true
}

fn projection(triangle: [Vec3; 3], axis: Vec3) -> [f64; 2] {
    let [a, b, c] = triangle.map(|v| v.dot(axis));
    [a.min(b).min(c), a.max(b).max(c)]
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SurfaceError {
    Empty,
    InvalidBounds,
    InvalidVertex(usize),
    InvalidPatch(usize),
    DegenerateFacet(usize),
    DuplicateFacet(usize),
    OpenEdge {
        first: usize,
        second: usize,
    },
    WindingOrNonManifold {
        first: usize,
        second: usize,
    },
    Disconnected,
    InvalidVolume,
    NonConvexFacet(usize),
    OverlappingFacets {
        first: usize,
        second: usize,
    },
    InvalidStl {
        line: Option<usize>,
        reason: &'static str,
    },
    PatchCount {
        facets: usize,
        names: usize,
    },
    FacetLimit {
        limit: usize,
    },
}

impl fmt::Display for SurfaceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("surface needs at least four facets"),
            Self::InvalidBounds => {
                f.write_str("surface needs finite, positive extents in three dimensions")
            }
            Self::InvalidVertex(id) => write!(f, "surface vertex {id} is invalid or unused"),
            Self::InvalidPatch(id) => {
                write!(f, "surface patch {id} has an invalid name or assignment")
            }
            Self::DegenerateFacet(id) => write!(f, "surface facet {id} is degenerate"),
            Self::DuplicateFacet(id) => write!(f, "surface facet {id} is duplicated"),
            Self::OpenEdge { first, second } => {
                write!(f, "surface has an open edge ({first}, {second})")
            }
            Self::WindingOrNonManifold { first, second } => write!(
                f,
                "surface edge ({first}, {second}) has inconsistent winding or non-manifold incidence"
            ),
            Self::Disconnected => f.write_str("surface must be one connected shell"),
            Self::InvalidVolume => f.write_str(
                "surface must enclose a positive representable volume with outward winding",
            ),
            Self::NonConvexFacet(id) => {
                write!(f, "surface facet {id} is not a convex supporting plane")
            }
            Self::OverlappingFacets { first, second } => {
                write!(f, "surface facets {first} and {second} overlap")
            }
            Self::InvalidStl {
                line: Some(line),
                reason,
            } => write!(f, "STL line {line}: {reason}"),
            Self::InvalidStl { line: None, reason } => write!(f, "STL: {reason}"),
            Self::PatchCount { facets, names } => {
                write!(f, "STL has {facets} facets and {names} patch assignments")
            }
            Self::FacetLimit { limit } => {
                write!(f, "STL exceeds the {limit}-facet limit for this operation")
            }
        }
    }
}

impl std::error::Error for SurfaceError {}
