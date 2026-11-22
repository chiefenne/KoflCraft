//! Exhaustive CPU reference construction for box-clipped Voronoi cells.

use std::fmt;

use crate::clipping::{ClipError, ConvexPolyhedron, FaceSource};
use mesh_core::geometry::{AxisAlignedBox, Plane, Vec3};
use mesh_geometry::ConvexSurface;

/// Signed-distance tolerance expressed as a fraction of the box's longest side.
pub const RELATIVE_CLIP_TOLERANCE: f64 = 64.0 * f64::EPSILON;

/// Generate one candidate polyhedron per site, in input order.
///
/// Each cell is clipped against every other site's perpendicular bisector, in input
/// order. Output vertices are in the input coordinate system. Indices are cell-local,
/// and bisector face tags refer to the original site indices.
///
/// Sites must be strictly interior and separated by more than twice the clipping
/// tolerance. Every box dimension must also exceed that threshold.
pub fn clipped_voronoi(
    domain: AxisAlignedBox,
    sites: &[Vec3],
) -> Result<Vec<ConvexPolyhedron>, VoronoiError> {
    let scale = domain.extent().max_abs_component();
    let extent = domain.extent() / scale;
    let normalized_sites = prepare_sites(domain, sites, |site| domain.contains_interior(site))?;
    let normalized_box = AxisAlignedBox::new(Vec3::ZERO, extent).ok_or(VoronoiError::ThinDomain)?;
    clipped_cells(
        domain,
        ConvexPolyhedron::from_box(normalized_box),
        &normalized_sites,
    )
}

/// Clip a convex triangulated shell by the same site bisectors as the analytic box.
pub fn clipped_surface_voronoi(
    surface: &ConvexSurface,
    sites: &[Vec3],
) -> Result<Vec<ConvexPolyhedron>, VoronoiError> {
    let bounds = surface.bounds();
    let sites = prepare_sites(bounds, sites, |site| surface.contains_interior(site))?;
    clipped_cells(bounds, ConvexPolyhedron::from_surface(surface), &sites)
}

fn clipped_cells(
    domain: AxisAlignedBox,
    initial: ConvexPolyhedron,
    normalized_sites: &[Vec3],
) -> Result<Vec<ConvexPolyhedron>, VoronoiError> {
    let origin = domain.min();
    let maximum = domain.max();
    let scale = domain.extent().max_abs_component();
    let extent = domain.extent() / scale;
    let mut cells = Vec::with_capacity(normalized_sites.len());
    for (index, &site) in normalized_sites.iter().enumerate() {
        let mut cell = initial.clone();
        for (other_index, &other_site) in normalized_sites.iter().enumerate() {
            if index == other_index {
                continue;
            }
            // Equidistant points lie on the plane through the midpoint, normal to
            // the segment joining the sites. The current site has negative distance.
            let plane = Plane::new((site + other_site) * 0.5, other_site - site).ok_or(
                VoronoiError::SitesTooClose {
                    first: index,
                    second: other_index,
                },
            )?;
            cell = cell
                .clip(
                    plane,
                    FaceSource::Bisector {
                        other_site: other_index,
                    },
                    RELATIVE_CLIP_TOLERANCE,
                )
                .map_err(|source| VoronoiError::Clip {
                    site: index,
                    other_site: other_index,
                    source,
                })?
                .ok_or(VoronoiError::EmptyCell { site: index })?;
        }
        cell.map_vertices(|point| {
            // Keep the six domain planes exact through the round trip.
            let restore = |t: f64, length: f64, min: f64, max: f64| {
                if t >= length {
                    max
                } else {
                    t.mul_add(scale, min).clamp(min, max)
                }
            };
            Vec3::new(
                restore(point.x, extent.x, origin.x, maximum.x),
                restore(point.y, extent.y, origin.y, maximum.y),
                restore(point.z, extent.z, origin.z, maximum.z),
            )
        });
        let vertices = cell.vertices();
        if vertices
            .iter()
            .enumerate()
            .any(|(i, vertex)| vertices[..i].contains(vertex))
        {
            return Err(VoronoiError::CollapsedVertices { site: index });
        }
        cells.push(cell);
    }
    Ok(cells)
}

/// Apply the input checks used by the clipping backend.
pub fn validate_sites(domain: AxisAlignedBox, sites: &[Vec3]) -> Result<(), VoronoiError> {
    prepare_sites(domain, sites, |site| domain.contains_interior(site)).map(|_| ())
}

pub fn validate_surface_sites(surface: &ConvexSurface, sites: &[Vec3]) -> Result<(), VoronoiError> {
    prepare_sites(surface.bounds(), sites, |site| {
        surface.contains_interior(site)
    })
    .map(|_| ())
}

fn prepare_sites(
    domain: AxisAlignedBox,
    sites: &[Vec3],
    contains: impl Fn(Vec3) -> bool,
) -> Result<Vec<Vec3>, VoronoiError> {
    if sites.is_empty() {
        return Err(VoronoiError::NoSites);
    }
    let origin = domain.min();
    let scale = domain.extent().max_abs_component();
    let extent = domain.extent() / scale;
    let resolution = 2.0 * RELATIVE_CLIP_TOLERANCE;
    if extent.x.min(extent.y).min(extent.z) <= resolution {
        return Err(VoronoiError::ThinDomain);
    }

    let mut normalized_sites = Vec::with_capacity(sites.len());
    for (index, &site) in sites.iter().enumerate() {
        if !contains(site) {
            return Err(VoronoiError::SiteOutside { index });
        }
        let point: Vec3 = (site - origin) / scale;
        if let Some(first) = normalized_sites
            .iter()
            .position(|&other| (point - other).length() <= resolution)
        {
            return Err(VoronoiError::SitesTooClose {
                first,
                second: index,
            });
        }
        normalized_sites.push(point);
    }
    Ok(normalized_sites)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum VoronoiError {
    NoSites,
    SiteOutside {
        index: usize,
    },
    SitesTooClose {
        first: usize,
        second: usize,
    },
    ThinDomain,
    Clip {
        site: usize,
        other_site: usize,
        source: ClipError,
    },
    EmptyCell {
        site: usize,
    },
    CollapsedVertices {
        site: usize,
    },
}

impl fmt::Display for VoronoiError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoSites => f.write_str("at least one site is required"),
            Self::SiteOutside { index } => {
                write!(f, "site {index} must be strictly inside the domain")
            }
            Self::SitesTooClose { first, second } => {
                write!(
                    f,
                    "sites {first} and {second} are closer than the clipping resolution"
                )
            }
            Self::ThinDomain => f.write_str("a box dimension is below the clipping resolution"),
            Self::Clip {
                site,
                other_site,
                source,
            } => {
                write!(
                    f,
                    "clipping site {site} against site {other_site} failed: {source}"
                )
            }
            Self::EmptyCell { site } => write!(f, "clipping removed the volume of cell {site}"),
            Self::CollapsedVertices { site } => {
                write!(
                    f,
                    "cell {site} has coincident vertices in the input coordinate system"
                )
            }
        }
    }
}

impl std::error::Error for VoronoiError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Clip { source, .. } => Some(source),
            _ => None,
        }
    }
}
