//! Clipped Voronoi cells in a box or a convex surface. Each cell is clipped
//! by the bisectors of its nearest sites until the security radius proves
//! it complete (ALG-04, ALG-13); the exhaustive all-pairs construction of
//! M1 stays as the reference ([`Strategy::AllPairs`]).

use std::fmt;

use crate::clipping::{ClipError, ConvexPolyhedron, FaceSource};
use crate::neighbours::KdTree;
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
    clipped_voronoi_with(domain, sites, Strategy::Nearest)
}

/// [`clipped_voronoi`] with a chosen clipping strategy.
pub fn clipped_voronoi_with(
    domain: AxisAlignedBox,
    sites: &[Vec3],
    strategy: Strategy,
) -> Result<Vec<ConvexPolyhedron>, VoronoiError> {
    let scale = domain.extent().max_abs_component();
    let extent = domain.extent() / scale;
    let normalized_sites = prepare_sites(domain, sites, |site| domain.contains_interior(site))?;
    let normalized_box = AxisAlignedBox::new(Vec3::ZERO, extent).ok_or(VoronoiError::ThinDomain)?;
    clipped_cells(
        domain,
        ConvexPolyhedron::from_box(normalized_box),
        &normalized_sites,
        strategy,
        &|source| match source {
            FaceSource::Boundary(side) => {
                let a = side.axis();
                let mut n = [0.0; 3];
                n[a] = 1.0;
                let c = if side.is_max() {
                    [extent.x, extent.y, extent.z][a]
                } else {
                    0.0
                };
                Some((Vec3::new(n[0], n[1], n[2]), c))
            }
            _ => None,
        },
    )
}

/// Clip a convex triangulated shell by the same site bisectors as the analytic box.
pub fn clipped_surface_voronoi(
    surface: &ConvexSurface,
    sites: &[Vec3],
) -> Result<Vec<ConvexPolyhedron>, VoronoiError> {
    clipped_surface_voronoi_with(surface, sites, Strategy::Nearest)
}

/// [`clipped_surface_voronoi`] with a chosen clipping strategy.
pub fn clipped_surface_voronoi_with(
    surface: &ConvexSurface,
    sites: &[Vec3],
    strategy: Strategy,
) -> Result<Vec<ConvexPolyhedron>, VoronoiError> {
    let bounds = surface.bounds();
    let sites = prepare_sites(bounds, sites, |site| surface.contains_interior(site))?;
    clipped_cells(
        bounds,
        ConvexPolyhedron::from_surface(surface),
        &sites,
        strategy,
        &|source| match source {
            FaceSource::Surface { facet, .. } => {
                let plane = surface.facets().get(facet.0)?.plane;
                Some((plane.normal(), -plane.signed_distance(Vec3::ZERO)))
            }
            _ => None,
        },
    )
}

/// Which bisectors clip a cell.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Strategy {
    /// The nearest sites in order of distance, until the security radius
    /// proves the cell complete (ALG-04, ALG-13).
    #[default]
    Nearest,
    /// Every other site in input order: the exhaustive reference of M1.
    AllPairs,
}

/// Sites queried first per cell; doubled until the cell is proven complete.
const INITIAL_NEIGHBOURS: usize = 32;
/// Relative slack on the security test, against rounding in the radius.
const SECURITY_MARGIN: f64 = 1.0 + 1e-9;

fn clipped_cells(
    domain: AxisAlignedBox,
    initial: ConvexPolyhedron,
    normalized_sites: &[Vec3],
    strategy: Strategy,
    boundary: &dyn Fn(FaceSource) -> Option<(Vec3, f64)>,
) -> Result<Vec<ConvexPolyhedron>, VoronoiError> {
    let origin = domain.min();
    let maximum = domain.max();
    let scale = domain.extent().max_abs_component();
    let extent = domain.extent() / scale;
    let mut cells = Vec::with_capacity(normalized_sites.len());
    let tree = (strategy == Strategy::Nearest).then(|| KdTree::new(normalized_sites));
    for (index, &site) in normalized_sites.iter().enumerate() {
        let mut cell = initial.clone();
        let clip = |cell: ConvexPolyhedron, other_index: usize| {
            let other_site = normalized_sites[other_index];
            // Equidistant points lie on the plane through the midpoint, normal to
            // the segment joining the sites. The current site has negative distance.
            let plane = Plane::new((site + other_site) * 0.5, other_site - site).ok_or(
                VoronoiError::SitesTooClose {
                    first: index,
                    second: other_index,
                },
            )?;
            cell.clip(
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
            .ok_or(VoronoiError::EmptyCell { site: index })
        };
        match &tree {
            None => {
                for other_index in 0..normalized_sites.len() {
                    if other_index != index {
                        cell = clip(cell, other_index)?;
                    }
                }
            }
            Some(tree) => {
                // Security radius (Ray, Sokolov, Lefebvre & Lévy 2018): once the
                // next nearest site is farther than twice the cell's radius
                // around its site (plus the clipping tolerance), no farther
                // bisector can touch a vertex, so the cell is complete.
                let others = normalized_sites.len() - 1;
                let mut k = INITIAL_NEIGHBOURS.min(others);
                let mut done = 0;
                'grow: loop {
                    let near = tree.nearest(index, k);
                    for &(squared, other_index) in &near[done..] {
                        let radius = cell
                            .vertices()
                            .iter()
                            .map(|&v| (v - site).length())
                            .fold(0.0, f64::max);
                        let reach = 2.0 * (radius + RELATIVE_CLIP_TOLERANCE);
                        if squared.sqrt() > reach * SECURITY_MARGIN {
                            break 'grow;
                        }
                        cell = clip(cell, other_index)?;
                        done += 1;
                    }
                    if k == others {
                        break;
                    }
                    k = (2 * k).min(others);
                }
            }
        }
        snap_vertices(&mut cell, index, normalized_sites, boundary);
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

/// Recompute every simple vertex of a cell (exactly three incident faces)
/// from the sites and boundary planes that define it, with the sites in
/// index order, so every cell sharing the vertex computes the same point to
/// the last bit, whatever order its planes were cut in. A vertex where
/// four or more faces meet (degenerate, as on a lattice), or whose system is
/// singular or moves the point more than the clipping tolerance, keeps its
/// clipped position.
fn snap_vertices(
    cell: &mut ConvexPolyhedron,
    site: usize,
    sites: &[Vec3],
    boundary: &dyn Fn(FaceSource) -> Option<(Vec3, f64)>,
) {
    let mut incident: Vec<Vec<FaceSource>> = vec![Vec::new(); cell.vertices().len()];
    for face in cell.faces() {
        for &v in face.vertices() {
            incident[v].push(face.source());
        }
    }
    let snapped: Vec<Option<Vec3>> = incident
        .iter()
        .zip(cell.vertices())
        .map(|(sources, &old)| {
            if sources.len() != 3 {
                return None;
            }
            let mut ids = vec![site];
            let mut rows: Vec<(Vec3, f64)> = Vec::with_capacity(3);
            for &source in sources {
                match source {
                    FaceSource::Bisector { other_site } => ids.push(other_site),
                    other => rows.push(boundary(other)?),
                }
            }
            ids.sort_unstable();
            let s0 = sites[ids[0]];
            let mut bisectors: Vec<(Vec3, f64)> = ids[1..]
                .iter()
                .map(|&m| {
                    let d = sites[m] - s0;
                    (d, d.dot(sites[m] + s0) * 0.5)
                })
                .collect();
            // boundary rows in a fixed order, after the bisectors
            rows.sort_by(|a, b| {
                [a.0.x, a.0.y, a.0.z, a.1]
                    .iter()
                    .zip([b.0.x, b.0.y, b.0.z, b.1].iter())
                    .map(|(x, y)| x.total_cmp(y))
                    .find(|o| o.is_ne())
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            bisectors.extend(rows);
            let point = solve3(&bisectors)?;
            ((point - old).length() <= RELATIVE_CLIP_TOLERANCE * SNAP_REACH).then_some(point)
        })
        .collect();
    cell.map_indexed_vertices(|i, v| snapped[i].unwrap_or(v));
}

/// How far, in clipping tolerances, a recomputed vertex may move.
const SNAP_REACH: f64 = 1024.0;

/// Solve the 3×3 system `n_i · x = c_i` by Gaussian elimination with
/// partial pivoting; `None` when it is singular.
fn solve3(rows: &[(Vec3, f64)]) -> Option<Vec3> {
    if rows.len() != 3 {
        return None;
    }
    let mut m: Vec<[f64; 4]> = rows.iter().map(|(n, c)| [n.x, n.y, n.z, *c]).collect();
    let scale = m
        .iter()
        .flat_map(|r| r[..3].iter())
        .fold(0.0f64, |a, &b| a.max(b.abs()));
    if scale == 0.0 || !scale.is_finite() {
        return None;
    }
    for col in 0..3 {
        let pivot = (col..3).max_by(|&a, &b| m[a][col].abs().total_cmp(&m[b][col].abs()))?;
        if m[pivot][col].abs() <= 1e-12 * scale {
            return None;
        }
        m.swap(col, pivot);
        let pivot_row = m[col];
        for row in m.iter_mut().skip(col + 1) {
            let f = row[col] / pivot_row[col];
            for (x, p) in row.iter_mut().zip(pivot_row).skip(col) {
                *x -= f * p;
            }
        }
    }
    let z = m[2][3] / m[2][2];
    let y = (m[1][3] - m[1][2] * z) / m[1][1];
    let x = (m[0][3] - m[0][1] * y - m[0][2] * z) / m[0][0];
    let p = Vec3::new(x, y, z);
    p.is_finite().then_some(p)
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
        normalized_sites.push((site - origin) / scale);
    }
    // as the exhaustive check: the first site with an earlier site within
    // the resolution, paired with the earliest such site
    let tree = KdTree::new(&normalized_sites);
    let others = normalized_sites.len() - 1;
    let mut worst = None;
    for second in 0..normalized_sites.len() {
        let mut k = 1.min(others);
        let first = loop {
            let near = tree.nearest(second, k);
            if k == others || near.last().is_none_or(|&(d, _)| d.sqrt() > resolution) {
                break near
                    .iter()
                    .filter(|&&(d, i)| d.sqrt() <= resolution && i < second)
                    .map(|&(_, i)| i)
                    .min();
            }
            k = (2 * k).min(others);
        };
        if let Some(first) = first {
            worst = Some((first, second));
            break;
        }
    }
    if let Some((first, second)) = worst {
        return Err(VoronoiError::SitesTooClose { first, second });
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
