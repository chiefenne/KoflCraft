//! Boundary-fitted tetrahedra on an explicitly partitioned shell, with optional
//! isoparametric wedge layers. Region interfaces follow existing STL edges.
mod columns;
mod convex;
mod front;
mod hex_join;
mod intersection;
mod polyhedral;
mod symbolic;
mod verify;
mod volume_dual;
mod voronoi_core;
mod wall_seeds;
mod yield_layers;
pub use volume_dual::{Scaffold, build_volume_dual, build_volume_dual_hex};
pub use voronoi_core::{
    DEFAULT_LLOYD_ITERATIONS, VoronoiOptions, VoronoiStats, WallJoin, build_voronoi_core,
};

use mesh_core::geometry::Vec3;
use mesh_geometry::ClosedSurface;
use std::collections::{BTreeMap, BTreeSet};

pub use polyhedral::Polyhedron;
pub use verify::Measurements;
pub use yield_layers::{YieldPlane, extrude as extrude_yielding, shorten as shorten_yielding};

pub struct Region {
    pub kernel: Vec3,
    pub facets: Vec<[usize; 2]>,
}

pub enum Element {
    Tet([usize; 4]),
    Wedge([usize; 6]),
    Polyhedron(Box<Polyhedron>),
}

impl Element {
    pub fn nodes(&self) -> &[usize] {
        match self {
            Self::Tet(n) => n,
            Self::Wedge(n) => n,
            Self::Polyhedron(p) => &p.nodes,
        }
    }

    pub(crate) fn faces(&self) -> Vec<Vec<usize>> {
        match *self {
            Self::Polyhedron(ref p) => p.faces.clone(),
            Self::Tet([a, b, c, d]) => {
                vec![vec![a, c, b], vec![a, b, d], vec![b, c, d], vec![c, a, d]]
            }
            Self::Wedge([a, b, c, d, e, f]) => vec![
                vec![a, c, b],
                vec![d, e, f],
                vec![a, b, e, d],
                vec![b, c, f, e],
                vec![c, a, d, f],
            ],
        }
    }
}

pub struct Cell {
    pub element: Element,
    pub layer: usize,
    pub region: usize,
}

pub struct BoundaryFace {
    pub nodes: Vec<usize>,
    /// 1 inlet, 2 outlet (port caps), 3 and above: wall layers (see
    /// [`is_wall`]).
    pub patch: usize,
    pub owner: usize,
}
pub struct DuctMesh {
    pub points: Vec<Vec3>,
    pub cells: Vec<Cell>,
    pub boundary: Vec<BoundaryFace>,
    pub measurements: Measurements,
    /// Name and role (`inlet`, `outlet` or `wall`) of every patch id,
    /// indexed by id - 1; the builders set [`classic_patches`].
    pub patches: Vec<(String, String)>,
}

/// The patches of a duct with an inlet cap, an outlet cap and one wall.
pub fn classic_patches() -> Vec<(String, String)> {
    [("inlet", "inlet"), ("outlet", "outlet"), ("walls", "wall")]
        .map(|(n, r)| (n.to_string(), r.to_string()))
        .to_vec()
}
type Edge = (usize, usize);

/// The patch id of the first wall patch. Ids 1 and 2 are the inlet and
/// outlet caps, which contract in their planes; every id from 3 up is a
/// patch the wall layers grow from (walls, and openings flush with a wall),
/// each a named boundary of its own.
pub const WALL: usize = 3;

/// Added to the id of a wall patch that grows no layers (`no_layers` in
/// the partition). Such a patch must be planar: where layered walls meet
/// it, their layers run straight into it, sliding along its plane, and its
/// own facets move in its plane to make room (as an inlet cap contracts).
pub const SLIP: usize = 1 << 16;

/// Whether a patch id grows wall layers (3 and above, below [`SLIP`]).
pub fn is_wall(patch: usize) -> bool {
    (WALL..SLIP).contains(&patch)
}

/// Whether a patch id is a wall patch without layers (see [`SLIP`]).
pub fn is_slip(patch: usize) -> bool {
    (SLIP..HEX).contains(&patch)
}

/// Scaffold boundary tag of the shell of a hex core ([`crate::hexcore`]):
/// an inner boundary of the scaffold that the cubes replace, no patch.
pub const HEX: usize = crate::hexcore::HEX;

/// Whether a scaffold boundary tag is the shell of a hex core.
pub fn is_hex(patch: usize) -> bool {
    patch == HEX
}

/// Whether a patch id is an inlet or outlet cap (1 or 2).
pub fn is_cap(patch: usize) -> bool {
    patch < WALL
}

/// The patch id without the [`SLIP`] mark: the index of its name.
pub fn patch_index(patch: usize) -> usize {
    if is_slip(patch) { patch - SLIP } else { patch }
}

fn edges([a, b, c]: [usize; 3]) -> [Edge; 3] {
    [(a, b), (b, c), (c, a)]
}

fn region_labels(count: usize, regions: &[Region]) -> Result<Vec<usize>, String> {
    let mut labels = vec![usize::MAX; count];
    for (id, r) in regions.iter().enumerate() {
        if !r.kernel.is_finite() || r.facets.is_empty() {
            return Err("invalid region or kernel".into());
        }
        for &[a, b] in &r.facets {
            let range = labels
                .get_mut(a..b)
                .filter(|v| !v.is_empty())
                .ok_or("invalid region facet range")?;
            for f in range {
                if *f != usize::MAX {
                    return Err("overlapping region facet ranges".into());
                }
                *f = id;
            }
        }
    }
    if labels.contains(&usize::MAX) {
        return Err("unassigned region facets".into());
    }
    Ok(labels)
}

/// Require a single oriented loop. This backend uses one disk per pair of adjacent regions.
fn loop_vertices(edges: &[Edge]) -> Result<Vec<usize>, String> {
    let next: BTreeMap<_, _> = edges.iter().copied().collect();
    if edges.len() < 3 || next.len() != edges.len() {
        return Err("branched or empty interface loop".into());
    }
    let mut ids = Vec::new();
    let mut seen = BTreeSet::new();
    let mut at = edges[0].0;
    while seen.insert(at) {
        ids.push(at);
        at = *next.get(&at).ok_or("open interface loop")?;
    }
    if at != ids[0] || ids.len() != edges.len() {
        return Err("interface requires one simple loop".into());
    }
    Ok(ids)
}

pub(super) fn mean(points: &[Vec3], ids: &[usize]) -> Vec3 {
    let origin = points[ids[0]];
    origin
        + ids
            .iter()
            .fold(Vec3::ZERO, |s, &i| s + (points[i] - origin))
            / ids.len() as f64
}

/// Ports are represented by their actual STL triangles and explicit patch IDs.
/// Thicknesses are individual increments, with an empty list selecting the bulk mesh.
pub fn build(
    surface: &ClosedSurface,
    patches: &[usize],
    regions: &[Region],
    thicknesses: &[f64],
) -> Result<DuctMesh, String> {
    let facets = surface.facets();
    if patches.len() != facets.len() || patches.contains(&0) {
        return Err("expected a patch on every STL facet".into());
    }
    let labels = region_labels(facets.len(), regions)?;
    let scale = surface.bounds().extent().max_abs_component();
    let outer: Vec<_> = facets.iter().map(|f| f.triangle.vertices).collect();
    intersection::embedded(surface.vertices(), &outer, scale)
        .map_err(|e| format!("input STL: {e}"))?;
    let mut pending = BTreeMap::new();
    let mut seams = BTreeMap::<(usize, usize), Vec<Edge>>::new();
    let mut rims = BTreeMap::<usize, Vec<Edge>>::new();
    for (i, f) in facets.iter().enumerate() {
        for (a, b) in edges(f.triangle.vertices) {
            if let Some(j) = pending.remove(&(b, a)) {
                let (left, right) = (labels[i], labels[j]);
                if left != right {
                    let (pair, edge) = if left < right {
                        ((left, right), (a, b))
                    } else {
                        ((right, left), (b, a))
                    };
                    seams.entry(pair).or_default().push(edge);
                }
                if patches[i] != patches[j]
                    && !(is_wall(patches[i]) && is_wall(patches[j]))
                    && !(is_slip(patches[i]) && is_slip(patches[j]))
                {
                    let (port, edge) = if is_wall(patches[i]) {
                        (patches[j], (a, b))
                    } else if is_wall(patches[j]) {
                        (patches[i], (b, a))
                    } else {
                        return Err(
                            "two patches without layers (inlet, outlet or no_layers) share an edge"
                                .into(),
                        );
                    };
                    rims.entry(port).or_default().push(edge); // wall winding
                }
            } else {
                pending.insert((a, b), i);
            }
        }
    }
    if !pending.is_empty() {
        return Err("unpaired surface edges".into());
    }
    let mut front = front::advance(surface, patches, &labels, &rims, thicknesses, None)?;
    let tol = surface.bounds().extent().max_abs_component() * 1e-10;
    let mut kernels = Vec::new();
    for r in regions {
        kernels.push(front.points.len());
        front.points.push(r.kernel);
    }
    for (i, f) in facets.iter().enumerate() {
        let [a, b, c] = f.triangle.vertices.map(|v| front.inner[v]);
        let owner = front.cells.len();
        front.cells.push(Cell {
            element: Element::Tet([kernels[labels[i]], a, b, c]),
            layer: 0,
            region: labels[i],
        });
        if !is_wall(patches[i]) || thicknesses.is_empty() {
            front.boundary.push(BoundaryFace {
                nodes: vec![a, b, c],
                patch: patches[i],
                owner,
            });
        }
    }
    for ((left, right), edge_loop) in &seams {
        let ids = loop_vertices(edge_loop)?;
        let center = mean(surface.vertices(), &ids);
        // The input disk must be planar and its fan consistently oriented.
        let delta = regions[*right].kernel - regions[*left].kernel;
        let mut first_normal = None;
        for &(a, b) in edge_loop {
            let v = (surface.vertices()[b] - center).cross(surface.vertices()[a] - center);
            if v.dot(delta) <= 0.0 {
                return Err(format!("invalid fan at region interface {left}/{right}"));
            }
            let normal = *first_normal.get_or_insert_with(|| v / v.length());
            if (surface.vertices()[a] - center).dot(normal).abs() > tol || v.dot(normal) <= 0.0 {
                return Err("source interface is nonplanar or has an inverted fan".into());
            }
        }
        let inner_ids: Vec<_> = ids.iter().map(|&v| front.inner[v]).collect();
        let center_id = front.points.len();
        front.points.push(mean(&front.points, &inner_ids));
        for &(a, b) in edge_loop {
            let (a, b) = (front.inner[a], front.inner[b]);
            for (region, ring) in [(*left, [b, a, center_id]), (*right, [a, b, center_id])] {
                front.cells.push(Cell {
                    element: Element::Tet([kernels[region], ring[0], ring[1], ring[2]]),
                    layer: 0,
                    region,
                });
            }
        }
    }
    let measurements = verify::check(
        surface,
        patches,
        &front.points,
        &front.cells,
        &front.boundary,
        thicknesses.len(),
    )?;
    Ok(DuctMesh {
        points: front.points,
        cells: front.cells,
        boundary: front.boundary,
        measurements,
        patches: classic_patches(),
    })
}

/// Build a boundary-fitted polyhedral core with the supplied wall layers.
pub fn build_polyhedral(
    surface: &ClosedSurface,
    patches: &[usize],
    regions: &[Region],
    thicknesses: &[f64],
) -> Result<DuctMesh, String> {
    let mut mesh = build(surface, patches, regions, thicknesses)?;
    let clearance = mesh.measurements.minimum_kernel_clearance;
    polyhedral::convert(&mut mesh)?;
    mesh.measurements = verify::check(
        surface,
        patches,
        &mesh.points,
        &mesh.cells,
        &mesh.boundary,
        thicknesses.len(),
    )?;
    mesh.measurements.minimum_kernel_clearance = clearance;
    Ok(mesh)
}
