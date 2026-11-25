//! Walls whose layers yield (`layers_yield` in the partition): a planar
//! wall keeps its layers, and every wall that crosses it runs its own layers
//! straight through them down to it, with no mitre at the junction.
//!
//! The duct is first built on a shortened domain: each yielding plane moves
//! into the fluid by the total layer height, and it is built as a wall
//! without layers ([`super::SLIP`]), so the crossing walls' layers end on the
//! moved plane. [`extrude`] then fills the slab between the moved plane and
//! the real wall by extruding every face on the moved plane along the
//! plane's normal through the layer heights: under core cells these are the
//! yielding wall's own layers, under the crossing walls' layer columns they
//! are blocks of (crossing layer × yielding layer) cells, conforming on all
//! sides. The crossing walls must meet the plane at a right angle, so their
//! extrusion through the slab is their own surface.

use super::*;
use mesh_core::geometry::PatchId;
use mesh_geometry::SurfaceTriangle;
use std::collections::HashMap;

/// A yielding plane: its outward unit normal and its offset on the real wall
/// (`normal · x = offset`).
#[derive(Clone, Copy, Debug)]
pub struct YieldPlane {
    pub normal: Vec3,
    pub offset: f64,
}

/// Largest |cos| between a crossing wall facet's normal and a yielding
/// plane's normal (about 6°): the facets of a triangulated tube that end on
/// the plane tilt by a few degrees (2.9° on the furnace tubes), and the
/// slab replaces their last layer height by a straight extrusion.
const RIGHT_ANGLE_SLACK: f64 = 0.1;

/// Width of the band in front of a yielding plane that moves with it, in
/// total layer heights.
const RAMP: f64 = 4.0;

/// Move the facets of the yielding patches (`yielding[f]`) into the fluid
/// by `total` along their plane's normal. Each edge-connected piece must be
/// planar, pieces must not touch each other, and every other facet at a
/// moved vertex must stand at a right angle on the plane and keep its
/// orientation. Returns the shortened surface and the planes.
pub fn shorten(
    surface: &ClosedSurface,
    yielding: &[bool],
    total: f64,
) -> Result<(ClosedSurface, Vec<YieldPlane>), String> {
    let original = surface.vertices();
    let facets = surface.facets();
    let scale = surface.bounds().extent().max_abs_component();
    let tol = scale * front::PORT_RELATIVE_TOLERANCE;
    // edge-connected pieces of the yielding facets
    let mut piece = vec![usize::MAX; facets.len()];
    let mut by_edge = HashMap::<Edge, usize>::new();
    for (f, facet) in facets.iter().enumerate() {
        for (a, b) in edges(facet.triangle.vertices) {
            by_edge.insert((a, b), f);
        }
    }
    let mut planes = Vec::new();
    let mut vertex_plane = vec![usize::MAX; original.len()];
    for start in 0..facets.len() {
        if !yielding[start] || piece[start] != usize::MAX {
            continue;
        }
        let id = planes.len();
        let normal = facets[start].plane.normal();
        let origin = original[facets[start].triangle.vertices[0]];
        let offset = normal.dot(origin);
        let mut stack = vec![start];
        piece[start] = id;
        while let Some(f) = stack.pop() {
            let facet = &facets[f];
            if facet.plane.normal().dot(normal) < 1.0 - front::PORT_NORMAL_SLACK
                || facet
                    .triangle
                    .vertices
                    .iter()
                    .any(|&v| (normal.dot(original[v]) - offset).abs() > tol)
            {
                return Err(
                    "a wall whose layers yield must be planar in each connected piece".into(),
                );
            }
            for (a, b) in edges(facet.triangle.vertices) {
                let g = by_edge[&(b, a)];
                if yielding[g] && piece[g] == usize::MAX {
                    piece[g] = id;
                    stack.push(g);
                }
            }
            for &v in &facet.triangle.vertices {
                if vertex_plane[v] != usize::MAX && vertex_plane[v] != id {
                    return Err(
                        "two walls whose layers yield meet; only one may yield at a junction"
                            .into(),
                    );
                }
                vertex_plane[v] = id;
            }
        }
        planes.push(YieldPlane { normal, offset });
    }
    // The plane moves in by `total`; the vertices in a band of `RAMP`
    // times that width in front of it move with a linear ramp, so a short
    // facet near the plane keeps its orientation. Walls crossing the plane
    // stand at a right angle on it, so their vertices slide in their own
    // planes and the domain keeps its shape outside the slab.
    let width = RAMP * total;
    let mut vertices = original.to_vec();
    let mut band = vec![usize::MAX; original.len()];
    for (v, &k) in vertex_plane.iter().enumerate() {
        if k != usize::MAX {
            vertices[v] = original[v] - planes[k].normal * total;
            band[v] = k;
            continue;
        }
        for (j, q) in planes.iter().enumerate() {
            let depth = q.offset - q.normal.dot(original[v]);
            if depth > -tol && depth < width {
                if band[v] != usize::MAX {
                    return Err(
                        "two walls whose layers yield are closer than four layer heights".into(),
                    );
                }
                band[v] = j;
                vertices[v] = original[v] - q.normal * (total * (1.0 - depth.max(0.0) / width));
            }
        }
    }
    for (f, facet) in facets.iter().enumerate() {
        if yielding[f] {
            continue;
        }
        let Some(&k) = facet
            .triangle
            .vertices
            .iter()
            .map(|&v| &band[v])
            .find(|&&k| k != usize::MAX)
        else {
            continue;
        };
        if facet.plane.normal().dot(planes[k].normal).abs() > RIGHT_ANGLE_SLACK {
            return Err(
                "a wall near a wall whose layers yield must meet it at a right angle".into(),
            );
        }
        let [a, b, c] = facet.triangle.vertices.map(|v| vertices[v]);
        let turn = (b - a).cross(c - a).dot(facet.plane.normal());
        if turn.is_nan() || turn <= 0.0 {
            return Err("a facet folds where a wall whose layers yield moves in".into());
        }
    }
    let triangles = facets.iter().map(|f| f.triangle).collect();
    let shortened = ClosedSurface::new(vertices, triangles, surface.patch_names().to_vec())
        .map_err(|e| format!("shortened surface: {e:?}"))?;
    Ok((shortened, planes))
}

/// Fill the slab of every yielding plane: extrude each boundary face of a
/// patch in `patches` that lies on a moved plane through the layer `steps`
/// (the first step at the real wall). `shell` is the shortened scaffold
/// surface the mesh was built on, with per-facet patch ids `shell_patches`
/// (without the [`SLIP`] mark). The mesh is checked again on the full
/// domain.
pub fn extrude(
    mesh: &mut DuctMesh,
    shell: &ClosedSurface,
    shell_patches: &[usize],
    planes: &[YieldPlane],
    patches: &[usize],
    steps: &[f64],
) -> Result<ClosedSurface, String> {
    let total: f64 = steps.iter().sum();
    let n = steps.len();
    // offsets of the levels from the moved plane: level n is the real wall
    let mut level = vec![0.0; n + 1];
    for k in 1..=n {
        level[k] = level[k - 1] + steps[n - k];
    }
    let scale = shell.bounds().extent().max_abs_component();
    let tol = scale * 1e-6;
    let on_plane = |x: Vec3| {
        planes
            .iter()
            .position(|q| (q.normal.dot(x) - (q.offset - total)).abs() <= tol)
    };
    let mut slab = Vec::new();
    let mut keep = Vec::new();
    for face in std::mem::take(&mut mesh.boundary) {
        let plane = patches
            .contains(&face.patch)
            .then(|| on_plane(mesh.points[face.nodes[0]]))
            .flatten()
            .filter(|&k| {
                face.nodes
                    .iter()
                    .all(|&v| on_plane(mesh.points[v]) == Some(k))
            });
        match plane {
            Some(k) => slab.push((face, k)),
            None => keep.push(face),
        }
    }
    if slab.is_empty() {
        return Err("no face lies on a wall whose layers yield".into());
    }
    // boundary faces along the rim of the slab: patch of each directed edge
    let mut edge_patch = HashMap::<Edge, usize>::new();
    for face in &keep {
        let m = face.nodes.len();
        for i in 0..m {
            edge_patch.insert((face.nodes[i], face.nodes[(i + 1) % m]), face.patch);
        }
    }
    let mut copies = HashMap::<usize, Vec<usize>>::new();
    let mut img = |points: &mut Vec<Vec3>, v: usize, k: usize, normal: Vec3| -> usize {
        if k == 0 {
            return v;
        }
        let list = copies.entry(v).or_insert_with(|| {
            (1..=n)
                .map(|j| {
                    points.push(points[v] + normal * level[j]);
                    points.len() - 1
                })
                .collect()
        });
        list[k - 1]
    };
    let mut interfaces = mesh.measurements.prism_core_faces;
    for (face, k) in &slab {
        let normal = planes[*k].normal;
        let ring = &face.nodes;
        let m = ring.len();
        if mesh.cells[face.owner].layer == 0 {
            interfaces += 1;
        }
        let levels: Vec<Vec<usize>> = (0..=n)
            .map(|j| {
                ring.iter()
                    .map(|&v| img(&mut mesh.points, v, j, normal))
                    .collect()
            })
            .collect();
        let flat: Vec<Vec3> = ring.iter().map(|&v| mesh.points[v]).collect();
        let fan = columns::ears(&flat);
        for j in 1..=n {
            let id = mesh.cells.len();
            let (low, high) = (&levels[j - 1], &levels[j]);
            let mut faces = vec![low.iter().rev().copied().collect::<Vec<_>>(), high.clone()];
            for i in 0..m {
                let (u, v) = (ring[i], ring[(i + 1) % m]);
                let quad = vec![low[i], low[(i + 1) % m], high[(i + 1) % m], high[i]];
                // a rim edge: the other side is a crossing wall
                if let Some(&patch) = edge_patch.get(&(v, u)) {
                    keep.push(BoundaryFace {
                        nodes: quad.clone(),
                        patch,
                        owner: id,
                    });
                }
                faces.push(quad);
            }
            if j == n {
                keep.push(BoundaryFace {
                    nodes: high.clone(),
                    patch: face.patch,
                    owner: id,
                });
            }
            let wedges: Vec<[usize; 6]> = fan
                .iter()
                .map(|&[a, b, c]| [low[a], low[b], low[c], high[a], high[b], high[c]])
                .collect();
            let expected: f64 = wedges
                .iter()
                .map(|w| verify::wedge(w.map(|i| mesh.points[i])).1)
                .sum();
            let interior = (0..m).fold(Vec3::ZERO, |s, i| {
                s + mesh.points[low[i]] + mesh.points[high[i]]
            }) / (2 * m) as f64;
            let mut poly = Polyhedron::new(faces, interior, None, expected);
            poly.prism = Some(wedges);
            mesh.cells.push(Cell {
                element: Element::Polyhedron(Box::new(poly)),
                // layer 1 at the real wall
                layer: n + 1 - j,
                region: 0,
            });
        }
    }
    mesh.boundary = keep;
    // the full domain: the shell without its moved planes, plus the slabs
    let full = full_surface(shell, shell_patches, planes, total, tol)?;
    let full_patches: Vec<usize> = full
        .facets()
        .iter()
        .map(|f| full.patch_names()[f.triangle.patch.0].parse().unwrap_or(0))
        .collect();
    mesh.measurements = verify::check_with_interface_count(
        &full,
        &full_patches,
        &mesh.points,
        &mesh.cells,
        &mesh.boundary,
        n,
        interfaces,
    )?;
    Ok(full)
}

/// The closed surface of the full domain: the shortened shell's facets on
/// the moved planes go back to the real walls, and each rim edge of a
/// moved plane gets a slab side (two triangles) with the patch beyond it.
fn full_surface(
    shell: &ClosedSurface,
    shell_patches: &[usize],
    planes: &[YieldPlane],
    total: f64,
    tol: f64,
) -> Result<ClosedSurface, String> {
    let names: Vec<String> = {
        let mut ids: Vec<usize> = shell_patches.to_vec();
        ids.sort_unstable();
        ids.dedup();
        ids.iter().map(|i| i.to_string()).collect()
    };
    let name_of = |p: usize| PatchId(names.iter().position(|s| *s == p.to_string()).unwrap_or(0));
    let mut vertices = shell.vertices().to_vec();
    let plane_of = |x: Vec3| {
        planes
            .iter()
            .position(|q| (q.normal.dot(x) - (q.offset - total)).abs() <= tol)
    };
    let facets = shell.facets();
    let on: Vec<Option<usize>> = facets
        .iter()
        .map(|f| {
            let k = plane_of(vertices[f.triangle.vertices[0]])?;
            (f.plane.normal().dot(planes[k].normal) >= 1.0 - 1e-6
                && f.triangle
                    .vertices
                    .iter()
                    .all(|&v| plane_of(vertices[v]) == Some(k)))
            .then_some(k)
        })
        .collect();
    let mut moved = HashMap::<usize, usize>::new();
    let mut lift = |v: usize, k: usize, vertices: &mut Vec<Vec3>| -> usize {
        *moved.entry(v).or_insert_with(|| {
            vertices.push(vertices[v] + planes[k].normal * total);
            vertices.len() - 1
        })
    };
    let mut owner = HashMap::<Edge, usize>::new();
    for (f, facet) in facets.iter().enumerate() {
        for (a, b) in edges(facet.triangle.vertices) {
            owner.insert((a, b), f);
        }
    }
    let mut triangles = Vec::new();
    for (f, facet) in facets.iter().enumerate() {
        let patch = name_of(shell_patches[f]);
        let Some(k) = on[f] else {
            triangles.push(SurfaceTriangle {
                vertices: facet.triangle.vertices,
                patch,
            });
            continue;
        };
        let t = facet.triangle.vertices;
        triangles.push(SurfaceTriangle {
            vertices: t.map(|v| lift(v, k, &mut vertices)),
            patch,
        });
        for (a, b) in edges(t) {
            let g = owner[&(b, a)];
            if on[g] == Some(k) {
                continue;
            }
            let (a2, b2) = (lift(a, k, &mut vertices), lift(b, k, &mut vertices));
            let side = name_of(shell_patches[g]);
            triangles.push(SurfaceTriangle {
                vertices: [a, b, b2],
                patch: side,
            });
            triangles.push(SurfaceTriangle {
                vertices: [a, b2, a2],
                patch: side,
            });
        }
    }
    // the moved planes' inner vertices are no longer used
    let mut index = vec![usize::MAX; vertices.len()];
    let mut used = Vec::new();
    for t in &mut triangles {
        for v in &mut t.vertices {
            if index[*v] == usize::MAX {
                index[*v] = used.len();
                used.push(vertices[*v]);
            }
            *v = index[*v];
        }
    }
    ClosedSurface::new(used, triangles, names).map_err(|e| format!("full domain surface: {e:?}"))
}
