use super::*;

pub(super) struct Front {
    pub points: Vec<Vec3>,
    pub cells: Vec<Cell>,
    pub boundary: Vec<BoundaryFace>,
    pub inner: Vec<usize>,
    /// The share of the layer height at each surface vertex (see `squeeze`).
    pub scale: Vec<f64>,
}

/// Largest distance of a port vertex from the port plane, relative to the
/// longest extent: a few f32 roundings of the coordinates.
pub(crate) const PORT_RELATIVE_TOLERANCE: f64 = 1e-6;
/// Largest 1 - cos of the angle between two port facet normals (0.08°).
pub(crate) const PORT_NORMAL_SLACK: f64 = 1e-6;

pub(super) fn advance(
    surface: &ClosedSurface,
    patches: &[usize],
    labels: &[usize],
    rims: &BTreeMap<usize, Vec<Edge>>,
    steps: &[f64],
    limit: Option<&[f64]>,
) -> Result<Front, String> {
    let original = surface.vertices();
    let mut points = original.to_vec();
    let total: f64 = steps.iter().sum();
    if !total.is_finite() || steps.iter().any(|&s| !s.is_finite() || s <= 0.0) {
        return Err("layer thicknesses must be finite and positive".into());
    }
    let tol = surface.bounds().extent().max_abs_component() * 1e-10;
    let mut direction = vec![Vec3::ZERO; original.len()];
    let mut on_wall = vec![false; original.len()];
    let mut port_vertex = vec![false; original.len()];
    let mut normals = vec![Vec::<Vec3>::new(); original.len()];
    for (f, &patch) in surface.facets().iter().zip(patches) {
        if !is_wall(patch) {
            continue;
        }
        let ids = f.triangle.vertices;
        for i in 0..3 {
            let p = original[ids[i]];
            let u = original[ids[(i + 1) % 3]] - p;
            let v = original[ids[(i + 2) % 3]] - p;
            let angle = u.cross(v).length().atan2(u.dot(v));
            let inward = f.plane.normal() * -1.0;
            let groups = &mut normals[ids[i]];
            if let Some(group) = groups.iter_mut().find(|v| v.dot(inward) > 0.9 * v.length()) {
                *group = *group + inward * angle;
            } else {
                groups.push(inward * angle);
            }
            on_wall[ids[i]] = true;
        }
    }
    if !steps.is_empty() {
        for (i, groups) in normals.iter().enumerate() {
            // Miter sharp features by satisfying a unit normal displacement for
            // each distinct wall sector. Orthogonalization keeps prior constraints.
            let mut basis = Vec::<Vec3>::new();
            for group in groups {
                let normal = *group / group.length();
                let tangent = basis.iter().fold(normal, |v, &axis| v - axis * v.dot(axis));
                let residual = 1.0 - direction[i].dot(normal);
                if tangent.length() > 1e-6 && basis.len() < 3 {
                    direction[i] = direction[i] + tangent * (residual / tangent.dot(tangent));
                    basis.push(tangent / tangent.length());
                } else if residual.abs() > 1e-3 {
                    return Err(format!("incompatible sharp wall normals at vertex {i}"));
                }
            }
        }
    }
    // single-precision STL: a tilted cap's corners are coplanar only to the
    // rounding of f32 coordinates
    let port_tol = surface.bounds().extent().max_abs_component() * PORT_RELATIVE_TOLERANCE;
    let mut slip_rim = vec![false; original.len()];
    let slip_inner = if steps.is_empty() {
        Vec::new()
    } else {
        slip_patches(
            surface,
            patches,
            port_tol,
            &on_wall,
            &mut direction,
            &mut slip_rim,
        )?
    };
    for (&id, edges) in rims {
        if edges.is_empty() || is_slip(id) {
            continue;
        }
        let index = id - 1;
        let ids = loop_vertices(edges)?;
        let center = mean(original, &ids);
        let mut caps = surface
            .facets()
            .iter()
            .zip(patches)
            .filter(|&(_, &p)| p == index + 1)
            .map(|(f, _)| f);
        let normal = caps.next().ok_or("missing port facets")?.plane.normal();
        let mut radius = f64::INFINITY;
        for &(a, b) in edges {
            let u = original[a] - center;
            let v = original[b] - center;
            radius = radius.min(v.cross(u).dot(normal) / (v - u).length());
        }
        if radius <= tol || total >= radius {
            return Err(format!("layer stack closes port {}", index + 1));
        }
        for (f, &patch) in surface.facets().iter().zip(patches) {
            if patch != index + 1 {
                continue;
            }
            if f.plane.normal().dot(normal) < 1.0 - PORT_NORMAL_SLACK {
                return Err("port patch must be planar and consistently oriented".into());
            }
            for &v in &f.triangle.vertices {
                let r = original[v] - center;
                if r.dot(normal).abs() > port_tol {
                    return Err("off-plane port vertex".into());
                }
                // A single in-plane similarity preserves every original cap triangle.
                direction[v] = (normal * r.dot(normal) - r) / radius;
                port_vertex[v] = true;
                if !on_wall[v] {
                    points[v] = original[v] + direction[v] * total;
                }
            }
        }
    }
    let mut inner: Vec<_> = (0..original.len()).collect();
    let mut cells = Vec::new();
    let mut boundary = Vec::new();
    if !steps.is_empty() {
        // At sharp edges the bisector must travel farther than either individual
        // face offset. Scale to the least normal advance of its incident walls.
        let mut projection = vec![f64::INFINITY; original.len()];
        for (f, &patch) in surface.facets().iter().zip(patches) {
            if !is_wall(patch) {
                continue;
            }
            for &v in &f.triangle.vertices {
                projection[v] = projection[v].min(-direction[v].dot(f.plane.normal()));
            }
        }
        for i in 0..original.len() {
            if on_wall[i] {
                if projection[i] <= 1e-10 || !projection[i].is_finite() {
                    return Err(format!("no common inward wall direction at vertex {i}"));
                }
                if !port_vertex[i] {
                    direction[i] = direction[i] / projection[i];
                }
            }
        }
        for (normal, inner_ids, faces) in &slip_inner {
            spread_in_plane(*normal, inner_ids, faces, &mut direction);
            for &v in inner_ids {
                points[v] = original[v] + direction[v] * total;
            }
            let moved = |v: usize| original[v] + direction[v] * total;
            let folded = faces
                .iter()
                .filter(|t| {
                    let [a, b, c] = t.map(moved);
                    (b - a).cross(c - a).dot(*normal) <= 0.0
                })
                .count();
            if folded > 0 {
                return Err(format!(
                    "{folded} facets of a patch without layers fold where the layers meet it"
                ));
            }
        }
    }
    let rim_patches: BTreeMap<_, _> = rims
        .iter()
        .flat_map(|(&id, r)| r.iter().map(move |&edge| (edge, id)))
        .collect();
    let triangles: Vec<_> = surface
        .facets()
        .iter()
        .map(|f| f.triangle.vertices)
        .collect();
    let extent = surface.bounds().extent().max_abs_component();
    let scale = squeeze(
        surface,
        patches,
        &triangles,
        &direction,
        &on_wall,
        &port_vertex,
        &slip_rim,
        steps,
        limit,
    )?;
    let mut distance = 0.0;
    for (level, &step) in steps.iter().enumerate() {
        distance += step;
        let displaced: Vec<_> = (0..original.len())
            .map(|i| original[i] + direction[i] * (distance * scale[i]))
            .collect();
        super::intersection::embedded(&displaced, &triangles, extent)
            .map_err(|e| format!("layer {} front: {e}", level + 1))?;
        let mut next = inner.clone();
        for (i, &wall) in on_wall.iter().enumerate() {
            if wall {
                next[i] = points.len();
                points.push(displaced[i]);
            }
        }
        for (i, f) in surface.facets().iter().enumerate() {
            if !is_wall(patches[i]) {
                continue;
            }
            let [a, b, c] = f.triangle.vertices;
            let owner = cells.len();
            cells.push(Cell {
                element: Element::Wedge([inner[a], inner[c], inner[b], next[a], next[c], next[b]]),
                layer: level + 1,
                region: labels[i],
            });
            if level == 0 {
                boundary.push(BoundaryFace {
                    nodes: vec![a, b, c],
                    patch: patches[i],
                    owner,
                });
            }
            for (a, b) in edges([a, b, c]) {
                if let Some(&patch) = rim_patches.get(&(a, b)) {
                    boundary.push(BoundaryFace {
                        nodes: vec![inner[b], inner[a], next[a], next[b]],
                        patch,
                        owner,
                    });
                }
            }
        }
        inner = next;
    }
    Ok(Front {
        points,
        cells,
        boundary,
        inner,
        scale,
    })
}

/// Plane normal, inner vertices and facets of a wall patch without layers.
type SlipPatch = (Vec3, Vec<usize>, Vec<[usize; 3]>);

/// Split the facets of the wall patches without layers into edge-connected
/// planes, check that each is planar, turn the layer direction of their rim
/// vertices into their plane (the layers of the walls around run straight
/// into them) and return each plane with its inner vertices and facets.
fn slip_patches(
    surface: &ClosedSurface,
    patches: &[usize],
    port_tol: f64,
    on_wall: &[bool],
    direction: &mut [Vec3],
    slip_rim: &mut [bool],
) -> Result<Vec<SlipPatch>, String> {
    let original = surface.vertices();
    let facets = surface.facets();
    let slip: Vec<usize> = (0..facets.len()).filter(|&f| is_slip(patches[f])).collect();
    if slip.is_empty() {
        return Ok(Vec::new());
    }
    // edge-connected components of the slip facets
    let mut parent: BTreeMap<usize, usize> = slip.iter().map(|&f| (f, f)).collect();
    fn root(parent: &mut BTreeMap<usize, usize>, mut f: usize) -> usize {
        while parent[&f] != f {
            let up = parent[&parent[&f]];
            parent.insert(f, up);
            f = up;
        }
        f
    }
    let mut by_edge = BTreeMap::<Edge, usize>::new();
    for &f in &slip {
        for (a, b) in edges(facets[f].triangle.vertices) {
            if let Some(g) = by_edge.remove(&(b, a)) {
                let (x, y) = (root(&mut parent, f), root(&mut parent, g));
                parent.insert(x, y);
            } else {
                by_edge.insert((a, b), f);
            }
        }
    }
    let mut groups = BTreeMap::<usize, Vec<usize>>::new();
    for &f in &slip {
        let r = root(&mut parent, f);
        groups.entry(r).or_default().push(f);
    }
    let mut result = Vec::new();
    for faces in groups.into_values() {
        let first = &facets[faces[0]];
        let normal = first.plane.normal();
        let origin = original[first.triangle.vertices[0]];
        let name = patch_index(patches[faces[0]]);
        let mut inner = Vec::new();
        for &f in &faces {
            let facet = &facets[f];
            let flat = facet.plane.normal().dot(normal) >= 1.0 - PORT_NORMAL_SLACK
                && facet
                    .triangle
                    .vertices
                    .iter()
                    .all(|&v| (original[v] - origin).dot(normal).abs() <= port_tol);
            if !flat {
                return Err(format!(
                    "patch {name} has no layers, so each connected piece of it must be planar"
                ));
            }
            for &v in &facet.triangle.vertices {
                if on_wall[v] {
                    if !slip_rim[v] {
                        slip_rim[v] = true;
                        direction[v] = direction[v] - normal * direction[v].dot(normal);
                    }
                } else {
                    inner.push(v);
                }
            }
        }
        // Vertices on a border between two patches without layers stay put,
        // so neither patch changes its area.
        let mut owner = BTreeMap::<usize, usize>::new();
        let mut pinned = std::collections::BTreeSet::new();
        for &f in &faces {
            for &v in &facets[f].triangle.vertices {
                if *owner.entry(v).or_insert(patches[f]) != patches[f] {
                    pinned.insert(v);
                }
            }
        }
        inner.retain(|v| !pinned.contains(v));
        inner.sort_unstable();
        inner.dedup();
        let triangles = faces.iter().map(|&f| facets[f].triangle.vertices).collect();
        result.push((normal, inner, triangles));
    }
    Ok(result)
}

/// Spread the in-plane layer offset of the rim of a patch without layers
/// over its inner vertices (a harmonic field, by over-relaxed Gauss-Seidel
/// sweeps), so its facets give way to the layers as a cap contracts.
fn spread_in_plane(normal: Vec3, inner: &[usize], faces: &[[usize; 3]], direction: &mut [Vec3]) {
    let mut neighbours = BTreeMap::<usize, Vec<usize>>::new();
    for &t in faces {
        for (a, b) in edges(t) {
            neighbours.entry(a).or_default().push(b);
            neighbours.entry(b).or_default().push(a);
        }
    }
    for list in neighbours.values_mut() {
        list.sort_unstable();
        list.dedup();
    }
    for &v in inner {
        direction[v] = Vec3::ZERO;
    }
    let largest = neighbours
        .keys()
        .map(|&v| direction[v].length())
        .fold(0.0, f64::max);
    for _ in 0..SPREAD_SWEEPS {
        let mut change = 0.0_f64;
        for &v in inner {
            let list = &neighbours[&v];
            let mean = list.iter().fold(Vec3::ZERO, |s, &u| s + direction[u]) / list.len() as f64;
            let step = (mean - direction[v]) * SPREAD_RELAXATION;
            direction[v] = direction[v] + step;
            change = change.max(step.length());
        }
        if change <= SPREAD_TOLERANCE * largest {
            break;
        }
    }
    for &v in inner {
        direction[v] = direction[v] - normal * direction[v].dot(normal);
    }
}

const SPREAD_SWEEPS: usize = 20_000;
const SPREAD_RELAXATION: f64 = 1.9;
const SPREAD_TOLERANCE: f64 = 1e-9;

/// Share of the local gap that the layers of one wall may fill.
const GAP_SHARE: f64 = 0.3;
/// Largest change of the layer height per unit length along the wall.
const HEIGHT_RATE: f64 = 0.3;
/// Factor on the height of the vertices of colliding front facets.
const COLLISION_FACTOR: f64 = 0.7;
const COLLISION_ROUNDS: usize = 60;
/// Layers are not squeezed below this share of their height.
const MIN_SCALE: f64 = 0.02;

/// The share of the layer height each wall vertex gets (1: all of it).
///
/// Layers are squeezed where they do not fit, with the same number of
/// layers everywhere (the height scaling of Garimella & Shephard 2000 and
/// Dyedov et al. 2009): first each wall vertex is limited to `GAP_SHARE`
/// of the gap in front of it (the medial-ball thickness of
/// `sizing::thickness`), then the vertices of front facets that still
/// intersect are scaled by `COLLISION_FACTOR` until every front is
/// embedded. The height changes by at most `HEIGHT_RATE` per unit length
/// along the wall. Vertices on the port rims keep the full height, so the
/// caps stay similar to themselves.
#[allow(clippy::too_many_arguments)]
fn squeeze(
    surface: &ClosedSurface,
    patches: &[usize],
    triangles: &[[usize; 3]],
    direction: &[Vec3],
    on_wall: &[bool],
    port_vertex: &[bool],
    slip_rim: &[bool],
    steps: &[f64],
    start: Option<&[f64]>,
) -> Result<Vec<f64>, String> {
    use crate::sizing::surface::SurfaceGeometry;
    use crate::sizing::thickness::OrientedSet;
    let original = surface.vertices();
    let mut scale = start.map_or_else(|| vec![1.0; original.len()], <[f64]>::to_vec);
    let total: f64 = steps.iter().sum();
    if steps.is_empty() {
        return Ok(scale);
    }
    let free = |v: usize| on_wall[v] && !port_vertex[v] && !slip_rim[v];
    let geometry = SurfaceGeometry::new(surface);
    let set = OrientedSet::new(&geometry);
    let wedge = 45f64.to_radians();
    let cap = 1.01 * total / GAP_SHARE;
    for (f, t) in triangles.iter().enumerate() {
        if !is_wall(patches[f]) {
            continue;
        }
        let n = geometry.normals[f];
        for &v in t {
            if !free(v) {
                continue;
            }
            if let Some(gap) = set.thickness(original[v], n, n * -1.0, wedge, cap) {
                scale[v] = scale[v].min(GAP_SHARE * gap / total);
            }
        }
    }
    let mut neighbours = vec![Vec::new(); original.len()];
    for (f, t) in triangles.iter().enumerate() {
        if !is_wall(patches[f]) {
            continue;
        }
        for (a, b) in edges(*t) {
            neighbours[a].push(b);
            neighbours[b].push(a);
        }
    }
    let limit = |scale: &mut Vec<f64>| {
        let mut queue: std::collections::VecDeque<usize> =
            (0..original.len()).filter(|&v| scale[v] < 1.0).collect();
        while let Some(v) = queue.pop_front() {
            for &u in &neighbours[v] {
                if !free(u) {
                    continue;
                }
                let bound = scale[v] + HEIGHT_RATE * (original[u] - original[v]).length() / total;
                if scale[u] > bound {
                    scale[u] = bound;
                    queue.push_back(u);
                }
            }
        }
    };
    limit(&mut scale);
    let extent = surface.bounds().extent().max_abs_component();
    let mut rounds = 0;
    loop {
        // every front must be embedded; fix the outermost first
        let mut hit = Vec::new();
        let mut distance = total;
        for level in (0..steps.len()).rev() {
            let displaced: Vec<_> = (0..original.len())
                .map(|i| original[i] + direction[i] * (distance * scale[i]))
                .collect();
            hit = super::intersection::collisions(&displaced, triangles, extent, usize::MAX);
            if !hit.is_empty() {
                break;
            }
            distance -= steps[level];
        }
        if hit.is_empty() {
            break;
        }
        rounds += 1;
        let mut touched = vec![false; original.len()];
        for (i, j) in hit {
            for v in triangles[i].into_iter().chain(triangles[j]) {
                touched[v] = true;
            }
        }
        let mut changed = false;
        for v in 0..original.len() {
            if touched[v] && free(v) && scale[v] * COLLISION_FACTOR >= MIN_SCALE {
                scale[v] *= COLLISION_FACTOR;
                changed = true;
            }
        }
        if !changed || rounds > COLLISION_ROUNDS {
            // the caller's embedding check reports the place
            break;
        }
        limit(&mut scale);
    }
    if crate::trace_enabled() {
        let squeezed = (0..original.len())
            .filter(|&v| free(v) && scale[v] < 1.0)
            .count();
        let least = scale.iter().copied().fold(1.0, f64::min);
        eprintln!(
            "kofl trace: layers squeezed at {squeezed} of {} wall vertices, least share {least:.3}, {rounds} collision rounds",
            on_wall.iter().filter(|&&w| w).count()
        );
    }
    Ok(scale)
}
