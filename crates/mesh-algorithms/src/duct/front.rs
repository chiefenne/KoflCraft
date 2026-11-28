use super::*;

pub(super) struct Front {
    pub points: Vec<Vec3>,
    pub cells: Vec<Cell>,
    pub boundary: Vec<BoundaryFace>,
    pub inner: Vec<usize>,
    /// The share of the layer height at each surface vertex (see `squeeze`).
    pub scale: Vec<f64>,
    /// The surface vertex each of the points above grows from (itself for
    /// a surface vertex).
    pub wall_of: Vec<usize>,
}

/// Largest distance of a port vertex from the port plane, relative to the
/// longest extent: a few f32 roundings of the coordinates.
pub(crate) const PORT_RELATIVE_TOLERANCE: f64 = 1e-6;
/// Largest 1 - cos of the angle between two port facet normals (0.08°).
pub(crate) const PORT_NORMAL_SLACK: f64 = 1e-6;

/// The triangles that may not collide: `triangles` with every vertex copy
/// replaced by its original (`twins`, see [`advance`]), each coincident
/// pair kept once.
pub(super) fn contact_triangles(triangles: &[[usize; 3]], twins: &[usize]) -> Vec<[usize; 3]> {
    if twins.is_empty() {
        return triangles.to_vec();
    }
    let mut seen = std::collections::BTreeSet::new();
    triangles
        .iter()
        .map(|t| t.map(|v| twins[v]))
        .filter(|t| {
            let mut k = *t;
            k.sort_unstable();
            seen.insert(k)
        })
        .collect()
}

/// The wall layers. `twins` (empty, or the original of every surface
/// vertex) names the vertex copies of a cracked zero-thickness shell
/// (SHL-01): a copy and its original lie at one place and move as one
/// vertex, with the wall normals of both sides.
pub(super) fn advance(
    surface: &ClosedSurface,
    patches: &[usize],
    labels: &[usize],
    rims: &BTreeMap<usize, Vec<Edge>>,
    steps: &[f64],
    limit: Option<&[f64]>,
    twins: &[usize],
) -> Result<Front, String> {
    let original = surface.vertices();
    if !twins.is_empty() && twins.len() != original.len() {
        return Err("invalid vertex copies".into());
    }
    let canon = |v: usize| twins.get(v).copied().unwrap_or(v);
    let unify_min = |values: &mut [f64]| {
        for v in 0..values.len() {
            let c = canon(v);
            if c != v {
                let m = values[v].min(values[c]);
                values[v] = m;
                values[c] = m;
            }
        }
        for v in 0..values.len() {
            values[v] = values[canon(v)];
        }
    };
    let start: Option<Vec<f64>> = limit.map(|l| {
        let mut l = l.to_vec();
        unify_min(&mut l);
        l
    });
    let limit = start.as_deref();
    let mut points = original.to_vec();
    let mut wall_of: Vec<usize> = (0..original.len()).collect();
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
            let groups = &mut normals[canon(ids[i])];
            if let Some(group) = groups.iter_mut().find(|v| v.dot(inward) > 0.9 * v.length()) {
                *group = *group + inward * angle;
            } else {
                groups.push(inward * angle);
            }
            on_wall[ids[i]] = true;
            on_wall[canon(ids[i])] = true;
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
        for v in 0..original.len() {
            direction[v] = direction[canon(v)];
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
    let mut frames: Vec<(Vec3, Vec3, f64, Vec<usize>)> = Vec::new();
    for (&id, edges) in rims {
        if edges.is_empty() || is_slip(id) {
            continue;
        }
        // one similarity per cap: a cap id may hold several planar caps
        // (the cuts of a sweep fill behind which the region has more than
        // two openings), each its own connected set of facets and rim loop
        for (cap, rim) in caps_of(surface, patches, id, edges)? {
            let ids = loop_vertices(&rim)?;
            let center = mean(original, &ids);
            let normal = surface.facets()[cap[0]].plane.normal();
            let mut radius = f64::INFINITY;
            for &(a, b) in &rim {
                let u = original[a] - center;
                let v = original[b] - center;
                radius = radius.min(v.cross(u).dot(normal) / (v - u).length());
            }
            if radius <= tol || total >= radius {
                return Err(format!("layer stack closes port {id}"));
            }
            frames.push((center, normal, radius, ids.clone()));
            for &f in &cap {
                let f = &surface.facets()[f];
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
    }
    if !steps.is_empty() {
        fade_cap_shift(
            surface,
            patches,
            &frames,
            total,
            &normals,
            &port_vertex,
            &mut direction,
        );
        for v in 0..original.len() {
            direction[v] = direction[canon(v)];
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
                let c = canon(v);
                projection[c] = projection[c].min(-direction[v].dot(f.plane.normal()));
            }
        }
        for i in 0..original.len() {
            if on_wall[i] {
                let projection = projection[canon(i)];
                if projection <= 1e-10 || !projection.is_finite() {
                    return Err(format!("no common inward wall direction at vertex {i}"));
                }
                if !port_vertex[i] {
                    direction[i] = direction[i] / projection;
                }
            }
        }
        for (normal, inner_ids, faces) in &slip_inner {
            spread_in_plane(*normal, inner_ids, faces, &mut direction);
        }
        for v in 0..original.len() {
            direction[v] = direction[canon(v)];
        }
        free_edges_keep_their_line(surface, patches, &on_wall, &mut direction)?;
        for (normal, inner_ids, faces) in &slip_inner {
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
    let contact = contact_triangles(&triangles, twins);
    let mut scale = squeeze(
        surface,
        patches,
        &triangles,
        &contact,
        &direction,
        &on_wall,
        &port_vertex,
        &slip_rim,
        steps,
        limit,
    )?;
    unify_min(&mut scale);
    let mut distance = 0.0;
    for (level, &step) in steps.iter().enumerate() {
        distance += step;
        let displaced: Vec<_> = (0..original.len())
            .map(|i| original[i] + direction[i] * (distance * scale[i]))
            .collect();
        super::intersection::embedded(&displaced, &contact, extent)
            .map_err(|e| format!("layer {} front: {e}", level + 1))?;
        let mut next = inner.clone();
        for (i, &wall) in on_wall.iter().enumerate() {
            if wall {
                next[i] = points.len();
                points.push(displaced[i]);
                wall_of.push(i);
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
        wall_of,
    })
}

/// One planar cap: its facets and its rim edges.
type Cap = (Vec<usize>, Vec<Edge>);

/// Distance from a cap plane, in layer stacks, over which the wall layers
/// take the cap's sideways shift.
const CAP_FADE: f64 = 4.0;

/// The wall layers next to a cap take its sideways shift, fading out.
/// The cap's similarity moves a rim vertex x by (c - x) / r per unit
/// height (in the cap plane): along the wall's normal by one, and, where
/// the rim is not a circle about c, also along the wall (up to one at a
/// square's corners). A wall vertex next to the rim moves along the
/// normal only, so the layer prisms between the rim and the next wall
/// vertices shear by that sideways shift over one wall edge, and thin
/// scaffold tetrahedra there turn over (the square cuts of a sweep fill).
/// Each wall vertex reached from the rim over wall edges within
/// `CAP_FADE` layer stacks of the plane, on walls across the plane, takes
/// the tangential part of the similarity at its place, times a smoothstep
/// from one at the plane to zero at that distance. The part along the wall
/// normals is taken out, so the layer heights do not change.
fn fade_cap_shift(
    surface: &ClosedSurface,
    patches: &[usize],
    frames: &[(Vec3, Vec3, f64, Vec<usize>)],
    total: f64,
    normals: &[Vec<Vec3>],
    port_vertex: &[bool],
    direction: &mut [Vec3],
) {
    let points = surface.vertices();
    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); points.len()];
    for (f, &patch) in surface.facets().iter().zip(patches) {
        if !is_wall(patch) {
            continue;
        }
        for (a, b) in edges(f.triangle.vertices) {
            neighbours[a].push(b);
            neighbours[b].push(a);
        }
    }
    let reach = CAP_FADE * total;
    let mut shift = vec![Vec3::ZERO; points.len()];
    for (center, normal, radius, rim) in frames {
        let (c, n, r) = (*center, *normal, *radius);
        let mut seen: BTreeSet<usize> = rim.iter().copied().collect();
        let mut stack: Vec<usize> = rim.clone();
        while let Some(v) = stack.pop() {
            for &w in &neighbours[v] {
                if port_vertex[w] || !seen.insert(w) {
                    continue;
                }
                let x = points[w];
                let s = (c - x).dot(n);
                if !(0.0..reach).contains(&s) {
                    continue;
                }
                // walls across the plane only (not a wall the plane runs along)
                let basis: Vec<Vec3> = normals[w].iter().map(|g| *g / g.length()).collect();
                if basis.is_empty()
                    || basis.len() >= 3
                    || basis.iter().any(|g| g.dot(n).abs() > 0.3)
                {
                    continue;
                }
                let d = x - c;
                let mut t = (n * d.dot(n) - d) / r;
                // keep the motion along every wall normal (orthonormalised)
                let mut axes: Vec<Vec3> = Vec::new();
                for g in &basis {
                    let a = axes.iter().fold(*g, |v, &e| v - e * v.dot(e));
                    if a.length() > 1e-6 {
                        axes.push(a / a.length());
                    }
                }
                for e in &axes {
                    t = t - *e * t.dot(*e);
                }
                let u = s / reach;
                let weight = 1.0 - u * u * (3.0 - 2.0 * u);
                shift[w] = shift[w] + t * weight;
                stack.push(w);
            }
        }
    }
    for v in 0..points.len() {
        direction[v] = direction[v] + shift[v];
    }
}

/// The caps of the cap id `id`: each connected set of its facets (over
/// shared edges) with the rim edges on it. A cap id usually holds one cap;
/// the region behind the cuts of a sweep fill can hold several under one
/// id. Two caps of one id that touch at a vertex are refused.
fn caps_of(
    surface: &ClosedSurface,
    patches: &[usize],
    id: usize,
    rim: &[Edge],
) -> Result<Vec<Cap>, String> {
    let facets: Vec<usize> = (0..patches.len()).filter(|&f| patches[f] == id).collect();
    if facets.is_empty() {
        return Err("missing port facets".into());
    }
    let mut parent: Vec<usize> = (0..facets.len()).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    let mut by_edge = BTreeMap::<Edge, usize>::new();
    for (k, &f) in facets.iter().enumerate() {
        for (a, b) in edges(surface.facets()[f].triangle.vertices) {
            let e = (a.min(b), a.max(b));
            if let Some(&j) = by_edge.get(&e) {
                let (ra, rb) = (root(&mut parent, k), root(&mut parent, j));
                parent[ra] = rb;
            } else {
                by_edge.insert(e, k);
            }
        }
    }
    let mut group = BTreeMap::<usize, usize>::new();
    let mut caps: Vec<Cap> = Vec::new();
    let mut cap_of_vertex = BTreeMap::<usize, usize>::new();
    for (k, &f) in facets.iter().enumerate() {
        let r = root(&mut parent, k);
        let next = caps.len();
        let c = *group.entry(r).or_insert(next);
        if c == next {
            caps.push((Vec::new(), Vec::new()));
        }
        caps[c].0.push(f);
        for &v in &surface.facets()[f].triangle.vertices {
            if *cap_of_vertex.entry(v).or_insert(c) != c {
                return Err(format!("two caps of patch {id} touch at a vertex"));
            }
        }
    }
    for &(a, b) in rim {
        let c = *cap_of_vertex
            .get(&a)
            .ok_or("a port rim edge is not on a port facet")?;
        caps[c].1.push((a, b));
    }
    Ok(caps)
}

/// The free edges of a zero-thickness shell (where its two sides meet back
/// to back) keep their line, so the shell keeps its outline and area: a
/// vertex inside a straight free edge slides along it, a free-edge corner
/// stays, and where a free edge meets a layered wall, the wall's layer
/// direction there is turned along the free edge.
fn free_edges_keep_their_line(
    surface: &ClosedSurface,
    patches: &[usize],
    on_wall: &[bool],
    direction: &mut [Vec3],
) -> Result<(), String> {
    let facets = surface.facets();
    let original = surface.vertices();
    let mut by_edge = BTreeMap::<Edge, usize>::new();
    let mut tangents = BTreeMap::<usize, Vec<Vec3>>::new();
    for (f, facet) in facets.iter().enumerate() {
        if !is_slip(patches[f]) {
            continue;
        }
        for (a, b) in edges(facet.triangle.vertices) {
            match by_edge.remove(&(b, a)) {
                Some(g)
                    if facet.plane.normal().dot(facets[g].plane.normal())
                        < -1.0 + PORT_NORMAL_SLACK =>
                {
                    let t = original[b] - original[a];
                    let t = t / t.length();
                    tangents.entry(a).or_default().push(t);
                    tangents.entry(b).or_default().push(t * -1.0);
                }
                Some(_) => {}
                None => {
                    by_edge.insert((a, b), f);
                }
            }
        }
    }
    for (v, t) in tangents {
        let d = direction[v];
        direction[v] = match t[..] {
            [t0] if on_wall[v] => {
                let along = d.dot(t0);
                if d.dot(d) == 0.0 {
                    d
                } else if along <= 0.3 * d.length() {
                    return Err(
                        "a free edge of a shell meets a layered wall at too shallow an angle for the layers to follow it"
                            .into(),
                    );
                } else {
                    t0 * (d.dot(d) / along)
                }
            }
            [t0, t1] if !on_wall[v] && t0.cross(t1).length() < 1e-6 => t0 * d.dot(t0),
            _ => Vec3::ZERO,
        };
    }
    Ok(())
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
                // the two sides of a zero-thickness shell meet back to back
                // at its free edges: each side is a plane of its own
                let back_to_back = facets[f].plane.normal().dot(facets[g].plane.normal())
                    < -1.0 + PORT_NORMAL_SLACK;
                if back_to_back {
                    continue;
                }
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
    contact: &[[usize; 3]],
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
            hit = super::intersection::collisions(&displaced, contact, extent, usize::MAX);
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
            for v in contact[i].into_iter().chain(contact[j]) {
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
