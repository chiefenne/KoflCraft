use super::*;

pub(super) struct Front {
    pub points: Vec<Vec3>,
    pub cells: Vec<Cell>,
    pub boundary: Vec<BoundaryFace>,
    pub inner: Vec<usize>,
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
    rims: &[Vec<Edge>; 2],
    steps: &[f64],
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
        if patch != 3 {
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
    for (index, edges) in rims.iter().enumerate() {
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
            if patch != 3 {
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
    }
    let rim_patches: BTreeMap<_, _> = rims
        .iter()
        .enumerate()
        .flat_map(|(i, r)| r.iter().map(move |&edge| (edge, i + 1)))
        .collect();
    let mut distance = 0.0;
    let triangles: Vec<_> = surface
        .facets()
        .iter()
        .map(|f| f.triangle.vertices)
        .collect();
    for (level, &step) in steps.iter().enumerate() {
        distance += step;
        let displaced: Vec<_> = original
            .iter()
            .zip(&direction)
            .map(|(&p, &d)| p + d * distance)
            .collect();
        super::intersection::embedded(
            &displaced,
            &triangles,
            surface.bounds().extent().max_abs_component(),
        )
        .map_err(|e| format!("layer {} front: {e}", level + 1))?;
        let mut next = inner.clone();
        for (i, &wall) in on_wall.iter().enumerate() {
            if wall {
                next[i] = points.len();
                points.push(displaced[i]);
            }
        }
        for (i, f) in surface.facets().iter().enumerate() {
            if patches[i] != 3 {
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
                    patch: 3,
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
    })
}
