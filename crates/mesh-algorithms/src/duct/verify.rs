use super::*;
use std::collections::HashMap;

pub struct Measurements {
    pub volume: f64,
    /// Minimum, fifth percentile and median over core polyhedra.
    pub core_sphericity: Option<[f64; 3]>,
    pub minimum_volume: f64,
    pub minimum_wedge_jacobian: Option<f64>,
    pub areas: [f64; 3],
    pub prism_core_faces: usize,
    pub port_layer_areas: [Vec<f64>; 2],
    pub minimum_tet_mean_ratio: f64,
    pub tetrahedra_below_mean_ratio_0_1: usize,
    pub normal_thickness_ranges: Vec<[f64; 2]>,
    pub minimum_kernel_clearance: f64,
}

fn det(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    a.cross(b).dot(c)
}

/// At fixed height the linear-wedge determinant is affine over the reference
/// triangle. On each of its three vertical edges it is quadratic in height.
/// Integrate that polynomial exactly to obtain the physical wedge volume.
fn wedge(p: [Vec3; 6]) -> (f64, f64, [f64; 2]) {
    let [a, b, c, d, e, f] = p;
    let (u, v) = (b - a, c - a);
    let (du, dv) = ((e - d) - u, (f - d) - v);
    let normal = u.cross(v) / u.cross(v).length();
    let mut minimum = f64::INFINITY;
    let mut volume = 0.0;
    let mut thickness = [f64::INFINITY, f64::NEG_INFINITY];
    for w in [d - a, e - b, f - c] {
        let (q, l, h) = (det(du, dv, w), det(du, v, w) + det(u, dv, w), det(u, v, w));
        minimum = minimum.min(h).min(q + l + h);
        if q > 0.0 {
            let t = -l / (2.0 * q);
            if (0.0..1.0).contains(&t) {
                minimum = minimum.min((q * t + l) * t + h);
            }
        }
        volume += (h + l * 0.5 + q / 3.0) / 6.0;
        thickness[0] = thickness[0].min(w.dot(normal));
        thickness[1] = thickness[1].max(w.dot(normal));
    }
    (minimum, volume, thickness)
}

fn key(ring: &[usize]) -> [usize; 4] {
    let mut k = [usize::MAX; 4];
    k[..ring.len()].copy_from_slice(ring);
    k.sort_unstable();
    k
}
fn same_ring(a: &[usize], b: &[usize], reverse: bool) -> bool {
    let n = a.len();
    b.len() == n
        && a.iter().position(|v| *v == b[0]).is_some_and(|start| {
            (0..n).all(|i| a[(start + if reverse { n - i } else { i }) % n] == b[i])
        })
}
/// Integrate the oriented shell from its certified interior point. The expected
/// volume is computed independently from the source tetrahedron partitions.
fn poly_volume(poly: &Polyhedron, points: &[Vec3]) -> Result<(f64, f64), &'static str> {
    let mut edges = HashMap::<(usize, usize), (usize, i32)>::new();
    let mut volume = 0.0;
    let mut area = 0.0;
    for face in &poly.faces {
        for (&a, &b) in face.iter().zip(face.iter().cycle().skip(1)) {
            let e = edges.entry((a.min(b), a.max(b))).or_default();
            e.0 += 1;
            e.1 += if a < b { 1 } else { -1 };
        }
        let a = points[face[0]];
        let first = (points[face[1]] - a).cross(points[face[2]] - a);
        let norm = first.length();
        for pair in face[1..].windows(2) {
            let cross = (points[pair[0]] - a).cross(points[pair[1]] - a);
            let signed = cross.dot(a - poly.interior);
            if !signed.is_finite() || signed <= 0.0 || cross.dot(first) <= 0.0 {
                return Err("inverted face or invalid polyhedron kernel");
            }
            // Face vertices are computed in world coordinates. Account for
            // cancellation in very slender faces near a large coordinate offset.
            let edge = (points[face[1]] - a)
                .length()
                .max((points[face[2]] - a).length());
            let coordinate_scale = a.max_abs_component().max(edge);
            let distance = (points[pair[1]] - a).dot(first).abs();
            if distance
                > 128.0
                    * f64::EPSILON
                    * coordinate_scale
                    * (norm + edge * (points[pair[1]] - a).length())
            {
                return Err("nonplanar polyhedron face");
            }
            volume += signed / 6.0;
            area += 0.5 * cross.length();
        }
    }
    if edges.values().any(|&use_count| use_count != (2, 0))
        || poly.nodes.len() + poly.faces.len() != edges.len() + 2
    {
        return Err("open or nonmanifold polyhedron shell");
    }
    if !volume.is_finite() || (volume / poly.expected_volume - 1.0).abs() > 1e-8 {
        return Err("polyhedron shell volume differs from its source partition");
    }
    Ok((
        volume,
        (36. * std::f64::consts::PI).cbrt() * volume.cbrt().powi(2) / area,
    ))
}

pub(super) fn check(
    surface: &ClosedSurface,
    patches: &[usize],
    p: &[Vec3],
    cells: &[Cell],
    boundary: &[BoundaryFace],
    layers: usize,
) -> Result<Measurements, String> {
    check_with_interface_count(
        surface,
        patches,
        p,
        cells,
        boundary,
        layers,
        patches.iter().filter(|&&p| p == 3).count(),
    )
}

pub(super) fn check_with_interface_count(
    surface: &ClosedSurface,
    patches: &[usize],
    p: &[Vec3],
    cells: &[Cell],
    boundary: &[BoundaryFace],
    layers: usize,
    expected_interfaces: usize,
) -> Result<Measurements, String> {
    let mut exterior = HashMap::new();
    let mut areas = [0.; 3];
    let mut port_layer_areas = [vec![0.; layers + 1], vec![0.; layers + 1]];
    for face in boundary {
        if exterior.insert(key(&face.nodes), face).is_some() {
            return Err("duplicate outer boundary face".into());
        }
        let a = p[face.nodes[0]];
        let area = face.nodes[1..]
            .windows(2)
            .map(|s| (p[s[0]] - a).cross(p[s[1]] - a).length() * 0.5)
            .sum::<f64>();
        areas[face.patch - 1] += area;
        if face.patch <= 2 {
            port_layer_areas[face.patch - 1][cells[face.owner].layer] += area;
        }
    }
    let mut expected = [0.; 3];
    let scale = surface.bounds().extent().max_abs_component();
    for (facet, &patch) in surface.facets().iter().zip(patches) {
        expected[patch - 1] += facet.normalized_area * scale * scale;
    }
    for (i, (&a, &e)) in areas.iter().zip(&expected).enumerate() {
        if !a.is_finite() || !e.is_finite() || e <= 0. || (a / e - 1.).abs() >= 1e-9 {
            return Err(format!(
                "patch {} area {a} differs from STL area {e}",
                i + 1
            ));
        }
    }
    let mut pending = HashMap::<[usize; 4], (Vec<usize>, usize)>::new();
    let mut prism_core_faces = 0;
    let mut complete = std::collections::HashSet::new();
    let mut minimum_volume = f64::INFINITY;
    let mut minimum_j = f64::INFINITY;
    let mut minimum_kernel_clearance = f64::INFINITY;
    let mut normal_thickness_ranges = vec![[f64::INFINITY, 0.0_f64]; layers];
    let mut volume = 0.;
    let mut minimum_tet_mean_ratio = 1.0_f64;
    let mut tetrahedra_below_mean_ratio_0_1 = 0;
    let mut sphericities = Vec::new();
    for (id, cell) in cells.iter().enumerate() {
        let value = match cell.element {
            Element::Polyhedron(ref poly) => {
                let (volume, sphericity) =
                    poly_volume(poly, p).map_err(|e| format!("cell {id}: {e}"))?;
                sphericities.push(sphericity);
                volume
            }
            Element::Tet(n) => {
                let [o, a, b, c] = n.map(|v| p[v]);
                let determinant = det(b - a, c - a, a - o);
                let clearance = determinant / (b - a).cross(c - a).length();
                if !clearance.is_finite() || clearance <= scale * 1e-12 {
                    return Err(format!(
                        "core collision or invalid kernel in region {} at cell {id}, clearance {clearance}",
                        cell.region
                    ));
                }
                minimum_kernel_clearance = minimum_kernel_clearance.min(clearance);
                let squared_edges = [b - a, c - b, a - c, a - o, b - o, c - o]
                    .iter()
                    .map(|e| e.dot(*e))
                    .sum::<f64>();
                let mean_ratio = 12.0 * (determinant / 2.0).cbrt().powi(2) / squared_edges;
                minimum_tet_mean_ratio = minimum_tet_mean_ratio.min(mean_ratio);
                tetrahedra_below_mean_ratio_0_1 += usize::from(mean_ratio < 0.1);
                determinant / 6.
            }
            Element::Wedge(n) => {
                let (jac, value, range) = wedge(n.map(|v| p[v]));
                if !jac.is_finite() || jac <= 0. {
                    return Err(format!(
                        "layer collision or inversion in region {}, layer {}, cell {id}, Jacobian {jac}",
                        cell.region, cell.layer
                    ));
                }
                minimum_j = minimum_j.min(jac);
                let r = &mut normal_thickness_ranges[cell.layer - 1];
                r[0] = r[0].min(range[0]);
                r[1] = r[1].max(range[1]);
                value
            }
        };
        if !value.is_finite() || value <= 0. {
            return Err(format!(
                "nonpositive or unrepresentable volume at cell {id}"
            ));
        }
        volume += value;
        minimum_volume = minimum_volume.min(value);
        for ring in cell.element.faces() {
            let k = key(&ring);
            if complete.contains(&k) {
                return Err("too many incidences on one face".into());
            }
            if let Some(face) = exterior.remove(&k) {
                complete.insert(k);
                if id != face.owner || !same_ring(&face.nodes, &ring, false) {
                    return Err("incorrect boundary owner or winding".into());
                }
            } else if let Some((previous, other_layer)) = pending.remove(&k) {
                if (other_layer == 0) != (cell.layer == 0) {
                    prism_core_faces += 1;
                }
                complete.insert(k);
                if !same_ring(&previous, &ring, true) {
                    return Err("interior face has inconsistent winding".into());
                }
            } else {
                pending.insert(k, (ring, cell.layer));
            }
        }
    }
    if !pending.is_empty() || !exterior.is_empty() {
        return Err("unpaired interior or missing exterior face".into());
    }
    if !volume.is_finite() || (volume / surface.volume() - 1.).abs() >= 1e-9 {
        return Err(format!(
            "mesh volume {volume} differs from STL volume {}",
            surface.volume()
        ));
    }
    if port_layer_areas.iter().flatten().any(|v| *v <= 0.) {
        return Err("a layer does not reach both ports".into());
    }
    if layers > 0 && prism_core_faces != expected_interfaces {
        return Err("incomplete prism-to-core interface".into());
    }
    sphericities.sort_unstable_by(f64::total_cmp);
    let core_sphericity = (!sphericities.is_empty()).then(|| {
        let n = sphericities.len() - 1;
        [sphericities[0], sphericities[n / 20], sphericities[n / 2]]
    });
    Ok(Measurements {
        core_sphericity,
        prism_core_faces,
        minimum_tet_mean_ratio,
        tetrahedra_below_mean_ratio_0_1,
        volume,
        minimum_volume,
        minimum_wedge_jacobian: minimum_j.is_finite().then_some(minimum_j),
        areas,
        port_layer_areas,
        normal_thickness_ranges,
        minimum_kernel_clearance,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wedge_with_positive_end_jacobians_can_invert_halfway() {
        let points = [
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [-2., 0., 1.],
            [0., -1., 1.],
        ]
        .map(|[x, y, z]| Vec3::new(x, y, z));
        // J(t) = (1 - 3t)(1 - 2t), with a negative minimum at t = 5/12.
        let (minimum, volume, _) = wedge(points);
        assert!((minimum + 1.0 / 24.0).abs() < 1e-14);
        assert!(volume > 0.0);
    }
}
