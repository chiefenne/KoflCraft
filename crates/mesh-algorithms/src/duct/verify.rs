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
pub(super) fn wedge(p: [Vec3; 6]) -> (f64, f64, [f64; 2]) {
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

fn key(ring: &[usize]) -> Vec<usize> {
    let mut k = ring.to_vec();
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
/// A face counts as planar when its corners are within this share of its
/// size of the fitted plane. Such a face is fanned from its vertex mean, as
/// FoehnCraft does with every face, so the checked volume is the solver's.
/// The core welds points up to 1e-10 of the model size apart, which moves a
/// corner of a small face by more than 1e-10 of that face's size.
const PLANAR_SHARE: f64 = 1e-7;

/// Two-point Gauss rule on [0, 1].
const GAUSS: [f64; 2] = [0.211_324_865_405_187_1, 0.788_675_134_594_812_9];

/// A face integrated from `origin`: (six times the cone volume, area), and
/// whether every part of it is seen strictly from the front.
///
/// Planar polygons are fanned from their vertex mean, as FoehnCraft does. A
/// quadrilateral that is not planar is the bilinear patch of a layer prism
/// side (as in a linear wedge); its cone volume is integrated exactly by the
/// 2×2 Gauss rule (the integrand has degree two per parameter).
fn face_integral(
    points: &[Vec3],
    face: &[usize],
    origin: Vec3,
) -> Result<(f64, f64, bool), &'static str> {
    let (newell, centre) = super::convex::ring_plane(points, face);
    let norm = newell.length();
    if norm.is_nan() || norm <= 0.0 {
        return Err("degenerate polyhedron face");
    }
    let unit = super::convex::fitted_normal(points, face, newell, centre);
    let diameter = face
        .iter()
        .map(|&v| (points[v] - centre).length())
        .fold(0.0, f64::max);
    let coordinate_scale = centre.max_abs_component().max(diameter);
    // a triangle is planar; on a sliver its Newell normal is too inaccurate
    // for the offset test
    let planar = face.len() == 3
        || face.iter().all(|&v| {
            (points[v] - centre).dot(unit).abs()
                <= PLANAR_SHARE * diameter + 256.0 * f64::EPSILON * coordinate_scale
        });
    if !planar {
        if face.len() != 4 {
            return Err("nonplanar polyhedron face");
        }
        let [p0, p1, p2, p3] = [0, 1, 2, 3].map(|i| points[face[i]]);
        let mut v6 = 0.0;
        let mut area = 0.0;
        let mut front = true;
        let at = |u: f64, v: f64| {
            let x = p0 * ((1. - u) * (1. - v))
                + p1 * (u * (1. - v))
                + p2 * (u * v)
                + p3 * ((1. - u) * v);
            let xu = (p1 - p0) * (1. - v) + (p2 - p3) * v;
            let xv = (p3 - p0) * (1. - u) + (p2 - p1) * u;
            (x, xu.cross(xv))
        };
        for &u in &GAUSS {
            for &v in &GAUSS {
                let (x, n) = at(u, v);
                // 2 = 6 / 3 for the cone, 1/4 for the Gauss weights
                v6 += 0.5 * (x - origin).dot(n);
                area += 0.25 * n.length();
                front &= (x - origin).dot(n) > 0.0;
            }
        }
        for (u, v) in [(0., 0.), (1., 0.), (1., 1.), (0., 1.)] {
            let (x, n) = at(u, v);
            front &= (x - origin).dot(n) > 0.0;
        }
        return Ok((v6, area, front));
    }
    let mut v6 = 0.0;
    let mut area = 0.0;
    let mut front = true;
    for k in 0..face.len() {
        let (a, b) = (points[face[k]], points[face[(k + 1) % face.len()]]);
        let cross = (a - centre).cross(b - centre);
        let signed = cross.dot(centre - origin);
        front &= signed > 0.0 && cross.dot(newell) > 0.0;
        v6 += signed;
        area += 0.5 * cross.length();
    }
    Ok((v6, area, front))
}

/// Integrate the oriented shell from its certified interior point. The expected
/// volume is computed independently from the source tetrahedron partitions.
/// Every edge used twice in opposite directions, Euler characteristic 2.
fn shell_closed(poly: &Polyhedron) -> Result<(), &'static str> {
    let mut edges = HashMap::<(usize, usize), (usize, i32)>::new();
    for face in &poly.faces {
        for (&a, &b) in face.iter().zip(face.iter().cycle().skip(1)) {
            let e = edges.entry((a.min(b), a.max(b))).or_default();
            e.0 += 1;
            e.1 += if a < b { 1 } else { -1 };
        }
    }
    if edges.values().any(|&use_count| use_count != (2, 0))
        || poly.nodes.len() + poly.faces.len() != edges.len() + 2
    {
        return Err("open or nonmanifold polyhedron shell");
    }
    Ok(())
}

fn poly_volume(poly: &Polyhedron, points: &[Vec3]) -> Result<(f64, f64), &'static str> {
    let mut volume = 0.0;
    let mut area = 0.0;
    for face in &poly.faces {
        let (v6, a, front) = face_integral(points, face, poly.interior)?;
        if !v6.is_finite() || !front {
            return Err("inverted face or invalid polyhedron kernel");
        }
        volume += v6 / 6.0;
        area += a;
    }
    shell_closed(poly)?;
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

/// Area of a boundary face from its vector area. For a planar face this is
/// the area also when the face is not convex, where a sum of unsigned fan
/// triangles counts a notch twice.
fn boundary_area(p: &[Vec3], nodes: &[usize]) -> f64 {
    let a = p[nodes[0]];
    nodes[1..]
        .windows(2)
        .fold(Vec3::ZERO, |s, w| s + (p[w[0]] - a).cross(p[w[1]] - a))
        .length()
        * 0.5
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
        let area = boundary_area(p, &face.nodes);
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
    let mut pending = HashMap::<Vec<usize>, (Vec<usize>, usize)>::new();
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
            Element::Polyhedron(ref poly) if poly.prism.is_some() => {
                let mut value = 0.0;
                for w in poly.prism.iter().flatten() {
                    let (jac, v, _) = wedge(w.map(|v| p[v]));
                    if !jac.is_finite() || jac <= 0. {
                        return Err(format!(
                            "layer collision or inversion in layer {}, cell {id}, Jacobian {jac}",
                            cell.layer
                        ));
                    }
                    minimum_j = minimum_j.min(jac);
                    value += v;
                }
                shell_closed(poly).map_err(|e| format!("cell {id}: {e}"))?;
                if !value.is_finite() || (value / poly.expected_volume - 1.0).abs() > 1e-8 {
                    return Err(format!(
                        "cell {id}: layer piece volume differs from its shell"
                    ));
                }
                value
            }
            Element::Polyhedron(ref poly) => {
                let (volume, sphericity) = poly_volume(poly, p)
                    .map_err(|e| format!("cell {id} (layer {}): {e}", cell.layer))?;
                if cell.layer == 0 {
                    sphericities.push(sphericity);
                }
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

    /// Keep rule: a long thin planar face far from the origin passes the
    /// planarity check. The ring is a core face of the Manifold as built
    /// with conforming wall seeds (2.6 long, 0.00007 wide, 263 from the
    /// origin): two wall cells that barely touch over a wall vertex. Its
    /// vertices are on one plane to 1e-13, but a Newell normal from
    /// absolute coordinates is tilted by 1e-10 over the length and the
    /// offset test rejected the whole mesh. Taking the normal from a plain
    /// Newell sum again, or dropping the least-squares plane, breaks this;
    /// the fixtures of the fast run have no such face, and the Manifold
    /// build is a full-run test.
    #[test]
    fn thin_face_far_from_the_origin_is_planar() {
        let points = [
            [6.414802894970417, 262.9597716879253, 18.751259654212134],
            [6.859459004309515, 262.6533452053385, 18.034802806686436],
            [6.983285982556163, 262.5680121477317, 17.835285281547556],
            [6.997576348350637, 262.5581642080918, 17.812259778890425],
            [7.023771642070697, 262.54011220892545, 17.770052335510094],
            [7.023766730368811, 262.54005614619194, 17.77009088156666],
            [6.997572759971385, 262.5581082208852, 17.812295683907582],
            [6.983282513930192, 262.56795661671845, 17.8353207160207],
            [6.859456505831417, 262.6532936738063, 18.034834272658962],
            [6.414835399147537, 262.9597127964614, 18.751226084952386],
            [5.4468412280837315, 263.62682445664535, 20.310898102481982],
        ]
        .map(|[x, y, z]| Vec3::new(x, y, z));
        let ring: Vec<usize> = (0..points.len()).collect();
        let (newell, centre) = crate::duct::convex::ring_plane(&points, &ring);
        let unit = crate::duct::convex::fitted_normal(&points, &ring, newell, centre);
        let flatness = points
            .iter()
            .map(|&p| (p - centre).dot(unit).abs())
            .fold(0.0, f64::max);
        assert!(
            flatness < 1e-12,
            "offset from the fitted plane {flatness:e}"
        );
        // seen from a point off the plane the face integrates as planar
        let origin = centre - unit;
        let result = face_integral(&points, &ring, origin);
        let (v6, area, _) = result.expect("the face is planar");
        assert!(area > 5e-5 && area < 2e-4, "area {area:e}");
        assert!((v6.abs() / 6.0 - area / 3.0).abs() < 1e-9 * area);
    }

    /// Keep rule: the fitted plane does not hide a face that is bent. A
    /// pentagon with one corner lifted by a millionth of its size is not
    /// planar, with or without the fit; a fit that absorbed the lift (a
    /// wrong sign or axis in the tilt) would accept warped polyhedron faces
    /// everywhere, and no mesh test builds one.
    #[test]
    fn bent_face_is_still_rejected() {
        let mut points = [[0., 0.], [1., 0.], [1.4, 0.8], [0.5, 1.5], [-0.4, 0.8]]
            .map(|[x, y]| Vec3::new(100.0 + x, 200.0 + y, 50.0));
        points[2].z += 1e-6;
        let ring: Vec<usize> = (0..points.len()).collect();
        let origin = Vec3::new(100.5, 200.7, 49.0);
        assert_eq!(
            face_integral(&points, &ring, origin).err(),
            Some("nonplanar polyhedron face")
        );
    }

    /// Invariant: the area of a notched planar boundary face is its true
    /// area. Summing unsigned fan triangles from the first corner (the code
    /// before) counts the notch twice, and the patch area check then fails
    /// a valid mesh; the fixtures only have convex boundary faces.
    #[test]
    fn notched_boundary_face_has_its_true_area() {
        let v = |x: f64, y: f64| Vec3::new(x, y, 0.0);
        let points = [v(0., 0.), v(4., 0.), v(4., 4.), v(2., 1.), v(0., 4.)];
        // square of 16 without the triangle (4,4) (2,1) (0,4) of area 6
        let area = boundary_area(&points, &[0, 1, 2, 3, 4]);
        assert!((area - 10.0).abs() < 1e-12);
    }
}
