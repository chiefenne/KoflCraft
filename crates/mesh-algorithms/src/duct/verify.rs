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
    /// The solver rules, counted beside the strict checks.
    pub solver_rules: SolverRules,
}

/// Cells counted against the validity rules of production finite-volume
/// codes (OpenFOAM's mesh description and checkMesh, the STAR-CCM+ and
/// Fluent cell checks; research/robust-polyhedral-meshing): a cell is
/// usable when its shell is closed, its volume is positive, every face
/// pyramid on the cell centre is positive and every face points away from
/// the cell (owner-to-neighbour vector, or centre-to-face vector at the
/// boundary, at less than 90° to the face normal). Faces need not be flat:
/// a face's area vector, centre and pyramids come from the triangle fan on
/// its vertex mean, as in those codes. Report only; the strict checks
/// below still decide.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SolverRules {
    pub cells: usize,
    /// Cells failing at least one rule.
    pub failing_cells: usize,
    /// Failing cells of the wall layers (layer > 0); the rest are core cells.
    pub failing_layer_cells: usize,
    /// Shell not closed: an edge not used by exactly two faces of the cell,
    /// or the face area vectors not summing to zero (1e-6 of their sum).
    pub open_cells: usize,
    pub nonpositive_volume: usize,
    /// Cells with a face pyramid of nonpositive volume.
    pub negative_pyramids: usize,
    /// Cells with a face at 90° or more (non-orthogonality).
    pub faces_turned_away: usize,
    /// Largest non-orthogonality over all faces, degrees.
    pub maximum_non_orthogonality: f64,
    /// Largest face warping, degrees: the angle between a triangle of a
    /// face's fan and the face's area vector (a graded quality measure, as
    /// the warpage angle of STAR-CCM+; no limit).
    pub maximum_warping: f64,
}

/// Area vector and centre of a face as finite-volume codes take them:
/// triangles fanned from the vertex mean, centres weighted by area.
fn fv_face(p: &[Vec3], ring: &[usize]) -> (Vec3, Vec3) {
    let mean = ring.iter().fold(Vec3::ZERO, |s, &v| s + p[v]) / ring.len() as f64;
    if ring.len() == 3 {
        let [a, b, c] = [p[ring[0]], p[ring[1]], p[ring[2]]];
        return ((b - a).cross(c - a) * 0.5, mean);
    }
    let mut area = Vec3::ZERO;
    let mut moment = Vec3::ZERO;
    let mut weight = 0.0;
    for (i, &v) in ring.iter().enumerate() {
        let (a, b) = (p[v], p[ring[(i + 1) % ring.len()]]);
        let n = (a - mean).cross(b - mean) * 0.5;
        area = area + n;
        let w = n.length();
        moment = moment + (a + b + mean) * (w / 3.0);
        weight += w;
    }
    let centre = if weight > 0.0 { moment / weight } else { mean };
    (area, centre)
}

/// Counts the cells of a mesh against [`SolverRules`].
#[cfg(test)]
fn solver_rules(p: &[Vec3], cells: &[Cell]) -> SolverRules {
    solver_rules_and_volumes(p, cells).0
}

/// Largest angle between a fan triangle of a face and its area vector.
fn warping(p: &[Vec3], ring: &[usize], area: Vec3) -> f64 {
    if ring.len() == 3 {
        return 0.0;
    }
    let mean = ring.iter().fold(Vec3::ZERO, |s, &v| s + p[v]) / ring.len() as f64;
    let mut worst = 0.0_f64;
    for (i, &v) in ring.iter().enumerate() {
        let n = (p[v] - mean).cross(p[ring[(i + 1) % ring.len()]] - mean);
        let cos = n.dot(area) / (n.length() * area.length());
        if cos.is_finite() {
            worst = worst.max(cos.clamp(-1.0, 1.0).acos().to_degrees());
        }
    }
    worst
}

/// A cell's shell as the solver sees it.
pub(super) struct Shell {
    pub open: bool,
    pub volume: f64,
    pub centre: Vec3,
    pub pyramids_ok: bool,
}

impl Shell {
    /// Closed, positive volume and every face pyramid on the centre positive
    /// (the solver rules that need no neighbour).
    pub fn valid(&self) -> bool {
        !self.open && self.volume > 0.0 && self.pyramids_ok
    }
}

/// The shell of a cell given by its outward face rings.
pub(super) fn cell_shell(p: &[Vec3], rings: &[Vec<usize>]) -> Shell {
    let geometry: Vec<(Vec3, Vec3)> = rings.iter().map(|r| fv_face(p, r)).collect();
    shell(rings, &geometry)
}

fn shell(rings: &[Vec<usize>], geometry: &[(Vec3, Vec3)]) -> Shell {
    // closed shell: every edge twice, area vectors summing to zero
    let mut edges: HashMap<(usize, usize), usize> = HashMap::new();
    for r in rings {
        for (i, &a) in r.iter().enumerate() {
            let b = r[(i + 1) % r.len()];
            *edges.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    let sum = geometry.iter().fold(Vec3::ZERO, |s, g| s + g.0);
    let total: f64 = geometry.iter().map(|g| g.0.length()).sum();
    let open = edges.values().any(|&n| n != 2)
        || sum.length().partial_cmp(&(1e-6 * total)) == Some(std::cmp::Ordering::Greater)
        || !total.is_finite();
    // volume and centre from pyramids on the mean of the face centres
    let estimate = geometry.iter().fold(Vec3::ZERO, |s, g| s + g.1) / geometry.len().max(1) as f64;
    let mut volume = 0.0_f64;
    let mut moment = Vec3::ZERO;
    for (area, centre) in geometry {
        let v = area.dot(*centre - estimate) / 3.0;
        volume += v;
        moment = moment + (*centre * 0.75 + estimate * 0.25) * v;
    }
    let centre = if volume.abs() > 0.0 && volume.is_finite() {
        moment / volume
    } else {
        estimate
    };
    let pyramids_ok = geometry.iter().all(|(area, c)| area.dot(*c - centre) > 0.0);
    Shell {
        open,
        volume,
        centre,
        pyramids_ok,
    }
}

/// The solver rules and each cell's volume as the solver takes it (the sum
/// of its face pyramids).
fn solver_rules_and_volumes(p: &[Vec3], cells: &[Cell]) -> (SolverRules, Vec<f64>) {
    let mut volumes = Vec::with_capacity(cells.len());
    let mut out = SolverRules {
        cells: cells.len(),
        ..SolverRules::default()
    };
    let mut centres = Vec::with_capacity(cells.len());
    let mut bad = vec![false; cells.len()];
    // face key -> (cell, area vector, face centre) of the first cell seen
    let mut faces: HashMap<Vec<usize>, (usize, Vec3, Vec3)> = HashMap::new();
    let mut shared: Vec<(usize, usize, Vec3)> = Vec::new();
    let mut open_faces: Vec<(usize, Vec3, Vec3)> = Vec::new();
    for (id, cell) in cells.iter().enumerate() {
        let rings = cell.element.faces();
        let geometry: Vec<(Vec3, Vec3)> = rings.iter().map(|r| fv_face(p, r)).collect();
        for (r, g) in rings.iter().zip(&geometry) {
            out.maximum_warping = out.maximum_warping.max(warping(p, r, g.0));
        }
        let Shell {
            open,
            volume,
            centre,
            pyramids_ok,
        } = shell(&rings, &geometry);
        out.open_cells += usize::from(open);
        let nonpositive = volume.is_nan() || volume <= 0.0;
        out.nonpositive_volume += usize::from(nonpositive);
        out.negative_pyramids += usize::from(!pyramids_ok);
        bad[id] = open || nonpositive || !pyramids_ok;
        centres.push(centre);
        volumes.push(volume);
        for (r, (area, c)) in rings.iter().zip(geometry) {
            match faces.remove(&key(r)) {
                Some((other, other_area, _)) => shared.push((other, id, other_area)),
                None => {
                    faces.insert(key(r), (id, area, c));
                }
            }
        }
    }
    for (_, (id, area, c)) in faces {
        open_faces.push((id, area, c));
    }
    let mut turned = vec![false; cells.len()];
    let mut worst = 0.0_f64;
    let mut check = |cells_of_face: &[usize], d: Vec3, area: Vec3| {
        let cos = d.dot(area) / (d.length() * area.length());
        let degrees = if cos.is_finite() {
            cos.clamp(-1.0, 1.0).acos().to_degrees()
        } else {
            180.0
        };
        worst = worst.max(degrees);
        if cos.is_nan() || cos <= 0.0 {
            for &c in cells_of_face {
                turned[c] = true;
            }
        }
    };
    for (owner, neighbour, area) in shared {
        check(
            &[owner, neighbour],
            centres[neighbour] - centres[owner],
            area,
        );
    }
    for (owner, area, c) in open_faces {
        check(&[owner], c - centres[owner], area);
    }
    out.maximum_non_orthogonality = worst;
    out.faces_turned_away = turned.iter().filter(|&&t| t).count();
    for ((b, t), cell) in bad.iter().zip(&turned).zip(cells) {
        if *b || *t {
            out.failing_cells += 1;
            out.failing_layer_cells += usize::from(cell.layer > 0);
        }
    }
    (out, volumes)
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
    let (solver_rules, solver_volumes) = solver_rules_and_volumes(p, cells);
    if crate::trace_enabled() {
        eprintln!(
            "kofl trace: solver rules: {} of {} cells fail ({} in layers; open {}, volume {}, pyramids {}, turned away {}), max non-orthogonality {:.2}°",
            solver_rules.failing_cells,
            solver_rules.cells,
            solver_rules.failing_layer_cells,
            solver_rules.open_cells,
            solver_rules.nonpositive_volume,
            solver_rules.negative_pyramids,
            solver_rules.faces_turned_away,
            solver_rules.maximum_non_orthogonality
        );
    }
    if solver_rules.failing_cells > 0 {
        return Err(format!(
            "{} of {} cells fail the solver rules ({} in layers; open {}, volume {}, pyramids {}, turned away {})",
            solver_rules.failing_cells,
            solver_rules.cells,
            solver_rules.failing_layer_cells,
            solver_rules.open_cells,
            solver_rules.nonpositive_volume,
            solver_rules.negative_pyramids,
            solver_rules.faces_turned_away
        ));
    }
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
        // The solver rules decided above; the element measures below are
        // quality reports only.
        match cell.element {
            Element::Polyhedron(ref poly) if poly.prism.is_some() => {
                for w in poly.prism.iter().flatten() {
                    let (jac, _, _) = wedge(w.map(|v| p[v]));
                    minimum_j = minimum_j.min(jac);
                }
            }
            Element::Polyhedron(ref poly) => {
                if cell.layer == 0
                    && let Ok((_, sphericity)) = poly_volume(poly, p)
                {
                    sphericities.push(sphericity);
                }
            }
            Element::Tet(n) => {
                let [o, a, b, c] = n.map(|v| p[v]);
                let determinant = det(b - a, c - a, a - o);
                let clearance = determinant / (b - a).cross(c - a).length();
                minimum_kernel_clearance = minimum_kernel_clearance.min(clearance);
                let squared_edges = [b - a, c - b, a - c, a - o, b - o, c - o]
                    .iter()
                    .map(|e| e.dot(*e))
                    .sum::<f64>();
                let mean_ratio = 12.0 * (determinant / 2.0).cbrt().powi(2) / squared_edges;
                minimum_tet_mean_ratio = minimum_tet_mean_ratio.min(mean_ratio);
                tetrahedra_below_mean_ratio_0_1 += usize::from(mean_ratio < 0.1);
            }
            Element::Wedge(n) => {
                let (jac, _, range) = wedge(n.map(|v| p[v]));
                minimum_j = minimum_j.min(jac);
                let r = &mut normal_thickness_ranges[cell.layer - 1];
                r[0] = r[0].min(range[0]);
                r[1] = r[1].max(range[1]);
            }
        }
        let value = solver_volumes[id];
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
        solver_rules,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cell(element: Element) -> Cell {
        Cell {
            element,
            layer: 0,
            region: 0,
        }
    }

    #[test]
    fn warped_faces_pass_the_solver_rules() {
        // a wedge whose top corner is lifted: its three side quads are not
        // flat, which the solver rules accept
        let p = [
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [1., 0., 1.3],
            [0., 1., 1.],
        ]
        .map(|[x, y, z]| Vec3::new(x, y, z));
        let r = solver_rules(&p, &[cell(Element::Wedge([0, 1, 2, 3, 4, 5]))]);
        assert_eq!(r.failing_cells, 0, "{r:?}");
        assert!(r.maximum_non_orthogonality < 90.0);
    }

    #[test]
    fn two_tetrahedra_pass_and_an_inverted_one_fails() {
        let p = [
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [0.6, 0.6, 0.6],
        ]
        .map(|[x, y, z]| Vec3::new(x, y, z));
        let good = [
            cell(Element::Tet([0, 1, 2, 3])),
            cell(Element::Tet([1, 2, 3, 4])),
        ];
        // the second tetrahedron must be oriented positively
        let r = solver_rules(&p, &good);
        assert_eq!((r.cells, r.failing_cells), (2, 0), "{r:?}");
        let inverted = [cell(Element::Tet([0, 2, 1, 3]))];
        let r = solver_rules(&p, &inverted);
        assert_eq!(r.failing_cells, 1);
        assert_eq!(r.nonpositive_volume, 1);
    }

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
