//! Vertex-centred median dual of a volume-filling tetrahedral scaffold.
use super::*;
use mesh_core::geometry::PatchId;
use mesh_geometry::SurfaceTriangle;
use std::collections::HashMap;

pub struct Scaffold {
    pub points: Vec<Vec3>,
    pub boundary: Vec<([usize; 3], usize)>,
    pub tetrahedra: Vec<[usize; 4]>,
}

#[derive(Default)]
struct Subdivision {
    edges: HashMap<[usize; 2], usize>,
    triangles: HashMap<[usize; 3], usize>,
}

fn sorted<const N: usize>(mut ids: [usize; N]) -> [usize; N] {
    ids.sort_unstable();
    ids
}

fn center(points: &mut Vec<Vec3>, ids: &[usize]) -> usize {
    let i = points.len();
    points.push(mean(points, ids));
    i
}

impl Subdivision {
    fn mid(&mut self, p: &mut Vec<Vec3>, a: usize, b: usize) -> usize {
        let edge = sorted([a, b]);
        *self.edges.entry(edge).or_insert_with(|| center(p, &edge))
    }
    fn face(&mut self, p: &mut Vec<Vec3>, ids: [usize; 3]) -> usize {
        let key = sorted(ids);
        *self.triangles.entry(key).or_insert_with(|| center(p, &key))
    }
    fn split(&mut self, p: &mut Vec<Vec3>, t: [usize; 3]) -> [[usize; 3]; 6] {
        let f = self.face(p, t);
        let m = [
            self.mid(p, t[0], t[1]),
            self.mid(p, t[1], t[2]),
            self.mid(p, t[2], t[0]),
        ];
        std::array::from_fn(|i| {
            let j = i / 2;
            if i % 2 == 0 {
                [t[j], m[j], f]
            } else {
                [t[j], f, m[(j + 2) % 3]]
            }
        })
    }
}

// The six subwedges exactly partition the original linear wedge, including
// its bilinear side surfaces. Their end triangles match the core dual pieces.
fn split_layers(front: &mut front::Front, cuts: &mut Subdivision) -> Result<(), String> {
    let mut patches = HashMap::new();
    let p = &mut front.points;
    for f in front.boundary.drain(..) {
        let mut add = |nodes: Vec<usize>| {
            let mut key = nodes;
            key.sort_unstable();
            patches.insert(key, f.patch);
        };
        match f.nodes[..] {
            [a, b, c] => {
                for t in cuts.split(p, [a, b, c]) {
                    add(t.to_vec())
                }
            }
            [a, b, c, d] => {
                let lo = cuts.mid(p, a, b);
                let hi = cuts.mid(p, c, d);
                add(vec![a, lo, hi, d]);
                add(vec![lo, b, c, hi]);
            }
            _ => return Err("unexpected layer boundary face".into()),
        }
    }
    let coarse = std::mem::take(&mut front.cells);
    for c in coarse {
        let Element::Wedge([a, b, c0, d, e, f]) = c.element else {
            return Err("expected layer prism".into());
        };
        let first = cuts.split(p, [a, b, c0]);
        let last = cuts.split(p, [d, e, f]);
        for ([a, b, c0], [d, e, f]) in first.into_iter().zip(last) {
            let element = Element::Wedge([a, b, c0, d, e, f]);
            for nodes in element.faces() {
                let mut key = nodes.clone();
                key.sort_unstable();
                if let Some(patch) = patches.remove(&key) {
                    front.boundary.push(BoundaryFace {
                        nodes,
                        patch,
                        owner: front.cells.len(),
                    });
                }
            }
            front.cells.push(Cell {
                element,
                layer: c.layer,
                region: 0,
            });
        }
    }
    if !patches.is_empty() {
        return Err("missing subdivided prism boundary".into());
    }
    Ok(())
}

#[derive(Default)]
struct Dual {
    faces: Vec<Vec<usize>>,
    volume: f64,
    // Sum of inward unit normals to find a strict witness at the boundary.
    inward: Vec3,
    boundary_faces: Vec<[usize; 3]>,
}

/// Extend the prescribed wall displacement through the scaffold. Weighting by
/// inverse squared edge length lets close neighbours track the same wall motion.
fn offset_interior(
    scaffold: &Scaffold,
    surface_ids: &BTreeMap<usize, usize>,
    front: &front::Front,
) -> Result<Vec<Vec3>, String> {
    let mut neighbors = vec![Vec::new(); scaffold.points.len()];
    for t in &scaffold.tetrahedra {
        for i in 0..4 {
            for j in i + 1..4 {
                neighbors[t[i]].push(t[j]);
                neighbors[t[j]].push(t[i]);
            }
        }
    }
    let mut displacement = vec![Vec3::ZERO; scaffold.points.len()];
    for (&id, &v) in surface_ids {
        displacement[id] = front.points[front.inner[v]] - scaffold.points[id];
    }
    let mut stencils = Vec::new();
    for (i, adj) in neighbors.iter_mut().enumerate() {
        if surface_ids.contains_key(&i) {
            continue;
        }
        adj.sort_unstable();
        adj.dedup();
        let mut row = Vec::new();
        let mut total = 0.;
        for &j in adj.iter() {
            let d = scaffold.points[i] - scaffold.points[j];
            let weight = 1. / d.dot(d);
            row.push((j, weight));
            total += weight;
        }
        if !total.is_finite() || total <= 0. {
            return Err("degenerate scaffold node neighbourhood".into());
        }
        for (_, w) in &mut row {
            *w /= total;
        }
        stencils.push((i, row));
    }
    for _ in 0..100 {
        let mut step = 0.0_f64;
        for (i, row) in &stencils {
            let d = row
                .iter()
                .fold(Vec3::ZERO, |a, &(j, w)| a + displacement[j] * w);
            step = step.max((displacement[*i] - d).max_abs_component());
            displacement[*i] = d;
        }
        if step <= f64::EPSILON {
            break;
        }
    }
    Ok(displacement)
}

// Interior scaffold vertices locate the dual cells, but only their surrounding
// midpoints and centroids belong to the final mesh connectivity.
fn compact_points(points: &mut Vec<Vec3>, cells: &mut [Cell], boundary: &mut [BoundaryFace]) {
    let mut ids = vec![usize::MAX; points.len()];
    for cell in cells.iter() {
        for &i in cell.element.nodes() {
            ids[i] = 0;
        }
    }
    let mut next = 0;
    for id in &mut ids {
        if *id == 0 {
            *id = next;
            next += 1;
        }
    }
    let mut cursor = 0;
    points.retain(|_| {
        let keep = ids[cursor] != usize::MAX;
        cursor += 1;
        keep
    });
    for cell in cells {
        match &mut cell.element {
            Element::Polyhedron(p) => {
                for n in &mut p.nodes {
                    *n = ids[*n];
                }
                for n in p.faces.iter_mut().flatten() {
                    *n = ids[*n];
                }
                p.cone_apex = p.cone_apex.map(|v| ids[v]);
            }
            Element::Wedge(n) => *n = n.map(|v| ids[v]),
            Element::Tet(n) => *n = n.map(|v| ids[v]),
        }
    }
    for n in boundary.iter_mut().flat_map(|f| &mut f.nodes) {
        *n = ids[*n];
    }
}

pub fn build_volume_dual(
    scaffold: Scaffold,
    steps: &[f64],
) -> Result<(DuctMesh, ClosedSurface), String> {
    if scaffold
        .boundary
        .iter()
        .any(|&(t, tag)| !(1..=3).contains(&tag) || t.iter().any(|&n| n >= scaffold.points.len()))
    {
        return Err("invalid scaffold boundary node or patch".into());
    }
    for t in &scaffold.tetrahedra {
        if t.iter().any(|&n| n >= scaffold.points.len()) {
            return Err("invalid scaffold tetrahedron node".into());
        }
        let [a, b, c, d] = t.map(|n| scaffold.points[n]);
        let det = (b - a).cross(c - a).dot(d - a);
        if !det.is_finite() || det <= 0. {
            return Err("nonpositive source scaffold tetrahedron".into());
        }
    }
    let mut surface_ids = BTreeMap::new();
    let mut surface_points = Vec::new();
    let mut triangles = Vec::new();
    let mut patches = Vec::new();
    for &(t, tag) in &scaffold.boundary {
        let vertices = t.map(|n| {
            *surface_ids.entry(n).or_insert_with(|| {
                let id = surface_points.len();
                surface_points.push(scaffold.points[n]);
                id
            })
        });
        triangles.push(SurfaceTriangle {
            vertices,
            patch: PatchId(tag - 1),
        });
        patches.push(tag);
    }
    let surface = ClosedSurface::new(
        surface_points,
        triangles,
        vec!["inlet".into(), "outlet".into(), "walls".into()],
    )
    .map_err(|e| e.to_string())?;
    let ts: Vec<_> = surface
        .facets()
        .iter()
        .map(|f| f.triangle.vertices)
        .collect();
    intersection::embedded(
        surface.vertices(),
        &ts,
        surface.bounds().extent().max_abs_component(),
    )?;
    let mut rims: [Vec<Edge>; 2] = Default::default();
    let mut port_edges = BTreeMap::new();
    for (&t, &tag) in ts.iter().zip(&patches) {
        if tag <= 2 {
            for (u, v) in edges(t) {
                if port_edges.remove(&(v, u)).is_none() {
                    port_edges.insert((u, v), tag);
                }
            }
        }
    }
    for ((a, b), tag) in port_edges {
        rims[tag - 1].push((b, a));
    }
    let mut front = front::advance(&surface, &patches, &vec![0; ts.len()], &rims, steps)?;
    let mut cuts = Subdivision::default();
    split_layers(&mut front, &mut cuts)?;
    let displacement = offset_interior(&scaffold, &surface_ids, &front)?;
    let mut mapped = Vec::with_capacity(scaffold.points.len());
    for (id, &pos) in scaffold.points.iter().enumerate() {
        mapped.push(if let Some(&v) = surface_ids.get(&id) {
            front.inner[v]
        } else {
            let i = front.points.len();
            front.points.push(pos + displacement[id]);
            i
        });
    }
    let p = &mut front.points;
    let mut dual: Vec<Dual> = (0..mapped.len()).map(|_| Dual::default()).collect();
    let owner_offset = front.cells.len();
    for &(t, patch) in &scaffold.boundary {
        let fine = cuts.split(p, t.map(|n| mapped[n]));
        for (j, nodes) in fine.into_iter().enumerate() {
            let id = t[j / 2];
            let [a, b, c] = nodes.map(|n| p[n]);
            let normal = (b - a).cross(c - a);
            dual[id].inward = dual[id].inward - normal / normal.length();
            dual[id].boundary_faces.push(nodes);
            dual[id].faces.push(nodes.to_vec());
            if patch <= 2 || steps.is_empty() {
                front.boundary.push(BoundaryFace {
                    nodes: nodes.to_vec(),
                    patch,
                    owner: owner_offset + id,
                });
            }
        }
    }
    for (index, t) in scaffold.tetrahedra.iter().enumerate() {
        let q = t.map(|n| mapped[n]);
        let [a, b, c, d] = q.map(|n| p[n]);
        let six = (b - a).cross(c - a).dot(d - a);
        if !six.is_finite() || six <= 0. {
            return Err(format!(
                "inverted scaffold tetrahedron {index} after wall offset, determinant {six}"
            ));
        }
        let centroid = center(p, &q);
        for &n in t {
            dual[n].volume += six / 24.;
        }
        // Each four-point face lies in the affine plane where the two
        // endpoint barycentric coordinates are equal.
        for [i, j, k, l] in [
            [0, 1, 2, 3],
            [0, 2, 3, 1],
            [0, 3, 1, 2],
            [1, 2, 0, 3],
            [1, 3, 2, 0],
            [2, 3, 0, 1],
        ] {
            let edge = cuts.mid(p, q[i], q[j]);
            let f = cuts.face(p, [q[i], q[j], q[k]]);
            let g = cuts.face(p, [q[i], q[j], q[l]]);
            let face = [edge, f, centroid, g];
            dual[t[i]].faces.push(face.to_vec());
            dual[t[j]].faces.push(face.into_iter().rev().collect());
        }
    }
    for (id, cell) in dual.into_iter().enumerate() {
        if cell.volume <= 0. {
            return Err(format!("unused volume node {id}"));
        }
        let site = p[mapped[id]];
        let mut witness = site;
        if !cell.boundary_faces.is_empty() {
            // A sufficiently short inward step preserves strict visibility of
            // every internal dual face, even in a slender boundary tetrahedron.
            let mut direction = cell.inward;
            for _ in 0..32 {
                let mut valid = true;
                for face in &cell.boundary_faces {
                    let [a, b, c] = face.map(|i| p[i]);
                    let normal = (b - a).cross(c - a);
                    let n = normal / normal.length();
                    if direction.dot(n) > -0.01 {
                        direction = direction - n * (direction.dot(n) + 0.02);
                        valid = false;
                    }
                }
                if valid {
                    break;
                }
            }
            let mut distance = cell.volume.cbrt() * 0.1;
            let direction = direction / direction.length();
            for face in &cell.faces {
                let a = p[face[0]];
                let n = (p[face[1]] - a).cross(p[face[2]] - a);
                if n.dot(direction) > 0. {
                    distance = distance.min(n.dot(a - site) / n.dot(direction) * 0.5);
                }
            }
            witness = site + direction * distance;
        }
        front.cells.push(Cell {
            element: Element::Polyhedron(Box::new(Polyhedron::new(
                cell.faces,
                witness,
                None,
                cell.volume,
            ))),
            layer: 0,
            region: 0,
        });
    }
    compact_points(p, &mut front.cells, &mut front.boundary);
    let measurements = verify::check_with_interface_count(
        &surface,
        &patches,
        p,
        &front.cells,
        &front.boundary,
        steps.len(),
        patches.iter().filter(|&&p| p == 3).count() * 6,
    )?;
    Ok((
        DuctMesh {
            points: front.points,
            cells: front.cells,
            boundary: front.boundary,
            measurements,
        },
        surface,
    ))
}
