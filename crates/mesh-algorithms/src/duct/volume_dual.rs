//! Vertex-centred median dual of a volume-filling tetrahedral scaffold.
use super::*;
use mesh_core::geometry::PatchId;
use mesh_geometry::SurfaceTriangle;
use std::collections::HashMap;

#[derive(Clone, Debug, PartialEq)]
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

/// The part of one scaffold tetrahedron nearest one of its vertices: the
/// vertex, three edge midpoints, three face centroids and the centroid.
#[derive(Default)]
struct Corner {
    /// The three median-dual quadrilaterals towards the other vertices.
    faces: Vec<Vec<usize>>,
    /// The tetrahedron faces through the vertex (sorted scaffold ids) with
    /// the quadrilateral of the corner piece on each, oriented outward.
    sides: Vec<([usize; 3], Vec<usize>)>,
    volume: f64,
    /// Mean of the eight corner points, used only to group pieces.
    centre: Vec3,
}

#[derive(Default)]
struct Dual {
    corners: Vec<Corner>,
    /// Fine boundary triangles, their patch and the scaffold face they lie on.
    boundary: Vec<([usize; 3], usize, [usize; 3])>,
}

/// The faces of the cell made of the corner pieces in `group`: their dual
/// quadrilaterals, the boundary triangles on their tetrahedra, and the corner
/// quadrilaterals on tetrahedron faces shared with a piece outside the group.
/// Also returns the faces through the vertex and the boundary triangles used.
fn part_faces(cell: &Dual, group: &[usize]) -> (Vec<Vec<usize>>, Vec<[usize; 3]>, Vec<usize>) {
    let owns = |key: &[usize; 3]| {
        group
            .iter()
            .any(|&c| cell.corners[c].sides.iter().any(|(k, _)| k == key))
    };
    let mut faces = Vec::new();
    let mut planes = Vec::new();
    for &c in group {
        let corner = &cell.corners[c];
        faces.extend(corner.faces.iter().cloned());
        for (key, quad) in &corner.sides {
            let inside = group
                .iter()
                .any(|&o| o != c && cell.corners[o].sides.iter().any(|(k, _)| k == key));
            let on_boundary = cell.boundary.iter().any(|(_, _, k)| k == key);
            if !inside && !on_boundary {
                faces.push(quad.clone());
                planes.push([quad[0], quad[1], quad[3]]);
            }
        }
    }
    let mut used = Vec::new();
    for (i, (nodes, _, key)) in cell.boundary.iter().enumerate() {
        if owns(key) {
            faces.push(nodes.to_vec());
            planes.push(*nodes);
            used.push(i);
        }
    }
    (faces, planes, used)
}

/// Centroid of a closed shell, integrated over a fan from each face's vertex
/// mean as FoehnCraft and the hand-off checks do.
fn shell_centroid(points: &[Vec3], faces: &[Vec<usize>]) -> Vec3 {
    let Some(&first) = faces.first().and_then(|f| f.first()) else {
        return Vec3::ZERO;
    };
    let origin = points[first];
    let mut v6 = 0.;
    let mut moment = Vec3::ZERO;
    for ring in faces {
        let c = mean(points, ring) - origin;
        for (k, &a) in ring.iter().enumerate() {
            let b = ring[(k + 1) % ring.len()];
            let (pa, pb) = (points[a] - origin, points[b] - origin);
            let t6 = c.dot(pa.cross(pb));
            v6 += t6;
            moment = moment + (c + pa + pb) * (t6 / 4.);
        }
    }
    origin + moment * (1. / v6)
}

/// Faces between a part and a wall prism whose area vector does not point
/// from the part's centroid towards the prism's (S·d ≤ 0).
fn wall_inversions(
    p: &[Vec3],
    cell: &Dual,
    group: &[usize],
    wall: &dyn Fn(&[usize; 3]) -> Option<Vec3>,
) -> usize {
    let (faces, _, used) = part_faces(cell, group);
    let centre = shell_centroid(p, &faces);
    used.iter()
        .filter(|&&i| {
            let nodes = cell.boundary[i].0;
            wall(&nodes).is_some_and(|prism| {
                let [a, b, c] = nodes.map(|n| p[n]);
                (b - a).cross(c - a).dot(prism - centre) <= 0.
            })
        })
        .count()
}

/// Feature angle that separates the wall planes of a split dual cell.
const SPLIT_FEATURE_COS: f64 = 0.866_025_403_784_438_6; // cos 30°

/// A median-dual cell around a vertex on a reentrant wall edge is L-shaped:
/// its centroid lies behind the wall planes next to the vertex, so the faces
/// to the wall prisms there point the wrong way (S·d ≤ 0, FoehnCraft calls
/// them inverted). Split such a cell along scaffold faces through the vertex,
/// one part per wall plane: each corner piece joins the plane it stands
/// farthest in front of. The split is kept only where it removes inversions
/// and every part is a connected, manifold shell.
fn split_reentrant(
    p: &[Vec3],
    site: usize,
    cell: &Dual,
    wall: &dyn Fn(&[usize; 3]) -> Option<Vec3>,
) -> Vec<Vec<usize>> {
    let all: Vec<usize> = (0..cell.corners.len()).collect();
    let before = wall_inversions(p, cell, &all, wall);
    if before == 0 {
        return vec![all];
    }
    let mut planes: Vec<Vec3> = Vec::new();
    for (nodes, _, _) in &cell.boundary {
        let [a, b, c] = nodes.map(|n| p[n]);
        let normal = (b - a).cross(c - a);
        let inward = normal * (-1. / normal.length());
        if !planes.iter().any(|n| n.dot(inward) > SPLIT_FEATURE_COS) {
            planes.push(inward);
        }
    }
    if planes.len() < 2 {
        return vec![all];
    }
    let mut groups = vec![Vec::new(); planes.len()];
    for (c, corner) in cell.corners.iter().enumerate() {
        let r = corner.centre - p[site];
        let mut best = 0;
        for (k, n) in planes.iter().enumerate() {
            if r.dot(*n) > r.dot(planes[best]) {
                best = k;
            }
        }
        groups[best].push(c);
    }
    // Pieces join through shared tetrahedron faces; each connected set is
    // its own part.
    let mut parts = Vec::new();
    for group in groups {
        let mut left = group;
        while let Some(seed) = left.pop() {
            let mut part = vec![seed];
            let mut i = 0;
            while i < part.len() {
                let keys: Vec<[usize; 3]> = cell.corners[part[i]]
                    .sides
                    .iter()
                    .map(|(k, _)| *k)
                    .collect();
                let (join, rest): (Vec<usize>, Vec<usize>) = left
                    .iter()
                    .partition(|&&o| cell.corners[o].sides.iter().any(|(k, _)| keys.contains(k)));
                part.extend(join);
                left = rest;
                i += 1;
            }
            part.sort_unstable();
            parts.push(part);
        }
    }
    parts.sort();
    if parts.len() < 2 {
        return vec![all];
    }
    for part in &parts {
        let (faces, _, _) = part_faces(cell, part);
        let mut edges = HashMap::<(usize, usize), i32>::new();
        for f in &faces {
            for (k, &a) in f.iter().enumerate() {
                let b = f[(k + 1) % f.len()];
                *edges.entry((a.min(b), a.max(b))).or_default() += if a < b { 1 } else { -1 };
            }
        }
        let mut uses = HashMap::<(usize, usize), usize>::new();
        for f in &faces {
            for (k, &a) in f.iter().enumerate() {
                let b = f[(k + 1) % f.len()];
                *uses.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        if edges.values().any(|&s| s != 0) || uses.values().any(|&u| u != 2) {
            return vec![all];
        }
    }
    let after: usize = parts
        .iter()
        .map(|part| wall_inversions(p, cell, part, wall))
        .sum();
    if after < before { parts } else { vec![all] }
}

/// A point strictly inside a boundary cell that sees every face: a short
/// step from the vertex into the cell, in front of every face through the
/// vertex (`planes`, oriented outward), and halfway to the nearest other face.
fn strict_witness(
    p: &[Vec3],
    site: Vec3,
    faces: &[Vec<usize>],
    planes: &[[usize; 3]],
    volume: f64,
) -> Vec3 {
    let unit = |[a, b, c]: [usize; 3]| {
        let n = (p[b] - p[a]).cross(p[c] - p[a]);
        n / n.length()
    };
    let mut direction = planes.iter().fold(Vec3::ZERO, |s, &f| s - unit(f));
    // A sufficiently short inward step preserves strict visibility of every
    // internal dual face, even in a slender boundary tetrahedron.
    for _ in 0..32 {
        let mut valid = true;
        for &face in planes {
            let n = unit(face);
            if direction.dot(n) > -0.01 {
                direction = direction - n * (direction.dot(n) + 0.02);
                valid = false;
            }
        }
        if valid {
            break;
        }
    }
    let mut distance = volume.cbrt() * 0.1;
    let direction = direction / direction.length();
    for face in faces {
        let a = p[face[0]];
        let n = (p[face[1]] - a).cross(p[face[2]] - a);
        if n.dot(direction) > 0. {
            distance = distance.min(n.dot(a - site) / n.dot(direction) * 0.5);
        }
    }
    site + direction * distance
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
pub(super) fn compact_points(
    points: &mut Vec<Vec3>,
    cells: &mut [Cell],
    boundary: &mut [BoundaryFace],
) {
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
                for w in p.prism.iter_mut().flatten() {
                    *w = w.map(|v| ids[v]);
                }
            }
            Element::Wedge(n) => *n = n.map(|v| ids[v]),
            Element::Tet(n) => *n = n.map(|v| ids[v]),
        }
    }
    for n in boundary.iter_mut().flat_map(|f| &mut f.nodes) {
        *n = ids[*n];
    }
}

pub(super) struct Prepared {
    pub surface: ClosedSurface,
    /// Patch tag (1 inlet, 2 outlet, 3 walls) of every scaffold boundary
    /// triangle, in scaffold order; the surface facets keep that order.
    pub patches: Vec<usize>,
    pub front: front::Front,
    /// Scaffold node id to surface vertex id for boundary nodes.
    pub surface_ids: BTreeMap<usize, usize>,
}

/// Check the scaffold, build its boundary surface and the wall layers.
pub(super) fn prepare(scaffold: &Scaffold, steps: &[f64]) -> Result<Prepared, String> {
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
    let front = front::advance(&surface, &patches, &vec![0; ts.len()], &rims, steps)?;
    Ok(Prepared {
        surface,
        patches,
        front,
        surface_ids,
    })
}

/// Move the interior scaffold nodes with the wall offset and return the
/// front point id of every scaffold node.
pub(super) fn map_scaffold(
    scaffold: &Scaffold,
    surface_ids: &BTreeMap<usize, usize>,
    front: &mut front::Front,
) -> Result<Vec<usize>, String> {
    let displacement = offset_interior(scaffold, surface_ids, front)?;
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
    Ok(mapped)
}

pub fn build_volume_dual(
    scaffold: Scaffold,
    steps: &[f64],
) -> Result<(DuctMesh, ClosedSurface), String> {
    let Prepared {
        surface,
        patches,
        mut front,
        surface_ids,
    } = prepare(&scaffold, steps)?;
    let mut cuts = Subdivision::default();
    split_layers(&mut front, &mut cuts)?;
    let mapped = map_scaffold(&scaffold, &surface_ids, &mut front)?;
    let p = &mut front.points;
    let mut dual: Vec<Dual> = (0..mapped.len()).map(|_| Dual::default()).collect();
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
        let corners: [usize; 4] = std::array::from_fn(|i| {
            let corner = dual[t[i]].corners.len();
            dual[t[i]].corners.push(Corner {
                volume: six / 24.,
                ..Corner::default()
            });
            corner
        });
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
            dual[t[i]].corners[corners[i]].faces.push(face.to_vec());
            dual[t[j]].corners[corners[j]]
                .faces
                .push(face.into_iter().rev().collect());
        }
        // The three tetrahedron faces through each vertex bound its corner
        // piece. They become cell faces only where a dual cell is split.
        for i in 0..4 {
            let v = q[i];
            let others: Vec<usize> = (0..4).filter(|&k| k != i).collect();
            let mut hint = p[v] + p[centroid];
            for (x, y, z) in [
                (others[0], others[1], others[2]),
                (others[1], others[2], others[0]),
                (others[2], others[0], others[1]),
            ] {
                let (mut ma, mut mb) = (cuts.mid(p, v, q[x]), cuts.mid(p, v, q[y]));
                let f = cuts.face(p, [v, q[x], q[y]]);
                if (p[ma] - p[v]).cross(p[mb] - p[v]).dot(p[q[z]] - p[v]) > 0. {
                    std::mem::swap(&mut ma, &mut mb);
                }
                hint = hint + p[ma] + p[f];
                dual[t[i]].corners[corners[i]]
                    .sides
                    .push((sorted([t[i], t[x], t[y]]), vec![v, ma, f, mb]));
            }
            dual[t[i]].corners[corners[i]].centre = hint / 8.;
        }
    }
    for &(t, patch) in &scaffold.boundary {
        let fine = cuts.split(p, t.map(|n| mapped[n]));
        for (j, nodes) in fine.into_iter().enumerate() {
            dual[t[j / 2]].boundary.push((nodes, patch, sorted(t)));
        }
    }
    // Wall prisms of the innermost layer by their inner end triangle.
    let mut prisms = HashMap::new();
    for (index, cell) in front.cells.iter().enumerate() {
        for face in cell.element.faces() {
            if let [a, b, c] = face[..] {
                prisms.insert(sorted([a, b, c]), index);
            }
        }
    }
    let prism_centres: Vec<Vec3> = front
        .cells
        .iter()
        .map(|c| shell_centroid(p, &c.element.faces()))
        .collect();
    for (id, cell) in dual.into_iter().enumerate() {
        let volume: f64 = cell.corners.iter().map(|c| c.volume).sum();
        if volume <= 0. {
            return Err(format!("unused volume node {id}"));
        }
        let site = mapped[id];
        let owner_offset = front.cells.len();
        let wall = |nodes: &[usize; 3]| prisms.get(&sorted(*nodes)).map(|&c| prism_centres[c]);
        let parts = split_reentrant(p, site, &cell, &wall);
        for (part, group) in parts.iter().enumerate() {
            let (faces, planes, used) = part_faces(&cell, group);
            let part_volume: f64 = group.iter().map(|&c| cell.corners[c].volume).sum();
            for i in used {
                let (nodes, patch, _) = cell.boundary[i];
                if patch <= 2 || steps.is_empty() {
                    front.boundary.push(BoundaryFace {
                        nodes: nodes.to_vec(),
                        patch,
                        owner: owner_offset + part,
                    });
                }
            }
            let witness = if planes.is_empty() {
                p[site]
            } else {
                strict_witness(p, p[site], &faces, &planes, part_volume)
            };
            front.cells.push(Cell {
                element: Element::Polyhedron(Box::new(Polyhedron::new(
                    faces,
                    witness,
                    None,
                    part_volume,
                ))),
                layer: 0,
                region: 0,
            });
        }
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
