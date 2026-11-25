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
    /// The scaffold vertex across each of `faces`.
    towards: Vec<usize>,
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

/// Splits a part that breaks the solver rules (verify::Shell::valid: closed,
/// positive volume, every face pyramid on the centroid positive) into parts
/// that keep them. Each corner piece is convex (the intersection of the
/// half-spaces where its vertex's barycentric coordinate is the largest),
/// so a single piece always passes. Parts grow greedily from the largest
/// remaining piece over pieces that share a tetrahedron face with them,
/// as long as the union still passes: the local split-and-keep fallback
/// of production meshers (research/robust-polyhedral-meshing).
fn split_invalid(
    p: &[Vec3],
    cell: &Dual,
    part: Vec<usize>,
    earlier: &[Vec<usize>],
    split_quads: &HashMap<[usize; 4], usize>,
    limit: Option<f64>,
) -> Vec<Vec<usize>> {
    let mut parts: Vec<Vec<usize>> = earlier.to_vec();
    // The part as the solver sees it after the hand-off: the faces it shares
    // with one neighbour merged into one polygon. Neighbours: the unsplit
    // dual cell across (`towards`), a part split before (`split_quads`),
    // a part of this cell made before (`done`) and the rest of this cell.
    let valid = |group: &[usize], done: &[Vec<usize>]| {
        let mut groups: HashMap<(u8, usize), Vec<&Vec<usize>>> = HashMap::new();
        let mut single: Vec<Vec<usize>> = Vec::new();
        for &c in group {
            let corner = &cell.corners[c];
            for (f, &across) in corner.faces.iter().zip(&corner.towards) {
                let key = match split_quads.get(&quad_key(f)) {
                    Some(&owner) => (0, owner),
                    None => (1, across),
                };
                groups.entry(key).or_default().push(f);
            }
            for (key, quad) in &corner.sides {
                let shares = |o: usize| cell.corners[o].sides.iter().any(|(k, _)| k == key);
                if group.iter().any(|&o| o != c && shares(o))
                    || cell.boundary.iter().any(|(_, _, k)| k == key)
                {
                    continue;
                }
                let to = done.iter().position(|d| d.iter().any(|&o| shares(o)));
                groups
                    .entry(to.map_or((3, 0), |j| (2, j)))
                    .or_default()
                    .push(quad);
            }
        }
        for (nodes, _, key) in &cell.boundary {
            if group
                .iter()
                .any(|&c| cell.corners[c].sides.iter().any(|(k, _)| k == key))
            {
                single.push(nodes.to_vec());
            }
        }
        for rings in groups.values() {
            match verify::merge_rings(rings) {
                Some(ring) => single.push(ring),
                None => single.extend(rings.iter().map(|r| r.to_vec())),
            }
        }
        verify::cell_shell(p, &single).valid()
    };
    let touches = |a: usize, b: usize| {
        cell.corners[a]
            .sides
            .iter()
            .any(|(k, _)| cell.corners[b].sides.iter().any(|(o, _)| o == k))
    };
    // A solver takes one face between two cells, so the faces a part shares
    // with each neighbouring cell must tile one disk: with an unsplit dual
    // cell (the scaffold vertex across, `towards`), with a
    // part of a dual cell split before (`split_quads`) and with each part of
    // this cell made before it. A single piece always does: its three side
    // quadrilaterals share its vertex, and its three dual quadrilaterals
    // face three different vertices.
    let simple = |group: &[usize], done: &[Vec<usize>]| {
        let mut by_edge: HashMap<(bool, usize), Vec<&Vec<usize>>> = HashMap::new();
        for &c in group {
            let corner = &cell.corners[c];
            for (f, &across) in corner.faces.iter().zip(&corner.towards) {
                let key = match split_quads.get(&quad_key(f)) {
                    Some(&owner) => (true, owner),
                    None => (false, across),
                };
                by_edge.entry(key).or_default().push(f);
            }
        }
        by_edge
            .values()
            .all(|rings| super::voronoi_core::one_loop(rings))
            && done.iter().all(|other| {
                let shared: Vec<&Vec<usize>> = group
                    .iter()
                    .flat_map(|&c| cell.corners[c].sides.iter())
                    .filter(|(k, _)| {
                        other
                            .iter()
                            .any(|&o| cell.corners[o].sides.iter().any(|(q, _)| q == k))
                    })
                    .map(|(_, quad)| quad)
                    .collect();
                super::voronoi_core::one_loop(&shared)
            })
    };
    // `limit`: a part that passes on its own but fails on the exported mesh
    // (against a neighbour) is split again into parts of at most this volume.
    if limit.is_none() && valid(&part, earlier) && simple(&part, earlier) {
        return vec![part];
    }
    let volume = |group: &[usize]| group.iter().map(|&c| cell.corners[c].volume).sum::<f64>();
    let mut left = part;
    left.sort_by(|&a, &b| {
        cell.corners[b]
            .volume
            .total_cmp(&cell.corners[a].volume)
            .then(a.cmp(&b))
    });
    while !left.is_empty() {
        let mut group = vec![left.remove(0)];
        loop {
            let next = left.iter().position(|&c| {
                group.iter().any(|&g| touches(g, c)) && {
                    let mut trial = group.clone();
                    trial.push(c);
                    limit.is_none_or(|l| volume(&trial) <= l)
                        && valid(&trial, &parts)
                        && simple(&trial, &parts)
                }
            });
            match next {
                Some(i) => group.push(left.remove(i)),
                None => break,
            }
        }
        group.sort_unstable();
        parts.push(group);
    }
    parts.split_off(earlier.len())
}

/// The cell made of the corner pieces `group` of a dual cell, and its
/// boundary faces (the walls only when no layers cover them).
fn build_part(
    p: &[Vec3],
    cell: &Dual,
    group: &[usize],
    site: usize,
    layered: bool,
    owner: usize,
) -> (Cell, Vec<BoundaryFace>) {
    let (faces, planes, used) = part_faces(cell, group);
    let part_volume: f64 = group.iter().map(|&c| cell.corners[c].volume).sum();
    let mut boundary = Vec::new();
    for i in used {
        let (nodes, patch, _) = cell.boundary[i];
        if !is_wall(patch) || !layered {
            boundary.push(BoundaryFace {
                nodes: nodes.to_vec(),
                patch,
                owner,
            });
        }
    }
    let witness = if planes.is_empty() {
        p[site]
    } else {
        strict_witness(p, p[site], &faces, &planes, part_volume)
    };
    let built = Cell {
        element: Element::Polyhedron(Box::new(Polyhedron::new(faces, witness, None, part_volume))),
        layer: 0,
        region: 0,
    };
    (built, boundary)
}

/// A quadrilateral's nodes, sorted: the same key from both sides.
fn quad_key(f: &[usize]) -> [usize; 4] {
    let mut k = [f[0], f[1], f[2], f[3]];
    k.sort_unstable();
    k
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
/// Rounds of squeezing the layers where the wall offset spoils the scaffold.
const OFFSET_ROUNDS: usize = 40;
/// A scaffold tetrahedron must keep this share of its volume.
const OFFSET_KEEP: f64 = 1e-3;
const OFFSET_FACTOR: f64 = 0.7;
const OFFSET_MIN: f64 = 0.02;

pub(super) fn prepare(scaffold: &Scaffold, steps: &[f64]) -> Result<Prepared, String> {
    if scaffold
        .boundary
        .iter()
        .any(|&(t, tag)| tag == 0 || t.iter().any(|&n| n >= scaffold.points.len()))
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
    // The shell names each patch by its id; ids no facet uses are left out.
    let mut ids: Vec<usize> = scaffold.boundary.iter().map(|b| b.1).collect();
    ids.sort_unstable();
    ids.dedup();
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
            patch: PatchId(ids.binary_search(&tag).unwrap_or(0)),
        });
        patches.push(tag);
    }
    let surface = ClosedSurface::new(
        surface_points,
        triangles,
        ids.iter().map(|id| id.to_string()).collect(),
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
    let mut rims = BTreeMap::<usize, Vec<Edge>>::new();
    let mut port_edges = BTreeMap::new();
    for (&t, &tag) in ts.iter().zip(&patches) {
        if !is_wall(tag) {
            for (u, v) in edges(t) {
                match port_edges.remove(&(v, u)) {
                    None => {
                        port_edges.insert((u, v), tag);
                    }
                    Some(other) if other != tag && !(is_slip(other) && is_slip(tag)) => {
                        return Err(
                            "two patches without layers (inlet, outlet or no_layers) share an edge"
                                .into(),
                        );
                    }
                    Some(_) => {}
                }
            }
        }
    }
    for ((a, b), tag) in port_edges {
        rims.entry(tag).or_default().push((b, a));
    }
    // Squeeze the layers further where the wall offset would turn a
    // scaffold tetrahedron over or flatten it (thin ribs, tight corners).
    let labels = vec![0; ts.len()];
    let mut front = front::advance(&surface, &patches, &labels, &rims, steps, None)?;
    if !steps.is_empty() {
        let mut neighbours = vec![Vec::new(); scaffold.points.len()];
        for t in &scaffold.tetrahedra {
            for i in 0..4 {
                for j in 0..4 {
                    if i != j {
                        neighbours[t[i]].push(t[j]);
                    }
                }
            }
        }
        for round in 0..OFFSET_ROUNDS {
            let displacement = offset_interior(scaffold, &surface_ids, &front)?;
            let mut scale = front.scale.clone();
            let mut bad = 0usize;
            for t in &scaffold.tetrahedra {
                let [a, b, c, d] = t.map(|n| scaffold.points[n]);
                let [e, f, g, h] = t.map(|n| scaffold.points[n] + displacement[n]);
                let before = (b - a).cross(c - a).dot(d - a);
                let after = (f - e).cross(g - e).dot(h - e);
                if after.is_finite() && after > OFFSET_KEEP * before {
                    continue;
                }
                bad += 1;
                // the wall vertices of the tetrahedron, or the nearest ones
                let mut ring: Vec<usize> = t.to_vec();
                for _ in 0..4 {
                    if ring.iter().any(|n| surface_ids.contains_key(n)) {
                        break;
                    }
                    ring = ring
                        .iter()
                        .flat_map(|&n| neighbours[n].iter().copied())
                        .collect();
                    ring.sort_unstable();
                    ring.dedup();
                }
                for n in ring {
                    if let Some(&v) = surface_ids.get(&n) {
                        scale[v] = (front.scale[v] * OFFSET_FACTOR).max(OFFSET_MIN);
                    }
                }
            }
            if crate::trace_enabled() {
                eprintln!(
                    "kofl trace: wall offset round {round}: {bad} scaffold tetrahedra turned over or flattened"
                );
            }
            if bad == 0 || scale == front.scale {
                break;
            }
            front = front::advance(&surface, &patches, &labels, &rims, steps, Some(&scale))?;
        }
    }
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
            dual[t[i]].corners[corners[i]].towards.push(t[j]);
            dual[t[j]].corners[corners[j]]
                .faces
                .push(face.into_iter().rev().collect());
            dual[t[j]].corners[corners[j]].towards.push(t[i]);
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
    // Dual quadrilaterals of cells already split, by the part that owns them.
    let mut split_quads: HashMap<[usize; 4], usize> = HashMap::new();
    let layer_cells = front.cells.len();
    // Per core cell: its dual cell and corner pieces; per dual cell: its cells.
    let mut origin: Vec<(usize, Vec<usize>)> = Vec::new();
    let mut cells_of: Vec<Vec<usize>> = vec![Vec::new(); dual.len()];
    let layered = !steps.is_empty();
    for (id, cell) in dual.iter().enumerate() {
        let volume: f64 = cell.corners.iter().map(|c| c.volume).sum();
        if volume <= 0. {
            return Err(format!("unused volume node {id}"));
        }
        let site = mapped[id];
        let owner_offset = front.cells.len();
        let wall = |nodes: &[usize; 3]| prisms.get(&sorted(*nodes)).map(|&c| prism_centres[c]);
        let mut parts: Vec<Vec<usize>> = Vec::new();
        for part in split_reentrant(p, site, cell, &wall) {
            let split = split_invalid(p, cell, part, &parts, &split_quads, None);
            parts.extend(split);
        }
        if parts.len() > 1 {
            for (k, group) in parts.iter().enumerate() {
                for &c in group {
                    for f in &cell.corners[c].faces {
                        split_quads.insert(quad_key(f), owner_offset + k);
                    }
                }
            }
        }
        for (part, group) in parts.into_iter().enumerate() {
            let (built, boundary) = build_part(p, cell, &group, site, layered, owner_offset + part);
            front.cells.push(built);
            front.boundary.extend(boundary);
            cells_of[id].push(owner_offset + part);
            origin.push((id, group));
        }
    }
    // Repeat the split until the mesh as exported passes the solver rules:
    // a core cell that fails (or whose layer neighbour fails) is split again
    // into halves by volume, until it passes or is a single corner piece.
    loop {
        let failing = verify::exported_failures(p, &front.cells);
        let mut targets = std::collections::BTreeSet::new();
        for (id, neighbours) in &failing {
            if *id >= layer_cells {
                targets.insert(*id);
            } else {
                targets.extend(neighbours.iter().filter(|&&n| n >= layer_cells));
            }
        }
        targets.retain(|&c| origin[c - layer_cells].1.len() > 1);
        if crate::trace_enabled() {
            eprintln!(
                "kofl trace: split again: {} cells fail as exported, {} core cells to split",
                failing.len(),
                targets.len()
            );
        }
        if targets.is_empty() {
            break;
        }
        for c in targets {
            let (d, group) = origin[c - layer_cells].clone();
            let cell = &dual[d];
            let earlier: Vec<Vec<usize>> = cells_of[d]
                .iter()
                .filter(|&&o| o != c)
                .map(|&o| origin[o - layer_cells].1.clone())
                .collect();
            let half = group.iter().map(|&k| cell.corners[k].volume).sum::<f64>() / 2.;
            let pieces = split_invalid(p, cell, group, &earlier, &split_quads, Some(half));
            front.boundary.retain(|b| b.owner != c);
            for (k, piece) in pieces.into_iter().enumerate() {
                let id = if k == 0 { c } else { front.cells.len() };
                for &corner in &piece {
                    for f in &cell.corners[corner].faces {
                        split_quads.insert(quad_key(f), id);
                    }
                }
                let (built, boundary) = build_part(p, cell, &piece, mapped[d], layered, id);
                front.boundary.extend(boundary);
                if k == 0 {
                    front.cells[c] = built;
                    origin[c - layer_cells].1 = piece;
                } else {
                    front.cells.push(built);
                    cells_of[d].push(id);
                    origin.push((d, piece));
                }
            }
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
        patches.iter().filter(|&&p| is_wall(p)).count() * 6,
    )?;
    Ok((
        DuctMesh {
            points: front.points,
            cells: front.cells,
            boundary: front.boundary,
            measurements,
            patches: classic_patches(),
        },
        surface,
    ))
}
