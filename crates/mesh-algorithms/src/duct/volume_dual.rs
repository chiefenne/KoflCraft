//! Vertex-centred median dual of a volume-filling tetrahedral scaffold.
use super::*;
use mesh_core::geometry::PatchId;
use mesh_geometry::SurfaceTriangle;
use std::collections::{HashMap, HashSet};

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Scaffold {
    pub points: Vec<Vec3>,
    pub boundary: Vec<([usize; 3], usize)>,
    pub tetrahedra: Vec<[usize; 4]>,
    /// Zero-thickness shell triangles inside the domain (SHL-01), each a
    /// face of two tetrahedra, counter-clockwise about its shell's normal,
    /// with the shell's index. [`crack`](super::crack_shells) turns them
    /// into boundary faces on both sides.
    pub shells: Vec<([usize; 3], usize)>,
    /// After cracking: (copy, original) for every node doubled at a shell.
    pub twins: Vec<(usize, usize)>,
    /// After cracking: the boundary tags of the two sides of each shell
    /// (front: the side its normal points into).
    pub shell_sides: Vec<(usize, usize)>,
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

/// Where a cell of the mesh comes from: a wall layer cell and its
/// sub-wedges, or a core cell, its dual cell and its corner pieces.
#[derive(Clone)]
enum Origin {
    /// Sub-wedges and column (see `merge_layer_prisms`).
    Layer(Vec<[usize; 6]>, usize),
    Core(usize, Vec<usize>),
    /// A cube of a hex core: no dual part, never split.
    Cube,
}

fn sorted_ring(ring: &[usize]) -> Vec<usize> {
    let mut key = ring.to_vec();
    key.sort_unstable();
    key
}

/// Wall triangles around a wall vertex whose normals differ by more than
/// this angle (cos 30°), or that lie on different patches, get separate
/// polygon prisms: the prism stays on one side of a feature line.
const PRISM_SECTOR_COS: f64 = 0.866_025_403_784_438_6;

/// One layer cell made of sub-wedges that share side faces: its faces are
/// the sub-wedge faces not shared inside it, and the sub-wedges stay its
/// exact prism decomposition.
fn prism_cell(p: &[Vec3], pieces: Vec<[usize; 6]>, layer: usize) -> Cell {
    if let [w] = pieces[..] {
        return Cell {
            element: Element::Wedge(w),
            layer,
            region: 0,
        };
    }
    // the side quadrilaterals two sub-wedges share are inside the cell (end
    // triangles never are); faces in the order of `Element::faces`
    let sides = |[a, b, c, d, e, f]: [usize; 6]| [[a, b, e, d], [b, c, f, e], [c, a, d, f]];
    let mut keys: Vec<[usize; 4]> = pieces
        .iter()
        .flat_map(|&w| sides(w))
        .map(|mut q| {
            q.sort_unstable();
            q
        })
        .collect();
    keys.sort_unstable();
    let inside = |q: &[usize; 4]| {
        let mut key = *q;
        key.sort_unstable();
        let at = keys.partition_point(|k| *k < key);
        keys.get(at + 1) == Some(&key)
    };
    let mut faces: Vec<Vec<usize>> = Vec::with_capacity(3 * pieces.len() + 2);
    for &w in &pieces {
        let [a, b, c, d, e, f] = w;
        faces.push(vec![a, c, b]);
        faces.push(vec![d, e, f]);
        faces.extend(sides(w).iter().filter(|q| !inside(q)).map(|q| q.to_vec()));
    }
    let volume: f64 = pieces
        .iter()
        .map(|w| verify::wedge(w.map(|i| p[i])).1)
        .sum();
    let interior =
        pieces.iter().flatten().fold(Vec3::ZERO, |s, &v| s + p[v]) / (6 * pieces.len()) as f64;
    let mut poly = Polyhedron::new(faces, interior, None, volume);
    poly.prism = Some(pieces);
    Cell {
        element: Element::Polyhedron(Box::new(poly)),
        layer,
        region: 0,
    }
}

/// The polygon prisms of one wall layer: the corners of the wall triangles
/// (corner `j` of triangle `k` is `3 k + j`) grouped by wall vertex into
/// sectors, so that the triangles of a sector lie on one patch and each
/// turns by less than 30° from its neighbour across their shared edge. A
/// triangle without a patch stays apart. Returns each corner's sector, named
/// by its first corner; the number of polygon prisms per layer is the
/// number of corners that name their own sector.
pub fn prism_sectors(
    points: &[Vec3],
    triangles: &[[usize; 3]],
    patches: &[Option<usize>],
) -> Vec<usize> {
    let normal = |t: [usize; 3]| {
        let [a, b, c] = t.map(|v| points[v]);
        let n = (b - a).cross(c - a);
        n / n.length()
    };
    // union-find over the corners
    let mut parent: Vec<usize> = (0..3 * triangles.len()).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    let mut edges: Vec<([usize; 2], usize)> = Vec::with_capacity(3 * triangles.len());
    for (k, t) in triangles.iter().enumerate() {
        for j in 0..3 {
            edges.push((sorted([t[j], t[(j + 1) % 3]]), k));
        }
    }
    edges.sort_unstable();
    for pair in edges.chunk_by(|a, b| a.0 == b.0) {
        let [(edge, k), (_, l)] = pair[..] else {
            continue;
        };
        let patch = patches.get(k).copied().flatten();
        let smooth = normal(triangles[k]).dot(normal(triangles[l])) >= PRISM_SECTOR_COS;
        if patch.is_none() || patch != patches.get(l).copied().flatten() || !smooth {
            continue;
        }
        let (tk, tl) = (triangles[k], triangles[l]);
        for v in edge {
            let (Some(i), Some(j)) = (
                tk.iter().position(|&x| x == v),
                tl.iter().position(|&x| x == v),
            ) else {
                continue;
            };
            let (a, b) = (find(&mut parent, 3 * k + i), find(&mut parent, 3 * l + j));
            parent[a.max(b)] = a.min(b);
        }
    }
    (0..parent.len()).map(|x| find(&mut parent, x)).collect()
}

/// Polygon prisms, as the prism layers of production dual meshers: the
/// sub-wedges of one layer around one wall vertex (two per wall triangle,
/// see `split_layers`) become one cell, one per smooth sector of the
/// vertex's wall triangles. Every wall layer column then matches the dual
/// cell of its wall vertex instead of six slivers per wall prism. A cell
/// keeps its sub-wedges as its exact prism decomposition (`prism`), so a
/// merged cell that fails the solver rules can be split back. Returns the
/// sub-wedges of every cell, in cell order, and its column (the sector of
/// wall triangles it stands on, the same in every layer).
fn merge_layer_prisms(front: &mut front::Front) -> Vec<(Vec<[usize; 6]>, usize)> {
    let wedges: Option<Vec<[usize; 6]>> = front
        .cells
        .iter()
        .map(|c| match c.element {
            Element::Wedge(n) => Some(n),
            _ => None,
        })
        .collect();
    let identity = |front: &front::Front| -> Vec<(Vec<[usize; 6]>, usize)> {
        front
            .cells
            .iter()
            .enumerate()
            .map(|(id, c)| match c.element {
                Element::Wedge(n) => (vec![n], usize::MAX - id),
                _ => (Vec::new(), usize::MAX - id),
            })
            .collect()
    };
    let Some(wedges) = wedges else {
        return identity(front);
    };
    // six sub-wedges per wall prism, the prisms layer after layer in the
    // same wall-triangle order (front::advance)
    let coarse = wedges.len() / 6;
    let layer_of = |k: usize| front.cells[6 * k].layer;
    let per_layer = (0..coarse).take_while(|&k| layer_of(k) == 1).count();
    if wedges.len() % 6 != 0
        || per_layer == 0
        || coarse % per_layer != 0
        || (0..coarse).any(|k| {
            layer_of(k) != k / per_layer + 1
                || (1..6).any(|i| front.cells[6 * k + i].layer != layer_of(k))
        })
    {
        return identity(front);
    }
    // corner j of wall prism k: sub-wedges 2j and 2j + 1 stand on it
    let corners = |k: usize| [0, 1, 2].map(|j| wedges[6 * k + 2 * j][0]);
    let p = &front.points;
    let mut patch = vec![None; per_layer];
    for face in &front.boundary {
        if face.nodes.len() == 3 && face.owner < 6 * per_layer {
            patch[face.owner / 6] = Some(face.patch);
        }
    }
    let triangles: Vec<[usize; 3]> = (0..per_layer).map(corners).collect();
    let sector_of = prism_sectors(p, &triangles, &patch);
    // one cell per (layer, sector), in the order the sub-wedges come
    let mut group_of: HashMap<(usize, usize), usize> = HashMap::new();
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut sectors: Vec<usize> = Vec::new();
    for k in 0..coarse {
        for j in 0..3 {
            let sector = sector_of[3 * (k % per_layer) + j];
            let g = *group_of.entry((layer_of(k), sector)).or_insert_with(|| {
                groups.push(Vec::new());
                sectors.push(sector);
                groups.len() - 1
            });
            groups[g].extend([6 * k + 2 * j, 6 * k + 2 * j + 1]);
        }
    }
    // Where the first-layer prism's centre lies behind one of its wall
    // triangles (a curved wall under a thin first layer), the merged cells
    // would fail the solver rules as exported: that column starts as the
    // prisms of each wall-triangle corner (see the split in
    // `build_volume_dual`, which still checks every cell).
    let behind = |group: &[usize]| {
        let (mut volume, mut moment) = (0.0, Vec3::ZERO);
        for &s in group {
            let w = wedges[s].map(|v| p[v]);
            let v = verify::wedge(w).1;
            volume += v;
            moment = moment + w.iter().fold(Vec3::ZERO, |m, &x| m + x) * (v / 6.0);
        }
        let centre = moment / volume;
        group.iter().any(|&s| {
            let [a, b, c] = [0, 1, 2].map(|i| p[wedges[s][i]]);
            // outward normal of the wall face [a, c, b]
            (centre - a).dot((c - a).cross(b - a)) >= 0.0
        })
    };
    let cornered: std::collections::BTreeSet<usize> = groups
        .iter()
        .zip(&sectors)
        .filter(|(group, _)| front.cells[group[0]].layer == 1 && behind(group))
        .map(|(_, &sector)| sector)
        .collect();
    let mut new_id = vec![0; wedges.len()];
    let mut cells = Vec::with_capacity(groups.len());
    let mut parts = Vec::with_capacity(groups.len());
    for (group, sector) in groups.iter().zip(&sectors) {
        let layer = front.cells[group[0]].layer;
        let size = if cornered.contains(sector) {
            2
        } else {
            group.len()
        };
        for part in group.chunks(size) {
            let pieces: Vec<[usize; 6]> = part.iter().map(|&s| wedges[s]).collect();
            for &s in part {
                new_id[s] = cells.len();
            }
            cells.push(prism_cell(p, pieces.clone(), layer));
            parts.push((pieces, *sector));
        }
    }
    for face in &mut front.boundary {
        face.owner = new_id[face.owner];
    }
    front.cells = cells;
    parts
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
    let mut parts = connected_parts(cell, groups);
    parts.sort();
    if parts.len() < 2 || !parts.iter().all(|part| manifold(cell, part)) {
        return vec![all];
    }
    let after: usize = parts
        .iter()
        .map(|part| wall_inversions(p, cell, part, wall))
        .sum();
    if after < before { parts } else { vec![all] }
}

/// The corner pieces of each group in connected parts: pieces join through
/// shared tetrahedron faces.
fn connected_parts(cell: &Dual, groups: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
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
    parts
}

/// Whether the cell made of the pieces `part` is a closed, consistently
/// oriented, manifold shell.
fn manifold(cell: &Dual, part: &[usize]) -> bool {
    let (faces, _, _) = part_faces(cell, part);
    let mut edges = HashMap::<(usize, usize), (i32, usize)>::new();
    for f in &faces {
        for (k, &a) in f.iter().enumerate() {
            let b = f[(k + 1) % f.len()];
            let e = edges.entry((a.min(b), a.max(b))).or_default();
            e.0 += if a < b { 1 } else { -1 };
            e.1 += 1;
        }
    }
    edges.values().all(|&(s, u)| s == 0 && u == 2)
}

/// The dual cell of a node on a free edge of a zero-thickness shell wraps
/// round the edge, so its centroid lies near the shell plane and its faces
/// on the shell are seen almost edge-on (non-orthogonality near 90°).
/// Split it into the pieces in front of the shell plane and those behind
/// (`normal` points to the front), so each part sits on one side as a wall
/// cell does; the faces between the parts lie on the scaffold faces
/// through the node, around the plane beyond the edge. Kept only where
/// every part is a closed manifold shell.
fn split_shell_edge(p: &[Vec3], site: usize, cell: &Dual, normal: Vec3) -> Vec<Vec<usize>> {
    let all: Vec<usize> = (0..cell.corners.len()).collect();
    let mut groups = vec![Vec::new(), Vec::new()];
    for (c, corner) in cell.corners.iter().enumerate() {
        groups[usize::from((corner.centre - p[site]).dot(normal) <= 0.)].push(c);
    }
    if groups.iter().any(|g| g.is_empty()) {
        return vec![all];
    }
    let mut parts = connected_parts(cell, groups);
    parts.sort();
    if parts.len() < 2 || !parts.iter().all(|part| manifold(cell, part)) {
        return vec![all];
    }
    parts
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
    // a site at a reflex corner of its part has no strict kernel point
    // next to it (the inward step is cut to nothing): take the mean of the
    // part's nodes when every face (fanned from its first node) sees it
    let witness = if faces.iter().flatten().any(|&v| p[v] == witness) {
        let mut nodes: Vec<usize> = faces.iter().flatten().copied().collect();
        nodes.sort_unstable();
        nodes.dedup();
        let mean = nodes.iter().fold(Vec3::ZERO, |s, &v| s + p[v]) / nodes.len() as f64;
        let sees = faces.iter().all(|f| {
            (1..f.len() - 1).all(|k| {
                let (a, b, c) = (p[f[0]], p[f[k]], p[f[k + 1]]);
                (b - a).cross(c - a).dot(mean - a) < 0.0
            })
        });
        if sees { mean } else { witness }
    } else {
        witness
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
    // the shell of a hex core stays where the cubes are
    let mut pinned = vec![false; scaffold.points.len()];
    for &(t, tag) in &scaffold.boundary {
        if is_hex(tag) {
            for n in t {
                pinned[n] = true;
            }
        }
    }
    let mut stencils = Vec::new();
    for (i, adj) in neighbors.iter_mut().enumerate() {
        if surface_ids.contains_key(&i) || pinned[i] {
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

/// What [`crack_shells`] did.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CrackStats {
    /// Shell triangles, now a boundary face on each side.
    pub shell_triangles: usize,
    /// Nodes given a copy (one per extra side).
    pub copies: usize,
    /// Shell nodes kept single: on free edges, where the fluid wraps round.
    pub single_nodes: usize,
}

/// Crack the scaffold along its zero-thickness shells (SHL-01, SHL-02):
/// every shell node gets one copy per connected set of its tetrahedra
/// (connected through faces that are not shell faces), so a node inside a
/// shell or on its line on the wall gets one copy per side, and a node on a
/// free edge, where the fluid wraps round the edge, stays single. Every
/// shell triangle becomes two boundary faces, oriented out of the fluid on
/// each side: `sides[k].0` for the side its normal points into (the
/// shell's name), `sides[k].1` for the other (`<name>_back`). The copies
/// are listed in `twins`.
pub fn crack_shells(
    mut scaffold: Scaffold,
    sides: &[(usize, usize)],
) -> Result<(Scaffold, CrackStats), String> {
    let mut stats = CrackStats::default();
    if scaffold.shells.is_empty() {
        return Ok((scaffold, stats));
    }
    let p = &scaffold.points;
    let mut shell_face: HashMap<[usize; 3], usize> = HashMap::new();
    for (i, (t, k)) in scaffold.shells.iter().enumerate() {
        if *k >= sides.len() || t.iter().any(|&v| v >= p.len()) {
            return Err("invalid shell triangle".into());
        }
        if shell_face.insert(sorted(*t), i).is_some() {
            return Err("a shell triangle is given twice".into());
        }
    }
    let on_shell: std::collections::BTreeSet<usize> =
        scaffold.shells.iter().flat_map(|(t, _)| *t).collect();
    // the two tetrahedra on each shell triangle: [front, back]
    let mut owners = vec![[usize::MAX; 2]; scaffold.shells.len()];
    let mut star: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut faces: HashMap<[usize; 3], Vec<usize>> = HashMap::new();
    for (ti, tet) in scaffold.tetrahedra.iter().enumerate() {
        if !tet.iter().any(|v| on_shell.contains(v)) {
            continue;
        }
        for &v in tet {
            if on_shell.contains(&v) {
                star.entry(v).or_default().push(ti);
            }
        }
        for o in 0..4 {
            let face = sorted(std::array::from_fn(|k| tet[(o + 1 + k) % 4]));
            if !face.iter().any(|v| on_shell.contains(v)) {
                continue;
            }
            faces.entry(face).or_default().push(ti);
            if let Some(&i) = shell_face.get(&face) {
                let [a, b, c] = scaffold.shells[i].0.map(|v| p[v]);
                let normal = (b - a).cross(c - a);
                let side = usize::from((p[tet[o]] - a).dot(normal) <= 0.0);
                if owners[i][side] != usize::MAX {
                    return Err(format!("shell triangle {i} has two tetrahedra on one side"));
                }
                owners[i][side] = ti;
            }
        }
    }
    if let Some(i) = owners.iter().position(|o| o.contains(&usize::MAX)) {
        return Err(format!(
            "shell triangle {i} is not inside the fluid (a tetrahedron on one side only)"
        ));
    }
    // one copy per extra connected set of tetrahedra around each node
    let mut copy: HashMap<(usize, usize), usize> = HashMap::new();
    let mut points = scaffold.points.clone();
    for &v in &on_shell {
        let tets = &star[&v];
        let mut parent: HashMap<usize, usize> = tets.iter().map(|&t| (t, t)).collect();
        fn root(parent: &mut HashMap<usize, usize>, mut t: usize) -> usize {
            while parent[&t] != t {
                let up = parent[&parent[&t]];
                parent.insert(t, up);
                t = up;
            }
            t
        }
        for &t in tets {
            let tet = scaffold.tetrahedra[t];
            for o in 0..4 {
                let face = sorted(std::array::from_fn(|k| tet[(o + 1 + k) % 4]));
                if !face.contains(&v) || shell_face.contains_key(&face) {
                    continue;
                }
                for &n in &faces[&face] {
                    if n != t {
                        let (x, y) = (root(&mut parent, t), root(&mut parent, n));
                        parent.insert(x.max(y), x.min(y));
                    }
                }
            }
        }
        let mut sets: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
        for &t in tets {
            let r = root(&mut parent, t);
            sets.entry(r).or_default().push(t);
        }
        if sets.len() == 1 {
            stats.single_nodes += 1;
            continue;
        }
        for group in sets.into_values().skip(1) {
            let id = points.len();
            points.push(points[v]);
            scaffold.twins.push((id, v));
            stats.copies += 1;
            for t in group {
                copy.insert((t, v), id);
            }
        }
    }
    let map = |t: usize, nodes: [usize; 3]| nodes.map(|v| copy.get(&(t, v)).copied().unwrap_or(v));
    for (wall, _) in &mut scaffold.boundary {
        if wall.iter().any(|v| on_shell.contains(v)) {
            let owner = faces
                .get(&sorted(*wall))
                .and_then(|o| (o.len() == 1).then_some(o[0]))
                .ok_or("a wall triangle at a shell has no single tetrahedron")?;
            *wall = map(owner, *wall);
        }
    }
    for (i, (t, k)) in scaffold.shells.iter().enumerate() {
        let [front, back] = owners[i];
        let f = map(front, [t[0], t[2], t[1]]);
        let b = map(back, *t);
        if sorted(f) == sorted(b) {
            return Err(format!(
                "shell triangle {i} has all its corners on free edges"
            ));
        }
        scaffold.boundary.push((f, sides[*k].0));
        scaffold.boundary.push((b, sides[*k].1));
    }
    for (ti, tet) in scaffold.tetrahedra.iter_mut().enumerate() {
        for v in tet.iter_mut() {
            if let Some(&c) = copy.get(&(ti, *v)) {
                *v = c;
            }
        }
    }
    stats.shell_triangles = scaffold.shells.len();
    scaffold.points = points;
    scaffold.shells.clear();
    scaffold.shell_sides = sides.to_vec();
    Ok((scaffold, stats))
}

/// Check the scaffold, build its boundary surface and the wall layers.
/// Rounds of squeezing the layers where the wall offset spoils the scaffold.
const OFFSET_ROUNDS: usize = 40;
/// A scaffold tetrahedron must keep this share of its volume.
const OFFSET_KEEP: f64 = 1e-3;
const OFFSET_FACTOR: f64 = 0.7;
const OFFSET_MIN: f64 = 0.02;

pub(super) fn prepare(scaffold: &Scaffold, steps: &[f64]) -> Result<Prepared, String> {
    prepare_within(scaffold, steps, None)
}

/// [`prepare`] with the layer height at each surface vertex limited to the
/// share `limit` of it from the start.
fn prepare_within(
    scaffold: &Scaffold,
    steps: &[f64],
    limit: Option<&[f64]>,
) -> Result<Prepared, String> {
    if !scaffold.shells.is_empty() {
        return Err("crack the scaffold along its shells before the dual is built".into());
    }
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
    let mut ids: Vec<usize> = scaffold
        .boundary
        .iter()
        .map(|b| b.1)
        .filter(|&t| !is_hex(t))
        .collect();
    ids.sort_unstable();
    ids.dedup();
    for &(t, tag) in &scaffold.boundary {
        if is_hex(tag) {
            continue;
        }
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
    let names = ids.iter().map(|id| id.to_string()).collect();
    // a shell that touches no wall is a closed piece of its own
    let surface = if scaffold.twins.is_empty() {
        ClosedSurface::new(surface_points, triangles, names)
    } else {
        ClosedSurface::with_pieces(surface_points, triangles, names)
    }
    .map_err(|e| e.to_string())?;
    let ts: Vec<_> = surface
        .facets()
        .iter()
        .map(|f| f.triangle.vertices)
        .collect();
    // the original surface vertex of every vertex copy at a shell
    let mut twins: Vec<usize> = Vec::new();
    if !scaffold.twins.is_empty() {
        twins = (0..surface.vertices().len()).collect();
        for &(c, o) in &scaffold.twins {
            match (surface_ids.get(&c), surface_ids.get(&o)) {
                (Some(&c), Some(&o)) => twins[c] = o,
                _ => return Err("a shell node copy is not on the boundary".into()),
            }
        }
    }
    intersection::embedded(
        surface.vertices(),
        &front::contact_triangles(&ts, &twins),
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
    let mut front = front::advance(&surface, &patches, &labels, &rims, steps, limit, &twins)?;
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
        // tetrahedra still turned over after their own wall vertices yielded
        let mut stuck: HashSet<usize> = HashSet::new();
        for round in 0..OFFSET_ROUNDS {
            let displacement = offset_interior(scaffold, &surface_ids, &front)?;
            let mut scale = front.scale.clone();
            let mut bad = 0usize;
            let mut now_bad: HashSet<usize> = HashSet::new();
            for (index, t) in scaffold.tetrahedra.iter().enumerate() {
                let [a, b, c, d] = t.map(|n| scaffold.points[n]);
                let [e, f, g, h] = t.map(|n| scaffold.points[n] + displacement[n]);
                let before = (b - a).cross(c - a).dot(d - a);
                let after = (f - e).cross(g - e).dot(h - e);
                if after.is_finite() && after > OFFSET_KEEP * before {
                    continue;
                }
                bad += 1;
                now_bad.insert(index);
                // the wall vertices of the tetrahedron, or the nearest ones
                let mut ring: Vec<usize> = t.to_vec();
                // Still turned over a round later: an interior corner moves
                // with walls the tetrahedron does not touch (the third wall
                // at a box corner), so those yield too.
                if stuck.contains(&index) {
                    for &n in t {
                        if !surface_ids.contains_key(&n) {
                            ring.extend(neighbours[n].iter().copied());
                        }
                    }
                }
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
            stuck = now_bad;
            front = front::advance(
                &surface,
                &patches,
                &labels,
                &rims,
                steps,
                Some(&scale),
                &twins,
            )?;
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
    build_dual(scaffold, steps, None)
}

/// [`build_volume_dual`] on a scaffold whose boundary includes the shell of
/// `core` (tag [`HEX`]): the dual fills the gap between the wall layers and
/// the shell, and every cube of the core becomes one cell, joined to the
/// dual cells through their triangles on the shell.
pub fn build_volume_dual_hex(
    scaffold: Scaffold,
    steps: &[f64],
    core: &crate::hexcore::HexCore,
) -> Result<(DuctMesh, ClosedSurface), String> {
    build_dual(scaffold, steps, Some(core))
}

/// Rebuilds of the dual with thinner layers where layer cells still fail
/// the solver rules after every split (see [`build_dual`]).
const LAYER_RETRIES: usize = 3;

/// One build of the dual: the mesh, or the surface vertices under layer
/// cells that still fail the solver rules after every split, with the
/// layer shares of that build.
enum Attempt {
    Built(Box<(DuctMesh, ClosedSurface)>),
    Folded { walls: Vec<usize>, scale: Vec<f64> },
}

/// The dual of the scaffold with the wall layers. A layer cell that still
/// fails the solver rules when no split is left (a sub-wedge whose own
/// faces fold: the top of the stack in a concave corner, as where a tube
/// passes through a plate) is not split or merged away; the layers under
/// it are squeezed (the share of the height at its surface vertices
/// times `OFFSET_FACTOR`) and the dual is built again, up to
/// [`LAYER_RETRIES`] times: the layers retract locally, as the prism
/// layers of production meshers do in tight corners, and keep their full
/// height everywhere else.
fn build_dual(
    scaffold: Scaffold,
    steps: &[f64],
    core: Option<&crate::hexcore::HexCore>,
) -> Result<(DuctMesh, ClosedSurface), String> {
    let mut limit: Option<Vec<f64>> = None;
    for retry in 0..=LAYER_RETRIES {
        match build_dual_once(
            &scaffold,
            steps,
            core,
            limit.as_deref(),
            retry < LAYER_RETRIES,
        )? {
            Attempt::Built(built) => return Ok(*built),
            Attempt::Folded { walls, mut scale } => {
                if crate::trace_enabled() {
                    eprintln!(
                        "kofl trace: dual retry {}: layers squeezed at {} wall vertices under folded layer cells",
                        retry + 1,
                        walls.len()
                    );
                }
                for v in walls {
                    scale[v] = (scale[v] * OFFSET_FACTOR).max(OFFSET_MIN);
                }
                limit = Some(scale);
            }
        }
    }
    Err("the layer cells still fold after the retries".into())
}

fn build_dual_once(
    scaffold: &Scaffold,
    steps: &[f64],
    core: Option<&crate::hexcore::HexCore>,
    limit: Option<&[f64]>,
    may_retry: bool,
) -> Result<Attempt, String> {
    let mut stage = std::time::Instant::now();
    let mut lap = |name: &str| {
        if crate::trace_enabled() {
            eprintln!(
                "kofl trace: dual {name}: {:.2} s",
                stage.elapsed().as_secs_f64()
            );
        }
        stage = std::time::Instant::now();
    };
    let Prepared {
        surface,
        patches,
        mut front,
        surface_ids,
    } = prepare_within(scaffold, steps, limit)?;
    lap("layers");
    let mut cuts = Subdivision::default();
    split_layers(&mut front, &mut cuts)?;
    let layer_parts = merge_layer_prisms(&mut front);
    lap("layer split and polygon prisms");
    let mapped = map_scaffold(scaffold, &surface_ids, &mut front)?;
    lap("interior offset");
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
    lap("dual pieces");
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
    // Per cell: its sub-wedges or its dual cell and corner pieces; per dual
    // cell: its cells.
    let mut columns: HashMap<usize, Vec<usize>> = HashMap::new();
    for (id, (_, column)) in layer_parts.iter().enumerate() {
        columns.entry(*column).or_default().push(id);
    }
    let mut origin: Vec<Origin> = layer_parts
        .into_iter()
        .map(|(pieces, column)| Origin::Layer(pieces, column))
        .collect();
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
        // a node on a shell's free edge: both sides of one shell
        let shell_normal = scaffold.shell_sides.iter().find_map(|&(front, back)| {
            let f = cell.boundary.iter().find(|b| b.1 == front)?;
            cell.boundary.iter().find(|b| b.1 == back)?;
            let [a, b, c] = f.0.map(|n| p[n]);
            // a front face points out of the fluid, against the normal
            Some((b - a).cross(c - a) * -1.)
        });
        let first = match shell_normal {
            Some(n) => split_shell_edge(p, site, cell, n),
            None => split_reentrant(p, site, cell, &wall),
        };
        for part in first {
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
            origin.push(Origin::Core(id, group));
        }
    }
    // Repeat the split until the mesh as exported passes the solver rules:
    // a core cell that fails (or whose layer neighbour fails) is split again
    // into halves by volume, until it passes or is a single corner piece.
    // After the first check of the whole mesh only the cells near a split
    // are checked again (verify::ExportedCheck), with the same verdicts.
    lap("dual cells");
    if let Some(core) = core {
        super::hex_join::join(core, p, &mut front.cells, &mut front.boundary)?;
    }
    // the cubes of a hex core sit between the first dual cells and those
    // added by the splits below
    origin.resize(front.cells.len(), Origin::Cube);
    let mut exported = verify::ExportedCheck::new(p, &front.cells)?;
    lap("first check");
    loop {
        let failing = exported.failing();
        let mut split_core = std::collections::BTreeSet::new();
        let mut split_layer = std::collections::BTreeSet::new();
        for (id, neighbours) in failing {
            match &origin[*id] {
                Origin::Core(..) => {
                    split_core.insert(*id);
                }
                // A polygon prism that fails splits its whole column into
                // the prisms of each wall-triangle corner (two sub-wedges,
                // whose faces keep the wall vertex); one of those that fails
                // goes back to its two sub-wedges.
                Origin::Layer(pieces, column) if pieces.len() > 2 => {
                    split_layer.extend(
                        columns[column]
                            .iter()
                            .filter(|&&c| matches!(&origin[c], Origin::Layer(p, _) if p.len() > 2)),
                    );
                }
                Origin::Layer(pieces, _) if pieces.len() == 2 => {
                    split_layer.insert(*id);
                }
                // a single sub-wedge or a cube of a hex core: its core
                // neighbours split
                Origin::Layer(..) | Origin::Cube => split_core.extend(
                    neighbours
                        .iter()
                        .filter(|&&n| matches!(origin[n], Origin::Core(..))),
                ),
            }
        }
        split_core.retain(|&c| matches!(&origin[c], Origin::Core(_, group) if group.len() > 1));
        if crate::trace_enabled() {
            eprintln!(
                "kofl trace: split again: {} cells fail as exported, {} core cells to split, {} layer prisms to split",
                failing.len(),
                split_core.len(),
                split_layer.len()
            );
        }
        if split_core.is_empty() && split_layer.is_empty() {
            // layer cells that fail with nothing left to split: squeeze the
            // layers under them and build again
            let mut walls: Vec<usize> = failing
                .keys()
                .filter(|&&id| matches!(origin[id], Origin::Layer(..)))
                .flat_map(|&id| front.cells[id].element.nodes().to_vec())
                .filter_map(|n| front.wall_of.get(n).copied())
                .collect();
            walls.sort_unstable();
            walls.dedup();
            if may_retry && !walls.is_empty() && !steps.is_empty() {
                return Ok(Attempt::Folded {
                    walls,
                    scale: front.scale.clone(),
                });
            }
            break;
        }
        let mut owned: HashMap<usize, Vec<usize>> = HashMap::new();
        for (i, b) in front.boundary.iter().enumerate() {
            if split_layer.contains(&b.owner) {
                owned.entry(b.owner).or_default().push(i);
            }
        }
        for c in split_layer {
            let Origin::Layer(pieces, column) = origin[c].clone() else {
                continue;
            };
            let groups: Vec<Vec<[usize; 6]>> = if pieces.len() > 2 {
                pieces.chunks(2).map(<[_]>::to_vec).collect()
            } else {
                pieces.iter().map(|&w| vec![w]).collect()
            };
            exported.remove(c, &front.cells[c].element);
            let level = front.cells[c].layer;
            let mut ids = Vec::with_capacity(groups.len());
            for (k, group) in groups.iter().enumerate() {
                let cell = prism_cell(p, group.clone(), level);
                let id = if k == 0 {
                    front.cells[c] = cell;
                    origin[c] = Origin::Layer(group.clone(), column);
                    c
                } else {
                    front.cells.push(cell);
                    origin.push(Origin::Layer(group.clone(), column));
                    let id = front.cells.len() - 1;
                    columns.entry(column).or_default().push(id);
                    id
                };
                exported.insert(id, &front.cells[id].element);
                ids.push(id);
            }
            // each boundary face of the prism goes to the part that has it
            for &i in owned.get(&c).into_iter().flatten() {
                let key = sorted_ring(&front.boundary[i].nodes);
                if let Some(k) = groups.iter().position(|group| {
                    group.iter().any(|&w| {
                        Element::Wedge(w)
                            .faces()
                            .iter()
                            .any(|f| sorted_ring(f) == key)
                    })
                }) {
                    front.boundary[i].owner = ids[k];
                }
            }
        }
        front.boundary.retain(|b| !split_core.contains(&b.owner));
        for c in split_core {
            let Origin::Core(d, group) = origin[c].clone() else {
                continue;
            };
            let cell = &dual[d];
            let earlier: Vec<Vec<usize>> = cells_of[d]
                .iter()
                .filter(|&&o| o != c)
                .filter_map(|&o| match &origin[o] {
                    Origin::Core(_, g) => Some(g.clone()),
                    Origin::Layer(..) | Origin::Cube => None,
                })
                .collect();
            let half = group.iter().map(|&k| cell.corners[k].volume).sum::<f64>() / 2.;
            let pieces = split_invalid(p, cell, group, &earlier, &split_quads, Some(half));
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
                    exported.remove(c, &front.cells[c].element);
                    front.cells[c] = built;
                    origin[c] = Origin::Core(d, piece);
                } else {
                    front.cells.push(built);
                    cells_of[d].push(id);
                    origin.push(Origin::Core(d, piece));
                }
                exported.insert(id, &front.cells[id].element);
            }
        }
        exported.recheck(p, &front.cells)?;
    }
    lap("split again");
    // the splits above re-add the dual triangles on the shell as boundary
    front.boundary.retain(|b| !is_hex(b.patch));
    compact_points(p, &mut front.cells, &mut front.boundary);
    lap("compact");
    let measurements = verify::check_with_interface_count(
        &surface,
        &patches,
        p,
        &front.cells,
        &front.boundary,
        steps.len(),
        patches.iter().filter(|&&p| is_wall(p)).count() * 6,
    )?;
    lap("final check");
    Ok(Attempt::Built(Box::new((
        DuctMesh {
            points: front.points,
            cells: front.cells,
            boundary: front.boundary,
            measurements,
            patches: classic_patches(),
        },
        surface,
    ))))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Keep rule (SW-18). Invariant: a polygon prism takes the corners of a
    /// wall vertex whose triangles lie on one patch and turn by less than
    /// 30°; a crease or a patch border separates them. Joining across a
    /// crease or a patch border breaks this (a layer cell would then reach
    /// round an edge, or carry two patches); only this test checks the
    /// sectors apart from a built mesh.
    #[test]
    fn sectors_part_at_creases_and_patch_borders() {
        // four triangles around vertex 0; vertex 4 is folded up by 45°
        let points = [
            Vec3::new(0., 0., 0.),
            Vec3::new(1., 0., 0.),
            Vec3::new(0., 1., 0.),
            Vec3::new(-1., 0., 0.),
            Vec3::new(0., -0.7, 0.7),
        ];
        let triangles = [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1]];
        // corner 0 of each triangle stands on vertex 0
        let at_zero = |patches: &[Option<usize>]| {
            let sectors = prism_sectors(&points, &triangles, patches);
            [0, 1, 2, 3].map(|k| sectors[3 * k])
        };
        let one = [Some(3); 4];
        assert_eq!(at_zero(&one), [0, 0, 6, 6]);
        let two = [Some(3), Some(4), Some(3), Some(3)];
        assert_eq!(at_zero(&two), [0, 3, 6, 6]);
        let none = [Some(3), None, Some(3), Some(3)];
        assert_eq!(at_zero(&none), [0, 3, 6, 6]);
        // flat: one sector all round
        let mut flat = points;
        flat[4] = Vec3::new(0., -1., 0.);
        let sectors = prism_sectors(&flat, &triangles, &one);
        assert_eq!([0, 1, 2, 3].map(|k| sectors[3 * k]), [0; 4]);
    }

    /// Keep rule (SW-18). Invariant: a polygon prism's faces are the faces of
    /// its sub-wedges except the side quadrilaterals two of them share, its
    /// volume is theirs, and the sub-wedges stay its prism decomposition.
    /// Keeping a shared quadrilateral leaves a face inside the cell; only
    /// this test checks the cell apart from the solver rules of a whole mesh.
    #[test]
    fn polygon_prism_drops_the_shared_sides() {
        // a unit square split into two triangles, extruded by 1
        let base = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]];
        let points: Vec<Vec3> = [0., 1.]
            .iter()
            .flat_map(|&z| base.map(|[x, y]| Vec3::new(x, y, z)))
            .collect();
        let pieces = vec![[0, 1, 2, 4, 5, 6], [0, 2, 3, 4, 6, 7]];
        let cell = prism_cell(&points, pieces.clone(), 1);
        let Element::Polyhedron(poly) = &cell.element else {
            panic!("two sub-wedges make a polyhedron");
        };
        assert_eq!(poly.faces.len(), 8);
        assert!(
            poly.faces
                .iter()
                .all(|f| sorted_ring(f) != vec![0, 2, 4, 6])
        );
        assert!((poly.expected_volume - 1.).abs() < 1e-12);
        assert_eq!(poly.prism.as_ref(), Some(&pieces));
        let shell = verify::cell_shell(&points, &poly.faces);
        assert!(!shell.open && (shell.volume - 1.).abs() < 1e-12);
        // one sub-wedge stays a wedge
        let single = prism_cell(&points, vec![pieces[0]], 1);
        assert!(matches!(single.element, Element::Wedge(_)));
    }

    /// Four tetrahedra round the edge ab; the shell triangle abc ends at
    /// the free edge ab.
    fn fan() -> Scaffold {
        let v = |x, y, z| Vec3::new(x, y, z);
        let (a, b, c, d, e, f) = (0, 1, 2, 3, 4, 5);
        let tetrahedra = vec![[a, b, c, d], [a, b, d, f], [a, b, f, e], [a, b, e, c]];
        let boundary = [
            [a, c, d],
            [b, c, d],
            [a, d, f],
            [b, d, f],
            [a, f, e],
            [b, f, e],
            [a, e, c],
            [b, e, c],
        ]
        .into_iter()
        .map(|t| (t, 1))
        .collect();
        Scaffold {
            points: vec![
                v(0., 0., 0.),
                v(1., 0., 0.),
                v(0.5, 1., 0.),
                v(0.5, 0., 1.),
                v(0.5, 0., -1.),
                v(0.5, -1., 0.),
            ],
            boundary,
            tetrahedra,
            shells: vec![([a, b, c], 0)],
            ..Scaffold::default()
        }
    }

    #[test]
    fn crack_doubles_the_shell_and_keeps_its_free_edge() {
        let (s, stats) = crack_shells(fan(), &[(7, 8)]).unwrap();
        assert_eq!(
            (stats.shell_triangles, stats.copies, stats.single_nodes),
            (1, 1, 2)
        );
        assert_eq!(s.twins, vec![(6, 2)]);
        assert_eq!(s.points[6], s.points[2]);
        assert!(s.shells.is_empty());
        assert_eq!(s.shell_sides, vec![(7, 8)]);
        // the tetrahedra above and below the shell now share the free edge only
        let shared = s.tetrahedra[0]
            .iter()
            .filter(|v| s.tetrahedra[3].contains(v))
            .count();
        assert_eq!(shared, 2);
        // every boundary face lies on one tetrahedron, both shell sides too
        assert_eq!(s.boundary.len(), 10);
        for (t, _) in &s.boundary {
            assert!(
                s.tetrahedra
                    .iter()
                    .any(|tet| t.iter().all(|v| tet.contains(v))),
                "{t:?}"
            );
        }
        let sides: Vec<_> = s.boundary[8..].iter().map(|b| b.1).collect();
        assert_eq!(sides, vec![7, 8]);
        // the front side faces away from the tetrahedron it bounds
        for (t, _) in &s.boundary[8..] {
            let owner = s
                .tetrahedra
                .iter()
                .find(|tet| t.iter().all(|v| tet.contains(v)))
                .unwrap();
            let apex = owner.iter().find(|v| !t.contains(v)).unwrap();
            let [a, b, c] = t.map(|n| s.points[n]);
            assert!((b - a).cross(c - a).dot(s.points[*apex] - a) < 0.0);
        }
    }

    #[test]
    fn crack_refuses_a_shell_on_the_fluid_boundary() {
        let mut s = fan();
        s.tetrahedra.pop();
        let e = crack_shells(s, &[(7, 8)]).unwrap_err();
        assert!(e.contains("not inside the fluid"), "{e}");
    }
}
