use crate::handoff::{CellFaces, HandoffError, HandoffMesh, patches_from_names};
use crate::msh::{GmshVersion, Model};
use crate::{end_array, end_vtu, start_array, vtu_points, write_array, write_ids};
use mesh_algorithms::duct::{DuctMesh, Element};
use mesh_core::geometry::Vec3;
use std::collections::{HashMap, HashSet};
use std::io::{self, Write};

pub fn volume(mesh: &DuctMesh, out: &mut dyn Write) -> io::Result<()> {
    vtu_points(&mesh.points, mesh.cells.len(), out)?;
    start_array(out, "Int64", "connectivity")?;
    let mut offsets = Vec::with_capacity(mesh.cells.len());
    let mut n = 0;
    for cell in &mesh.cells {
        let nodes = cell.element.nodes();
        write_ids(out, nodes, false)?;
        writeln!(out)?;
        n += nodes.len();
        offsets.push(n);
    }
    end_array(out)?;
    write_array(out, "offsets", &offsets)?;
    start_array(out, "UInt8", "types")?;
    for cell in &mesh.cells {
        writeln!(
            out,
            "{}",
            match cell.element {
                Element::Tet(_) => 10,
                Element::Wedge(_) => 13,
                Element::Polyhedron(_) => 42,
            }
        )?;
    }
    end_array(out)?;
    if mesh
        .cells
        .iter()
        .any(|c| matches!(c.element, Element::Polyhedron(_)))
    {
        start_array(out, "Int64", "faces")?;
        let mut offset = 0_i64;
        let mut face_offsets = Vec::with_capacity(mesh.cells.len());
        for cell in &mesh.cells {
            let Element::Polyhedron(ref p) = cell.element else {
                face_offsets.push(-1);
                continue;
            };
            write!(out, "{}", p.faces.len())?;
            offset += 1;
            for f in &p.faces {
                write!(out, " {} ", f.len())?;
                write_ids(out, f, false)?;
                offset += 1 + f.len() as i64;
            }
            writeln!(out)?;
            face_offsets.push(offset);
        }
        end_array(out)?;
        start_array(out, "Int64", "faceoffsets")?;
        for end in face_offsets {
            writeln!(out, "{end}")?;
        }
        end_array(out)?;
    }
    writeln!(out, "</Cells>\n<CellData Scalars=\"LayerIndex\">")?;
    write_array(
        out,
        "LayerIndex",
        &mesh.cells.iter().map(|c| c.layer).collect::<Vec<_>>(),
    )?;
    write_array(
        out,
        "RegionId",
        &mesh.cells.iter().map(|c| c.region).collect::<Vec<_>>(),
    )?;
    writeln!(out, "</CellData>")?;
    end_vtu(out)
}

pub fn patch(mesh: &DuctMesh, tag: usize, out: &mut dyn Write) -> io::Result<()> {
    let faces = mesh
        .boundary
        .iter()
        .filter(|f| f.patch == tag)
        .map(|f| (f.nodes.as_slice(), f.owner, mesh.cells[f.owner].layer));
    crate::channel::write_boundary_vtu(&mesh.points, tag, faces, out)
}

/// Triangles of a polygon face for Gmsh, in the ring's orientation: the
/// end triangles of the poly-prism wedges where the face is a prism cap (so
/// both sides of the face split it alike), otherwise a fan from a centre
/// point shared by both sides.
struct Splits {
    caps: HashMap<Vec<usize>, Vec<[usize; 3]>>,
    /// Faces fanned from a centre although they have at most four corners:
    /// one of each pair of a cell's faces that share more than one edge.
    /// Cones over both faces would be two elements sharing two faces.
    folds: HashSet<Vec<usize>>,
    centres: HashMap<Vec<usize>, usize>,
}

impl Splits {
    fn new(mesh: &DuctMesh) -> Self {
        let mut caps = HashMap::new();
        for cell in &mesh.cells {
            if let Element::Polyhedron(p) = &cell.element
                && native(p)
                && let Some(wedges) = &p.prism
            {
                let ends: Vec<[usize; 3]> = wedges
                    .iter()
                    .flat_map(|w| [[w[0], w[1], w[2]], [w[3], w[4], w[5]]])
                    .collect();
                for face in &p.faces {
                    let tris = end_triangles(&ends, face);
                    if !tris.is_empty() {
                        caps.insert(sorted(face), tris);
                    }
                }
            }
        }
        // faces of tetrahedra, wedges and wall columns keep their shape
        let mut fixed = HashSet::new();
        for cell in &mesh.cells {
            match &cell.element {
                Element::Tet(n) => {
                    for f in [[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]] {
                        fixed.insert(sorted(&f.map(|i| n[i])));
                    }
                }
                Element::Wedge(n) => {
                    for f in [[0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]] {
                        fixed.insert(sorted(&f.map(|i| n[i])));
                    }
                }
                Element::Polyhedron(p) if native(p) => {
                    for f in &p.faces {
                        fixed.insert(sorted(f));
                    }
                }
                Element::Polyhedron(_) => {}
            }
        }
        let mut folds = HashSet::new();
        for cell in &mesh.cells {
            let Element::Polyhedron(p) = &cell.element else {
                continue;
            };
            if native(p) {
                continue;
            }
            for (i, j) in folded_pairs(&p.faces) {
                let (fi, fj) = (&p.faces[i], &p.faces[j]);
                let fanned = |f: &Vec<usize>| {
                    f.len() > 4 || caps.contains_key(&sorted(f)) || folds.contains(&sorted(f))
                };
                if fanned(fi) || fanned(fj) {
                    continue;
                }
                if let Some(f) = [fi, fj].into_iter().find(|f| !fixed.contains(&sorted(f))) {
                    folds.insert(sorted(f));
                }
            }
        }
        Self {
            caps,
            folds,
            centres: HashMap::new(),
        }
    }

    /// Whether a face is split: tiled by the end triangles of a native layer cell
    /// (a quadrilateral cap of a wall column is two triangles, not a quad), or
    /// a fold (see `folds`).
    fn is_cap(&self, ring: &[usize]) -> bool {
        let key = sorted(ring);
        self.caps.contains_key(&key) || self.folds.contains(&key)
    }

    fn split(&mut self, points: &mut Vec<Vec3>, ring: &[usize]) -> Vec<[usize; 3]> {
        let key = sorted(ring);
        if let Some(tris) = self.caps.get(&key) {
            let n = newell(points, ring);
            return tris
                .iter()
                .map(|&[a, b, c]| {
                    if (points[b] - points[a]).cross(points[c] - points[a]).dot(n) < 0.0 {
                        [a, c, b]
                    } else {
                        [a, b, c]
                    }
                })
                .collect();
        }
        let centre = *self.centres.entry(key).or_insert_with(|| {
            let c = ring.iter().fold(Vec3::ZERO, |s, &v| s + points[v]) / ring.len() as f64;
            points.push(c);
            points.len() - 1
        });
        (0..ring.len())
            .map(|k| [centre, ring[k], ring[(k + 1) % ring.len()]])
            .collect()
    }
}

/// The pairs of a cell's faces that share more than one edge (in face
/// order): a vertex between them belongs to no other face of the cell.
fn folded_pairs(faces: &[Vec<usize>]) -> Vec<(usize, usize)> {
    let mut first: HashMap<(usize, usize), usize> = HashMap::new();
    let mut shared: std::collections::BTreeMap<(usize, usize), usize> =
        std::collections::BTreeMap::new();
    for (i, f) in faces.iter().enumerate() {
        for k in 0..f.len() {
            let (a, b) = (f[k], f[(k + 1) % f.len()]);
            let edge = (a.min(b), a.max(b));
            if let Some(&j) = first.get(&edge) {
                *shared.entry((j, i)).or_insert(0) += 1;
            } else {
                first.insert(edge, i);
            }
        }
    }
    shared
        .into_iter()
        .filter(|&(_, n)| n > 1)
        .map(|(pair, _)| pair)
        .collect()
}

/// The end triangles of a cell's wedges that lie in one of its faces.
fn end_triangles(ends: &[[usize; 3]], face: &[usize]) -> Vec<[usize; 3]> {
    ends.iter()
        .filter(|t| t.iter().all(|v| face.contains(v)))
        .copied()
        .collect()
}

/// A layer cell that can be written as its native Gmsh prisms: every face
/// is a wedge side (a quadrilateral) or exactly tiled by wedge end
/// triangles. A wall column qualifies; a split cap, whose faces carry core
/// vertices that are no wedge corners, does not.
fn native(p: &mesh_algorithms::duct::Polyhedron) -> bool {
    let Some(wedges) = &p.prism else {
        return false;
    };
    let ends: Vec<[usize; 3]> = wedges
        .iter()
        .flat_map(|w| [[w[0], w[1], w[2]], [w[3], w[4], w[5]]])
        .collect();
    p.faces.iter().all(|f| {
        let tris = end_triangles(&ends, f).len();
        tris + 2 == f.len() || (f.len() == 4 && tris == 0)
    })
}

fn sorted(ring: &[usize]) -> Vec<usize> {
    let mut k = ring.to_vec();
    k.sort_unstable();
    k
}

fn newell(points: &[Vec3], ring: &[usize]) -> Vec3 {
    let mut n = Vec3::ZERO;
    for k in 0..ring.len() {
        let (a, b) = (points[ring[k]], points[ring[(k + 1) % ring.len()]]);
        n = n + Vec3::new(
            (a.y - b.y) * (a.z + b.z),
            (a.z - b.z) * (a.x + b.x),
            (a.x - b.x) * (a.y + b.y),
        );
    }
    n
}

pub fn gmsh(mesh: &DuctMesh, version: GmshVersion, out: &mut dyn Write) -> io::Result<()> {
    let mut points = mesh.points.clone();
    let mut apexes = Vec::with_capacity(mesh.cells.len());
    for cell in &mesh.cells {
        apexes.push(match &cell.element {
            Element::Polyhedron(p) if !native(p) => Some(p.cone_apex.unwrap_or_else(|| {
                points.push(p.interior);
                points.len() - 1
            })),
            _ => None,
        });
    }
    // polygons beyond quadrilaterals are split into triangles first
    let mut splits = Splits::new(mesh);
    let mut boundary_split = Vec::with_capacity(mesh.boundary.len());
    for face in &mesh.boundary {
        let split = face.nodes.len() > 4 || splits.is_cap(&face.nodes);
        boundary_split.push(split.then(|| splits.split(&mut points, &face.nodes)));
    }
    let mut cell_split = Vec::with_capacity(mesh.cells.len());
    for cell in &mesh.cells {
        cell_split.push(match &cell.element {
            Element::Polyhedron(p) if !native(p) => p
                .faces
                .iter()
                .map(|f| (f.len() > 4 || splits.is_cap(f)).then(|| splits.split(&mut points, f)))
                .collect(),
            _ => Vec::new(),
        });
    }
    let mut model = Model::new(points);
    for (tag, name) in [(1, "inlet"), (2, "outlet"), (3, "walls")] {
        model.group(2, tag, name);
    }
    model.group(3, 4, "fluid");
    for (face, split) in mesh.boundary.iter().zip(&boundary_split) {
        match split {
            Some(tris) => {
                for t in tris {
                    model.element(2, face.patch, face.patch, t)?;
                }
            }
            None => {
                let kind = if face.nodes.len() == 3 { 2 } else { 3 };
                model.element(kind, face.patch, face.patch, &face.nodes)?;
            }
        }
    }
    let polyhedral = mesh
        .cells
        .iter()
        .any(|c| matches!(c.element, Element::Polyhedron(_)));
    for (parent, c) in mesh.cells.iter().enumerate() {
        let entity = c.region + 1;
        let parent_id = polyhedral.then_some(parent);
        match &c.element {
            Element::Tet(n) => model.element_of(4, 4, entity, n, parent_id)?,
            Element::Wedge(n) => model.element_of(6, 4, entity, n, parent_id)?,
            Element::Polyhedron(p) if native(p) => {
                for w in p.prism.iter().flatten() {
                    model.element_of(6, 4, entity, w, parent_id)?;
                }
            }
            Element::Polyhedron(p) => {
                let apex = apexes[parent].expect("apex of a general polyhedron");
                for (face, split) in p.faces.iter().zip(&cell_split[parent]) {
                    if face.contains(&apex) {
                        continue;
                    }
                    match split {
                        Some(tris) => {
                            for t in tris {
                                model.element_of(
                                    4,
                                    4,
                                    entity,
                                    &[t[2], t[1], t[0], apex],
                                    parent_id,
                                )?;
                            }
                        }
                        None => {
                            let mut nodes: Vec<usize> = face.iter().rev().copied().collect();
                            nodes.push(apex);
                            let kind = if face.len() == 3 { 4 } else { 7 };
                            model.element_of(kind, 4, entity, &nodes, parent_id)?;
                        }
                    }
                }
            }
        }
    }
    model.write(version, out)
}

/// The duct as the FoehnCraft hand-off mesh: tetrahedra, wedges and
/// polyhedra as face rings, patches inlet, outlet, walls; one zone `fluid`
/// (the regions of the partition are a construction aid, not cell zones).
pub fn handoff(mesh: &DuctMesh) -> Result<HandoffMesh, HandoffError> {
    let cells: Vec<CellFaces> = mesh
        .cells
        .iter()
        .map(|cell| {
            let faces: Vec<Vec<usize>> = match &cell.element {
                Element::Tet(n) => [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]]
                    .iter()
                    .map(|f| f.iter().map(|&i| n[i]).collect())
                    .collect(),
                Element::Wedge(n) => {
                    let tri = [[0, 2, 1], [3, 4, 5]];
                    let quad = [[0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]];
                    tri.iter()
                        .map(|f| f.iter().map(|&i| n[i]).collect())
                        .chain(quad.iter().map(|f| f.iter().map(|&i| n[i]).collect()))
                        .collect()
                }
                Element::Polyhedron(p) => p.faces.clone(),
            };
            CellFaces { faces, zone: 0 }.outward(&mesh.points)
        })
        .collect();
    let mut cells: Vec<CellFaces> = cells;
    let boundary = merge_boundary(&mut cells, mesh);
    HandoffMesh::from_cells(
        mesh.points.clone(),
        cells,
        vec!["fluid".to_string()],
        &boundary,
        patches_from_names(&["inlet", "outlet", "walls"]),
    )
}

/// Boundary faces of one cell on one patch, merged into one polygon where
/// they form a single loop, as the hand-off already does for the faces two
/// cells share. A wide wall-layer cell then has one wall face whose centroid
/// lies over the cell's, instead of many small ones beside it.
fn merge_boundary(cells: &mut [CellFaces], mesh: &DuctMesh) -> Vec<(Vec<usize>, usize)> {
    let mut groups: std::collections::BTreeMap<(usize, usize), Vec<&Vec<usize>>> =
        std::collections::BTreeMap::new();
    for f in &mesh.boundary {
        groups.entry((f.owner, f.patch)).or_default().push(&f.nodes);
    }
    let mut out = Vec::with_capacity(mesh.boundary.len());
    for ((owner, patch), rings) in groups {
        let tag = patch.saturating_sub(1);
        let merged = (rings.len() > 1)
            .then(|| single_loop(&rings))
            .flatten()
            .filter(|ring| planar(&mesh.points, ring, &rings));
        let Some(ring) = merged else {
            out.extend(rings.into_iter().map(|r| (r.clone(), tag)));
            continue;
        };
        let keys: Vec<Vec<usize>> = rings.iter().map(|r| sorted(r)).collect();
        let cell = &mut cells[owner];
        // the cell lists the same faces, possibly reversed by `outward`
        let reversed = cell
            .faces
            .iter()
            .find(|f| sorted(f) == keys[0])
            .is_some_and(|f| !same_cycle(f, rings[0]));
        cell.faces.retain(|f| !keys.contains(&sorted(f)));
        let mut face = ring.clone();
        if reversed {
            face.reverse();
        }
        cell.faces.push(face.clone());
        out.push((face, tag));
    }
    out
}

/// Whether a merged boundary polygon is planar to roundoff. FoehnCraft fans
/// a polygon about its node average, so a bent wall face would move the wall
/// and change the volume; such faces stay as their planar pieces. Every
/// vertex of the pieces is tested, not only the merged loop: pieces with a
/// flat rim around a raised inner vertex are no flat polygon.
fn planar(points: &[Vec3], ring: &[usize], pieces: &[&Vec<usize>]) -> bool {
    let centre = ring.iter().fold(Vec3::ZERO, |s, &v| s + points[v]) / ring.len() as f64;
    let n = newell(points, ring);
    let len = n.length();
    if len.is_nan() || len <= 0.0 {
        return false;
    }
    let unit = n / len;
    let size = ring
        .iter()
        .map(|&v| (points[v] - centre).length())
        .fold(0.0, f64::max);
    let tol = 1e-10 * size + 256.0 * f64::EPSILON * centre.max_abs_component().max(size);
    pieces
        .iter()
        .flat_map(|r| r.iter())
        .chain(ring)
        .all(|&v| (points[v] - centre).dot(unit).abs() <= tol)
}

fn same_cycle(a: &[usize], b: &[usize]) -> bool {
    let n = a.len();
    n == b.len()
        && a.iter()
            .position(|&v| v == b[0])
            .is_some_and(|s| (0..n).all(|i| a[(s + i) % n] == b[i]))
}

fn single_loop(rings: &[&Vec<usize>]) -> Option<Vec<usize>> {
    let mut edges: HashMap<(usize, usize), i32> = HashMap::new();
    for ring in rings {
        for k in 0..ring.len() {
            let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
            if let Some(n) = edges.get_mut(&(b, a)) {
                *n -= 1;
                if *n == 0 {
                    edges.remove(&(b, a));
                }
            } else {
                *edges.entry((a, b)).or_insert(0) += 1;
            }
        }
    }
    let mut next = HashMap::new();
    for (&(a, b), &n) in &edges {
        if n != 1 || next.insert(a, b).is_some() {
            return None;
        }
    }
    let &start = next.keys().min()?;
    let mut ring = vec![start];
    let mut at = next[&start];
    while at != start {
        if ring.len() >= next.len() {
            return None;
        }
        ring.push(at);
        at = *next.get(&at)?;
    }
    (ring.len() == next.len()).then_some(ring)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Keep: the hand-off merges the boundary faces one cell has on a patch
    /// only when the result is flat. A bent wall face would be fanned about
    /// its node average by FoehnCraft and move the wall (0.14% of the
    /// Manifold volume with merged wall columns). Removing the test in
    /// `merge_boundary` breaks this; the Manifold volume check of the
    /// workflow test catches it too, but only after a long build.
    #[test]
    fn only_flat_boundary_polygons_merge() {
        let p = [
            [0., 0., 0.],
            [1., 0., 0.],
            [1., 1., 0.],
            [0., 1., 0.],
            [0., 1., 0.01],
            [0.5, 0.5, 0.01],
        ]
        .map(|[x, y, z]| Vec3::new(x, y, z));
        assert!(planar(&p, &[0, 1, 2, 3], &[]));
        assert!(!planar(&p, &[0, 1, 2, 4], &[]));
        // a flat rim around a raised inner vertex (a tent) is not flat
        let tent = [0, 1, 2, 3].map(|k| vec![k, (k + 1) % 4, 5]);
        let refs: Vec<&Vec<usize>> = tent.iter().collect();
        assert!(!planar(&p, &[0, 1, 2, 3], &refs));
    }

    /// Keep: two quadrilateral faces of one polyhedron that share the edges
    /// a-v and v-b (v belongs to no other face) are found as a fold, so the
    /// Gmsh export fans one of them about its centre. Pyramids over both
    /// would share two triangles, which FoehnCraft refuses ("cells ... share
    /// more than one face"); this happened on Manifold after the wall flips.
    #[test]
    fn faces_sharing_two_edges_are_a_fold() {
        // a = 0, v = 1, b = 2; the faces close a wedge-like cell
        let faces = vec![
            vec![0, 1, 2, 3],
            vec![2, 1, 0, 4],
            vec![0, 3, 4],
            vec![2, 4, 3],
        ];
        assert_eq!(folded_pairs(&faces), vec![(0, 1)]);
        let cube = vec![
            vec![0, 1, 2, 3],
            vec![4, 7, 6, 5],
            vec![0, 4, 5, 1],
            vec![1, 5, 6, 2],
            vec![2, 6, 7, 3],
            vec![3, 7, 4, 0],
        ];
        assert!(folded_pairs(&cube).is_empty());
    }
}
