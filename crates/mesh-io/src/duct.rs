use crate::handoff::{CellFaces, HandoffError, HandoffMesh, patches_from_names};
use crate::msh::{GmshVersion, Model};
use crate::{end_array, end_vtu, start_array, vtu_points, write_array, write_ids};
use mesh_algorithms::duct::{DuctMesh, Element};
use mesh_core::geometry::Vec3;
use std::collections::HashMap;
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
                for face in p.faces.iter().filter(|f| f.len() > 4 || f.len() == 3) {
                    let tris: Vec<[usize; 3]> = ends
                        .iter()
                        .filter(|t| t.iter().all(|v| face.contains(v)))
                        .copied()
                        .collect();
                    caps.insert(sorted(face), tris);
                }
            }
        }
        Self {
            caps,
            centres: HashMap::new(),
        }
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

/// A one-wedge layer cell that is still a plain wedge (5 faces of at most
/// 4 nodes) and so can be written as a native Gmsh prism.
fn native(p: &mesh_algorithms::duct::Polyhedron) -> bool {
    p.prism.as_ref().is_some_and(|w| w.len() == 1)
        && p.faces.len() == 5
        && p.faces.iter().all(|f| f.len() <= 4)
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
        boundary_split.push((face.nodes.len() > 4).then(|| splits.split(&mut points, &face.nodes)));
    }
    let mut cell_split = Vec::with_capacity(mesh.cells.len());
    for cell in &mesh.cells {
        cell_split.push(match &cell.element {
            Element::Polyhedron(p) if !native(p) => p
                .faces
                .iter()
                .map(|f| (f.len() > 4).then(|| splits.split(&mut points, f)))
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
        let merged = (rings.len() > 1).then(|| single_loop(&rings)).flatten();
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
