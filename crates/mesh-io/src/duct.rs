use crate::msh::{GmshVersion, Model};
use crate::{end_array, end_vtu, start_array, vtu_points, write_array, write_ids};
use mesh_algorithms::duct::{DuctMesh, Element};
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

pub fn gmsh(mesh: &DuctMesh, version: GmshVersion, out: &mut dyn Write) -> io::Result<()> {
    let mut points = mesh.points.clone();
    for cell in &mesh.cells {
        if let Element::Polyhedron(p) = &cell.element
            && p.cone_apex.is_none()
        {
            points.push(p.interior);
        }
    }
    let mut model = Model::new(points);
    for (tag, name) in [(1, "inlet"), (2, "outlet"), (3, "walls")] {
        model.group(2, tag, name);
    }
    model.group(3, 4, "fluid");
    for face in &mesh.boundary {
        let kind = if face.nodes.len() == 3 { 2 } else { 3 };
        model.element(kind, face.patch, face.patch, &face.nodes)?;
    }
    let polyhedral = mesh
        .cells
        .iter()
        .any(|c| matches!(c.element, Element::Polyhedron(_)));
    let mut next_center = mesh.points.len();
    for (parent, c) in mesh.cells.iter().enumerate() {
        let entity = c.region + 1;
        let parent = polyhedral.then_some(parent);
        match &c.element {
            Element::Tet(n) => model.element_of(4, 4, entity, n, parent)?,
            Element::Wedge(n) => model.element_of(6, 4, entity, n, parent)?,
            Element::Polyhedron(p) => {
                let apex = p.cone_apex.unwrap_or_else(|| {
                    let id = next_center;
                    next_center += 1;
                    id
                });
                for face in &p.faces {
                    if face.contains(&apex) {
                        continue;
                    }
                    let mut nodes: Vec<usize> = face.iter().rev().copied().collect();
                    nodes.push(apex);
                    let kind = if face.len() == 3 { 4 } else { 7 };
                    model.element_of(kind, 4, entity, &nodes, parent)?;
                }
            }
        }
    }
    model.write(version, out)
}
