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

pub fn gmsh(mesh: &DuctMesh, out: &mut dyn Write) -> io::Result<()> {
    let mut points = mesh.points.clone();
    let mut piece_counts = Vec::with_capacity(mesh.cells.len());
    for cell in &mesh.cells {
        let count = match &cell.element {
            Element::Polyhedron(p) => match p.cone_apex {
                Some(apex) => p.faces.iter().filter(|f| !f.contains(&apex)).count(),
                None => {
                    points.push(p.interior);
                    p.faces.len()
                }
            },
            _ => 1,
        };
        piece_counts.push(count);
    }
    let volume_elements = piece_counts.iter().sum::<usize>();
    crate::gmsh::write_header(
        out,
        &points,
        &["inlet", "outlet", "walls"],
        mesh.boundary.len() + volume_elements,
    )?;
    let mut id = 0;
    let mut elem = |kind, patch, entity, nodes: &[usize]| -> io::Result<()> {
        id += 1;
        write!(out, "{id} {kind} 2 {patch} {entity}")?;
        for &n in nodes {
            write!(out, " {}", n + 1)?;
        }
        writeln!(out)
    };
    for face in &mesh.boundary {
        elem(
            if face.nodes.len() == 3 { 2 } else { 3 },
            face.patch,
            face.patch,
            &face.nodes,
        )?;
    }
    let mut next_center = mesh.points.len();
    for c in &mesh.cells {
        let entity = c.region + 1;
        match &c.element {
            Element::Tet(n) => elem(4, 4, entity, n)?,
            Element::Wedge(n) => elem(6, 4, entity, n)?,
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
                    let mut nodes = [0; 5];
                    for (dest, src) in nodes.iter_mut().zip(face.iter().rev()) {
                        *dest = *src;
                    }
                    nodes[face.len()] = apex;
                    elem(
                        if face.len() == 3 { 4 } else { 7 },
                        4,
                        entity,
                        &nodes[..=face.len()],
                    )?;
                }
            }
        }
    }
    writeln!(out, "$EndElements")?;
    if mesh
        .cells
        .iter()
        .any(|c| matches!(c.element, Element::Polyhedron(_)))
    {
        writeln!(
            out,
            "$ElementData\n1\n\"ParentCellId\"\n1\n0\n3\n0\n1\n{volume_elements}"
        )?;
        let mut id = mesh.boundary.len() + 1;
        for (parent, count) in piece_counts.into_iter().enumerate() {
            for _ in 0..count {
                writeln!(out, "{id} {parent}")?;
                id += 1;
            }
        }
        writeln!(out, "$EndElementData")?;
    }
    Ok(())
}
