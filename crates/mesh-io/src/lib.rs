//! Serializers for VTK XML polyhedra and Gmsh meshes.

pub mod block;
pub mod channel;
mod directory;
mod gmsh;
pub mod handoff;
mod msh;

pub use directory::{BLOCK, Bundle, MeshDirectory, StreamError};
pub use gmsh::{write_gmsh, write_gmsh_domain};
pub use msh::GmshVersion;

use std::collections::BTreeSet;
use std::fmt;
use std::io::{self, Write};

use mesh_algorithms::channel::ChannelCell;
use mesh_core::geometry::{AxisAlignedBox, BoxSide, Vec3};
use mesh_core::mesh::{FaceAdjacency, Mesh};
use mesh_quality::{MeshDomain, ValidationError, validate_domain, validate_mesh};

/// The five output streams for an OpenFOAM `constant/polyMesh` directory.
pub struct PolyMeshWriters<'a> {
    pub points: &'a mut dyn Write,
    pub faces: &'a mut dyn Write,
    pub owner: &'a mut dyn Write,
    pub neighbour: &'a mut dyn Write,
    pub boundary: &'a mut dyn Write,
}

#[derive(Debug)]
pub enum ExportError {
    Validation(ValidationError),
    InvalidSubdivision { cell: usize, face: Option<usize> },
    Io(io::Error),
}

impl fmt::Display for ExportError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Validation(error) => error.fmt(f),
            Self::InvalidSubdivision {
                cell,
                face: Some(face),
            } => {
                write!(
                    f,
                    "cell {cell}, face {face}: Gmsh subdivision has a nonpositive tetrahedron"
                )
            }
            Self::InvalidSubdivision { cell, face: None } => {
                write!(
                    f,
                    "cell {cell}: Gmsh subdivision does not conserve cell volume"
                )
            }
            Self::Io(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for ExportError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Validation(error) => Some(error),
            Self::Io(error) => Some(error),
            Self::InvalidSubdivision { .. } => None,
        }
    }
}

/// Validate the mesh and serialize it to the five caller-provided streams.
pub fn write_openfoam(
    domain: AxisAlignedBox,
    mesh: &Mesh,
    out: PolyMeshWriters<'_>,
) -> Result<(), ExportError> {
    validate_mesh(domain, mesh).map_err(ExportError::Validation)?;
    serialize_openfoam(mesh, out).map_err(ExportError::Io)
}

fn serialize_openfoam(mesh: &Mesh, out: PolyMeshWriters<'_>) -> io::Result<()> {
    let mut internal = Vec::new();
    let mut patches: [Vec<_>; BoxSide::ALL.len()] = Default::default();
    for face in 0..mesh.face_count() {
        match mesh.face_adjacency(face).unwrap() {
            FaceAdjacency::Interior { owner, neighbour } => {
                internal.push((owner.min(neighbour), owner.max(neighbour), face));
            }
            FaceAdjacency::Boundary { patch, .. } => {
                patches[patch as usize].push(face);
            }
            FaceAdjacency::Surface { .. } => unreachable!("box validation checks boundary labels"),
        }
    }
    internal.sort_unstable();

    foam_list(out.points, "vectorField", "points", mesh.vertices().len())?;
    for point in mesh.vertices() {
        writeln!(out.points, "({:e} {:e} {:e})", point.x, point.y, point.z)?;
    }
    writeln!(out.points, ")")?;

    foam_list(out.faces, "faceList", "faces", mesh.face_count())?;
    foam_list(out.owner, "labelList", "owner", mesh.face_count())?;
    foam_list(out.neighbour, "labelList", "neighbour", internal.len())?;
    for &(owner, neighbour, face) in &internal {
        let vertices = mesh.face_vertices(face).unwrap();
        let reverse = mesh.face_adjacency(face).unwrap().owner() != owner;
        foam_face(out.faces, vertices, reverse)?;
        writeln!(out.owner, "{owner}")?;
        writeln!(out.neighbour, "{neighbour}")?;
    }
    writeln!(out.neighbour, ")")?;

    foam_list(out.boundary, "polyBoundaryMesh", "boundary", patches.len())?;
    let mut start = internal.len();
    for (patch, faces) in BoxSide::ALL.into_iter().zip(&patches) {
        let name = patch.name();
        writeln!(
            out.boundary,
            "{name}\n{{\n    type patch;\n    nFaces {};\n    startFace {start};\n}}",
            faces.len()
        )?;
        start += faces.len();
        for &face in faces {
            foam_face(out.faces, mesh.face_vertices(face).unwrap(), false)?;
            writeln!(out.owner, "{}", mesh.face_adjacency(face).unwrap().owner())?;
        }
    }
    writeln!(out.faces, ")")?;
    writeln!(out.owner, ")")?;
    writeln!(out.boundary, ")")
}

fn foam_list(out: &mut dyn Write, class: &str, object: &str, size: usize) -> io::Result<()> {
    writeln!(
        out,
        "FoamFile\n{{\n    version 2.0;\n    format ascii;\n    class {class};\n    location \"constant/polyMesh\";\n    object {object};\n}}\n\n{size}\n("
    )
}

fn foam_face(out: &mut dyn Write, vertices: &[usize], reverse: bool) -> io::Result<()> {
    write!(out, "{}(", vertices.len())?;
    write_ids(out, vertices, reverse)?;
    writeln!(out, ")")
}

/// Validate the mesh and write a VTK XML grid with outward polygon faces for every cell.
pub fn write_vtu(
    domain: AxisAlignedBox,
    mesh: &Mesh,
    out: &mut dyn Write,
) -> Result<(), ExportError> {
    write_vtu_domain(MeshDomain::Box(domain), mesh, out)
}

pub fn write_vtu_domain(
    domain: MeshDomain<'_>,
    mesh: &Mesh,
    out: &mut dyn Write,
) -> Result<(), ExportError> {
    validate_domain(domain, mesh).map_err(ExportError::Validation)?;
    serialize_vtu(mesh, None, out).map_err(ExportError::Io)
}

fn vtu_points(points: &[Vec3], cells: usize, out: &mut dyn Write) -> io::Result<()> {
    writeln!(
        out,
        r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="{}" NumberOfCells="{}">
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">"#,
        points.len(),
        cells
    )?;
    for point in points {
        writeln!(out, "          {:e} {:e} {:e}", point.x, point.y, point.z)?;
    }
    writeln!(out, "        </DataArray>\n      </Points>\n      <Cells>")?;

    Ok(())
}

fn serialize_vtu(
    mesh: &Mesh,
    shapes: Option<&[ChannelCell]>,
    out: &mut dyn Write,
) -> io::Result<()> {
    vtu_points(mesh.vertices(), mesh.cell_count(), out)?;
    start_array(out, "Int64", "connectivity")?;
    let mut offsets = Vec::with_capacity(mesh.cell_count());
    let mut point_ids = BTreeSet::new();
    let mut offset = 0;
    for cell in 0..mesh.cell_count() {
        if let Some(shapes) = shapes {
            let nodes: &[usize] = match &shapes[cell] {
                ChannelCell::Prism { nodes, .. } => nodes,
                ChannelCell::Core { nodes } => nodes,
            };
            write_ids(out, nodes, false)?;
            writeln!(out)?;
            offset += nodes.len();
            offsets.push(offset);
            continue;
        }
        for &face in mesh.cell_faces(cell).unwrap() {
            point_ids.extend(mesh.face_vertices(face).unwrap().iter().copied());
        }
        write!(out, "         ")?;
        for point in &point_ids {
            write!(out, " {point}")?;
        }
        writeln!(out)?;
        offset += point_ids.len();
        offsets.push(offset);
        point_ids.clear();
    }
    end_array(out)?;
    write_array(out, "offsets", &offsets)?;

    start_array(out, "UInt8", "types")?;
    let is_wedge =
        |cell: usize| shapes.is_some_and(|s| matches!(s[cell], ChannelCell::Prism { .. }));
    for cell in 0..mesh.cell_count() {
        writeln!(out, "          {}", if is_wedge(cell) { 13 } else { 42 })?;
    }
    end_array(out)?;

    start_array(out, "Int64", "faces")?;
    let mut face_offsets = Vec::with_capacity(mesh.cell_count());
    offset = 0;
    for cell in 0..mesh.cell_count() {
        if is_wedge(cell) {
            face_offsets.push(-1);
            continue;
        }
        let faces = mesh.cell_faces(cell).unwrap();
        write!(out, "          {}", faces.len())?;
        offset += 1;
        for &face in faces {
            let vertices = mesh.face_vertices(face).unwrap();
            let reverse = mesh.face_adjacency(face).unwrap().owner() != cell;
            write!(out, " {} ", vertices.len())?;
            write_ids(out, vertices, reverse)?;
            offset += 1 + vertices.len();
        }
        writeln!(out)?;
        face_offsets.push(offset as i64);
    }
    end_array(out)?;
    start_array(out, "Int64", "faceoffsets")?;
    for offset in face_offsets {
        writeln!(out, "          {offset}")?;
    }
    end_array(out)?;
    writeln!(out, "      </Cells>")?;
    if let Some(shapes) = shapes {
        writeln!(out, "      <CellData Scalars=\"LayerIndex\">")?;
        write_array(
            out,
            "LayerIndex",
            &shapes.iter().map(ChannelCell::layer).collect::<Vec<_>>(),
        )?;
        writeln!(out, "      </CellData>")?;
    }
    end_vtu(out)
}

fn end_vtu(out: &mut dyn Write) -> io::Result<()> {
    writeln!(out, "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>")
}

fn start_array(out: &mut dyn Write, value_type: &str, name: &str) -> io::Result<()> {
    writeln!(
        out,
        r#"        <DataArray type="{value_type}" Name="{name}" format="ascii">"#
    )
}

fn end_array(out: &mut dyn Write) -> io::Result<()> {
    writeln!(out, "        </DataArray>")
}

fn write_array(out: &mut dyn Write, name: &str, values: &[usize]) -> io::Result<()> {
    start_array(out, "Int64", name)?;
    write!(out, "          ")?;
    write_ids(out, values, false)?;
    writeln!(out)?;
    end_array(out)
}

fn write_ids(out: &mut dyn Write, values: &[usize], reverse: bool) -> io::Result<()> {
    let mut ids = values.iter();
    for index in 0..values.len() {
        if index != 0 {
            write!(out, " ")?;
        }
        let id = if reverse { ids.next_back() } else { ids.next() };
        write!(out, "{}", id.unwrap())?;
    }
    Ok(())
}

pub mod duct;

pub mod sweep;
