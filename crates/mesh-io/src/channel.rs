//! Straight-channel export. Both formats keep the wall cells as six-node prisms.

use std::collections::BTreeMap;
use std::io::{self, Write};

use mesh_algorithms::channel::{ChannelCell, ChannelMesh, ChannelPatch};

use crate::msh::{GmshVersion, Model};
use crate::{end_vtu, serialize_vtu, vtu_points, write_array};

pub fn write_volume(channel: &ChannelMesh, out: &mut dyn Write) -> io::Result<()> {
    serialize_vtu(channel.mesh(), Some(channel.cells()), out)
}

pub fn write_patch(
    channel: &ChannelMesh,
    patch: ChannelPatch,
    out: &mut dyn Write,
) -> io::Result<()> {
    let mesh = channel.mesh();
    let faces = channel
        .boundary()
        .filter(|&(_, _, name)| name == patch)
        .map(|(face, owner, _)| {
            (
                mesh.face_vertices(face).unwrap(),
                owner,
                channel.cells()[owner].layer(),
            )
        });
    write_boundary_vtu(mesh.vertices(), patch as usize, faces, out)
}

pub(super) fn write_boundary_vtu<'a>(
    vertices: &[mesh_core::geometry::Vec3],
    patch: usize,
    faces: impl Iterator<Item = (&'a [usize], usize, usize)>,
    out: &mut dyn Write,
) -> io::Result<()> {
    let mut points = Vec::new();
    let mut local_ids = BTreeMap::new();
    let mut connectivity = Vec::new();
    let mut offsets = Vec::new();
    let mut types = Vec::new();
    let mut owners = Vec::new();
    let mut layers = Vec::new();
    for (ring, owner, layer) in faces {
        for &global in ring {
            let local = *local_ids.entry(global).or_insert_with(|| {
                points.push(vertices[global]);
                points.len() - 1
            });
            connectivity.push(local);
        }
        offsets.push(connectivity.len());
        types.push(match ring.len() {
            3 => 5,
            4 => 9,
            _ => 7,
        }); // VTK triangle, quad or polygon
        owners.push(owner);
        layers.push(layer);
    }
    vtu_points(&points, owners.len(), out)?;
    write_array(out, "connectivity", &connectivity)?;
    write_array(out, "offsets", &offsets)?;
    crate::start_array(out, "UInt8", "types")?;
    crate::write_ids(out, &types, false)?;
    crate::end_array(out)?;
    writeln!(
        out,
        "      </Cells>\n      <CellData Scalars=\"LayerIndex\">"
    )?;
    write_array(out, "PatchId", &vec![patch; owners.len()])?;
    write_array(out, "OwnerCellId", &owners)?;
    write_array(out, "LayerIndex", &layers)?;
    writeln!(out, "      </CellData>")?;
    end_vtu(out)
}

pub const MULTIBLOCK: &[u8] = br#"<?xml version="1.0"?>
<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">
  <vtkMultiBlockDataSet>
    <DataSet index="0" name="volume" file="volume.vtu"/>
    <DataSet index="1" name="inlet" file="inlet.vtu"/>
    <DataSet index="2" name="outlet" file="outlet.vtu"/>
    <DataSet index="3" name="walls" file="walls.vtu"/>
  </vtkMultiBlockDataSet>
</VTKFile>
"#;

/// A core cell is coned to its centre, retaining its original triangular and
/// quadrilateral faces. All pieces use the same elementary volume tag as their parent.
pub fn write_gmsh(
    channel: &ChannelMesh,
    version: GmshVersion,
    out: &mut dyn Write,
) -> io::Result<()> {
    let mesh = channel.mesh();
    let mut points = mesh.vertices().to_vec();
    for cell in channel.cells() {
        if let ChannelCell::Core { nodes } = cell {
            let a = points[nodes[0]];
            points.push(a + (points[nodes[6]] - a) * 0.5);
        }
    }
    let names = ChannelPatch::ALL.map(ChannelPatch::name);
    let fluid = names.len() + 1;
    let mut model = Model::new(points);
    for (patch, name) in names.iter().enumerate() {
        model.group(2, patch + 1, name);
    }
    model.group(3, fluid, "fluid");
    for (face, _, patch) in channel.boundary() {
        let ring = mesh.face_vertices(face).unwrap();
        let kind = if ring.len() == 3 { 2 } else { 3 };
        model.element(kind, patch as usize, patch as usize, ring)?;
    }
    let mut centre = mesh.vertices().len();
    for (cell, shape) in channel.cells().iter().enumerate() {
        match shape {
            ChannelCell::Prism { nodes, .. } => model.element(6, fluid, cell + 1, nodes)?,
            ChannelCell::Core { .. } => {
                for &face in mesh.cell_faces(cell).unwrap() {
                    let mut base = mesh.face_vertices(face).unwrap().to_vec();
                    if mesh.face_adjacency(face).unwrap().owner() == cell {
                        base.reverse();
                    }
                    let kind = if base.len() == 3 { 4 } else { 7 };
                    base.push(centre);
                    model.element(kind, fluid, cell + 1, &base)?;
                }
                centre += 1;
            }
        }
    }
    model.write(version, out)
}
