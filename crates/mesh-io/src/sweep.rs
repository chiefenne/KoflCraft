//! Mixed polygonal prisms on a section swept along a straight axis.
use crate::{channel::write_boundary_vtu, end_array, end_vtu, start_array, vtu_points, write_ids};
use mesh_algorithms::sweep::Sweep;
use mesh_core::geometry::Vec3;
use std::io::{self, Write};

pub const BLIND_MULTIBLOCK: &[u8] = br#"<?xml version="1.0"?>
<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">
  <vtkMultiBlockDataSet>
    <DataSet index="0" name="volume" file="volume.vtu"/>
    <DataSet index="1" name="inlet" file="inlet.vtu"/>
    <DataSet index="2" name="walls" file="walls.vtu"/>
  </vtkMultiBlockDataSet>
</VTKFile>
"#;

fn vtk_type(s: &Sweep, k: usize, t: usize) -> usize {
    match (s.layer(k, t), s.tiles[t].ring.len()) {
        (1.., 4) => 12,
        (1.., 3) => 13,
        _ => 42,
    }
}
fn cap(s: &Sweep, k: usize, t: usize) -> Vec<usize> {
    s.tiles[t].ring.iter().map(|&v| s.point_id(k, v)).collect()
}
fn side(s: &Sweep, k: usize, [a, b]: [usize; 2]) -> [usize; 4] {
    [
        s.point_id(k, a),
        s.point_id(k, b),
        s.point_id(k + 1, b),
        s.point_id(k + 1, a),
    ]
}
fn visit(s: &Sweep, mut write: impl FnMut(usize, usize) -> io::Result<()>) -> io::Result<()> {
    for k in 0..s.axial_layers.len() {
        for t in 0..s.tiles.len() {
            write(k, t)?;
        }
    }
    Ok(())
}

pub fn volume(s: &Sweep, points: &[Vec3], out: &mut dyn Write) -> io::Result<()> {
    vtu_points(points, s.cells(), out)?;
    start_array(out, "Int64", "connectivity")?;
    visit(s, |k, t| {
        for z in [k, k + 1] {
            write_ids(out, &cap(s, z, t), false)?;
            write!(out, " ")?;
        }
        writeln!(out)
    })?;
    end_array(out)?;
    start_array(out, "Int64", "offsets")?;
    let mut offset = 0;
    visit(s, |_, t| {
        offset += 2 * s.tiles[t].ring.len();
        writeln!(out, "{offset}")
    })?;
    end_array(out)?;
    start_array(out, "UInt8", "types")?;
    visit(s, |k, t| writeln!(out, "{}", vtk_type(s, k, t)))?;
    end_array(out)?;
    start_array(out, "Int64", "faces")?;
    let mut ends = Vec::with_capacity(s.cells());
    let mut n = 0_i64;
    visit(s, |k, t| {
        if vtk_type(s, k, t) != 42 {
            ends.push(-1);
            return Ok(());
        }
        let ring = &s.tiles[t].ring;
        let count = ring.len();
        write!(out, "{}", count + 2)?;
        for z in [k, k + 1] {
            write!(out, " {count} ")?;
            write_ids(out, &cap(s, z, t), z == k)?;
        }
        for (&a, &b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
            write!(out, " 4 ")?;
            write_ids(out, &side(s, k, [a, b]), false)?;
        }
        n += (3 + 7 * count) as i64;
        ends.push(n);
        writeln!(out)
    })?;
    end_array(out)?;
    start_array(out, "Int64", "faceoffsets")?;
    for v in ends {
        writeln!(out, "{v}")?;
    }
    end_array(out)?;
    writeln!(out, "</Cells>\n<CellData Scalars=\"LayerIndex\">")?;
    start_array(out, "Int64", "LayerIndex")?;
    visit(s, |k, t| writeln!(out, "{}", s.layer(k, t)))?;
    end_array(out)?;
    for (name, metric) in [
        (
            "CellSizeMm",
            (|s: &Sweep, k, t| s.volume(k, t).cbrt()) as fn(&Sweep, usize, usize) -> f64,
        ),
        ("MaxNeighbourSizeRatio", Sweep::neighbour_ratio),
        ("MaxInternalNonOrthogonalityDeg", |s, _, t| {
            s.tiles[t].non_orthogonality
        }),
        ("MaxInternalSkewness", |s, _, t| s.tiles[t].skewness),
    ] {
        start_array(out, "Float64", name)?;
        visit(s, |k, t| writeln!(out, "{:e}", metric(s, k, t)))?;
        end_array(out)?;
    }
    writeln!(out, "</CellData>")?;
    end_vtu(out)
}

struct Boundary {
    ring: Vec<usize>,
    owner: usize,
    tag: usize,
    layer: usize,
}
fn boundary(s: &Sweep) -> Vec<Boundary> {
    let last = s.axial_layers.len();
    let mut faces = Vec::new();
    for t in 0..s.tiles.len() {
        let mut inlet = cap(s, 0, t);
        inlet.reverse();
        faces.push(Boundary {
            ring: inlet,
            owner: t,
            tag: 1,
            layer: s.layer(0, t),
        });
        faces.push(Boundary {
            ring: cap(s, last, t),
            owner: (last - 1) * s.tiles.len() + t,
            tag: if s.profile.blind() { 3 } else { 2 },
            layer: s.layer(last - 1, t),
        });
    }
    for k in 0..last {
        for &(edge, t) in &s.outside {
            faces.push(Boundary {
                ring: side(s, k, edge).to_vec(),
                tag: 3,
                owner: k * s.tiles.len() + t,
                layer: s.layer(k, t),
            });
        }
    }
    faces
}

pub fn patches(s: &Sweep, points: &[Vec3]) -> io::Result<[Vec<u8>; 3]> {
    let faces = boundary(s);
    let mut streams: [Vec<u8>; 3] = Default::default();
    for (i, stream) in streams.iter_mut().enumerate() {
        if i == 1 && s.profile.blind() {
            continue;
        }
        let tag = i + 1;
        write_boundary_vtu(
            points,
            tag,
            faces
                .iter()
                .filter(|f| f.tag == tag)
                .map(|f| (f.ring.as_slice(), f.owner, f.layer)),
            stream,
        )?;
    }
    Ok(streams)
}

// Retain a triangle or quad natively, or fan a convex polygon into triangles.
// Both axial incidences use exactly the same fan and diagonals.
fn pieces(ring: &[usize]) -> Vec<Vec<usize>> {
    if ring.len() <= 4 {
        return vec![ring.to_vec()];
    }
    ring[1..]
        .windows(2)
        .map(|p| vec![ring[0], p[0], p[1]])
        .collect()
}

pub fn gmsh(s: &Sweep, points: &[Vec3], out: &mut dyn Write) -> io::Result<()> {
    let pieces: Vec<_> = s.tiles.iter().map(|t| pieces(&t.ring)).collect();
    let per_section: usize = pieces.iter().map(Vec::len).sum();
    let last = s.axial_layers.len();
    let volume_count = per_section * last;
    let surface_count = 2 * per_section + s.outside.len() * last;
    let tags: &[(usize, &str)] = if s.profile.blind() {
        &[(1, "inlet"), (3, "walls")]
    } else {
        &[(1, "inlet"), (2, "outlet"), (3, "walls")]
    };
    crate::gmsh::write_tagged_header(out, points, tags, 4, surface_count + volume_count)?;
    let mut id = 0;
    let mut element = |kind: usize, tag: usize, nodes: &[usize]| -> io::Result<()> {
        id += 1;
        write!(out, "{id} {kind} 2 {tag} {tag}")?;
        for &n in nodes {
            write!(out, " {}", n + 1)?;
        }
        writeln!(out)
    };
    for (k, tag) in [(0, 1), (last, if s.profile.blind() { 3 } else { 2 })] {
        for base in pieces.iter().flatten() {
            let mut nodes: Vec<_> = base.iter().map(|&v| s.point_id(k, v)).collect();
            if k == 0 {
                nodes.reverse();
            }
            element(if base.len() == 4 { 3 } else { 2 }, tag, &nodes)?;
        }
    }
    for k in 0..last {
        for &(edge, _) in &s.outside {
            element(3, 3, &side(s, k, edge))?;
        }
    }
    for k in 0..last {
        for base in pieces.iter().flatten() {
            let nodes: Vec<_> = [k, k + 1]
                .into_iter()
                .flat_map(|z| base.iter().map(move |&v| s.point_id(z, v)))
                .collect();
            element(if base.len() == 4 { 5 } else { 6 }, 4, &nodes)?;
        }
    }
    writeln!(
        out,
        "$EndElements\n$ElementData\n1\n\"ParentCellId\"\n1\n0\n3\n0\n1\n{volume_count}"
    )?;
    let mut id = surface_count;
    for k in 0..last {
        for (t, p) in pieces.iter().enumerate() {
            for _ in p {
                id += 1;
                writeln!(out, "{id} {}", k * s.tiles.len() + t)?;
            }
        }
    }
    writeln!(out, "$EndElementData")
}
