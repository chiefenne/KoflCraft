use std::collections::{BTreeMap, BTreeSet};

use koflcraft::clipping::{ConvexPolyhedron, FaceSource};
use koflcraft::generate_mesh;
use koflcraft::geometry::{AxisAlignedBox, Vec3};
use koflcraft::mesh::{FaceAdjacency, Mesh};
use koflcraft::mesh_export::{ExportError, PolyMeshWriters, write_gmsh, write_openfoam, write_vtu};

fn unit_box() -> AxisAlignedBox {
    AxisAlignedBox::new(Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap()
}

// These widely spaced seeds produce internal faces and polygons with varying vertex counts.
fn sites() -> [Vec3; 4] {
    [
        Vec3::new(0.15, 0.22, 0.33),
        Vec3::new(0.72, 0.28, 0.11),
        Vec3::new(0.60, 0.72, 0.54),
        Vec3::new(0.30, 0.50, 0.82),
    ]
}

#[test]
fn two_cell_box_has_analytic_metrics() {
    let mesh = generate_mesh(
        unit_box(),
        &[Vec3::new(0.2, 0.5, 0.5), Vec3::new(0.6, 0.5, 0.5)],
    )
    .unwrap();
    let metrics = mesh.metrics();
    assert_eq!(metrics.boundary_face_count, 10);
    assert_eq!(metrics.faces.len(), 11);
    for (cell, (volume, x)) in metrics.cells.iter().zip([(0.4, 0.2), (0.6, 0.7)]) {
        assert!((cell.volume - volume).abs() < 1e-14);
        assert!((cell.centroid - Vec3::new(x, 0.5, 0.5)).length() < 1e-14);
        assert_eq!(cell.face_count, 6);
        assert_eq!(cell.neighbour_count, 1);
    }
    let interior = (0..mesh.mesh().face_count())
        .find(|&face| {
            mesh.mesh()
                .face_adjacency(face)
                .unwrap()
                .neighbour()
                .is_some()
        })
        .unwrap();
    let face = metrics.faces[interior];
    assert!((face.centroid - Vec3::new(0.4, 0.5, 0.5)).length() < 1e-14);
    assert!((face.normal - Vec3::new(1.0, 0.0, 0.0)).length() < 1e-14);
    assert!((face.area - 1.0).abs() < 1e-14);
    assert!(face.non_orthogonality_degrees < 1e-6);
    assert!(face.skewness < 1e-14);
}

#[test]
#[ignore = "outside the VTU and Gmsh acceptance scope"]
fn openfoam_files_parse_as_a_conforming_partition() {
    let domain = unit_box();
    let generated = generate_mesh(domain, &sites()).unwrap();
    let mesh = generated.mesh();
    let [
        mut points,
        mut faces,
        mut owner,
        mut neighbour,
        mut boundary,
    ] = std::array::from_fn(|_| Vec::new());
    write_openfoam(
        domain,
        mesh,
        PolyMeshWriters {
            points: &mut points,
            faces: &mut faces,
            owner: &mut owner,
            neighbour: &mut neighbour,
            boundary: &mut boundary,
        },
    )
    .unwrap();
    let points = read_list(&points, |tokens| {
        tokens.expect("(");
        let point = Vec3::new(tokens.read(), tokens.read(), tokens.read());
        tokens.expect(")");
        point
    });
    let faces = read_list(&faces, |tokens| {
        let len: usize = tokens.read();
        tokens.expect("(");
        let ids: Vec<usize> = (0..len).map(|_| tokens.read()).collect();
        tokens.expect(")");
        ids
    });
    let owner: Vec<usize> = read_list(&owner, |tokens| tokens.read());
    let neighbour: Vec<usize> = read_list(&neighbour, |tokens| tokens.read());
    let patches = read_list(&boundary, |tokens| {
        let name = tokens.next().to_owned();
        tokens.expect("{");
        tokens.expect("type");
        tokens.expect("patch");
        tokens.expect(";");
        tokens.expect("nFaces");
        let count: usize = tokens.read();
        tokens.expect(";");
        tokens.expect("startFace");
        let start: usize = tokens.read();
        tokens.expect(";");
        tokens.expect("}");
        (name, start, count)
    });

    assert_eq!(points, mesh.vertices());
    assert_eq!(faces.len(), mesh.face_count());
    assert_eq!(owner.len(), faces.len());
    assert_eq!(
        faces.len() - neighbour.len(),
        generated.metrics().boundary_face_count
    );
    assert!(!neighbour.is_empty());
    let mut cell_faces = vec![Vec::new(); mesh.cell_count()];
    for (face, vertices) in faces.iter().enumerate() {
        cell_faces[owner[face]].push(vertices.clone());
        if let Some(&other) = neighbour.get(face) {
            assert!(owner[face] < other);
            cell_faces[other].push(vertices.iter().rev().copied().collect());
        }
    }
    assert_eq!(patches.len(), 6);
    let mut previous_patch_end = neighbour.len();
    for (axis, expected) in ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]
        .into_iter()
        .enumerate()
    {
        let (name, start, count) = &patches[axis];
        assert_eq!(name, expected);
        assert_eq!(*start, previous_patch_end);
        assert!(*count > 0);
        previous_patch_end = start + count;
        for vertices in &faces[*start..previous_patch_end] {
            for &id in vertices {
                let point = points[id];
                let coordinate = [point.x, point.y, point.z][axis / 2];
                assert!((coordinate - (axis % 2) as f64).abs() < 1e-13);
            }
        }
    }
    assert_eq!(previous_patch_end, faces.len());
    let total_volume: f64 = cell_faces
        .iter()
        .zip(sites())
        .map(|(faces, site)| shell_volume(faces, &points, site))
        .sum();
    assert!((total_volume - 1.0).abs() < 1e-13);
}

#[test]
fn vtu_contains_complete_outward_polyhedron_faces() {
    let domain = unit_box();
    let generated = generate_mesh(domain, &sites()).unwrap();
    let mut output = Vec::new();
    write_vtu(domain, generated.mesh(), &mut output).unwrap();
    let output = String::from_utf8(output).unwrap();
    let arrays = |name: &str| {
        let attribute = format!("Name=\"{name}\"");
        let (_, array) = output.split_once(&attribute).unwrap();
        let (_, data) = array.split_once('>').unwrap();
        let (values, _) = data.split_once("</DataArray>").unwrap();
        values
            .split_whitespace()
            .map(|id| id.parse::<usize>().unwrap())
            .collect::<Vec<_>>()
    };
    let connectivity = arrays("connectivity");
    let offsets = arrays("offsets");
    let types = arrays("types");
    let faces = arrays("faces");
    let face_offsets = arrays("faceoffsets");
    assert_eq!(types, vec![42; sites().len()]);
    assert_eq!(offsets.len(), types.len());
    assert_eq!(face_offsets.len(), types.len());
    let mut face_start = 0;
    let mut connectivity_start = 0;
    let mut volume = 0.0;
    for (cell, site) in sites().into_iter().enumerate() {
        let mut stream = faces[face_start..face_offsets[cell]].iter().copied();
        let face_count = stream.next().unwrap();
        let mut cell_faces = Vec::new();
        for _ in 0..face_count {
            let count = stream.next().unwrap();
            let face: Vec<_> = stream.by_ref().take(count).collect();
            assert_eq!(face.len(), count);
            cell_faces.push(face);
        }
        assert_eq!(stream.next(), None);
        let points: BTreeSet<_> = cell_faces.iter().flatten().copied().collect();
        assert!(
            points
                .into_iter()
                .eq(connectivity[connectivity_start..offsets[cell]]
                    .iter()
                    .copied())
        );
        volume += shell_volume(&cell_faces, generated.mesh().vertices(), site);
        face_start = face_offsets[cell];
        connectivity_start = offsets[cell];
    }
    assert_eq!(face_start, faces.len());
    assert_eq!(connectivity_start, connectivity.len());
    assert!((volume - 1.0).abs() < 1e-13);
}

#[test]
fn failed_validation_writes_no_vtu_data() {
    let generated = generate_mesh(unit_box(), &sites()).unwrap();
    let other_domain = AxisAlignedBox::new(Vec3::ZERO, Vec3::new(2.0, 2.0, 2.0)).unwrap();
    let mut output = Vec::new();
    assert!(matches!(
        write_vtu(other_domain, generated.mesh(), &mut output),
        Err(ExportError::Validation(_))
    ));
    assert!(output.is_empty());
}

#[test]
fn gmsh_triangulation_keeps_vertices_on_split_boundary_edges() {
    let cube = ConvexPolyhedron::from_box(unit_box());
    let mut vertices = cube.vertices().to_vec();
    let midpoint = vertices.len();
    vertices.push(Vec3::new(0.0, 0.0, 0.5));
    let mut offsets = vec![0];
    let mut connectivity = Vec::new();
    let mut adjacency = Vec::new();
    for face in cube.faces() {
        let FaceSource::Boundary(patch) = face.source() else {
            unreachable!()
        };
        let ring = face.vertices();
        for (&a, &b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
            connectivity.push(a);
            if matches!((a, b), (0, 4) | (4, 0)) {
                connectivity.push(midpoint);
            }
        }
        offsets.push(connectivity.len());
        adjacency.push(FaceAdjacency::Boundary { owner: 0, patch });
    }
    let mesh = Mesh::from_faces(vertices, offsets, connectivity, adjacency, 1).unwrap();
    let mut msh = Vec::new();
    write_gmsh(unit_box(), &mesh, &mut msh).unwrap();
    let text = String::from_utf8(msh).unwrap();
    assert!(text.contains("$Nodes\n12\n"));
    assert!(text.contains("9 0e0 0e0 5e-1\n"));
}

struct Tokens<'a>(std::str::SplitWhitespace<'a>);

impl<'a> Tokens<'a> {
    fn next(&mut self) -> &'a str {
        self.0.next().unwrap()
    }

    fn expect(&mut self, s: &str) {
        assert_eq!(self.next(), s);
    }

    fn read<T: std::str::FromStr>(&mut self) -> T
    where
        T::Err: std::fmt::Debug,
    {
        self.next().parse().unwrap()
    }
}

fn read_list<T>(data: &[u8], mut element: impl FnMut(&mut Tokens<'_>) -> T) -> Vec<T> {
    let text = std::str::from_utf8(data).unwrap();
    let (header, body) = text.split_once('}').unwrap();
    assert!(header.contains("FoamFile"));
    assert!(header.contains("format ascii;"));
    let separated = ['(', ')', '{', '}', ';']
        .into_iter()
        .fold(body.to_owned(), |body, c| {
            body.replace(c, &format!(" {c} "))
        });
    let mut tokens = Tokens(separated.split_whitespace());
    let count: usize = tokens.read();
    tokens.expect("(");
    let values = (0..count).map(|_| element(&mut tokens)).collect();
    tokens.expect(")");
    assert_eq!(tokens.0.next(), None);
    values
}

fn shell_volume(faces: &[Vec<usize>], points: &[Vec3], interior: Vec3) -> f64 {
    let mut edges = BTreeMap::<(usize, usize), [u8; 2]>::new();
    let mut volume = 0.0;
    for face in faces {
        assert!(face.len() >= 3);
        let a = points[face[0]] - interior;
        let face_volume: f64 = face[1..]
            .windows(2)
            .map(|pair| a.dot((points[pair[0]] - interior).cross(points[pair[1]] - interior)) / 6.0)
            .sum();
        assert!(face_volume > 0.0);
        volume += face_volume;
        for (&a, &b) in face.iter().zip(face.iter().cycle().skip(1)) {
            edges.entry((a.min(b), a.max(b))).or_default()[usize::from(a > b)] += 1;
        }
    }
    assert!(edges.values().all(|&uses| uses == [1, 1]));
    volume
}
