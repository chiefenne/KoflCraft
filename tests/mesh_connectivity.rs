use koflcraft::geometry::{BoxSide, Vec3};
use koflcraft::mesh::{ConnectivityError, FaceAdjacency, Mesh};

fn two_tetrahedra() -> (Vec<Vec3>, Vec<usize>, Vec<usize>, Vec<FaceAdjacency>) {
    let vertices = vec![
        Vec3::ZERO,
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 0.0, -1.0),
    ];
    let adjacency = [0, 0, 0, 1, 1, 1]
        .into_iter()
        .map(|owner| FaceAdjacency::Boundary {
            owner,
            patch: BoxSide::XMin,
        });
    (
        vertices,
        (0..=21).step_by(3).collect(),
        vec![
            0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3, 0, 4, 1, 0, 2, 4, 1, 4, 2,
        ],
        [FaceAdjacency::Interior {
            owner: 0,
            neighbour: 1,
        }]
        .into_iter()
        .chain(adjacency)
        .collect(),
    )
}

#[test]
fn shared_face_is_stored_once_and_indexed_by_both_cells() {
    let (vertices, offsets, ids, adjacency) = two_tetrahedra();
    let mesh = Mesh::from_faces(vertices, offsets, ids, adjacency, 2).unwrap();
    assert_eq!(mesh.face_count(), 7);
    assert_eq!(mesh.cell_count(), 2);
    assert_eq!(mesh.face_vertices(0), Some(&[0, 2, 1][..]));
    assert_eq!(
        mesh.face_adjacency(0),
        Some(FaceAdjacency::Interior {
            owner: 0,
            neighbour: 1
        })
    );
    assert_eq!(mesh.cell_faces(0), Some(&[0, 1, 2, 3][..]));
    assert_eq!(mesh.cell_faces(1), Some(&[0, 4, 5, 6][..]));
    assert_eq!(mesh.cell_faces(2), None);
    assert_eq!(mesh.face_vertices(7), None);
}

#[test]
fn invalid_offsets_and_incidence_are_rejected() {
    let (vertices, mut offsets, ids, mut adjacency) = two_tetrahedra();
    offsets[1] = 50;
    assert_eq!(
        Mesh::from_faces(vertices.clone(), offsets, ids.clone(), adjacency.clone(), 2),
        Err(ConnectivityError::FaceVertices { face: 0 })
    );

    adjacency[0] = FaceAdjacency::Interior {
        owner: 0,
        neighbour: 0,
    };
    let offsets = (0..=21).step_by(3).collect();
    assert_eq!(
        Mesh::from_faces(vertices, offsets, ids, adjacency, 2),
        Err(ConnectivityError::FaceCells { face: 0 })
    );
}
