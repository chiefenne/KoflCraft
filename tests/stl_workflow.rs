use koflcraft::certification::{MeshDomain, validate_domain};
use koflcraft::geometry::{AxisAlignedBox, Vec3};
use koflcraft::mesh::FaceAdjacency;
use koflcraft::mesh_export::{write_gmsh_domain, write_vtu_domain};
use koflcraft::surface_geometry::{ConvexSurface, StlOptions, SurfaceError, import_stl};
use koflcraft::voronoi::VoronoiError;
use koflcraft::{GenerationError, generate_mesh, generate_surface_mesh};

const UNIT_BOX: &[u8] = include_bytes!("../test-fixtures/unit-box.stl");
const SKEWED_BOX: &[u8] = include_bytes!("../test-fixtures/skewed-box.stl");
const PATCHES: [&str; 6] = [
    "west inlet",
    "east outlet",
    "south wall",
    "north wall",
    "floor",
    "ceiling",
];
const SITES: [Vec3; 4] = [
    Vec3::new(0.15, 0.22, 0.33),
    Vec3::new(0.72, 0.28, 0.11),
    Vec3::new(0.60, 0.72, 0.54),
    Vec3::new(0.30, 0.50, 0.82),
];

#[test]
fn triangulated_box_matches_the_analytic_box() {
    let surface = import_stl(UNIT_BOX, StlOptions::default()).unwrap();
    let bounds = AxisAlignedBox::new(Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    assert_eq!(surface.bounds(), bounds);
    assert_eq!(surface.vertices().len(), 8);
    assert_eq!(surface.facets().len(), 12);
    assert_eq!(surface.patch_names(), PATCHES);
    let mesh = generate_surface_mesh(&surface, &SITES).unwrap();
    let reference = generate_mesh(bounds, &SITES).unwrap();
    assert_eq!(mesh, generate_surface_mesh(&surface, &SITES).unwrap());
    assert_eq!(mesh.metrics().domain_volume, 1.0);
    for (actual, reference) in mesh.metrics().cells.iter().zip(&reference.metrics().cells) {
        assert!((actual.volume - reference.volume).abs() < 1e-12);
        assert!((actual.centroid - reference.centroid).length() < 1e-12);
        assert_eq!(actual.neighbour_count, reference.neighbour_count);
    }
    let mut areas = [0.0; 6];
    for (face, metrics) in mesh.metrics().faces.iter().enumerate() {
        if let FaceAdjacency::Surface { patch, .. } = mesh.mesh().face_adjacency(face).unwrap() {
            areas[patch.0] += metrics.area;
        }
    }
    assert!(areas.iter().all(|area| (area - 1.0).abs() < 1e-12));
}

#[test]
fn binary_detection_and_explicit_patches_preserve_coplanar_labels() {
    let ascii = import_stl(UNIT_BOX, StlOptions::default()).unwrap();
    let triangles = triangle_points(&ascii);
    let bytes = binary(&triangles);
    let unnamed = import_stl(&bytes, StlOptions::default()).unwrap();
    assert_eq!(unnamed.patch_names(), ["surface"]);
    let mut names: Vec<_> = ascii
        .facets()
        .iter()
        .map(|facet| ascii.patch_names()[facet.triangle.patch.0].clone())
        .collect();
    names[0] = "lower inlet".to_owned();
    names[1] = "upper inlet".to_owned();
    let surface = import_stl(
        &bytes,
        StlOptions {
            facet_patches: Some(&names),
            facet_limit: None,
        },
    )
    .unwrap();
    assert_eq!(&surface.patch_names()[..2], ["lower inlet", "upper inlet"]);
    assert_ne!(
        surface.facets()[0].triangle.patch,
        surface.facets()[1].triangle.patch
    );
    let mesh = generate_surface_mesh(&surface, &SITES).unwrap();
    let mut areas = [0.0; 2];
    for (face, metrics) in mesh.metrics().faces.iter().enumerate() {
        if let FaceAdjacency::Surface { patch, .. } = mesh.mesh().face_adjacency(face).unwrap()
            && let Some(area) = areas.get_mut(patch.0)
        {
            *area += metrics.area;
        }
    }
    assert!(areas.iter().all(|area| (area - 0.5).abs() < 1e-12));
    let mut msh = Vec::new();
    write_gmsh_domain(MeshDomain::Surface(&surface), mesh.mesh(), &mut msh).unwrap();
    let text = String::from_utf8(msh).unwrap();
    assert!(text.contains("2 1 \"lower inlet\""));
    assert!(text.contains("2 2 \"upper inlet\""));
}

#[test]
fn skewed_box_clips_against_the_surface_and_exports() {
    let surface = import_stl(SKEWED_BOX, StlOptions::default()).unwrap();
    let exterior = Vec3::new(0.1, 0.9, 0.5);
    assert!(surface.bounds().contains_interior(exterior));
    assert!(matches!(
        generate_surface_mesh(&surface, &[exterior]),
        Err(GenerationError::Clip(VoronoiError::SiteOutside {
            index: 0
        }))
    ));
    let result = generate_surface_mesh(&surface, &SITES).unwrap();
    assert_eq!(result.mesh().cell_count(), 4);
    let domain = MeshDomain::Surface(&surface);
    let metrics = validate_domain(domain, result.mesh()).unwrap();
    assert!((metrics.domain_volume - 1.0).abs() < 1e-12);
    assert!((metrics.total_volume - 1.0).abs() < 1e-12);
    assert!(metrics.cells.iter().all(|cell| cell.volume > 0.0));
    let mut vtu = Vec::new();
    write_vtu_domain(domain, result.mesh(), &mut vtu).unwrap();
    assert!(
        String::from_utf8(vtu)
            .unwrap()
            .contains("NumberOfCells=\"4\"")
    );
    write_gmsh_domain(domain, result.mesh(), &mut std::io::sink()).unwrap();
}

#[test]
fn invalid_stl_shells_are_rejected() {
    let surface = import_stl(UNIT_BOX, StlOptions::default()).unwrap();
    let triangles = triangle_points(&surface);
    let parse = |triangles: &[[Vec3; 3]]| import_stl(&binary(triangles), StlOptions::default());
    assert!(matches!(
        parse(&triangles[1..]),
        Err(SurfaceError::OpenEdge { .. })
    ));
    let mut modified = triangles.clone();
    modified[0].swap(1, 2);
    assert!(matches!(
        parse(&modified),
        Err(SurfaceError::WindingOrNonManifold { .. })
    ));
    let mut modified = triangles.clone();
    for triangle in &mut modified {
        triangle.swap(1, 2);
    }
    assert_eq!(parse(&modified), Err(SurfaceError::InvalidVolume));
    let mut modified = triangles.clone();
    modified[0][1] = modified[0][0];
    assert_eq!(parse(&modified), Err(SurfaceError::DegenerateFacet(0)));
    let mut modified = triangles.clone();
    modified.push(modified[0]);
    assert_eq!(parse(&modified), Err(SurfaceError::DuplicateFacet(12)));
    let mut modified = triangles.clone();
    for point in modified.iter_mut().flatten() {
        if *point == Vec3::new(1.0, 1.0, 1.0) {
            *point = Vec3::new(0.75, 0.75, 0.75);
        }
    }
    assert!(matches!(
        parse(&modified),
        Err(SurfaceError::NonConvexFacet(_))
    ));
    let bytes = binary(&triangles);
    assert!(matches!(
        import_stl(&bytes[..bytes.len() - 1], StlOptions::default()),
        Err(SurfaceError::InvalidStl { .. })
    ));
    assert_eq!(
        import_stl(
            &bytes,
            StlOptions {
                facet_patches: Some(&[]),
                facet_limit: None
            }
        ),
        Err(SurfaceError::PatchCount {
            facets: 12,
            names: 0
        })
    );
    assert_eq!(
        import_stl(
            &bytes,
            StlOptions {
                facet_patches: None,
                facet_limit: Some(11)
            }
        ),
        Err(SurfaceError::FacetLimit { limit: 11 })
    );
}

fn triangle_points(surface: &ConvexSurface) -> Vec<[Vec3; 3]> {
    surface
        .facets()
        .iter()
        .map(|facet| facet.triangle.vertices.map(|id| surface.vertices()[id]))
        .collect()
}

fn binary(triangles: &[[Vec3; 3]]) -> Vec<u8> {
    let mut bytes = vec![0; 80];
    // A header starting with "solid" still has binary records.
    bytes[..5].copy_from_slice(b"solid");
    bytes.extend(u32::try_from(triangles.len()).unwrap().to_le_bytes());
    for points in triangles {
        bytes.extend([0; 12]);
        for &Vec3 { x, y, z } in points {
            for coordinate in [x, y, z] {
                bytes.extend((coordinate as f32).to_le_bytes());
            }
        }
        bytes.extend([0; 2]);
    }
    bytes
}
