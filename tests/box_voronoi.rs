use std::collections::BTreeMap;

use koflcraft::certification::validate_mesh;
use koflcraft::clipping::{ClipError, ConvexPolyhedron, FaceSource};
use koflcraft::geometry::{AxisAlignedBox, Plane, Vec3};
use koflcraft::mesh_export::{GmshVersion, write_gmsh};
use koflcraft::voronoi::{RELATIVE_CLIP_TOLERANCE, VoronoiError, clipped_voronoi};
use koflcraft::{GenerationError, generate_mesh};

#[path = "../test-fixtures/sites.rs"]
mod fixtures;

fn unit_box() -> AxisAlignedBox {
    AxisAlignedBox::new(Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap()
}

// Independently integrate outward-wound triangles about a point inside the cell.
fn volume_and_closedness(cell: &ConvexPolyhedron, inside: Vec3) -> f64 {
    let points = cell.vertices();
    let mut edges = BTreeMap::<_, (usize, i32)>::new();
    let mut volume = 0.0;
    for face in cell.faces() {
        let ids = face.vertices();
        assert!(ids.len() >= 3);
        let a = points[ids[0]] - inside;
        let mut face_volume = 0.0;
        for pair in ids[1..].windows(2) {
            let b = points[pair[0]] - inside;
            let c = points[pair[1]] - inside;
            face_volume += a.dot(b.cross(c)) / 6.0;
        }
        assert!(
            face_volume > 0.0,
            "face must be outward-wound with positive area"
        );
        volume += face_volume;
        for (&a, &b) in ids.iter().zip(ids.iter().cycle().skip(1)) {
            let (uses, balance) = edges.entry((a.min(b), a.max(b))).or_default();
            *uses += 1;
            *balance += if a < b { 1 } else { -1 };
        }
    }
    assert!(edges.values().all(|count| *count == (2, 0)));
    volume
}

#[test]
fn box_clipping_handles_crossings_coplanar_edges_and_tangency() {
    let cube = ConvexPolyhedron::from_box(unit_box());
    let plane = |point, normal| Plane::new(point, normal).unwrap();
    let pos_x = Vec3::new(1.0, 0.0, 0.0);
    let cut = |point, normal| {
        cube.clone()
            .clip(
                plane(point, normal),
                FaceSource::Cut,
                RELATIVE_CLIP_TOLERANCE,
            )
            .unwrap()
    };

    let measured = volume_and_closedness(&cube, Vec3::new(0.5, 0.5, 0.5));
    assert!((measured - 1.0).abs() < 1e-14);
    assert_eq!(cut(pos_x, pos_x), Some(cube.clone()));
    assert_eq!(cut(Vec3::ZERO, pos_x), None);
    assert_eq!(
        cube.clone()
            .clip(plane(pos_x, pos_x), FaceSource::Cut, f64::NAN),
        Err(ClipError::InvalidTolerance)
    );

    let cases = [
        (Vec3::new(0.25, 0.0, 0.0), pos_x, 4, 0.25),
        (pos_x, Vec3::new(1.0, 1.0, 0.0), 4, 0.5),
        (pos_x, Vec3::new(1.0, 1.0, 1.0), 3, 1.0 / 6.0),
    ];
    for (point, normal, cap_vertices, volume) in cases {
        let clipped = cut(point, normal).unwrap();
        let cap = clipped.faces().last().unwrap();
        assert_eq!(cap.source(), FaceSource::Cut);
        assert_eq!(cap.vertices().len(), cap_vertices);
        for &id in cap.vertices() {
            assert!(
                plane(point, normal)
                    .signed_distance(clipped.vertices()[id])
                    .abs()
                    < 1e-14
            );
        }
        let measured = volume_and_closedness(&clipped, Vec3::new(0.125, 0.125, 0.125));
        assert!((measured - volume).abs() < 1e-14);
    }
}

#[test]
fn rectangular_partition_survives_rescaling_and_translation() {
    for (origin, scale) in [(0.0, 1e-30), (1e9, 2.0), (-1e12, 1e12)] {
        let min = Vec3::new(origin, -origin, origin);
        let max = min + Vec3::new(4.0, 2.0, 1.0) * scale;
        let domain = AxisAlignedBox::new(min, max).unwrap();
        let sites = [
            min + Vec3::new(0.5, 1.0, 0.5) * scale,
            min + Vec3::new(1.5, 1.0, 0.5) * scale,
        ];
        let cells = clipped_voronoi(domain, &sites).unwrap();
        for (id, (cell, unscaled_volume)) in cells.iter().zip([2.0, 6.0]).enumerate() {
            let expected = unscaled_volume * scale.powi(3);
            let actual = volume_and_closedness(cell, sites[id]);
            assert!((actual / expected - 1.0).abs() < 1e-13);
            let shared = cell
                .faces()
                .iter()
                .find(|face| face.source() == FaceSource::Bisector { other_site: 1 - id })
                .unwrap();
            assert!(
                shared
                    .vertices()
                    .iter()
                    .all(|&v| { (cell.vertices()[v].x - (origin + scale)).abs() <= 1e-13 * scale })
            );
        }
    }
}

#[test]
fn deterministic_lattices_and_nonuniform_sites_partition_the_cube() {
    for dimensions in fixtures::LATTICES {
        let count = dimensions.into_iter().product::<u32>() as usize;
        for sites in [fixtures::lattice(dimensions), fixtures::clustered(count)] {
            let cells = clipped_voronoi(unit_box(), &sites).unwrap();
            assert_eq!(cells.len(), count);
            assert_eq!(cells, clipped_voronoi(unit_box(), &sites).unwrap());
            let mut volume = 0.0;
            for (cell, &site) in cells.iter().zip(&sites) {
                volume += volume_and_closedness(cell, site);
                for &vertex in cell.vertices() {
                    for axis in [vertex.x, vertex.y, vertex.z] {
                        assert!((0.0..=1.0).contains(&axis));
                    }
                    let own = vertex - site;
                    for &other in &sites {
                        let alternative = vertex - other;
                        assert!(own.dot(own) <= alternative.dot(alternative) + 1e-12);
                    }
                }
            }
            assert!(
                (volume - 1.0).abs() < 1e-10,
                "{count} sites measured {volume}"
            );
            let generated = generate_mesh(unit_box(), &sites).unwrap();
            assert_eq!(generated.mesh().cell_count(), count);
            assert_eq!(generated, generate_mesh(unit_box(), &sites).unwrap());
            assert_eq!(
                *generated.metrics(),
                validate_mesh(unit_box(), generated.mesh()).unwrap()
            );
            assert!((generated.metrics().total_volume - 1.0).abs() < 1e-10);
            assert!(
                generated
                    .metrics()
                    .cells
                    .iter()
                    .all(|cell| cell.volume > 0.0)
            );
        }
    }
}

#[test]
fn detects_invalid_domains_and_sites() {
    assert!(AxisAlignedBox::new(Vec3::ZERO, Vec3::new(0.0, 1.0, 1.0)).is_none());
    assert!(Plane::new(Vec3::ZERO, Vec3::ZERO).is_none());
    assert_eq!(clipped_voronoi(unit_box(), &[]), Err(VoronoiError::NoSites));
    for point in [Vec3::ZERO, Vec3::new(0.5, 0.5, f64::NAN)] {
        assert_eq!(
            clipped_voronoi(unit_box(), &[point]),
            Err(VoronoiError::SiteOutside { index: 0 })
        );
    }
    let site = Vec3::new(0.5, 0.5, 0.5);
    for other in [site, site + Vec3::new(f64::EPSILON, 0.0, 0.0)] {
        assert_eq!(
            clipped_voronoi(unit_box(), &[site, other]),
            Err(VoronoiError::SitesTooClose {
                first: 0,
                second: 1
            })
        );
    }
    let thin = AxisAlignedBox::new(Vec3::ZERO, Vec3::new(1.0, 1.0, f64::EPSILON)).unwrap();
    assert_eq!(
        clipped_voronoi(thin, &[site]),
        Err(VoronoiError::ThinDomain)
    );
}

#[test]
fn nearby_distinct_sites_keep_separate_closed_cells() {
    let sites = [Vec3::new(0.5, 0.5, 0.5), Vec3::new(0.500_000_01, 0.5, 0.5)];
    let mesh = generate_mesh(unit_box(), &sites).unwrap();
    assert_eq!(mesh.mesh().cell_count(), sites.len());
    assert!(mesh.metrics().cells.iter().all(|cell| cell.volume > 0.0));
    assert!((mesh.metrics().total_volume - 1.0).abs() < 1e-12);
}

#[test]
fn incompatible_near_degenerate_face_loops_are_rejected() {
    let sites = [
        Vec3::new(0.5, 0.5, 0.5),
        Vec3::new(0.500_000_01, 0.500_000_01, 0.5),
        Vec3::new(0.25, 0.25, 0.75),
        Vec3::new(0.75, 0.75, 0.25),
    ];
    assert!(matches!(
        generate_mesh(unit_box(), &sites),
        Err(GenerationError::UnmatchedFace { .. })
    ));
}

#[test]
fn seeded_near_wall_and_nearly_coincident_sites_remain_resolvable() {
    let mut seed = fixtures::SEED;
    for axis in 0..3 {
        for high in [false, true] {
            for case in 0..8 {
                let mut position = [0.0; 3].map(|_| 0.2 + 0.6 * fixtures::sample(&mut seed));
                let gap = 1e-7 * (1.0 + fixtures::sample(&mut seed));
                let clearance = 1e-9 * (1.0 + fixtures::sample(&mut seed));
                let distances = [clearance, clearance + gap, 0.7];
                let sites = distances.map(|distance| {
                    position[axis] = if high { 1.0 - distance } else { distance };
                    Vec3::new(position[0], position[1], position[2])
                });
                let mesh = generate_mesh(unit_box(), &sites).unwrap_or_else(|error| {
                    panic!("axis {axis}, high wall {high}, case {case}: {error}")
                });
                assert_eq!(mesh.mesh().cell_count(), 3);
                let metrics = mesh.metrics();
                assert!((metrics.cells[0].volume - (clearance + 0.5 * gap)).abs() < 2e-14);
                assert!(metrics.cells.iter().all(|cell| cell.volume > 0.0));
                assert!((metrics.total_volume - 1.0).abs() < 1e-12);
                write_gmsh(
                    unit_box(),
                    mesh.mesh(),
                    GmshVersion::V4_1,
                    &mut std::io::sink(),
                )
                .unwrap();
            }
        }
    }
}

/// The bisector neighbours of a cell, with each face's area and the cell's
/// volume, independent of the order its faces were cut in.
fn signature(cell: &ConvexPolyhedron, site: Vec3) -> (f64, Vec<(String, f64)>) {
    let points = cell.vertices();
    let mut faces: Vec<(String, f64)> = cell
        .faces()
        .iter()
        .map(|face| {
            let ids = face.vertices();
            let a = points[ids[0]];
            let area: f64 = ids[1..]
                .windows(2)
                .map(|p| (points[p[0]] - a).cross(points[p[1]] - a).length() / 2.0)
                .sum();
            (format!("{:?}", face.source()), area)
        })
        .collect();
    faces.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
    (volume_and_closedness(cell, site), faces)
}

/// Keep rule (SW-18). Invariant (ALG-04, ALG-13): clipping each cell by its
/// nearest sites until the security radius proves it complete gives the
/// same cells as the exhaustive all-pairs reference: the same bisector
/// neighbours and boundary faces, face areas and volumes to round-off, on
/// degenerate lattices (many cospherical sites), clustered sites and 2,000
/// scattered sites with a 100:1 density range. A security test that stops
/// too early drops a neighbour and fails here; the mesh validation would
/// catch only a cell that no longer tiles the box.
#[test]
fn nearest_site_clipping_equals_the_all_pairs_reference() {
    use koflcraft::voronoi::{Strategy, clipped_voronoi_with};
    let mut cases: Vec<Vec<Vec3>> = fixtures::LATTICES
        .iter()
        .map(|&l| fixtures::lattice(l))
        .collect();
    cases.push(fixtures::lattice([8, 8, 8]));
    cases.push(fixtures::clustered(256));
    let mut seed = fixtures::SEED;
    cases.push(
        (0..2000)
            .map(|_| {
                // dense near the x = 0 face: x = u^3 spans two decades
                let u = fixtures::sample(&mut seed);
                Vec3::new(
                    (u * u * u).clamp(1e-4, 1.0 - 1e-4),
                    fixtures::sample(&mut seed).clamp(1e-4, 1.0 - 1e-4),
                    fixtures::sample(&mut seed).clamp(1e-4, 1.0 - 1e-4),
                )
            })
            .collect(),
    );
    for sites in &cases {
        let fast = clipped_voronoi_with(unit_box(), sites, Strategy::Nearest).unwrap();
        let slow = clipped_voronoi_with(unit_box(), sites, Strategy::AllPairs).unwrap();
        assert_eq!(fast.len(), slow.len());
        for (i, (a, b)) in fast.iter().zip(&slow).enumerate() {
            let (va, fa) = signature(a, sites[i]);
            let (vb, fb) = signature(b, sites[i]);
            assert!(
                (va - vb).abs() <= 1e-12,
                "{} sites, cell {i}: volume {va} vs {vb}",
                sites.len()
            );
            let names = |f: &[(String, f64)]| f.iter().map(|x| x.0.clone()).collect::<Vec<_>>();
            assert_eq!(names(&fa), names(&fb), "{} sites, cell {i}", sites.len());
            for (x, y) in fa.iter().zip(&fb) {
                assert!(
                    (x.1 - y.1).abs() <= 1e-12,
                    "{} sites, cell {i}, face {}",
                    sites.len(),
                    x.0
                );
            }
        }
        generate_mesh(unit_box(), sites).unwrap();
    }
}
