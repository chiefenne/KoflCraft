use std::collections::BTreeMap;

use mesh_algorithms::channel::straight_channel;
use mesh_core::geometry::Vec3;
use mesh_io::channel::write_gmsh;
use mesh_quality::validate_mesh;

/// Read the emitted MSH elements independently. Match every external facet to
/// one volume element, every internal facet to two, and measure signed volumes.
#[test]
fn native_channel_prisms_and_core_subdivision_share_the_tagged_boundary() {
    for (stations, layers, growth) in [(1, 1, 1.0), (3, 2, 1.4)] {
        let channel =
            straight_channel(Vec3::new(3.0, 2.0, 1.0), stations, layers, 0.05, growth).unwrap();
        let metrics = validate_mesh(channel.bounds(), channel.mesh()).unwrap();
        let mut output = Vec::new();
        write_gmsh(&channel, &mut output).unwrap();
        let file = std::str::from_utf8(&output).unwrap();
        assert!(file.contains("2 1 \"inlet\"\n2 2 \"outlet\"\n2 3 \"walls\"\n3 4 \"fluid\""));
        let mut points = vec![Vec3::ZERO]; // MSH uses one-based node IDs.
        let nodes = file
            .split_once("$Nodes\n")
            .unwrap()
            .1
            .split_once("$EndNodes")
            .unwrap()
            .0;
        for line in nodes.lines().skip(1) {
            let mut words = line.split_whitespace();
            assert_eq!(
                words.next().unwrap().parse::<usize>().unwrap(),
                points.len()
            );
            let mut coord = || words.next().unwrap().parse().unwrap();
            points.push(Vec3::new(coord(), coord(), coord()));
        }
        let mut uses = BTreeMap::new();
        let mut boundary = BTreeMap::new();
        let mut areas = [0.; 3];
        let mut volumes = vec![0.0; channel.cells().len()];
        let mut prisms = 0;
        let stream = file
            .split_once("$Elements\n")
            .unwrap()
            .1
            .split_once("$EndElements")
            .unwrap()
            .0;
        for line in stream.lines().skip(1) {
            let e: Vec<usize> = line
                .split_whitespace()
                .map(|s| s.parse().unwrap())
                .collect();
            let n = &e[5..];
            let phys = e[3];
            if matches!(e[1], 2 | 3) {
                let mut key = n.to_vec();
                key.sort_unstable();
                assert!(boundary.insert(key, phys).is_none());
                let origin = points[n[0]];
                let normal = n[1..].windows(2).fold(Vec3::ZERO, |sum, p| {
                    sum + (points[p[0]] - origin).cross(points[p[1]] - origin) * 0.5
                });
                areas[phys - 1] += normal.length();
                match phys {
                    1 | 2 => {
                        assert!(
                            n.iter()
                                .all(|&id| points[id].x == if phys == 1 { 0.0 } else { 3.0 })
                        );
                        assert!(if phys == 1 {
                            normal.x < 0.0
                        } else {
                            normal.x > 0.0
                        });
                    }
                    3 => assert!(normal.x.abs() < 1e-12),
                    _ => panic!("unexpected boundary group"),
                }
                continue;
            }
            assert_eq!(phys, 4);
            let (facets, tets): (&[&[usize]], &[[usize; 4]]) = match e[1] {
                4 => (
                    &[&[0, 1, 2], &[0, 1, 3], &[0, 2, 3], &[1, 2, 3]],
                    &[[0, 1, 2, 3]],
                ),
                6 => {
                    prisms += 1;
                    (
                        &[
                            &[0, 1, 2],
                            &[3, 4, 5],
                            &[0, 1, 3, 4],
                            &[1, 2, 4, 5],
                            &[0, 2, 3, 5],
                        ],
                        &[[0, 1, 2, 3], [1, 4, 2, 3], [2, 4, 5, 3]],
                    )
                }
                7 => (
                    &[
                        &[0, 1, 2, 3],
                        &[0, 1, 4],
                        &[1, 2, 4],
                        &[2, 3, 4],
                        &[0, 3, 4],
                    ],
                    &[[0, 1, 2, 4], [0, 2, 3, 4]],
                ),
                kind => panic!("unexpected volume element {kind}"),
            };
            for facet in facets {
                let mut key: Vec<_> = facet.iter().map(|&i| n[i]).collect();
                key.sort_unstable();
                *uses.entry(key).or_insert(0) += 1;
            }
            for tet in tets {
                let [a, b, c, d] = tet.map(|i| points[n[i]]);
                let volume = (b - a).dot((c - a).cross(d - a)) / 6.0;
                assert!(volume > 0.0);
                volumes[e[4] - 1] += volume;
            }
        }
        assert_eq!(prisms, channel.cells().len() - stations);
        for (key, count) in uses {
            assert_eq!(
                count,
                if boundary.remove(&key).is_some() {
                    1
                } else {
                    2
                }
            );
        }
        assert!(boundary.is_empty());
        for (area, expected) in areas.into_iter().zip([2., 2., 18.]) {
            assert!((area - expected).abs() < 1e-12);
        }
        for (volume, cell) in volumes.iter().zip(metrics.cells) {
            assert!((volume / cell.volume - 1.0).abs() < 1e-12);
        }
    }
}
