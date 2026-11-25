use std::collections::HashMap;

use mesh_core::block::{BlockSpec, corner_slot};

use super::grid::HEX_FACES;
use super::templates::{AirfoilCgrid, Naca4, Section, airfoil_cgrid, channel, elbow, pipe_ogrid};
use super::*;

/// Every hex face is shared by two hexes in opposite winding or is one
/// boundary quad, wound outward.
fn assert_conforming(g: &BlockGrid) {
    let mut seen: HashMap<Vec<usize>, Vec<Vec<usize>>> = HashMap::new();
    for h in &g.hexes {
        for f in HEX_FACES {
            let ring: Vec<usize> = f.iter().map(|&s| h[s]).collect();
            let mut k = ring.clone();
            k.sort_unstable();
            seen.entry(k).or_default().push(ring);
        }
    }
    let mut boundary: HashMap<Vec<usize>, Vec<usize>> = HashMap::new();
    for (q, _, _) in &g.boundary {
        let mut k = q.to_vec();
        k.sort_unstable();
        assert!(
            boundary.insert(k, q.to_vec()).is_none(),
            "boundary quad twice"
        );
    }
    for (k, rings) in &seen {
        match rings.len() {
            1 => assert_eq!(
                boundary.get(k),
                Some(&rings[0]),
                "unpatched or misoriented boundary face"
            ),
            2 => {
                assert!(!boundary.contains_key(k));
                let (a, b) = (&rings[0], &rings[1]);
                let s = b.iter().position(|&v| v == a[0]).unwrap();
                assert!(
                    (0..4).all(|i| b[(s + 4 - i) % 4] == a[i]),
                    "same winding twice"
                );
            }
            n => panic!("a face is shared by {n} cells"),
        }
    }
    assert_eq!(
        boundary.len(),
        seen.values().filter(|r| r.len() == 1).count()
    );
}

fn total(g: &BlockGrid) -> f64 {
    g.volumes().iter().sum()
}

/// Invariant: a box channel is reproduced exactly (volume, unit scaled
/// Jacobians) and every face conforms. A wrong TFI weight or hex winding
/// breaks it.
#[test]
fn channel_is_exact() {
    let t = channel([2.0, 1.0, 0.5], [6, 4, 3], None).unwrap();
    let g = build(&t.topology).unwrap();
    assert_eq!(g.hexes.len(), 72);
    assert_eq!(g.points.len(), 7 * 5 * 4);
    assert_conforming(&g);
    assert!((total(&g) - 1.0).abs() < 1e-13);
    assert!((g.boundary_volume() - 1.0).abs() < 1e-13);
    for q in g.quality() {
        assert!((q.min_scaled_jacobian - 1.0).abs() < 1e-12);
    }
}

/// Invariant: two-sided clustering keeps a tensor grid (TFI of straight
/// edges with equal laws) and reaches the requested wall spacing within
/// Vinokur's slope approximation.
#[test]
fn channel_clusters_at_walls() {
    let t = channel([1.0, 0.2, 0.1], [4, 20, 10], Some(0.002)).unwrap();
    let g = build(&t.topology).unwrap();
    assert_conforming(&g);
    assert!((total(&g) - 0.02).abs() < 1e-15);
    for q in g.quality() {
        assert!((q.min_scaled_jacobian - 1.0).abs() < 1e-9);
    }
    for e in &g.edges {
        if e.cells == 20 {
            assert!((e.first_spacing / 0.002 - 1.0).abs() < 0.25, "{e:?}");
        }
    }
}

/// Invariant: the pipe O-grid conforms, its wall nodes lie on the
/// cylinder, and its volume equals the inscribed polygon prism exactly and
/// πR²L within the chord error. A wrong face orientation map breaks the
/// conformity or the node count.
#[test]
fn pipe_ogrid_matches_its_cylinder() {
    let (r, l, core) = (0.5, 3.0, 6);
    let t = pipe_ogrid(r, l, 0.5, [core, 5, 8], Some(0.01)).unwrap();
    let g = build(&t.topology).unwrap();
    assert_eq!(g.hexes.len(), 8 * (core * core + 4 * 5 * core));
    assert_conforming(&g);
    let n = 4 * core;
    let polygon = 0.5 * r * r * n as f64 * (std::f64::consts::TAU / n as f64).sin();
    assert!((total(&g) / (polygon * l) - 1.0).abs() < 1e-12);
    assert!((g.boundary_volume() / total(&g) - 1.0).abs() < 1e-12);
    assert!((total(&g) / t.analytic_volume - 1.0).abs() < 0.02);
    let walls = 2;
    for (q, p, _) in &g.boundary {
        if *p == walls {
            for &v in q {
                let pt = g.points[v];
                assert!((pt.y.hypot(pt.z) - r).abs() < 1e-12);
            }
        }
    }
    for q in g.quality() {
        assert!(q.min_jacobian > 0.0 && q.min_scaled_jacobian > 0.3, "{q:?}");
    }
    for e in g.edges.iter().filter(|e| e.cells == 5) {
        assert!((e.last_spacing / 0.01 - 1.0).abs() < 1e-9, "{e:?}");
    }
}

/// A block with its local axes permuted cyclically (still right-handed).
fn permuted(b: &BlockSpec, perm: [usize; 3]) -> BlockSpec {
    let mut vertices = [0; 8];
    for (slot, v) in vertices.iter_mut().enumerate() {
        let a = slot_corner(slot);
        let mut c = [0; 3];
        for d in 0..3 {
            c[perm[d]] = a[d];
        }
        *v = b.vertices[corner_slot(c[0], c[1], c[2])];
    }
    BlockSpec {
        vertices,
        cells: std::array::from_fn(|d| b.cells[perm[d]]),
        domain: b.domain.clone(),
    }
}

/// Invariant: the grid does not depend on how each block's local frame is
/// chosen. Rotated and permuted frames give the same nodes and volume; this
/// exercises all face orientation maps of the shared-face reader.
#[test]
fn local_frames_do_not_change_the_grid() {
    let t = pipe_ogrid(1.0, 2.0, 0.55, [4, 3, 5], None).unwrap();
    let reference = build(&t.topology).unwrap();
    let mut topo = t.topology.clone();
    topo.blocks[1] = permuted(&topo.blocks[1], [1, 2, 0]);
    topo.blocks[2] = permuted(&topo.blocks[2], [2, 0, 1]);
    // rotate block 3 a quarter turn about k
    let b = topo.blocks[3].clone();
    let v = b.vertices;
    topo.blocks[3] = BlockSpec {
        vertices: [v[1], v[2], v[3], v[0], v[5], v[6], v[7], v[4]],
        cells: [b.cells[1], b.cells[0], b.cells[2]],
        domain: b.domain,
    };
    let g = build(&topo).unwrap();
    assert_conforming(&g);
    assert_eq!(g.points.len(), reference.points.len());
    assert!((total(&g) / total(&reference) - 1.0).abs() < 1e-13);
}

/// Invariant: an inconsistent chain is refused before any grid is built,
/// naming the edge and both blocks (BS-05).
#[test]
fn inconsistent_counts_are_named() {
    let mut t = pipe_ogrid(1.0, 2.0, 0.55, [4, 3, 5], None)
        .unwrap()
        .topology;
    t.blocks[1].cells[1] = 5;
    let issues = check(&t).unwrap_err();
    assert!(
        issues
            .iter()
            .any(|i| i.kind == "count_chain" && i.message.contains("block"))
    );
}

/// Invariant: a boundary face outside every patch and a left-handed block
/// are refused, never meshed.
#[test]
fn open_faces_and_left_handed_blocks_are_refused() {
    let mut t = channel([1.0, 1.0, 1.0], [2, 2, 2], None).unwrap().topology;
    t.patches[2].faces.pop();
    assert!(check(&t).unwrap_err().iter().any(|i| i.kind == "patch"));
    let mut t = channel([1.0, 1.0, 1.0], [2, 2, 2], None).unwrap().topology;
    t.blocks[0].vertices = [0, 3, 2, 1, 4, 7, 6, 5];
    assert!(check(&t).unwrap_err().iter().any(|i| i.kind == "block"));
}

mod elliptic {
    use super::super::smooth::{SmoothOptions, smooth};
    use super::*;

    fn walls(t: &mesh_core::block::BlockTopology) -> SmoothOptions {
        SmoothOptions {
            wall_patches: t
                .patches
                .iter()
                .enumerate()
                .filter(|(_, p)| p.name == "walls")
                .map(|(i, _)| i)
                .collect(),
            ..SmoothOptions::default()
        }
    }

    /// Invariant: a channel of straight, orthogonal, stretched lines
    /// already solves the system with its distribution controls, so
    /// smoothing leaves it in place. A sign error in a control or in the
    /// SOR update moves it.
    #[test]
    fn smoothing_keeps_a_stretched_channel() {
        let t = channel([1.0, 0.2, 0.1], [6, 12, 8], Some(0.002)).unwrap();
        let mut g = build(&t.topology).unwrap();
        let tfi = g.points.clone();
        let r = smooth(&mut g, &walls(&t.topology)).unwrap();
        assert!(r.converged && r.wall_faces == 4, "{r:?}");
        for (a, b) in g.points.iter().zip(&tfi) {
            assert!((*a - *b).length() < 1e-9, "{a:?} {b:?}");
        }
        assert!(r.wall_angle_max_after < 1e-6);
    }

    /// Invariant: on the pipe O-grid, nodes on shared faces and inner
    /// block edges move, boundary nodes do not (BS-08), the first grid
    /// line leaves the wall at a right angle with the first spacing kept,
    /// cells stay valid and conforming, and the volume is unchanged.
    #[test]
    fn pipe_walls_become_orthogonal() {
        let t = pipe_ogrid(0.05, 0.2, 0.55, [8, 10, 6], Some(0.0005)).unwrap();
        let mut g = build(&t.topology).unwrap();
        let tfi = g.points.clone();
        let volume = total(&g);
        let r = smooth(&mut g, &walls(&t.topology)).unwrap();
        assert!(r.converged, "{r:?}");
        assert!(r.face_nodes > 0 && r.edge_nodes == 4 * 5, "{r:?}");
        assert!(
            r.wall_angle_max_before > 2.0 && r.wall_angle_max_after < 0.5,
            "{r:?}"
        );
        assert!(r.first_spacing_change < 1e-4, "{r:?}");
        for &(q, _, _) in &g.boundary {
            for v in q {
                assert_eq!(g.points[v], tfi[v]);
            }
        }
        let moved = g.points.iter().zip(&tfi).filter(|(a, b)| a != b).count();
        assert!(moved > r.block_nodes, "shared-face nodes must move too");
        assert_conforming(&g);
        assert!((total(&g) / volume - 1.0).abs() < 1e-12);
        for q in g.quality() {
            assert!(
                q.min_jacobian > 0.0 && q.min_scaled_jacobian > 0.69,
                "{q:?}"
            );
        }
    }

    fn patch(t: &mesh_core::block::BlockTopology, name: &str) -> usize {
        t.patches.iter().position(|p| p.name == name).unwrap()
    }

    /// Invariant: on the elbow, the interior block vertices and the bend
    /// joins move (frozen TFI residual, so only the wall controls act),
    /// the wall nodes stay on the pipe, the wall angle drops and the worst
    /// scaled Jacobian stays within 1 % of the TFI grid's (measured: 0.704
    /// against 0.706; FoehnCraft max non-orthogonality 23.45° against
    /// 23.7° for TFI and 34.5° for the MB1b smoother).
    #[test]
    fn elbow_smoothing_is_not_worse_than_tfi() {
        let (r, big) = (0.5, 1.5);
        let t = elbow(r, big, 90.0, 1.0, 1.0, 0.5, [6, 6, 36], Some(0.01)).unwrap();
        let mut g = build(&t.topology).unwrap();
        let worst = |g: &BlockGrid| {
            g.quality()
                .iter()
                .map(|q| q.min_scaled_jacobian)
                .fold(f64::INFINITY, f64::min)
        };
        let before = worst(&g);
        let res = smooth(&mut g, &walls(&t.topology)).unwrap();
        assert!(res.converged, "{res:?}");
        assert!(res.vertex_nodes > 0 && res.face_nodes > 0, "{res:?}");
        assert!(
            res.wall_angle_max_after < 0.5 * res.wall_angle_max_before,
            "{res:?}"
        );
        assert!(elbow_wall_error(&g, r, big, 1.0) < 1e-12);
        assert_conforming(&g);
        assert!(worst(&g) >= 0.99 * before, "{} < {before}", worst(&g));
    }

    /// Invariant: on projection surfaces (pipe cylinder, elbow torus) the
    /// block-edge wall nodes slide and the face nodes follow; they stay on
    /// the surface, the smoother converges, and the grid is no worse than
    /// with fixed walls (already orthogonal to 0.4° to 0.6°).
    #[test]
    fn surface_walls_slide_and_stay() {
        let worst = |g: &BlockGrid| {
            g.quality()
                .iter()
                .map(|q| q.min_scaled_jacobian)
                .fold(f64::INFINITY, f64::min)
        };
        let (r, big) = (0.5, 1.5);
        let cases = [
            elbow(r, big, 90.0, 1.0, 1.0, 0.5, [6, 6, 36], Some(0.01)).unwrap(),
            pipe_ogrid(0.05, 0.2, 0.55, [8, 10, 6], Some(0.0005)).unwrap(),
        ];
        for (c, t) in cases.iter().enumerate() {
            let mut fixed = build(&t.topology).unwrap();
            let mut slid = fixed.clone();
            let a = smooth(&mut fixed, &walls(&t.topology)).unwrap();
            let o = SmoothOptions {
                slide_walls: true,
                ..walls(&t.topology)
            };
            let b = smooth(&mut slid, &o).unwrap();
            assert!(b.converged && b.contour_nodes > 0, "{b:?}");
            assert!(b.contour_slide_max > 0.0, "{b:?}");
            assert!(
                b.wall_angle_max_after <= a.wall_angle_max_after + 0.05,
                "{b:?} {a:?}"
            );
            assert!(worst(&slid) >= 0.99 * worst(&fixed));
            if c == 0 {
                assert!(elbow_wall_error(&slid, r, big, 1.0) < 1e-12);
            }
            assert_conforming(&slid);
        }
    }

    /// Invariant: nodes of a planar sliding patch move inside their plane
    /// and stay on it; a patch that is not planar (the pipe wall) keeps its
    /// nodes fixed even when it is named as sliding.
    #[test]
    fn sliding_nodes_stay_on_their_plane() {
        let t = pipe_ogrid(0.05, 0.2, 0.55, [8, 10, 6], Some(0.0005)).unwrap();
        let mut g = build(&t.topology).unwrap();
        let tfi = g.points.clone();
        let (inlet, outlet, wall) = (
            patch(&t.topology, "inlet"),
            patch(&t.topology, "outlet"),
            patch(&t.topology, "walls"),
        );
        let options = SmoothOptions {
            sliding_patches: vec![inlet, outlet, wall],
            ..walls(&t.topology)
        };
        let res = smooth(&mut g, &options).unwrap();
        assert!(res.converged && res.sliding_nodes > 0, "{res:?}");
        let mut moved = 0;
        for &(q, p, _) in &g.boundary {
            for v in q {
                let (a, b) = (g.points[v], tfi[v]);
                if p == wall {
                    // rim nodes of the discs lie on the wall too
                    assert!((a.y.hypot(a.z) - 0.05).abs() < 1e-12, "{a:?}");
                    continue;
                }
                // the discs are the planes x = 0 and x = length
                assert!((a.x - b.x).abs() < 1e-12, "{a:?} left its plane");
                if a != b {
                    moved += 1;
                }
            }
        }
        assert!(moved > 0, "no disc node slid");
        assert_conforming(&g);
        for q in g.quality() {
            assert!(q.min_jacobian > 0.0, "{q:?}");
        }
    }
}

/// Largest distance of a wall node of a 90° elbow (bend in the x-y plane
/// about (0, R, ·), legs along x and y) from the pipe surface.
fn elbow_wall_error(g: &BlockGrid, r: f64, big: f64, inlet: f64) -> f64 {
    let walls = 2;
    let mut worst: f64 = 0.0;
    for (q, p, _) in &g.boundary {
        if *p != walls {
            continue;
        }
        for &v in q {
            let x = g.points[v];
            // distance to the centreline
            let d = if x.x < 0.0 {
                debug_assert!(x.x >= -inlet - 1e-12);
                (x.y * x.y + x.z * x.z).sqrt()
            } else if x.y > big {
                ((x.x - big).powi(2) + x.z * x.z).sqrt()
            } else {
                let rho = (x.x * x.x + (x.y - big).powi(2)).sqrt();
                ((rho - big).powi(2) + x.z * x.z).sqrt()
            };
            worst = worst.max((d - r).abs());
        }
    }
    worst
}

/// Invariant: the elbow conforms, has positive Jacobians, its wall nodes
/// lie on the cylinders and the torus (projection, BS-07), and its volume
/// approaches the Pappus value from below at second order as the
/// cross-section is refined.
#[test]
fn elbow_matches_its_torus() {
    let (r, big) = (0.5, 1.5);
    let mut last = f64::NAN;
    for core in [4, 8, 16] {
        let t = elbow(
            r,
            big,
            90.0,
            1.0,
            1.0,
            0.5,
            [core, core, 6 * core],
            Some(0.01),
        )
        .unwrap();
        let g = build(&t.topology).unwrap();
        assert_eq!(g.blocks.len(), 25);
        assert_conforming(&g);
        assert!(g.quality().iter().all(|q| q.min_jacobian > 0.0));
        let v = total(&g);
        assert!((v - g.boundary_volume()).abs() < 1e-10 * v);
        assert!(elbow_wall_error(&g, r, big, 1.0) < 1e-12);
        let d = (v - t.analytic_volume) / t.analytic_volume;
        assert!(d < 0.0, "chords lie inside the circle: {d}");
        if last.is_finite() {
            let ratio = last / d;
            assert!(ratio > 3.5 && ratio < 4.5, "order of convergence: {ratio}");
        }
        last = d;
    }
}

/// Invariant: a face whose corners are off its surface, or a surface on a
/// face between two blocks, is refused before any grid is built.
#[test]
fn misplaced_surfaces_are_refused() {
    let mut t = pipe_ogrid(1.0, 2.0, 0.5, [4, 4, 4], None).unwrap().topology;
    let mut bad = t.clone();
    if let mesh_core::block::Surface::Cylinder { radius, .. } = &mut bad.surfaces[0].surface {
        *radius = 1.1;
    }
    let issues = check(&bad).unwrap_err();
    assert!(issues.iter().all(|i| i.kind == "surface"));
    assert!(issues[0].message.contains("not on the surface"));
    // the core block's i = 0 face lies between the core and an outer block
    let core = t.blocks[0].vertices;
    t.surfaces[0].faces.push(face_vertices(&core, 0));
    let issues = check(&t).unwrap_err();
    assert!(issues[0].message.contains("no boundary block face"));
}

fn naca_case(digits: &str, side: usize) -> (BlockGrid, f64) {
    let t = airfoil_cgrid(
        &Naca4::parse(digits).unwrap().section(),
        &AirfoilCgrid {
            chord: 1.0,
            span: 0.1,
            farfield: 5.0,
            wake: 5.0,
            cells: [side, side / 2, side / 2, 1],
            wall_spacing: 1e-3,
            le_spacing: None,
            te_spacing: None,
        },
    )
    .unwrap();
    (build(&t.topology).unwrap(), t.analytic_volume)
}

/// Largest distance of a wall node from the section, measured against the
/// segments of a fine polygon of the analytic shape.
fn section_error(g: &BlockGrid, section: &Naca4) -> f64 {
    let n = 20_000;
    let mut ring = Vec::new();
    for side in [1.0, -1.0] {
        for i in 0..=n {
            let xc = 0.5 * (1.0 - (std::f64::consts::PI * i as f64 / n as f64).cos());
            ring.push(section.point(xc, side));
        }
    }
    let mut worst: f64 = 0.0;
    for (q, p, _) in &g.boundary {
        if *p != 0 {
            continue;
        }
        for &v in q {
            let x = g.points[v];
            let d = ring
                .windows(2)
                .map(|w| {
                    let (a, b) = (w[0], w[1]);
                    let (ex, ey) = (b[0] - a[0], b[1] - a[1]);
                    let l2 = ex * ex + ey * ey;
                    let t = if l2 > 0.0 {
                        (((x.x - a[0]) * ex + (x.y - a[1]) * ey) / l2).clamp(0.0, 1.0)
                    } else {
                        0.0
                    };
                    (a[0] + t * ex - x.x).hypot(a[1] + t * ey - x.y)
                })
                .fold(f64::INFINITY, f64::min);
            worst = worst.max(d);
        }
    }
    worst
}

/// Invariant: the airfoil C-grid conforms, has positive Jacobians, keeps
/// the wall nodes on the section (spline through 160 points per side), and
/// its volume converges at second order to the analytic one (C, wake and
/// section areas) as the cells double.
#[test]
fn airfoil_cgrid_wraps_its_section() {
    for digits in ["0012", "2412"] {
        let section = Naca4::parse(digits).unwrap();
        let mut last = f64::NAN;
        for side in [24, 48, 96] {
            let (g, exact) = naca_case(digits, side);
            assert_eq!(g.blocks.len(), 6);
            assert_conforming(&g);
            let jmin = g
                .quality()
                .iter()
                .map(|q| q.min_scaled_jacobian)
                .fold(1.0, f64::min);
            assert!(jmin > 0.4, "{digits} {side}: {jmin}");
            let v = total(&g);
            assert!((v - g.boundary_volume()).abs() < 1e-10 * v);
            let d = (v - exact) / exact;
            assert!(d < 0.0, "chords lie inside the C: {d}");
            if last.is_finite() {
                let ratio = last / d;
                assert!(ratio > 3.0 && ratio < 4.5, "{digits}: order {ratio}");
            }
            last = d;
            if side == 48 {
                let e = section_error(&g, &section);
                assert!(e < 1e-6, "{digits}: wall nodes {e} off the section");
            }
        }
    }
    assert!(Naca4::parse("12").is_err() && Naca4::parse("2012").is_err());
}

/// Invariant: a Selig and a Lednicer file of the same section, moved,
/// turned, scaled and with a blunt trailing edge, read to the same unit
/// section with a closed trailing edge; the C-grid builds from it. (The
/// section is symmetric so that its data leading edge is the point
/// farthest from the trailing edge, as the Selig reader defines it.)
#[test]
fn coordinate_files_read_to_one_section() {
    let naca = Naca4::parse("0012").unwrap().section();
    // to chord 2, turned by 5°, moved, with a gap of 0.01 chords
    let (c, s) = (5f64.to_radians().cos(), 5f64.to_radians().sin());
    let place = |p: [f64; 2], sign: f64| {
        let (x, y) = (2.0 * p[0], 2.0 * (p[1] + sign * 0.005 * p[0]));
        [3.0 + c * x - s * y, -1.0 + s * x + c * y]
    };
    let mut selig = String::from("NACA 0012 blunt\n");
    for p in naca.upper.iter().rev() {
        let q = place(*p, 1.0);
        selig += &format!("{:.12} {:.12}\n", q[0], q[1]);
    }
    for p in &naca.lower[1..] {
        let q = place(*p, -1.0);
        selig += &format!("{:.12} {:.12}\n", q[0], q[1]);
    }
    let mut lednicer = format!(
        "NACA 0012 blunt\n{}. {}.\n\n",
        naca.upper.len(),
        naca.lower.len()
    );
    for side in [(&naca.upper, 1.0), (&naca.lower, -1.0)] {
        for p in side.0 {
            let q = place(*p, side.1);
            lednicer += &format!("{:.12} {:.12}\n", q[0], q[1]);
        }
        lednicer += "\n";
    }
    for text in [&selig, &lednicer] {
        let sec = Section::from_coordinates(text).unwrap();
        assert!(
            (sec.trailing_edge_gap - 0.01).abs() < 1e-6,
            "{}",
            sec.trailing_edge_gap
        );
        for (a, b) in sec
            .upper
            .iter()
            .zip(&naca.upper)
            .chain(sec.lower.iter().zip(&naca.lower))
        {
            assert!((a[0] - b[0]).abs() < 1e-6 && (a[1] - b[1]).abs() < 1e-6);
        }
        assert!((sec.area - naca.area).abs() < 1e-4 * naca.area);
        let t = airfoil_cgrid(
            &sec,
            &AirfoilCgrid {
                chord: 1.0,
                span: 0.1,
                farfield: 5.0,
                wake: 5.0,
                cells: [32, 16, 16, 1],
                wall_spacing: 1e-3,
                le_spacing: None,
                te_spacing: None,
            },
        )
        .unwrap();
        assert_conforming(&build(&t.topology).unwrap());
    }
    assert!(Section::from_coordinates("x\n1 0\n0 0\n").is_err());
}

/// Invariant: on the NACA 2412 C-grid the wall nodes slide along the
/// section (they stay on it), the smoother still converges, and the grid
/// is no worse than with fixed wall nodes: smaller largest wall angle, no
/// smaller scaled Jacobian.
#[test]
fn wall_nodes_slide_along_the_section() {
    use super::smooth::{SmoothOptions, smooth};
    let section = Naca4::parse("2412").unwrap();
    let mut runs = Vec::new();
    for slide in [false, true] {
        let (mut g, _) = naca_case("2412", 48);
        let cg = airfoil_cgrid(
            &section.section(),
            &AirfoilCgrid {
                chord: 1.0,
                span: 0.1,
                farfield: 5.0,
                wake: 5.0,
                cells: [48, 24, 24, 1],
                wall_spacing: 1e-3,
                le_spacing: None,
                te_spacing: None,
            },
        )
        .unwrap();
        let named = |f: &dyn Fn(&str) -> bool| -> Vec<usize> {
            cg.topology
                .patches
                .iter()
                .enumerate()
                .filter(|(_, p)| f(&p.name))
                .map(|(i, _)| i)
                .collect()
        };
        let o = SmoothOptions {
            wall_patches: named(&|n| n.contains("wall")),
            sliding_patches: named(&|n| n.starts_with("symmetry")),
            slide_walls: slide,
            ..SmoothOptions::default()
        };
        let res = smooth(&mut g, &o).unwrap();
        assert!(res.converged, "{res:?}");
        let jmin = g
            .quality()
            .iter()
            .map(|q| q.min_scaled_jacobian)
            .fold(1.0, f64::min);
        runs.push((res, jmin, section_error(&g, &section)));
    }
    let (fixed, slid) = (&runs[0], &runs[1]);
    assert_eq!(fixed.0.contour_nodes, 0);
    // two contours (upper, lower) on each of the two symmetry planes
    assert!(slid.0.contour_nodes > 100, "{:?}", slid.0);
    assert!(slid.0.contour_slide_max > 0.0);
    assert!(
        slid.0.wall_angle_max_after < fixed.0.wall_angle_max_after,
        "{:?} {:?}",
        slid.0,
        fixed.0
    );
    assert!(slid.1 >= fixed.1, "{} < {}", slid.1, fixed.1);
    assert!(slid.2 < 2e-5, "wall nodes {} off the section", slid.2);
}

fn multi_case(
    slat: Option<(f64, f64, f64)>,
    flap: Option<(f64, f64, f64)>,
) -> Result<templates::Template, String> {
    use super::multi_element::{MultiElement, Placement, multi_element_airfoil};
    let place = |digits: &str, chord: f64, (deflection, gap, overlap): (f64, f64, f64)| Placement {
        section: Naca4::parse(digits).unwrap().section(),
        chord,
        deflection,
        gap,
        overlap,
    };
    let slat = slat.map(|s| place("4410", 0.16, s));
    let flap = flap.map(|f| place("4412", 0.3, f));
    multi_element_airfoil(
        &Naca4::parse("4412").unwrap().section(),
        slat.as_ref(),
        flap.as_ref(),
        &MultiElement {
            chord: 1.0,
            span: 0.1,
            farfield: 5.0,
            wake: 5.0,
            surface_cells: 32,
            normal_cells: 16,
            wake_cells: 16,
            span_cells: 1,
            thickness_cells: 8,
            gap_cells: 8,
            upstream_cells: 16,
            wall_spacing: 1e-3,
            trailing_edge_thickness: 0.005,
        },
    )
}

/// Invariant: the multi-element template (slat, main, flap) conforms, has
/// positive Jacobians, gets the box area minus the three sections (the
/// spline walls are second-order close), places the flap at the asked
/// gap from the main trailing edge, and works with the slat or the flap
/// alone; its layout optimisation leaves no cell over 70° of
/// non-orthogonality and reports what the written topology grids to; a
/// slat whose trailing edge sits well above the main nose (the case that
/// folded before the wall rings and the layout optimisation) meshes; it
/// refuses overlapping elements and names the gap whose blocks still fold.
/// A change that dropped the rings or the optimisation (max 86-90°, up to
/// 433 cells over 70°) breaks it; no other test grids the template.
#[test]
fn multi_element_airfoil_places_and_conforms() {
    for (slat, flap, blocks) in [
        (Some((25.0, 0.05, -0.04)), Some((30.0, 0.03, -0.02)), 40),
        (None, Some((20.0, 0.025, -0.015)), 28),
        (Some((15.0, 0.03, -0.025)), None, 28),
        (Some((20.0, 0.04, -0.03)), None, 28),
    ] {
        let t = multi_case(slat, flap).unwrap();
        let g = grid::build(&t.topology).unwrap();
        assert_eq!(g.blocks.len(), blocks);
        let l = t.layout.as_ref().unwrap();
        assert_eq!(l.after.cells_over_70, 0, "{slat:?} {flap:?} {l:?}");
        assert_eq!(l.after.folded_corners, 0, "{l:?}");
        assert!(
            l.after.max_deg < 70.0 && l.after.max_deg < l.before.max_deg,
            "{l:?}"
        );
        let again = super::layout::orthogonality(&t.topology).unwrap();
        assert!(
            (again.max_deg - l.after.max_deg).abs() < 1e-9,
            "{again:?} {l:?}"
        );
        assert_conforming(&g);
        let jmin = g
            .quality()
            .iter()
            .map(|q| q.min_jacobian)
            .fold(f64::INFINITY, f64::min);
        assert!(jmin > 0.0, "{jmin}");
        let v = total(&g);
        assert!((v - g.boundary_volume()).abs() < 1e-9 * v);
        assert!(((v - t.analytic_volume) / t.analytic_volume).abs() < 1e-4);
        let names: Vec<&str> = t.topology.patches.iter().map(|p| p.name.as_str()).collect();
        assert_eq!(names.contains(&"slat_wall"), slat.is_some());
        assert_eq!(names.contains(&"flap_wall"), flap.is_some());
        assert!(names.contains(&"main_wall") && names.contains(&"inlet"));
    }
    // the flap gap: the main element's trailing edge (middle vertex of its
    // blunt edge) to the nearest flap wall vertex on the section
    let t = multi_case(None, Some((20.0, 0.025, -0.015))).unwrap();
    let g = grid::build(&t.topology).unwrap();
    let wall = |name: &str| {
        let p = t
            .topology
            .patches
            .iter()
            .position(|p| p.name == name)
            .unwrap();
        g.boundary
            .iter()
            .filter(move |(_, patch, _)| *patch == p)
            .flat_map(|(q, _, _)| *q)
    };
    let te = Vec3::new(1.0, 0.0, 0.0);
    let gap = wall("flap_wall")
        .map(|i| {
            let p = g.points[i];
            (Vec3::new(p.x, p.y, 0.0) - te).length()
        })
        .fold(f64::INFINITY, f64::min);
    assert!((gap - 0.025).abs() < 2e-3, "flap gap {gap}");
    let e = multi_case(None, Some((20.0, 0.025, 0.01))).err().unwrap();
    assert!(e.contains("negative overlap"), "{e}");
    let e = multi_case(Some((10.0, 0.02, -0.005)), None).err().unwrap();
    assert!(e.contains("between the slat and the main element"), "{e}");
}
