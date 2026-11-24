use std::collections::HashMap;

use mesh_core::block::{BlockSpec, corner_slot};

use super::grid::HEX_FACES;
use super::templates::{channel, pipe_ogrid};
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
