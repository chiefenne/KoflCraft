//! Finite-volume quality of a block grid in 3D (BS-13, BS-20): the
//! non-orthogonality and growth across every inner face of the hexahedra
//! and the scaled Jacobian of every corner, for the whole grid and per
//! block. The 2D layout model of MB1g measures extruded topologies only;
//! this works on any grid, so the 3D trigger can judge a split.

use std::collections::HashMap;

use mesh_core::geometry::Vec3;

use super::corner_edges;
use super::grid::{BlockGrid, HEX_FACES};

/// The quality of a whole grid.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct GridQuality {
    /// Largest and mean angle (degrees) between an inner face's normal and
    /// the line through its two cell centroids.
    pub max_deg: f64,
    pub mean_deg: f64,
    /// Cells with an inner face above 70°.
    pub cells_over_70: usize,
    /// Largest ratio of the two centroid distances from an inner face.
    pub max_growth: f64,
    /// Hexahedron corners with a determinant ≤ 0.
    pub folded_corners: usize,
    pub min_scaled_jacobian: f64,
}

/// The quality of the cells of one block.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BlockStats {
    pub max_deg: f64,
    pub max_growth: f64,
    pub min_scaled_jacobian: f64,
}

/// Measure a grid: the whole and every block, in the order of
/// `grid.blocks`.
pub fn measure(grid: &BlockGrid) -> (GridQuality, Vec<BlockStats>) {
    let x = &grid.points;
    let mut block_of = vec![0usize; grid.hexes.len()];
    for (b, r) in grid.blocks.iter().enumerate() {
        let count = r.cells.iter().product::<usize>();
        block_of[r.first_cell..r.first_cell + count].fill(b);
    }
    let mut per: Vec<BlockStats> = grid
        .blocks
        .iter()
        .map(|_| BlockStats {
            max_deg: 0.0,
            max_growth: 1.0,
            min_scaled_jacobian: f64::INFINITY,
        })
        .collect();
    let mut q = GridQuality {
        max_growth: 1.0,
        min_scaled_jacobian: f64::INFINITY,
        ..GridQuality::default()
    };
    let centroid = |h: &[usize; 8]| h.iter().fold(Vec3::ZERO, |a, &v| a + x[v]) / 8.0;
    let mut faces: HashMap<[usize; 4], (usize, Vec3, Vec3)> = HashMap::new();
    let mut over = vec![false; grid.hexes.len()];
    let (mut sum, mut inner) = (0.0, 0usize);
    for (c, h) in grid.hexes.iter().enumerate() {
        let p = h.map(|v| x[v]);
        for slot in 0..8 {
            let e = corner_edges(&p, slot);
            let det = e[0].cross(e[1]).dot(e[2]);
            let l = e.map(Vec3::length);
            if det.is_nan() || det <= 0.0 {
                q.folded_corners += 1;
            }
            let sj = det / (l[0] * l[1] * l[2]);
            q.min_scaled_jacobian = q.min_scaled_jacobian.min(sj);
            let s = &mut per[block_of[c]];
            s.min_scaled_jacobian = s.min_scaled_jacobian.min(sj);
        }
        let cc = centroid(h);
        for f in HEX_FACES {
            let ring = f.map(|s| h[s]);
            let mut key = ring;
            key.sort_unstable();
            let fc = ring.iter().fold(Vec3::ZERO, |a, &v| a + x[v]) / 4.0;
            let n = (x[ring[2]] - x[ring[0]]).cross(x[ring[3]] - x[ring[1]]);
            match faces.remove(&key) {
                None => {
                    faces.insert(key, (c, fc, n));
                }
                Some((o, fo, no)) => {
                    let co = centroid(&grid.hexes[o]);
                    let d = cc - co;
                    let deg = (no.dot(d) / (no.length() * d.length()))
                        .clamp(-1.0, 1.0)
                        .acos()
                        .to_degrees();
                    let (da, db) = ((fo - co).length(), (cc - fc).length());
                    let growth = (da / db).max(db / da);
                    sum += deg;
                    inner += 1;
                    q.max_deg = q.max_deg.max(deg);
                    q.max_growth = q.max_growth.max(growth);
                    if deg > 70.0 {
                        over[c] = true;
                        over[o] = true;
                    }
                    for cell in [c, o] {
                        let s = &mut per[block_of[cell]];
                        s.max_deg = s.max_deg.max(deg);
                        s.max_growth = s.max_growth.max(growth);
                    }
                }
            }
        }
    }
    q.mean_deg = if inner > 0 { sum / inner as f64 } else { 0.0 };
    q.cells_over_70 = over.iter().filter(|&&o| o).count();
    (q, per)
}
