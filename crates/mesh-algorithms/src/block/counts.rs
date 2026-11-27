//! The count solve of a block topology (BS-05, BS-21): one integer per
//! sheet, close to target counts, with the sum rule at every T-junction.
//!
//! A *sheet* is a class of block axes that must carry one count: the four
//! parallel edges of a block, and every block that shares one of them
//! whole (the chains of BS-05). Without T-junctions each sheet is free and
//! its best count is its rounded target. A T-junction ties sheets by a sum,
//! count(whole edge) = Σ count(pieces), so the counts of the tied sheets
//! lie on an integer lattice: the kernel of the integer matrix of the sums.
//! The solve takes an integer basis of that kernel (column Hermite
//! reduction), the least-squares point of the relative deviations on it,
//! rounds it and improves it by lattice moves until no move helps, with
//! every count at least 1.
//!
//! This is the quantisation step of quad and hex layouts as integer
//! programs over the layout's arcs or sheets (Lyon, Campen and Kobbelt
//! 2021; Brückler, Bommes and Campen 2022; the sheet form follows the
//! integer-sheet-pump of Brückler et al., SGP 2024), here on block axes
//! with a small exact lattice search instead of a general integer solver.
//! No code was read.

use std::collections::{BTreeMap, HashMap};

use mesh_core::block::BlockTopology;

use super::{Junctions, edge_key, edge_slots};

/// What the count solve did.
#[derive(Clone, Debug, PartialEq)]
pub struct CountFit {
    /// The new counts of every block.
    pub cells: Vec<[usize; 3]>,
    /// Sheets (classes of block axes with one count) and how many of them
    /// the T-junction sums tie together.
    pub sheets: usize,
    pub tied_sheets: usize,
    /// Sum rules (one per T-junction edge).
    pub sums: usize,
    /// Largest |count − target| / target over the sheets.
    pub max_deviation: f64,
}

struct Classes(Vec<usize>);

impl Classes {
    fn find(&mut self, mut x: usize) -> usize {
        while self.0[x] != x {
            self.0[x] = self.0[self.0[x]];
            x = self.0[x];
        }
        x
    }
    fn union(&mut self, a: usize, b: usize) {
        let (a, b) = (self.find(a), self.find(b));
        self.0[a] = b;
    }
}

/// Fit the counts of `t` to `targets` (cells per block and axis, > 0):
/// every chain of opposite edges one count, every T-junction edge the sum
/// of its pieces, every count at least 1, the relative deviations from
/// the targets as small as the integers allow (least squares).
pub fn fit(
    t: &BlockTopology,
    junctions: &Junctions,
    targets: &[[f64; 3]],
) -> Result<CountFit, String> {
    let nb = t.blocks.len();
    if targets.len() != nb {
        return Err("one target per block is needed".into());
    }
    let mut cl = Classes((0..3 * nb).collect());
    let mut on_edge: HashMap<(usize, usize), usize> = HashMap::new();
    for (b, block) in t.blocks.iter().enumerate() {
        for (axis, s, e) in edge_slots() {
            let key = edge_key(block.vertices[s], block.vertices[e]);
            let id = 3 * b + axis;
            match on_edge.get(&key) {
                Some(&o) => cl.union(id, o),
                None => {
                    on_edge.insert(key, id);
                }
            }
        }
    }
    // sheets and their targets (mean of the axes' targets)
    let mut sheet_of: BTreeMap<usize, usize> = BTreeMap::new();
    let mut members: Vec<Vec<usize>> = Vec::new();
    for id in 0..3 * nb {
        let r = cl.find(id);
        let s = *sheet_of.entry(r).or_insert_with(|| {
            members.push(Vec::new());
            members.len() - 1
        });
        members[s].push(id);
    }
    let ns = members.len();
    let sheet = |cl: &mut Classes, id: usize| sheet_of[&cl.find(id)];
    let mut target = vec![0.0; ns];
    let mut weight = vec![0.0; ns];
    for (s, m) in members.iter().enumerate() {
        for &id in m {
            let v = targets[id / 3][id % 3];
            if !(v.is_finite() && v > 0.0) {
                return Err(format!("block {} has no positive target count", id / 3));
            }
            target[s] += v;
        }
        weight[s] = m.len() as f64;
        target[s] /= weight[s];
    }
    // sum rules over sheets
    let mut rows: Vec<BTreeMap<usize, i128>> = Vec::new();
    let mut keys: Vec<_> = junctions.through.keys().copied().collect();
    keys.sort_unstable();
    for key in keys {
        let chain = &junctions.through[&key];
        let Some(&w) = on_edge.get(&key) else {
            return Err(format!(
                "the junction edge {}–{} is no block edge",
                key.0, key.1
            ));
        };
        let mut row: BTreeMap<usize, i128> = BTreeMap::new();
        *row.entry(sheet(&mut cl, w)).or_default() += 1;
        for p in chain.windows(2) {
            let Some(&id) = on_edge.get(&edge_key(p[0], p[1])) else {
                return Err(format!("the piece {}–{} is no block edge", p[0], p[1]));
            };
            *row.entry(sheet(&mut cl, id)).or_default() -= 1;
        }
        row.retain(|_, v| *v != 0);
        if row.is_empty() {
            return Err(format!(
                "edge {}–{} and its pieces are the same sheet; no count can satisfy the T-junction",
                key.0, key.1
            ));
        }
        rows.push(row);
    }
    let mut x: Vec<i128> = target.iter().map(|&v| (v.round() as i128).max(1)).collect();
    // the tied sheets: columns of the sum rules
    let mut tied: Vec<usize> = rows.iter().flat_map(|r| r.keys().copied()).collect();
    tied.sort_unstable();
    tied.dedup();
    if !tied.is_empty() {
        let col: HashMap<usize, usize> = tied.iter().enumerate().map(|(i, &s)| (s, i)).collect();
        let a: Vec<Vec<i128>> = rows
            .iter()
            .map(|r| {
                let mut v = vec![0; tied.len()];
                for (&s, &c) in r {
                    v[col[&s]] = c;
                }
                v
            })
            .collect();
        let basis = kernel(&a, tied.len())?;
        if basis.is_empty() {
            return Err("the T-junction sums leave no free count".into());
        }
        let tt: Vec<f64> = tied.iter().map(|&s| target[s]).collect();
        let ww: Vec<f64> = tied
            .iter()
            .map(|&s| weight[s] / (target[s] * target[s]))
            .collect();
        let y = lattice_fit(&basis, &tt, &ww)?;
        for (i, &s) in tied.iter().enumerate() {
            x[s] = y[i];
        }
    }
    if let Some(s) = (0..ns).find(|&s| x[s] < 1) {
        return Err(format!(
            "no counts of at least 1 satisfy the T-junctions (sheet of block {} axis {})",
            members[s][0] / 3,
            members[s][0] % 3
        ));
    }
    let mut cells = vec![[0usize; 3]; nb];
    for (s, m) in members.iter().enumerate() {
        for &id in m {
            cells[id / 3][id % 3] = x[s] as usize;
        }
    }
    let max_deviation = (0..ns)
        .map(|s| (x[s] as f64 - target[s]).abs() / target[s])
        .fold(0.0, f64::max);
    Ok(CountFit {
        cells,
        sheets: ns,
        tied_sheets: tied.len(),
        sums: rows.len(),
        max_deviation,
    })
}

/// An integer basis of { x ∈ ℤⁿ : A x = 0 } by column Hermite reduction:
/// unimodular column operations bring A to echelon form; the columns of
/// the transform past the pivots span the kernel.
fn kernel(a: &[Vec<i128>], n: usize) -> Result<Vec<Vec<i128>>, String> {
    let mut m: Vec<Vec<i128>> = a.to_vec();
    let mut u: Vec<Vec<i128>> = (0..n)
        .map(|i| (0..n).map(|j| i128::from(i == j)).collect())
        .collect();
    let col_op =
        |m: &mut Vec<Vec<i128>>, u: &mut Vec<Vec<i128>>, dst: usize, src: usize, q: i128| {
            for row in m.iter_mut() {
                row[dst] -= q * row[src];
            }
            for row in u.iter_mut() {
                row[dst] -= q * row[src];
            }
        };
    let swap = |m: &mut Vec<Vec<i128>>, u: &mut Vec<Vec<i128>>, i: usize, j: usize| {
        for row in m.iter_mut() {
            row.swap(i, j);
        }
        for row in u.iter_mut() {
            row.swap(i, j);
        }
    };
    let mut p = 0;
    for r in 0..m.len() {
        if p == n {
            break;
        }
        while let Some(best) = (p..n)
            .filter(|&c| m[r][c] != 0)
            .min_by_key(|&c| m[r][c].abs())
        {
            swap(&mut m, &mut u, p, best);
            let mut done = true;
            for c in p + 1..n {
                if m[r][c] != 0 {
                    let q = m[r][c].div_euclid(m[r][p]);
                    col_op(&mut m, &mut u, c, p, q);
                    if m[r][c] != 0 {
                        done = false;
                    }
                }
            }
            if u.iter().flatten().any(|v| v.abs() > 1 << 40) {
                return Err("the count solve's lattice grew too large".into());
            }
            if done {
                p += 1;
                break;
            }
        }
    }
    let mut basis: Vec<Vec<i128>> = (p..n)
        .map(|c| u.iter().map(|row| row[c]).collect())
        .collect();
    // pairwise size reduction keeps the moves short
    for _ in 0..8 {
        let mut changed = false;
        for i in 0..basis.len() {
            for j in 0..basis.len() {
                if i == j {
                    continue;
                }
                let dot: i128 = basis[i].iter().zip(&basis[j]).map(|(a, b)| a * b).sum();
                let nn: i128 = basis[j].iter().map(|b| b * b).sum();
                let q = (dot as f64 / nn as f64).round() as i128;
                if q != 0 {
                    let bj = basis[j].clone();
                    for (x, y) in basis[i].iter_mut().zip(&bj) {
                        *x -= q * y;
                    }
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    Ok(basis)
}

/// The integer point y = K z (K: the basis vectors as columns) nearest to
/// the targets in the weighted norm, every entry at least 1: least squares
/// in z, rounded, then single and paired basis moves while they help.
fn lattice_fit(basis: &[Vec<i128>], t: &[f64], w: &[f64]) -> Result<Vec<i128>, String> {
    let (n, d) = (t.len(), basis.len());
    // normal equations (Kᵀ W K) z = Kᵀ W t
    let mut m = vec![vec![0.0; d + 1]; d];
    for i in 0..d {
        for j in 0..d {
            m[i][j] = (0..n)
                .map(|c| basis[i][c] as f64 * w[c] * basis[j][c] as f64)
                .sum();
        }
        m[i][d] = (0..n).map(|c| basis[i][c] as f64 * w[c] * t[c]).sum();
    }
    for col in 0..d {
        let piv = (col..d)
            .max_by(|&a, &b| m[a][col].abs().total_cmp(&m[b][col].abs()))
            .unwrap_or(col);
        m.swap(col, piv);
        if m[col][col].abs() < 1e-300 {
            return Err("the count solve is singular".into());
        }
        let pivot = m[col].clone();
        for (r, row) in m.iter_mut().enumerate() {
            if r != col {
                let f = row[col] / pivot[col];
                for (x, y) in row.iter_mut().zip(&pivot).skip(col) {
                    *x -= f * y;
                }
            }
        }
    }
    let z: Vec<i128> = (0..d)
        .map(|i| (m[i][d] / m[i][i]).round() as i128)
        .collect();
    let point = |z: &[i128]| -> Vec<i128> {
        (0..n)
            .map(|c| (0..d).map(|i| basis[i][c] * z[i]).sum())
            .collect()
    };
    let cost = |y: &[i128]| -> f64 {
        (0..n)
            .map(|c| {
                let v = y[c] as f64;
                let dev = w[c] * (v - t[c]) * (v - t[c]);
                let low = (1.0 - v).max(0.0);
                dev + 1e9 * low * low
            })
            .sum()
    };
    let mut z = z;
    let mut best = cost(&point(&z));
    let mut moves: Vec<Vec<i128>> = Vec::new();
    for i in 0..d {
        for s in [-1, 1] {
            let mut mv = vec![0; d];
            mv[i] = s;
            moves.push(mv);
        }
    }
    if d <= 40 {
        for i in 0..d {
            for j in i + 1..d {
                for (si, sj) in [(1, 1), (1, -1), (-1, 1), (-1, -1)] {
                    let mut mv = vec![0; d];
                    mv[i] = si;
                    mv[j] = sj;
                    moves.push(mv);
                }
            }
        }
    }
    for _ in 0..10_000 {
        let mut improved = false;
        for mv in &moves {
            let cand: Vec<i128> = z.iter().zip(mv).map(|(a, b)| a + b).collect();
            let c = cost(&point(&cand));
            if c < best - 1e-12 {
                best = c;
                z = cand;
                improved = true;
            }
        }
        if !improved {
            break;
        }
    }
    Ok(point(&z))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Invariant: the kernel basis spans exactly the solutions of the
    /// sums; a wrong column operation leaves a vector with A x ≠ 0.
    #[test]
    fn kernel_vectors_solve_the_sums() {
        // x0 = x1 + x2, x3 = x1 + x4, x0 = x3
        let a = vec![
            vec![1, -1, -1, 0, 0],
            vec![0, -1, 0, 1, -1],
            vec![1, 0, 0, -1, 0],
        ];
        let k = kernel(&a, 5).unwrap();
        assert_eq!(k.len(), 2);
        for v in &k {
            for row in &a {
                assert_eq!(row.iter().zip(v).map(|(a, b)| a * b).sum::<i128>(), 0);
            }
        }
    }

    /// Invariant: the fit meets the sums exactly and lands near the
    /// targets; a target that the sum cannot meet is shared out.
    #[test]
    fn lattice_fit_shares_out_the_mismatch() {
        // x0 = x1 + x2 with targets 4, 1, 4: the best integers share the
        // gap between 4 and 5
        let a = vec![vec![1, -1, -1]];
        let k = kernel(&a, 3).unwrap();
        let t = [4.0, 1.0, 4.0];
        let w: Vec<f64> = t.iter().map(|v| 1.0 / (v * v)).collect();
        let y = lattice_fit(&k, &t, &w).unwrap();
        assert_eq!(y[0], y[1] + y[2]);
        assert!(y.iter().all(|&v| v >= 1));
        assert_eq!(y[1], 1);
        assert!((4..=5).contains(&y[0]), "{y:?}");
    }
}
