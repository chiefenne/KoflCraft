//! A compressed-row sparse matrix and a Jacobi-preconditioned conjugate
//! gradient solver with Dirichlet values (Hestenes & Stiefel 1952; the
//! textbook form of Shewchuk 1994, *An Introduction to the Conjugate
//! Gradient Method Without the Agonizing Pain*).

/// A symmetric matrix in compressed rows.
#[derive(Clone, Debug, Default)]
pub struct Csr {
    pub n: usize,
    pub start: Vec<usize>,
    pub col: Vec<usize>,
    pub val: Vec<f64>,
}

/// What a solve did.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SolveStats {
    pub iterations: usize,
    /// Final residual norm over the norm of the right-hand side.
    pub residual: f64,
    pub converged: bool,
}

impl Csr {
    /// Sum duplicate (row, col) entries.
    pub fn from_triplets(n: usize, mut t: Vec<(usize, usize, f64)>) -> Self {
        t.sort_unstable_by_key(|&(i, j, _)| (i, j));
        let mut start = vec![0; n + 1];
        let mut col = Vec::with_capacity(t.len());
        let mut val: Vec<f64> = Vec::with_capacity(t.len());
        let mut last: Option<(usize, usize)> = None;
        for (i, j, v) in t {
            if last == Some((i, j)) {
                if let Some(x) = val.last_mut() {
                    *x += v;
                }
                continue;
            }
            last = Some((i, j));
            start[i + 1] += 1;
            col.push(j);
            val.push(v);
        }
        for i in 0..n {
            start[i + 1] += start[i];
        }
        Self { n, start, col, val }
    }

    pub fn diagonal(&self) -> Vec<f64> {
        (0..self.n)
            .map(|i| {
                (self.start[i]..self.start[i + 1])
                    .find(|&k| self.col[k] == i)
                    .map_or(0.0, |k| self.val[k])
            })
            .collect()
    }

    /// `s * self + diag(d)`.
    pub fn scaled_plus_diagonal(&self, s: f64, d: &[f64]) -> Self {
        let mut out = self.clone();
        for v in &mut out.val {
            *v *= s;
        }
        for i in 0..self.n {
            if let Some(k) = (self.start[i]..self.start[i + 1]).find(|&k| self.col[k] == i) {
                out.val[k] += d[i];
            }
        }
        out
    }

    pub fn mul(&self, x: &[f64], y: &mut [f64]) {
        for i in 0..self.n {
            let mut s = 0.0;
            for k in self.start[i]..self.start[i + 1] {
                s += self.val[k] * x[self.col[k]];
            }
            y[i] = s;
        }
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Solve `a x = b` for the rows without a fixed value; fixed rows take
/// their value, and their columns move to the right-hand side, so the
/// reduced system stays symmetric positive definite. `x0` is the start
/// (free entries only).
pub fn solve(
    a: &Csr,
    b: &[f64],
    fixed: &[Option<f64>],
    x0: Option<&[f64]>,
    tolerance: f64,
    max_iterations: usize,
) -> (Vec<f64>, SolveStats) {
    let n = a.n;
    let free: Vec<bool> = fixed.iter().map(Option::is_none).collect();
    let mut x: Vec<f64> = (0..n)
        .map(|i| match fixed[i] {
            Some(v) => v,
            None => x0.map_or(0.0, |s| s[i]),
        })
        .collect();
    // right-hand side of the free rows: b - A_fixed x_fixed
    let mut rhs = vec![0.0; n];
    for i in 0..n {
        if !free[i] {
            continue;
        }
        let mut s = b[i];
        for k in a.start[i]..a.start[i + 1] {
            let j = a.col[k];
            if !free[j] {
                s -= a.val[k] * x[j];
            }
        }
        rhs[i] = s;
    }
    let apply = |v: &[f64], out: &mut [f64]| {
        for i in 0..n {
            if !free[i] {
                out[i] = 0.0;
                continue;
            }
            let mut s = 0.0;
            for k in a.start[i]..a.start[i + 1] {
                let j = a.col[k];
                if free[j] {
                    s += a.val[k] * v[j];
                }
            }
            out[i] = s;
        }
    };
    let diag = a.diagonal();
    let inv: Vec<f64> = (0..n)
        .map(|i| {
            if free[i] && diag[i] > 0.0 {
                1.0 / diag[i]
            } else {
                0.0
            }
        })
        .collect();
    let norm_b = dot(&rhs, &rhs).sqrt();
    if norm_b == 0.0 {
        for i in 0..n {
            if free[i] {
                x[i] = 0.0;
            }
        }
        return (
            x,
            SolveStats {
                iterations: 0,
                residual: 0.0,
                converged: true,
            },
        );
    }
    let mut xf: Vec<f64> = (0..n).map(|i| if free[i] { x[i] } else { 0.0 }).collect();
    let mut ax = vec![0.0; n];
    apply(&xf, &mut ax);
    let mut r: Vec<f64> = (0..n).map(|i| rhs[i] - ax[i]).collect();
    let mut z: Vec<f64> = (0..n).map(|i| inv[i] * r[i]).collect();
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    let mut q = vec![0.0; n];
    let mut stats = SolveStats::default();
    for it in 0..max_iterations {
        let res = dot(&r, &r).sqrt() / norm_b;
        stats.iterations = it;
        stats.residual = res;
        if res <= tolerance {
            stats.converged = true;
            break;
        }
        apply(&p, &mut q);
        let pq = dot(&p, &q);
        if pq <= 0.0 || !pq.is_finite() {
            break;
        }
        let alpha = rz / pq;
        for i in 0..n {
            xf[i] += alpha * p[i];
            r[i] -= alpha * q[i];
        }
        for i in 0..n {
            z[i] = inv[i] * r[i];
        }
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
    }
    if !stats.converged {
        let res = dot(&r, &r).sqrt() / norm_b;
        stats.residual = res;
        stats.converged = res <= tolerance;
        stats.iterations = max_iterations;
    }
    for i in 0..n {
        if free[i] {
            x[i] = xf[i];
        }
    }
    (x, stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Keep: the 1D Laplacian with both ends fixed must give the straight
    /// line; a sign or elimination error in `solve` bends it, and the
    /// field tests above it would only show a vaguer symptom.
    #[test]
    fn line_between_fixed_ends() {
        let n = 11;
        let mut t = Vec::new();
        for i in 0..n - 1 {
            t.extend([
                (i, i, 1.0),
                (i + 1, i + 1, 1.0),
                (i, i + 1, -1.0),
                (i + 1, i, -1.0),
            ]);
        }
        let a = Csr::from_triplets(n, t);
        let mut fixed = vec![None; n];
        fixed[0] = Some(0.0);
        fixed[n - 1] = Some(1.0);
        let (x, s) = solve(&a, &vec![0.0; n], &fixed, None, 1e-12, 100);
        assert!(s.converged);
        for (i, v) in x.iter().enumerate() {
            assert!((v - i as f64 / 10.0).abs() < 1e-9);
        }
    }
}
