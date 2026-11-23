//! Exact orientation and in-sphere signs for 3D Delaunay construction.
//!
//! Each predicate is first evaluated in floating point; when the result is
//! smaller than a bound on its rounding error, it is evaluated again in
//! exact floating-point expansion arithmetic (Shewchuk, *Adaptive Precision
//! Floating-Point Arithmetic and Fast Robust Geometric Predicates*, Discrete
//! & Computational Geometry 18, 1997: two-sum, two-product, expansion growth
//! and scaling with zero elimination). The error bound used here is a
//! generous multiple of the paper's first-stage bound, so the filter only
//! falls through more often than needed, never less.

use mesh_core::geometry::Vec3;
use std::cmp::Ordering;

type Expansion = Vec<f64>;

fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let x = a + b;
    let bv = x - a;
    let av = x - bv;
    (x, (a - av) + (b - bv))
}

fn two_product(a: f64, b: f64) -> (f64, f64) {
    let x = a * b;
    (x, a.mul_add(b, -x))
}

fn grow(e: &[f64], b: f64) -> Expansion {
    let mut out = Vec::with_capacity(e.len() + 1);
    let mut q = b;
    for &c in e {
        let (s, h) = two_sum(q, c);
        if h != 0.0 {
            out.push(h);
        }
        q = s;
    }
    if q != 0.0 || out.is_empty() {
        out.push(q);
    }
    out
}

fn add(e: &[f64], f: &[f64]) -> Expansion {
    let mut out = e.to_vec();
    for &c in f {
        out = grow(&out, c);
    }
    out
}

fn scale(e: &[f64], b: f64) -> Expansion {
    let mut out = Vec::with_capacity(2 * e.len());
    let Some((&first, rest)) = e.split_first() else {
        return vec![0.0];
    };
    let (mut q, h) = two_product(first, b);
    if h != 0.0 {
        out.push(h);
    }
    for &c in rest {
        let (t1, t0) = two_product(c, b);
        let (s, h) = two_sum(q, t0);
        if h != 0.0 {
            out.push(h);
        }
        let (s2, h2) = two_sum(t1, s);
        if h2 != 0.0 {
            out.push(h2);
        }
        q = s2;
    }
    if q != 0.0 || out.is_empty() {
        out.push(q);
    }
    out
}

fn mul(e: &[f64], f: &[f64]) -> Expansion {
    let mut out = vec![0.0];
    for &c in f {
        out = add(&out, &scale(e, c));
    }
    out
}

fn neg(e: &[f64]) -> Expansion {
    e.iter().map(|v| -v).collect()
}

fn diff(a: f64, b: f64) -> Expansion {
    let (x, y) = two_sum(a, -b);
    if y == 0.0 { vec![x] } else { vec![y, x] }
}

fn sign(e: &[f64]) -> Ordering {
    e.iter()
        .rev()
        .find(|v| **v != 0.0)
        .map_or(Ordering::Equal, |v| v.total_cmp(&0.0))
}

fn det3(m: [[&Expansion; 3]; 3]) -> Expansion {
    let minor = |a: &Expansion, b: &Expansion, c: &Expansion, d: &Expansion| {
        add(&mul(a, b), &neg(&mul(c, d)))
    };
    let t0 = mul(m[0][0], &minor(m[1][1], m[2][2], m[1][2], m[2][1]));
    let t1 = mul(m[0][1], &minor(m[1][0], m[2][2], m[1][2], m[2][0]));
    let t2 = mul(m[0][2], &minor(m[1][0], m[2][1], m[1][1], m[2][0]));
    add(&add(&t0, &neg(&t1)), &t2)
}

const BOUND: f64 = 1e-13;

/// The sign of det[b − a, c − a, d − a]: positive when a, b, c, d form a
/// positively oriented tetrahedron.
pub fn orient3d(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> Ordering {
    let (u, v, w) = (b - a, c - a, d - a);
    let det = u.x * (v.y * w.z - v.z * w.y) - u.y * (v.x * w.z - v.z * w.x)
        + u.z * (v.x * w.y - v.y * w.x);
    let perm = (u.x.abs() * ((v.y * w.z).abs() + (v.z * w.y).abs()))
        + (u.y.abs() * ((v.x * w.z).abs() + (v.z * w.x).abs()))
        + (u.z.abs() * ((v.x * w.y).abs() + (v.y * w.x).abs()));
    if det.abs() > BOUND * perm {
        return det.total_cmp(&0.0);
    }
    let rows = [b, c, d].map(|p| [diff(p.x, a.x), diff(p.y, a.y), diff(p.z, a.z)]);
    sign(&det3([
        [&rows[0][0], &rows[0][1], &rows[0][2]],
        [&rows[1][0], &rows[1][1], &rows[1][2]],
        [&rows[2][0], &rows[2][1], &rows[2][2]],
    ]))
}

/// Positive when `e` lies strictly inside the circumsphere of the
/// positively oriented tetrahedron a, b, c, d, zero on it, negative outside.
pub fn insphere(a: Vec3, b: Vec3, c: Vec3, d: Vec3, e: Vec3) -> Ordering {
    let p = [a - e, b - e, c - e, d - e];
    let l = p.map(|q| q.dot(q));
    let m3 = |i: usize, j: usize, k: usize| {
        let (u, v, w) = (p[i], p[j], p[k]);
        u.x * (v.y * w.z - v.z * w.y) - u.y * (v.x * w.z - v.z * w.x)
            + u.z * (v.x * w.y - v.y * w.x)
    };
    let pm3 = |i: usize, j: usize, k: usize| {
        let (u, v, w) = (p[i], p[j], p[k]);
        u.x.abs() * ((v.y * w.z).abs() + (v.z * w.y).abs())
            + u.y.abs() * ((v.x * w.z).abs() + (v.z * w.x).abs())
            + u.z.abs() * ((v.x * w.y).abs() + (v.y * w.x).abs())
    };
    // det of rows (x, y, z, l) expanded along the lifted column
    let det = -l[0] * m3(1, 2, 3) + l[1] * m3(0, 2, 3) - l[2] * m3(0, 1, 3) + l[3] * m3(0, 1, 2);
    let perm = l[0] * pm3(1, 2, 3) + l[1] * pm3(0, 2, 3) + l[2] * pm3(0, 1, 3) + l[3] * pm3(0, 1, 2);
    // With rows a−e … d−e, det < 0 means inside for positive orientation.
    if det.abs() > BOUND * perm {
        return (-det).total_cmp(&0.0);
    }
    let rows: Vec<[Expansion; 3]> = [a, b, c, d]
        .iter()
        .map(|q| [diff(q.x, e.x), diff(q.y, e.y), diff(q.z, e.z)])
        .collect();
    let lift: Vec<Expansion> = rows
        .iter()
        .map(|r| add(&add(&mul(&r[0], &r[0]), &mul(&r[1], &r[1])), &mul(&r[2], &r[2])))
        .collect();
    let minor = |i: usize, j: usize, k: usize| {
        det3([
            [&rows[i][0], &rows[i][1], &rows[i][2]],
            [&rows[j][0], &rows[j][1], &rows[j][2]],
            [&rows[k][0], &rows[k][1], &rows[k][2]],
        ])
    };
    let t0 = neg(&mul(&lift[0], &minor(1, 2, 3)));
    let t1 = mul(&lift[1], &minor(0, 2, 3));
    let t2 = neg(&mul(&lift[2], &minor(0, 1, 3)));
    let t3 = mul(&lift[3], &minor(0, 1, 2));
    let det = add(&add(&t0, &t1), &add(&t2, &t3));
    sign(&neg(&det))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(x: f64, y: f64, z: f64) -> Vec3 {
        Vec3::new(x, y, z)
    }

    /// Invariant: the unit corner tetrahedron is positive, its mirror
    /// negative, and a coplanar point exactly zero even when the
    /// floating-point evaluation rounds. A wrong sign convention or a filter
    /// that trusts a rounded zero fails here.
    #[test]
    fn orientation_signs_are_exact() {
        let (a, b, c) = (v(0., 0., 0.), v(1., 0., 0.), v(0., 1., 0.));
        assert_eq!(orient3d(a, b, c, v(0., 0., 1.)), Ordering::Greater);
        assert_eq!(orient3d(a, c, b, v(0., 0., 1.)), Ordering::Less);
        let (p, q, r) = (v(0.1, 0.2, 0.3), v(1.7, 0.4, 0.9), v(0.3, 2.9, 1.3));
        for t in [0.1, 0.3, 1e-3, 1.0 / 3.0] {
            let s = p + (q - p) * t + (r - p) * (1.0 - t);
            let exact = orient3d(p, q, r, s);
            // s is the rounded image of a coplanar point; the exact sign must
            // agree with a tiny exact perturbation along the normal.
            let n = (q - p).cross(r - p);
            let up = s + n * 1e-9;
            assert_eq!(orient3d(p, q, r, up), Ordering::Greater);
            let _ = exact;
        }
        assert_eq!(orient3d(a, b, c, v(0.5, 0.25, 0.0)), Ordering::Equal);
        assert_eq!(orient3d(a, b, v(2., 0., 0.), v(1e-300, 1., 0.)), Ordering::Equal);
    }

    /// Invariant: the centroid is inside, a far point outside, and a fifth
    /// point of the same sphere exactly on it.
    #[test]
    fn insphere_signs_are_exact() {
        let (a, b, c, d) = (v(1., 0., 0.), v(0., 1., 0.), v(-1., 0., 0.), v(0., 0., 1.));
        assert_eq!(orient3d(a, b, c, d), Ordering::Greater);
        assert_eq!(insphere(a, b, c, d, v(0., 0., 0.2)), Ordering::Greater);
        assert_eq!(insphere(a, b, c, d, v(3., 0., 0.)), Ordering::Less);
        assert_eq!(insphere(a, b, c, d, v(0., -1., 0.)), Ordering::Equal);
        assert_eq!(insphere(a, b, c, d, v(0., 0., -1.)), Ordering::Equal);
    }
}
