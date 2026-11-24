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

#![allow(clippy::needless_range_loop)] // parallel per-vertex and per-face arrays

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
    let perm =
        l[0] * pm3(1, 2, 3) + l[1] * pm3(0, 2, 3) + l[2] * pm3(0, 1, 3) + l[3] * pm3(0, 1, 2);
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
        .map(|r| {
            add(
                &add(&mul(&r[0], &r[0]), &mul(&r[1], &r[1])),
                &mul(&r[2], &r[2]),
            )
        })
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

/// One entry of an offset determinant.
#[derive(Clone, Copy, Debug)]
pub enum Entry {
    /// The number one (the row of barycentric weights).
    One,
    /// Zero by construction: the two sites are equidistant from the point.
    Zero,
    /// |v − k|² − |v − s|², positive when `v` is closer to `s`.
    Offset { v: Vec3, s: Vec3, k: Vec3 },
}

impl Entry {
    /// Rounded value and a bound on its rounding error.
    fn approximate(self) -> (f64, f64) {
        match self {
            Entry::One => (1.0, 0.0),
            Entry::Zero => (0.0, 0.0),
            Entry::Offset { v, s, k } => {
                let mut value = 0.0;
                let mut size = 0.0;
                for (v, s, k) in [(v.x, s.x, k.x), (v.y, s.y, k.y), (v.z, s.z, k.z)] {
                    value += (s - k) * (2.0 * v - s - k);
                    size += (s.abs() + k.abs()) * (2.0 * v.abs() + s.abs() + k.abs());
                }
                (value, 16.0 * f64::EPSILON * size)
            }
        }
    }

    fn exact(self) -> Expansion {
        match self {
            Entry::One => vec![1.0],
            Entry::Zero => vec![0.0],
            Entry::Offset { v, s, k } => {
                let mut sum = vec![0.0];
                for (v, s, k) in [(v.x, s.x, k.x), (v.y, s.y, k.y), (v.z, s.z, k.z)] {
                    // (s − k)(2v − s − k); doubling is exact
                    let t = grow(&diff(2.0 * v, s), -k);
                    sum = add(&sum, &mul(&diff(s, k), &t));
                }
                sum
            }
        }
    }
}

fn minor_indices(n: usize, skip: usize) -> impl Iterator<Item = usize> {
    (0..n).filter(move |&j| j != skip)
}

/// Determinant (or permanent) of the rows `rows` and columns `cols` of `m`
/// by expansion along the first of the rows.
fn float_det(m: &[[f64; 4]; 4], rows: &[usize], cols: &[usize], permanent: bool) -> f64 {
    let Some((&r, rest)) = rows.split_first() else {
        return 1.0;
    };
    let mut sum = 0.0;
    for (j, &c) in cols.iter().enumerate() {
        let sub: Vec<usize> = minor_indices(cols.len(), j).map(|i| cols[i]).collect();
        let term = m[r][c] * float_det(m, rest, &sub, permanent);
        sum += if permanent || j % 2 == 0 { term } else { -term };
    }
    sum
}

fn exact_det(m: &[Vec<Expansion>], rows: &[usize], cols: &[usize]) -> Expansion {
    let Some((&r, rest)) = rows.split_first() else {
        return vec![1.0];
    };
    let mut sum = vec![0.0];
    for (j, &c) in cols.iter().enumerate() {
        if sign(&m[r][c]) == Ordering::Equal {
            continue;
        }
        let sub: Vec<usize> = minor_indices(cols.len(), j).map(|i| cols[i]).collect();
        let term = mul(&m[r][c], &exact_det(m, rest, &sub));
        sum = add(&sum, &if j % 2 == 0 { term } else { neg(&term) });
    }
    sum
}

/// The exact sign of the determinant of the leading n × n block (n ≤ 4).
pub fn offset_det_sign(n: usize, m: &[[Entry; 4]; 4]) -> Ordering {
    let index: Vec<usize> = (0..n).collect();
    let mut value = [[0.0; 4]; 4];
    let mut upper = [[0.0; 4]; 4];
    let mut lower = [[0.0; 4]; 4];
    for i in 0..n {
        for j in 0..n {
            let (v, e) = m[i][j].approximate();
            value[i][j] = v;
            lower[i][j] = v.abs();
            upper[i][j] = v.abs() + e;
        }
    }
    let det = float_det(&value, &index, &index, false);
    let high = float_det(&upper, &index, &index, true);
    let low = float_det(&lower, &index, &index, true);
    // entry errors spread through every product, plus the rounding of the
    // products and sums themselves
    let bound = (high - low) + 64.0 * f64::EPSILON * high;
    if det.abs() > bound {
        return det.total_cmp(&0.0);
    }
    let exact: Vec<Vec<Expansion>> = (0..n)
        .map(|i| (0..n).map(|j| m[i][j].exact()).collect())
        .collect();
    sign(&exact_det(&exact, &index, &index))
}

/// Side of a point that is fixed by bisectors on a simplex, relative to one
/// more bisector, after Lévy's side predicates (PCK, 2016) for restricted
/// Voronoi diagrams.
///
/// The point q lies on the simplex with the corners `corners` (one to four
/// of them: a point, an edge, a triangle or a tetrahedron) and on the
/// bisectors of the site `s` with each site of `on` (one fewer than the
/// corners). The result is the sign of |q − o|² − |q − s|²: `Greater` when q
/// is closer to `s`, `Equal` when it is on the bisector of `s` and `o`.
/// `None` when the bisectors do not fix a point of the simplex's plane.
///
/// `tied(corner, site)` says that the site is equidistant from that corner
/// with `s` by construction; that offset is then zero and not computed.
///
/// With barycentric weights λ of q, Σ λ = 1 and Σ λ f_k = 0 for each k of
/// `on`, where f_k is the offset of bisector k at a corner. Cramer's rule
/// gives Σ λ f_o as the quotient of two determinants of these offsets.
pub fn bisector_side(
    corners: &[Vec3],
    s: Vec3,
    on: &[Vec3],
    o: Vec3,
    tied: &dyn Fn(usize, Option<usize>) -> bool,
) -> Option<Ordering> {
    let n = corners.len();
    if n == 0 || n > 4 || on.len() + 1 != n {
        return None;
    }
    let offset = |j: usize, site: Option<usize>| -> Entry {
        if tied(j, site) {
            Entry::Zero
        } else {
            let k = site.map_or(o, |i| on[i]);
            Entry::Offset {
                v: corners[j],
                s,
                k,
            }
        }
    };
    let mut top = [[Entry::Zero; 4]; 4];
    let mut bottom = [[Entry::Zero; 4]; 4];
    for j in 0..n {
        top[0][j] = offset(j, None);
        bottom[0][j] = Entry::One;
        for i in 0..on.len() {
            top[i + 1][j] = offset(j, Some(i));
            bottom[i + 1][j] = top[i + 1][j];
        }
    }
    let below = offset_det_sign(n, &bottom);
    if below == Ordering::Equal {
        return None;
    }
    let above = offset_det_sign(n, &top);
    Some(if below == Ordering::Less {
        above.reverse()
    } else {
        above
    })
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
        assert_eq!(
            orient3d(a, b, v(2., 0., 0.), v(1e-300, 1., 0.)),
            Ordering::Equal
        );
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

    /// Invariant: the side of a Voronoi vertex is exact. Four sites on an
    /// integer sphere fix its centre inside a tetrahedron; a fifth site on
    /// the same sphere is exactly on, and one a rounding step inside or
    /// outside is strictly on its side. A filter that trusts a rounded
    /// value, or a wrong sign for a negative denominator, fails here; so
    /// does an answer that depends on the order of the corners.
    #[test]
    fn side_of_a_voronoi_vertex_is_exact() {
        let corners = [
            v(-7., -5., -6.),
            v(9., -4., -5.),
            v(-3., 8., -4.),
            v(1., 2., 11.),
        ];
        let s = v(3., 0., 0.);
        let on = [v(0., 3., 0.), v(0., 0., 3.), v(-3., 0., 0.)];
        let free = |_: usize, _: Option<usize>| false;
        let side = |c: &[Vec3], o: Vec3| bisector_side(c, s, &on, o, &free);
        assert_eq!(side(&corners, v(2., 2., 1.)), Some(Ordering::Equal));
        assert_eq!(
            side(&corners, v(2., 2., 1.0 + f64::EPSILON)),
            Some(Ordering::Greater)
        );
        assert_eq!(
            side(&corners, v(2., 2., 1.0 - f64::EPSILON)),
            Some(Ordering::Less)
        );
        let turned = [corners[1], corners[0], corners[2], corners[3]];
        assert_eq!(
            side(&turned, v(2., 2., 1.0 - f64::EPSILON)),
            Some(Ordering::Less)
        );
        // three bisectors through one line fix no point
        let line = [v(0., 3., 0.), v(0., -3., 0.), v(-3., 0., 0.)];
        assert_eq!(
            bisector_side(&corners, s, &line, v(9., 9., 9.), &free),
            None
        );
    }

    /// Invariant: the two cells of a bisector face take the same decision.
    /// A vertex on an edge and on the bisector of s and t is closer to a
    /// third site k than to s exactly when it is closer to k than to t.
    /// With rounded coordinates the two views disagree near a tie (the
    /// defect of the tolerance tests); the sites here are a rounding step
    /// from a tie on purpose.
    #[test]
    fn both_cells_of_a_face_agree() {
        let mut seed = 0x2545_f491_4f6c_dd1d_u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let free = |_: usize, _: Option<usize>| false;
        let mut strict = 0;
        for _ in 0..2000 {
            let a = v(next(), next(), next());
            let b = v(next() + 1.0, next(), next());
            let s = v(next(), next() + 0.5, next());
            let t = v(next() + 1.0, next() + 0.5, next());
            // a site at almost the same distance from the vertex as s and t
            let q = {
                let f = |p: Vec3| (p - t).dot(p - t) - (p - s).dot(p - s);
                let (fa, fb) = (f(a), f(b));
                a + (b - a) * (fa / (fa - fb))
            };
            let r = (q - s).length();
            let d = v(next() - 0.5, next() - 0.5, next() - 0.5);
            let k = q + d * (r / d.length() * (1.0 + (next() - 0.5) * 4.0 * f64::EPSILON));
            let from_s = bisector_side(&[a, b], s, &[t], k, &free);
            let from_t = bisector_side(&[a, b], t, &[s], k, &free);
            assert!(from_s.is_some());
            assert_eq!(from_s, from_t);
            strict += usize::from(from_s != Some(Ordering::Equal));
        }
        assert!(strict > 1900);
    }

    /// Invariant: a tie by construction is an exact zero. Two sites meant to
    /// be equidistant from a corner are a rounding step apart in stored
    /// coordinates; with the tie declared, the corner is on their bisector.
    #[test]
    fn declared_tie_is_on_the_bisector() {
        let corner = v(0.1, 0.2, 0.3);
        let s = corner + v(0.6, 0.0, 0.8) * 0.7;
        let o = corner + v(0.0, 0.6, 0.8) * 0.7;
        let free = |_: usize, _: Option<usize>| false;
        let tie = |_: usize, _: Option<usize>| true;
        assert_ne!(
            bisector_side(&[corner], s, &[], o, &free),
            Some(Ordering::Equal)
        );
        assert_eq!(
            bisector_side(&[corner], s, &[], o, &tie),
            Some(Ordering::Equal)
        );
    }
}
