//! Edge curves and node distributions (BS-06).
//!
//! The two-sided stretching is Vinokur's (M. Vinokur, "On one-dimensional
//! stretching functions for finite-difference calculations", J. Comput.
//! Phys. 50, 1983): u(ξ) = ½ (1 + tanh(δ (ξ − ½)) / tanh(δ/2)) and
//! s = u / (A + (1 − A) u), with sinh δ / δ = 1 / (N √(Δs₀ Δs₁)) and
//! A = √(Δs₁ / Δs₀). The slopes match the requested spacings; the achieved
//! first and last cell lengths are reported, not assumed.

use mesh_core::block::{Curve, Distribution, EdgeEnd};
use mesh_core::geometry::Vec3;

/// A curve between two points, evaluated at the fraction of its length.
#[derive(Clone, Debug)]
pub enum Evaluator {
    Line {
        a: Vec3,
        b: Vec3,
    },
    Arc {
        centre: Vec3,
        radius: f64,
        e1: Vec3,
        e2: Vec3,
        angle: f64,
    },
    Spline(Spline),
}

/// A natural cubic spline in the chord-length parameter, with a table of
/// arc length against parameter for evaluation by length.
#[derive(Clone, Debug)]
pub struct Spline {
    knots: Vec<f64>,
    points: Vec<Vec3>,
    /// Second derivatives at the knots.
    second: Vec<Vec3>,
    /// (parameter, arc length) samples, increasing.
    table: Vec<(f64, f64)>,
}

/// Samples of the arc-length table per spline segment.
const SPLINE_SAMPLES: usize = 64;

impl Spline {
    fn new(points: Vec<Vec3>) -> Result<Self, String> {
        let n = points.len();
        let mut knots = vec![0.0; n];
        for i in 1..n {
            let h = (points[i] - points[i - 1]).length();
            if h.is_nan() || h <= 0.0 {
                return Err("a spline repeats a point".into());
            }
            knots[i] = knots[i - 1] + h;
        }
        // natural spline: tridiagonal system for the second derivatives
        let mut second = vec![Vec3::ZERO; n];
        if n > 2 {
            let m = n - 2;
            let (mut diag, mut rhs) = (vec![0.0; m], vec![Vec3::ZERO; m]);
            let mut upper = vec![0.0; m];
            for r in 0..m {
                let i = r + 1;
                let (h0, h1) = (knots[i] - knots[i - 1], knots[i + 1] - knots[i]);
                diag[r] = 2.0 * (h0 + h1);
                upper[r] = h1;
                rhs[r] =
                    ((points[i + 1] - points[i]) / h1 - (points[i] - points[i - 1]) / h0) * 6.0;
            }
            // Thomas algorithm (the matrix is symmetric, lower = upper shifted)
            for r in 1..m {
                let w = upper[r - 1] / diag[r - 1];
                diag[r] -= w * upper[r - 1];
                rhs[r] = rhs[r] - rhs[r - 1] * w;
            }
            let mut x = vec![Vec3::ZERO; m];
            x[m - 1] = rhs[m - 1] / diag[m - 1];
            for r in (0..m - 1).rev() {
                x[r] = (rhs[r] - x[r + 1] * upper[r]) / diag[r];
            }
            second[1..n - 1].copy_from_slice(&x);
        }
        let mut s = Self {
            knots,
            points,
            second,
            table: Vec::new(),
        };
        let total = s.knots[n - 1];
        let samples = SPLINE_SAMPLES * (n - 1);
        let mut prev = s.eval(0.0);
        let mut length = 0.0;
        s.table.push((0.0, 0.0));
        for q in 1..=samples {
            let t = total * q as f64 / samples as f64;
            let p = s.eval(t);
            length += (p - prev).length();
            prev = p;
            s.table.push((t, length));
        }
        Ok(s)
    }

    fn eval(&self, t: f64) -> Vec3 {
        let n = self.knots.len();
        let i = match self.knots.partition_point(|&k| k <= t) {
            0 => 0,
            k if k >= n => n - 2,
            k => k - 1,
        };
        let (t0, t1) = (self.knots[i], self.knots[i + 1]);
        let h = t1 - t0;
        let (a, b) = ((t1 - t) / h, (t - t0) / h);
        self.points[i] * a
            + self.points[i + 1] * b
            + (self.second[i] * (a * a * a - a) + self.second[i + 1] * (b * b * b - b))
                * (h * h / 6.0)
    }

    fn length(&self) -> f64 {
        self.table[self.table.len() - 1].1
    }

    fn at(&self, s: f64) -> Vec3 {
        let target = s * self.length();
        let k = self
            .table
            .partition_point(|&(_, l)| l < target)
            .clamp(1, self.table.len() - 1);
        let ((t0, l0), (t1, l1)) = (self.table[k - 1], self.table[k]);
        let t = t0 + (t1 - t0) * ((target - l0) / (l1 - l0)).clamp(0.0, 1.0);
        self.eval(t)
    }
}

fn unit(v: Vec3) -> Option<Vec3> {
    let l = v.length();
    (l > 0.0 && l.is_finite()).then(|| v / l)
}

impl Evaluator {
    pub fn new(a: Vec3, b: Vec3, curve: &Curve) -> Result<Self, String> {
        if (b - a).length() <= 0.0 {
            return Err("an edge joins two coincident vertices".into());
        }
        match curve {
            Curve::Line => Ok(Self::Line { a, b }),
            Curve::Spline { through } => {
                if through.iter().any(|p| !p.is_finite()) {
                    return Err("a spline point is not finite".into());
                }
                let mut points = Vec::with_capacity(through.len() + 2);
                points.push(a);
                points.extend_from_slice(through);
                points.push(b);
                Ok(Self::Spline(Spline::new(points)?))
            }
            &Curve::Arc { through: p } => {
                let ab = b - a;
                let ap = p - a;
                let n = ab.cross(ap);
                let n2 = n.dot(n);
                if n2.is_nan() || n2 <= 1e-24 * ab.dot(ab) * ap.dot(ap) {
                    return Err("an arc's through-point is on the line of its ends".into());
                }
                // circumcentre relative to a
                let rel = (n.cross(ab) * ap.dot(ap) + ap.cross(n) * ab.dot(ab)) / (2.0 * n2);
                let centre = a + rel;
                let radius = rel.length();
                let e1 = unit(a - centre).ok_or("degenerate arc")?;
                // a → p → b turns positively about (p − a) × (b − a)
                let axis = unit(ap.cross(ab)).ok_or("degenerate arc")?;
                let e2 = axis.cross(e1);
                let angle_of = |q: Vec3| {
                    let d = q - centre;
                    let t = d.dot(e2).atan2(d.dot(e1));
                    if t < 0.0 {
                        t + std::f64::consts::TAU
                    } else {
                        t
                    }
                };
                let (tp, tb) = (angle_of(p), angle_of(b));
                if !(tp > 0.0 && tp < tb) {
                    return Err("an arc's through-point is not between its ends".into());
                }
                Ok(Self::Arc {
                    centre,
                    radius,
                    e1,
                    e2,
                    angle: tb,
                })
            }
        }
    }

    pub fn length(&self) -> f64 {
        match *self {
            Self::Line { a, b } => (b - a).length(),
            Self::Arc { radius, angle, .. } => radius * angle,
            Self::Spline(ref s) => s.length(),
        }
    }

    /// The point at fraction `s` of the length; exact at 0 and 1 only up to
    /// rounding, so callers place the end vertices themselves.
    pub fn at(&self, s: f64) -> Vec3 {
        match *self {
            Self::Line { a, b } => a + (b - a) * s,
            Self::Arc {
                centre,
                radius,
                e1,
                e2,
                angle,
            } => {
                let t = s * angle;
                centre + (e1 * t.cos() + e2 * t.sin()) * radius
            }
            Self::Spline(ref sp) => sp.at(s),
        }
    }

    /// The unit tangent at the start (`s = 0`) pointing along the curve,
    /// or at the end pointing back along it.
    pub fn end_tangent(&self, at_start: bool) -> Vec3 {
        let eps = 1e-6;
        let (p, q) = if at_start {
            (self.at(0.0), self.at(eps))
        } else {
            (self.at(1.0), self.at(1.0 - eps))
        };
        let d = q - p;
        d / d.length()
    }
}

/// Fractions s₀ = 0 < s₁ < … < s_n = 1 of the edge length for `n` cells.
pub fn fractions(law: Distribution, n: usize, length: f64) -> Result<Vec<f64>, String> {
    if n == 0 {
        return Err("an edge needs at least one cell".into());
    }
    let uniform = || (0..=n).map(|k| k as f64 / n as f64).collect::<Vec<_>>();
    let s = match law {
        Distribution::Uniform => uniform(),
        Distribution::Geometric { spacing, at } => {
            let h = positive(spacing, "geometric spacing")? / length;
            if h >= 1.0 {
                return Err(format!(
                    "the first spacing {spacing:e} is not shorter than the edge ({length:e})"
                ));
            }
            let r = growth_ratio(h, n)?;
            let mut s = Vec::with_capacity(n + 1);
            let (mut x, mut d) = (0.0, h);
            s.push(0.0);
            for _ in 0..n {
                x += d;
                d *= r;
                s.push(x);
            }
            let total = x;
            let mut s: Vec<f64> = s.into_iter().map(|v| v / total).collect();
            if at == EdgeEnd::End {
                s = s.iter().rev().map(|v| 1.0 - v).collect();
            }
            s
        }
        Distribution::Tanh { first, last } => {
            let d0 = positive(first, "first spacing")? / length;
            let d1 = positive(last, "last spacing")? / length;
            if d0 + d1 >= 1.0 {
                return Err("the end spacings are longer than the edge together".into());
            }
            let b = 1.0 / (n as f64 * (d0 * d1).sqrt());
            let a = (d1 / d0).sqrt();
            let u: Box<dyn Fn(f64) -> f64> = if (b - 1.0).abs() < 1e-9 {
                Box::new(|x| x)
            } else if b > 1.0 {
                let delta = solve_increasing(|d| d.sinh() / d, b, 1e-12, 1.0)?;
                Box::new(move |x| 0.5 * (1.0 + (delta * (x - 0.5)).tanh() / (0.5 * delta).tanh()))
            } else {
                // sin δ / δ decreases from 1 to 0 on (0, π)
                let delta =
                    solve_increasing(|d| -(d.sin() / d), -b, 1e-12, std::f64::consts::PI - 1e-12)?;
                Box::new(move |x| 0.5 * (1.0 + (delta * (x - 0.5)).tan() / (0.5 * delta).tan()))
            };
            (0..=n)
                .map(|k| {
                    if k == 0 {
                        0.0
                    } else if k == n {
                        1.0
                    } else {
                        let u = u(k as f64 / n as f64);
                        u / (a + (1.0 - a) * u)
                    }
                })
                .collect()
        }
    };
    if s.windows(2).any(|w| w[1] <= w[0]) || s.iter().any(|v| !v.is_finite()) {
        return Err("the distribution is not strictly increasing at this precision".into());
    }
    let mut s = s;
    s[0] = 0.0;
    s[n] = 1.0;
    Ok(s)
}

fn positive(v: f64, what: &str) -> Result<f64, String> {
    if v.is_finite() && v > 0.0 {
        Ok(v)
    } else {
        Err(format!("{what} must be positive and finite"))
    }
}

/// The ratio r with h (rⁿ − 1)/(r − 1) = 1 (fractions of the length).
fn growth_ratio(h: f64, n: usize) -> Result<f64, String> {
    let sum = |r: f64| (0..n).fold((0.0, 1.0), |(s, p), _| (s + p, p * r)).0 * h;
    if (sum(1.0) - 1.0).abs() < 1e-14 {
        return Ok(1.0);
    }
    let mut hi = 2.0;
    while sum(hi) < 1.0 {
        hi *= 2.0;
        if hi > 1e6 {
            return Err("the spacing is too small for the cell count".into());
        }
    }
    solve_increasing(sum, 1.0, 1e-300, hi)
}

/// Bisection for f(x) = target with f increasing on [lo, hi].
fn solve_increasing(f: impl Fn(f64) -> f64, target: f64, lo: f64, hi: f64) -> Result<f64, String> {
    let (mut lo, mut hi) = (lo, hi);
    let mut hi_v = f(hi);
    while hi_v < target {
        hi *= 2.0;
        hi_v = f(hi);
        if !hi_v.is_finite() {
            return Err("no stretching parameter matches the spacings".into());
        }
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Invariant: a geometric law gives the requested first spacing and a
    /// constant ratio. A wrong ratio solve breaks it; the grid tests would
    /// only see a different but valid grid.
    #[test]
    fn geometric_hits_the_first_spacing() {
        let s = fractions(
            Distribution::Geometric {
                spacing: 0.01,
                at: EdgeEnd::Start,
            },
            20,
            2.0,
        )
        .unwrap();
        assert!((s[1] * 2.0 - 0.01).abs() < 1e-12);
        let r0 = (s[2] - s[1]) / s[1];
        let r1 = (s[20] - s[19]) / (s[19] - s[18]);
        assert!((r0 - r1).abs() < 1e-9);
        let e = fractions(
            Distribution::Geometric {
                spacing: 0.01,
                at: EdgeEnd::End,
            },
            20,
            2.0,
        )
        .unwrap();
        assert!(((1.0 - e[19]) * 2.0 - 0.01).abs() < 1e-12);
    }

    /// Invariant: Vinokur's law clusters at both ends close to the requested
    /// spacings and stays monotone, including the weak-stretching branch.
    #[test]
    fn tanh_clusters_both_ends() {
        let s = fractions(
            Distribution::Tanh {
                first: 0.001,
                last: 0.004,
            },
            40,
            1.0,
        )
        .unwrap();
        assert!((s[1] / 0.001 - 1.0).abs() < 0.2, "first {}", s[1]);
        assert!(
            ((1.0 - s[39]) / 0.004 - 1.0).abs() < 0.2,
            "last {}",
            1.0 - s[39]
        );
        let weak = fractions(
            Distribution::Tanh {
                first: 0.04,
                last: 0.04,
            },
            10,
            1.0,
        )
        .unwrap();
        assert!(weak[1] > 0.04 && weak[1] < 0.1);
    }

    #[test]
    fn arc_through_a_quarter() {
        let e = Evaluator::new(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            &Curve::Arc {
                through: Vec3::new(0.5f64.sqrt(), 0.5f64.sqrt(), 0.0),
            },
        )
        .unwrap();
        assert!((e.length() - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
        let m = e.at(0.5);
        assert!((m.length() - 1.0).abs() < 1e-12 && (m.x - m.y).abs() < 1e-12);
        assert!(
            Evaluator::new(
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                &Curve::Arc {
                    through: Vec3::new(2.0, -1.0, 0.0)
                }
            )
            .is_err()
        );
    }

    /// Invariant: a spline through points of a circle follows the circle
    /// closely and is evaluated by arc length.
    #[test]
    fn spline_through_a_half_circle() {
        let through: Vec<Vec3> = (1..40)
            .map(|i| {
                let t = std::f64::consts::PI * i as f64 / 40.0;
                Vec3::new(t.cos(), t.sin(), 0.0)
            })
            .collect();
        let e = Evaluator::new(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
            &Curve::Spline { through },
        )
        .unwrap();
        assert!((e.length() - std::f64::consts::PI).abs() < 1e-3);
        for k in 1..10 {
            let p = e.at(k as f64 / 10.0);
            assert!((p.length() - 1.0).abs() < 1e-4, "{k}: {}", p.length());
            let angle = p.y.atan2(p.x) / std::f64::consts::PI;
            assert!((angle - k as f64 / 10.0).abs() < 1e-3);
        }
        let t = e.end_tangent(true);
        assert!(t.y > 0.999);
    }
}
