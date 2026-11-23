//! Small convex polytopes with tagged faces for the restricted Voronoi core:
//! a tetrahedron (or box) cut by half-spaces, each face remembering the
//! plane it came from.
use mesh_core::geometry::Vec3;

/// The plane a face of a piece lies on.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Src {
    /// Face of the scaffold tetrahedron opposite its local vertex.
    Tet(u8),
    /// Bisector towards another site.
    Bis(u32),
    /// Any other cutting plane (kernels).
    Cut,
}

#[derive(Clone, Debug)]
pub(super) struct Convex {
    pub pts: Vec<Vec3>,
    /// Outward counter-clockwise rings of local vertex indices.
    pub faces: Vec<(Src, Vec<usize>)>,
}

pub(super) enum Cut {
    Inside,
    Clipped,
    Empty,
}

impl Convex {
    /// A positively oriented tetrahedron; face `Tet(k)` is opposite vertex k.
    pub fn tet(p: [Vec3; 4]) -> Self {
        Self {
            pts: p.to_vec(),
            faces: vec![
                (Src::Tet(3), vec![0, 2, 1]),
                (Src::Tet(2), vec![0, 1, 3]),
                (Src::Tet(0), vec![1, 2, 3]),
                (Src::Tet(1), vec![2, 0, 3]),
            ],
        }
    }

    pub fn cube(min: Vec3, max: Vec3) -> Self {
        let pts = vec![
            min,
            Vec3::new(max.x, min.y, min.z),
            Vec3::new(max.x, max.y, min.z),
            Vec3::new(min.x, max.y, min.z),
            Vec3::new(min.x, min.y, max.z),
            Vec3::new(max.x, min.y, max.z),
            max,
            Vec3::new(min.x, max.y, max.z),
        ];
        let faces = [
            [0, 4, 7, 3],
            [1, 2, 6, 5],
            [0, 1, 5, 4],
            [3, 7, 6, 2],
            [0, 3, 2, 1],
            [4, 5, 6, 7],
        ]
        .into_iter()
        .map(|f| (Src::Cut, f.to_vec()))
        .collect();
        Self { pts, faces }
    }

    /// Keep `n·x ≤ c`; vertices within `tol` count as on the plane.
    pub fn clip(&mut self, n: Vec3, c: f64, src: Src, tol: f64) -> Result<Cut, &'static str> {
        let d: Vec<f64> = self.pts.iter().map(|&p| n.dot(p) - c).collect();
        if d.iter().any(|x| !x.is_finite()) {
            return Err("nonfinite clipping distance");
        }
        if d.iter().all(|&x| x <= tol) {
            return Ok(Cut::Inside);
        }
        if d.iter().all(|&x| x >= -tol) {
            return Ok(Cut::Empty);
        }
        let mut pts = Vec::with_capacity(self.pts.len() + 4);
        let mut on = Vec::with_capacity(self.pts.len() + 4);
        let mut map = vec![usize::MAX; self.pts.len()];
        for (i, &x) in d.iter().enumerate() {
            if x <= tol {
                map[i] = pts.len();
                pts.push(self.pts[i]);
                on.push(x >= -tol);
            }
        }
        let mut crossings: Vec<((usize, usize), usize)> = Vec::new();
        let mut faces = Vec::with_capacity(self.faces.len() + 1);
        for (s, ring) in &self.faces {
            if !ring.iter().any(|&v| d[v] < -tol) {
                continue;
            }
            let mut out = Vec::with_capacity(ring.len() + 1);
            for k in 0..ring.len() {
                let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
                if d[a] <= tol {
                    out.push(map[a]);
                }
                if (d[a] < -tol && d[b] > tol) || (d[a] > tol && d[b] < -tol) {
                    let key = (a.min(b), a.max(b));
                    let id = if let Some(&(_, id)) = crossings.iter().find(|(k, _)| *k == key) {
                        id
                    } else {
                        let (lo, hi) = key;
                        let t = d[lo] / (d[lo] - d[hi]);
                        let id = pts.len();
                        pts.push(self.pts[lo] + (self.pts[hi] - self.pts[lo]) * t);
                        on.push(true);
                        crossings.push((key, id));
                        id
                    };
                    out.push(id);
                }
            }
            if out.len() >= 3 {
                faces.push((*s, out));
            }
        }
        // Cap: in-plane edges of the kept faces, reversed.
        let mut exposed: Vec<(usize, usize)> = Vec::new();
        for (_, ring) in &faces {
            for k in 0..ring.len() {
                let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
                if on[a] && on[b] {
                    if let Some(i) = exposed.iter().position(|&e| e == (a, b)) {
                        exposed.swap_remove(i);
                    } else {
                        exposed.push((b, a));
                    }
                }
            }
        }
        // remove opposite pairs (an in-plane edge shared by two kept faces)
        let mut i = 0;
        while i < exposed.len() {
            let (a, b) = exposed[i];
            if let Some(j) = exposed.iter().position(|&e| e == (b, a)) {
                exposed.swap_remove(j.max(i));
                exposed.swap_remove(j.min(i));
            } else {
                i += 1;
            }
        }
        if exposed.len() < 3 {
            return Err("degenerate clipping cap");
        }
        let start = exposed[0].0;
        let mut cap = Vec::with_capacity(exposed.len());
        let mut at = start;
        for _ in 0..exposed.len() {
            cap.push(at);
            let Some(&(_, next)) = exposed.iter().find(|e| e.0 == at) else {
                return Err("open clipping cap");
            };
            at = next;
            if at == start {
                break;
            }
        }
        if at != start || cap.len() != exposed.len() {
            return Err("clipping cap is not one loop");
        }
        faces.push((src, cap));
        // drop unreferenced vertices
        let mut used = vec![usize::MAX; pts.len()];
        let mut kept = Vec::with_capacity(pts.len());
        for (_, ring) in &mut faces {
            for v in ring.iter_mut() {
                if used[*v] == usize::MAX {
                    used[*v] = kept.len();
                    kept.push(pts[*v]);
                }
                *v = used[*v];
            }
        }
        self.pts = kept;
        self.faces = faces;
        Ok(Cut::Clipped)
    }

    /// Volume, density-weighted mass and first moment for a density that is
    /// linear over the piece (`rho` evaluated at points).
    pub fn moments(&self, rho: &dyn Fn(Vec3) -> f64) -> (f64, f64, Vec3) {
        let o = self.pts[0];
        let ro = rho(o);
        let (mut vol, mut mass, mut moment) = (0.0, 0.0, Vec3::ZERO);
        for (_, ring) in &self.faces {
            if ring.contains(&0) {
                continue;
            }
            let a = self.pts[ring[0]];
            let ra = rho(a);
            for k in 1..ring.len() - 1 {
                let (b, c) = (self.pts[ring[k]], self.pts[ring[k + 1]]);
                let v = (a - o).cross(b - o).dot(c - o) / 6.0;
                let (rb, rc) = (rho(b), rho(c));
                let sr = ro + ra + rb + rc;
                vol += v;
                mass += v * sr / 4.0;
                let weighted = o * ro + a * ra + b * rb + c * rc;
                moment = moment + (weighted + (o + a + b + c) * sr) * (v / 20.0);
            }
        }
        (vol, mass, moment)
    }

    pub fn volume(&self) -> f64 {
        self.moments(&|_| 1.0).0
    }
}

/// Newell normal (area-weighted, twice the area) and vertex mean of a ring.
pub(super) fn ring_plane(points: &[Vec3], ring: &[usize]) -> (Vec3, Vec3) {
    // About the first vertex: far from the origin the products of absolute
    // coordinates lose the normal of a long thin face to rounding.
    let origin = ring.first().map_or(Vec3::ZERO, |&v| points[v]);
    let mut n = Vec3::ZERO;
    let mut c = Vec3::ZERO;
    for k in 0..ring.len() {
        let a = points[ring[k]] - origin;
        let b = points[ring[(k + 1) % ring.len()]] - origin;
        n = n + Vec3::new(
            (a.y - b.y) * (a.z + b.z),
            (a.z - b.z) * (a.x + b.x),
            (a.x - b.x) * (a.y + b.y),
        );
        c = c + a;
    }
    (n, origin + c / ring.len().max(1) as f64)
}

/// Unit normal of the least-squares plane through the vertex mean `centre`
/// of a ring, starting from its Newell normal `newell`.
///
/// The Newell normal of a long thin face (two cells that barely touch) is
/// tilted by rounding over the face's width, which shows as offsets at the
/// far ends although the vertices are on one plane. The fit removes the
/// tilt: offsets from the fitted plane are the face's true flatness.
pub(super) fn fitted_normal(points: &[Vec3], ring: &[usize], newell: Vec3, centre: Vec3) -> Vec3 {
    let unit = newell / newell.length();
    // in-plane axes, the first along the longest reach from the centre
    let far = ring
        .iter()
        .map(|&v| points[v] - centre)
        .map(|r| r - unit * r.dot(unit))
        .fold(
            Vec3::ZERO,
            |best, r| if r.dot(r) > best.dot(best) { r } else { best },
        );
    let reach = far.length();
    if reach.is_nan() || reach <= 0.0 {
        return unit;
    }
    let ex = far / reach;
    let ey = unit.cross(ex);
    let (mut sxx, mut sxy, mut syy, mut sxh, mut syh) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for &v in ring {
        let r = points[v] - centre;
        let (x, y, h) = (r.dot(ex), r.dot(ey), r.dot(unit));
        sxx += x * x;
        sxy += x * y;
        syy += y * y;
        sxh += x * h;
        syh += y * h;
    }
    let det = sxx * syy - sxy * sxy;
    // a ring without width fixes the tilt along its length only
    let (tx, ty) = if det > 1e-24 * sxx * sxx {
        ((sxh * syy - syh * sxy) / det, (syh * sxx - sxh * sxy) / det)
    } else {
        (sxh / sxx, 0.0)
    };
    let fitted = unit - ex * tx - ey * ty;
    let length = fitted.length();
    if length.is_finite() && length > 0.0 {
        fitted / length
    } else {
        unit
    }
}

/// A point strictly inside every face plane of a closed shell (its kernel),
/// or `None` when the shell is not star-shaped.
pub(super) fn kernel_point(points: &[Vec3], faces: &[Vec<usize>]) -> Option<Vec3> {
    let mut min = Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut max = min * -1.0;
    for &v in faces.iter().flatten() {
        let p = points[v];
        min = Vec3::new(min.x.min(p.x), min.y.min(p.y), min.z.min(p.z));
        max = Vec3::new(max.x.max(p.x), max.y.max(p.y), max.z.max(p.z));
    }
    let size = (max - min).max_abs_component();
    let pad = Vec3::new(size, size, size) * 0.01;
    let mut k = Convex::cube(min - pad, max + pad);
    let margin = size * 1e-7;
    for ring in faces {
        let (n, c) = ring_plane(points, ring);
        let len = n.length();
        if len.is_nan() || len <= 0.0 {
            continue;
        }
        let n = n / len;
        match k.clip(n, n.dot(c) - margin, Src::Cut, size * 1e-12).ok()? {
            Cut::Empty => return None,
            Cut::Inside | Cut::Clipped => {}
        }
    }
    let (vol, _, m) = k.moments(&|_| 1.0);
    (vol > size * size * size * 1e-18).then(|| m / vol)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn halving_a_cube_keeps_half_the_volume_and_tags_the_cut() {
        let mut c = Convex::cube(Vec3::ZERO, Vec3::new(1., 1., 1.));
        let r = c
            .clip(Vec3::new(1., 0., 0.), 0.5, Src::Bis(7), 1e-12)
            .unwrap();
        assert!(matches!(r, Cut::Clipped));
        assert!((c.volume() - 0.5).abs() < 1e-14);
        assert!(c.faces.iter().any(|(s, _)| *s == Src::Bis(7)));
        let (_, m, x) = c.moments(&|_| 1.0);
        assert!(((x / m).x - 0.25).abs() < 1e-14);
    }

    #[test]
    fn linear_density_moves_the_centroid() {
        let c = Convex::tet([
            Vec3::ZERO,
            Vec3::new(1., 0., 0.),
            Vec3::new(0., 1., 0.),
            Vec3::new(0., 0., 1.),
        ]);
        assert!((c.volume() - 1. / 6.).abs() < 1e-15);
        let (_, mass, m) = c.moments(&|p| p.x);
        // ∫x = 1/24, ∫x² = 1/60 over the unit tetrahedron
        assert!((mass - 1. / 24.).abs() < 1e-15);
        assert!((m.x - 1. / 60.).abs() < 1e-15);
    }
}
