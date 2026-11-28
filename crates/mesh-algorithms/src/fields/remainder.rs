//! The closed surface of the region behind the cuts of a sweep fill: the
//! STL without the swept ducts, closed at every cut by the cut's cap
//! ([`CutCap`]), so the polyhedral mesher can fill it (MD-25, BS-17).
//!
//! At each cut the STL facets that cross the cut plane next to the cut's
//! rim form a band. The band and everything the band separates on the
//! sweep's side (found by a flood over facet edges) are dropped. The hole
//! left in the STL is joined to the cap's rim (the wall nodes of the cut,
//! which lie on the STL in the cut plane) by a strip of triangles, zipped
//! along both loops by the shorter diagonal (Fuchs, Kedem & Uselton 1977,
//! *Optimal surface reconstruction from planar contours*, CACM 20(10), in
//! its greedy form).

use super::sweep::CutCap;
use mesh_core::geometry::Vec3;
use std::collections::HashMap;

/// Where a triangle of the remainder surface comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Source {
    /// An STL facet (its index).
    Facet(usize),
    /// The strip at a cut, next to this STL facet.
    Strip(usize),
    /// The cap of a cut (its index in the caps given).
    Cap(usize),
}

/// A closed surface: the STL behind the cuts and the caps.
#[derive(Clone, Debug)]
pub struct RemainderSurface {
    pub points: Vec<Vec3>,
    pub triangles: Vec<[usize; 3]>,
    pub source: Vec<Source>,
    /// The STL facets dropped (the swept ducts and the bands).
    pub dropped: usize,
}

fn signed(c: Vec3, n: Vec3, x: Vec3) -> f64 {
    (x - c).dot(n)
}

/// Cut the swept ducts off a closed, outward-oriented triangle surface and
/// close it with the caps.
pub fn remainder_surface(
    points: &[Vec3],
    triangles: &[[usize; 3]],
    caps: &[CutCap],
) -> Result<RemainderSurface, String> {
    if caps.is_empty() {
        return Err("no cut to close".into());
    }
    // facets by directed edge
    let mut by_edge: HashMap<(usize, usize), usize> = HashMap::with_capacity(3 * triangles.len());
    for (f, t) in triangles.iter().enumerate() {
        for e in 0..3 {
            if by_edge.insert((t[e], t[(e + 1) % 3]), f).is_some() {
                return Err("the surface has a directed edge twice".into());
            }
        }
    }
    let across = |f: usize, e: usize| -> Option<usize> {
        let t = triangles[f];
        by_edge.get(&(t[(e + 1) % 3], t[e])).copied()
    };
    // every cap's band first, so the flood of a duct cut at both ends (a
    // tube between two headers) stops at the band of its other end
    let mut band = vec![false; triangles.len()];
    let mut rings: Vec<(Vec<usize>, Vec3, Vec3, f64)> = Vec::with_capacity(caps.len());
    for cap in caps {
        let (c, n) = (cap.center, cap.normal);
        let reach = cap.points[..cap.rim]
            .iter()
            .map(|&p| (p - c).length())
            .fold(0.0, f64::max);
        // the band: the ring of facets across the plane around the rim. The
        // facets across the plane, joined over the edges across it, are the
        // curves where the plane meets the STL; the band is the one closest
        // to the cut's centre (another duct may cross the plane further out)
        // an STL node within a quarter of the rim spacing in front of the
        // plane counts as the sweep's side: a ring of STL nodes on (or
        // across, the cut being a level set) the plane would leave a hole
        // that zigzags through it and a strip of needles to the rim, whose
        // boundary recovery in the scaffold splits without end
        let rim_spacing = (0..cap.rim)
            .map(|k| (cap.points[k] - cap.points[(k + 1) % cap.rim]).length())
            .sum::<f64>()
            / cap.rim.max(1) as f64;
        let slack = 0.25 * rim_spacing;
        let ahead = |x: Vec3| signed(c, n, x) > -slack;
        let side = |v: usize| ahead(points[v]);
        let crossing = |f: usize| {
            let t = triangles[f];
            let s0 = side(t[0]);
            s0 != side(t[1]) || s0 != side(t[2])
        };
        let mut component = vec![usize::MAX; triangles.len()];
        let mut best: Option<(f64, Vec<usize>)> = None;
        for f0 in 0..triangles.len() {
            if component[f0] != usize::MAX || !crossing(f0) {
                continue;
            }
            component[f0] = f0;
            let mut members = vec![f0];
            let mut k = 0;
            while k < members.len() {
                let f = members[k];
                k += 1;
                let t = triangles[f];
                for e in 0..3 {
                    if side(t[e]) != side(t[(e + 1) % 3])
                        && let Some(g) = across(f, e)
                        && component[g] == usize::MAX
                    {
                        component[g] = f0;
                        members.push(g);
                    }
                }
            }
            // the curve's distance from the axis, and whether it stays next to the rim
            let mut far: f64 = 0.0;
            let mut fits = true;
            for &f in &members {
                let t = triangles[f];
                let longest = (0..3)
                    .map(|e| (points[t[e]] - points[t[(e + 1) % 3]]).length())
                    .fold(0.0, f64::max);
                let g = (points[t[0]] + points[t[1]] + points[t[2]]) / 3.0;
                let r = g - c;
                let radial = (r - n * r.dot(n)).length();
                far = far.max(radial);
                fits &= radial <= reach + longest;
            }
            if fits && best.as_ref().is_none_or(|(b, _)| far < *b) {
                best = Some((far, members));
            }
        }
        let Some((_, ring)) = best else {
            return Err("no STL facet crosses a cut next to its rim".into());
        };
        for &f in &ring {
            band[f] = true;
        }
        rings.push((ring, c, n, slack));
    }
    let mut dropped = vec![false; triangles.len()];
    for (ring, c, n, slack) in &rings {
        let ahead = |x: Vec3| signed(*c, *n, x) > -slack;
        // flood the sweep's side from the band
        let mut stack: Vec<usize> = Vec::new();
        for &f in ring {
            dropped[f] = true;
            for e in 0..3 {
                if let Some(g) = across(f, e)
                    && !band[g]
                    && !dropped[g]
                    && triangles[g].iter().all(|&v| ahead(points[v]))
                {
                    dropped[g] = true;
                    stack.push(g);
                }
            }
        }
        while let Some(f) = stack.pop() {
            for e in 0..3 {
                if let Some(g) = across(f, e)
                    && !band[g]
                    && !dropped[g]
                {
                    dropped[g] = true;
                    stack.push(g);
                }
            }
        }
    }
    // the holes: directed edges of kept facets with a dropped neighbour
    let mut next: HashMap<usize, (usize, usize)> = HashMap::new();
    for (f, t) in triangles.iter().enumerate() {
        if dropped[f] {
            continue;
        }
        for e in 0..3 {
            match across(f, e) {
                Some(g) if !dropped[g] => {}
                _ => {
                    // the strip runs this edge backwards
                    if next.insert(t[(e + 1) % 3], (t[e], f)).is_some() {
                        return Err("a hole left by the cut is not a simple loop".into());
                    }
                }
            }
        }
    }
    let mut holes: Vec<Vec<(usize, usize)>> = Vec::new();
    let mut seen: HashMap<usize, bool> = HashMap::new();
    let mut starts: Vec<usize> = next.keys().copied().collect();
    starts.sort_unstable();
    for s in starts {
        if seen.contains_key(&s) {
            continue;
        }
        let mut hole = Vec::new();
        let mut v = s;
        loop {
            seen.insert(v, true);
            let Some(&(w, f)) = next.get(&v) else {
                return Err("a hole left by the cut is open".into());
            };
            hole.push((v, f));
            v = w;
            if v == s {
                break;
            }
            if hole.len() > next.len() {
                return Err("a hole left by the cut does not close".into());
            }
        }
        holes.push(hole);
    }
    if holes.len() != caps.len() {
        return Err(format!(
            "the cuts leave {} holes in the STL, expected {}",
            holes.len(),
            caps.len()
        ));
    }
    // the output: kept STL points, then the caps' points
    let mut map = vec![usize::MAX; points.len()];
    let mut out_points: Vec<Vec3> = Vec::new();
    let mut out_triangles: Vec<[usize; 3]> = Vec::new();
    let mut source: Vec<Source> = Vec::new();
    for (f, t) in triangles.iter().enumerate() {
        if dropped[f] {
            continue;
        }
        let ids = t.map(|v| {
            if map[v] == usize::MAX {
                map[v] = out_points.len();
                out_points.push(points[v]);
            }
            map[v]
        });
        out_triangles.push(ids);
        source.push(Source::Facet(f));
    }
    let mut used = vec![false; holes.len()];
    for (k, cap) in caps.iter().enumerate() {
        let base = out_points.len();
        out_points.extend_from_slice(&cap.points);
        for t in &cap.triangles {
            out_triangles.push(t.map(|v| base + v));
            source.push(Source::Cap(k));
        }
        // the cap's rim as the strip runs it (against the cap's own edges)
        let mut rim_next: HashMap<usize, usize> = HashMap::new();
        let mut cap_edges: HashMap<(usize, usize), ()> = HashMap::new();
        for t in &cap.triangles {
            for e in 0..3 {
                cap_edges.insert((t[e], t[(e + 1) % 3]), ());
            }
        }
        for t in &cap.triangles {
            for e in 0..3 {
                let (a, b) = (t[e], t[(e + 1) % 3]);
                if !cap_edges.contains_key(&(b, a)) {
                    if a >= cap.rim || b >= cap.rim {
                        return Err("a cap's boundary is not its rim".into());
                    }
                    rim_next.insert(b, a);
                }
            }
        }
        if rim_next.len() != cap.rim {
            return Err("a cap's rim is not one loop".into());
        }
        let mut q_loop = Vec::with_capacity(cap.rim);
        let mut v = 0;
        for _ in 0..cap.rim {
            q_loop.push(v);
            v = rim_next[&v];
        }
        if v != 0 {
            return Err("a cap's rim is not one loop".into());
        }
        // the hole of this cap: the one closest to its plane and centre
        let (c, n) = (cap.center, cap.normal);
        let h = (0..holes.len())
            .filter(|&h| !used[h])
            .min_by(|&a, &b| {
                let score = |h: usize| {
                    (holes[h].iter().fold(Vec3::ZERO, |s, &(v, _)| s + points[v])
                        / holes[h].len() as f64
                        - c)
                        .length()
                        + signed(c, n, points[holes[h][0].0]).abs()
                };
                score(a).total_cmp(&score(b))
            })
            .ok_or("a cut has no hole left")?;
        used[h] = true;
        // zip: P runs the hole as the strip runs it, Q the rim the other way
        // round, so both turn the same way
        let p: Vec<(usize, usize)> = holes[h].iter().map(|&(v, f)| (map[v], f)).collect();
        let q: Vec<usize> = q_loop.iter().rev().map(|&v| base + v).collect();
        for (t, f) in strip(&out_points, &p, &q, c, n) {
            out_triangles.push(t);
            source.push(Source::Strip(f));
        }
    }
    Ok(RemainderSurface {
        points: out_points,
        triangles: out_triangles,
        source,
        dropped: dropped.iter().filter(|&&d| d).count(),
    })
}

/// A rim node where the rim turns by more than this (degrees) is a corner
/// of the cut's section (a square duct's edge), matched to the hole node
/// on that edge.
const CORNER_TURN: f64 = 30.0;

/// The strip between a hole `p` (nodes with the STL facet beside each) and
/// a cap's rim `q`, both running the same way round the cut's axis
/// `(c, n)`: each triangle with the STL facet it lies next to.
///
/// The zip runs by the shorter diagonal (Fuchs, Kedem & Uselton 1977), and
/// a strip triangle must face away from the axis, as the duct wall: the
/// zip never takes the triangle that folds when the other does not (a hole
/// edge running steeply along the duct, from a coarse STL, puts the rim
/// node behind the line of that edge).
///
/// A section with corners (a square duct) is zipped side by side: each
/// rim corner is matched to the hole node on its edge, and each side runs
/// from corner to corner, so no triangle spans an edge of the duct (one
/// that did cut the corner off, and the layer stack squeezed the cap
/// unevenly there until its facets collided). A round section, or corners
/// that do not match, are zipped round the loop from the hole's first node
/// and the rim node nearest to it; its last steps are forced, so the first
/// start whose strip has no fold is taken (the hole's first node, then the
/// others in turn), else the start with the fewest folds.
fn strip(
    points: &[Vec3],
    p: &[(usize, usize)],
    q: &[usize],
    c: Vec3,
    n: Vec3,
) -> Vec<([usize; 3], usize)> {
    let at = |v: usize| points[v];
    let faces_out = |t: [usize; 3]| {
        let [a, b, d] = t.map(at);
        let g = (a + b + d) / 3.0 - c;
        (b - a).cross(d - a).dot(g - n * g.dot(n)) > 0.0
    };
    // zip a hole run onto a rim run, first and last nodes paired
    let zip = |pr: &[(usize, usize)], qr: &[usize]| -> (Vec<([usize; 3], usize)>, usize) {
        let (m, l) = (pr.len() - 1, qr.len() - 1);
        let mut out = Vec::with_capacity(m + l);
        let mut folds = 0;
        let (mut i, mut j) = (0, 0);
        while i < m || j < l {
            let (pi, qj) = (pr[i].0, qr[j]);
            let advance_p = if i == m {
                false
            } else if j == l {
                true
            } else {
                let (pn, qn) = (pr[i + 1].0, qr[j + 1]);
                match (faces_out([pi, pn, qj]), faces_out([pi, qn, qj])) {
                    (true, false) => true,
                    (false, true) => false,
                    _ => (at(pn) - at(qj)).length() <= (at(pi) - at(qn)).length(),
                }
            };
            let t = if advance_p {
                i += 1;
                [pi, pr[i].0, qj]
            } else {
                j += 1;
                [pi, qr[j], qj]
            };
            folds += usize::from(!faces_out(t));
            out.push((t, pr[i.min(m) - usize::from(advance_p)].1));
        }
        (out, folds)
    };
    let (np, nq) = (p.len(), q.len());
    let p_run = |s: usize, len: usize| -> Vec<(usize, usize)> {
        (0..=len).map(|i| p[(s + i) % np]).collect()
    };
    let q_run =
        |s: usize, len: usize| -> Vec<usize> { (0..=len).map(|j| q[(s + j) % nq]).collect() };
    if let Some(pairs) = corners(points, p, q, n) {
        let mut out = Vec::with_capacity(np + nq);
        for k in 0..pairs.len() {
            let ((pa, qa), (pb, qb)) = (pairs[k], pairs[(k + 1) % pairs.len()]);
            let (lp, lq) = ((pb + np - pa) % np, (qb + nq - qa) % nq);
            let (lp, lq) = if pairs.len() == 1 { (np, nq) } else { (lp, lq) };
            out.extend(zip(&p_run(pa, lp), &q_run(qa, lq)).0);
        }
        return out;
    }
    let nearest = |v: usize| {
        (0..nq)
            .min_by(|&a, &b| {
                (at(q[a]) - at(v))
                    .length()
                    .total_cmp(&(at(q[b]) - at(v)).length())
            })
            .unwrap_or(0)
    };
    let around = |s: usize| zip(&p_run(s, np), &q_run(nearest(p[s].0), nq));
    let mut best = around(0);
    for s in 1..np {
        if best.1 == 0 {
            break;
        }
        let other = around(s);
        if other.1 < best.1 {
            best = other;
        }
    }
    best.0
}

/// The corners of a cut's section paired with the hole: `(hole index, rim
/// index)` for every rim node where the rim turns by more than
/// [`CORNER_TURN`], with the hole node on that duct edge (in the plane's
/// view within a quarter of the rim spacing of the corner, the nearest to
/// the plane of those), in the order of both loops. `None` when the rim
/// has no corner or a corner has no hole node on its edge, or the pairs do
/// not run round both loops the same way once.
fn corners(
    points: &[Vec3],
    p: &[(usize, usize)],
    q: &[usize],
    n: Vec3,
) -> Option<Vec<(usize, usize)>> {
    let (np, nq) = (p.len(), q.len());
    if np < 3 || nq < 3 {
        return None;
    }
    let at = |v: usize| points[v];
    let flat = |v: Vec3| v - n * v.dot(n);
    let spacing = (0..nq)
        .map(|j| (at(q[(j + 1) % nq]) - at(q[j])).length())
        .sum::<f64>()
        / nq as f64;
    let cos_turn = CORNER_TURN.to_radians().cos();
    let mut pairs = Vec::new();
    for j in 0..nq {
        let (a, b, d) = (at(q[(j + nq - 1) % nq]), at(q[j]), at(q[(j + 1) % nq]));
        let (u, w) = (flat(b - a), flat(d - b));
        if u.length() == 0.0 || w.length() == 0.0 {
            return None;
        }
        if u.dot(w) / (u.length() * w.length()) > cos_turn {
            continue;
        }
        let i = (0..np)
            .filter(|&i| flat(at(p[i].0) - b).length() < 0.25 * spacing)
            .min_by(|&x, &y| {
                let axial = |i: usize| (at(p[i].0) - b).dot(n).abs();
                axial(x).total_cmp(&axial(y))
            })?;
        pairs.push((i, j));
    }
    if pairs.is_empty() {
        return None;
    }
    // the hole indices run round once, in the rim's order
    let mut turn = 0;
    for k in 0..pairs.len() {
        let (a, b) = (pairs[k].0, pairs[(k + 1) % pairs.len()].0);
        let step = (b + np - a) % np;
        if step == 0 && pairs.len() > 1 {
            return None;
        }
        turn += step;
    }
    (turn == np || pairs.len() == 1).then_some(pairs)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A closed box `lo`..`hi`, outward, each side two triangles.
    fn add_box(points: &mut Vec<Vec3>, triangles: &mut Vec<[usize; 3]>, lo: Vec3, hi: Vec3) {
        let b = points.len();
        for k in 0..8 {
            points.push(Vec3::new(
                if k & 1 == 0 { lo.x } else { hi.x },
                if k & 2 == 0 { lo.y } else { hi.y },
                if k & 4 == 0 { lo.z } else { hi.z },
            ));
        }
        for q in [
            [0, 2, 3, 1],
            [4, 5, 7, 6],
            [0, 1, 5, 4],
            [2, 6, 7, 3],
            [0, 4, 6, 2],
            [1, 3, 7, 5],
        ] {
            triangles.push([b + q[0], b + q[1], b + q[2]]);
            triangles.push([b + q[0], b + q[2], b + q[3]]);
        }
    }

    /// The cap of the square duct x, y in [0, 2] at z = 5, the sweep above:
    /// eight rim points (corners and edge midpoints) and a centre.
    fn cap() -> CutCap {
        let rim = [
            (0.0, 0.0),
            (1.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (2.0, 2.0),
            (1.0, 2.0),
            (0.0, 2.0),
            (0.0, 1.0),
        ];
        let mut points: Vec<Vec3> = rim.iter().map(|&(x, y)| Vec3::new(x, y, 5.0)).collect();
        points.push(Vec3::new(1.0, 1.0, 5.0));
        let triangles = (0..8).map(|i| [8, i, (i + 1) % 8]).collect();
        CutCap {
            end: 1,
            normal: Vec3::new(0.0, 0.0, 1.0),
            center: Vec3::new(1.0, 1.0, 5.0),
            radius: 1.0,
            points,
            rim: 8,
            triangles,
            node: (0..9).collect(),
        }
    }

    fn volume(s: &RemainderSurface) -> f64 {
        s.triangles
            .iter()
            .map(|t| s.points[t[0]].dot(s.points[t[1]].cross(s.points[t[2]])) / 6.0)
            .sum()
    }

    fn closed(s: &RemainderSurface) -> bool {
        let mut edges: HashMap<(usize, usize), usize> = HashMap::new();
        for t in &s.triangles {
            for e in 0..3 {
                *edges.entry((t[e], t[(e + 1) % 3])).or_default() += 1;
            }
        }
        edges
            .iter()
            .all(|(&(a, b), &n)| n == 1 && edges.get(&(b, a)) == Some(&1))
    }

    /// Keep: the cut drops the duct on the sweep's side of the cap and
    /// closes the rest with the cap and a strip, outward and watertight, so
    /// the polyhedral mesher can fill it; a wrong strip direction or a
    /// missed band facet leaves an open or folded surface, which the
    /// Manifold test only sees as a failed scaffold.
    #[test]
    fn cut_duct_is_closed_by_its_cap() {
        let (mut points, mut triangles) = (Vec::new(), Vec::new());
        add_box(
            &mut points,
            &mut triangles,
            Vec3::ZERO,
            Vec3::new(2.0, 2.0, 10.0),
        );
        let s = remainder_surface(&points, &triangles, &[cap()]).unwrap();
        assert!(closed(&s));
        assert!((volume(&s) - 20.0).abs() < 1e-12, "{}", volume(&s));
        // the top and the eight sides go, the bottom stays
        assert_eq!(s.dropped, 10);
        let caps = s
            .source
            .iter()
            .filter(|x| matches!(x, Source::Cap(0)))
            .count();
        assert_eq!(caps, 8);
    }

    /// Keep: another duct that crosses the cut plane next to the rim (here
    /// within one of its long facets) is not part of the cut: only the ring
    /// of crossing facets around the rim is the band. Taking every
    /// crossing facet within reach dropped the Manifold's second pipe and
    /// left more holes than cuts.
    #[test]
    fn a_duct_beside_the_cut_stays() {
        let (mut points, mut triangles) = (Vec::new(), Vec::new());
        add_box(
            &mut points,
            &mut triangles,
            Vec3::ZERO,
            Vec3::new(2.0, 2.0, 10.0),
        );
        add_box(
            &mut points,
            &mut triangles,
            Vec3::new(4.0, 0.0, 0.0),
            Vec3::new(6.0, 2.0, 10.0),
        );
        let s = remainder_surface(&points, &triangles, &[cap()]).unwrap();
        assert!(closed(&s));
        assert!((volume(&s) - 60.0).abs() < 1e-12, "{}", volume(&s));
        assert_eq!(s.dropped, 10);
    }

    /// Keep: a ring of STL nodes just behind the cut plane (0.01 here, a
    /// hundredth of the rim spacing) is taken with the sweep's side, so the
    /// strip reaches back to the ring before it and has no needle. The
    /// three-pipe box had its cut level set 0.003 across a ring of pipe
    /// nodes: the hole zigzagged through that ring, the strip held 179°
    /// triangles and the scaffold's boundary recovery split without end.
    #[test]
    fn a_ring_next_to_the_plane_goes_with_the_sweep() {
        let square = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)];
        let levels = [0.0, 4.99, 10.0];
        let mut points: Vec<Vec3> = Vec::new();
        for &z in &levels {
            for &(x, y) in &square {
                points.push(Vec3::new(x, y, z));
            }
        }
        let mut triangles: Vec<[usize; 3]> = vec![[0, 2, 1], [0, 3, 2], [8, 9, 10], [8, 10, 11]];
        for l in 0..2 {
            for i in 0..4 {
                let (a, b) = (4 * l + i, 4 * l + (i + 1) % 4);
                triangles.push([a, b, b + 4]);
                triangles.push([a, b + 4, a + 4]);
            }
        }
        let s = remainder_surface(&points, &triangles, &[cap()]).unwrap();
        assert!(closed(&s));
        assert!((volume(&s) - 20.0).abs() < 1e-12, "{}", volume(&s));
        for (t, src) in s.triangles.iter().zip(&s.source) {
            if matches!(src, Source::Strip(_)) {
                let low = t
                    .iter()
                    .map(|&v| s.points[v].z)
                    .fold(f64::INFINITY, f64::min);
                assert!(low < 1.0, "a strip triangle ends at z = {low}");
            }
        }
    }

    /// Keep: the strip between a coarse hole and a fine rim has no
    /// triangle folded inward. Here the hole's last edge runs steeply down
    /// the tube to its first node, and the rim node nearest that first node
    /// lies behind the edge: the zip from the first node ends in a forced
    /// triangle that faces the axis. A 5-tube header's polyhedral layers
    /// failed on such a fold ("no common inward wall direction").
    #[test]
    fn a_steep_hole_edge_does_not_fold_the_strip() {
        let (c, n) = (Vec3::new(0.0, 0.0, 5.0), Vec3::new(0.0, 0.0, 1.0));
        let at = |deg: f64, depth: f64| {
            let a = deg.to_radians();
            Vec3::new(a.cos(), a.sin(), 5.0 - depth)
        };
        let mut points = vec![at(0.0, 0.6)];
        for k in 1..12 {
            points.push(at(30.0 * k as f64, if k % 2 == 0 { 0.2 } else { 0.4 }));
        }
        points.push(at(359.0, 0.05));
        let np = points.len();
        for k in 0..16 {
            points.push(at(-8.0 + 22.5 * k as f64, 0.0));
        }
        let p: Vec<(usize, usize)> = (0..np).map(|i| (i, i)).collect();
        let q: Vec<usize> = (np..points.len()).collect();
        let tris = strip(&points, &p, &q, c, n);
        assert_eq!(tris.len(), np + q.len());
        for (t, _) in &tris {
            let [a, b, d] = t.map(|v| points[v]);
            let g = (a + b + d) / 3.0 - c;
            let out = (b - a).cross(d - a).dot(g - n * g.dot(n));
            assert!(out > 0.0, "the strip triangle {t:?} faces the axis");
        }
        // every hole edge and every rim edge once
        let mut edges: HashMap<(usize, usize), usize> = HashMap::new();
        for (t, _) in &tris {
            for e in 0..3 {
                *edges.entry((t[e], t[(e + 1) % 3])).or_default() += 1;
            }
        }
        for i in 0..np {
            assert_eq!(edges.get(&(i, (i + 1) % np)), Some(&1), "hole edge {i}");
        }
    }

    /// Keep: at a square cut no strip triangle spans an edge of the duct:
    /// each side of the rim is zipped to the hole from corner to corner.
    /// Here the hole's corner nodes lie deeper than their neighbours (as
    /// a Gmsh STL along a box edge), and the zip round the loop joined a
    /// hole node on one wall to rim nodes on the next. The cut corner was
    /// cut off, the layer stack squeezed the cap unevenly and its facets
    /// collided (square tee with a round branch, "surface collision").
    #[test]
    fn a_square_strip_keeps_the_duct_edges() {
        let (c, n) = (Vec3::new(1.0, 1.0, 5.0), Vec3::new(0.0, 0.0, 1.0));
        // a point on the square's perimeter, counter-clockwise from (0, 0)
        let around = |t: f64, z: f64| {
            let (side, u) = ((t / 2.0).floor() as i32 % 4, t % 2.0);
            let (x, y) = match side {
                0 => (u, 0.0),
                1 => (2.0, u),
                2 => (2.0 - u, 2.0),
                _ => (0.0, 2.0 - u),
            };
            Vec3::new(x, y, z)
        };
        let mut points = Vec::new();
        for k in 0..4 {
            let t = 2.0 * k as f64;
            points.push(around(t, 4.4));
            points.push(around(t + 0.1, 4.9));
            points.push(around(t + 1.0, 4.7));
            points.push(around(t + 1.9, 4.9));
        }
        let np = points.len();
        for k in 0..24 {
            points.push(around(8.0 * k as f64 / 24.0, 5.0));
        }
        let p: Vec<(usize, usize)> = (0..np).map(|i| (i, i)).collect();
        let q: Vec<usize> = (np..points.len()).collect();
        let tris = strip(&points, &p, &q, c, n);
        assert_eq!(tris.len(), np + q.len());
        let wall = |x: Vec3| -> Vec<u8> {
            let mut w = Vec::new();
            for (k, on) in [x.y < 1e-9, x.x > 2.0 - 1e-9, x.y > 2.0 - 1e-9, x.x < 1e-9]
                .into_iter()
                .enumerate()
            {
                if on {
                    w.push(k as u8);
                }
            }
            w
        };
        for (t, _) in &tris {
            let walls: Vec<Vec<u8>> = t.iter().map(|&v| wall(points[v])).collect();
            let shared = (0..4u8).any(|k| walls.iter().all(|w| w.contains(&k)));
            assert!(shared, "the strip triangle {t:?} spans a duct edge");
        }
    }
}
