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
    let mut dropped = vec![false; triangles.len()];
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
        let mut band = vec![false; triangles.len()];
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
        // flood the sweep's side from the band
        let mut stack: Vec<usize> = Vec::new();
        for f in 0..triangles.len() {
            if !band[f] {
                continue;
            }
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
        let at = |v: usize| out_points[v];
        let q0 = (0..q.len())
            .min_by(|&a, &b| {
                (at(q[a]) - at(p[0].0))
                    .length()
                    .total_cmp(&(at(q[b]) - at(p[0].0)).length())
            })
            .unwrap_or(0);
        let (np, nq) = (p.len(), q.len());
        let (mut i, mut j) = (0, 0);
        while i < np || j < nq {
            let pi = p[i % np].0;
            let qj = q[(q0 + j) % nq];
            let advance_p = if i == np {
                false
            } else if j == nq {
                true
            } else {
                let pn = p[(i + 1) % np].0;
                let qn = q[(q0 + j + 1) % nq];
                (at(pn) - at(qj)).length() <= (at(pi) - at(qn)).length()
            };
            if advance_p {
                let pn = p[(i + 1) % np].0;
                out_triangles.push([pi, pn, qj]);
                source.push(Source::Strip(p[i % np].1));
                i += 1;
            } else {
                let qn = q[(q0 + j + 1) % nq];
                out_triangles.push([pi, qn, qj]);
                source.push(Source::Strip(p[i % np].1));
                j += 1;
            }
        }
    }
    Ok(RemainderSurface {
        points: out_points,
        triangles: out_triangles,
        source,
        dropped: dropped.iter().filter(|&&d| d).count(),
    })
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
}
