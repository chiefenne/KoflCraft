//! Voronoi pieces with symbolic vertices.
//!
//! A piece is a scaffold tetrahedron cut by bisector planes of one site.
//! Cutting each tetrahedron on its own and welding the vertices by their
//! coordinates lets two tetrahedra (or two sites) disagree about a thin
//! piece: one side keeps a needle that the other side has collapsed
//! (research/voronoi-symbolic.md). Here a vertex is named by the three
//! planes that meet in it, as in the restricted Voronoi diagrams of Lévy &
//! Liu (2010, Lp centroidal Voronoi tessellation) and Yan, Wang, Lévy & Liu
//! (2013): tetrahedron faces and bisectors, written with the global ids of
//! scaffold vertices and sites. Its coordinates and its side of any further
//! bisector are functions of that name alone, so every tetrahedron and
//! every site that meets the vertex takes the same decision.
//!
//! The side of a vertex is the exact sign of a quotient of two determinants
//! of squared distances from the scaffold vertices of the name to the sites
//! (`predicates::offset_side`, after the side predicates of Lévy's PCK,
//! 2016). It has a third answer, "on the bisector". Conforming wall seeds
//! put bisectors through wall edges on purpose: seeds on the ball of one
//! wall vertex are equidistant from it, and their distance from that vertex
//! is taken to one stand-in seed, so the tie is an exact zero for every
//! site that asks.
use super::convex::{Convex, Cut, Src};
use crate::predicates::{Entry, offset_side};
use mesh_core::geometry::Vec3;
use std::cmp::Ordering;

/// Name of a piece vertex. Arrays are sorted.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(super) enum Key {
    /// A scaffold vertex.
    V(usize),
    /// Scaffold edge and the bisector of two sites.
    E([usize; 2], [u32; 2]),
    /// Scaffold face and the line equidistant from three sites.
    F([usize; 3], [u32; 3]),
    /// The point equidistant from four sites.
    X([u32; 4]),
}

/// One scaffold tetrahedron and the sites.
pub(super) struct Frame<'a> {
    /// Global scaffold vertex ids of the corners.
    pub nodes: [usize; 4],
    pub corners: [Vec3; 4],
    pub sites: &'a [Vec3],
    /// Distance scale for the guard on computed positions.
    pub tol: f64,
    /// The site whose distance from a scaffold vertex stands for that of
    /// the given site: all wall seeds on the ball of the vertex are
    /// equidistant from it by construction and share one stand-in; any
    /// other site stands for itself.
    pub stand_in: &'a (dyn Fn(usize, u32) -> u32 + Sync),
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Class {
    In,
    On,
    Out,
}

fn det3(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    a.cross(b).dot(c)
}

impl Frame<'_> {
    fn corner(&self, node: usize) -> Option<Vec3> {
        self.nodes
            .iter()
            .position(|&n| n == node)
            .map(|i| self.corners[i])
    }

    fn site(&self, s: u32) -> Vec3 {
        self.sites[s as usize]
    }

    /// Signed offset of `p` from the bisector of sites `a < b` along b − a
    /// (unscaled): negative on a's side.
    fn offset(&self, p: Vec3, a: u32, b: u32) -> f64 {
        let (pa, pb) = (self.site(a), self.site(b));
        (p - (pa + pb) * 0.5).dot(pb - pa)
    }

    /// Coordinates of a named vertex, from the name alone.
    fn position(&self, key: Key) -> Option<Vec3> {
        let p = match key {
            Key::V(n) => self.corner(n)?,
            Key::E([n1, n2], [a, b]) => {
                let (v1, v2) = (self.corner(n1)?, self.corner(n2)?);
                let (f1, f2) = (self.offset(v1, a, b), self.offset(v2, a, b));
                let t = f1 / (f1 - f2);
                if !t.is_finite() {
                    return None;
                }
                v1 + (v2 - v1) * t.clamp(0.0, 1.0)
            }
            Key::F([n1, n2, n3], [a, b, c]) => {
                let v1 = self.corner(n1)?;
                let (e1, e2) = (self.corner(n2)? - v1, self.corner(n3)? - v1);
                let (pa, pb, pc) = (self.site(a), self.site(b), self.site(c));
                let (w1, w2) = (pb - pa, pc - pa);
                let (m11, m12, m21, m22) = (e1.dot(w1), e2.dot(w1), e1.dot(w2), e2.dot(w2));
                let r1 = ((pa + pb) * 0.5 - v1).dot(w1);
                let r2 = ((pa + pc) * 0.5 - v1).dot(w2);
                let det = m11 * m22 - m12 * m21;
                let (u, w) = ((r1 * m22 - r2 * m12) / det, (m11 * r2 - m21 * r1) / det);
                if !u.is_finite() || !w.is_finite() {
                    return None;
                }
                v1 + e1 * u + e2 * w
            }
            Key::X([a, b, c, d]) => {
                let pa = self.site(a);
                let (w1, w2, w3) = (self.site(b) - pa, self.site(c) - pa, self.site(d) - pa);
                let (r1, r2, r3) = (w1.dot(w1) * 0.5, w2.dot(w2) * 0.5, w3.dot(w3) * 0.5);
                let det = det3(w1, w2, w3);
                let x = (w2.cross(w3) * r1 + w3.cross(w1) * r2 + w1.cross(w2) * r3) / det;
                if !(x.x.is_finite() && x.y.is_finite() && x.z.is_finite()) {
                    return None;
                }
                pa + x
            }
        };
        Some(p)
    }

    /// Side of vertex `key` of the bisector between site `s`, whose piece
    /// it belongs to, and site `o`. `In` is s's side. The answer is the
    /// exact sign of |q − o|² − |q − s|² at the point q that the name
    /// defines (`predicates::bisector_side`), so it is the same from every
    /// tetrahedron and from every site of the name.
    fn side(&self, key: Key, s: u32, o: u32) -> Class {
        let mut corners = [Vec3::ZERO; 4];
        let mut nodes = [usize::MAX; 4];
        let mut others = [0u32; 3];
        let (count, sites): (usize, &[u32]) = match &key {
            Key::V(n) => {
                nodes[0] = *n;
                (1, &[])
            }
            Key::E(n, k) => {
                nodes[..2].copy_from_slice(n);
                (2, k)
            }
            Key::F(n, k) => {
                nodes[..3].copy_from_slice(n);
                (3, k)
            }
            Key::X(k) => {
                nodes = self.nodes;
                nodes.sort_unstable();
                (4, k)
            }
        };
        if sites.contains(&o) {
            return Class::On;
        }
        // the vertex is equidistant from its sites; s is one of them
        let base = if sites.is_empty() || sites.contains(&s) {
            s
        } else {
            sites[0]
        };
        let mut m = 0;
        for &k in sites {
            if k != base && m < 3 {
                others[m] = k;
                m += 1;
            }
        }
        if m + 1 != count {
            return Class::On;
        }
        for j in 0..count {
            let Some(c) = self.corner(nodes[j]) else {
                return Class::On;
            };
            corners[j] = c;
        }
        // Each squared distance is taken to the site that stands for the
        // given one at that scaffold vertex, so all offsets come from one
        // table and every site of the vertex gets the same answer.
        let offset = |j: usize, site: Option<usize>| -> Entry {
            let k = (self.stand_in)(nodes[j], site.map_or(o, |i| others[i]));
            let b = (self.stand_in)(nodes[j], base);
            if k == b {
                Entry::Zero
            } else {
                Entry::Offset {
                    v: corners[j],
                    s: self.site(b),
                    k: self.site(k),
                }
            }
        };
        match offset_side(count, &offset) {
            Some(Ordering::Greater) => Class::In,
            Some(Ordering::Less) => Class::Out,
            Some(Ordering::Equal) | None => Class::On,
        }
    }

    fn name(&self, planes: [Src; 3], s: u32) -> Option<Key> {
        let mut faces: Vec<u8> = Vec::with_capacity(3);
        let mut sites: Vec<u32> = vec![s];
        for p in planes {
            match p {
                Src::Tet(k) => faces.push(k),
                Src::Bis(o) => sites.push(o),
                Src::Cut => return None,
            }
        }
        sites.sort_unstable();
        sites.dedup();
        faces.sort_unstable();
        faces.dedup();
        if faces.len() + sites.len() != 4 {
            return None;
        }
        let mut nodes: Vec<usize> = (0..4u8)
            .filter(|k| !faces.contains(k))
            .map(|k| self.nodes[k as usize])
            .collect();
        nodes.sort_unstable();
        match faces.len() {
            2 => Some(Key::E([nodes[0], nodes[1]], [sites[0], sites[1]])),
            1 => Some(Key::F(
                [nodes[0], nodes[1], nodes[2]],
                [sites[0], sites[1], sites[2]],
            )),
            0 => Some(Key::X([sites[0], sites[1], sites[2], sites[3]])),
            _ => None,
        }
    }
}

/// A convex piece whose vertices carry their names.
#[derive(Clone, Debug)]
pub(super) struct SymConvex {
    pub keys: Vec<Key>,
    pub pts: Vec<Vec3>,
    /// Outward counter-clockwise rings of local vertex indices.
    pub faces: Vec<(Src, Vec<usize>)>,
    /// An upper bound on the distance of each vertex from the piece's site.
    /// The coordinates of a vertex where three planes come close to sharing
    /// a line are not accurate; the security radius must not trust them.
    pub reach: Vec<f64>,
}

impl SymConvex {
    /// The whole tetrahedron, as the first piece of site `s`.
    pub fn tet(frame: &Frame, s: u32) -> Self {
        let c = Convex::tet(frame.corners);
        let site = frame.site(s);
        Self {
            keys: frame.nodes.iter().map(|&n| Key::V(n)).collect(),
            reach: c
                .pts
                .iter()
                .map(|&p| (p - site).length() * (1.0 + 1e-12))
                .collect(),
            pts: c.pts,
            faces: c.faces,
        }
    }

    pub fn into_convex(self) -> (Vec<Key>, Convex) {
        (
            self.keys,
            Convex {
                pts: self.pts,
                faces: self.faces,
            },
        )
    }

    /// Keep the side of site `s` of its bisector with site `o`.
    pub fn clip(&mut self, frame: &Frame, s: u32, o: u32) -> Result<Cut, &'static str> {
        let class: Vec<Class> = self
            .keys
            .iter()
            .zip(&self.pts)
            .map(|(&k, _)| frame.side(k, s, o))
            .collect();
        if !class.contains(&Class::Out) {
            return Ok(Cut::Inside);
        }
        if !class.contains(&Class::In) {
            return Ok(Cut::Empty);
        }
        let mut keys = Vec::with_capacity(self.keys.len() + 4);
        let mut pts = Vec::with_capacity(self.pts.len() + 4);
        let mut on = Vec::with_capacity(self.pts.len() + 4);
        let mut reach = Vec::with_capacity(self.pts.len() + 4);
        let mut map = vec![usize::MAX; self.pts.len()];
        let own = frame.site(s);
        for (i, &c) in class.iter().enumerate() {
            if c != Class::Out {
                map[i] = pts.len();
                keys.push(self.keys[i]);
                pts.push(self.pts[i]);
                reach.push(self.reach[i]);
                on.push(c == Class::On);
            }
        }
        // the face on the other side of each directed edge
        let mut twin: Vec<((usize, usize), usize)> = Vec::new();
        for (fi, (_, ring)) in self.faces.iter().enumerate() {
            for k in 0..ring.len() {
                twin.push(((ring[(k + 1) % ring.len()], ring[k]), fi));
            }
        }
        let mut crossings: Vec<((usize, usize), usize)> = Vec::new();
        let mut faces = Vec::with_capacity(self.faces.len() + 1);
        for (src, ring) in &self.faces {
            if !ring.iter().any(|&v| class[v] == Class::In) {
                continue;
            }
            let mut out = Vec::with_capacity(ring.len() + 1);
            for k in 0..ring.len() {
                let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
                if class[a] != Class::Out {
                    out.push(map[a]);
                }
                let crosses = matches!(
                    (class[a], class[b]),
                    (Class::In, Class::Out) | (Class::Out, Class::In)
                );
                if !crosses {
                    continue;
                }
                let edge = (a.min(b), a.max(b));
                let id = if let Some(&(_, id)) = crossings.iter().find(|(e, _)| *e == edge) {
                    id
                } else {
                    let Some(&(_, other)) = twin.iter().find(|(e, _)| *e == (a, b)) else {
                        return Err("piece edge with one face");
                    };
                    let key = frame
                        .name([*src, self.faces[other].0, Src::Bis(o)], s)
                        .ok_or("piece edge on fewer than two planes")?;
                    // The point where this edge crosses, and the point of
                    // the name. The name's formula is the same for every
                    // piece but loses its footing when the three planes
                    // come close to sharing a line (four sites on a
                    // circle); the crossing is always on the edge.
                    let (lo, hi) = edge;
                    let (sa, sb) = (s.min(o), s.max(o));
                    let (d1, d2) = (
                        frame.offset(self.pts[lo], sa, sb),
                        frame.offset(self.pts[hi], sa, sb),
                    );
                    let t = d1 / (d1 - d2);
                    let t = if t.is_finite() {
                        t.clamp(0.0, 1.0)
                    } else {
                        0.5
                    };
                    let crossing = self.pts[lo] + (self.pts[hi] - self.pts[lo]) * t;
                    let p = match frame.position(key) {
                        Some(p) if (p - crossing).length() <= frame.tol * 1e4 => p,
                        _ => crossing,
                    };
                    // The true vertex is on the edge, so no farther from
                    // the site than the farther end. Its coordinates are
                    // good to the rounding of t along the edge, plus what
                    // the ends themselves are uncertain by.
                    let ends = self.reach[lo].max(self.reach[hi]);
                    let doubt = |v: usize| (self.reach[v] - (self.pts[v] - own).length()).max(0.0);
                    let size = [self.pts[lo], self.pts[hi], frame.site(sa), frame.site(sb)]
                        .iter()
                        .map(|q| q.max_abs_component())
                        .fold(0.0, f64::max);
                    let slip = 128.0 * f64::EPSILON * size * size / (d1 - d2).abs();
                    let near = (p - own).length()
                        + (p - crossing).length()
                        + slip * (self.pts[hi] - self.pts[lo]).length()
                        + doubt(lo).max(doubt(hi));
                    let id = pts.len();
                    keys.push(key);
                    pts.push(p);
                    reach.push(if near.is_nan() { ends } else { ends.min(near) });
                    on.push(true);
                    crossings.push((edge, id));
                    id
                };
                out.push(id);
            }
            if out.len() >= 3 {
                faces.push((*src, out));
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
        faces.push((Src::Bis(o), cap));
        let mut used = vec![usize::MAX; pts.len()];
        let mut kept_keys = Vec::with_capacity(pts.len());
        let mut kept = Vec::with_capacity(pts.len());
        let mut kept_reach = Vec::with_capacity(pts.len());
        for (_, ring) in &mut faces {
            for v in ring.iter_mut() {
                if used[*v] == usize::MAX {
                    used[*v] = kept.len();
                    kept.push(pts[*v]);
                    kept_reach.push(reach[*v]);
                    kept_keys.push(keys[*v]);
                }
                *v = used[*v];
            }
        }
        self.keys = kept_keys;
        self.pts = kept;
        self.reach = kept_reach;
        self.faces = faces;
        Ok(Cut::Clipped)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(x: f64, y: f64, z: f64) -> Vec3 {
        Vec3::new(x, y, z)
    }

    /// Two tetrahedra on the face 0-1-2, with apexes 3 and 4.
    fn pair() -> ([Vec3; 5], [usize; 4], [usize; 4]) {
        let p = [
            v(0., 0., 0.),
            v(1., 0., 0.),
            v(0., 1., 0.),
            v(0.2, 0.2, 1.),
            v(0.3, 0.1, -1.),
        ];
        (p, [0, 1, 2, 3], [0, 2, 1, 4])
    }

    fn piece(p: &[Vec3; 5], nodes: [usize; 4], sites: &[Vec3], s: u32) -> Option<SymConvex> {
        let frame = Frame {
            nodes,
            corners: nodes.map(|n| p[n]),
            sites,
            tol: 1e-12,
            stand_in: &|_, k| k,
        };
        let mut c = SymConvex::tet(&frame, s);
        for o in 0..sites.len() as u32 {
            if o != s && matches!(c.clip(&frame, s, o).unwrap(), Cut::Empty) {
                return None;
            }
        }
        Some(c)
    }

    /// The vertices of a piece on the scaffold face 0-1-2, by name.
    fn on_shared_face(c: &SymConvex) -> Vec<(Key, [u64; 3])> {
        let mut out: Vec<(Key, [u64; 3])> = c
            .keys
            .iter()
            .zip(&c.pts)
            .filter(|(k, _)| match k {
                Key::V(n) => *n < 3,
                Key::E(n, _) => n.iter().all(|n| *n < 3),
                Key::F(n, _) => n.iter().all(|n| *n < 3),
                Key::X(_) => false,
            })
            .map(|(k, p)| (*k, [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]))
            .collect();
        out.sort();
        out
    }

    /// Keep: the two tetrahedra on a scaffold face give a site's cell the
    /// same vertices on that face, by name and bit for bit, although each
    /// cuts on its own and meets the vertices along different edges. With
    /// positions taken from each tetrahedron's own edge crossing (as before
    /// this module) the bits differ, and a weld by distance had to decide
    /// whether they are one vertex: the jacket failure "piece has no
    /// continuation". No fixture in the fast run has a piece thin enough to
    /// show that.
    #[test]
    fn neighbouring_tetrahedra_share_names_and_coordinates() {
        let (p, up, down) = pair();
        let sites = [
            v(0.21, 0.33, 0.07),
            v(0.63, 0.18, -0.11),
            v(0.17, 0.71, 0.05),
            v(0.45, 0.42, 0.52),
        ];
        let mut seen = 0;
        for s in 0..sites.len() as u32 {
            let (Some(a), Some(b)) = (piece(&p, up, &sites, s), piece(&p, down, &sites, s)) else {
                continue;
            };
            let (fa, fb) = (on_shared_face(&a), on_shared_face(&b));
            assert_eq!(fa, fb, "site {s}");
            seen += fa.iter().filter(|(k, _)| !matches!(k, Key::V(_))).count();
        }
        assert!(seen >= 4, "the bisectors must cross the shared face");
    }

    /// Keep: two sites that share a bisector face name its vertices alike,
    /// each from its own side (the sets of sites are sorted, and the side
    /// test asks about the least site of a name). A name that depended on
    /// the asking site would give the two cells different vertex ids on
    /// their common face.
    #[test]
    fn two_sites_name_their_common_face_alike() {
        let (p, up, _) = pair();
        let sites = [v(0.2, 0.3, 0.3), v(0.6, 0.2, 0.4), v(0.2, 0.6, 0.2)];
        let cap = |s: u32, o: u32| -> Vec<Key> {
            let c = piece(&p, up, &sites, s).unwrap();
            let (_, ring) = c
                .faces
                .iter()
                .find(|(src, _)| *src == Src::Bis(o))
                .expect("the two cells touch");
            let mut keys: Vec<Key> = ring.iter().map(|&i| c.keys[i]).collect();
            keys.sort();
            keys
        };
        assert_eq!(cap(0, 1), cap(1, 0));
        assert_eq!(cap(0, 2), cap(2, 0));
        assert_eq!(cap(1, 2), cap(2, 1));
    }

    /// Keep: a scaffold vertex on the bisector stays a vertex of both
    /// pieces, and no vertex appears beside it. Conforming wall seeds put
    /// bisectors through wall vertices on purpose; a sign taken from
    /// rounding there cuts a needle off one cell.
    #[test]
    fn a_scaffold_vertex_on_the_bisector_belongs_to_both_pieces() {
        let (p, up, _) = pair();
        // mirrored about the plane x = y, which holds vertices 0 and 3
        let sites = [v(0.4, 0.1, 0.3), v(0.1, 0.4, 0.3)];
        for s in 0..2 {
            let c = piece(&p, up, &sites, s).unwrap();
            assert!(c.keys.contains(&Key::V(0)) && c.keys.contains(&Key::V(3)));
            // the cut runs through 0 and 3 and crosses the edge 1-2 only
            let made: Vec<&Key> = c.keys.iter().filter(|k| !matches!(k, Key::V(_))).collect();
            assert_eq!(made, [&Key::E([1, 2], [0, 1])]);
        }
    }
}
