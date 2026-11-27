//! A multi-element airfoil (slat, main element, flap) in a box, extruded
//! along +z (BS-03, MB1g).
//!
//! **Layout.** An H-topology of 2 E + 1 block columns (E elements: a front
//! column, one column per element, a gap column between two elements and
//! a wake column) and four block rows (below, lower middle, upper middle,
//! above). Each element sits in a wall ring (an O-grid one block thick,
//! `RING` of its chord, at most `RING_GAP_SHARE` of a gap) of six blocks
//! whose vertices lie on the wall: the leading edge (the point farthest
//! forward), a lower and an upper corner (the lowest and highest point of
//! the front `CORNER_REACH` of each side), the two trailing-edge corners
//! and the middle of the blunt trailing edge. The ring's outer loop is the
//! wall offset along its normals, so the near-wall cells are orthogonal
//! by construction and the H-blocks meet a smooth loop instead of the
//! wall's corners. Where the next element's nose sits below the previous
//! element's chord (a flap under the main element) the lower trailing-edge
//! vertex slides forward along the lower side by up to `COVE_REACH` of it
//! (the cove), so the gap blocks face the nose across the slot instead of
//! wrapping round the trailing edge. Between two elements the row lines
//! are splines that pass in front of the next element's nose.
//!
//! **Layout optimisation.** The H-blocks off the rings then move toward
//! orthogonal cells (`super::layout`, after Sauer & Morsbach 2023): their
//! corners slide along the box or move freely, their edges bend and their
//! end spacings change, minimising the face non-orthogonality of the TFI
//! grid; the rings stay. On the cases of the tests this takes the largest
//! non-orthogonality from 86-90° to 47-68° with no cell over 70°.
//!
//! The Python prototype of the contour sliding animation (2026-11-24) used
//! one middle row with the corners at the leading edge; that folds as soon
//! as a flap is deflected.
//!
//! **Placement** follows the usual high-lift conventions, in fractions of
//! the main chord: the flap gap is the distance from the main element's
//! trailing edge to the flap, the slat gap the distance from the slat's
//! trailing edge to the main element; the overlap is positive where the
//! next element reaches forward under (over) the previous one. This layout
//! needs the elements to follow each other along x, so the overlaps must
//! be negative; overlapping elements need a slot topology (later).

use mesh_core::block::{
    BlockSpec, BlockTopology, Curve, Distribution, EdgeEnd, EdgeSpec, PatchSpec,
};
use mesh_core::geometry::Vec3;

use super::layout::LayoutOptions;
use super::templates::{Section, Template, counts, positive};

/// The lower and upper corners are the extreme points of this front part
/// of each side (fraction of the element's chord).
const CORNER_REACH: f64 = 0.3;
/// Least angle between the row lines that leave a trailing edge into a
/// gap.
const FAN: f64 = 40.0;
/// Largest angle between a trailing edge's chord and the middle line of
/// the fan that leaves it.
const AIM: f64 = 60.0;
/// Fractions of a gap edge's chord at which it leaves the trailing edge
/// along the fan and reaches the next element along its normal.
const DEPART: f64 = 0.1;
const ARRIVE: f64 = 0.2;
/// Half height of the middle rows where they meet the box front and end,
/// in farfield distances.
const END_SPREAD: f64 = 0.3;
/// How far the column lines move toward even spacing where they meet the
/// box bottom and top.
const BOX_SPREAD: f64 = 0.5;
/// How far (in main chords) a column line runs from a trailing-edge corner
/// in its own direction before it turns toward the box.
const TE_COLUMN_REACH: f64 = 0.05;
/// Default trailing-edge thickness over the element's chord; the sections
/// are closed when read, and the block layout needs a blunt edge.
pub const TRAILING_EDGE_THICKNESS: f64 = 0.005;
/// Thickness of the wall ring around each element, over its chord.
const RING: f64 = 0.05;
/// The ring takes at most this share of the gap on either side of an
/// element, so two rings leave the gap's middle to the gap blocks.
const RING_GAP_SHARE: f64 = 0.3;
/// Largest growth ratio of the cells across the ring.
const RING_GROWTH: f64 = 1.2;
const RING_MIN_CELLS: usize = 4;
/// Longest cove, over the lower side of its element.
const COVE_REACH: f64 = 0.4;

/// The six vertices of an element's wall and of its ring, clockwise: the
/// lower trailing-edge corner, the lower corner, the leading edge, the
/// upper corner, the upper trailing-edge corner, the middle of the blunt
/// trailing edge.
const TEL: usize = 0;
const LOW: usize = 1;
const LE: usize = 2;
const UP: usize = 3;
const TEU: usize = 4;
const TEM: usize = 5;

/// A slat or flap and where it sits relative to the main element.
#[derive(Clone, Debug)]
pub struct Placement {
    pub section: Section,
    /// Chord over the main element's chord.
    pub chord: f64,
    /// Degrees; slat: nose down, flap: trailing edge down.
    pub deflection: f64,
    /// Over the main chord.
    pub gap: f64,
    /// Over the main chord; positive is overlapping, must be negative here.
    pub overlap: f64,
}

/// Parameters of the multi-element template.
#[derive(Clone, Copy, Debug)]
pub struct MultiElement {
    /// Chord of the main element; the leading edge sits at the origin.
    pub chord: f64,
    pub span: f64,
    /// Distance of the box sides (upstream, below, above) from the main
    /// element's leading edge, in main chords.
    pub farfield: f64,
    /// Length of the wake column behind the last trailing edge, in main
    /// chords.
    pub wake: f64,
    /// Cells along the main element's lower side; the slat and flap get
    /// `surface · √(chord ratio)`, at least eight.
    pub surface_cells: usize,
    /// Cells from the elements to the box, above and below.
    pub normal_cells: usize,
    pub wake_cells: usize,
    pub span_cells: usize,
    /// Cells across each middle row: along each nose half and trailing-edge
    /// half.
    pub thickness_cells: usize,
    /// Cells along each gap.
    pub gap_cells: usize,
    /// Cells in front of the first element.
    pub upstream_cells: usize,
    /// First cell height at the walls.
    pub wall_spacing: f64,
    /// Trailing-edge thickness over each element's chord.
    pub trailing_edge_thickness: f64,
}

type P2 = [f64; 2];

fn dist(a: P2, b: P2) -> f64 {
    (b[0] - a[0]).hypot(b[1] - a[1])
}

/// One element placed in the plane: its wall as one ring from the lower
/// trailing-edge corner over the lower side, the leading edge and the
/// upper side to the upper trailing-edge corner.
struct Element {
    name: &'static str,
    ring: Vec<P2>,
    /// Local chord position of each ring point.
    xc: Vec<f64>,
    /// Ring indices of the lower corner, the leading edge and the upper
    /// corner.
    lower: usize,
    le: usize,
    upper: usize,
    area: f64,
    cells: usize,
    chord: f64,
    /// The ring around the wall: its thickness, its cell count across,
    /// the outer curve as a closed clockwise polyline that starts at the
    /// lower trailing-edge corner, and the polyline index of each of the
    /// six vertices.
    delta: f64,
    ring_cells: usize,
    outer: Vec<P2>,
    ov: [usize; 6],
    /// Ring index of the cove vertex on the lower side, where the next
    /// element sits below this one's trailing edge: vertex TEL moves there
    /// and TEM becomes the lower trailing-edge corner, so the gap blocks
    /// meet the cove instead of the trailing edge alone.
    cove: Option<usize>,
}

impl Element {
    /// The section opened to a blunt trailing edge of `thickness` (over the
    /// chord), scaled to `chord`, turned by `angle` (radians, positive
    /// counter-clockwise) about its leading edge, which goes to the origin.
    fn place(
        name: &'static str,
        s: &Section,
        thickness: f64,
        chord: f64,
        angle: f64,
        cells: usize,
    ) -> Result<Self, String> {
        let (c, sn) = (angle.cos(), angle.sin());
        let map = |p: P2, sign: f64| {
            let (x, y) = (p[0], p[1] + sign * 0.5 * thickness * p[0]);
            [chord * (c * x - sn * y), chord * (sn * x + c * y)]
        };
        let mut ring: Vec<P2> = s.lower.iter().rev().map(|&p| map(p, -1.0)).collect();
        let mut xc: Vec<f64> = s.lower.iter().rev().map(|p| p[0]).collect();
        ring.extend(s.upper.iter().skip(1).map(|&p| map(p, 1.0)));
        xc.extend(s.upper.iter().skip(1).map(|p| p[0]));
        let mut e = Self {
            name,
            ring,
            xc,
            lower: 0,
            le: 0,
            upper: 0,
            area: chord * chord * (s.area + 0.5 * thickness),
            cells,
            chord,
            delta: 0.0,
            ring_cells: 0,
            outer: Vec::new(),
            ov: [0; 6],
            cove: None,
        };
        e.find_corners(None)?;
        Ok(e)
    }

    /// The corners on the nose: the middle one is the point nearest to
    /// `facing` (the previous element's trailing edge), or the point
    /// farthest forward; the lower and upper corners are where the wall
    /// normal has turned at least 45° from the middle one's, down and up,
    /// and points within 22.5° of the vertical where it can.
    fn find_corners(&mut self, facing: Option<P2>) -> Result<(), String> {
        let n = self.ring.len();
        let front: Vec<usize> = (2..n - 2).filter(|&i| self.xc[i] <= CORNER_REACH).collect();
        let key = |i: usize| match facing {
            Some(f) => dist(self.ring[i], f),
            None => self.ring[i][0],
        };
        let le = front
            .iter()
            .copied()
            .min_by(|&a, &b| key(a).total_cmp(&key(b)))
            .ok_or("the section has no front")?;
        let angle = |i: usize| {
            let m = self.normal(i);
            m[1].atan2(m[0])
        };
        let off = |i: usize, target: f64| {
            let d = (angle(i) - target).rem_euclid(std::f64::consts::TAU);
            d.min(std::f64::consts::TAU - d)
        };
        let pick = |range: &mut dyn Iterator<Item = usize>, target: f64| {
            range.min_by(|&a, &b| off(a, target).total_cmp(&off(b, target)))
        };
        // normal directions in [0, 2π): lower side 3π/2, nose π, upper π/2;
        // the column lines leave the corners vertically, so their normals
        // point at most 22.5° off the vertical, and 45° away from the middle
        let mid = angle(le).rem_euclid(std::f64::consts::TAU);
        let eighth = std::f64::consts::FRAC_PI_4;
        let lower = pick(
            &mut front.iter().copied().filter(|&i| i + 2 <= le),
            (mid + eighth).max(1.375 * std::f64::consts::PI),
        );
        let upper = pick(
            &mut (le + 2..n - 2).filter(|&i| self.xc[i] <= 2.0 * CORNER_REACH),
            (mid - eighth).min(0.625 * std::f64::consts::PI),
        );
        match (lower, upper) {
            (Some(lo), Some(up)) => {
                (self.lower, self.le, self.upper) = (lo, le, up);
                Ok(())
            }
            _ => Err(format!(
                "the {} section's nose cannot be split into a lower and an upper half \
                 (is it turned too far?)",
                self.name
            )),
        }
    }

    /// The outward unit normal at ring point `i`.
    fn normal(&self, i: usize) -> P2 {
        let (a, b) = (self.ring[i - 1], self.ring[i + 1]);
        let (tx, ty) = (b[0] - a[0], b[1] - a[1]);
        let l = tx.hypot(ty);
        [-ty / l, tx / l]
    }

    fn shifted(mut self, d: P2) -> Self {
        for p in self.ring.iter_mut() {
            p[0] += d[0];
            p[1] += d[1];
        }
        self
    }

    fn te_lower(&self) -> P2 {
        self.ring[0]
    }
    fn te_upper(&self) -> P2 {
        self.ring[self.ring.len() - 1]
    }
    fn te_mid(&self) -> P2 {
        let (a, b) = (self.te_upper(), self.te_lower());
        [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0]
    }

    /// The wall point of vertex `v`.
    fn wall_point(&self, v: usize) -> P2 {
        match (v, self.cove) {
            (TEL, Some(c)) => self.ring[c],
            (TEL, None) => self.te_lower(),
            (LOW, _) => self.ring[self.lower],
            (LE, _) => self.ring[self.le],
            (UP, _) => self.ring[self.upper],
            (TEU, _) => self.te_upper(),
            (_, Some(_)) => self.te_lower(),
            (_, None) => self.te_mid(),
        }
    }

    /// The ring point of vertex `v`.
    fn outer_point(&self, v: usize) -> P2 {
        self.outer[self.ov[v]]
    }

    /// Lay the ring of thickness `delta` around the wall: the offset of
    /// the sides along their normals, with circular arcs about the two
    /// trailing-edge corners, which the ring leaves along the bisector.
    fn set_ring(&mut self, delta: f64, cells: usize) {
        let n = self.ring.len();
        let unit = |v: P2| {
            let l = v[0].hypot(v[1]);
            [v[0] / l, v[1] / l]
        };
        let rot = |t: P2| unit([-t[1], t[0]]);
        let (c0, c1) = (self.ring[0], self.ring[n - 1]);
        let side0 = rot([self.ring[1][0] - c0[0], self.ring[1][1] - c0[1]]);
        let side1 = rot([c1[0] - self.ring[n - 2][0], c1[1] - self.ring[n - 2][1]]);
        let te = rot([c0[0] - c1[0], c0[1] - c1[1]]);
        let bis0 = unit([side0[0] + te[0], side0[1] + te[1]]);
        let bis1 = unit([side1[0] + te[0], side1[1] + te[1]]);
        let at = |c: P2, d: P2| [c[0] + delta * d[0], c[1] + delta * d[1]];
        // points on the arc about c from direction a to b, without a,
        // with b
        let arc = |c: P2, a: P2, b: P2| -> Vec<P2> {
            let phi = (a[0] * b[1] - a[1] * b[0]).atan2(a[0] * b[0] + a[1] * b[1]);
            let m = ((phi.abs() / (std::f64::consts::PI / 12.0)).ceil() as usize).max(1);
            (1..=m)
                .map(|k| {
                    let t = phi * k as f64 / m as f64;
                    let (cs, sn) = (t.cos(), t.sin());
                    at(c, [cs * a[0] - sn * a[1], sn * a[0] + cs * a[1]])
                })
                .collect()
        };
        let mut o = vec![at(c0, bis0)];
        o.extend(arc(c0, bis0, side0));
        let mut ov = [0; 6];
        for i in 1..n - 1 {
            if Some(i) == self.cove {
                ov[TEL] = o.len();
            } else if i == self.lower {
                ov[LOW] = o.len();
            } else if i == self.le {
                ov[LE] = o.len();
            } else if i == self.upper {
                ov[UP] = o.len();
            }
            o.push(at(self.ring[i], self.normal(i)));
        }
        o.push(at(c1, side1));
        o.extend(arc(c1, side1, bis1));
        ov[TEU] = o.len() - 1;
        o.extend(arc(c1, bis1, te));
        if self.cove.is_none() {
            ov[TEM] = o.len();
        }
        o.push(at(self.te_mid(), te));
        o.push(at(c0, te));
        let mut back = arc(c0, te, bis0);
        back.pop();
        o.extend(back);
        self.outer = o;
        self.ov = ov;
        self.delta = delta;
        self.ring_cells = cells;
    }

    /// The vertex after `v`, clockwise.
    fn next_vertex(v: usize) -> usize {
        match v {
            TEL => LOW,
            LOW => LE,
            LE => UP,
            UP => TEU,
            TEU => TEM,
            _ => TEL,
        }
    }

    /// Ring points strictly between the neighbouring vertices `a` and `b`
    /// (either order), from `a` to `b`.
    fn outer_between(&self, a: usize, b: usize) -> Vec<P2> {
        let m = self.outer.len();
        let cw = |x: usize, y: usize| {
            let (px, py) = (self.ov[x], self.ov[y]);
            let end = if py == 0 { m } else { py };
            self.outer[px + 1..end].to_vec()
        };
        if Self::next_vertex(a) == b {
            cw(a, b)
        } else {
            let mut v = cw(b, a);
            v.reverse();
            v
        }
    }

    /// Length of the ring curve between two neighbouring vertices.
    fn outer_length(&self, a: usize, b: usize) -> f64 {
        let mut pts = vec![self.outer_point(a)];
        pts.extend(self.outer_between(a, b));
        pts.push(self.outer_point(b));
        pts.windows(2).map(|w| dist(w[0], w[1])).sum()
    }

    /// How much longer the ring is than the wall about ring point `i`.
    fn stretch(&self, i: usize) -> f64 {
        let off = |j: usize| {
            let m = self.normal(j);
            [
                self.ring[j][0] + self.delta * m[0],
                self.ring[j][1] + self.delta * m[1],
            ]
        };
        dist(off(i - 1), off(i + 1)) / dist(self.ring[i - 1], self.ring[i + 1])
    }

    /// Distance from `p` to the closed wall.
    fn distance(&self, p: P2) -> f64 {
        let n = self.ring.len();
        (0..n)
            .map(|i| {
                let (a, b) = (self.ring[i], self.ring[(i + 1) % n]);
                let d = [b[0] - a[0], b[1] - a[1]];
                let dd = d[0] * d[0] + d[1] * d[1];
                let t = if dd > 0.0 {
                    (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / dd).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                (p[0] - a[0] - t * d[0]).hypot(p[1] - a[1] - t * d[1])
            })
            .fold(f64::INFINITY, f64::min)
    }

    /// Polyline length between two ring indices.
    fn length(&self, a: usize, b: usize) -> f64 {
        let (lo, hi) = (a.min(b), a.max(b));
        (lo..hi).map(|i| dist(self.ring[i], self.ring[i + 1])).sum()
    }
}

/// The parameter t from `t0` toward `t1` at which `distance(t)` first falls
/// to `target`, by a scan and bisection.
fn approach(distance: impl Fn(f64) -> f64, target: f64, t0: f64, t1: f64) -> Option<f64> {
    let steps = 4000;
    if distance(t0) <= target {
        return None;
    }
    let mut prev = t0;
    for k in 1..=steps {
        let t = t0 + (t1 - t0) * k as f64 / steps as f64;
        if distance(t) <= target {
            let (mut far, mut near) = (prev, t);
            for _ in 0..60 {
                let mid = 0.5 * (far + near);
                if distance(mid) <= target {
                    near = mid;
                } else {
                    far = mid;
                }
            }
            return Some(far);
        }
        prev = t;
    }
    None
}

fn v(p: P2, z: f64) -> Vec3 {
    Vec3::new(p[0], p[1], z)
}

fn check_deflection(what: &str, d: f64) -> Result<(), String> {
    if d.is_finite() && d.abs() <= 80.0 {
        Ok(())
    } else {
        Err(format!("the {what} deflection must lie within ±80°"))
    }
}

/// The elements placed one behind the other.
fn place_all(
    main: &Section,
    slat: Option<&Placement>,
    flap: Option<&Placement>,
    p: &MultiElement,
    c: f64,
    far: f64,
) -> Result<Vec<Element>, String> {
    let tte = p.trailing_edge_thickness;
    let cells_for = |ratio: f64| ((p.surface_cells as f64 * ratio.sqrt()).round() as usize).max(8);
    let main_el = Element::place("main", main, tte, c, 0.0, p.surface_cells)?;
    let mut elements = Vec::new();
    if let Some(s) = slat {
        let cs = positive(s.chord, "slat chord")?;
        positive(s.gap, "slat gap")?;
        check_deflection("slat", s.deflection)?;
        if s.overlap >= 0.0 {
            return Err(format!(
                "slat overlap {}: this block layout needs the slat's trailing edge ahead of the \
                 main element's leading edge (a negative overlap); an overlapping slat needs \
                 the slot topology, which is not built yet",
                s.overlap
            ));
        }
        // nose down: the trailing edge turns up about the leading edge
        let el = Element::place(
            "slat",
            &s.section,
            tte,
            cs * c,
            s.deflection.to_radians(),
            cells_for(cs),
        )?;
        let te = el.te_mid();
        let x = s.overlap * c;
        let y = approach(|y| main_el.distance([x, y]), s.gap * c, 0.5 * far, -0.5 * c).ok_or_else(
            || {
                format!(
                    "no slat position has gap {} with overlap {}; the gap must be larger than \
                     the overlap's size",
                    s.gap, s.overlap
                )
            },
        )?;
        elements.push(el.shifted([x - te[0], y - te[1]]));
    }
    elements.push(main_el);
    if let Some(f) = flap {
        let cf = positive(f.chord, "flap chord")?;
        positive(f.gap, "flap gap")?;
        check_deflection("flap", f.deflection)?;
        if f.overlap >= 0.0 {
            return Err(format!(
                "flap overlap {}: this block layout needs the flap's leading edge behind the \
                 main element's trailing edge (a negative overlap); an overlapping flap needs \
                 the slot topology, which is not built yet",
                f.overlap
            ));
        }
        let el = Element::place(
            "flap",
            &f.section,
            tte,
            cf * c,
            -f.deflection.to_radians(),
            cells_for(cf),
        )?;
        let main_te = elements[elements.len() - 1].te_mid();
        let le = el.ring[el.le];
        let x = c - f.overlap * c - le[0];
        let y = approach(
            |y| el.distance([main_te[0] - x, main_te[1] - y]),
            f.gap * c,
            -0.5 * far,
            0.5 * c,
        )
        .ok_or_else(|| format!("no flap position has gap {}", f.gap))?;
        elements.push(el.shifted([x, y]));
    }
    for k in 1..elements.len() {
        let te = elements[k - 1].te_mid();
        elements[k].find_corners(Some(te))?;
    }
    Ok(elements)
}

/// The multi-element airfoil template: `main` with an optional slat in
/// front and an optional flap behind. Patches: `slat_wall`, `main_wall`,
/// `flap_wall`, `inlet` (front, top and bottom of the box), `outlet` (the
/// box end behind the wake), `symmetry_low` and `symmetry_high`.
pub fn multi_element_airfoil(
    main: &Section,
    slat: Option<&Placement>,
    flap: Option<&Placement>,
    p: &MultiElement,
) -> Result<Template, String> {
    let c = positive(p.chord, "chord")?;
    let span = positive(p.span, "span")?;
    let far = positive(p.farfield, "farfield distance")? * c;
    let wake = positive(p.wake, "wake length")? * c;
    if far < 2.0 * c {
        return Err("the farfield distance must be at least two chords".into());
    }
    counts(&[
        p.surface_cells,
        p.normal_cells,
        p.wake_cells,
        p.span_cells,
        p.thickness_cells,
        p.gap_cells,
        p.upstream_cells,
    ])?;
    if p.surface_cells < 8 {
        return Err("the main element needs at least eight surface cells".into());
    }
    let h = positive(p.wall_spacing, "wall spacing")?;
    let tte = positive(p.trailing_edge_thickness, "trailing-edge thickness")?;
    if tte > 0.05 {
        return Err("the trailing-edge thickness must be at most 5 % of the chord".into());
    }
    let mut elements = place_all(main, slat, flap, p, c, far)?;
    let ne = elements.len();
    // a cove where the next element's leading edge lies below this one's
    // chord: as long as the next element's lower nose half
    for k in 0..ne.saturating_sub(1) {
        let (prev, next) = (&elements[k], &elements[k + 1]);
        let (le, te, q) = (prev.ring[prev.le], prev.te_mid(), next.ring[next.le]);
        let below = (te[0] - le[0]) * (q[1] - le[1]) - (te[1] - le[1]) * (q[0] - le[0]) < 0.0;
        if !below {
            continue;
        }
        let want = next
            .length(next.lower, next.le)
            .min(COVE_REACH * prev.length(prev.lower, 0));
        let (mut c, mut acc) = (1, 0.0);
        while c + 2 < prev.lower && acc + dist(prev.ring[c - 1], prev.ring[c]) < want {
            acc += dist(prev.ring[c - 1], prev.ring[c]);
            c += 1;
        }
        if c >= 2 {
            elements[k].cove = Some(c);
        }
    }
    // the ring of each element: a share of its chord, at most a share of
    // the gap on either side, with geometric cells from the wall spacing
    let gaps: Vec<f64> = (0..ne.saturating_sub(1))
        .map(|k| {
            let (a, b) = (&elements[k], &elements[k + 1]);
            let cove = a.ring[..=a.cove.unwrap_or(0)].iter();
            cove.chain(&[a.te_mid(), a.te_upper()])
                .map(|&q| b.distance(q))
                .fold(f64::INFINITY, f64::min)
        })
        .collect();
    for k in 0..ne {
        // thick enough for a few cells from a coarse wall spacing
        let mut delta = (RING * elements[k].chord).max(2.0 * h * RING_MIN_CELLS as f64);
        if k > 0 {
            delta = delta.min(RING_GAP_SHARE * gaps[k - 1]);
        }
        if k + 1 < ne {
            delta = delta.min(RING_GAP_SHARE * gaps[k]);
        }
        if h * RING_MIN_CELLS as f64 >= delta {
            return Err(format!(
                "the wall spacing {h} is too large for the {} ring ({delta:.3e} thick); use a \
                 smaller wall spacing",
                elements[k].name
            ));
        }
        let cells = (((1.0 + delta * (RING_GROWTH - 1.0) / h).ln() / RING_GROWTH.ln()).ceil()
            as usize)
            .max(RING_MIN_CELLS);
        elements[k].set_ring(delta, cells);
    }
    let elements = elements;
    // column lines 0..nl (element k: lines 2k+1 and 2k+2), row lines 0..=4
    let nl = 2 * ne + 2;
    let (x_min, y_min, y_max) = (-far, -far, far);
    let last = &elements[ne - 1];
    let x_max = last.te_upper()[0].max(last.te_lower()[0]) + wake;
    for e in &elements {
        for q in e.ring.iter().chain(&e.outer) {
            if q[0] <= x_min + 0.5 * c || q[1] <= y_min + 0.5 * c || q[1] >= y_max - 0.5 * c {
                return Err("the farfield distance is too short for the elements".into());
            }
        }
    }
    // the element vertex at a column line and row
    let vertex_of = |line: usize, row: usize| -> usize {
        match (line % 2 == 1, row) {
            (true, 1) => LOW,
            (true, 2) => LE,
            (true, _) => UP,
            (false, 1) => TEL,
            (false, 2) => TEM,
            (false, _) => TEU,
        }
    };
    // the background blocks meet the elements' rings
    let hole = |line: usize, row: usize| -> P2 {
        elements[(line - 1) / 2].outer_point(vertex_of(line, row))
    };
    let corner = |line: usize, row: usize| -> P2 {
        match row {
            // at the box ends the middle rows open up like the front of a
            // C-grid, so their blocks are wedges, not long thin strips
            1..=3 if line == 0 => [x_min, hole(1, 2)[1] + (row as f64 - 2.0) * END_SPREAD * far],
            1..=3 if line == nl - 1 => [
                x_max,
                hole(nl - 2, 2)[1] + (row as f64 - 2.0) * END_SPREAD * far,
            ],
            1..=3 => hole(line, row),
            _ => {
                let x = if line == 0 {
                    x_min
                } else if line == nl - 1 {
                    x_max
                } else {
                    // halfway to even spacing between the first and the
                    // last element line, so the narrow gap columns widen
                    // toward the box
                    let r = if row == 0 { 1 } else { 3 };
                    let (a, b) = (hole(1, r)[0], hole(nl - 2, r)[0]);
                    let even = a + (b - a) * (line - 1) as f64 / (nl - 3) as f64;
                    hole(line, r)[0] + BOX_SPREAD * (even - hole(line, r)[0])
                };
                [x, if row == 0 { y_min } else { y_max }]
            }
        }
    };
    // the elements must follow each other along x on every corner row
    for row in 1..=3 {
        for line in 0..nl - 1 {
            if corner(line + 1, row)[0] <= corner(line, row)[0] {
                return Err(format!(
                    "the elements overlap along the chord (row {row}, column line {}); move \
                     the slat forward or the flap back (a more negative overlap)",
                    line + 1
                ));
            }
        }
    }
    let id = |line: usize, row: usize, layer: usize| line + nl * row + 5 * nl * layer;
    // the wall vertices of element k, inside its ring
    let wid = |k: usize, v: usize, layer: usize| 10 * nl + 6 * k + v + 6 * ne * layer;
    // the ring vertex of element k's vertex v
    let oid = |k: usize, v: usize, layer: usize| -> usize {
        let (line, row) = match v {
            LOW => (2 * k + 1, 1),
            LE => (2 * k + 1, 2),
            UP => (2 * k + 1, 3),
            TEL => (2 * k + 2, 1),
            TEM => (2 * k + 2, 2),
            _ => (2 * k + 2, 3),
        };
        id(line, row, layer)
    };
    let mut vertices = Vec::with_capacity(10 * nl + 12 * ne);
    for z in [0.0, span] {
        for row in 0..5 {
            for line in 0..nl {
                vertices.push(v(corner(line, row), z));
            }
        }
    }
    for z in [0.0, span] {
        for e in &elements {
            for vtx in 0..6 {
                vertices.push(v(e.wall_point(vtx), z));
            }
        }
    }
    let mut col_cells = vec![p.upstream_cells];
    for (k, e) in elements.iter().enumerate() {
        col_cells.push(e.cells);
        col_cells.push(if k + 1 < ne {
            p.gap_cells
        } else {
            p.wake_cells
        });
    }
    let row_cells = [
        p.normal_cells,
        p.thickness_cells,
        p.thickness_cells,
        p.normal_cells,
    ];
    // cell lengths along each element at the leading edge, the corners and
    // the trailing edge, on the wall and on its ring
    #[derive(Clone, Copy)]
    struct Spacing {
        le: f64,
        corner: [f64; 2],
        te: f64,
        /// At the cove vertex, and at the trailing-edge corner along the
        /// cove.
        cove: f64,
        cove_te: f64,
    }
    let nt = p.thickness_cells as f64;
    let wall_spacing: Vec<Spacing> = elements
        .iter()
        .map(|e| {
            let mean = e.length(e.lower, 0) / e.cells as f64;
            let (nose_lo, nose_up) = (e.length(e.lower, e.le), e.length(e.le, e.upper));
            let te_half = 0.5 * dist(e.te_upper(), e.te_lower());
            let le = (0.2 * mean).min(0.3 * nose_lo.min(nose_up) / nt);
            let te = (0.2 * mean).max(te_half / nt);
            Spacing {
                le,
                corner: [
                    (nose_lo / nt).min(mean).max(1.5 * le),
                    (nose_up / nt).min(mean).max(1.5 * le),
                ],
                te,
                cove: te,
                cove_te: te,
            }
        })
        .collect();
    // the cove takes the law of the nose half it faces, scaled to its
    // length, so the gap block's sides and the ring's lines across agree
    let mut wall_spacing = wall_spacing;
    for k in 0..ne.saturating_sub(1) {
        if let Some(c) = elements[k].cove {
            let next = &elements[k + 1];
            let ratio = elements[k].length(0, c) / next.length(next.lower, next.le);
            let n = wall_spacing[k + 1];
            wall_spacing[k].cove = n.corner[0] * ratio;
            wall_spacing[k].cove_te = n.le * ratio;
        }
    }
    let wall_spacing = wall_spacing;
    // on the ring the cells are longer where the wall is convex
    let spacing: Vec<Spacing> = elements
        .iter()
        .zip(&wall_spacing)
        .map(|(e, s)| Spacing {
            le: s.le * e.stretch(e.le),
            corner: [
                s.corner[0] * e.stretch(e.lower),
                s.corner[1] * e.stretch(e.upper),
            ],
            te: s.te,
            cove: e.cove.map_or(s.cove, |c| s.cove * e.stretch(c)),
            cove_te: s.cove_te,
        })
        .collect();
    // the last cell across each ring
    let ring_out: Vec<f64> = elements
        .iter()
        .map(|e| last_cell(h, e.ring_cells, e.delta))
        .collect();
    // the directions in which the three row lines leave each trailing edge
    // into the gap behind it: the middle one toward the next element's
    // middle corner (at most AIM off the chord), the others toward their
    // corners but at least FAN beside the middle one; behind the last
    // element, toward the box end
    let fans: Vec<[f64; 3]> = (0..ne)
        .map(|k| {
            let prev = &elements[k];
            let te = [
                prev.outer_point(TEL),
                prev.outer_point(TEM),
                prev.outer_point(TEU),
            ];
            let toward = |from: P2, to: P2| (to[1] - from[1]).atan2(to[0] - from[0]);
            if k + 1 == ne {
                return [1, 2, 3].map(|row| toward(te[row - 1], corner(nl - 1, row)));
            }
            let next = &elements[k + 1];
            let base = toward(prev.ring[prev.le], prev.te_mid());
            let off = |t: f64| {
                (t - base + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
                    - std::f64::consts::PI
            };
            let (aim, fan, limit) = (AIM.to_radians(), FAN.to_radians(), 80f64.to_radians());
            let mid = off(toward(te[1], next.outer_point(LE))).clamp(-aim, aim);
            let lower = off(toward(te[0], next.outer_point(LOW)));
            let upper = off(toward(te[2], next.outer_point(UP)));
            [
                base + lower.min(mid - fan).max(-limit),
                base + mid,
                base + upper.max(mid + fan).min(limit),
            ]
        })
        .collect();
    // the column line leaving a trailing-edge corner of the ring halves the
    // angle between the ring and the gap's lower (upper) row line
    let te_column = |k: usize, upper: bool| -> f64 {
        let e = &elements[k];
        if let (false, Some(c)) = (upper, e.cove) {
            let m = e.normal(c);
            return m[1].atan2(m[0]);
        }
        let (corner, along) = if upper {
            (e.outer_point(TEU), e.outer[e.ov[TEU] - 1])
        } else {
            (e.outer_point(TEL), e.outer[e.ov[TEL] + 1])
        };
        let wall = (along[1] - corner[1]).atan2(along[0] - corner[0]);
        let tau = std::f64::consts::TAU;
        if upper {
            let row = fans[k][2];
            row + 0.5 * (wall - row).rem_euclid(tau)
        } else {
            let row = fans[k][0];
            row - 0.5 * (row - wall).rem_euclid(tau)
        }
    };
    // the ring curve of element k's nose half below (row 1) or above
    // (row 2) the leading edge, its law and its length
    let nose = |k: usize, row: usize| {
        let (e, s) = (&elements[k], &spacing[k]);
        if row == 1 {
            (
                e.outer_between(LOW, LE),
                Distribution::Tanh {
                    first: s.corner[0],
                    last: s.le,
                },
                e.outer_length(LOW, LE),
            )
        } else {
            (
                e.outer_between(LE, UP),
                Distribution::Tanh {
                    first: s.le,
                    last: s.corner[1],
                },
                e.outer_length(LE, UP),
            )
        }
    };
    let scale = |law: Distribution, ratio: f64| match law {
        Distribution::Geometric { spacing, at } => Distribution::Geometric {
            spacing: spacing * ratio,
            at,
        },
        Distribution::Tanh { first, last } => Distribution::Tanh {
            first: first * ratio,
            last: last * ratio,
        },
        Distribution::Uniform => Distribution::Uniform,
    };
    // the spacing at the far end of the column lines from the box to the
    // middle rows: the last ring cell at an element, the spacing of the
    // box side next to it at the box front and end
    let to_box = |line: usize, row: usize| -> f64 {
        if line == 0 {
            let (_, law, l) = nose(0, if row == 1 { 1 } else { 2 });
            let ratio = dist(corner(0, 2), corner(0, row)) / l;
            match law {
                Distribution::Tanh { first, last } => {
                    if row == 1 {
                        first * ratio
                    } else {
                        last * ratio
                    }
                }
                _ => ratio * l / nt,
            }
        } else if line == nl - 1 {
            dist(corner(line, 2), corner(line, row)) / nt
        } else {
            ring_out[(line - 1) / 2]
        }
    };
    let mut edges = Vec::new();
    for layer in 0..2 {
        let z = if layer == 0 { 0.0 } else { span };
        let spline = |pts: Vec<P2>| Curve::Spline {
            through: pts.into_iter().map(|q| v(q, z)).collect(),
        };
        let curve_length = |a: P2, curve: &Curve, b: P2| -> f64 {
            match curve {
                Curve::Spline { through } => {
                    let mut l = 0.0;
                    let mut prev = v(a, z);
                    for &q in through.iter().chain(std::iter::once(&v(b, z))) {
                        l += (q - prev).length();
                        prev = q;
                    }
                    l
                }
                _ => dist(a, b),
            }
        };
        // along the row lines 1–3, column by column, and the box edges
        for col in 0..nl - 1 {
            for row in 1..=3 {
                let (a, b) = (corner(col, row), corner(col + 1, row));
                // the row lines leave the rings across them: their first
                // cells continue the rings' last ones
                let (curve, law) = if col == 0 {
                    (
                        Curve::Line,
                        Distribution::Geometric {
                            spacing: ring_out[0],
                            at: EdgeEnd::End,
                        },
                    )
                } else if col == nl - 2 {
                    (
                        Curve::Line,
                        Distribution::Geometric {
                            spacing: ring_out[ne - 1],
                            at: EdgeEnd::Start,
                        },
                    )
                } else if col % 2 == 1 {
                    if row == 2 {
                        continue; // inside the element
                    }
                    let k = (col - 1) / 2;
                    let (e, s) = (&elements[k], &wall_spacing[k]);
                    let (from, to, wall_from, wall_to, end) = if row == 1 {
                        (LOW, TEL, e.lower, e.cove.unwrap_or(0), s.cove)
                    } else {
                        (UP, TEU, e.upper, e.ring.len() - 1, s.te)
                    };
                    // the wall's law on the ring's length: the same
                    // fractions, so the lines across the ring stay normal
                    let ratio = e.outer_length(from, to) / e.length(wall_from, wall_to);
                    (
                        spline(e.outer_between(from, to)),
                        scale(
                            Distribution::Tanh {
                                first: s.corner[if row == 1 { 0 } else { 1 }],
                                last: end,
                            },
                            ratio,
                        ),
                    )
                } else {
                    // a gap: pass in front of the next element's nose
                    let k = (col - 2) / 2;
                    let next = &elements[k + 1];
                    // leave the trailing edge in the fan and reach the next
                    // element's ring across it
                    let (vtx, i) = match row {
                        1 => (LOW, next.lower),
                        2 => (LE, next.le),
                        _ => (UP, next.upper),
                    };
                    let nrm = next.normal(i);
                    let q = next.outer_point(vtx);
                    let d = dist(a, b);
                    // halfway between the wall normal and the way back to
                    // the trailing edge, so the line stays on the gap's side
                    // of the corner's column line
                    let nrm = {
                        let w = [nrm[0] + (a[0] - q[0]) / d, nrm[1] + (a[1] - q[1]) / d];
                        let l = w[0].hypot(w[1]);
                        [w[0] / l, w[1] / l]
                    };
                    let turn = fans[k][row - 1];
                    let mut out = [turn.cos(), turn.sin()];
                    if let (1, Some(c)) = (row, elements[k].cove) {
                        // from the cove: between its normal and the target
                        let m = elements[k].normal(c);
                        let w = [m[0] + (q[0] - a[0]) / d, m[1] + (q[1] - a[1]) / d];
                        let l = w[0].hypot(w[1]);
                        out = [w[0] / l, w[1] / l];
                    }
                    let (d1, d2) = (DEPART * d, ARRIVE * d);
                    let curve = spline(vec![
                        [a[0] + d1 * out[0], a[1] + d1 * out[1]],
                        [q[0] + d2 * nrm[0], q[1] + d2 * nrm[1]],
                    ]);
                    (
                        curve,
                        Distribution::Tanh {
                            first: ring_out[k],
                            last: ring_out[k + 1],
                        },
                    )
                };
                let length = curve_length(a, &curve, b);
                edges.push(EdgeSpec {
                    vertices: [id(col, row, layer), id(col + 1, row, layer)],
                    curve,
                    distribution: law,
                });
                if row == 2 {
                    continue;
                }
                // the box edge on this side, with the spacings scaled
                let box_row = if row == 1 { 0 } else { 4 };
                let ratio = dist(corner(col, box_row), corner(col + 1, box_row)) / length;
                edges.push(EdgeSpec {
                    vertices: [id(col, box_row, layer), id(col + 1, box_row, layer)],
                    curve: Curve::Line,
                    distribution: scale(law, ratio),
                });
            }
        }
        // along the column lines
        for line in 0..nl {
            // from a trailing-edge corner the line leaves between the ring
            // and the gap's row line, then turns toward the box
            let leave = |row: usize| -> Curve {
                if line == 0 || line == nl - 1 || line % 2 == 1 {
                    return Curve::Line;
                }
                let k = line / 2 - 1;
                let t = te_column(k, row == 3);
                let p = corner(line, row);
                let q = corner(line, if row == 3 { 4 } else { 0 });
                // a quadratic Bézier from the corner, sampled densely near it
                let ctrl = [
                    p[0] + TE_COLUMN_REACH * c * t.cos(),
                    p[1] + TE_COLUMN_REACH * c * t.sin(),
                ];
                let bezier = |u: f64| {
                    let (a, b, d) = ((1.0 - u) * (1.0 - u), 2.0 * u * (1.0 - u), u * u);
                    [
                        a * p[0] + b * ctrl[0] + d * q[0],
                        a * p[1] + b * ctrl[1] + d * q[1],
                    ]
                };
                spline(
                    [0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7]
                        .iter()
                        .map(|&u| bezier(u))
                        .collect(),
                )
            };
            let down = match leave(1) {
                Curve::Spline { mut through } => {
                    through.reverse();
                    Curve::Spline { through }
                }
                other => other,
            };
            edges.push(EdgeSpec {
                vertices: [id(line, 0, layer), id(line, 1, layer)],
                curve: down,
                distribution: Distribution::Geometric {
                    spacing: to_box(line, 1),
                    at: EdgeEnd::End,
                },
            });
            edges.push(EdgeSpec {
                vertices: [id(line, 3, layer), id(line, 4, layer)],
                curve: leave(3),
                distribution: Distribution::Geometric {
                    spacing: to_box(line, 3),
                    at: EdgeEnd::Start,
                },
            });
            for row in [1, 2] {
                // the trailing-edge half of element k's ring below (row 1)
                // or above (row 2) the middle
                let te_half = |k: usize| {
                    let e = &elements[k];
                    let (a, b) = if row == 1 { (TEL, TEM) } else { (TEM, TEU) };
                    (e.outer_between(a, b), e.outer_length(a, b))
                };
                let (curve, law) = if line == nl - 1 {
                    (Curve::Line, Distribution::Uniform)
                } else if line == nl - 2 {
                    (spline(te_half(ne - 1).0), Distribution::Uniform)
                } else if line % 2 == 1 {
                    let (pts, law, _) = nose((line - 1) / 2, row);
                    (spline(pts), law)
                } else {
                    // the box front or a trailing edge: the law of the nose
                    // across the block, scaled to this edge's length, so that
                    // the interpolation inside the block stays untangled
                    let next = if line == 0 { 0 } else { line / 2 };
                    let (_, law, l_nose) = nose(next, row);
                    if line == 0 {
                        let ratio = dist(corner(line, row), corner(line, row + 1)) / l_nose;
                        (Curve::Line, scale(law, ratio))
                    } else {
                        let (pts, l) = te_half(line / 2 - 1);
                        (spline(pts), scale(law, l / l_nose))
                    }
                };
                edges.push(EdgeSpec {
                    vertices: [id(line, row, layer), id(line, row + 1, layer)],
                    curve,
                    distribution: law,
                });
            }
        }
        // the rings: the wall edges and the edges across
        for (k, e) in elements.iter().enumerate() {
            let s = &wall_spacing[k];
            let n = e.ring.len();
            let c = e.cove.unwrap_or(0);
            let walls = [
                (
                    TEL,
                    LOW,
                    e.ring[c + 1..e.lower].to_vec(),
                    s.cove,
                    s.corner[0],
                ),
                (
                    LOW,
                    LE,
                    e.ring[e.lower + 1..e.le].to_vec(),
                    s.corner[0],
                    s.le,
                ),
                (
                    LE,
                    UP,
                    e.ring[e.le + 1..e.upper].to_vec(),
                    s.le,
                    s.corner[1],
                ),
                (
                    UP,
                    TEU,
                    e.ring[e.upper + 1..n - 1].to_vec(),
                    s.corner[1],
                    s.te,
                ),
            ];
            for (a, b, pts, first, last) in walls {
                edges.push(EdgeSpec {
                    vertices: [wid(k, a, layer), wid(k, b, layer)],
                    curve: spline(pts),
                    distribution: Distribution::Tanh { first, last },
                });
            }
            if e.cove.is_some() {
                // the blunt trailing edge, then the cove from its corner
                edges.push(EdgeSpec {
                    vertices: [wid(k, TEU, layer), wid(k, TEM, layer)],
                    curve: Curve::Line,
                    distribution: Distribution::Uniform,
                });
                edges.push(EdgeSpec {
                    vertices: [wid(k, TEM, layer), wid(k, TEL, layer)],
                    curve: spline(e.ring[1..c].to_vec()),
                    distribution: Distribution::Tanh {
                        first: s.cove_te,
                        last: s.cove,
                    },
                });
            } else {
                for (a, b) in [(TEU, TEM), (TEM, TEL)] {
                    edges.push(EdgeSpec {
                        vertices: [wid(k, a, layer), wid(k, b, layer)],
                        curve: Curve::Line,
                        distribution: Distribution::Uniform,
                    });
                }
            }
            for vtx in 0..6 {
                edges.push(EdgeSpec {
                    vertices: [wid(k, vtx, layer), oid(k, vtx, layer)],
                    curve: Curve::Line,
                    distribution: Distribution::Geometric {
                        spacing: h,
                        at: EdgeEnd::Start,
                    },
                });
            }
        }
    }
    let quad = |col: usize, row: usize| -> [usize; 8] {
        let q = [
            (col, row),
            (col + 1, row),
            (col + 1, row + 1),
            (col, row + 1),
        ];
        let mut out = [0; 8];
        for layer in 0..2 {
            for (s, &(l, r)) in q.iter().enumerate() {
                out[4 * layer + s] = id(l, r, layer);
            }
        }
        out
    };
    let mut blocks = Vec::new();
    // the column of each background block, None for a ring block
    let mut block_col = Vec::new();
    for (row, &rc) in row_cells.iter().enumerate() {
        for (col, &cc) in col_cells.iter().enumerate().take(nl - 1) {
            if (row == 1 || row == 2) && col % 2 == 1 {
                continue; // an element
            }
            block_col.push(Some(col));
            blocks.push(BlockSpec {
                vertices: quad(col, row),
                cells: [cc, rc, p.span_cells],
                domain: "fluid".into(),
            });
        }
    }
    // ring blocks: i along the wall clockwise, j away from it
    let ring_blocks = |k: usize| {
        let (e, nt) = (&elements[k], p.thickness_cells);
        [
            (TEL, LOW, e.cells),
            (LOW, LE, nt),
            (LE, UP, nt),
            (UP, TEU, e.cells),
            (TEU, TEM, nt),
            (TEM, TEL, nt),
        ]
    };
    let mut ring_owner = Vec::new();
    for (k, e) in elements.iter().enumerate() {
        for (a, b, cells) in ring_blocks(k) {
            block_col.push(None);
            ring_owner.push(k);
            blocks.push(BlockSpec {
                vertices: [
                    wid(k, a, 0),
                    wid(k, b, 0),
                    oid(k, b, 0),
                    oid(k, a, 0),
                    wid(k, a, 1),
                    wid(k, b, 1),
                    oid(k, b, 1),
                    oid(k, a, 1),
                ],
                cells: [cells, e.ring_cells, p.span_cells],
                domain: "fluid".into(),
            });
        }
    }
    let side = |a: (usize, usize), b: (usize, usize)| {
        [
            id(a.0, a.1, 0),
            id(b.0, b.1, 0),
            id(b.0, b.1, 1),
            id(a.0, a.1, 1),
        ]
    };
    let mut patches = Vec::new();
    for (k, e) in elements.iter().enumerate() {
        patches.push(PatchSpec {
            name: format!("{}_wall", e.name),
            faces: ring_blocks(k)
                .iter()
                .map(|&(a, b, _)| [wid(k, a, 0), wid(k, b, 0), wid(k, b, 1), wid(k, a, 1)])
                .collect(),
        });
    }
    let mut inlet: Vec<[usize; 4]> = (0..4).map(|row| side((0, row), (0, row + 1))).collect();
    for col in 0..nl - 1 {
        inlet.push(side((col, 0), (col + 1, 0)));
        inlet.push(side((col, 4), (col + 1, 4)));
    }
    patches.push(PatchSpec {
        name: "inlet".into(),
        faces: inlet,
    });
    patches.push(PatchSpec {
        name: "outlet".into(),
        faces: (0..4)
            .map(|row| side((nl - 1, row), (nl - 1, row + 1)))
            .collect(),
    });
    for (layer, name) in [(0, "symmetry_low"), (1, "symmetry_high")] {
        patches.push(PatchSpec {
            name: name.into(),
            faces: blocks
                .iter()
                .map(|b: &BlockSpec| {
                    let o = 4 * layer;
                    [
                        b.vertices[o],
                        b.vertices[o + 1],
                        b.vertices[o + 2],
                        b.vertices[o + 3],
                    ]
                })
                .collect(),
        });
    }
    let mut topology = BlockTopology {
        junctions: Vec::new(),
        vertices,
        edges,
        blocks,
        patches,
        surfaces: Vec::new(),
    };
    let unoptimised = topology.clone();
    // the corners, shapes and spacings of the blocks off the rings move
    // toward orthogonal cells; the rings (on the wall patches) stay
    let layout = super::layout::optimize(
        &mut topology,
        &LayoutOptions {
            rigid_patches: (0..ne).collect(),
            ..LayoutOptions::default()
        },
    )?;
    // the block edges between the elements must not cut through them: the
    // TFI grid of the layout must not fold
    let grid = super::grid::build(&topology).map_err(|issues| {
        issues
            .iter()
            .map(|i| i.message.clone())
            .collect::<Vec<_>>()
            .join("; ")
    })?;
    if let Some(q) = grid.quality().iter().find(|q| q.min_jacobian <= 0.0) {
        let names: Vec<&str> = elements
            .iter()
            .map(|e| {
                if e.name == "main" {
                    "main element"
                } else {
                    e.name
                }
            })
            .collect();
        let place = match block_col[q.block] {
            None => {
                let k = ring_owner[q.block - (block_col.len() - ring_owner.len())];
                format!(
                    "in the wall ring around the {}: a smaller wall spacing or more surface \
                     cells may help",
                    names[k]
                )
            }
            Some(col) if col % 2 == 0 && col > 0 && col < nl - 2 => {
                let k = (col - 2) / 2;
                format!(
                    "in the gap between the {} and the {}: this layout needs the {}'s trailing \
                     edge about level with the {}'s nose; make that overlap more negative or the \
                     gap smaller",
                    names[k],
                    names[k + 1],
                    names[k],
                    names[k + 1]
                )
            }
            Some(_) => {
                "away from the gaps; a larger farfield distance or other cell counts may help"
                    .into()
            }
        };
        return Err(format!(
            "the grid of block {} folds {place} (a slot topology for overlapping elements is \
             not built yet)",
            q.block
        ));
    }
    let area = (x_max - x_min) * (y_max - y_min) - elements.iter().map(|e| e.area).sum::<f64>();
    Ok(Template {
        topology,
        analytic_volume: area * span,
        layout: Some(layout),
        unoptimised: Some(unoptimised),
    })
}

/// The last of `n` geometric cells from `h` that fill `length`.
fn last_cell(h: f64, n: usize, length: f64) -> f64 {
    let sum = |r: f64| (0..n).fold((0.0, 1.0), |(s, q), _| (s + q, q * r)).0 * h;
    let (mut lo, mut hi) = (1e-3, 2.0);
    while sum(hi) < length && hi < 1e6 {
        hi *= 2.0;
    }
    for _ in 0..100 {
        let m = 0.5 * (lo + hi);
        if sum(m) < length {
            lo = m;
        } else {
            hi = m;
        }
    }
    h * lo.powi(n as i32 - 1)
}
