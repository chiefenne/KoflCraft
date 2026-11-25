//! A multi-element airfoil (slat, main element, flap) in a box, extruded
//! along +z (BS-03, MB1g).
//!
//! **Layout.** An H-topology of 2 E + 1 block columns (E elements: a front
//! column, one column per element, a gap column between two elements and
//! a wake column) and four block rows (below, lower middle, upper middle,
//! above). Each element is a hole two blocks tall with six vertices on its
//! wall: the leading edge (the point farthest forward), a lower and an
//! upper corner (the lowest and highest point of the front `CORNER_REACH`
//! of each side), the two trailing-edge corners and the middle of the
//! blunt trailing edge. So the leading edge is a vertex of two blocks only
//! and the middle row line leaves it along the wall normal; three blocks
//! meet at the lower and upper corners, where the column line leaves the
//! wall vertically, normal to it. Between two elements the row lines are
//! splines that pass in front of the next element's nose where a straight
//! line would cut it. The Python prototype of the contour sliding
//! animation (2026-11-24) used one middle row with the corners at the
//! leading edge; that folds as soon as a flap is deflected.
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

    /// Ring points strictly between indices `from` and `to` (either order),
    /// from `from` to `to`.
    fn between(&self, from: usize, to: usize) -> Vec<P2> {
        if from < to {
            self.ring[from + 1..to].to_vec()
        } else {
            self.ring[to + 1..from].iter().rev().copied().collect()
        }
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
    let elements = place_all(main, slat, flap, p, c, far)?;
    let ne = elements.len();
    // column lines 0..nl (element k: lines 2k+1 and 2k+2), row lines 0..=4
    let nl = 2 * ne + 2;
    let (x_min, y_min, y_max) = (-far, -far, far);
    let last = &elements[ne - 1];
    let x_max = last.te_upper()[0].max(last.te_lower()[0]) + wake;
    for e in &elements {
        for q in &e.ring {
            if q[0] <= x_min + 0.5 * c || q[1] <= y_min + 0.5 * c || q[1] >= y_max - 0.5 * c {
                return Err("the farfield distance is too short for the elements".into());
            }
        }
    }
    let wall = |line: usize, row: usize| -> P2 {
        let e = &elements[(line - 1) / 2];
        match (line % 2 == 1, row) {
            (true, 1) => e.ring[e.lower],
            (true, 2) => e.ring[e.le],
            (true, _) => e.ring[e.upper],
            (false, 1) => e.te_lower(),
            (false, 2) => e.te_mid(),
            (false, _) => e.te_upper(),
        }
    };
    let corner = |line: usize, row: usize| -> P2 {
        match row {
            // at the box ends the middle rows open up like the front of a
            // C-grid, so their blocks are wedges, not long thin strips
            1..=3 if line == 0 => [x_min, wall(1, 2)[1] + (row as f64 - 2.0) * END_SPREAD * far],
            1..=3 if line == nl - 1 => [
                x_max,
                wall(nl - 2, 2)[1] + (row as f64 - 2.0) * END_SPREAD * far,
            ],
            1..=3 => wall(line, row),
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
                    let (a, b) = (wall(1, r)[0], wall(nl - 2, r)[0]);
                    let even = a + (b - a) * (line - 1) as f64 / (nl - 3) as f64;
                    wall(line, r)[0] + BOX_SPREAD * (even - wall(line, r)[0])
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
    let mut vertices = Vec::with_capacity(10 * nl);
    for z in [0.0, span] {
        for row in 0..5 {
            for line in 0..nl {
                vertices.push(v(corner(line, row), z));
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
    // the trailing edge
    struct Spacing {
        le: f64,
        corner: [f64; 2],
        te: f64,
    }
    let nt = p.thickness_cells as f64;
    let spacing: Vec<Spacing> = elements
        .iter()
        .map(|e| {
            let mean = e.length(e.lower, 0) / e.cells as f64;
            let (nose_lo, nose_up) = (e.length(e.lower, e.le), e.length(e.le, e.upper));
            let te_half = 0.5 * dist(e.te_upper(), e.te_lower());
            let le = (0.2 * mean).min(0.3 * nose_lo.min(nose_up) / nt);
            Spacing {
                le,
                corner: [
                    (nose_lo / nt).min(mean).max(1.5 * le),
                    (nose_up / nt).min(mean).max(1.5 * le),
                ],
                te: (0.2 * mean).max(te_half / nt),
            }
        })
        .collect();
    // the directions in which the three row lines leave each trailing edge
    // into the gap behind it: the middle one toward the next element's
    // middle corner (at most AIM off the chord), the others toward their
    // corners but at least FAN beside the middle one; behind the last
    // element, toward the box end
    let fans: Vec<[f64; 3]> = (0..ne)
        .map(|k| {
            let prev = &elements[k];
            let te = [prev.te_lower(), prev.te_mid(), prev.te_upper()];
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
            let mid = off(toward(te[1], next.ring[next.le])).clamp(-aim, aim);
            let lower = off(toward(te[0], next.ring[next.lower]));
            let upper = off(toward(te[2], next.ring[next.upper]));
            [
                base + lower.min(mid - fan).max(-limit),
                base + mid,
                base + upper.max(mid + fan).min(limit),
            ]
        })
        .collect();
    // the column line leaving a trailing-edge corner halves the angle
    // between the wall and the gap's lower (upper) row line
    let te_column = |k: usize, upper: bool| -> f64 {
        let e = &elements[k];
        let n = e.ring.len();
        let (corner, along) = if upper {
            (e.te_upper(), e.ring[n - 2])
        } else {
            (e.te_lower(), e.ring[1])
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
                let (curve, law) = if col == 0 {
                    let s = &spacing[0];
                    let at = match row {
                        1 => s.corner[0],
                        2 => s.le,
                        _ => s.corner[1],
                    };
                    (
                        Curve::Line,
                        Distribution::Geometric {
                            spacing: at,
                            at: EdgeEnd::End,
                        },
                    )
                } else if col == nl - 2 {
                    (
                        Curve::Line,
                        Distribution::Geometric {
                            spacing: spacing[ne - 1].te,
                            at: EdgeEnd::Start,
                        },
                    )
                } else if col % 2 == 1 {
                    if row == 2 {
                        continue; // inside the element
                    }
                    let k = (col - 1) / 2;
                    let (e, s) = (&elements[k], &spacing[k]);
                    let (from, to) = if row == 1 {
                        (e.lower, 0)
                    } else {
                        (e.upper, e.ring.len() - 1)
                    };
                    (
                        spline(e.between(from, to)),
                        Distribution::Tanh {
                            first: s.corner[if row == 1 { 0 } else { 1 }],
                            last: s.te,
                        },
                    )
                } else {
                    // a gap: pass in front of the next element's nose
                    let k = (col - 2) / 2;
                    let (next, s) = (&elements[k + 1], &spacing[k + 1]);
                    // leave the trailing edge in the fan and reach the next
                    // element across its wall
                    let (i, nrm) = match row {
                        1 => (next.lower, next.normal(next.lower)),
                        2 => (next.le, next.normal(next.le)),
                        _ => (next.upper, next.normal(next.upper)),
                    };
                    let q = next.ring[i];
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
                    let out = [turn.cos(), turn.sin()];
                    let (d1, d2) = (DEPART * d, ARRIVE * d);
                    let curve = spline(vec![
                        [a[0] + d1 * out[0], a[1] + d1 * out[1]],
                        [q[0] + d2 * nrm[0], q[1] + d2 * nrm[1]],
                    ]);
                    (
                        curve,
                        Distribution::Tanh {
                            first: spacing[k].te,
                            last: match row {
                                1 => s.corner[0],
                                2 => s.le,
                                _ => s.corner[1],
                            },
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
                let scaled = match law {
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
                edges.push(EdgeSpec {
                    vertices: [id(col, box_row, layer), id(col + 1, box_row, layer)],
                    curve: Curve::Line,
                    distribution: scaled,
                });
            }
        }
        // along the column lines
        for line in 0..nl {
            // from a trailing-edge corner the line leaves between the wall
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
                    spacing: h,
                    at: EdgeEnd::End,
                },
            });
            edges.push(EdgeSpec {
                vertices: [id(line, 3, layer), id(line, 4, layer)],
                curve: leave(3),
                distribution: Distribution::Geometric {
                    spacing: h,
                    at: EdgeEnd::Start,
                },
            });
            for row in [1, 2] {
                // the nose half of element k, its law and length
                let nose = |k: usize| {
                    let (e, s) = (&elements[k], &spacing[k]);
                    if row == 1 {
                        (
                            e.between(e.lower, e.le),
                            Distribution::Tanh {
                                first: s.corner[0],
                                last: s.le,
                            },
                            e.length(e.lower, e.le),
                        )
                    } else {
                        (
                            e.between(e.le, e.upper),
                            Distribution::Tanh {
                                first: s.le,
                                last: s.corner[1],
                            },
                            e.length(e.le, e.upper),
                        )
                    }
                };
                let (curve, law) = if line == nl - 1 || line == nl - 2 {
                    // the last trailing edge and the box end behind it
                    (Curve::Line, Distribution::Uniform)
                } else if line % 2 == 1 {
                    let (pts, law, _) = nose((line - 1) / 2);
                    (spline(pts), law)
                } else {
                    // the box front or a trailing edge: the law of the nose
                    // across the block, scaled to this edge's length, so that
                    // the interpolation inside the block stays untangled
                    let next = if line == 0 { 0 } else { line / 2 };
                    let (_, law, l_nose) = nose(next);
                    let ratio = dist(corner(line, row), corner(line, row + 1)) / l_nose;
                    let Distribution::Tanh { first, last } = law else {
                        unreachable!()
                    };
                    (
                        Curve::Line,
                        Distribution::Tanh {
                            first: first * ratio,
                            last: last * ratio,
                        },
                    )
                };
                edges.push(EdgeSpec {
                    vertices: [id(line, row, layer), id(line, row + 1, layer)],
                    curve,
                    distribution: law,
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
    let mut block_col = Vec::new();
    for (row, &rc) in row_cells.iter().enumerate() {
        for (col, &cc) in col_cells.iter().enumerate().take(nl - 1) {
            if (row == 1 || row == 2) && col % 2 == 1 {
                continue; // an element
            }
            block_col.push(col);
            blocks.push(BlockSpec {
                vertices: quad(col, row),
                cells: [cc, rc, p.span_cells],
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
        let l = 2 * k + 1;
        patches.push(PatchSpec {
            name: format!("{}_wall", e.name),
            faces: vec![
                side((l, 1), (l + 1, 1)),
                side((l, 3), (l + 1, 3)),
                side((l, 1), (l, 2)),
                side((l, 2), (l, 3)),
                side((l + 1, 1), (l + 1, 2)),
                side((l + 1, 2), (l + 1, 3)),
            ],
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
    let topology = BlockTopology {
        vertices,
        edges,
        blocks,
        patches,
        surfaces: Vec::new(),
    };
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
        let col = block_col[q.block];
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
        let place = if col % 2 == 0 && col > 0 && col < nl - 2 {
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
        } else {
            "away from the gaps; a larger farfield distance or other cell counts may help".into()
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
    })
}
