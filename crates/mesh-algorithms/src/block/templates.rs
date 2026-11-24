//! Parametric block topologies for standard cases (BS-03). A template only
//! writes a topology; the grid is built from it like from any other.
//! MB1a: the straight channel and the pipe O-grid with an H-core.

use mesh_core::block::{
    BlockSpec, BlockTopology, Curve, Distribution, EdgeEnd, EdgeSpec, PatchSpec, Surface,
    SurfaceSpec,
};
use mesh_core::geometry::Vec3;

/// A template's topology with the analytic volume of the shape it meshes.
pub struct Template {
    pub topology: BlockTopology,
    pub analytic_volume: f64,
}

fn positive(v: f64, what: &str) -> Result<f64, String> {
    if v.is_finite() && v > 0.0 {
        Ok(v)
    } else {
        Err(format!("{what} must be positive and finite"))
    }
}

fn counts(cells: &[usize]) -> Result<(), String> {
    if cells.contains(&0) {
        Err("cell counts must be at least 1".into())
    } else {
        Ok(())
    }
}

/// A straight rectangular channel along +x: one block, inlet at x = 0,
/// outlet at x = length, the four other faces `walls`. With `wall_spacing`
/// the cross-section edges cluster at both walls (two-sided tanh).
pub fn channel(
    size: [f64; 3],
    cells: [usize; 3],
    wall_spacing: Option<f64>,
) -> Result<Template, String> {
    let [l, w, h] = [
        positive(size[0], "length")?,
        positive(size[1], "width")?,
        positive(size[2], "height")?,
    ];
    counts(&cells)?;
    let vertices = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(l, 0.0, 0.0),
        Vec3::new(l, w, 0.0),
        Vec3::new(0.0, w, 0.0),
        Vec3::new(0.0, 0.0, h),
        Vec3::new(l, 0.0, h),
        Vec3::new(l, w, h),
        Vec3::new(0.0, w, h),
    ];
    let mut edges = Vec::new();
    if let Some(s) = wall_spacing {
        positive(s, "wall spacing")?;
        let law = Distribution::Tanh { first: s, last: s };
        for v in [
            [0, 3],
            [1, 2],
            [4, 7],
            [5, 6],
            [0, 4],
            [1, 5],
            [3, 7],
            [2, 6],
        ] {
            edges.push(EdgeSpec {
                vertices: v,
                curve: Curve::Line,
                distribution: law,
            });
        }
    }
    Ok(Template {
        topology: BlockTopology {
            vertices,
            edges,
            blocks: vec![BlockSpec {
                vertices: [0, 1, 2, 3, 4, 5, 6, 7],
                cells,
                domain: "fluid".into(),
            }],
            patches: vec![
                PatchSpec {
                    name: "inlet".into(),
                    faces: vec![[0, 3, 7, 4]],
                },
                PatchSpec {
                    name: "outlet".into(),
                    faces: vec![[1, 2, 6, 5]],
                },
                PatchSpec {
                    name: "walls".into(),
                    faces: vec![[0, 1, 5, 4], [3, 2, 6, 7], [0, 1, 2, 3], [4, 5, 6, 7]],
                },
            ],
            surfaces: Vec::new(),
        },
        analytic_volume: l * w * h,
    })
}

/// A straight circular pipe along +x: an H-core block and four outer
/// blocks (the O-grid "butterfly"), inlet at x = 0, outlet at x = length,
/// wall on the cylinder. `core_fraction` is the radius of the core's
/// corners over the pipe radius; `core_cells` the cells along a core side
/// (a quarter of the circumference); `radial_cells` from core to wall.
/// With `wall_spacing` the radial edges grow geometrically from the wall.
pub fn pipe_ogrid(
    radius: f64,
    length: f64,
    core_fraction: f64,
    cells: [usize; 3],
    wall_spacing: Option<f64>,
) -> Result<Template, String> {
    let r = positive(radius, "radius")?;
    let l = positive(length, "length")?;
    let stations = [
        Station::straight(Vec3::ZERO),
        Station::straight(Vec3::new(l, 0.0, 0.0)),
    ];
    let topology = swept_ogrid(r, core_fraction, &stations, &[None], cells, wall_spacing)?;
    Ok(Template {
        topology,
        analytic_volume: std::f64::consts::PI * r * r * l,
    })
}

/// A circular pipe elbow: an optional straight inlet leg along +x ending at
/// the origin, a bend of `angle_deg` about the axis through (0, R, ·) along
/// +z with centreline radius `bend_radius`, and an optional straight outlet
/// leg. The cross-section is the pipe O-grid; the bend is split into
/// segments of at most 30° so that every longitudinal block edge is an
/// exact circle on the torus. `cells` is [core, radial, axial] with the
/// axial cells shared out by centreline length.
#[allow(clippy::too_many_arguments)]
pub fn elbow(
    radius: f64,
    bend_radius: f64,
    angle_deg: f64,
    inlet_length: f64,
    outlet_length: f64,
    core_fraction: f64,
    cells: [usize; 3],
    wall_spacing: Option<f64>,
) -> Result<Template, String> {
    let r = positive(radius, "radius")?;
    let big = positive(bend_radius, "bend radius")?;
    if big <= r * 1.05 {
        return Err("bend radius must exceed 1.05 pipe radii".into());
    }
    let theta = positive(angle_deg, "bend angle")?;
    if theta > 180.0 {
        return Err("bend angle must be at most 180°".into());
    }
    for (v, what) in [
        (inlet_length, "inlet length"),
        (outlet_length, "outlet length"),
    ] {
        if !(v.is_finite() && v >= 0.0) {
            return Err(format!("{what} must be zero or positive and finite"));
        }
    }
    let segments = (theta / 30.0).ceil() as usize;
    let station_at = |phi: f64| {
        let (s, c) = phi.to_radians().sin_cos();
        Station {
            origin: Vec3::new(big * s, big * (1.0 - c), 0.0),
            normal: Vec3::new(-s, c, 0.0),
            tangent: Vec3::new(c, s, 0.0),
        }
    };
    let mut stations = Vec::new();
    let mut bends = Vec::new();
    if inlet_length > 0.0 {
        stations.push(Station::straight(Vec3::new(-inlet_length, 0.0, 0.0)));
        bends.push(None);
    }
    for m in 0..=segments {
        stations.push(station_at(theta * m as f64 / segments as f64));
        if m < segments {
            bends.push(Some(BendAxis {
                centre: Vec3::new(0.0, big, 0.0),
                radius: big,
                mid_angle: theta * (m as f64 + 0.5) / segments as f64,
            }));
        }
    }
    if outlet_length > 0.0 {
        let end = station_at(theta);
        stations.push(Station {
            origin: end.origin + end.tangent * outlet_length,
            ..end
        });
        bends.push(None);
    }
    let topology = swept_ogrid(r, core_fraction, &stations, &bends, cells, wall_spacing)?;
    Ok(Template {
        topology,
        analytic_volume: std::f64::consts::PI
            * r
            * r
            * (inlet_length + outlet_length + big * theta.to_radians()),
    })
}

/// A cross-section plane of a swept O-grid: centre, the in-plane direction
/// of the cross-section's local y (angle 0) and the sweep direction; the
/// local z is `tangent × normal`.
#[derive(Clone, Copy)]
struct Station {
    origin: Vec3,
    normal: Vec3,
    tangent: Vec3,
}

impl Station {
    fn straight(origin: Vec3) -> Self {
        Self {
            origin,
            normal: Vec3::new(0.0, 1.0, 0.0),
            tangent: Vec3::new(1.0, 0.0, 0.0),
        }
    }

    /// The point at radius `rad` and angle `deg` in the cross-section.
    fn at(&self, rad: f64, deg: f64) -> Vec3 {
        let t = deg.to_radians();
        let binormal = self.tangent.cross(self.normal);
        self.origin + (self.normal * t.cos() + binormal * t.sin()) * rad
    }
}

/// A bend segment between two stations: the bend axis passes through
/// `centre` along +z with centreline radius `radius`, and `mid_angle` is the bend angle (degrees) halfway
/// along the segment.
#[derive(Clone, Copy)]
struct BendAxis {
    centre: Vec3,
    radius: f64,
    mid_angle: f64,
}

/// An O-grid with an H-core swept through `stations`; segment s joins
/// station s and s + 1 straight (`None`) or along a bend. Inlet is the
/// first station, outlet the last, walls the swept outer ring.
fn swept_ogrid(
    r: f64,
    core_fraction: f64,
    stations: &[Station],
    bends: &[Option<BendAxis>],
    cells: [usize; 3],
    wall_spacing: Option<f64>,
) -> Result<BlockTopology, String> {
    if !(core_fraction > 0.1 && core_fraction < 0.9) {
        return Err("core_fraction must lie between 0.1 and 0.9".into());
    }
    counts(&cells)?;
    let spacing = wall_spacing
        .map(|s| positive(s, "wall spacing"))
        .transpose()?;
    let [core, radial, axial] = cells;
    let rc = core_fraction * r;
    let segments = stations.len() - 1;
    // the axial cells shared out by centreline length, at least one each
    let lengths: Vec<f64> = (0..segments)
        .map(|s| match bends[s] {
            None => (stations[s + 1].origin - stations[s].origin).length(),
            Some(b) => {
                let a = stations[s].origin - b.centre;
                let c = stations[s + 1].origin - b.centre;
                a.length() * a.cross(c).length().atan2(a.dot(c))
            }
        })
        .collect();
    if axial < segments {
        return Err(format!(
            "axial cells must be at least the number of segments ({segments})"
        ));
    }
    let seg_cells = share(axial, &lengths);
    // vertex 8 p + m: inner corner m at 45° + 90° m; 8 p + 4 + m: outer
    let mut vertices = Vec::with_capacity(8 * stations.len());
    for st in stations {
        for m in 0..4 {
            vertices.push(st.at(rc, 45.0 + 90.0 * m as f64));
        }
        for m in 0..4 {
            vertices.push(st.at(r, 45.0 + 90.0 * m as f64));
        }
    }
    let c = |p: usize, m: usize| 8 * p + m % 4;
    let o = |p: usize, m: usize| 8 * p + 4 + m % 4;
    let mut blocks = Vec::new();
    let mut edges = Vec::new();
    let mut walls = Vec::new();
    let mut surfaces: Vec<SurfaceSpec> = Vec::new();
    for s in 0..segments {
        let (p0, p1) = (s, s + 1);
        blocks.push(BlockSpec {
            // counter-clockwise in the cross-section from 225°, k along the sweep
            vertices: [
                c(p0, 2),
                c(p0, 3),
                c(p0, 0),
                c(p0, 1),
                c(p1, 2),
                c(p1, 3),
                c(p1, 0),
                c(p1, 1),
            ],
            cells: [core, core, seg_cells[s]],
            domain: "fluid".into(),
        });
        for q in 0..4 {
            let m1 = q + 3; // the corner at 90° q − 45°
            blocks.push(BlockSpec {
                vertices: [
                    c(p0, m1),
                    o(p0, m1),
                    o(p0, q),
                    c(p0, q),
                    c(p1, m1),
                    o(p1, m1),
                    o(p1, q),
                    c(p1, q),
                ],
                cells: [radial, core, seg_cells[s]],
                domain: "fluid".into(),
            });
            walls.push([o(p0, m1), o(p0, q), o(p1, q), o(p1, m1)]);
        }
        let ring = walls[walls.len() - 4..].to_vec();
        match bends[s] {
            None => surfaces.push(SurfaceSpec {
                name: format!("cylinder_{s}"),
                surface: Surface::Cylinder {
                    point: stations[s].origin,
                    axis: stations[s].tangent,
                    radius: r,
                },
                faces: ring,
            }),
            Some(b) => match surfaces.iter_mut().find(|x| x.name == "torus") {
                Some(t) => t.faces.extend(ring),
                None => surfaces.push(SurfaceSpec {
                    name: "torus".into(),
                    surface: Surface::Torus {
                        centre: b.centre,
                        axis: Vec3::new(0.0, 0.0, 1.0),
                        major: b.radius,
                        minor: r,
                    },
                    faces: ring,
                }),
            },
        }
        if let Some(b) = bends[s] {
            // longitudinal edges: circles about the bend axis
            let (sn, cs) = b.mid_angle.to_radians().sin_cos();
            let mid = Station {
                origin: b.centre + Vec3::new(sn, -cs, 0.0) * b.radius,
                normal: Vec3::new(-sn, cs, 0.0),
                tangent: Vec3::new(cs, sn, 0.0),
            };
            for m in 0..4 {
                let deg = 45.0 + 90.0 * m as f64;
                for (v0, v1, rad) in [(c(p0, m), c(p1, m), rc), (o(p0, m), o(p1, m), r)] {
                    edges.push(EdgeSpec {
                        vertices: [v0, v1],
                        curve: Curve::Arc {
                            through: mid.at(rad, deg),
                        },
                        distribution: Distribution::Uniform,
                    });
                }
            }
        }
    }
    for (p, st) in stations.iter().enumerate() {
        for q in 0..4 {
            let m1 = q + 3;
            edges.push(EdgeSpec {
                vertices: [o(p, m1), o(p, q)],
                curve: Curve::Arc {
                    through: st.at(r, 90.0 * q as f64),
                },
                distribution: Distribution::Uniform,
            });
            if let Some(h) = spacing {
                edges.push(EdgeSpec {
                    vertices: [c(p, q), o(p, q)],
                    curve: Curve::Line,
                    distribution: Distribution::Geometric {
                        spacing: h,
                        at: EdgeEnd::End,
                    },
                });
            }
        }
    }
    let last = segments;
    let ring = |p: usize| {
        let mut f = vec![[c(p, 2), c(p, 3), c(p, 0), c(p, 1)]];
        for q in 0..4 {
            f.push([c(p, q + 3), o(p, q + 3), o(p, q), c(p, q)]);
        }
        f
    };
    Ok(BlockTopology {
        vertices,
        edges,
        blocks,
        patches: vec![
            PatchSpec {
                name: "inlet".into(),
                faces: ring(0),
            },
            PatchSpec {
                name: "outlet".into(),
                faces: ring(last),
            },
            PatchSpec {
                name: "walls".into(),
                faces: walls,
            },
        ],
        surfaces,
    })
}

/// `total` cells shared out in proportion to `weights`, at least one each,
/// by largest remainder.
fn share(total: usize, weights: &[f64]) -> Vec<usize> {
    let n = weights.len();
    let sum: f64 = weights.iter().sum();
    let spare = (total - n) as f64;
    let exact: Vec<f64> = weights.iter().map(|w| spare * w / sum).collect();
    let mut out: Vec<usize> = exact.iter().map(|e| 1 + e.floor() as usize).collect();
    let mut left = total - out.iter().sum::<usize>();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        let fa = exact[a] - exact[a].floor();
        let fb = exact[b] - exact[b].floor();
        fb.total_cmp(&fa).then(a.cmp(&b))
    });
    for &i in order.iter().cycle() {
        if left == 0 {
            break;
        }
        out[i] += 1;
        left -= 1;
    }
    out
}

/// A NACA 4-digit section with the closed trailing edge (last thickness
/// coefficient −0.1036), chord 1 from the leading edge at the origin to
/// the trailing edge at (1, 0). `side` +1 is the upper surface, −1 the
/// lower; `xc` ∈ [0, 1] is the position along the camber line.
#[derive(Clone, Copy, Debug)]
pub struct Naca4 {
    camber: f64,
    position: f64,
    thickness: f64,
}

impl Naca4 {
    pub fn parse(digits: &str) -> Result<Self, String> {
        let d: Vec<u32> = digits.chars().filter_map(|c| c.to_digit(10)).collect();
        if digits.len() != 4 || d.len() != 4 {
            return Err(format!("NACA designation {digits:?} is not four digits"));
        }
        let camber = d[0] as f64 / 100.0;
        let position = d[1] as f64 / 10.0;
        let thickness = (10 * d[2] + d[3]) as f64 / 100.0;
        if thickness <= 0.0 || thickness > 0.4 {
            return Err("NACA thickness must lie between 1 and 40 % of the chord".into());
        }
        if (camber > 0.0) != (position > 0.0) {
            return Err("a cambered NACA section needs a camber position, and the reverse".into());
        }
        Ok(Self {
            camber,
            position,
            thickness,
        })
    }

    /// The surface point of `side` at camber-line position `xc`.
    pub fn point(&self, xc: f64, side: f64) -> [f64; 2] {
        let t = self.thickness;
        let yt = 5.0
            * t
            * (0.2969 * xc.sqrt() - 0.1260 * xc - 0.3516 * xc * xc + 0.2843 * xc.powi(3)
                - 0.1036 * xc.powi(4));
        let (m, p) = (self.camber, self.position);
        let (yc, slope) = if m == 0.0 {
            (0.0, 0.0)
        } else if xc < p {
            (
                m / (p * p) * (2.0 * p * xc - xc * xc),
                2.0 * m / (p * p) * (p - xc),
            )
        } else {
            (
                m / ((1.0 - p) * (1.0 - p)) * (1.0 - 2.0 * p + 2.0 * p * xc - xc * xc),
                2.0 * m / ((1.0 - p) * (1.0 - p)) * (p - xc),
            )
        };
        let th = slope.atan();
        [xc - side * yt * th.sin(), yc + side * yt * th.cos()]
    }

    /// The area of the section (chord 1), by the shoelace formula on a
    /// fine cosine-spaced polygon.
    pub fn area(&self) -> f64 {
        let n = 20_000;
        let mut ring = Vec::with_capacity(2 * n);
        for i in 0..=n {
            let xc = 0.5 * (1.0 - (std::f64::consts::PI * i as f64 / n as f64).cos());
            ring.push(self.point(xc, -1.0));
        }
        for i in (1..n).rev() {
            let xc = 0.5 * (1.0 - (std::f64::consts::PI * i as f64 / n as f64).cos());
            ring.push(self.point(xc, 1.0));
        }
        let mut a = 0.0;
        for k in 0..ring.len() {
            let (p, q) = (ring[k], ring[(k + 1) % ring.len()]);
            a += p[0] * q[1] - q[0] * p[1];
        }
        0.5 * a
    }
}

impl Naca4 {
    /// The section as points of both sides, cosine-spaced along the
    /// camber line, with the area of the analytic shape.
    pub fn section(&self) -> Section {
        let side = |sign: f64| -> Vec<[f64; 2]> {
            (0..=AIRFOIL_POINTS)
                .map(|i| {
                    let xc = 0.5
                        * (1.0 - (std::f64::consts::PI * i as f64 / AIRFOIL_POINTS as f64).cos());
                    self.point(xc, sign)
                })
                .collect()
        };
        Section {
            upper: side(1.0),
            lower: side(-1.0),
            area: self.area(),
            trailing_edge_gap: 0.0,
        }
    }
}

/// An airfoil section of unit chord: the leading edge at the origin and
/// the closed trailing edge at (1, 0); each side runs from the leading to
/// the trailing edge. The wall of the C-grid is the cubic spline through
/// these points.
#[derive(Clone, Debug, PartialEq)]
pub struct Section {
    pub upper: Vec<[f64; 2]>,
    pub lower: Vec<[f64; 2]>,
    /// The area of the shape (analytic for NACA, of the polygon of the
    /// points for a coordinate file).
    pub area: f64,
    /// The trailing-edge gap of the file over the chord before it was
    /// closed (0 for NACA).
    pub trailing_edge_gap: f64,
}

impl Section {
    /// Read a coordinate file in Selig format (a name line, then x y from
    /// the trailing edge over the upper side to the leading edge and back
    /// over the lower side) or Lednicer format (a name line, the point
    /// counts of the two sides, then each side from the leading edge).
    /// The section is moved to the leading edge at the origin, turned and
    /// scaled to unit chord, and a blunt trailing edge is closed by a
    /// linear thickness correction.
    pub fn from_coordinates(text: &str) -> Result<Self, String> {
        let rows: Vec<Option<[f64; 2]>> = text
            .lines()
            .map(|l| {
                let v: Vec<f64> = l
                    .split_whitespace()
                    .map_while(|w| w.parse::<f64>().ok())
                    .collect();
                (v.len() >= 2 && l.split_whitespace().count() == 2).then(|| [v[0], v[1]])
            })
            .collect();
        let numbers: Vec<(usize, [f64; 2])> = rows
            .iter()
            .enumerate()
            .filter_map(|(i, r)| r.map(|p| (i, p)))
            .collect();
        if numbers.len() < 8 {
            return Err("a coordinate file needs at least eight points".into());
        }
        let (first, rest) = (numbers[0].1, &numbers[1..]);
        let points: Vec<[f64; 2]> = rest.iter().map(|&(_, p)| p).collect();
        let (upper, lower) = if first[0] > 1.5 && first[1] > 1.5 {
            // Lednicer: counts, then upper and lower from the leading edge
            let (nu, nl) = (first[0] as usize, first[1] as usize);
            if first[0].fract() != 0.0 || first[1].fract() != 0.0 || nu + nl != points.len() {
                return Err(format!(
                    "the Lednicer counts {nu} + {nl} do not match {} points",
                    points.len()
                ));
            }
            (points[..nu].to_vec(), points[nu..].to_vec())
        } else {
            // Selig: one loop from the trailing edge; split at the leading
            // edge, the point farthest from the trailing-edge midpoint
            let all: Vec<[f64; 2]> = numbers.iter().map(|&(_, p)| p).collect();
            let (a, b) = (all[0], all[all.len() - 1]);
            let te = [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0];
            let le = (0..all.len())
                .max_by(|&i, &j| {
                    let d = |p: [f64; 2]| (p[0] - te[0]).hypot(p[1] - te[1]);
                    d(all[i]).total_cmp(&d(all[j]))
                })
                .unwrap_or(0);
            let mut upper = all[..=le].to_vec();
            upper.reverse();
            (upper, all[le..].to_vec())
        };
        Self::normalise(upper, lower)
    }

    /// Both sides from the leading edge: to unit chord and closed.
    fn normalise(mut upper: Vec<[f64; 2]>, mut lower: Vec<[f64; 2]>) -> Result<Self, String> {
        if upper.len() < 4 || lower.len() < 4 {
            return Err("each side of the section needs at least four points".into());
        }
        if upper[0] != lower[0] {
            return Err("the two sides do not start at the same leading-edge point".into());
        }
        if upper
            .iter()
            .chain(&lower)
            .any(|p| !(p[0].is_finite() && p[1].is_finite()))
        {
            return Err("a coordinate is not finite".into());
        }
        let le = upper[0];
        let (tu, tl) = (upper[upper.len() - 1], lower[lower.len() - 1]);
        let te = [(tu[0] + tl[0]) / 2.0, (tu[1] + tl[1]) / 2.0];
        let chord = (te[0] - le[0]).hypot(te[1] - le[1]);
        if chord.is_nan() || chord <= 0.0 {
            return Err("the section has no chord".into());
        }
        let (c, s) = ((te[0] - le[0]) / chord, (te[1] - le[1]) / chord);
        let map = |p: [f64; 2]| {
            let (dx, dy) = (p[0] - le[0], p[1] - le[1]);
            [(c * dx + s * dy) / chord, (-s * dx + c * dy) / chord]
        };
        upper.iter_mut().for_each(|p| *p = map(*p));
        lower.iter_mut().for_each(|p| *p = map(*p));
        let gap = upper[upper.len() - 1][1] - lower[lower.len() - 1][1];
        if gap < -1e-9 {
            return Err(
                "the upper side ends below the lower side; is the file upside down?".into(),
            );
        }
        for (side, sign) in [(&mut upper, -1.0), (&mut lower, 1.0)] {
            for p in side.iter_mut() {
                p[1] += sign * 0.5 * gap * p[0];
            }
            let n = side.len();
            side[n - 1] = [1.0, 0.0];
            side[0] = [0.0, 0.0];
            if side
                .windows(2)
                .any(|w| (w[1][0] - w[0][0]).hypot(w[1][1] - w[0][1]) <= 0.0)
            {
                return Err("the section repeats a point".into());
            }
        }
        let ring: Vec<[f64; 2]> = lower
            .iter()
            .chain(upper.iter().rev().skip(1).take(upper.len() - 2))
            .copied()
            .collect();
        let mut area = 0.0;
        for k in 0..ring.len() {
            let (p, q) = (ring[k], ring[(k + 1) % ring.len()]);
            area += p[0] * q[1] - q[0] * p[1];
        }
        if area <= 0.0 {
            return Err("the section's sides are crossed or in the wrong order".into());
        }
        Ok(Self {
            upper,
            lower,
            area: 0.5 * area,
            trailing_edge_gap: gap,
        })
    }
}

/// Parameters of the extruded airfoil C-grid.
#[derive(Clone, Copy, Debug)]
pub struct AirfoilCgrid {
    pub chord: f64,
    pub span: f64,
    /// Radius of the C's half circle about the trailing edge.
    pub farfield: f64,
    /// Length of the wake blocks behind the trailing edge.
    pub wake: f64,
    /// Cells along each side of the airfoil, normal to it, along the wake
    /// and along the span.
    pub cells: [usize; 4],
    /// First cell height at the airfoil and on the wake cut.
    pub wall_spacing: f64,
    /// Cell lengths along the surface at the leading and trailing edges;
    /// a fifth of the mean cell length when absent.
    pub le_spacing: Option<f64>,
    pub te_spacing: Option<f64>,
}

/// Where along the chord the C-grid splits each side into a front
/// and a rear block (the airfoil edges between two vertices must differ
/// on the two sides, so the leading edge cannot be the only split).
const AIRFOIL_SPLIT: f64 = 0.3;
/// Spline points per airfoil side.
const AIRFOIL_POINTS: usize = 160;

/// An airfoil C-grid extruded along +z: the section (chord along +x from
/// the leading edge at the origin) inside a C of a half circle about the
/// trailing edge and a wake to `trailing edge + wake`; eight blocks, a
/// front and a rear block on each side and a wake block above and below
/// the wake cut. Patches: `airfoil_wall`, `inlet` (the C and the lines
/// above and below the wake), `outlet` (the wake end), `symmetry_low`
/// and `symmetry_high` (the two span ends).
pub fn airfoil_cgrid(section: &Section, p: &AirfoilCgrid) -> Result<Template, String> {
    let c = positive(p.chord, "chord")?;
    let span = positive(p.span, "span")?;
    let big = positive(p.farfield, "farfield radius")?;
    let wake = positive(p.wake, "wake length")?;
    if big < 2.0 * c {
        return Err("the farfield radius must be at least two chords".into());
    }
    counts(&p.cells)?;
    let [n_side, n_normal, n_wake, n_span] = p.cells;
    if n_side < 2 {
        return Err("each airfoil side needs at least two cells".into());
    }
    let h = positive(p.wall_spacing, "wall spacing")?;
    let side_points = |side: f64| {
        if side > 0.0 {
            &section.upper
        } else {
            &section.lower
        }
    };
    // the split vertex: the point of each side nearest AIRFOIL_SPLIT
    let split_of = |side: f64| {
        let pts = side_points(side);
        (1..pts.len() - 1)
            .min_by(|&i, &j| {
                (pts[i][0] - AIRFOIL_SPLIT)
                    .abs()
                    .total_cmp(&(pts[j][0] - AIRFOIL_SPLIT).abs())
            })
            .unwrap_or(1)
    };
    let pt = |side: f64, i: usize, z: f64| {
        let [x, y] = side_points(side)[i];
        Vec3::new(c * x, c * y, z)
    };
    // 2D vertices: 0 LE, 1 TE, 2 upper split, 3 lower split, 4 front,
    // 5 C top above TE, 6 C bottom, 7 C above the upper split, 8 C below
    // the lower split, 9 wake end, 10 its top, 11 its bottom
    let te_x = c;
    let on_c = |x: f64, sign: f64| {
        let dx = te_x - x;
        Vec3::new(x, sign * (big * big - dx * dx).sqrt(), 0.0)
    };
    let (ku, kl) = (split_of(1.0), split_of(-1.0));
    if ku < 2 || kl < 2 {
        return Err("the section has too few points before the split".into());
    }
    let su = pt(1.0, ku, 0.0);
    let sl = pt(-1.0, kl, 0.0);
    let plane = [
        Vec3::ZERO,
        Vec3::new(te_x, 0.0, 0.0),
        su,
        sl,
        Vec3::new(te_x - big, 0.0, 0.0),
        Vec3::new(te_x, big, 0.0),
        Vec3::new(te_x, -big, 0.0),
        on_c(su.x, 1.0),
        on_c(sl.x, -1.0),
        Vec3::new(te_x + wake, 0.0, 0.0),
        Vec3::new(te_x + wake, big, 0.0),
        Vec3::new(te_x + wake, -big, 0.0),
    ];
    let nv = plane.len();
    let mut vertices = Vec::with_capacity(2 * nv);
    for z in [0.0, span] {
        vertices.extend(plane.iter().map(|&v| Vec3::new(v.x, v.y, z)));
    }
    let at = |v: usize, layer: usize| v + nv * layer;
    let quad = |q: [usize; 4]| -> [usize; 8] {
        [
            at(q[0], 0),
            at(q[1], 0),
            at(q[2], 0),
            at(q[3], 0),
            at(q[0], 1),
            at(q[1], 1),
            at(q[2], 1),
            at(q[3], 1),
        ]
    };
    // spline points of each side strictly between two of its indices
    let through = |side: f64, lo: usize, hi: usize, z: f64| -> Vec<Vec3> {
        (lo + 1..hi).map(|i| pt(side, i, z)).collect()
    };
    let poly_len = |side: f64, lo: usize, hi: usize| {
        (lo..hi)
            .map(|i| (pt(side, i + 1, 0.0) - pt(side, i, 0.0)).length())
            .sum::<f64>()
    };
    let mut edges = Vec::new();
    let mut blocks = Vec::new();
    let fluid = || "fluid".to_string();
    for (s, side) in [1.0, -1.0].into_iter().enumerate() {
        let k = if s == 0 { ku } else { kl };
        let last = side_points(side).len() - 1;
        let (front, rear) = (poly_len(side, 0, k), poly_len(side, k, last));
        let split = share(n_side, &[front, rear]);
        let mean = (front + rear) / n_side as f64;
        let le = p.le_spacing.unwrap_or(0.2 * mean);
        let te = p.te_spacing.unwrap_or(0.2 * mean);
        positive(le, "leading-edge spacing")?;
        positive(te, "trailing-edge spacing")?;
        let (split_v, c_split) = if s == 0 { (2, 7) } else { (3, 8) };
        let c_end = if s == 0 { 5 } else { 6 };
        // the C's arcs, with spacings scaled by the length ratio so that
        // they carry the airfoil's relative distribution
        // angles about the trailing edge in [0, 2π): the C runs from π/2
        // over π to 3π/2
        let angle = |a: Vec3| a.y.atan2(a.x - te_x).rem_euclid(std::f64::consts::TAU);
        let arc_mid = |a: Vec3, b: Vec3| {
            let t = 0.5 * (angle(a) + angle(b));
            Vec3::new(te_x + big * t.cos(), big * t.sin(), 0.0)
        };
        let arc_length = |a: Vec3, b: Vec3| big * (angle(b) - angle(a)).abs();
        let lf = arc_length(plane[4], plane[c_split]);
        let lr = arc_length(plane[c_split], plane[c_end]);
        let mid_spacing = mean;
        for layer in 0..2 {
            let z = if layer == 0 { 0.0 } else { span };
            let lift = |v: Vec3| Vec3::new(v.x, v.y, z);
            edges.push(EdgeSpec {
                vertices: [at(0, layer), at(split_v, layer)],
                curve: Curve::Spline {
                    through: through(side, 0, k, z),
                },
                distribution: Distribution::Tanh {
                    first: le,
                    last: mid_spacing,
                },
            });
            edges.push(EdgeSpec {
                vertices: [at(split_v, layer), at(1, layer)],
                curve: Curve::Spline {
                    through: through(side, k, last, z),
                },
                distribution: Distribution::Tanh {
                    first: mid_spacing,
                    last: te,
                },
            });
            edges.push(EdgeSpec {
                vertices: [at(4, layer), at(c_split, layer)],
                curve: Curve::Arc {
                    through: lift(arc_mid(plane[4], plane[c_split])),
                },
                distribution: Distribution::Tanh {
                    first: le * lf / front,
                    last: mid_spacing * lf / front,
                },
            });
            edges.push(EdgeSpec {
                vertices: [at(c_split, layer), at(c_end, layer)],
                curve: Curve::Arc {
                    through: lift(arc_mid(plane[c_split], plane[c_end])),
                },
                distribution: Distribution::Tanh {
                    first: mid_spacing * lr / rear,
                    last: te * lr / rear,
                },
            });
            // normal edges from the wall
            for (a, b) in [
                (split_v, c_split),
                (1, c_end),
                (9, if s == 0 { 10 } else { 11 }),
            ] {
                edges.push(EdgeSpec {
                    vertices: [at(a, layer), at(b, layer)],
                    curve: Curve::Line,
                    distribution: Distribution::Geometric {
                        spacing: h,
                        at: EdgeEnd::Start,
                    },
                });
            }
            if s == 0 {
                edges.push(EdgeSpec {
                    vertices: [at(0, layer), at(4, layer)],
                    curve: Curve::Line,
                    distribution: Distribution::Geometric {
                        spacing: h,
                        at: EdgeEnd::Start,
                    },
                });
                // the wake cut, clustered like the trailing edge
                edges.push(EdgeSpec {
                    vertices: [at(1, layer), at(9, layer)],
                    curve: Curve::Line,
                    distribution: Distribution::Geometric {
                        spacing: te,
                        at: EdgeEnd::Start,
                    },
                });
            }
            // the C's straight line on this side starts with the cell length
            // that the arc ends with, so the cells match across the vertical
            // line above (below) the trailing edge
            edges.push(EdgeSpec {
                vertices: [at(c_end, layer), at(if s == 0 { 10 } else { 11 }, layer)],
                curve: Curve::Line,
                distribution: Distribution::Geometric {
                    spacing: te * lr / rear,
                    at: EdgeEnd::Start,
                },
            });
        }
        if s == 0 {
            // i along the surface from the LE, j away from the wall
            blocks.push(BlockSpec {
                vertices: quad([0, 2, 7, 4]),
                cells: [split[0], n_normal, n_span],
                domain: fluid(),
            });
            blocks.push(BlockSpec {
                vertices: quad([2, 1, 5, 7]),
                cells: [split[1], n_normal, n_span],
                domain: fluid(),
            });
        } else {
            // i towards the LE, j away from the wall
            blocks.push(BlockSpec {
                vertices: quad([1, 3, 8, 6]),
                cells: [split[1], n_normal, n_span],
                domain: fluid(),
            });
            blocks.push(BlockSpec {
                vertices: quad([3, 0, 4, 8]),
                cells: [split[0], n_normal, n_span],
                domain: fluid(),
            });
        }
    }
    blocks.push(BlockSpec {
        vertices: quad([1, 9, 10, 5]),
        cells: [n_wake, n_normal, n_span],
        domain: fluid(),
    });
    blocks.push(BlockSpec {
        vertices: quad([9, 1, 6, 11]),
        cells: [n_wake, n_normal, n_span],
        domain: fluid(),
    });
    let side_face = |a: usize, b: usize| [at(a, 0), at(b, 0), at(b, 1), at(a, 1)];
    let span_faces = |layer: usize| -> Vec<[usize; 4]> {
        blocks
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
            .collect()
    };
    let patches = vec![
        PatchSpec {
            name: "airfoil_wall".into(),
            faces: vec![
                side_face(0, 2),
                side_face(2, 1),
                side_face(1, 3),
                side_face(3, 0),
            ],
        },
        PatchSpec {
            name: "inlet".into(),
            faces: vec![
                side_face(4, 7),
                side_face(7, 5),
                side_face(5, 10),
                side_face(6, 8),
                side_face(8, 4),
                side_face(11, 6),
            ],
        },
        PatchSpec {
            name: "outlet".into(),
            faces: vec![side_face(9, 10), side_face(11, 9)],
        },
        PatchSpec {
            name: "symmetry_low".into(),
            faces: span_faces(0),
        },
        PatchSpec {
            name: "symmetry_high".into(),
            faces: span_faces(1),
        },
    ];
    let area = std::f64::consts::PI * big * big / 2.0 + 2.0 * big * wake - c * c * section.area;
    Ok(Template {
        topology: BlockTopology {
            vertices,
            edges,
            blocks,
            patches,
            surfaces: Vec::new(),
        },
        analytic_volume: area * span,
    })
}
