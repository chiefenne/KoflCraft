//! The scalar size field h(x) of GRD-01 to GRD-03 and GRD-05.
//!
//! Every source is a shape with a target size. A source's size holds within
//! its influence distance of the shape and grows linearly with the distance
//! beyond it. The field is the minimum over all sources, capped by the
//! maximum size and floored by the minimum size.
//!
//! Gradation (GRD-05) holds by construction instead of by smoothing: each
//! distance function is 1-Lipschitz, so a source that grows at rate `r` is
//! `r`-Lipschitz; the rate is capped at `gradation - 1`, and the minimum of
//! such functions and the clamps keep that bound. This is the closed form of
//! gradient limiting (Persson 2006) for sources given as distance fields:
//! the limited field of a source with size `h_i` on a set `S_i` is exactly
//! `h_i + β·dist(x, S_i)`, so no background grid is needed. Distances are
//! Euclidean, which equals the distance through the domain for the convex
//! domains of M2.

use std::fmt;
use std::sync::Arc;

use mesh_core::geometry::Vec3;

use super::surface::WeightedSet;

/// A right-handed orthonormal coordinate system (GRD-03).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Frame {
    origin: Vec3,
    axes: [Vec3; 3],
}

impl Frame {
    pub const GLOBAL: Self = Self {
        origin: Vec3::ZERO,
        axes: [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ],
    };

    /// The frame with `x_axis` and the part of `y_axis` normal to it; `z = x × y`.
    pub fn new(origin: Vec3, x_axis: Vec3, y_axis: Vec3) -> Option<Self> {
        if !origin.is_finite() || !x_axis.is_finite() || !y_axis.is_finite() {
            return None;
        }
        let x = unit(x_axis)?;
        let normal_part = y_axis - x * y_axis.dot(x);
        // Reject nearly parallel axes, which would make the frame ill-defined.
        if normal_part.length() <= 1e-9 * y_axis.length() {
            return None;
        }
        let y = unit(normal_part)?;
        Some(Self {
            origin,
            axes: [x, y, x.cross(y)],
        })
    }

    /// The point's coordinates in this frame.
    pub fn local(&self, point: Vec3) -> Vec3 {
        let d = point - self.origin;
        Vec3::new(
            d.dot(self.axes[0]),
            d.dot(self.axes[1]),
            d.dot(self.axes[2]),
        )
    }
}

fn unit(v: Vec3) -> Option<Vec3> {
    let length = v.length();
    (length > 0.0 && length.is_finite()).then(|| v / length)
}

/// The shape of a size source, in its frame's coordinates. Bodies are solid:
/// the distance is zero inside them.
#[derive(Clone, Debug, PartialEq)]
pub enum Shape {
    Point(Vec3),
    Segment {
        start: Vec3,
        end: Vec3,
    },
    /// Two-sided: the distance to the plane.
    Plane {
        point: Vec3,
        normal: Vec3,
    },
    /// Axis-aligned in the frame.
    Box {
        min: Vec3,
        max: Vec3,
    },
    /// All points with `min <= x[axis] <= max`.
    Slab {
        axis: usize,
        min: f64,
        max: f64,
    },
    Sphere {
        centre: Vec3,
        radius: f64,
    },
    /// A solid frustum from `start` to `end`; a cylinder has equal radii and a
    /// cone a zero radius at one end.
    Frustum {
        start: Vec3,
        end: Vec3,
        start_radius: f64,
        end_radius: f64,
    },
    /// Primitives with their own sizes and reaches, derived from the
    /// surface (GRD-04); the distance is the reach distance.
    Weighted(Arc<WeightedSet>),
}

impl Shape {
    fn check(&self) -> Result<(), &'static str> {
        let finite = |v: Vec3| v.is_finite();
        let ok = match *self {
            Self::Weighted(ref set) => set.is_valid(),
            Self::Point(p) => finite(p),
            Self::Segment { start, end } => {
                finite(start) && finite(end) && (end - start).length() > 0.0
            }
            Self::Plane { point, normal } => {
                finite(point) && finite(normal) && normal.length() > 0.0
            }
            Self::Box { min, max } => {
                finite(min) && finite(max) && min.x < max.x && min.y < max.y && min.z < max.z
            }
            Self::Slab { axis, min, max } => {
                axis < 3 && min.is_finite() && max.is_finite() && min < max
            }
            Self::Sphere { centre, radius } => finite(centre) && radius.is_finite() && radius > 0.0,
            Self::Frustum {
                start,
                end,
                start_radius,
                end_radius,
            } => {
                finite(start)
                    && finite(end)
                    && (end - start).length() > 0.0
                    && start_radius.is_finite()
                    && end_radius.is_finite()
                    && start_radius >= 0.0
                    && end_radius >= 0.0
                    && start_radius.max(end_radius) > 0.0
            }
        };
        if ok {
            Ok(())
        } else {
            Err("the shape has nonfinite, empty or degenerate dimensions")
        }
    }

    /// The Euclidean distance from `p` to the shape (zero inside a body).
    pub fn distance(&self, p: Vec3) -> f64 {
        match *self {
            Self::Weighted(ref set) => set.reach_distance(p),
            Self::Point(q) => (p - q).length(),
            Self::Segment { start, end } => segment_distance(p, start, end),
            Self::Plane { point, normal } => ((p - point).dot(normal) / normal.length()).abs(),
            Self::Box { min, max } => {
                let excess = |v: f64, lo: f64, hi: f64| (lo - v).max(v - hi).max(0.0);
                Vec3::new(
                    excess(p.x, min.x, max.x),
                    excess(p.y, min.y, max.y),
                    excess(p.z, min.z, max.z),
                )
                .length()
            }
            Self::Slab { axis, min, max } => {
                let v = [p.x, p.y, p.z][axis];
                (min - v).max(v - max).max(0.0)
            }
            Self::Sphere { centre, radius } => ((p - centre).length() - radius).max(0.0),
            Self::Frustum {
                start,
                end,
                start_radius,
                end_radius,
            } => {
                let axis = end - start;
                let length = axis.length();
                let along = (p - start).dot(axis) / length;
                let radial = ((p - start) - axis * (along / length)).length();
                frustum_distance(along, radial, length, start_radius, end_radius)
            }
        }
    }
}

fn segment_distance(p: Vec3, a: Vec3, b: Vec3) -> f64 {
    let ab = b - a;
    let t = ((p - a).dot(ab) / ab.dot(ab)).clamp(0.0, 1.0);
    (p - (a + ab * t)).length()
}

/// Distance in the meridian half-plane from `(t, rho)` to the trapezoid
/// `0 <= t <= length, 0 <= rho <= r(t)`; it equals the 3D distance to the
/// solid of revolution.
fn frustum_distance(t: f64, rho: f64, length: f64, ra: f64, rb: f64) -> f64 {
    if (0.0..=length).contains(&t) && rho <= ra + (rb - ra) * (t / length) {
        return 0.0;
    }
    let d2 = |px: f64, py: f64, ax: f64, ay: f64, bx: f64, by: f64| {
        let (ex, ey) = (bx - ax, by - ay);
        let s = (((px - ax) * ex + (py - ay) * ey) / (ex * ex + ey * ey)).clamp(0.0, 1.0);
        let (dx, dy) = (px - ax - s * ex, py - ay - s * ey);
        dx * dx + dy * dy
    };
    let mut best = d2(t, rho, 0.0, ra, length, rb);
    // A cap of zero radius has no area; the slant edge covers its point.
    if ra > 0.0 {
        best = best.min(d2(t, rho, 0.0, 0.0, 0.0, ra));
    }
    if rb > 0.0 {
        best = best.min(d2(t, rho, length, 0.0, length, rb));
    }
    best.sqrt()
}

/// Whether a source is a datum object (GRD-02), a coordinate region
/// (GRD-03) or derived from the surface (GRD-04).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SourceKind {
    Datum,
    Region,
    Curvature,
    Gap,
    FeatureEdges,
}

impl SourceKind {
    pub fn name(self) -> &'static str {
        match self {
            Self::Datum => "datum",
            Self::Region => "region",
            Self::Curvature => "curvature",
            Self::Gap => "gap",
            Self::FeatureEdges => "feature_edges",
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SizeSource {
    pub name: String,
    pub kind: SourceKind,
    pub shape: Shape,
    pub frame: Frame,
    /// Target size within `distance` of the shape; for weighted primitives
    /// the smallest of their sizes (`max_size` when there are none).
    pub size: f64,
    /// Influence distance; zero for regions.
    pub distance: f64,
    /// Size ratio per unit of size beyond the influence distance; `None`
    /// takes the gradation limit. Capped at the gradation limit.
    pub growth: Option<f64>,
}

/// The combined, gradation-limited size field.
#[derive(Clone, Debug, PartialEq)]
pub struct SizeField {
    sources: Vec<SizeSource>,
    rates: Vec<f64>,
    min_size: f64,
    max_size: f64,
    gradation: f64,
}

impl SizeField {
    /// `max_size` is also the background size where no source applies.
    /// `gradation` is the largest size ratio of neighbouring cells (> 1).
    pub fn new(
        sources: Vec<SizeSource>,
        min_size: f64,
        max_size: f64,
        gradation: f64,
    ) -> Result<Self, SizingError> {
        if !(min_size.is_finite() && min_size > 0.0) {
            return Err(SizingError::Field(
                "min_size must be positive and finite".into(),
            ));
        }
        if !(max_size.is_finite() && max_size >= min_size) {
            return Err(SizingError::Field(
                "max_size must be finite and at least min_size".into(),
            ));
        }
        if !(gradation.is_finite() && gradation > 1.0) {
            return Err(SizingError::Field(
                "gradation must be finite and greater than 1".into(),
            ));
        }
        let mut rates = Vec::with_capacity(sources.len());
        for (index, source) in sources.iter().enumerate() {
            let fail = |why: &str| SizingError::Source {
                name: source.name.clone(),
                message: why.to_owned(),
            };
            if source.name.is_empty() || source.name.chars().any(char::is_control) {
                return Err(fail(
                    "the name must be nonempty, without control characters",
                ));
            }
            if sources[..index].iter().any(|s| s.name == source.name) {
                return Err(fail("the name is used twice"));
            }
            source.shape.check().map_err(fail)?;
            if !(source.size.is_finite() && min_size <= source.size && source.size <= max_size) {
                return Err(fail("size must lie between min_size and max_size"));
            }
            if let Shape::Weighted(set) = &source.shape {
                let items = set.items();
                if items
                    .iter()
                    .any(|w| !(min_size <= w.size && w.size <= max_size))
                    || set.smallest().unwrap_or(max_size) != source.size
                {
                    return Err(fail(
                        "derived sizes must lie between min_size and max_size, the smallest being the size",
                    ));
                }
            }
            if !(source.distance.is_finite() && source.distance >= 0.0) {
                return Err(fail("distance must be finite and not negative"));
            }
            if source.kind == SourceKind::Region && source.distance != 0.0 {
                return Err(fail("a region has no influence distance"));
            }
            let growth = source.growth.unwrap_or(gradation);
            if !(growth.is_finite() && growth > 1.0) {
                return Err(fail("growth must be finite and greater than 1"));
            }
            rates.push((growth - 1.0).min(gradation - 1.0));
        }
        Ok(Self {
            sources,
            rates,
            min_size,
            max_size,
            gradation,
        })
    }

    pub fn sources(&self) -> &[SizeSource] {
        &self.sources
    }

    pub fn min_size(&self) -> f64 {
        self.min_size
    }

    pub fn max_size(&self) -> f64 {
        self.max_size
    }

    pub fn gradation(&self) -> f64 {
        self.gradation
    }

    /// The growth rate (|∇h| bound) source `index` uses after the cap.
    pub fn rate(&self, index: usize) -> f64 {
        self.rates[index]
    }

    /// The smallest size the field can take anywhere.
    pub fn lower_bound(&self) -> f64 {
        self.sources
            .iter()
            .map(|s| s.size)
            .fold(self.max_size, f64::min)
            .max(self.min_size)
    }

    /// The distance from `p` to source `index`'s shape.
    pub fn source_distance(&self, index: usize, p: Vec3) -> f64 {
        let source = &self.sources[index];
        source.shape.distance(source.frame.local(p))
    }

    /// The size source `index` alone asks for at `p`.
    pub fn source_size(&self, index: usize, p: Vec3) -> f64 {
        let source = &self.sources[index];
        if let Shape::Weighted(set) = &source.shape {
            return set.envelope(
                source.frame.local(p),
                self.rates[index],
                source.distance,
                f64::INFINITY,
            );
        }
        let beyond = (self.source_distance(index, p) - source.distance).max(0.0);
        source.size + self.rates[index] * beyond
    }

    /// h(p).
    pub fn size(&self, p: Vec3) -> f64 {
        let mut best = self.max_size;
        for (i, source) in self.sources.iter().enumerate() {
            let size = match &source.shape {
                // The current minimum bounds the search of the hierarchy.
                Shape::Weighted(set) => {
                    set.envelope(source.frame.local(p), self.rates[i], source.distance, best)
                }
                _ => self.source_size(i, p),
            };
            best = best.min(size);
        }
        best.clamp(self.min_size, self.max_size)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum SizingError {
    Field(String),
    Source {
        name: String,
        message: String,
    },
    /// Rejection sampling did not place the seeds within its budget.
    SamplingBudget {
        placed: usize,
        target: usize,
    },
    TooManySeeds {
        target: usize,
        limit: usize,
    },
    Voronoi(crate::voronoi::VoronoiError),
}

impl fmt::Display for SizingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Field(message) => write!(f, "size field: {message}"),
            Self::Source { name, message } => write!(f, "size source {name}: {message}"),
            Self::SamplingBudget { placed, target } => write!(
                f,
                "seed sampling placed {placed} of {target} seeds within its budget; \
                 the ratio of max_size to the smallest size is too large"
            ),
            Self::TooManySeeds { target, limit } => write!(
                f,
                "the size field asks for {target} seeds, more than the limit of {limit}"
            ),
            Self::Voronoi(error) => write!(f, "seed relaxation: {error}"),
        }
    }
}

impl std::error::Error for SizingError {}

#[cfg(test)]
mod tests {
    use super::*;

    fn datum(name: &str, shape: Shape, size: f64, distance: f64) -> SizeSource {
        SizeSource {
            name: name.into(),
            kind: SourceKind::Datum,
            shape,
            frame: Frame::GLOBAL,
            size,
            distance,
            growth: None,
        }
    }

    // Keep rule (SW-18): the exact distances of every shape. A sign or axis
    // error in one distance function breaks this; the seeding tests only see
    // the combined field.
    #[test]
    fn shape_distances_are_exact() {
        let p = Vec3::new(3.0, 4.0, 0.0);
        assert_eq!(Shape::Point(Vec3::ZERO).distance(p), 5.0);
        let seg = Shape::Segment {
            start: Vec3::ZERO,
            end: Vec3::new(1.0, 0.0, 0.0),
        };
        assert!((seg.distance(p) - (4.0f64 + 16.0).sqrt()).abs() < 1e-15);
        let plane = Shape::Plane {
            point: Vec3::ZERO,
            normal: Vec3::new(0.0, 2.0, 0.0),
        };
        assert_eq!(plane.distance(p), 4.0);
        assert_eq!(plane.distance(Vec3::new(0.0, -4.0, 0.0)), 4.0);
        let cube = Shape::Box {
            min: Vec3::ZERO,
            max: Vec3::new(1.0, 1.0, 1.0),
        };
        assert_eq!(cube.distance(Vec3::new(0.5, 0.5, 0.5)), 0.0);
        assert!((cube.distance(Vec3::new(4.0, 5.0, 0.5)) - 5.0).abs() < 1e-15);
        let slab = Shape::Slab {
            axis: 2,
            min: -1.0,
            max: 1.0,
        };
        assert_eq!(slab.distance(Vec3::new(9.0, 9.0, 3.0)), 2.0);
        let ball = Shape::Sphere {
            centre: Vec3::ZERO,
            radius: 1.0,
        };
        assert_eq!(ball.distance(p), 4.0);
        let cylinder = Shape::Frustum {
            start: Vec3::ZERO,
            end: Vec3::new(0.0, 0.0, 2.0),
            start_radius: 1.0,
            end_radius: 1.0,
        };
        assert_eq!(cylinder.distance(Vec3::new(0.5, 0.0, 1.0)), 0.0);
        assert!((cylinder.distance(Vec3::new(3.0, 0.0, 1.0)) - 2.0).abs() < 1e-15);
        assert!((cylinder.distance(Vec3::new(0.0, 0.0, -1.5)) - 1.5).abs() < 1e-15);
        assert!((cylinder.distance(Vec3::new(4.0, 0.0, 6.0)) - 5.0).abs() < 1e-15);
        // A cone of half-angle 45° with its tip at z = 1: the distance to the
        // slant from (2, 0, 0) is 2 / sqrt(2) beyond the base rim's normal.
        let cone = Shape::Frustum {
            start: Vec3::ZERO,
            end: Vec3::new(0.0, 0.0, 1.0),
            start_radius: 1.0,
            end_radius: 0.0,
        };
        assert_eq!(cone.distance(Vec3::new(0.2, 0.0, 0.5)), 0.0);
        assert!((cone.distance(Vec3::new(1.0, 0.0, 1.0)) - 0.5f64.sqrt()).abs() < 1e-15);
        assert!((cone.distance(Vec3::new(0.0, 0.0, 3.0)) - 2.0).abs() < 1e-15);
    }

    #[test]
    fn frames_move_regions() {
        let frame = Frame::new(
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
        )
        .unwrap();
        assert_eq!(
            frame.local(Vec3::new(1.0, 3.0, 1.0)),
            Vec3::new(2.0, 0.0, 0.0)
        );
        assert_eq!(
            frame.local(Vec3::new(1.0, 1.0, 4.0)),
            Vec3::new(0.0, 0.0, 3.0)
        );
        assert!(
            Frame::new(
                Vec3::ZERO,
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(2.0, 0.0, 0.0)
            )
            .is_none()
        );
    }

    // Keep rule (SW-18): the field is the clamped minimum and its slope never
    // exceeds gradation - 1 (GRD-05). Dropping the rate cap or the clamp
    // breaks it; no other test probes the field between seeds.
    #[test]
    fn the_field_is_a_clamped_gradation_limited_minimum() {
        let field = SizeField::new(
            vec![
                datum("tip", Shape::Point(Vec3::ZERO), 0.1, 0.2),
                SizeSource {
                    growth: Some(3.0),
                    ..datum(
                        "wall",
                        Shape::Slab {
                            axis: 0,
                            min: 2.0,
                            max: 3.0,
                        },
                        0.2,
                        0.0,
                    )
                },
            ],
            0.05,
            1.0,
            1.25,
        )
        .unwrap();
        assert_eq!(field.size(Vec3::new(0.1, 0.0, 0.0)), 0.1);
        assert!((field.size(Vec3::new(0.6, 0.0, 0.0)) - (0.1 + 0.25 * 0.4)).abs() < 1e-15);
        assert_eq!(field.size(Vec3::new(2.5, 0.0, 0.0)), 0.2);
        assert_eq!(field.size(Vec3::new(50.0, 0.0, 0.0)), 1.0);
        assert_eq!(field.rate(1), 0.25);
        assert_eq!(field.lower_bound(), 0.1);
        let mut previous = field.size(Vec3::new(0.0, 0.3, 0.1));
        let step = 1e-3;
        for i in 1..6000 {
            let size = field.size(Vec3::new(i as f64 * step, 0.3, 0.1));
            assert!((size - previous).abs() <= 0.25 * step * (1.0 + 1e-9));
            previous = size;
        }
    }

    #[test]
    fn invalid_fields_and_sources_are_rejected() {
        let good = datum("a", Shape::Point(Vec3::ZERO), 0.1, 0.0);
        assert!(SizeField::new(vec![], 0.0, 1.0, 1.2).is_err());
        assert!(SizeField::new(vec![], 0.5, 0.4, 1.2).is_err());
        assert!(SizeField::new(vec![], 0.1, 1.0, 1.0).is_err());
        assert!(SizeField::new(vec![good.clone(), good.clone()], 0.1, 1.0, 1.2).is_err());
        let small = SizeSource {
            size: 0.01,
            ..good.clone()
        };
        assert!(SizeField::new(vec![small], 0.1, 1.0, 1.2).is_err());
        let region = SizeSource {
            kind: SourceKind::Region,
            distance: 0.1,
            ..good.clone()
        };
        assert!(SizeField::new(vec![region], 0.1, 1.0, 1.2).is_err());
        let flat = SizeSource {
            shape: Shape::Sphere {
                centre: Vec3::ZERO,
                radius: 0.0,
            },
            ..good
        };
        assert!(SizeField::new(vec![flat], 0.1, 1.0, 1.2).is_err());
    }
}
