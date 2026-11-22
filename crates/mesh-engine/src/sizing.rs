//! The size field as a typed, named input of the STL mesher (GRD-01 to
//! GRD-03, GRD-05, GRD-08) and the achieved sizes it is reported with.

use std::sync::Arc;

use mesh_algorithms::neighbours::KdTree;
use mesh_algorithms::sizing::surface::{
    ConvexPlanes, GapOptions, Primitive, SurfaceGeometry, Thickness, WeightedSet,
    curvature_primitives, feature_primitives, gap_primitives,
};
use mesh_algorithms::sizing::{
    Frame, SeedOptions, SeedSet, Shape, SizeField, SizeSource, SizingError, SourceKind, seed_sites,
    target_count,
};
use mesh_core::geometry::Vec3;
use mesh_core::mesh::{FaceAdjacency, Mesh};
use mesh_geometry::ConvexSurface;
use mesh_quality::MeshMetrics;
use serde::{Deserialize, Serialize};

use crate::artifact::{Point, vector};
use crate::{EngineError, ErrorCode};

/// Lloyd iterations unless the config says otherwise (research/m2-size-control.md).
pub const DEFAULT_RELAXATION_ITERATIONS: usize = 20;
/// The largest seed count a size field may ask for while clipping is all-pairs (ALG-04).
pub const MAX_SIZED_SEEDS: usize = 200_000;
/// Edges whose facet normals turn more than this are feature edges, not
/// curvature (degrees, GRD-04).
pub const DEFAULT_FEATURE_ANGLE: f64 = 30.0;
/// Walls meeting at an interior angle below this form a gap (degrees, GRD-04).
pub const DEFAULT_WEDGE_ANGLE: f64 = 45.0;
/// Wall samples at most at which a gap source's cells across are counted.
const CELLS_ACROSS_SAMPLES: usize = 2_000;

#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SizeFieldConfig {
    /// The smallest cell size anywhere.
    pub min_size: f64,
    /// The largest cell size, also the size where no source applies.
    pub max_size: f64,
    /// The largest size ratio of neighbouring cells the field allows (> 1);
    /// it bounds the field's slope by `gradation - 1`.
    pub gradation: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub relaxation_iterations: Option<usize>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub frames: Vec<FrameConfig>,
    #[serde(default)]
    pub sources: Vec<SourceConfig>,
}

/// A named coordinate system: `z = x_axis × y_axis`, with `y_axis` made
/// normal to `x_axis`.
#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct FrameConfig {
    pub name: String,
    pub origin: Point,
    pub x_axis: Point,
    pub y_axis: Point,
}

#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum SourceConfig {
    /// GRD-02: a size near an object, growing away from it.
    Datum {
        name: String,
        shape: ShapeConfig,
        size: f64,
        #[serde(default)]
        distance: f64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        growth: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        frame: Option<String>,
    },
    /// GRD-03: a size inside coordinate bounds (box, slab or cylinder).
    Region {
        name: String,
        shape: ShapeConfig,
        size: f64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        frame: Option<String>,
    },
    /// GRD-04: near curved walls a cell spans at most `angle` degrees of
    /// the wall's normal turn (size = angle / largest principal curvature).
    Curvature {
        name: String,
        angle: f64,
        #[serde(default)]
        distance: f64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        growth: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        feature_angle: Option<f64>,
    },
    /// GRD-04: at least `cells_across` cells across the gap between
    /// opposite walls (size = local thickness / cells_across).
    Gap {
        name: String,
        cells_across: f64,
        #[serde(default)]
        distance: f64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        growth: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        wedge_angle: Option<f64>,
    },
    /// GRD-04: a size along the sharp edges of the surface.
    FeatureEdges {
        name: String,
        size: f64,
        #[serde(default)]
        distance: f64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        growth: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        feature_angle: Option<f64>,
    },
}

impl SourceConfig {
    fn derived(&self) -> bool {
        matches!(
            self,
            Self::Curvature { .. } | Self::Gap { .. } | Self::FeatureEdges { .. }
        )
    }
}

/// An angle in degrees within `(0, max]`, in radians.
fn angle(name: &str, what: &str, degrees: f64, max: f64) -> Result<f64, EngineError> {
    if degrees.is_finite() && degrees > 0.0 && degrees <= max {
        Ok(degrees.to_radians())
    } else {
        Err(EngineError::new(
            ErrorCode::InvalidSizeField,
            format!("size source {name}: {what} must be above 0 and at most {max} degrees"),
        ))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Axis {
    X,
    Y,
    Z,
}

#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum ShapeConfig {
    Point {
        point: Point,
    },
    Segment {
        start: Point,
        end: Point,
    },
    Plane {
        point: Point,
        normal: Point,
    },
    Box {
        min: Point,
        max: Point,
    },
    Slab {
        axis: Axis,
        min: f64,
        max: f64,
    },
    Sphere {
        centre: Point,
        radius: f64,
    },
    Cylinder {
        start: Point,
        end: Point,
        radius: f64,
    },
    Cone {
        start: Point,
        end: Point,
        start_radius: f64,
        end_radius: f64,
    },
}

impl ShapeConfig {
    fn shape(&self) -> Shape {
        match *self {
            Self::Point { point } => Shape::Point(vector(point)),
            Self::Segment { start, end } => Shape::Segment {
                start: vector(start),
                end: vector(end),
            },
            Self::Plane { point, normal } => Shape::Plane {
                point: vector(point),
                normal: vector(normal),
            },
            Self::Box { min, max } => Shape::Box {
                min: vector(min),
                max: vector(max),
            },
            Self::Slab { axis, min, max } => Shape::Slab {
                axis: axis as usize,
                min,
                max,
            },
            Self::Sphere { centre, radius } => Shape::Sphere {
                centre: vector(centre),
                radius,
            },
            Self::Cylinder { start, end, radius } => Shape::Frustum {
                start: vector(start),
                end: vector(end),
                start_radius: radius,
                end_radius: radius,
            },
            Self::Cone {
                start,
                end,
                start_radius,
                end_radius,
            } => Shape::Frustum {
                start: vector(start),
                end: vector(end),
                start_radius,
                end_radius,
            },
        }
    }
}

impl SizeFieldConfig {
    /// The field on `surface`; curvature, gap and feature-edge sources are
    /// derived from its facets (GRD-04).
    pub fn field(&self, surface: &ConvexSurface) -> Result<SizeField, EngineError> {
        let invalid = |message: String| EngineError::new(ErrorCode::InvalidSizeField, message);
        let mut frames: Vec<(&str, Frame)> = vec![("global", Frame::GLOBAL)];
        for frame in &self.frames {
            if frames.iter().any(|(name, _)| *name == frame.name) {
                return Err(invalid(format!("frame {} is defined twice", frame.name)));
            }
            let built = Frame::new(
                vector(frame.origin),
                vector(frame.x_axis),
                vector(frame.y_axis),
            )
            .ok_or_else(|| {
                invalid(format!(
                    "frame {} needs finite, nonparallel axes",
                    frame.name
                ))
            })?;
            frames.push((&frame.name, built));
        }
        let find = |name: &Option<String>| {
            let name = name.as_deref().unwrap_or("global");
            frames
                .iter()
                .find(|(n, _)| *n == name)
                .map(|(_, f)| *f)
                .ok_or_else(|| invalid(format!("unknown frame {name}")))
        };
        let geometry = self
            .sources
            .iter()
            .any(SourceConfig::derived)
            .then(|| SurfaceGeometry::new(surface));
        let (min, max) = (self.min_size, self.max_size);
        let derived = |name: &String, kind, items: Vec<_>, distance: f64, growth: Option<f64>| {
            let set = WeightedSet::new(items);
            SizeSource {
                name: name.clone(),
                kind,
                size: set.smallest().unwrap_or(max),
                shape: Shape::Weighted(Arc::new(set)),
                frame: Frame::GLOBAL,
                distance,
                growth,
            }
        };
        let mut sources = Vec::with_capacity(self.sources.len());
        for source in &self.sources {
            sources.push(match source {
                SourceConfig::Datum {
                    name,
                    shape,
                    size,
                    distance,
                    growth,
                    frame,
                } => SizeSource {
                    name: name.clone(),
                    kind: SourceKind::Datum,
                    shape: shape.shape(),
                    frame: find(frame)?,
                    size: *size,
                    distance: *distance,
                    growth: *growth,
                },
                SourceConfig::Region {
                    name,
                    shape,
                    size,
                    frame,
                } => {
                    if !matches!(
                        shape,
                        ShapeConfig::Box { .. }
                            | ShapeConfig::Slab { .. }
                            | ShapeConfig::Cylinder { .. }
                    ) {
                        return Err(invalid(format!(
                            "region {name} must be a box, slab or cylinder; use a datum for other shapes"
                        )));
                    }
                    SizeSource {
                        name: name.clone(),
                        kind: SourceKind::Region,
                        shape: shape.shape(),
                        frame: find(frame)?,
                        size: *size,
                        distance: 0.0,
                        growth: None,
                    }
                }
                SourceConfig::Curvature {
                    name,
                    angle: turn,
                    distance,
                    growth,
                    feature_angle,
                } => {
                    let turn = angle(name, "angle", *turn, 90.0)?;
                    let feature = angle(
                        name,
                        "feature_angle",
                        feature_angle.unwrap_or(DEFAULT_FEATURE_ANGLE),
                        180.0,
                    )?;
                    let items = geometry
                        .as_ref()
                        .map(|g| curvature_primitives(g, turn, feature, min, max))
                        .unwrap_or_default();
                    derived(name, SourceKind::Curvature, items, *distance, *growth)
                }
                SourceConfig::Gap {
                    name,
                    cells_across,
                    distance,
                    growth,
                    wedge_angle,
                } => {
                    if !(cells_across.is_finite() && *cells_across >= 1.0) {
                        return Err(invalid(format!(
                            "size source {name}: cells_across must be at least 1"
                        )));
                    }
                    let wedge = angle(
                        name,
                        "wedge_angle",
                        wedge_angle.unwrap_or(DEFAULT_WEDGE_ANGLE),
                        90.0,
                    )?;
                    let items = geometry
                        .as_ref()
                        .map(|g| {
                            gap_primitives(
                                g,
                                &ConvexPlanes::new(g),
                                GapOptions {
                                    cells_across: *cells_across,
                                    wedge_angle: wedge,
                                    method: Thickness::MedialBall,
                                    min_size: min,
                                    max_size: max,
                                },
                            )
                        })
                        .unwrap_or_default();
                    derived(name, SourceKind::Gap, items, *distance, *growth)
                }
                SourceConfig::FeatureEdges {
                    name,
                    size,
                    distance,
                    growth,
                    feature_angle,
                } => {
                    if !(size.is_finite() && min <= *size && *size <= max) {
                        return Err(invalid(format!(
                            "size source {name}: size must lie between min_size and max_size"
                        )));
                    }
                    let feature = angle(
                        name,
                        "feature_angle",
                        feature_angle.unwrap_or(DEFAULT_FEATURE_ANGLE),
                        180.0,
                    )?;
                    let items = geometry
                        .as_ref()
                        .map(|g| feature_primitives(g, feature, *size))
                        .unwrap_or_default();
                    derived(name, SourceKind::FeatureEdges, items, *distance, *growth)
                }
            });
        }
        SizeField::new(sources, self.min_size, self.max_size, self.gradation)
            .map_err(|error| invalid(error.to_string()))
    }

    pub(crate) fn seeds(
        &self,
        surface: &ConvexSurface,
        seed: u64,
    ) -> Result<(SizeField, SeedSet), EngineError> {
        let field = self.field(surface)?;
        let set = seed_sites(
            surface,
            &field,
            SeedOptions {
                seed,
                relaxation_iterations: self.iterations(),
                max_seeds: MAX_SIZED_SEEDS,
            },
        )
        .map_err(|error| {
            let code = match error {
                SizingError::Voronoi(_) => ErrorCode::GenerationFailed,
                _ => ErrorCode::InvalidSizeField,
            };
            EngineError::new(code, error)
        })?;
        Ok((field, set))
    }

    pub(crate) fn count(&self, surface: &ConvexSurface) -> Result<usize, EngineError> {
        Ok(target_count(surface, &self.field(surface)?))
    }

    fn iterations(&self) -> usize {
        self.relaxation_iterations
            .unwrap_or(DEFAULT_RELAXATION_ITERATIONS)
    }
}

/// What the seeding did and the sizes the mesh achieved (GRD-08). Sizes are
/// equivalent-cube edges, cbrt(cell volume); a cell belongs to a source when
/// its centroid lies within the source's influence distance and the source
/// sets h there (for curvature, gap and feature-edge sources: within half a
/// cell size beyond the distance, since they lie on the wall).
#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SizingEvidence {
    pub method: String,
    pub size_field: SizeFieldConfig,
    pub seed: u64,
    pub seeds: usize,
    pub relaxation_iterations: usize,
    /// The largest seed move of the last Lloyd iteration, relative to h.
    pub last_relative_shift: Option<f64>,
    /// cbrt(V) / h(centroid) over all cells.
    pub size_ratio: Quantiles,
    /// Larger over smaller cbrt(V) over interior faces.
    pub neighbour_ratio: Quantiles,
    /// Interior faces whose neighbour ratio exceeds the field's gradation.
    pub faces_over_gradation: usize,
    pub interior_faces: usize,
    pub sources: Vec<SourceEvidence>,
}

#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Quantiles {
    pub min: f64,
    pub p10: f64,
    pub median: f64,
    pub p90: f64,
    pub p99: f64,
    pub max: f64,
}

#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SourceEvidence {
    pub name: String,
    pub kind: String,
    pub target_size: f64,
    /// Cells with their centroid within the influence distance (inside a
    /// region) where this source sets h.
    pub cells: usize,
    /// Achieved cbrt(V) of those cells; absent when no centroid is that close.
    pub achieved: Option<Quantiles>,
    /// Cells where this source sets h (the smallest of all sources there).
    pub governs: usize,
    /// For sources derived from the surface (GRD-04).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub derived: Option<DerivedEvidence>,
}

/// What a curvature, gap or feature-edge source found on the surface.
#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct DerivedEvidence {
    pub method: String,
    /// Facets, facet pieces or edges that ask for a size below `max_size`.
    pub primitives: usize,
    /// Of those, the ones held at `min_size`: there the source asks for
    /// less than `min_size` and its target (angle, cells across) is not met.
    pub at_min_size: usize,
    /// Gap sources: cells crossed along the inward wall normal to the
    /// opposite wall, counted at wall samples whose size is above `min_size`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cells_across: Option<Quantiles>,
    /// Gap sources: samples counted, and those with fewer cells than asked.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub samples: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub samples_below_target: Option<usize>,
}

pub(crate) fn evidence(
    config: &SizeFieldConfig,
    field: &SizeField,
    set: &SeedSet,
    seed: u64,
    surface: &ConvexSurface,
    mesh: &Mesh,
    metrics: &MeshMetrics,
) -> SizingEvidence {
    let mut across = Across::new(surface, &set.sites);
    let size = |c: usize| metrics.cells[c].volume.cbrt();
    let centroids: Vec<Vec3> = metrics.cells.iter().map(|c| c.centroid).collect();
    let ratios = (0..centroids.len())
        .map(|c| size(c) / field.size(centroids[c]))
        .collect();
    let mut neighbour = Vec::new();
    for face in 0..mesh.face_count() {
        if let Some(FaceAdjacency::Interior {
            owner,
            neighbour: other,
        }) = mesh.face_adjacency(face)
        {
            let (a, b) = (size(owner), size(other));
            neighbour.push(a.max(b) / a.min(b));
        }
    }
    let over = neighbour.iter().filter(|&&r| r > field.gradation()).count();
    let sources = field
        .sources()
        .iter()
        .enumerate()
        .map(|(index, source)| {
            let governed: Vec<bool> = centroids
                .iter()
                .map(|&c| {
                    let own = field.source_size(index, c);
                    own < field.max_size()
                        && (0..field.sources().len())
                            .all(|other| own <= field.source_size(other, c))
                })
                .collect();
            let governs = governed.iter().filter(|&&g| g).count();
            // A surface-derived source sits on the wall, where no centroid
            // lies: its cells are those within half a cell size of its reach.
            let derived = matches!(source.shape, Shape::Weighted(_));
            let near: Vec<usize> = (0..centroids.len())
                .filter(|&c| {
                    let slack = if derived { 0.5 * size(c) } else { 0.0 };
                    governed[c]
                        && field.source_distance(index, centroids[c]) <= source.distance + slack
                })
                .collect();
            let derived = match (&source.shape, &config.sources[index]) {
                (Shape::Weighted(items), kind) => {
                    Some(derived_evidence(items, kind, field, &mut across))
                }
                _ => None,
            };
            SourceEvidence {
                name: source.name.clone(),
                kind: source.kind.name().into(),
                target_size: source.size,
                cells: near.len(),
                achieved: quantiles(near.iter().map(|&c| size(c)).collect()),
                governs,
                derived,
            }
        })
        .collect();
    SizingEvidence {
        method: "density_weighted_cvt_v1".into(),
        size_field: config.clone(),
        seed,
        seeds: set.sites.len(),
        relaxation_iterations: set.iterations,
        last_relative_shift: set.last_shift.is_finite().then_some(set.last_shift),
        size_ratio: quantiles(ratios).unwrap_or(Quantiles::ONE),
        interior_faces: neighbour.len(),
        neighbour_ratio: quantiles(neighbour).unwrap_or(Quantiles::ONE),
        faces_over_gradation: over,
        sources,
    }
}

fn derived_evidence(
    items: &WeightedSet,
    kind: &SourceConfig,
    field: &SizeField,
    across: &mut Across,
) -> DerivedEvidence {
    let primitives = items.items().len();
    let at_min_size = items
        .items()
        .iter()
        .filter(|w| w.size <= field.min_size())
        .count();
    let mut evidence = DerivedEvidence {
        method: match kind {
            SourceConfig::Curvature { .. } => "normal_cycle_tensor_v1",
            SourceConfig::Gap { .. } => "medial_ball_thickness_v1",
            _ => "dihedral_angle_v1",
        }
        .into(),
        primitives,
        at_min_size,
        cells_across: None,
        samples: None,
        samples_below_target: None,
    };
    if let SourceConfig::Gap { cells_across, .. } = kind {
        let stride = primitives.div_ceil(CELLS_ACROSS_SAMPLES).max(1);
        let counts: Vec<f64> = items
            .items()
            .iter()
            .step_by(stride)
            .filter(|w| w.size > field.min_size())
            .filter_map(|w| match w.primitive {
                Primitive::Triangle(t) => across.count(t),
                Primitive::Segment(_) => None,
            })
            .map(|n| n as f64)
            .collect();
        evidence.samples = Some(counts.len());
        evidence.samples_below_target = Some(counts.iter().filter(|&&n| n < *cells_across).count());
        evidence.cells_across = quantiles(counts);
    }
    evidence
}

/// Counts the cells a segment across the domain crosses: in a clipped
/// Voronoi mesh the cell containing a point is its nearest site's cell.
struct Across {
    planes: ConvexPlanes,
    sites: KdTree,
    surface: ConvexSurface,
}

impl Across {
    fn new(surface: &ConvexSurface, sites: &[Vec3]) -> Self {
        Self {
            planes: ConvexPlanes::new(&SurfaceGeometry::new(surface)),
            sites: KdTree::new(sites),
            surface: surface.clone(),
        }
    }

    /// Cells crossed from the centroid of wall piece `t` along its inward
    /// normal to the wall where the ray leaves the domain.
    fn count(&self, t: [Vec3; 3]) -> Option<usize> {
        let normal = (t[1] - t[0]).cross(t[2] - t[0]);
        let n = normal / normal.length();
        let start = (t[0] + t[1] + t[2]) / 3.0;
        // Every wall the ray can leave through (interior angle below 90°).
        let length =
            self.planes
                .thickness(start, n, std::f64::consts::FRAC_PI_2, Thickness::NormalRay)?;
        const STEPS: usize = 512;
        let mut cells = 0;
        let mut last = None;
        for k in 0..STEPS {
            let p = start - n * (length * (k as f64 + 0.5) / STEPS as f64);
            if !self.surface.contains_interior(p) {
                continue;
            }
            let cell = self.sites.closest(p);
            if cell != last {
                cells += 1;
                last = cell;
            }
        }
        (cells > 0).then_some(cells)
    }
}

impl Quantiles {
    const ONE: Self = Self {
        min: 1.0,
        p10: 1.0,
        median: 1.0,
        p90: 1.0,
        p99: 1.0,
        max: 1.0,
    };
}

fn quantiles(mut values: Vec<f64>) -> Option<Quantiles> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(f64::total_cmp);
    let at = |q: f64| values[((values.len() - 1) as f64 * q).round() as usize];
    Some(Quantiles {
        min: values[0],
        p10: at(0.1),
        median: at(0.5),
        p90: at(0.9),
        p99: at(0.99),
        max: values[values.len() - 1],
    })
}
