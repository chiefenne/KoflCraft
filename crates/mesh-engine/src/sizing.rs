//! The size field as a typed, named input of the STL mesher (GRD-01 to
//! GRD-03, GRD-05, GRD-08) and the achieved sizes it is reported with.

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
pub const MAX_SIZED_SEEDS: usize = 20_000;

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
    pub fn field(&self) -> Result<SizeField, EngineError> {
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
        let field = self.field()?;
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
        Ok(target_count(surface, &self.field()?))
    }

    fn iterations(&self) -> usize {
        self.relaxation_iterations
            .unwrap_or(DEFAULT_RELAXATION_ITERATIONS)
    }
}

/// What the seeding did and the sizes the mesh achieved (GRD-08). Sizes are
/// equivalent-cube edges, cbrt(cell volume); a cell belongs to a source when
/// its centroid lies within the source's influence distance and the source
/// sets h there.
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
}

pub(crate) fn evidence(
    config: &SizeFieldConfig,
    field: &SizeField,
    set: &SeedSet,
    seed: u64,
    mesh: &Mesh,
    metrics: &MeshMetrics,
) -> SizingEvidence {
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
            let near: Vec<usize> = (0..centroids.len())
                .filter(|&c| {
                    governed[c] && field.source_distance(index, centroids[c]) <= source.distance
                })
                .collect();
            SourceEvidence {
                name: source.name.clone(),
                kind: match source.kind {
                    SourceKind::Datum => "datum",
                    SourceKind::Region => "region",
                }
                .into(),
                target_size: source.size,
                cells: near.len(),
                achieved: quantiles(near.iter().map(|&c| size(c)).collect()),
                governs,
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
