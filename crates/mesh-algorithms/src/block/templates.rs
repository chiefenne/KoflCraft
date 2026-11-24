//! Parametric block topologies for standard cases (BS-03). A template only
//! writes a topology; the grid is built from it like from any other.
//! MB1a: the straight channel and the pipe O-grid with an H-core.

use mesh_core::block::{
    BlockSpec, BlockTopology, Curve, Distribution, EdgeEnd, EdgeSpec, PatchSpec,
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
    if !(core_fraction > 0.1 && core_fraction < 0.9) {
        return Err("core_fraction must lie between 0.1 and 0.9".into());
    }
    counts(&cells)?;
    let [core, radial, axial] = cells;
    let rc = core_fraction * r;
    let at = |x: f64, rad: f64, deg: f64| {
        let t = deg.to_radians();
        Vec3::new(x, rad * t.cos(), rad * t.sin())
    };
    // vertex 8 p + m: inner corner m at 45° + 90° m; 8 p + 4 + m: outer
    let mut vertices = Vec::with_capacity(16);
    for x in [0.0, l] {
        for m in 0..4 {
            vertices.push(at(x, rc, 45.0 + 90.0 * m as f64));
        }
        for m in 0..4 {
            vertices.push(at(x, r, 45.0 + 90.0 * m as f64));
        }
    }
    let c = |p: usize, m: usize| 8 * p + m % 4;
    let o = |p: usize, m: usize| 8 * p + 4 + m % 4;
    let mut blocks = vec![BlockSpec {
        // counter-clockwise in (y, z) from 225°, k along +x
        vertices: [
            c(0, 2),
            c(0, 3),
            c(0, 0),
            c(0, 1),
            c(1, 2),
            c(1, 3),
            c(1, 0),
            c(1, 1),
        ],
        cells: [core, core, axial],
        domain: "fluid".into(),
    }];
    let mut edges = Vec::new();
    let mut inlet = vec![[c(0, 2), c(0, 3), c(0, 0), c(0, 1)]];
    let mut outlet = vec![[c(1, 2), c(1, 3), c(1, 0), c(1, 1)]];
    let mut walls = Vec::new();
    for q in 0..4 {
        let m1 = q + 3; // the corner at 90° q − 45°
        blocks.push(BlockSpec {
            vertices: [
                c(0, m1),
                o(0, m1),
                o(0, q),
                c(0, q),
                c(1, m1),
                o(1, m1),
                o(1, q),
                c(1, q),
            ],
            cells: [radial, core, axial],
            domain: "fluid".into(),
        });
        for p in 0..2 {
            let x = if p == 0 { 0.0 } else { l };
            edges.push(EdgeSpec {
                vertices: [o(p, m1), o(p, q)],
                curve: Curve::Arc {
                    through: at(x, r, 90.0 * q as f64),
                },
                distribution: Distribution::Uniform,
            });
            if let Some(s) = wall_spacing {
                edges.push(EdgeSpec {
                    vertices: [c(p, q), o(p, q)],
                    curve: Curve::Line,
                    distribution: Distribution::Geometric {
                        spacing: positive(s, "wall spacing")?,
                        at: EdgeEnd::End,
                    },
                });
            }
        }
        inlet.push([c(0, m1), o(0, m1), o(0, q), c(0, q)]);
        outlet.push([c(1, m1), o(1, m1), o(1, q), c(1, q)]);
        walls.push([o(0, m1), o(0, q), o(1, q), o(1, m1)]);
    }
    Ok(Template {
        topology: BlockTopology {
            vertices,
            edges,
            blocks,
            patches: vec![
                PatchSpec {
                    name: "inlet".into(),
                    faces: inlet,
                },
                PatchSpec {
                    name: "outlet".into(),
                    faces: outlet,
                },
                PatchSpec {
                    name: "walls".into(),
                    faces: walls,
                },
            ],
        },
        analytic_volume: std::f64::consts::PI * r * r * l,
    })
}
