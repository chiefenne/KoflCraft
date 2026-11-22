//! ALG-16 benchmark of the geometry-derived size sources of GRD-04
//! (research/grd-04-geometry-sizing.md).
//!
//! Curvature: three estimators of the largest principal curvature per facet
//! against the analytic value on spheres, an ellipsoid and CAD-style
//! cylinders:
//! - `normal-cycle`: the edge-based normal-cycle tensor (Cohen-Steiner &
//!   Morvan 2003) over the smooth vertex ring (the method chosen, in
//!   `mesh_algorithms::sizing::surface`);
//! - `cotan`: the cotangent mean curvature and angle-defect Gaussian
//!   curvature per vertex (Meyer, Desbrun, Schröder & Barr 2003),
//!   κ_max = H + sqrt(H² − K), averaged over the facet's corners; written
//!   here from the paper for the benchmark only;
//! - `normal-deviation`: the classical turn of the facet normal over the
//!   centroid distance to each edge neighbour, the largest; benchmark only.
//!
//! Gaps: the medial-ball thickness (graphics; shrinking-ball method of Ma,
//! Bae & Choi 2012) against the inward-normal ray cast (classical), each
//! meshed by the production path and measured by counting the cells crossed
//! along the wall normals at fixed wall samples.
//!
//! Run: `cargo run --release -p mesh-engine --example geometry_sizing_benchmark`

#[path = "../../../test-fixtures/shapes.rs"]
mod shapes;

use std::sync::Arc;
use std::time::Instant;

use mesh_algorithms::neighbours::KdTree;
use mesh_algorithms::sizing::surface::{
    ConvexPlanes, GapOptions, SurfaceGeometry, Thickness, WeightedSet, gap_primitives,
    normal_cycle_curvature,
};
use mesh_algorithms::sizing::{
    Frame, SeedOptions, Shape, SizeField, SizeSource, SourceKind, seed_sites,
};
use mesh_core::geometry::{PatchId, Vec3};
use mesh_engine::{StlMeshConfig, build_stl_mesh, check_mesh, generate_surface_mesh};
use mesh_geometry::{ConvexSurface, SurfaceTriangle};

const FEATURE: f64 = 30.0;

fn main() {
    curvature_table();
    gap_table();
    mesh_table();
}

// ---------------------------------------------------------------- curvature

/// The analytic largest principal curvature at a facet centroid.
type Truth = Box<dyn Fn(Vec3) -> f64>;

fn curvature_table() {
    println!("## Curvature estimators: |κ_est/κ − 1| over curved facets\n");
    println!("| Surface | Facets | Estimator | median | p90 | max | max κ·R on flat facets |");
    println!("|---|---|---|---|---|---|---|");
    let o = Vec3::ZERO;
    let axes = Vec3::new(0.5, 0.35, 0.25);
    let cases: Vec<(&str, shapes::Shape, Truth)> = vec![
        (
            "icosphere R=0.5",
            shapes::icosphere(o, 0.5, 3),
            Box::new(|_| 2.0),
        ),
        (
            "UV sphere R=0.5, 24×48",
            shapes::uv_sphere(o, 0.5, 24, 48),
            Box::new(|_| 2.0),
        ),
        (
            "ellipsoid 0.5/0.35/0.25",
            shapes::ellipsoid(o, axes, 3),
            Box::new(move |c| ellipsoid_kmax(c, axes)),
        ),
        (
            "strip cylinder R=0.25, n=48",
            shapes::strip_cylinder(0.0, 0.0, 0.25, 1.0, 48),
            Box::new(|c: Vec3| {
                if c.z > 1e-9 && c.z < 1.0 - 1e-9 {
                    4.0
                } else {
                    0.0
                }
            }),
        ),
        (
            "strip cylinder R=0.25, n=16",
            shapes::strip_cylinder(0.0, 0.0, 0.25, 1.0, 16),
            Box::new(|c: Vec3| {
                if c.z > 1e-9 && c.z < 1.0 - 1e-9 {
                    4.0
                } else {
                    0.0
                }
            }),
        ),
    ];
    for (name, (v, f), truth) in cases {
        let g = SurfaceGeometry::from_triangles(v, f);
        let feature = FEATURE.to_radians();
        let estimators: [(&str, Vec<f64>); 3] = [
            ("normal-cycle", normal_cycle_curvature(&g, feature)),
            ("cotan", cotan(&g)),
            ("normal-deviation", normal_deviation(&g, feature)),
        ];
        for (label, k) in estimators {
            let mut errors = Vec::new();
            let mut flat: f64 = 0.0;
            for (facet, &value) in k.iter().enumerate() {
                let c = g.centroid(facet);
                // The analytic curvature where the side facets' centroids
                // lie (cylinder sides: z strictly inside; caps: 0).
                let exact = truth(c);
                if exact > 0.0 {
                    errors.push((value / exact - 1.0).abs());
                } else {
                    flat = flat.max(value * 0.25);
                }
            }
            errors.sort_by(f64::total_cmp);
            let q = |p: f64| errors[((errors.len() - 1) as f64 * p).round() as usize];
            println!(
                "| {name} | {} | {label} | {:.3} | {:.3} | {:.3} | {:.2} |",
                g.triangles.len(),
                q(0.5),
                q(0.9),
                errors[errors.len() - 1],
                flat
            );
        }
    }
    println!();
}

/// κ_max of the ellipsoid x²/a² + y²/b² + z²/c² = 1 at the radial
/// projection of `p`: H and K in closed form, κ_max = H + sqrt(H² − K).
fn ellipsoid_kmax(p: Vec3, a: Vec3) -> f64 {
    let s = (p.x * p.x / (a.x * a.x) + p.y * p.y / (a.y * a.y) + p.z * p.z / (a.z * a.z)).sqrt();
    let q = p / s;
    let (a2, b2, c2) = (a.x * a.x, a.y * a.y, a.z * a.z);
    let w = q.x * q.x / (a2 * a2) + q.y * q.y / (b2 * b2) + q.z * q.z / (c2 * c2);
    let k = 1.0 / (a2 * b2 * c2 * w * w);
    let h = (a2 + b2 + c2 - q.dot(q)) / (2.0 * a.x * a.y * a.z * a.x * a.y * a.z * w.powf(1.5));
    h + (h * h - k).max(0.0).sqrt()
}

/// Meyer et al. (2003): per vertex, the cotangent Laplacian of the position
/// gives the mean curvature normal, the angle defect the Gaussian
/// curvature, both over the mixed Voronoi area.
fn cotan(g: &SurfaceGeometry) -> Vec<f64> {
    let n = g.vertices.len();
    let mut laplace = vec![Vec3::ZERO; n];
    let mut area = vec![0.0; n];
    let mut angle = vec![0.0; n];
    for t in &g.triangles {
        let p = t.map(|i| g.vertices[i]);
        let mut corner = [0.0; 3];
        for k in 0..3 {
            let (a, b) = (p[(k + 1) % 3] - p[k], p[(k + 2) % 3] - p[k]);
            corner[k] = a.cross(b).length().atan2(a.dot(b));
        }
        let tri_area = (p[1] - p[0]).cross(p[2] - p[0]).length() / 2.0;
        let obtuse = corner.iter().position(|&c| c > std::f64::consts::FRAC_PI_2);
        for k in 0..3 {
            let (i, j, l) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
            angle[i] += corner[k];
            // The edge (i, j) is opposite corner l; (i, l) opposite j.
            let cot = |c: f64| c.cos() / c.sin();
            let (cl, cj) = (cot(corner[(k + 2) % 3]), cot(corner[(k + 1) % 3]));
            laplace[i] = laplace[i]
                + (g.vertices[i] - g.vertices[j]) * cl
                + (g.vertices[i] - g.vertices[l]) * cj;
            area[i] += match obtuse {
                None => {
                    ((g.vertices[i] - g.vertices[j]).dot(g.vertices[i] - g.vertices[j]) * cl
                        + (g.vertices[i] - g.vertices[l]).dot(g.vertices[i] - g.vertices[l]) * cj)
                        / 8.0
                }
                Some(o) if o == k => tri_area / 2.0,
                Some(_) => tri_area / 4.0,
            };
        }
    }
    let vertex: Vec<f64> = (0..n)
        .map(|i| {
            let h = laplace[i].length() / (4.0 * area[i]);
            let k = (std::f64::consts::TAU - angle[i]) / area[i];
            h + (h * h - k).max(0.0).sqrt()
        })
        .collect();
    g.triangles
        .iter()
        .map(|t| t.iter().map(|&i| vertex[i]).sum::<f64>() / 3.0)
        .collect()
}

/// The classical estimate: the facet normal's turn over the centroid
/// distance, the largest over the edge neighbours (features excluded).
fn normal_deviation(g: &SurfaceGeometry, feature: f64) -> Vec<f64> {
    (0..g.triangles.len())
        .map(|f| {
            g.facet_edges[f]
                .iter()
                .map(|&e| {
                    let edge = &g.edges[e];
                    let other = if edge.facets[0] == f {
                        edge.facets[1]
                    } else {
                        edge.facets[0]
                    };
                    if edge.angle > feature {
                        0.0
                    } else {
                        edge.angle / (g.centroid(f) - g.centroid(other)).length()
                    }
                })
                .fold(0.0, f64::max)
        })
        .collect()
}

// --------------------------------------------------------------------- gaps

fn surface((v, f): &shapes::Shape) -> ConvexSurface {
    let triangles = f
        .iter()
        .map(|&vertices| SurfaceTriangle {
            vertices,
            patch: PatchId(0),
        })
        .collect();
    ConvexSurface::new(v.clone(), triangles, vec!["surface".into()]).expect("convex shape")
}

fn gap_table() {
    const ACROSS: f64 = 4.0;
    let (min_size, max_size, gradation) = (0.008, 0.15, 1.3);
    println!(
        "## Gaps: {ACROSS} cells across asked, min_size {min_size}, max_size {max_size}, gradation {gradation}\n"
    );
    println!(
        "| Shape | Thickness | Valid | Seeds | Samples | Cells across min / p10 / median | Below {ACROSS} | Time (s) |"
    );
    println!("|---|---|---|---|---|---|---|---|");
    // (name, shape, wedge_angle in degrees)
    let cases: [(&str, shapes::Shape, f64); 5] = [
        (
            "slab 1×1×0.1",
            shapes::cuboid(Vec3::new(1.0, 1.0, 0.1)),
            45.0,
        ),
        (
            "wedge 30°, 1 long",
            shapes::wedge(30f64.to_radians(), 1.0, 0.4),
            45.0,
        ),
        (
            "wedge 50°, 0.6 long, wedge_angle 60°",
            shapes::wedge(50f64.to_radians(), 0.6, 0.4),
            60.0,
        ),
        (
            "disc R=0.5, 0.08 thick, n=48",
            shapes::strip_cylinder(0.5, 0.5, 0.5, 0.08, 48),
            45.0,
        ),
        (
            "icosphere R=0.2",
            shapes::icosphere(Vec3::new(0.5, 0.5, 0.5), 0.2, 2),
            45.0,
        ),
    ];
    for (name, shape, wedge) in &cases {
        let s = surface(shape);
        let g = SurfaceGeometry::new(&s);
        let planes = ConvexPlanes::new(&g);
        for (label, method) in [
            ("medial ball", Thickness::MedialBall),
            ("normal ray", Thickness::NormalRay),
        ] {
            let started = Instant::now();
            let items = gap_primitives(
                &g,
                &planes,
                GapOptions {
                    cells_across: ACROSS,
                    wedge_angle: (*wedge).to_radians(),
                    method,
                    min_size,
                    max_size,
                },
            );
            let set = WeightedSet::new(items);
            let field = SizeField::new(
                vec![SizeSource {
                    name: "gap".into(),
                    kind: SourceKind::Gap,
                    size: set.smallest().unwrap_or(max_size),
                    shape: Shape::Weighted(Arc::new(set)),
                    frame: Frame::GLOBAL,
                    distance: 0.0,
                    growth: None,
                }],
                min_size,
                max_size,
                gradation,
            )
            .expect("field");
            let seeds = match seed_sites(
                &s,
                &field,
                SeedOptions {
                    seed: 1,
                    relaxation_iterations: 20,
                    max_seeds: 200_000,
                },
            ) {
                Ok(set) => set,
                Err(error) => {
                    println!("| {name} | {label} | seeding failed: {error} | | | | | |");
                    continue;
                }
            };
            let valid = generate_surface_mesh(&s, &seeds.sites).is_ok();
            let counts = cells_across(&g, &planes, &s, &seeds.sites, ACROSS, min_size);
            let below = counts.iter().filter(|&&c| (c as f64) < ACROSS).count();
            let mut sorted = counts.clone();
            sorted.sort_unstable();
            let q = |p: f64| sorted[((sorted.len() - 1) as f64 * p).round() as usize];
            println!(
                "| {name} | {label} | {} | {} | {} | {} / {} / {} | {} ({:.1} %) | {:.2} |",
                if valid { "yes" } else { "NO" },
                seeds.sites.len(),
                counts.len(),
                q(0.0),
                q(0.1),
                q(0.5),
                below,
                100.0 * below as f64 / counts.len() as f64,
                started.elapsed().as_secs_f64()
            );
        }
    }
    println!();
}

/// Cells crossed along the inward normal from a fixed grid of wall samples
/// (barycentric, 8 per edge) of every facet, where the path is long enough
/// that `across` cells of at least `min_size` fit.
fn cells_across(
    g: &SurfaceGeometry,
    planes: &ConvexPlanes,
    s: &ConvexSurface,
    sites: &[Vec3],
    across: f64,
    min_size: f64,
) -> Vec<usize> {
    let tree = KdTree::new(sites);
    let mut out = Vec::new();
    const GRID: usize = 8;
    for f in 0..g.triangles.len() {
        let [a, b, c] = g.corners(f);
        let n = g.normals[f];
        for i in 1..GRID {
            for j in 1..GRID - i {
                let start =
                    a + (b - a) * (i as f64 / GRID as f64) + (c - a) * (j as f64 / GRID as f64);
                let Some(length) =
                    planes.thickness(start, n, std::f64::consts::FRAC_PI_2, Thickness::NormalRay)
                else {
                    continue;
                };
                if length / across < min_size {
                    continue;
                }
                let mut cells = 0;
                let mut last = None;
                for k in 0..512 {
                    let p = start - n * (length * (k as f64 + 0.5) / 512.0);
                    if !s.contains_interior(p) {
                        continue;
                    }
                    let cell = tree.closest(p);
                    if cell != last {
                        cells += 1;
                        last = cell;
                    }
                }
                if cells > 0 {
                    out.push(cells);
                }
            }
        }
    }
    out
}

// ------------------------------------------------------------------- meshes

fn mesh_table() {
    println!("## Production builds (curvature, feature edges, all three)\n");
    println!(
        "| Case | Valid | Seeds | Source | Target | Cells near | Achieved median (p10–p90) | Time (s) |"
    );
    println!("|---|---|---|---|---|---|---|---|");
    let cylinder = shapes::ascii_stl(&shapes::strip_cylinder(0.5, 0.5, 0.25, 0.6, 48));
    let cube = shapes::ascii_stl(&shapes::cuboid(Vec3::new(1.0, 1.0, 1.0)));
    let wedge = shapes::ascii_stl(&shapes::wedge(30f64.to_radians(), 1.0, 0.4));
    let cases = [
        (
            "cylinder, curvature 15°",
            &cylinder,
            serde_json::json!([{"kind":"curvature","name":"wall","angle":15.0}]),
        ),
        (
            "unit box, feature edges 0.04",
            &cube,
            serde_json::json!([{"kind":"feature_edges","name":"edges","size":0.04}]),
        ),
        (
            "wedge 30°, all three",
            &wedge,
            serde_json::json!([
                {"kind":"curvature","name":"wall","angle":15.0},
                {"kind":"gap","name":"gap","cells_across":4.0},
                {"kind":"feature_edges","name":"edges","size":0.03}
            ]),
        ),
    ];
    for (name, stl, sources) in cases {
        let config: StlMeshConfig = serde_json::from_value(serde_json::json!({
            "schema_version": 1,
            "size_field": {"min_size": 0.01, "max_size": 0.15, "gradation": 1.3, "sources": sources},
            "seed": 1
        }))
        .expect("config");
        let started = Instant::now();
        let artifact = match build_stl_mesh(stl, config, None, None) {
            Ok(a) => a,
            Err(error) => {
                println!("| {name} | NO: {error} | | | | | | |");
                continue;
            }
        };
        let valid = check_mesh(&artifact).is_ok();
        let elapsed = started.elapsed().as_secs_f64();
        let sizing = artifact.sizing().expect("sizing");
        for source in &sizing.sources {
            let achieved = source
                .achieved
                .as_ref()
                .map(|a| format!("{:.4} ({:.4}–{:.4})", a.median, a.p10, a.p90))
                .unwrap_or_else(|| "–".into());
            println!(
                "| {name} | {} | {} | {} | {:.4} | {} | {achieved} | {elapsed:.2} |",
                if valid { "yes" } else { "NO" },
                sizing.seeds,
                source.name,
                source.target_size,
                source.cells,
            );
            if let Some(d) = &source.derived
                && let Some(c) = &d.cells_across
            {
                println!(
                    "| | | | {} cells across | | {} samples | min {} / median {} ; {} below |",
                    source.name,
                    d.samples.unwrap_or(0),
                    c.min,
                    c.median,
                    d.samples_below_target.unwrap_or(0)
                );
            }
        }
    }
}
