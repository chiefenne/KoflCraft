//! ALG-16 benchmark of seed placement for the size field (research/m2-size-control.md).
//!
//! Three placements of the same seed count on the same size field:
//! - `sampled`: rejection sampling only (step 2 of mesh_algorithms::sizing);
//! - `springs`: the classical size-function method of DistMesh (Persson &
//!   Strang 2004): repulsive bar forces toward the target lengths, written
//!   here from the paper for the benchmark only;
//! - `cvt`: density-weighted Lloyd relaxation, ρ = h⁻⁵ (the method chosen).
//!
//! Each mesh is built and validated by the production path first (validity),
//! then measured: size ratio cbrt(V)/h(centroid), neighbour size ratio,
//! non-orthogonality, sphericity and time.
//!
//! Run: `cargo run --release -p mesh-engine --example sizing_benchmark`

use std::time::Instant;

use mesh_algorithms::clipping::FaceSource;
use mesh_algorithms::sizing::{
    Frame, SeedOptions, Shape, SizeField, SizeSource, SourceKind, sample, seed_sites, target_count,
};
use mesh_algorithms::voronoi::clipped_surface_voronoi;
use mesh_core::geometry::Vec3;
use mesh_core::mesh::FaceAdjacency;
use mesh_engine::generate_surface_mesh;
use mesh_geometry::{ConvexSurface, StlOptions, import_stl};

const ITERATIONS: usize = 20;
/// Lloyd iteration counts of the convergence sweep.
const SWEEP: [usize; 5] = [0, 5, 10, 20, 40];

fn main() {
    let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../test-fixtures/");
    let load = |name: &str| {
        let bytes = std::fs::read(format!("{root}{name}")).expect("fixture");
        import_stl(
            &bytes,
            StlOptions {
                facet_patches: None,
                facet_limit: None,
            },
        )
        .expect("surface")
    };
    let cube = load("unit-box.stl");
    let skewed = load("skewed-box.stl");
    let datum = |name: &str, shape, size, distance| SizeSource {
        name: name.into(),
        kind: SourceKind::Datum,
        shape,
        frame: Frame::GLOBAL,
        size,
        distance,
        growth: None,
    };
    let cases: Vec<(&str, &ConvexSurface, SizeField)> = vec![
        (
            "uniform h=0.125, unit box",
            &cube,
            SizeField::new(vec![], 0.1, 0.125, 1.2).unwrap(),
        ),
        (
            "point datum 0.04 → 0.25, unit box",
            &cube,
            SizeField::new(
                vec![datum(
                    "tip",
                    Shape::Point(Vec3::new(0.25, 0.25, 0.25)),
                    0.04,
                    0.05,
                )],
                0.04,
                0.25,
                1.3,
            )
            .unwrap(),
        ),
        (
            "slab region 0.06 → 0.2, skewed box",
            &skewed,
            SizeField::new(
                vec![SizeSource {
                    kind: SourceKind::Region,
                    ..datum(
                        "band",
                        Shape::Slab {
                            axis: 2,
                            min: 0.45,
                            max: 0.55,
                        },
                        0.06,
                        0.0,
                    )
                }],
                0.06,
                0.2,
                1.25,
            )
            .unwrap(),
        ),
        (
            "cylinder datum 0.05 → 0.2, unit box",
            &cube,
            SizeField::new(
                vec![datum(
                    "wake",
                    Shape::Frustum {
                        start: Vec3::new(0.0, 0.5, 0.5),
                        end: Vec3::new(1.0, 0.5, 0.5),
                        start_radius: 0.05,
                        end_radius: 0.05,
                    },
                    0.05,
                    0.0,
                )],
                0.05,
                0.2,
                1.3,
            )
            .unwrap(),
        ),
    ];
    if std::env::args().any(|a| a == "--gradation") {
        gradation_sweep(&cube);
        return;
    }
    println!(
        "| Case | Method | Seeds | Valid | Size ratio median | p90/p10 | Neighbour ratio p99 / max (max, interior cells) | Non-orth. mean / max (°) | Sphericity min / mean | Time (s) |"
    );
    println!("|---|---|---|---|---|---|---|---|---|---|");
    for (name, surface, field) in &cases {
        let n = target_count(surface, field);
        for method in ["sampled", "springs", "cvt"] {
            if method == "cvt" && std::env::args().any(|a| a == "--sweep") {
                for iterations in SWEEP {
                    let started = Instant::now();
                    let set = seed_sites(
                        surface,
                        field,
                        SeedOptions {
                            seed: 1,
                            relaxation_iterations: iterations,
                            max_seeds: usize::MAX,
                        },
                    )
                    .unwrap();
                    let row = measure(surface, field, &set.sites);
                    println!(
                        "| {name} | cvt {iterations} (ran {}, shift {:.4}) | {n} | {} | {:.3} | {:.2} | {:.2} / {:.2} ({:.2}) | {:.1} / {:.1} | {:.3} / {:.3} | {:.2} |",
                        set.iterations,
                        set.last_shift,
                        row.valid,
                        row.median,
                        row.spread,
                        row.neighbour_p99,
                        row.neighbour_max,
                        row.interior_max,
                        row.nonorth_mean,
                        row.nonorth_max,
                        row.sphericity_min,
                        row.sphericity_mean,
                        started.elapsed().as_secs_f64()
                    );
                }
                continue;
            }
            let started = Instant::now();
            let sites = match method {
                "sampled" => sample(surface, field, n, 1).unwrap(),
                "springs" => springs(surface, field, sample(surface, field, n, 1).unwrap()),
                _ => {
                    seed_sites(
                        surface,
                        field,
                        SeedOptions {
                            seed: 1,
                            relaxation_iterations: ITERATIONS,
                            max_seeds: usize::MAX,
                        },
                    )
                    .unwrap()
                    .sites
                }
            };
            let row = measure(surface, field, &sites);
            println!(
                "| {name} | {method} | {n} | {} | {:.3} | {:.2} | {:.2} / {:.2} ({:.2}) | {:.1} / {:.1} | {:.3} / {:.3} | {:.2} |",
                row.valid,
                row.median,
                row.spread,
                row.neighbour_p99,
                row.neighbour_max,
                row.interior_max,
                row.nonorth_mean,
                row.nonorth_max,
                row.sphericity_min,
                row.sphericity_mean,
                started.elapsed().as_secs_f64()
            );
        }
    }
}

struct Row {
    valid: &'static str,
    median: f64,
    spread: f64,
    neighbour_median: f64,
    neighbour_p99: f64,
    neighbour_max: f64,
    interior_max: f64,
    nonorth_mean: f64,
    nonorth_max: f64,
    sphericity_min: f64,
    sphericity_mean: f64,
}

fn measure(surface: &ConvexSurface, field: &SizeField, sites: &[Vec3]) -> Row {
    let Ok(generated) = generate_surface_mesh(surface, sites) else {
        return Row {
            valid: "no",
            median: f64::NAN,
            spread: f64::NAN,
            neighbour_median: f64::NAN,
            neighbour_p99: f64::NAN,
            neighbour_max: f64::NAN,
            interior_max: f64::NAN,
            nonorth_mean: f64::NAN,
            nonorth_max: f64::NAN,
            sphericity_min: f64::NAN,
            sphericity_mean: f64::NAN,
        };
    };
    let (mesh, metrics) = (generated.mesh(), generated.metrics());
    let size = |c: usize| metrics.cells[c].volume.cbrt();
    let mut ratios: Vec<f64> = (0..mesh.cell_count())
        .map(|c| size(c) / field.size(metrics.cells[c].centroid))
        .collect();
    ratios.sort_by(f64::total_cmp);
    let pick = |v: &[f64], q: f64| v[((v.len() - 1) as f64 * q).round() as usize];
    let mut neighbour = Vec::new();
    let mut nonorth = Vec::new();
    let mut area = vec![0.0; mesh.cell_count()];
    let mut on_boundary = vec![false; mesh.cell_count()];
    for (face, face_metrics) in metrics.faces.iter().enumerate() {
        let adjacency = mesh.face_adjacency(face).unwrap();
        area[adjacency.owner()] += face_metrics.area;
        if adjacency.neighbour().is_none() {
            on_boundary[adjacency.owner()] = true;
        }
    }
    let mut interior_max: f64 = 1.0;
    for (face, face_metrics) in metrics.faces.iter().enumerate() {
        if let FaceAdjacency::Interior {
            owner,
            neighbour: other,
        } = mesh.face_adjacency(face).unwrap()
        {
            area[other] += face_metrics.area;
            let (a, b) = (size(owner), size(other));
            neighbour.push(a.max(b) / a.min(b));
            if !on_boundary[owner] && !on_boundary[other] {
                interior_max = interior_max.max(a.max(b) / a.min(b));
            }
            nonorth.push(face_metrics.non_orthogonality_degrees);
        }
    }
    neighbour.sort_by(f64::total_cmp);
    let sphericity: Vec<f64> = (0..mesh.cell_count())
        .map(|c| 36.0 * std::f64::consts::PI * metrics.cells[c].volume.powi(2) / area[c].powi(3))
        .collect();
    Row {
        valid: "yes",
        median: pick(&ratios, 0.5),
        spread: pick(&ratios, 0.9) / pick(&ratios, 0.1),
        neighbour_median: pick(&neighbour, 0.5),
        neighbour_p99: pick(&neighbour, 0.99),
        neighbour_max: *neighbour.last().unwrap_or(&1.0),
        interior_max,
        nonorth_mean: nonorth.iter().sum::<f64>() / nonorth.len().max(1) as f64,
        nonorth_max: nonorth.iter().copied().fold(0.0, f64::max),
        sphericity_min: sphericity.iter().copied().fold(f64::INFINITY, f64::min),
        sphericity_mean: sphericity.iter().sum::<f64>() / sphericity.len() as f64,
    }
}

/// DistMesh-style relaxation (Persson & Strang 2004, section 2): bars are the
/// Delaunay edges (Voronoi neighbours), the target length is the size field
/// at the bar midpoint scaled so that the bars are slightly compressed
/// (F_scale = 1.2), forces are repulsive only, and a forward Euler step
/// (Δt = 0.2) moves the points. Points that would leave the convex domain
/// take a halved step instead (the paper projects them onto the boundary,
/// which a Voronoi seed must not reach).
fn springs(surface: &ConvexSurface, field: &SizeField, mut points: Vec<Vec3>) -> Vec<Vec3> {
    for _ in 0..ITERATIONS {
        let cells = clipped_surface_voronoi(surface, &points).unwrap();
        let mut bars = Vec::new();
        for (i, cell) in cells.iter().enumerate() {
            for face in cell.faces() {
                if let FaceSource::Bisector { other_site } = face.source()
                    && i < other_site
                {
                    bars.push((i, other_site));
                }
            }
        }
        let lengths: Vec<f64> = bars
            .iter()
            .map(|&(i, j)| (points[i] - points[j]).length())
            .collect();
        let targets: Vec<f64> = bars
            .iter()
            .map(|&(i, j)| field.size((points[i] + points[j]) * 0.5))
            .collect();
        let scale = 1.2
            * (lengths.iter().map(|l| l.powi(3)).sum::<f64>()
                / targets.iter().map(|h| h.powi(3)).sum::<f64>())
            .cbrt();
        let mut moves = vec![Vec3::ZERO; points.len()];
        for (k, &(i, j)) in bars.iter().enumerate() {
            let force = (targets[k] * scale - lengths[k]).max(0.0);
            let direction = (points[i] - points[j]) / lengths[k];
            moves[i] = moves[i] + direction * (0.2 * force);
            moves[j] = moves[j] - direction * (0.2 * force);
        }
        for (p, m) in points.iter_mut().zip(moves) {
            let mut step = m;
            for _ in 0..20 {
                if surface.contains_interior(*p + step) {
                    *p = *p + step;
                    break;
                }
                step = step * 0.5;
            }
        }
    }
    points
}

/// Achieved neighbour size ratios against the field's slope bound, for the
/// point datum case with CVT (ITERATIONS iterations).
fn gradation_sweep(cube: &ConvexSurface) {
    println!("| Field gradation | Seeds | Neighbour ratio median / p99 / max | Interior max |");
    println!("|---|---|---|---|");
    for gradation in [1.1, 1.2, 1.3, 1.5, 2.0] {
        let field = SizeField::new(
            vec![SizeSource {
                name: "tip".into(),
                kind: SourceKind::Datum,
                shape: Shape::Point(Vec3::new(0.25, 0.25, 0.25)),
                frame: Frame::GLOBAL,
                size: 0.04,
                distance: 0.05,
                growth: None,
            }],
            0.04,
            0.25,
            gradation,
        )
        .unwrap();
        let set = seed_sites(
            cube,
            &field,
            SeedOptions {
                seed: 1,
                relaxation_iterations: ITERATIONS,
                max_seeds: usize::MAX,
            },
        )
        .unwrap();
        let row = measure(cube, &field, &set.sites);
        println!(
            "| {gradation} | {} | {:.2} / {:.2} / {:.2} | {:.2} |",
            set.sites.len(),
            row.neighbour_median,
            row.neighbour_p99,
            row.neighbour_max,
            row.interior_max
        );
    }
}
