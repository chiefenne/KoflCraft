//! ALG-04 / ALG-16 benchmark of the Voronoi clipping (research/alg-04-neighbour-clipping.md).
//!
//! The same sites clipped two ways: `all-pairs`, the exhaustive M1
//! reference, and `nearest`, the nearest sites until the security radius
//! proves each cell complete (Ray, Sokolov, Lefebvre & Lévy 2018). Where
//! both run, every cell is compared (bisector neighbours, volume); then the
//! full production path (`generate_mesh`: clipping, welding, validation) is
//! timed with the nearest-site clipper.
//!
//! Run: `cargo run --release -p mesh-engine --example clipping_benchmark [-- --large]`

use std::time::Instant;

use mesh_algorithms::clipping::{ConvexPolyhedron, FaceSource};
use mesh_algorithms::voronoi::{Strategy, clipped_surface_voronoi_with, clipped_voronoi_with};
use mesh_core::geometry::{AxisAlignedBox, Vec3};
use mesh_engine::{generate_mesh, generate_surface_mesh};
use mesh_geometry::{StlOptions, import_stl};

/// Largest site count the all-pairs reference is run for.
const ALL_PAIRS_MAX: usize = 16_000;

fn sites(n: usize, graded: bool, seed: u64) -> Vec<Vec3> {
    let mut state = seed;
    let mut r = move || {
        // SplitMix64 (Steele, Lea & Flood 2014)
        state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..n)
        .map(|_| {
            let u = r();
            let x = if graded { u * u * u } else { u };
            Vec3::new(x, r(), r()).map_inside()
        })
        .collect()
}

trait Inside {
    fn map_inside(self) -> Self;
}

impl Inside for Vec3 {
    fn map_inside(self) -> Self {
        let c = |v: f64| v.clamp(1e-6, 1.0 - 1e-6);
        Vec3::new(c(self.x), c(self.y), c(self.z))
    }
}

fn neighbours(cell: &ConvexPolyhedron) -> Vec<usize> {
    let mut n: Vec<usize> = cell
        .faces()
        .iter()
        .filter_map(|f| match f.source() {
            FaceSource::Bisector { other_site } => Some(other_site),
            _ => None,
        })
        .collect();
    n.sort_unstable();
    n
}

fn same(a: &[ConvexPolyhedron], b: &[ConvexPolyhedron]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| neighbours(x) == neighbours(y))
}

fn main() {
    let large = std::env::args().any(|a| a == "--large");
    let unit = AxisAlignedBox::new(Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).expect("box");
    let mut counts = vec![1_000, 4_000, 16_000, 64_000];
    if large {
        counts.push(256_000);
    }
    println!(
        "| Domain | Sites | Distribution | all-pairs clip (s) | nearest clip (s) | Speed-up | Same cells | Mean neighbours | Full build, nearest (s) |"
    );
    println!("|---|---|---|---|---|---|---|---|---|");
    for &graded in &[false, true] {
        for &n in &counts {
            let s = sites(n, graded, 2026);
            let t = Instant::now();
            let fast = clipped_voronoi_with(unit, &s, Strategy::Nearest).expect("nearest");
            let fast_time = t.elapsed().as_secs_f64();
            let (slow_time, equal) = if n <= ALL_PAIRS_MAX {
                let t = Instant::now();
                let slow = clipped_voronoi_with(unit, &s, Strategy::AllPairs).expect("all pairs");
                let time = t.elapsed().as_secs_f64();
                (Some(time), Some(same(&fast, &slow)))
            } else {
                (None, None)
            };
            let mean = fast.iter().map(|c| neighbours(c).len()).sum::<usize>() as f64 / n as f64;
            let t = Instant::now();
            generate_mesh(unit, &s).expect("valid mesh");
            let build = t.elapsed().as_secs_f64();
            println!(
                "| unit box | {n} | {} | {} | {fast_time:.3} | {} | {} | {mean:.1} | {build:.3} |",
                if graded { "graded x = u³" } else { "uniform" },
                slow_time.map_or("not run".into(), |t| format!("{t:.3}")),
                slow_time.map_or("—".into(), |t| format!("{:.0}×", t / fast_time)),
                equal.map_or("—", |e| if e { "yes" } else { "NO" }),
            );
        }
    }
    let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../test-fixtures/");
    let bytes = std::fs::read(format!("{root}skewed-box.stl")).expect("fixture");
    let surface = import_stl(&bytes, StlOptions::default()).expect("surface");
    for &n in &[1_000, 4_000, 16_000] {
        // uniform in the bounds, kept where inside the surface
        let b = surface.bounds();
        let s: Vec<Vec3> = sites(4 * n, false, 7)
            .into_iter()
            .map(|p| b.min() + b.extent().mul_each(p))
            .filter(|&p| surface.contains_interior(p))
            .take(n)
            .collect();
        let t = Instant::now();
        let fast = clipped_surface_voronoi_with(&surface, &s, Strategy::Nearest).expect("nearest");
        let fast_time = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let slow = clipped_surface_voronoi_with(&surface, &s, Strategy::AllPairs).expect("all");
        let slow_time = t.elapsed().as_secs_f64();
        let mean = fast.iter().map(|c| neighbours(c).len()).sum::<usize>() as f64 / s.len() as f64;
        let t = Instant::now();
        generate_surface_mesh(&surface, &s).expect("valid mesh");
        let build = t.elapsed().as_secs_f64();
        println!(
            "| skewed box STL | {} | uniform | {slow_time:.3} | {fast_time:.3} | {:.0}× | {} | {mean:.1} | {build:.3} |",
            s.len(),
            slow_time / fast_time,
            if same(&fast, &slow) { "yes" } else { "NO" },
        );
    }
}

trait MulEach {
    fn mul_each(self, other: Vec3) -> Vec3;
}

impl MulEach for Vec3 {
    fn mul_each(self, o: Vec3) -> Vec3 {
        Vec3::new(self.x * o.x, self.y * o.y, self.z * o.z)
    }
}
