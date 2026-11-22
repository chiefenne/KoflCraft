//! Size control for the convex STL mesher (ALG-05, GRD-01 to GRD-03, GRD-05).
//!
//! The size field becomes a seed set in three steps (design note
//! `research/m2-size-control.md`):
//!
//! 1. The seed count is `N = round(V · mean(h⁻³))`, so a cell of size h has
//!    volume h³ on average. The mean comes from a fixed Halton sequence in the
//!    bounding box, so the count is a pure function of the surface and the
//!    field, and `mesh.estimate` reports it exactly.
//! 2. N seeds are drawn by rejection sampling with acceptance `(h_lo/h)³`
//!    from a seeded generator (M1-11 determinism).
//! 3. Density-weighted Lloyd relaxation of the clipped Voronoi cells (a
//!    centroidal Voronoi tessellation, CVT) with density ρ = h⁻⁵: in 3D the
//!    CVT cell diameter scales as ρ^(-1/5) (Gersho; Du, Faber & Gunzburger
//!    1999), so this density makes the cell size follow h. Each new seed is
//!    the ρ-weighted centroid of its convex cell, so it stays strictly inside
//!    the convex domain.

mod field;

pub use field::{Frame, Shape, SizeField, SizeSource, SizingError, SourceKind};

use mesh_core::geometry::Vec3;
use mesh_geometry::ConvexSurface;

use crate::clipping::ConvexPolyhedron;
use crate::voronoi::clipped_surface_voronoi;

/// Halton points used to integrate h⁻³ for the seed count.
pub const COUNT_SAMPLES: usize = 16_384;
/// Rejection-sampling candidates allowed per seed before giving up.
const CANDIDATES_PER_SEED: usize = 20_000;
/// Lloyd stops early once no seed moves by more than this fraction of h.
pub const RELAXATION_TOLERANCE: f64 = 1e-3;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SeedOptions {
    /// Random seed of the sampler.
    pub seed: u64,
    /// Lloyd iterations at most (0 keeps the sampled seeds).
    pub relaxation_iterations: usize,
    /// The largest seed count accepted.
    pub max_seeds: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SeedSet {
    pub sites: Vec<Vec3>,
    /// The seed count step 1 asked for (equal to `sites.len()`).
    pub target_count: usize,
    pub iterations: usize,
    /// The largest seed move of the last iteration, relative to h there.
    pub last_shift: f64,
}

/// The seed count of step 1.
pub fn target_count(surface: &ConvexSurface, field: &SizeField) -> usize {
    let bounds = surface.bounds();
    let (min, extent) = (bounds.min(), bounds.extent());
    let mut sum = 0.0;
    let mut inside = 0usize;
    for i in 0..COUNT_SAMPLES {
        let index = i as u64 + 1;
        let p = min
            + Vec3::new(
                extent.x * halton(index, 2),
                extent.y * halton(index, 3),
                extent.z * halton(index, 5),
            );
        if surface.contains_interior(p) {
            inside += 1;
            sum += field.size(p).powi(-3);
        }
    }
    if inside == 0 {
        return 1;
    }
    let count = surface.volume() * sum / inside as f64;
    // Saturating float-to-int conversion; NaN cannot occur for a valid field.
    (count.round() as usize).max(1)
}

/// Steps 1 to 3: seeds whose Voronoi cells follow the size field.
pub fn seed_sites(
    surface: &ConvexSurface,
    field: &SizeField,
    options: SeedOptions,
) -> Result<SeedSet, SizingError> {
    let target = target_count(surface, field);
    if target > options.max_seeds {
        return Err(SizingError::TooManySeeds {
            target,
            limit: options.max_seeds,
        });
    }
    let mut sites = sample(surface, field, target, options.seed)?;
    let mut iterations = 0;
    let mut last_shift = f64::INFINITY;
    while iterations < options.relaxation_iterations {
        let cells = clipped_surface_voronoi(surface, &sites).map_err(SizingError::Voronoi)?;
        last_shift = 0.0;
        for (site, cell) in sites.iter_mut().zip(&cells) {
            let next = weighted_centroid(cell, field).unwrap_or(*site);
            last_shift = f64::max(last_shift, (next - *site).length() / field.size(*site));
            *site = next;
        }
        iterations += 1;
        if last_shift < RELAXATION_TOLERANCE {
            break;
        }
    }
    Ok(SeedSet {
        sites,
        target_count: target,
        iterations,
        last_shift,
    })
}

/// Step 2: rejection sampling, in a fixed order for a fixed seed.
pub fn sample(
    surface: &ConvexSurface,
    field: &SizeField,
    count: usize,
    seed: u64,
) -> Result<Vec<Vec3>, SizingError> {
    let bounds = surface.bounds();
    let (min, extent) = (bounds.min(), bounds.extent());
    let lowest = field.lower_bound();
    let mut random = SplitMix64(seed);
    let mut sites = Vec::with_capacity(count);
    let budget = count.saturating_mul(CANDIDATES_PER_SEED).max(1_000_000);
    for _ in 0..budget {
        if sites.len() == count {
            break;
        }
        let p = min
            + Vec3::new(
                extent.x * random.unit(),
                extent.y * random.unit(),
                extent.z * random.unit(),
            );
        let accept = (lowest / field.size(p)).powi(3);
        if random.unit() < accept && surface.contains_interior(p) {
            sites.push(p);
        }
    }
    if sites.len() < count {
        return Err(SizingError::SamplingBudget {
            placed: sites.len(),
            target: count,
        });
    }
    Ok(sites)
}

/// The ρ-weighted centroid of a convex cell, ρ = h⁻⁵, by a degree-2
/// four-point rule on the tetrahedra fanned from the vertex mean.
pub fn weighted_centroid(cell: &ConvexPolyhedron, field: &SizeField) -> Option<Vec3> {
    const A: f64 = 0.585_410_196_624_968_5;
    const B: f64 = 0.138_196_601_125_010_5;
    let vertices = cell.vertices();
    let centre = vertices.iter().fold(Vec3::ZERO, |sum, &v| sum + v) / vertices.len() as f64;
    let mut mass = 0.0;
    let mut moment = Vec3::ZERO;
    for face in cell.faces() {
        let ring = face.vertices();
        for k in 1..ring.len().saturating_sub(1) {
            let tet = [
                centre,
                vertices[ring[0]],
                vertices[ring[k]],
                vertices[ring[k + 1]],
            ];
            let volume = (tet[1] - tet[0])
                .dot((tet[2] - tet[0]).cross(tet[3] - tet[0]))
                .abs()
                / 6.0;
            let sum = tet[0] + tet[1] + tet[2] + tet[3];
            for corner in tet {
                let q = corner * (A - B) + sum * B;
                let weight = volume * 0.25 * field.size(q).powi(-5);
                mass += weight;
                moment = moment + q * weight;
            }
        }
    }
    (mass > 0.0 && mass.is_finite()).then(|| moment / mass)
}

/// The radical-inverse (Halton) coordinate of `index` in `base`.
fn halton(mut index: u64, base: u64) -> f64 {
    let mut result = 0.0;
    let mut fraction = 1.0 / base as f64;
    while index > 0 {
        result += fraction * (index % base) as f64;
        index /= base;
        fraction /= base as f64;
    }
    result
}

/// SplitMix64 (Steele, Lea & Flood 2014): a small, fixed, portable generator.
pub struct SplitMix64(pub u64);

impl SplitMix64 {
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in [0, 1) with 53 random bits.
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use mesh_core::geometry::PatchId;
    use mesh_geometry::SurfaceTriangle;

    pub(crate) fn unit_cube() -> ConvexSurface {
        let v = |x, y, z| Vec3::new(x, y, z);
        let vertices = vec![
            v(0.0, 0.0, 0.0),
            v(1.0, 0.0, 0.0),
            v(1.0, 1.0, 0.0),
            v(0.0, 1.0, 0.0),
            v(0.0, 0.0, 1.0),
            v(1.0, 0.0, 1.0),
            v(1.0, 1.0, 1.0),
            v(0.0, 1.0, 1.0),
        ];
        let quads = [
            [0, 3, 2, 1],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [2, 3, 7, 6],
            [0, 4, 7, 3],
            [1, 2, 6, 5],
        ];
        let triangles = quads
            .iter()
            .flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]])
            .map(|vertices| SurfaceTriangle {
                vertices,
                patch: PatchId(0),
            })
            .collect();
        ConvexSurface::new(vertices, triangles, vec!["surface".into()]).unwrap()
    }

    fn uniform(size: f64) -> SizeField {
        SizeField::new(vec![], size / 2.0, size, 1.2).unwrap()
    }

    // Keep rule (SW-18): a uniform field of size h in a unit cube asks for
    // round(1 / h³) seeds and the same inputs give the same seeds (M1-11).
    // A change to the count rule or an unseeded generator breaks this; the
    // estimate test only compares the count with the build.
    #[test]
    fn uniform_count_and_determinism() {
        let cube = unit_cube();
        assert_eq!(target_count(&cube, &uniform(0.25)), 64);
        let options = SeedOptions {
            seed: 7,
            relaxation_iterations: 3,
            max_seeds: 1_000,
        };
        let first = seed_sites(&cube, &uniform(0.25), options).unwrap();
        let second = seed_sites(&cube, &uniform(0.25), options).unwrap();
        assert_eq!(first, second);
        assert_eq!(first.sites.len(), 64);
        assert!(first.sites.iter().all(|&p| cube.contains_interior(p)));
        let other = seed_sites(&cube, &uniform(0.25), SeedOptions { seed: 8, ..options }).unwrap();
        assert_ne!(first.sites, other.sites);
        assert!(matches!(
            seed_sites(&cube, &uniform(0.05), options),
            Err(SizingError::TooManySeeds { target: 8000, .. })
        ));
    }

    #[test]
    fn weighted_centroid_of_a_uniform_box_is_its_centre() {
        let cell = ConvexPolyhedron::from_box(
            mesh_core::geometry::AxisAlignedBox::new(Vec3::ZERO, Vec3::new(2.0, 1.0, 1.0)).unwrap(),
        );
        let c = weighted_centroid(&cell, &uniform(0.5)).unwrap();
        assert!((c - Vec3::new(1.0, 0.5, 0.5)).length() < 1e-14);
    }

    #[test]
    fn halton_and_splitmix_are_fixed() {
        assert_eq!(halton(1, 2), 0.5);
        assert_eq!(halton(3, 2), 0.75);
        assert!((halton(5, 3) - 7.0 / 9.0).abs() < 1e-15);
        let mut random = SplitMix64(0);
        assert_eq!(random.next_u64(), 0xE220_A839_7B1D_CDAF);
    }
}
