use mesh_core::geometry::Vec3;

pub const SEED: u64 = 2026;
pub const LATTICES: [[u32; 3]; 3] = [[2, 2, 4], [4, 4, 4], [4, 8, 8]];

pub fn lattice([nx, ny, nz]: [u32; 3]) -> Vec<Vec3> {
    let coordinate = |i, n| (f64::from(i) + 0.5) / f64::from(n);
    let mut sites = Vec::new();
    for x in 0..nx {
        for y in 0..ny {
            for z in 0..nz {
                sites.push(Vec3::new(
                    coordinate(x, nx),
                    coordinate(y, ny),
                    coordinate(z, nz),
                ));
            }
        }
    }
    sites
}

pub fn clustered(count: usize) -> Vec<Vec3> {
    let mut seed = SEED;
    (0..count)
        .map(|_| {
            Vec3::new(
                sample(&mut seed).powi(3),
                sample(&mut seed),
                sample(&mut seed),
            )
        })
        .collect()
}

pub fn sample(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
    (f64::from((*seed >> 32) as u32) + 0.5) / 4_294_967_296.0
}
