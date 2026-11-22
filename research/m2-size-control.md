# M2 size control: size field and seed placement (ALG-05, GRD-01 to GRD-05)

Design note of ALG-16 for the first grading of the convex STL mesher. Question: how does a size field h(x) from datum objects (GRD-02) and coordinate regions (GRD-03) become Voronoi seeds whose cells follow it, with a gradation limit (GRD-05)? The outcome is built in [`mesh_algorithms::sizing`](../crates/mesh-algorithms/src/sizing/mod.rs) and [`mesh_engine::sizing`](../crates/mesh-engine/src/sizing.rs); usage is in [m2-stl.md](../docs/m2-stl.md).

## The size field

Each source is a shape with a target size `h_i`, an influence distance `d_i` and a growth ratio `g_i`. With `β = gradation − 1` and the rate `r_i = min(g_i − 1, β)`:

h(x) = clamp( min( h_max, min_i [ h_i + r_i · max(0, dist_i(x) − d_i) ] ), h_min, h_max )

`dist_i` is the exact Euclidean distance to the shape in the source's frame: point, segment, two-sided plane, box, slab, sphere, and a solid frustum that covers cylinders and cones (the distance in the meridian half-plane to a trapezoid). Bodies are solid, so the distance is zero inside them. Regions are datum sources with `d_i = 0` and the gradation as growth, limited to boxes, slabs and cylinders in a named frame (GRD-03).

**Gradation without smoothing.** The classical method limits the gradient of a sampled size function on a background grid by solving `|∇h| ≤ β` (Persson 2006). For sources given as distance fields the limited field has a closed form: the smallest β-Lipschitz function above `h_i` on a set `S_i` is `h_i + β·dist(x, S_i)`, every distance function is 1-Lipschitz, and the minimum and the clamps keep the bound. So the field above satisfies `|∇h| ≤ β` exactly, everywhere, with no grid, no iterations and no Monte Carlo noise. The Walk-on-Spheres diffusion of the screening note (Sawhney & Crane 2020) is not needed for these sources; it was kept as a candidate for the geometry-derived sources of GRD-04, which turned out not to need it either: written as weighted surface primitives they keep the closed form ([note](grd-04-geometry-sizing.md)). Distances are Euclidean, which equals the geodesic distance through the domain only for convex domains, as in M2; nonconvex domains (M4) need the geodesic distance.

## Seed placement: candidates

| | Method | Source |
|---|---|---|
| Graphics | **Density-weighted centroidal Voronoi tessellation** by Lloyd iterations on the clipped cells, ρ = h⁻⁵ | Du, Faber & Gunzburger (1999); Lévy's group: Liu, Wang, Lévy et al. (2009), Yan, Wang, Lévy & Liu (2010) |
| Graphics, deferred | Semi-discrete optimal transport: a power diagram with cell volumes h³ | Lévy (2015); Lévy (2022) |
| Classical | **DistMesh** bar forces toward target lengths h(midpoint) | Persson & Strang (2004) |
| Baseline | Rejection sampling with acceptance (h_lo/h)³, no relaxation | — |

The density exponent: in 3D the cell diameter of a CVT scales as ρ^(−1/(d+2)) = ρ^(−1/5) (Gersho's asymptotics, used by Du, Faber & Gunzburger), so ρ = h⁻⁵ gives cell sizes ∝ h. The seed count is `N = round(V · mean(h⁻³))` over a fixed Halton sequence of 16,384 points in the bounding box, so a cell of size h has volume h³ on average, the count is a pure function of the inputs, and `mesh.estimate` reports it exactly. All three placements start from the same N seeds drawn by rejection sampling from a SplitMix64 generator (Steele, Lea & Flood 2014) seeded by the config's `seed`. The CVT centroid integral uses a degree-2 four-point rule on tetrahedra fanned from each cell's vertex mean. A CVT seed is the weighted centroid of its convex cell, so it stays strictly inside the convex domain; DistMesh steps that would leave the domain are halved (the paper projects onto the boundary, which a Voronoi seed must not reach). DistMesh is written in the [benchmark](../crates/mesh-engine/examples/sizing_benchmark.rs) only, from the paper, with its constants F_scale = 1.2 and Δt = 0.2; neighbour pairs are the Voronoi neighbours.

Optimal transport is deferred, not rejected: it needs a weighted (power) clipper and a Newton solve over all seeds, and the CVT below already meets the size targets to a few per cent; it is the candidate when cell volumes must match h³ exactly (GRD-06, MD-03).

## Benchmark

`cargo run --release -p mesh-engine --example sizing_benchmark [-- --sweep | --gradation]`, Linux x86_64, release build, 2026-11-22, 20 iterations for springs and CVT. Every mesh was built and validated by the production path (`generate_surface_mesh`: closure, orientation, volume and patch-area conservation) before it was measured. Size ratio is cbrt(V)/h(centroid); neighbour ratio is the larger over the smaller cbrt(V) across an interior face ("interior cells" excludes cells with a boundary face); sphericity is 36πV²/A³ (1 for a ball).

| Case | Method | Seeds | Valid | Size ratio median | p90/p10 | Neighbour ratio p99 / max (interior cells) | Non-orth. mean / max (°) | Sphericity min / mean | Time (s) |
|---|---|---|---|---|---|---|---|---|---|
| uniform h=0.125, unit box | sampled | 512 | yes | 0.975 | 1.55 | 1.69 / 2.23 (2.09) | 10.3 / 37.7 | 0.199 / 0.481 | 0.14 |
| | springs | 512 | yes | 0.948 | 1.40 | 1.49 / 1.86 (1.12) | 5.1 / 27.6 | 0.254 / 0.530 | 2.46 |
| | **cvt** | 512 | yes | **0.998** | **1.09** | **1.10 / 1.16 (1.10)** | **0.5 / 2.3** | **0.517 / 0.651** | 2.37 |
| point datum 0.04 → 0.25, unit box | sampled | 323 | yes | 0.971 | 1.51 | 2.10 / 3.21 (2.49) | 11.6 / 44.2 | 0.186 / 0.478 | 0.08 |
| | springs | 323 | yes | 1.045 | 1.49 | 1.55 / 1.79 (1.65) | 5.6 / 28.2 | 0.280 / 0.539 | 1.23 |
| | **cvt** | 323 | yes | **0.994** | **1.08** | **1.54 / 1.67 (1.64)** | **1.3 / 7.4** | **0.495 / 0.656** | 1.31 |
| slab region 0.06 → 0.2, skewed box | sampled | 1440 | yes | 0.975 | 1.51 | 1.83 / 3.17 (3.17) | 10.4 / 45.3 | 0.124 / 0.503 | 0.64 |
| | springs | 1440 | yes | 1.023 | 1.44 | 1.54 / 1.92 (1.49) | 4.2 / 23.6 | 0.275 / 0.576 | 12.50 |
| | **cvt** | 1440 | yes | **0.997** | **1.10** | **1.40 / 1.59 (1.54)** | **1.4 / 6.2** | **0.474 / 0.670** | 12.42 |
| cylinder datum 0.05 → 0.2, unit box | sampled | 658 | yes | 0.967 | 1.47 | 2.08 / 3.26 (3.15) | 11.3 / 47.8 | 0.138 / 0.495 | 0.20 |
| | springs | 658 | yes | 1.022 | 1.50 | 1.58 / 1.87 (1.67) | 4.4 / 25.5 | 0.288 / 0.576 | 3.61 |
| | **cvt** | 658 | yes | **0.993** | **1.08** | **1.52 / 1.71 (1.68)** | **1.7 / 6.4** | **0.520 / 0.666** | 3.68 |

Gradations: 1.2 (uniform), 1.3 (point, cylinder), 1.25 (slab).

**Validity** is equal: all twelve meshes pass. **Size control:** CVT follows the field to a 1.08 to 1.10 p90/p10 spread with the median within 1 %; springs and plain sampling spread 1.40 to 1.55. **Quality:** CVT has 4 to 10 times lower mean non-orthogonality (the CFD-relevant measure) and about twice the minimum sphericity. **Speed:** springs and CVT cost the same per iteration (both rebuild the all-pairs Voronoi diagram, ALG-04), so time is not a discriminator. **Decision: density-weighted CVT**, the graphics candidate, wins on every measure after validity.

### Iterations

| Case | 5 | 10 | 20 | 40 |
|---|---|---|---|---|
| uniform: p90/p10, max non-orth. (°) | 1.17, 7.8 | 1.12, 4.3 | 1.09, 2.3 | 1.06, 1.8 |
| point datum | 1.14, 16.2 | 1.10, 7.4 | 1.08, 7.4 | 1.06, 5.5 |
| slab region | 1.18, 11.6 | 1.14, 8.7 | 1.10, 6.2 | 1.08, 5.4 |
| cylinder datum | 1.16, 15.7 | 1.11, 10.1 | 1.08, 6.4 | 1.07, 5.8 |

Lloyd converges slowly (the last relative move after 20 iterations is about 0.03 h), as is known; 20 iterations is the default (`relaxation_iterations` changes it) because the gain from 20 to 40 is small and the cost doubles. Faster CVT solvers (L-BFGS, Liu et al. 2009) become worth it with ALG-04.

### Gradation limit on cells

The field obeys `|∇h| ≤ gradation − 1` exactly, but cells are not the field: even in a uniform field, CVT cell sizes vary by up to 16 % between neighbours after 20 iterations. Point datum, CVT 20 iterations:

| Field gradation | Seeds | Neighbour ratio median / p99 / max | Interior max |
|---|---|---|---|
| 1.1 | 1676 | 1.06 / 1.18 / 1.27 | 1.25 |
| 1.2 | 612 | 1.12 / 1.34 / 1.46 | 1.43 |
| 1.3 | 323 | 1.18 / 1.54 / 1.67 | 1.64 |
| 1.5 | 164 | 1.20 / 1.86 / 2.01 | 2.01 |
| 2.0 | 96 | 1.07 / 3.19 / 3.51 | 3.30 |

The median neighbour ratio stays below the gradation, the tail does not: the maximum is about 1 + 2.2·(gradation − 1). So the gradation bounds the field, and the build reports the achieved neighbour ratios and the number of faces above the limit in its evidence and limitations (AG-17) instead of claiming the limit holds per cell. The maintainer decided on 2026-11-22 to keep this: the gradation bounds the field, and the faces over it are reported, not enforced.

## What is reported (GRD-08)

`mesh.generate_stl` with a `size_field` stores, in the artifact and in its config digest, the field as given, the seed, the seed count, the iterations and the last move, the quantiles of the size ratio and of the neighbour ratio, the faces above the gradation, and per source the target size, the cells near it where it sets h, their achieved sizes, and the cells it governs.

## Limits

- Convex STL domains only, Euclidean distances (exact there).
- STL surface and closed-body datum objects are open; curvature, gap and feature-edge sizes (GRD-04) are in [grd-04-geometry-sizing.md](grd-04-geometry-sizing.md).
- All-pairs clipping (ALG-04) made each Lloyd iteration O(N²). Since the nearest-site clipper ([note](alg-04-neighbour-clipping.md)) an iteration is about O(N log N); the engine refuses fields that ask for more than 200,000 seeds (decided 2026-11-22), and the MCP limits apply to the seed count.
- Anisotropy (GRD-06) is open; Lp-CVT (Lévy & Liu 2010) is the lead.

## Clean room

Built from the papers below. No source code of any meshing or CVT library was opened (LIC-13); nothing was added to the licence register.

## Sources

- Q. Du, V. Faber, M. Gunzburger (1999), *Centroidal Voronoi tessellations: applications and algorithms*, SIAM Review 41(4).
- Y. Liu, W. Wang, B. Lévy, F. Sun, D.-M. Yan, L. Lu, C. Yang (2009), *On centroidal Voronoi tessellation: energy smoothness and fast computation*, ACM Transactions on Graphics 28(4).
- D.-M. Yan, W. Wang, B. Lévy, Y. Liu (2010), *Efficient computation of 3D clipped Voronoi diagram*, Geometric Modeling and Processing, LNCS 6130; extended as *Efficient computation of clipped Voronoi diagram for mesh generation*, Computer-Aided Design 45(4), 2013.
- B. Lévy (2015), *A numerical algorithm for L2 semi-discrete optimal transport in 3D*, ESAIM: M2AN 49(6).
- B. Lévy, Y. Liu (2010), *Lp centroidal Voronoi tessellation and its applications*, ACM Transactions on Graphics 29(4).
- A. Gersho (1979), *Asymptotically optimal block quantization*, IEEE Transactions on Information Theory 25(4).
- P.-O. Persson, G. Strang (2004), *A simple mesh generator in MATLAB*, SIAM Review 46(2).
- P.-O. Persson (2006), *Mesh size functions for implicit geometries and PDE-based gradient limiting*, Engineering with Computers 22.
- R. Sawhney, K. Crane (2020), *Monte Carlo geometry processing*, ACM Transactions on Graphics 39(4).
- J. H. Halton (1960), *On the efficiency of certain quasi-random sequences of points in evaluating multi-dimensional integrals*, Numerische Mathematik 2.
- G. L. Steele, D. Lea, C. H. Flood (2014), *Fast splittable pseudorandom number generators*, OOPSLA.
