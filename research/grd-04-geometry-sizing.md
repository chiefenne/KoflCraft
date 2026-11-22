# GRD-04 geometry-derived sizes: curvature, gaps and feature edges

Design note of ALG-16 for the last open piece of M2 size control. Question: how does the STL itself make the size field finer near curved walls, across narrow gaps between opposite walls and along sharp edges, on top of the datum objects and regions of [m2-size-control.md](m2-size-control.md)? The outcome is built in [`mesh_algorithms::sizing::surface`](../crates/mesh-algorithms/src/sizing/surface.rs) and [`mesh_engine::sizing`](../crates/mesh-engine/src/sizing.rs); usage is in [m2-stl.md](../docs/m2-stl.md#size-control).

## From the surface to the field

Each derived source turns the surface into **weighted primitives**: facets, facet pieces or edges, each with its own size `h_k` and a reach `ρ_k` over which that size holds. The source asks for

h(x) = min_k [ h_k + r · max(0, dist(x, P_k) − ρ_k − d) ]

with the exact Euclidean distance to the triangle (closest point by Voronoi regions, Ericson 2005) or segment, the source's optional `distance` d and the rate r = min(growth, gradation) − 1. This is the closed form of gradient limiting already used for datum objects, applied per primitive: each term is r-Lipschitz, so the field keeps `|∇h| ≤ gradation − 1` with no grid and no smoothing. The Walk-on-Spheres diffusion named as a candidate in the M2 note is therefore not needed for these sources either. A bounding-volume hierarchy stores the smallest size and the largest reach of each node; the bound `h_min(node) + r·max(0, dist(x, box) − ρ_max(node) − d)` prunes the minimum exactly (a unit test compares it with brute force on 300 mixed primitives), and the running minimum over the other sources bounds the search.

Derived sizes below `min_size` are held at it and counted; there the target is not met, and the result says so.

## Curvature: candidates

| | Method | Source |
|---|---|---|
| Graphics, **chosen** | **Normal-cycle curvature tensor**: T = (1/area) Σ_e β_e·|e| ê êᵀ over a region, with β_e the angle between the facet normals at edge e; each edge counts half for each of its facets, the region is the facet's vertex ring without facets across a feature edge, and κ_max is the largest eigenvalue of T (convex surface: β ≥ 0) | Cohen-Steiner & Morvan (2003); used for remeshing by Alliez, Cohen-Steiner, Devillers, Lévy & Desbrun (2003) |
| Graphics | **Cotangent mean curvature and angle-defect Gaussian curvature** per vertex over the mixed Voronoi area, κ_max = H + sqrt(H² − K), averaged over the facet's corners | Meyer, Desbrun, Schröder & Barr (2003); the operators of Crane's discrete differential geometry notes (2018) |
| Classical | **Normal deviation**: the turn of the facet normal over the centroid distance, largest over the edge neighbours | the common engineering estimate, written from its definition |

The size is h = angle / κ_max: one cell spans at most `angle` of the wall normal's turn. Edges whose normals turn more than `feature_angle` (default 30°) are features and are left to the feature-edge source; without that rule a cylinder's rim would count as extreme curvature.

## Curvature: benchmark

`cargo run --release -p mesh-engine --example geometry_sizing_benchmark`, Linux x86_64, release, 2026-11-22. Error is |κ_est/κ − 1| per facet against the analytic κ_max at the facet (ellipsoid: closed-form H and K at the centroid's radial projection). The strip cylinders are written as CAD exporters do: one long triangle pair per side, fan caps.

| Surface | Facets | Estimator | median | p90 | max | max κ·R on flat facets |
|---|---|---|---|---|---|---|
| icosphere R=0.5 | 1280 | **normal-cycle** | 0.008 | 0.019 | 0.023 | – |
| | | cotan | 0.000 | 0.000 | 0.000 | – |
| | | normal-deviation | 0.048 | 0.081 | 0.147 | – |
| UV sphere R=0.5, 24×48 | 2208 | **normal-cycle** | 0.004 | 0.008 | 0.053 | – |
| | | cotan | 0.001 | 0.002 | 0.002 | – |
| | | normal-deviation | 0.421 | 0.498 | 0.507 | – |
| ellipsoid 0.5/0.35/0.25 | 1280 | **normal-cycle** | 0.014 | 0.030 | 0.078 | – |
| | | cotan | 0.007 | 0.018 | 0.047 | – |
| | | normal-deviation | 0.073 | 0.159 | 0.273 | – |
| strip cylinder R=0.25, n=48 | 192 | **normal-cycle** | 0.001 | 0.001 | 0.001 | 0.00 |
| | | cotan | 0.188 | 0.188 | 0.188 | 0.54 |
| | | normal-deviation | 0.902 | 0.902 | 0.902 | 0.00 |
| strip cylinder R=0.25, n=16 | 64 | **normal-cycle** | 0.006 | 0.006 | 0.006 | 0.00 |
| | | cotan | 0.206 | 0.206 | 0.206 | 0.53 |
| | | normal-deviation | 0.711 | 0.711 | 0.711 | 0.00 |

**Decision: the normal-cycle tensor.** The cotangent operators are the most accurate on smooth, well-shaped spheres, where every vertex lies on the surface, but every vertex of a CAD cylinder sits on a rim, so their vertex areas and angle defects mix in the caps: 19 to 21 % error on the wall and a false curvature of about 0.5/R on the flat caps, which would refine the caps for nothing. The classical normal deviation fails on anisotropic triangles (the centroid distance follows the long side: 71 to 90 % too low on strips, 42 to 51 % on the UV sphere). The normal cycle is within 2.5 % on the spheres, 8 % at worst on the ellipsoid, 0.1 to 0.6 % on the strip cylinders (the residual is the chord geometry, (π/n)/sin(π/n)), and exactly zero on flat facets, because curvature lives on the edges and features are cut out of the sum. For sizing, an 8 % error in κ is an 8 % error in h, well inside the spread of CVT cells (about ±10 %, see the M2 note).

## Gaps: candidates

A gap source asks for h = t / `cells_across` on each wall piece, held over the reach t/2 so the size spans the gap instead of growing towards its middle (without the reach, a slab with 4 cells asked would get cells of 0.4·t, 2.5 across, in its middle at gradation 1.3). Walls count as opposite when their normals make an angle above 180° − `wedge_angle` (default 45°): walls that face each other or meet at a sharp interior angle form a gap, a box's 90° edge does not. Without this rule every convex corner would be a gap of thickness zero, the well-known collapse of the local feature size at sharp edges.

| | Thickness t at a wall point p with outward normal n | Source |
|---|---|---|
| Graphics, **chosen (pending the card)** | **Medial ball**: the diameter of the largest ball tangent to the wall at p that the opposite walls leave free. The shrinking-ball method finds it by nearest-neighbour queries on surface samples; for a convex domain it is exact in closed form, t = min_g 2·D_g / (1 − n·n_g) over the opposite planes g, with D_g the distance of p from plane g | Ma, Bae & Choi (2012); medial balls as in Amenta, Choi & Kolluri (2001) |
| Classical | **Normal ray**: the distance along −n to the opposite wall, t = min_g D_g / (−n·n_g) | ray casting along the wall normal |

Both thicknesses are minima of affine functions of p, hence concave along a facet, so their minimum over a triangle is at a corner. Facets are split in four until the thickness varies by at most 10 % over a piece (at most five levels), and each piece asks for its corner minimum: exact at the corners, conservative inside.

## Gaps: benchmark

Same run. Four cells across asked, min_size 0.008, max_size 0.15, gradation 1.3, 20 Lloyd iterations; every mesh built and validated by `generate_surface_mesh`. Cells across are counted on the built mesh along the inward normal from a fixed barycentric grid on every facet (21 points per facet) to where the ray leaves the domain, as the number of nearest-site changes along 512 points (the cell containing a point of a clipped Voronoi mesh is its nearest site's), wherever four cells of min_size fit. The median is high because samples on the large side walls count cells along the slab, not across it.

GAP_TABLE

**Validity** is equal: all meshes pass. **Cells across** are equal sample for sample on every shape: no measure misses the target where the other meets it. The 5 % below target on the 30° wedge are samples whose path crosses one cell fewer than asked, a discretisation effect of CVT cells on a tapering gap, not of the thickness (both measures miss the same samples). **Cost**: the medial ball is never wider than the ray (for a wedge of opening α it is 2·tan(α/2) against tan α per unit of distance from the apex), so it asks for up to 11 % more seeds on the wedges and the sphere and the same on parallel walls. **Recommendation: the medial ball**, for two reasons beyond the tie. It errs on the safe side, and it is the measure that carries to the nonconvex domains of M4, where the shrinking-ball search replaces the closed form, while a ray from a wall can leave through an opening or graze a corner and report a gap that is not there or miss one that is. The choice is on a decision card to the maintainer; the build uses the medial ball until it is answered.

## Feature edges

An edge is a feature when the normals of its two facets turn more than `feature_angle` (default 30°), the dihedral-angle rule of the standard texts (Botsch, Kobbelt, Pauly, Alliez & Lévy 2010, ch. 2) and of classical meshers alike; there is no competing graphics candidate to benchmark for detection. The twelve edges of a box are found and its six face diagonals are not; a 48-sided strip cylinder has 96 rim edges and no other. Each feature edge asks for the source's `size` along the segment.

## Results on built meshes

The integration tests ([geometry_sizing.rs](../crates/mesh-engine/tests/geometry_sizing.rs)) build and validate: a 1 × 1 × 0.1 slab with 3 cells asked gets at least 3 on every sample, where the background alone gives one; a strip cylinder with `angle` 15° gets the wall cells at angle·R within 10 % (median) and nothing on the caps; a box with a feature size of 0.05 gets it along the edges within 10 %. A source's cells in the evidence are those whose centroid lies within half their own size of the source's reach (plus `distance`), because a derived source sits on the wall where no centroid lies.

The benchmark's production builds:

MESH_TABLE

## Limits

- Convex domains only: Euclidean distances equal distances through the domain and the thickness has a closed form. Nonconvex domains (M4) need the shrinking-ball search on surface samples and geodesic or visibility-aware distances.
- Curvature is per facet and isotropic: the cell size follows κ_max in every direction. Anisotropic cells along a cylinder axis wait for GRD-06 (Lp-CVT).
- Cells across a tapering gap are discrete; samples below the target are reported, not enforced.
- Cost: the field evaluation is a pruned search; builds near 8,000 seeds on a 192-facet disc take about 30 s, most of it in clipping each cell against all facets (an open ALG-04 limit).

## Clean room

Built from the papers and book sections below. No source code of any meshing, geometry-processing or medial-axis library was opened (LIC-13); nothing was added to the licence register.

## Sources

- D. Cohen-Steiner, J.-M. Morvan (2003), *Restricted Delaunay triangulations and normal cycle*, Symposium on Computational Geometry.
- P. Alliez, D. Cohen-Steiner, O. Devillers, B. Lévy, M. Desbrun (2003), *Anisotropic polygonal remeshing*, ACM Transactions on Graphics 22(3).
- M. Meyer, M. Desbrun, P. Schröder, A. H. Barr (2003), *Discrete differential-geometry operators for triangulated 2-manifolds*, Visualization and Mathematics III.
- K. Crane (2018), *Discrete differential geometry: an applied introduction*, course notes, Carnegie Mellon University.
- J. Ma, S. W. Bae, S. Choi (2012), *3D medial axis point approximation using nearest neighbors and the normal field*, The Visual Computer 28(1).
- N. Amenta, S. Choi, R. K. Kolluri (2001), *The power crust*, ACM Symposium on Solid Modeling and Applications.
- M. Botsch, L. Kobbelt, M. Pauly, P. Alliez, B. Lévy (2010), *Polygon mesh processing*, A K Peters.
- C. Ericson (2005), *Real-time collision detection*, Morgan Kaufmann, section 5.1.5.
- O. K. Smith (1961), *Eigenvalues of a symmetric 3 × 3 matrix*, Communications of the ACM 4(4).
- P.-O. Persson (2006), *Mesh size functions for implicit geometries and PDE-based gradient limiting*, Engineering with Computers 22.
