# ALG-04: nearest-site Voronoi clipping with a completeness proof

Design note of ALG-16 for replacing the all-pairs clipping of M1. Question: how can each clipped Voronoi cell be built from a few nearby sites instead of all N, without the heuristic cutoffs ALG-13 forbids? The outcome is built in [`mesh_algorithms::voronoi`](../crates/mesh-algorithms/src/voronoi.rs) and [`mesh_algorithms::neighbours`](../crates/mesh-algorithms/src/neighbours.rs).

## Method

The candidate comes from the screening of OBJ-19 (requirements section 9a): Ray, Sokolov, Lefebvre & Lévy (2018), *Meshless Voronoi on the GPU*, from Lévy's group. Built from the paper; no source code was read (LIC-13). The `meshless_voronoi` crate in the licence register (section 14a) was not opened.

1. **Nearest sites.** A static k-d tree over the sites (Bentley 1975; Friedman, Bentley & Finkel 1977), split at the median of the widest axis, leaves of 8. A query returns the k nearest sites of a site sorted by distance, ties by index, so the result is a pure function of the input. k starts at 32 and doubles when needed.
2. **Clipping in order of distance.** A cell starts as the domain (box or convex surface) and is clipped by the bisector of each neighbour in turn.
3. **Security radius.** Let R be the largest distance from the site to a vertex of the current cell. A site farther than 2R has its bisector beyond the ball of radius R around the site, so it cannot cut the cell. Clipping stops at the first neighbour farther than 2(R + τ)(1 + 10⁻⁹), with τ the clipping tolerance, so no farther bisector can even come within τ of a vertex. If the k neighbours run out first, k doubles; with all N − 1 sites the cell is exhaustive. The stop is a proof, not a cutoff (ALG-13).
4. **Shared vertices, bit for bit.** Cells sharing a vertex clip their planes in different orders, so their clipped vertices can differ by rounding. On 64,000 strongly graded sites two neighbours computed a shared vertex 5.8·10⁻¹³ apart, beyond the weld tolerance, and the mesh did not close (the all-pairs order happened to stay within it). Each simple vertex (exactly three incident faces) is therefore recomputed from its defining equations: the sites of its bisectors in index order and the boundary planes in a fixed order, solved by Gaussian elimination with partial pivoting. Every cell around the vertex solves the same system the same way, so the coordinates agree exactly. Degenerate vertices (four or more faces, as on lattices), singular systems and moves larger than 1024 τ keep the clipped position. This is the combinatorial vertex idea of the same paper line: a vertex is named by the planes that define it.

The exhaustive all-pairs clipper stays as `Strategy::AllPairs`, the reference that ALG-13 asks for.

## Benchmark

`cargo run --release -p mesh-engine --example clipping_benchmark -- --large`, Linux x86_64 (Xeon 2.1 GHz, one thread), release build, 2026-11-22. "Same cells": every cell has the same bisector neighbours in both strategies. "Full build" is the production path (`generate_mesh`: clipping, welding and validation) with the nearest-site clipper. Graded sites have x = u³ (density spans several decades near x = 0).

| Domain | Sites | Distribution | All-pairs clip (s) | Nearest clip (s) | Speed-up | Same cells | Mean neighbours | Full build (s) |
|---|---|---|---|---|---|---|---|---|
| unit box | 1,000 | uniform | 0.274 | 0.070 | 4× | yes | 13.4 | 0.104 |
| unit box | 4,000 | uniform | 2.682 | 0.330 | 8× | yes | 14.2 | 0.514 |
| unit box | 16,000 | uniform | 37.519 | 1.419 | 26× | yes | 14.6 | 2.184 |
| unit box | 64,000 | uniform | not run | 6.199 | — | — | 15.0 | 11.966 |
| unit box | 256,000 | uniform | not run | 29.387 | — | — | 15.2 | 58.316 |
| unit box | 1,000 | graded | 0.261 | 0.067 | 4× | yes | 12.4 | 0.103 |
| unit box | 4,000 | graded | 2.715 | 0.303 | 9× | yes | 13.2 | 0.466 |
| unit box | 16,000 | graded | 38.230 | 1.509 | 25× | yes | 13.9 | 2.513 |
| unit box | 64,000 | graded | not run | 6.752 | — | — | 14.3 | 12.588 |
| unit box | 256,000 | graded | not run | 28.968 | — | — | 14.6 | 55.133 |
| skewed box STL | 1,000 | uniform | 0.272 | 0.071 | 4× | yes | 13.4 | 0.122 |
| skewed box STL | 4,000 | uniform | 2.878 | 0.319 | 9× | yes | 14.1 | 0.627 |
| skewed box STL | 16,000 | uniform | 36.525 | 1.527 | 24× | yes | 14.6 | 2.805 |

All-pairs grows as N² (×14 from 4,000 to 16,000 sites); the nearest-site clipper grows as about N log N (×4.3 to ×4.7 per fourfold). All-pairs at 64,000 sites was not run; extrapolated it would take about 10 minutes for clipping alone.

`nearest_site_clipping_equals_the_all_pairs_reference` checks the same on the M1 lattices (2×2×4 to 8×8×8, full of cospherical sites), 256 clustered sites and 2,000 graded sites, also comparing face areas and volumes to 10⁻¹². Shortening the security radius to 0.6 of its value makes that test fail on the first clustered cell, so the test can see a missing neighbour.

## Decision and consequences

The nearest-site clipper is the default for box and STL meshes and for every CVT iteration of the size field. The build algorithm ids became `cpu_box_clipped_voronoi_v2` and `cpu_convex_surface_clipped_voronoi_v2`, because vertices can differ from v1 in the last bits. With the maintainer's card of 2026-11-22 the size-field cap went from 20,000 to 200,000 seeds and the MCP site limit from 256 to 20,000; builds over MCP above 500,000 estimated cells still ask for approval (AG-15).

## Limits

- Serial. Cells are independent, so threads (CORE-06) or the GPU form of the paper (GPU-02) are the next step.
- The cost of a cell grows with the sites within twice its radius: a cell next to a large empty region (sites only in part of the domain) needs many neighbours. Graded fields keep neighbouring sizes close, so this does not arise from a size field.
- A cell inside a convex STL starts from all its facets, so large STLs cost O(facets) per cell; clipping the surface by a bounding box of the cell first is open.
