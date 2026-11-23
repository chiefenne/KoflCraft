# Restricted Voronoi core for polyhedral ducts (M4)

Design note of ALG-16 for the polyhedral core of nonconvex ducts. Question: can a clipped Voronoi diagram, relaxed towards a centroidal Voronoi tessellation (CVT), replace the median dual of the Gmsh scaffold on the nonconvex elbow and Manifold fixtures, and does it remove the near-parallel faces that the median dual left (87° and 89.6° non-orthogonality after the fix of [m4-reentrant-dual-split](m4-reentrant-dual-split.md))? The code is [`duct/voronoi_core.rs`](../crates/mesh-algorithms/src/duct/voronoi_core.rs) and [`duct/convex.rs`](../crates/mesh-algorithms/src/duct/convex.rs), algorithm id `volume_clipped_voronoi_cvt_v1`. It is the default core; `"core": "median_dual"` in `volume_scaffold` keeps `volume_median_dual_v2`.

## Candidates

| | Method | Source |
|---|---|---|
| Classical | Vertex-centred median dual of the scaffold tetrahedra, cells split at reentrant edges (v2) | Standard finite-volume dual; the repair of the earlier note |
| Graphics (OBJ-19) | Voronoi diagram of sites restricted to the scaffold, computed per convex tetrahedron, with density-weighted Lloyd relaxation | Yan, Wang, Lévy & Liu (2013), *Efficient computation of clipped Voronoi diagram for mesh generation*; Lévy & Liu (2010) for the per-element restriction; Ray, Sokolov, Lefebvre & Lévy (2018), *Meshless Voronoi on the GPU*, for the security-radius stop; Du, Faber & Gunzburger (1999) for CVT |

No source code was read. Everything below is built from the papers and from first principles (LIC-12 to LIC-15).

## Construction

1. **Domain.** The scaffold tetrahedra after the wall offset are the convex pieces of the nonconvex core. The density is ρ = (h_min/h)^5 from the mean scaffold edge length h at each node, linear over each tetrahedron, so the CVT cell size follows the scaffold spacing.
2. **Restricted cells.** Each tetrahedron is clipped by the bisector planes of its nearest sites, in order of distance from the site nearest its centroid. The walk stops when the next site is farther than twice the largest vertex distance of the piece (the security radius). Each face of a piece remembers its plane: a tetrahedron face or the bisector towards one site.
3. **Lloyd.** Each site moves to the density-weighted centroid of its restricted cell. Ten iterations are the default (`lloyd_iterations`, at most 200). The last relative shift is reported.
4. **Final diagram.** Vertices are welded by position on a grid hash. A site's pieces are joined through tetrahedron faces into connected parts. Each part must have a kernel point, found by clipping a box with every face plane. A part without one is split greedily. Coplanar bisector polygons between two parts are merged when the result is planar, a single loop and star-shaped from its vertex mean.
5. **Reentrant edges.** A cell can wrap around a reentrant wall edge. Its centroid then lies behind one of the two walls, so the face to the wall layer there would have S·d ≤ 0. The site of such a cell is replaced by two sites mirrored across the plane that halves the angle between the walls, at 0.3 of the cell size, and the diagram is rebuilt. This repeats for up to three rounds. Then every face is still on a bisector. Cells that still wrap after that are repaired locally. A wall piece goes to a neighbouring cell if that cell stays star-shaped. Otherwise it is detached as a cell of its own, grown over neighbouring pieces up to 30% of the cell. As a last resort the cell is split per wall plane, as in v2.
6. **Robustness.** A repaired cell must keep a single polygon of contact with every neighbouring cell and with the wall cell over each wall triangle, because FoehnCraft allows one face per cell pair. The local repairs run for up to three passes. A part smaller than 1% of its neighbours' mean volume is merged into the neighbour with the largest contact. If the hand-off would still be rejected (an inverted face or a split contact), the engine rebuilds with one more Lloyd iteration, at most twice, and reports `retries`. The repairs depend on the exact site layout, which rounding can shift between platforms. The first CI run built a different Manifold mesh than the local build did, and its export failed. Lloyd counts 7 to 13 on Manifold were used as perturbed layouts. All seven are valid after these measures: maximum non-orthogonality 86.1° to 88.8°, maximum volume ratio 340 to 427. One of them (11 iterations) needed a retry.
7. **Wall layers ("split cap", now `"wall_join": "split_cap"`; the default is the wall columns of [m4-wall-columns](m4-wall-columns.md)).** The coarse wall wedges of the offset are all kept. An innermost wedge whose wall triangle touches one core cell stays a native prism. One whose triangle is shared by several core cells becomes a polyhedron. Its faces are the top triangle, the three side quads with the core vertices on their inner edges inserted, and the core wall polygons, reversed. Such a cell is validated through its wedge: the wedge Jacobian, a closed shell and the wedge volume.

## Benchmark

Release build, Linux x86_64. Meshes were read back with FoehnCraft `foehn mesh` 0.7.0 (main of 2026-11-22). Both columns use the hand-off of this PR, which merges the boundary faces one cell has on one patch into one polygon. This is why the median-dual Manifold count of faces above 70° is lower than in the earlier note (7,012 there).

| Case | Core | Cells | Inverted | Max / mean non-orth. | Faces > 70° | Max skewness | Max volume ratio | Core sphericity min / p05 / median | Build |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|
| Elbow | median dual v2 | 7,005 | 0 | 87.0° / 20.3° | 37 | 2.72 | 80.9 | 0.476 / 0.566 / 0.762 | 0.04 s |
| Elbow | Voronoi CVT | 1,355 | 0 | 85.9° / 21.1° | 44 | 0.76 | 11.4 | 0.731 / 0.783 / 0.824 | 0.2 s |
| Manifold | median dual v2 | 417,457 | 0 | 89.6° / 11.3° | 3,038 | 5.12 | 1,341 | 0.344 / 0.695 / 0.833 | 7 s |
| Manifold | Voronoi CVT | 92,058 | 0 | 87.7° / 11.2° | 605 | 0.99 | 355 | 0.352 / 0.809 / 0.877 | 75 s |

Volumes and patch areas match the scaffold to 1e-10 in both cores. All MSH 4.1 subdivisions load valid as well. On Manifold 11 sites were split in pairs over three rounds (2, 3 and 6). After that, two wall pieces were reassigned, two cells detached (2.2% and 4.4% of their cells) and one cell split. The least round cell, 0.35, is one of the detached ones. Without it the minimum is 0.70.

Almost every face above 60° lies between a core cell and an innermost wall cell: all 350 on the elbow and 6,769 of 6,787 on Manifold. The rest are 17 faces between layer cells and a single core-to-core face. Core-to-core faces are bisectors. The bad faces come from the wall join: a split-cap cell has one centroid over its whole wall triangle, while a core cell that touches only a corner of that triangle sees it far to the side.

### Wall join variants tried

| Variant | Elbow max / > 70° | Manifold | Why not |
|---|---|---|---|
| Split cap (kept) | 85.9° / 44 | valid, 87.7° / 605 | |
| Poly-prism column per core cell (wall polygon extruded through all layers, merged up to 30° feature angle) | 80.0° / 2 | 2 inverted faces, 129°, FoehnCraft volume 0.14% low | Over curved walls the columns bend. Their merged side faces are not planar, so FoehnCraft's face centroids and volumes are off |
| Columns over planar charts only | | 11 inverted faces, 131 above 70° | Same cause at the chart borders |

Outcome: the Voronoi CVT core. It is valid on both fixtures, has 5 times (elbow) and 4.5 times (Manifold) fewer cells, and gives rounder cells with lower skewness and volume ratio. Its worst face is slightly better on both fixtures. On Manifold it has 5 times fewer faces above 70°. On the elbow it has 7 more, all at the wall join. The price is build time: 75 s against 7 s on Manifold. Ten Lloyd iterations take about 2 s each, and each rebuild for a site split about 12 s, with welding the largest part.

## Limits

- Columns per core cell are the better wall join on flat walls (80° on the elbow). They are valid on curved walls since [m4-wall-columns](m4-wall-columns.md), which made them the default join.
- The mirrored site pairs at reentrant edges are a seed-placement rule inside the domain, not mirroring across the boundary. They are still close to the mirrored-seed family, so the LIC-08 patent check should cover them.
- The fallback repairs can leave small cells (2 to 4% of a neighbour) at reentrant edges. The validity of the reentrant repair rests on local rules plus a rebuild, not on a construction that is valid by design. Manifold's counts can therefore vary slightly between platforms, and the regression checks them within a band. A retry costs a full rebuild.
- `mesh.estimate` counts one cell per scaffold node and one per wall triangle and layer. Cells added at reentrant edges come later (none on the elbow, 14 on Manifold), and `cells_basis` says so.
- The scaffold still comes from Gmsh and a Python script (OPEN-07). A Rust tetrahedralization is its own step.

## Regression

`nonconvex_elbow_and_manifold_preserve_volume_and_expose_each_port_layer` asserts the Voronoi counts, sphericity floors and `inverted_faces = 0` for all duct meshes. `median_dual_core_stays_available_and_valid_on_the_elbow` keeps v2 covered. `cargo m1 --require-foehn` reads every duct back with `foehn mesh`.
