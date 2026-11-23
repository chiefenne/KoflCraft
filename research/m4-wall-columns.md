# Wall columns under the Voronoi core (M4)

Design note of ALG-16 for the join between the prism wall layers and the Voronoi core of [m4-voronoi-core](m4-voronoi-core.md). Question: can the layer cells follow the core cells (one column per core cell) on curved walls too, with exact volumes and no inverted faces, and does that lower the worst non-orthogonality, which sat at the join (85.9° on the elbow, 87.7° on Manifold)? The code is [`duct/columns.rs`](../crates/mesh-algorithms/src/duct/columns.rs). It is the default join; `"wall_join": "split_cap"` in `volume_scaffold` keeps the split cap of PR #22.

## Why the first columns failed on curved walls

The first column variant of PR #22 extruded each core cell's wall polygon through all layers and merged it over the wall triangles up to a 30° feature angle. FoehnCraft integrates every face as a fan about its node average (`polygon_area_centroid` and `cell_volume_centroid_of` in `foehncraft-mesh/src/geometry.rs`, read as documentation of the consumer). A face merged over a curved wall is bent, and its fan surface lies on the chord, off the wall by about w²/8R. On Manifold, with core cells about 6 units wide on pipes of radius 25, that is about 0.18, against a first layer of 0.08. Two things follow. The merged wall face moves the wall, which costs 0.14% of the volume. The fans of the two caps of a thin layer cell cross, which inverts faces. Any cell whose cap spans non-coplanar triangles has this problem, so on curved walls no layer cell may span two wall triangles.

## Construction

1. **Units.** The core cells cut the innermost layer surface into wall polygons. A unit is the polygon (or polygons) of one core cell in one wall triangle. These are the facets of the restricted Voronoi diagram on the wall triangles (Yan, Lévy, Liu, Sun & Wang 2009, *Isotropic remeshing with fast and exact computation of restricted Voronoi diagram*), here taken from the clipped volume diagram.
2. **Extrusion.** A polygon vertex with barycentric coordinates λ on the innermost triangle is placed at λ on the same triangle of every layer surface. A layer prism is affine in λ at each height, so each extruded unit is a union of sub-wedges (a fan of its polygon), these sub-wedges partition the coarse prism exactly, and every cap is planar. A vertex on a triangle edge gets the same image from both triangles. Side faces are bilinear quads, as on a native prism.
3. **Columns.** Units of one core cell over adjacent triangles are merged when the triangles are coplanar on every layer surface: on the wall to 1e-10 of the model size, inside to 1e-6 of the thinnest layer. On a flat wall a core cell then sits on one layer cell per layer and wall plane, with its centroid under the core cell's. On a curved wall each unit stays its own column.
4. **One face per pair.** FoehnCraft allows one face between two cells. A merged column must be one disk, and its contact with each neighbouring column one open chain. A column that fails this is split back into its units. A unit with one polygon is convex, so two such units meet in one chain; if two units still meet in more, the build fails and the engine rebuilds with one more Lloyd iteration.
   A merged column must also pass FoehnCraft's orthogonality test against each neighbouring column on every layer: the merged side face's area vector S and the line d between the two centroids must satisfy S·d > cos 85° |S||d|. The centroids are the volume-weighted centroids of the sub-wedges. A core cell built from merged or reassigned pieces can have a nonconvex wall footprint, and the merged column's centroid can then lie behind a side face. Such a column is split back into its units, like a column that fails the contact test.
5. **Native prisms.** A unit that covers its whole triangle and merges with nothing keeps its native wedges.
6. **Gmsh export.** A layer cell whose faces are all wedge sides or tiled by wedge end triangles is written as its native sub-prisms; the core cells split their faces to match. FoehnCraft reads both MSH 4.1 files as valid.
7. **Hand-off.** The hand-off merges the boundary faces one cell has on one patch only when the merged polygon is flat to roundoff (relative 1e-10). A merged column's wall face on a flat wall becomes one face. Nothing bent is merged.

## Benchmark

Release build, Linux x86_64, ten Lloyd iterations. Meshes were read back with FoehnCraft `foehn mesh` 0.7.0 (main of 2026-11-23). Volumes are FoehnCraft's, compared with KoflCraft's exact sum.

| Case | Join | Cells | Inverted | Volume error | Max / mean non-orth. | Faces > 70° | Max skewness | Max volume ratio | Smallest cell |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Elbow | split cap | 1,355 | 0 | 3e-16 | 85.9° / 21.1° | 44 | 0.76 | 11 | 9.1e-4 |
| Elbow | columns | 1,082 | 0 | 2e-16 | 80.1° / 7.6° | 2 | 0.91 | 727 | 5.4e-6 |
| Elbow | units only (no merge) | 4,313 | 0 | 7e-16 | 86.9° / 15.9° | 9 | 1.96 | 2.7e7 | |
| Manifold | split cap | 92,058 | 0 | 7e-16 | 87.7° / 11.2° | 605 | 0.99 | 355 | 2.6e-3 |
| Manifold | columns | 162,435 | 0 | 0 | 82.8° / 7.2° | 67 | 2.34 | 1.5e10 | 4.7e-11 |
| Manifold | first columns (PR #22) | | 2 | −1.4e-3 | 129° | | | | |

The elbow's 285 columns cover 376 wall triangles. On Manifold 45,153 columns cover 21,694 wall triangles. Lloyd counts 7, 8, 9, 11, 12 and 13 on Manifold were used as perturbed layouts. All six are valid with 0 inverted faces: maximum non-orthogonality 83.3° to 85.9°, 60 to 73 faces above 70°. One of them (11) needed a retry, as with the split cap. The core cells and their faces are unchanged. Manifold builds in about 80 s, as before.

The worst faces on Manifold now lie between a curved-wall unit and its core cell when the unit sits at the edge of the core cell's footprint (82.8° to 81°), and between units of one core cell in two adjacent triangles (about 80°).

### Slivers

Where a core cell's boundary clips a wall triangle near a corner, its unit is a sliver. On Manifold 2,030 units cover less than 1% of their triangle, 206 less than 1e-4 and 3 less than 1e-8. Each gives one tiny cell per layer, hence the volume ratio and the skewness (the largest skewness, 2.3 to 8 on the perturbed layouts, is between two such slivers). Variants tried:

| Variant | Manifold | Why not |
|---|---|---|
| Sliver joins the unit of another core cell in its triangle (below 10%, 1%, 0.1% of the triangle) | 87.0° / 350, 85.9° / 232, 82.8° / 199 faces > 70°; up to 3 inverted lateral faces; rebuilds needed | The sliver's face to its own core cell is seen from the other cell's centroid. A joined column can be nonconvex |
| The same with a predicted-angle (70°) and convexity check | 82.8° / 253 | Still more faces above 70°, and the smallest slivers fail the check |
| Sliver joins its own core cell's unit across the crease | not built | The merged wall face is bent: kept as pieces, the sliver's wall face is seen from the column's centroid at about 88°; merged, the wall moves |

Removing slivers needs a change in the core: no core cell should touch a wall triangle in a sliver. That change is the conforming wall seeds below.

### Manifold on the Rust scaffold

On the finer Rust-built Manifold scaffold of PR #26 (56k wall triangles), the columns join first failed every layout. At Lloyd 10 and 11, nine lateral faces between a nonconvex merged column on a flat wall and a neighbouring unit were inverted, and two faces at sliver units were degenerate or seen from behind. The engine retried. At Lloyd 12 the core failed (MESH006, shell volume off by 1.3e-8). The core face merge (`merge_rings`) had checked only that the outer loop of the merged face was flat. Three tiny pieces formed a tent 8e-3 high whose rim lay in a flat face, and the merge flattened it. The merge now tests every vertex of every ring, to the shell check's tolerance (1e-10 of the face size; before 1e-9, which let a face through that the shell check then called nonplanar). The hand-off's boundary merge tests every vertex too.

With the S·d split of step 4, the layout at Lloyd 11 is valid (FoehnCraft 0.7.0: 505,937 cells, 0 inverted, volume error 2e-15, max non-orthogonality 87.9°, 412 faces above 70°, skewness 2.6). The layout at Lloyd 10 still has the two sliver faces. Slivers remained the limit: the volume ratio was about 1e13, and the Manifold fixture kept `"wall_join": "split_cap"` (89.0°, volume ratio 271) until the conforming wall seeds below. A full-run test builds Manifold with columns and requires 0 inverted faces.

A first try at the core-side sliver fix handed the pieces of a cell in the scaffold tetrahedron on a wall triangle to the cell covering most of the triangle when the cell covered below 1% of it, with the checks of the reentrant repairs. It moved 106 slivers on the elbow, but the smallest ones (below 1e-5 of the triangle) failed the kernel and in-front checks, and the elbow's smallest cell fell from 5.4e-6 to 1.5e-7. Not kept.

### Other candidates

A planar polygonal wall, built by intersecting tangent planes at the sites as for the planar hexagonal meshes of architectural geometry, would let a column have a flat wall face on curved walls. It changes the wall itself, so the scaffold's volume and areas would no longer be the reference. Not tried.

## Conforming wall seeds

The sliver fix. Code: [`duct/wall_seeds.rs`](../crates/mesh-algorithms/src/duct/wall_seeds.rs), used by the core in `voronoi_core.rs`. Built from the papers alone (LIC-13).

### Method

A wall unit is a sliver when the boundary of a core cell crosses a wall triangle near a corner or an edge. The fix is to let the core cells end on the wall edges: the Voronoi diagram is made to conform to the wall triangulation of the innermost layer surface.

1. **Balls.** Every wall vertex carries a ball. Its radius is 1.05 times the largest circumradius of the triangles at the vertex, so the balls of a triangle's corners cover it. The radius falls off along the wall by at most 0.3 per edge length. A deterministic spread of 1% on the squared radii keeps four corners of two triangles from lying on one sphere.
2. **Seeds.** The seed of a wall triangle is the point on the fluid side that lies on the balls of its three corners (the foot has equal power to the three balls; the height follows from one of them). Two seeds over triangles that share an edge are both at the ball radius from each end of the edge, so their bisector holds the edge. If no other site lies inside a ball, the Voronoi cell of a seed cuts the wall in exactly its triangle: the seed is the nearest site at the three corners, and a bisector half-space that holds no corner holds no point of the triangle. This is the seed placement of VoroCrust (Abdelkader, Bajaj, Ebeida, Mahmoud, Mitchell, Owens & Rushdi 2020, *VoroCrust: Voronoi meshing without clipping*, ACM Trans. Graph. 39(3); sampling conditions in the SoCG 2018 paper of the same group), used on one side of a given triangulation instead of on a fresh sampling.
3. **Regular triangulation.** A seed must lie outside every other ball. The triangles for which that holds are those of the weighted Delaunay (regular) triangulation of the balls. It is reached from the wall triangulation by edge flips with the empty-ball test in space (Edelsbrunner & Shah 1996, *Incremental topological flipping works for regular triangulations*, Algorithmica 15). The flips change the seeds only, not the wall mesh.
4. **Close seeds.** Two seeds on one edge closer than 0.2% of the edge leave the direction of their bisector to rounding. The ball opposite the edge shrinks by 3% and the flips are repeated, up to 12 rounds. A seed that is still inside a foreign ball, lower than 3% of its triangle size, or too close to a neighbour is dropped. A vertex no kept seed uses has no ball.
5. **Free sites.** A wall vertex with a ball has no site of its own. Every other site stays 1.5 radii from each ball centre, at the start and after each Lloyd step (a constrained centroidal Voronoi tessellation; Du, Gunzburger & Ju 2003, *Constrained centroidal Voronoi tessellations for surfaces*). The seeds do not move.
6. **Where.** Seeds are placed on the wall triangles that have an edge at which columns end (the neighbour is not coplanar on every layer surface, or is a port) and no corner on a sharp convex wall edge (dihedral above 40° on the innermost surface). On flat walls the units merge into columns anyway. At a sharp convex edge the seeds of the two sides would each have to stay under the plane that halves the corner, and most do not (784 of 903 crease triangles on Manifold had their seed inside a ball of the other side).
7. **Anchored vertex sites.** The site of a wall vertex without a ball at an edge where columns end is kept on the line its layer nodes were marched along, in the Lloyd steps too. The cells then cut the wall triangles there like the Voronoi dual of the wall mesh, each triangle in up to three parts around a common point.
8. **Small columns.** A column whose footprint is below 5% of its wall triangle joins the neighbouring column in the same plane with which it shares the longest edge, as in cell merging for cut cells (Clarke, Salas & Hassan 1986; Quirk 1994; Ye, Mittal, Udaykumar & Shyy 1999). The merged column's innermost cell touches two core cells. A merge that fails the contact checks of step 4 above is undone. This is the safety net for the triangles without a seed.
9. **Torn contacts.** Two cells that both reach round a reentrant wall edge can touch on either side of it with a third cell between them at the edge. FoehnCraft allows one face per cell pair. Such a cell is now treated like a cell behind its wall: its site is replaced by two (root repair), or the part is split along the wall planes.

### Benchmark

Release build, Linux x86_64, ten Lloyd iterations, FoehnCraft `foehn mesh` (main of 2026-11-23). All rows have 0 negative volumes and the exact volume. "Columns before" on Manifold is the retried layout of PR #28.

| Case | Variant | Cells | Inverted | Max / mean non-orth. | Faces > 70° | Max skewness | Max volume ratio | Smallest cell |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Elbow | columns before | 1,082 | 0 | 80.1° / 7.6° | 2 | 0.91 | 727 | 5.4e-6 |
| Elbow | wall snapping (rejected) | 1,079 | 0 | 80.5° / 7.5° | 2 | 0.87 | 277 | 3.7e-5 |
| Elbow | anchored vertex sites only | 1,100 | 0 | 64.8° / 8.9° | 0 | 0.56 | 225 | 4.2e-5 |
| Elbow | seeds on every wall triangle | 1,431 | 0 | 84.9° / 15.3° | 12 | 5.52 | 5.3e5 | 1.2e-8 |
| Elbow | seeds on flat walls too, not at convex corners | 911 | 0 | 81.9° / 21.5° | 18 | 1.22 | 8.6e4 | 1.7e-7 |
| Elbow | chosen rule, free sites at 1.02 radii, no merge | 1,037 | 0 | 68.4° / 8.6° | 0 | 1.09 | 3,111 | 3.6e-6 |
| Elbow | chosen rule, no merge | 1,013 | 0 | 56.6° / 9.4° | 0 | 0.73 | 519 | 2.5e-5 |
| Elbow | **chosen: merge below 5%** | 1,004 | 0 | 56.6° / 9.4° | 0 | 0.73 | 77 | 8.6e-5 |
| Elbow | merge below 15% | 995 | 0 | 56.6° / 9.2° | 0 | 0.73 | 28 | 2.3e-4 |
| Manifold | split cap | 224,028 | 0 | 89.1° / 13.4° | 1,931 | 2.91 | 285 | 5.6e-3 |
| Manifold | columns before | 505,937 | 0 | 87.9° | 412 | 2.59 | 1.1e13 | 8.9e-13 |
| Manifold | wall snapping (rejected) | 474,077 | 0 | 89.4° / 9.1° | 190 | 10.9 | 7.0e11 | 1.4e-11 |
| Manifold | small-column merge only (5%) | 433,533 | 0 | 88.6° / 9.5° | 1,118 | 2.92 | 5.1e5 | 5.2e-6 |
| Manifold | anchored vertex sites only | 502,099 | 1 | 96.3° / 9.8° | 541 | 2.55 | 5.7e12 | 2.2e-12 |
| Manifold | seeds on every wall triangle, no merge | 262,401 | 0 | 78.4° / 6.4° | 14 | 4.81 | 9.9e6 | 8.8e-8 |
| Manifold | chosen rule, no merge | 239,551 | 0 | 73.6° / 5.8° | 2 | 1.44 | 8.3e6 | 3.2e-7 |
| Manifold | **chosen: merge below 5%** | 239,347 | 0 | 74.5° / 5.8° | 3 | 1.44 | 1,882 | 1.1e-3 |
| Manifold | merge below 15% | 239,116 | 0 | 74.4° / 5.8° | 3 | 1.44 | 1,077 | 1.5e-3 |

Per cell on Manifold (chosen): the volume ratio is 2.7 at the median, 27 at the 90th percentile and 170 at the 99.9th; 18 cells are above 1,000. Before, the median was 32 and 55,900 cells were above 1,000. With split caps the median is 1.5, the 99.9th percentile 150. Six cells are above 70° (split cap: 2,698 cells, 100 of them above 85°).

On Manifold 46,635 wall triangles are offered, 46,571 get a seed after 1,383 flips, and 43,774 triangles keep their native prisms (before: none of the curved ones). The build is valid at the first layout (no retry) and takes 130 s (before about 9 minutes with two retries). The maintainer chose the 5% merge threshold and switched the Manifold fixture to columns on 2026-11-23.

### Candidates measured and not kept

- **Wall snapping.** Near incidences of the diagram with wall vertices and edges were made exact by cutting in a copy of the core tetrahedra with displaced wall vertices and mapping the pieces back. It removed 9,638 of 12,486 slivers on Manifold, but 5,269 incidences could not be made exact, the map stretched cells (skewness 10.9) and the worst angle rose. About 850 lines; removed.
- **Anchored vertex sites alone.** Good on the elbow's right angles. On Manifold the three vertex cells of a wall triangle meet outside the triangle in 55% of the triangles (obtuse triangles, unequal heights), which leaves the slivers.
- **Seeds at sharp convex corners and on flat walls.** See the rows above. On flat walls a seed stops the merge of units into wide columns, and the needle triangles at box edges give skewed prisms (4.81).
- **Free sites just outside the balls** (1.02 radii). At a wall vertex where only some triangles have a seed, a free site next to the ball leaves the seeds a slice of the unseeded triangles. 1.25 radii gave 52.0° but volume ratio 2.2e6 on the elbow, 2.0 gave 60.4° and 4,910; 1.5 is kept.
- **Radii from a weight solve.** Vertex weights fitted by projection onto "seed high enough" and "seed outside the ball across each edge" stalled with about 2% of the triangles violated (Python prototype); simultaneous updates diverged. Radii from the circumradii and flips need no solve.

### The planarity check and thin faces

Seeds round a wall vertex are all on its ball. Where they are nearly on one circle, cells that are not neighbours on the wall touch in a long thin face over the vertex (2.6 long and 7e-5 wide in the case that stopped the Manifold build). Its vertices were on one plane to 1e-13, but the Newell normal summed from absolute coordinates (263 from the origin) was tilted by 1e-10, and the offset test called the face nonplanar. `ring_plane` now sums about the first vertex, and the planarity tests of the shell check and of the core face merge measure offsets from the least-squares plane through the vertex mean (`fitted_normal`). A bent face is rejected as before.

## Limits

- Sharp convex wall edges and flat walls have no seeds. The remaining tail sits there and at wall vertices where only some triangles have a seed: 18 cells above volume ratio 1,000 on Manifold.
- The flips change the seeds, not the wall mesh: a flipped pair of wall triangles is cut along the other diagonal (two halves each; 1,383 pairs on Manifold). A Delaunay flip pass in the wall remesher would give these whole cells too.
- A seed is dropped where two walls are closer than about the wall mesh size (the ball of a vertex across the gap holds it). Those triangles are cut as before.
- The cells of the seeds are low cells over one wall triangle. The core sphericity on Manifold is 0.74 at the median (split cap: 0.85).
- The repair of torn contacts is not reached on the two fixtures with the chosen layout (it was needed by the anchored-only layout on Manifold); no test covers it.
- `mesh.estimate` counts two layer cells per wall triangle and layer for columns (built: 0.70 on the elbow, 0.98 on Manifold), with `cells_basis` "approximate".
- On curved walls a layer cell covers at most one wall triangle, so the join there is only as good as the scaffold's wall triangulation is fine against the core.
