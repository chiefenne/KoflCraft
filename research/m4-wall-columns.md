# Wall columns under the Voronoi core (M4)

Design note of ALG-16 for the join between the prism wall layers and the Voronoi core of [m4-voronoi-core](m4-voronoi-core.md). Question: can the layer cells follow the core cells (one column per core cell) on curved walls too, with exact volumes and no inverted faces, and does that lower the worst non-orthogonality, which sat at the join (85.9° on the elbow, 87.7° on Manifold)? The code is [`duct/columns.rs`](../crates/mesh-algorithms/src/duct/columns.rs). It is the default join; `"wall_join": "split_cap"` in `volume_scaffold` keeps the split cap of PR #22.

## Why the first columns failed on curved walls

The first column variant of PR #22 extruded each core cell's wall polygon through all layers and merged it over the wall triangles up to a 30° feature angle. FoehnCraft integrates every face as a fan about its node average (`polygon_area_centroid` and `cell_volume_centroid_of` in `foehncraft-mesh/src/geometry.rs`, read as documentation of the consumer). A face merged over a curved wall is bent, and its fan surface lies on the chord, off the wall by about w²/8R. On Manifold, with core cells about 6 units wide on pipes of radius 25, that is about 0.18, against a first layer of 0.08. Two things follow. The merged wall face moves the wall, which costs 0.14% of the volume. The fans of the two caps of a thin layer cell cross, which inverts faces. Any cell whose cap spans non-coplanar triangles has this problem, so on curved walls no layer cell may span two wall triangles.

## Construction

1. **Units.** The core cells cut the innermost layer surface into wall polygons. A unit is the polygon (or polygons) of one core cell in one wall triangle. These are the facets of the restricted Voronoi diagram on the wall triangles (Yan, Lévy, Liu, Sun & Wang 2009, *Isotropic remeshing with fast and exact computation of restricted Voronoi diagram*), here taken from the clipped volume diagram.
2. **Extrusion.** A polygon vertex with barycentric coordinates λ on the innermost triangle is placed at λ on the same triangle of every layer surface. A layer prism is affine in λ at each height, so each extruded unit is a union of sub-wedges (a fan of its polygon), these sub-wedges partition the coarse prism exactly, and every cap is planar. A vertex on a triangle edge gets the same image from both triangles. Side faces are bilinear quads, as on a native prism.
3. **Columns.** Units of one core cell over adjacent triangles are merged when the triangles are coplanar on every layer surface: on the wall to 1e-10 of the model size, inside to 1e-6 of the thinnest layer. On a flat wall a core cell then sits on one layer cell per layer and wall plane, with its centroid under the core cell's. On a curved wall each unit stays its own column.
4. **One face per pair.** FoehnCraft allows one face between two cells. A merged column must be one disk, and its contact with each neighbouring column one open chain. A column that fails this is split back into its units. A unit with one polygon is convex, so two such units meet in one chain; if two units still meet in more, the build fails and the engine rebuilds with one more Lloyd iteration.
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

Removing slivers needs a change in the core: no core cell should touch a wall triangle in a sliver. That is follow-up work.

### Other candidates

A planar polygonal wall, built by intersecting tangent planes at the sites as for the planar hexagonal meshes of architectural geometry, would let a column have a flat wall face on curved walls. It changes the wall itself, so the scaffold's volume and areas would no longer be the reference. Not tried.

## Limits

- Slivers on curved walls, as above.
- `mesh.estimate` counts two layer cells per wall triangle and layer for columns (built: 0.76 on the elbow, 2.08 on Manifold), with `cells_basis` "approximate".
- On curved walls a layer cell covers at most one wall triangle, so the join there is only as good as the scaffold's wall triangulation is fine against the core.
