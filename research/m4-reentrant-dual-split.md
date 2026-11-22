# Inverted faces of the polyhedral duct core at reentrant wall edges

Design note of ALG-16 for the bug that made FoehnCraft reject the polyhedral elbow (M1b, FC-06). Question: why did `foehn mesh` find interior faces with S·d ≤ 0 in the volume median dual, and what removes them? The fix is in [`duct/volume_dual.rs`](../crates/mesh-algorithms/src/duct/volume_dual.rs) (`split_reentrant`), algorithm id `volume_median_dual_v2`.

## Root cause

FoehnCraft calls an interior face inverted when its area vector S does not point from the owner's centroid to the neighbour's (S·d ≤ 0, d between the two cell centroids). Before the fix the elbow had 10 such faces and SolidCraft 18. All 28 have the same shape: a thin wall prism of the innermost layer on one side and the median-dual cell of a wall vertex on the other, and that vertex sits on a **reentrant** edge of the core boundary, where the offset walls meet at more than 180° seen from the core (the inner corner of the elbow at x = 2, y = 1, and the box-to-pipe junction of SolidCraft at x ≈ 180.36).

The median dual of a vertex is the union of its corner pieces of the incident tetrahedra. At a reentrant edge the coarse scaffold makes this cell L-shaped: one arm along each wall, about 0.33 long and 0.1 thick on the elbow. Its centroid lies between the arms, on the solid side of both wall planes near the vertex. For a prism end face on one arm, d then points back across the face. Measured on the elbow (cell 6774): its centroid is at y = 0.919 while the wall faces lie on y = 0.869 with the core below; the worst face had cos(S, d) = −0.41.

The dual faces themselves are not at fault: each lies in the plane where two barycentric coordinates are equal, so S·(x_b − x_a) > 0 always holds for the two vertex positions. Only the centroids, which FoehnCraft uses, break it. Cells at convex edges and interior cells were never affected.

## Candidates

| | Method | Source |
|---|---|---|
| Classical | Split the vertex-centred dual cell along scaffold faces through the vertex, one part per wall plane | Standard repair of vertex-centred finite-volume duals at sharp boundary features (cells split at corners so that each part has one boundary normal) |
| Graphics (OBJ-19) | Replace the median dual by a clipped Voronoi diagram of the scaffold vertices, relaxed to a centroidal Voronoi tessellation | Yan, Wang, Lévy & Liu (2013), *Efficient computation of clipped Voronoi diagram for mesh generation*; Du, Faber & Gunzburger (1999) for CVT. In a CVT each site is its cell's centroid and every face is the bisector plane of its two sites, so S·d > 0 holds by construction |

Lévy's clipped Voronoi with CVT is the direction OBJ-20 already sets for the polyhedral core, and it would make the property structural instead of repaired. It cannot be benchmarked on these fixtures yet: KoflCraft's clipped Voronoi (ALG-04) clips by a convex shell only, and clipping by the nonconvex elbow and SolidCraft walls is the open M4 work itself. Crane's and Kobbelt's work offered no method specific to this repair. The classical split is local, keeps the scaffold, the wall prisms and the exact prism-to-core face matching, and was built here. No source code was read; both methods are from the papers and from first principles.

## The split

For every wall vertex cell with an inverted wall face (S·d ≤ 0 against the innermost prism, computed with the same centroid formula as FoehnCraft):

1. Group the cell's boundary triangles into wall planes: a triangle joins the first plane whose inward normal is within 30° of its own.
2. Give each corner piece (one per incident tetrahedron) to the plane it stands farthest in front of: the largest r·n, with r from the vertex to the mean of the piece's eight points and n the plane's inward normal.
3. Split each group into connected sets through shared tetrahedron faces. The faces between two parts are the corner quadrilaterals on the tetrahedron faces through the vertex, which are planar.
4. Keep the split only if every part is a closed manifold shell and the number of inverted wall faces falls; otherwise the cell stays whole. A part's kernel point is a short step from the vertex in front of every face through it, as for the whole cell before.

Volumes, patch areas and the matching of all prism end faces are unchanged; the shell checks of `verify` run on the parts as on any cell.

## Benchmark

Fixtures as in `docs/nonconvex-duct.md`, release build, read back with FoehnCraft `foehn mesh` 0.7.0 (main of 2026-11-22), Linux x86_64. Non-orthogonality is the angle between S and d as FoehnCraft reports it.

| Case | Version | Cells | Inverted faces | `foehn mesh` | Max non-orthogonality | Faces > 70° | Max skewness | Min / p05 / median sphericity |
|---|---|---:|---:|---|---:|---:|---:|---|
| Elbow polyhedral | before (v1) | 6,995 | 10 | invalid | 114.0° | 57 | 74.8 | 0.522 / 0.570 / 0.764 |
| Elbow polyhedral | split (v2) | 7,005 | 0 | valid | 87.0° | 37 | 2.72 | 0.476 / 0.566 / 0.762 |
| SolidCraft polyhedral | before (v1) | 417,454 | 18 | invalid | 103.2° | 7,035 | 5.12 | 0.344 / 0.695 / 0.833 |
| SolidCraft polyhedral | split (v2) | 417,457 | 0 | valid | 89.6° | 7,012 | 5.12 | 0.344 / 0.695 / 0.833 |

The elbow splits 10 cells into two parts each, SolidCraft 3. The MSH 4.1 subdivisions of both load valid as well. Build time of the SolidCraft polyhedral mesh stays about 19 s. The bulk and layered modes have no dual cells and were valid before.

Outcome: the classical split, since the graphics candidate cannot run on nonconvex domains yet. The cost is roundness: the elbow's least round core cell drops from 0.522 to 0.476 sphericity, and its test floor from 0.50 to 0.45 (maintainer decision, 2026-11-22).

## Limits

- The split repairs centroids; it does not make the cells orthogonal. Non-orthogonality still reaches 87° and 89.6°, so some faces pass S·d > 0 with little margin. A finer scaffold at reentrant edges, or the CVT core above, would raise it.
- Only faces to wall prisms trigger a split. Faces between two dual cells are counted in the hand-off evidence (`inverted_faces`) and checked by CI, and are zero on both fixtures, but no rule repairs them.
- `mesh.estimate` counts the cells from the scaffold without building the wall layers, so it no longer sees the split: it reports the count before it (10 below the built elbow, 3 below SolidCraft) and says so in `cells_basis`. Building the dual for an exact count would cost about 7 s on SolidCraft.
- A cell on a reentrant edge with a single smooth (concave) wall plane cannot be split by planes; none occurs on the fixtures.

## Regression

`nonconvex_elbow_and_solidcraft_preserve_volume_and_expose_each_port_layer` asserts `inverted_faces = 0` in `mesh-evidence.json` for all six duct meshes, and `cargo m1 --require-foehn` now reads the SolidCraft ducts as well as the elbow back with `foehn mesh` and fails on any export FoehnCraft calls invalid.
