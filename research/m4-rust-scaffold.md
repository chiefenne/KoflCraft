# Wall fidelity and a Rust tetrahedral scaffold (M4, OPEN-07)

Design note of ALG-16 for the volume scaffold of the polyhedral duct. Manifold is the duct test part formerly named SolidCraft (`test-fixtures/solidcraft-vs-foehncraft`). Question from the maintainer (2026-11-23): the Manifold STL is fine enough, so why does the polyhedral mesh look as if it were built on a coarser surface, although its cells are finer? The maintainer then asked that the scaffold no longer depend on Gmsh. Code: [`remesh.rs`](../crates/mesh-algorithms/src/remesh.rs) (wall), [`predicates.rs`](../crates/mesh-algorithms/src/predicates.rs) and [`delaunay.rs`](../crates/mesh-algorithms/src/delaunay.rs) (tetrahedralization), [`scaffold.rs`](../crates/mesh-algorithms/src/scaffold.rs) (scaffold), configured by `volume_scaffold.build` in the partition. The scaffold source is reported as `rust_remesh_delaunay_v1`.

## Diagnosis

The polyhedral mode takes its wall from the scaffold's boundary triangles, not from the STL. `walls.vtu` of the earlier mesh is exactly the scaffold wall: 21,694 triangles, 11,327 nodes. The Gmsh preparation (`CreateGeometry` on discrete charts, `MeshSizeMax` 6.25, curvature sizing off) had remeshed the 1-unit STL wall (median edge 1.0) into triangles with a median edge of 5.9 units. The finer core and layer cells subdivide those flat triangles.

Measured with sampled closest-point distances (each facet at its centroid, edge midpoints and three interior points):

| Wall | Triangles | Scaffold → STL max / p99 | STL → scaffold p99 |
|---|---:|---:|---:|
| Gmsh scaffold (previous fixture) | 21,694 | 0.543 / 0.227 | 0.411 |
| Gmsh, `MeshSizeFromCurvature = 48` | 22,500 | 0.631 / 0.214 | 0.320 |
| Gmsh, chord size field (0.05) as a background view | 30,222 | 0.723 / 0.183 | 0.181 |
| STL triangles as the wall | 89,872 | 0 | 0 |
| Rust remesh, tolerance 0.05 | 56,606 | 0.050 / – | 0.050 (max) |

The nodes of every Gmsh wall lie on the STL; the error is the chord between them. The largest deviations are on the pipe near the box junction (x ≈ 179) and on the inlet pipe, about seven times the first prism layer (0.08). Gmsh's curvature sizing does not act on these discrete charts, and with a chord size field the reparametrisation of the long pipe charts still stretches triangles to 12 units. The STL itself cannot be the wall: its 28-unit needle triangles invert a scaffold tetrahedron at the wall offset (MESH006).

## Wall: adaptive isotropic remeshing

Candidates: the Gmsh surface remesh (classical, above) and the remeshing loop of Botsch & Kobbelt, *A Remeshing Approach to Multiresolution Modeling* (SGP 2004), with the curvature-adaptive target length of Dunyach, Vanderhaeghe, Barthe & Botsch, *Adaptive Remeshing for Real-Time Mesh Deformation* (Eurographics 2013) (OBJ-19: Kobbelt). No source code was read (LIC-12 to LIC-15).

- Target length L = sqrt(6ε/κ − 3ε²): an equilateral triangle of edge L on a sphere of curvature κ deviates from it by ε. κ is the largest principal curvature of the normal-cycle tensor of the STL (GRD-04, `sizing::surface`), clamped to `[min_size, max_size]`, plus `rim_size` at the port rims, gradient limited with `rate` (branch-and-bound envelope of GRD-04). Targets are cached per vertex.
- Ten iterations of: split edges above 4/3 L at the midpoint, collapse edges below 4/5 L (link condition, no normal flip, no new edge above 4/3 L), flip edges towards valence six (four at a rim), tangential area-centroid smoothing, projection onto the STL wall. A final split pass removes long diagonals the last flips left.
- Constraints: port facets are kept as they are, with their rim vertices and edges. STL edges sharper than `feature_angle` (30°) are feature lines: split along the line, never flipped, collapsed only along the line; feature vertices move along their line and project onto the STL feature edges; corners (feature valence other than two, or a turn above 45°) stay. New facets may not be worse in shape (4√3·area/Σl²) than 0.1 or the facets they replace.
- Deviation pass: every facet whose sampled distance to the STL exceeds ε has its longest edge split, up to eight rounds, so the reported bound is measured, not estimated.

Manifold, ε = 0.05: 56,606 wall triangles, edges 0.25 to 8.0 (median 2.7), largest distance 0.050 both ways, 252 feature edges and 8 corners kept, 11 s.

## Volume: Delaunay tetrahedralization with recovery and refinement

Candidates: Gmsh (classical, separate process, LIC-09) and a Rust pipeline from the papers:

1. Bowyer–Watson insertion (Bowyer 1981; Watson 1981) of the wall vertices in Morton order, located by a visibility walk. The cavity grows until every boundary face sees the new point strictly from the front, so cospherical input gives no flat tetrahedra. Orientation and in-sphere signs are exact: floating point first, exact expansion arithmetic when the result is within a generous error bound (Shewchuk 1997).
2. Boundary recovery by conforming Steiner points (Murphy, Mount & Gable 2001; Cohen-Steiner, Colin de Verdière & Yvinec 2002): a missing wall edge is split at its midpoint, a missing triangle has its longest edge split. The points lie on the wall triangles, so the wall shape and its deviation do not change. Manifold needs 582 such points.
3. Tetrahedra reachable from the bounding box without crossing a wall triangle are outside.
4. Delaunay refinement of the interior (Shewchuk 1998): the circumcentre of a tetrahedron larger than the size field (wall spacing grown at `rate`, capped at `max_size`) or with a radius-to-shortest-edge ratio above 2 is inserted if it is inside and at least half a spacing from the wall. Insertions never cross a wall triangle.
5. Slivers. Most of the flat tetrahedra have all four corners on a curved wall: two wall triangles and the opposite diagonal. Their circumspheres bulge outward, so no interior point removes them, and points just inside the wall only add points (8,000 tried, the slivers stayed). They are removed by edge removal around their diagonal (Shewchuk 2002, *Two Discrete Optimization Algorithms for the Topological Improvement of Tetrahedral Meshes*) with the best triangulation of the ring found by dynamic programming (Klincsek 1980). Interior slivers get their circumcentre or a random point near it (Li & Teng 2001).

| Manifold scaffold | Points | Tetrahedra | Worst / p1 / median mean ratio | Below 0.2 | Volume vs STL | Time |
|---|---:|---:|---|---:|---:|---:|
| Gmsh 4.12.1 (previous fixture, optimised) | 26,962 | 119,124 | 0.257 / 0.528 / 0.886 | 0 | −0.183 % | separate step |
| Rust, before sliver removal | 49,763 | 219,116 | 0.002 / 0.124 / 0.771 | 2,851 | −0.028 % | 8 s |
| Rust, with edge removal | 50,384 | 219,169 | 0.052 / 0.32 / 0.77 | 48 | −0.028 % | 11 s |

The Gmsh scaffold is smoother (Gmsh also relocates vertices); the Rust scaffold is adequate for the Voronoi core, which moves the sites anyway (ten Lloyd iterations), and has twice the points because its wall is 2.6 times finer. Vertex smoothing (ODT, Alliez, Cohen-Steiner, Yvinec & Desbrun 2005) is the obvious next step if the scaffold quality matters for the median-dual core.

## Result on the polyhedral mesh

Manifold, polyhedral mode, Voronoi core, split-cap join, release build, Linux x86_64:

| Scaffold | Cells | Wall faces | Inverted (hand-off) | Wall → STL max | Volume vs STL | Build |
|---|---:|---:|---:|---:|---:|---:|
| Gmsh fixture | 92,058 | 21,694 | 0 | 0.543 | −0.183 % | 75 s |
| Rust, ε = 0.05 | 223,216 | 57,630 | 0 | 0.050 | −0.028 % | 188 s |

The finer wall needed one fix outside the scaffold: the verifier treated a sliver triangle face of a core cell (two corners 1.3·10⁻⁸ apart) as nonplanar, because its Newell normal is too inaccurate for the offset test. A triangle is planar, so the test now skips triangles.

## Limits

- The ports are kept as they are. A port with few, large triangles (the elbow's unit squares) keeps coarse cells at its rim; remeshing ports in their plane is a follow-up. The elbow keeps its Gmsh scaffold for now.
- The mean-ratio floor of the Rust scaffold (0.05) is lower than Gmsh's (0.26).
- Feature lines are found by the dihedral angle only.
