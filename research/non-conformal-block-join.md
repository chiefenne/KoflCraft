# Non-conformal block join: face overlay (BS-16, step 3)

Step 3 of the 3D tetrahedron and prism insertion (plan doc "KoflCraft: 3D
tetrahedron and prism block insertion (plan)", section "Step 3 design",
2026-11-28). A split piece keeps its own cell counts, its neighbours stay as
they are, and KoflCraft joins the two grids on each face between them into
one conformal polyhedral mesh (plan card 2, option B1).

## Decisions (maintainer cards, 2026-11-28)

| Card | Answer |
| --- | --- |
| A When non-conformal is used | automatic fallback: conformal where it works; when a neighbour cannot take prism columns the split joins non-conformally and the report says so; `conformal_only` forbids it |
| B Sliver faces | snap within 2 % of the shortest edge at a point; the largest snap is reported |
| C Count n0 of the piece | match the neighbours' cell size along the interface |
| E Count when size match fails the bar | coarsest n0 from half the largest count up that reaches scaled Jacobian 0.4 |
| D `mesh.msh` | polyhedra split into pyramids and tetrahedra from the centroid, with parents, as the duct mesher writes them |

Reading of card C: the tet keeps the counts a, b, c of its fourth face, so
its edge P0–P1 carries a + n0 cells against the neighbour's a, and the piece
is always finer. The coarsest n0 its shape allows is half its largest count
(the new vertex then sits at most 2/3 of the way from P0, spacing equal on
both sides of it). Measured on the tet on a box it is too coarse for the
quality bar: n0 = 3 gives the tet's cells a worst scaled Jacobian of 0.35,
4 gives 0.43, 5 (the mean) 0.49. So the default (card E, maintainer 2026-11-28 08:25 UTC) is the coarsest n0 from half the largest count up
whose four hexahedra reach 0.4.

## Method

Both sides of a block face are grids on the same transfinite patch, so the
overlay is done in the patch's own parameters (s, t) with no projection.
The kernel (`mesh-algorithms/src/block/overlay.rs`) builds the planar
arrangement of the two edge sets:

1. points of the two tilings closer than the snap distance merge (the more
   constrained one keeps its place: region corner, then border);
2. a point on (or within the snap distance of) an edge of the other tiling
   is put into that edge;
3. edges of the two tilings that cross properly are split at the crossing;
4. the faces of the planar graph are traced around each point by angle;
5. each face takes the cell of each tiling on its left; a face bounded only
   by the other tiling's edges inherits the cell across them.

Every topological decision (collinear, crossing) uses the exact planar
orientation predicate (`predicates::orient2d`, Shewchuk 1997 expansions);
only the coordinates of crossing points are rounded. Each overlap polygon is
one face seen from both cells, so the interface flux is conserved to
round-off. The kernel also returns, for every edge of each tiling, the
overlay points inside it, which the 3D assembly puts into the other faces at
that edge so every cell closes (as `join.rs` does for the sweep cut).

## Sources (papers only, no code read)

- P. E. Farrell, M. D. Piggott, C. C. Pain, G. J. Gorman, C. R. Wilson,
  "Conservative interpolation between unstructured meshes via supermesh
  construction", CMAME 198, 2009: the supermesh, also behind FoehnCraft's
  sliding interface.
- M. Cherchi, M. Livesu, R. Scateni, M. Attene, "Fast and robust mesh
  arrangements using floating-point arithmetic", ACM TOG 39(6), 2020:
  arrangements with exact predicates and rounded constructions.
- J. R. Shewchuk, "Adaptive precision floating-point arithmetic and fast
  robust geometric predicates", DCG 18, 1997.

## Why not the alternatives

- Hanging nodes 2:1 (B3) fix the ratio, so the piece loses its free count.
- Paired patches (B2) need a general static imprint in the solver, which
  FoehnCraft does not have; it can come later as an export option from the
  same overlay.
- A 3D supermesh by projection (FoehnCraft D58) carries a projection error
  that a shared parameter domain avoids.

## Results of the kernel (tests in `overlay.rs`)

| Case | Faces | Area error |
| --- | --- | --- |
| 3 × 3 against 4 × 5 grid | 42 (6 columns × 7 rows) | < 1e-12 |
| 2 × 2 against 4 × 4 (shared lines, exact) | 16, no crossing points | 0 |
| sheared 5 × 4 against 7 × 3 | per cell conserved | < 1e-12 |
| tri-center quads against 7 × 6, snap 2 % | per cell conserved | < 1e-9 |
| lines 1e-5 apart, exact / snap 2 % | 9 / 6, smallest 2-area < 1e-4 / > 0.1 | < 1e-12 |

## The 3D assembly (PR 3b, `mesh-io/src/block_join.rs`)

An interface in the topology names two patches: `coarse` is one whole block
face (the neighbour's), `fine` the tet's three parts of it. Each side is
gridded on its own as boundary faces. The coarse face's grid gives the
parameters (node (i, j) at (i / n_i, j / n_j)); fine nodes are placed there
by inverting the coarse bilinear quads (Gauss–Newton per quad, nearest
kept), fine nodes on the border are put exactly on it. A crossing point is
put on its coarse quad, so the coarse cells keep their shape. Each
overlap polygon replaces the two cell faces it lies in; the points inside
every edge are put into every other face of the mesh at that edge (side
faces of both layers, boundary faces along the face's border). Merged fine
nodes are aliased to their coarse twins everywhere.

| Fixture | Cells | Polyhedra | Interface faces | Inverted | Volume (poly / hex − 1) |
| --- | --- | --- | --- | --- | --- |
| tet on a box (bottom not flat, so no prism column) | 462 | 57 | 113 (20 coarse, 47 fine, 45 crossings) | 0 | < 1e-9 |

FoehnCraft (`foehn mesh`, main of 2026-11-28) reads `foehncraft.vtu` as
valid: 462 cells, 0 inverted faces, closure 1.3e-16. It reads `mesh.msh` as
1,676 cells (the pieces), same volume, 0 inverted faces.

## Next

PR 3b uses it in `block.split_tet` and `block.split_triangle` (the
interface record, nodes on the parent face by inverse TFI, polyhedral
hand-off and writers, FoehnCraft round trip); PR 3c the trigger in the
automatic-fallback mode, the report and requirements BS-16.
