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

## Next

PR 3b uses it in `block.split_tet` and `block.split_triangle` (the
interface record, nodes on the parent face by inverse TFI, polyhedral
hand-off and writers, FoehnCraft round trip); PR 3c the trigger in the
automatic-fallback mode, the report and requirements BS-16.
