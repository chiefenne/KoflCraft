# Triangle split (tri-center) for block topologies (MB1h)

Date: 2026-11-25. Requirements: BS-03, BS-05, BS-13. Follows
[`block-layout-optimisation.md`](block-layout-optimisation.md) (MB1g).

## Task

The maintainer asked for a "tri-center" split: a block whose cross-section is
a triangle becomes three quadrilaterals around a centre of valence 3, and the
same exists in 3D (a prism becomes three hexahedra, a tetrahedron four). It is
meant for the block-structured mesher in general, not only for the
multi-element airfoil, and an automatic criterion should decide where a split
goes in, so that the cell-size variation and the quality stay under control.
On the scope card the maintainer chose "Also tetrahedron to 4 hexes": a general
split with an automatic trigger for 2D and prism blocks, plus the tetrahedron
split as a building block for 3D templates, in one pull request.

## Sources

The constructions are standard, and no source code was read for them.

- The split of an n-gon into n quadrilaterals through its edge midpoints and
  a face point is the subdivision step of Catmull and Clark ("Recursively
  generated B-spline surfaces on arbitrary topological meshes", CAD 10(6),
  1978). For a triangle it is the tri-center map.
- The split of a tetrahedron into four hexahedra through its edge midpoints,
  face centroids and body centroid is the simplest all-hex conversion. It is
  described as the baseline in the survey of Pietroni, Campen, Sheffer,
  Cherchi, Bommes, Gao, Scateni, Ledoux, Remacle and Livesu ("Hex-mesh
  generation and processing: a survey", ACM TOG 42(2), 2022).
- The new count of the triangle enters as a new row of cells that runs through
  the neighbouring blocks to the boundary. In hexahedral meshing this is sheet
  insertion (Merkley, Ernst, Shepherd and Borden, "Methods and applications
  of generalized sheet insertion for hexahedral meshing", 16th International
  Meshing Roundtable, 2007). Here it is limited to one layer of blocks
  (extruded topologies) and to chains that `block.split_block` already
  handles.
- The trigger judges a split by the layout optimiser of MB1g (after Sauer and
  Morsbach, IJNME 2023, doi:10.1002/nme.7308) and by the same measure: the
  finite-volume non-orthogonality and growth of the TFI grid.

## Counts

Each corner i of the triangle (or tetrahedron) carries a count n_i. The piece
at corner i has n_j cells along its side toward corner j, so the edge (i, j)
gets n_i + n_j cells. The three side counts of a triangle therefore sum to an
even number, and each one is below the sum of the other two. Because the counts
sit on the corners, neighbouring triangles or tetrahedra with the same corner
counts meet node to node.

In an existing block with its flat corner M on side AB, the two edges at M
fix n_B (A–M) and n_A (M–B). The opposite edges B–C and C–A keep their
counts. n_C is new. By default it is the distance from M to the centroid
divided by the mean spacing at M, rounded, and at least 1.

## What was built

- `block.split_triangle` (engine `block_triangle::split_triangle`, edit
  `SplitTriangle`) splits a block at a corner of at least 120°. The row of n_C
  cells leaves through side BC next to B and side CA next to A, and cuts every
  block on its way in two. New vertices sit at node n_C of each edge's own
  curve and law, refined to n + n_C cells. The halves get fitted curves and
  laws, and the cut faces follow the old grid. Patches and projection
  surfaces keep their faces. A split is refused in these cases: a neighbour
  across the k faces (more than one layer), a row that would enter a block from
  both ends, cross it twice or come back into the triangle, and a row through a
  block with a face on a projection surface.
- `block.optimize_layout` runs the MB1g optimiser on any extruded topology as
  a tool, with the trigger (`triangle_splits`, default true).
- The trigger (`inject_triangle_splits`) takes the block corners of at least
  150° in the layout before optimisation as candidates, flattest first, at
  most four per round. Each trial splits one candidate, optimises the layout
  and compares it with the layout optimised without the split. It keeps the
  split only if all of these hold:
  - the worst angle of the triangle's blocks drops by at least 2°;
  - the whole grid's worst angle rises by no more than 0.05°;
  - no more cells pass 70°;
  - nothing folds;
  - the growth stays within the larger of 2 and the start growth, plus 5 %.

  The best trial of a round is kept, the next round starts from it, and it
  stops after four splits. Every corner it tries and leaves is reported with
  the reason.
- Growth hold (`LayoutOptions::hold_growth`, MB1g optimiser). Without it, the
  optimiser traded a cell size ratio of 3.3 across the new row for the angles
  in the first test. With it, each face's growth may reach the larger of 2 and
  its start growth. Past that limit the cost rises steeply, about 4,000 at 25 %
  over, which is as much as a folded corner. The hold is on for the trigger's
  trials and off for a start layout that folds (unfolding needs free sizes).
  The multi-element template keeps the MB1g objective.
- Building blocks `block::tri_center` (mesh-algorithms):
  - `TriCenter::prism` gives three hexahedra from a triangle and its top.
  - `TriCenter::tetrahedron` gives four hexahedra; edge and face points are
    shared between neighbours, and either corner orientation works.
  - Templates `triangle` and `tetrahedron`, for viewing and tests.

## Measurements

| Case | Max non-orthogonality | Cells | Notes |
|---|---|---|---|
| Wedge, a 171° triangle beside a quad, open boundaries | 64.0° raw, 61.9° optimised → **26.2°** with the split | 144 → 274 | Growth 1.76 → 2.00; the second round tried a 164° corner of a new block and left it ("less than 2° better"); 0.2 s |
| Same wedge, without the growth hold | 26.1° | 274 | Growth 3.27 across the row: refused by the growth guard |
| Same wedge with walls (blocks rigid) | 64.0° → 58.0° | 274 | Kept: the walls hold the triangle's blocks |
| Nine multi-element cases of MB1g | 47.4–68.2°, unchanged | unchanged | No candidate of 150° in any case; in cases A and B the flattest corners are 127–135° inside the wall rings and at most 125° elsewhere |

The multi-element template from MB1g already has no triangle blocks. The rings
and the cove turned the 174° gap corner of the first layout into ordinary
quadrilaterals. So the trigger runs there, finds nothing and keeps the layout.
`block.template` reports that under `layout.triangle_splits`. The overlapping
slat that folds (10°, gap 0.02, overlap −0.005) has no flat corner either. Its
fold needs a slot topology.

Two variants were tried and dropped:

- A steep growth cost on every face, past growth 2, in the default objective.
  Six of nine airfoil cases then folded, because their start layouts have
  large growth next to the rings.
- The growth hold as the default. Case F rose from 64.6° to 84.4° (4 cells
  over 70°), because its start layout needs sizes to change.

## Limits

- A split works on one layer of blocks (extruded topologies), like the
  layout optimiser. In 3D the tetrahedron split is a template building block,
  not an edit.
- A corner of exactly 180° cannot be stored: the topology check refuses the
  block as not right-handed. A triangle block needs a slight bend at M.
- The trigger tries at most four candidates per round. Each trial costs one
  layout optimisation: 0.2 s for the wedge, 0.5–2.5 s for an airfoil
  layout.

## 3D: the split through a column (2026-11-27)

Step 1 of the 3D insertion plan (maintainer cards of 2026-11-27: conformal
first, integer count solve, geometry plus measured trial as trigger).

- **The flat corner is a block edge.** A prism-shaped block has an edge whose
  cross-section corner is flat (150° or more at both ends, `flat_edges`). The
  split takes that edge as the prism axis (`along`, or the flattest edge at
  `vertex`). The block is turned into a local frame with the edge along k
  (one of the 24 right-handed rotations of the hex frame), so the 2D code of
  MB1h applies unchanged in the cross-section.
- **The column.** Every block stacked on it across the faces at the two ends
  of the edge is turned to the same frame (face corner by face corner) and
  split the same way, so the faces between the layers stay point-matched
  (BS-05). A column that closes into a ring is refused for now.
- **The sheet.** The new count n_C enters on the edges B–C and C–A of every
  layer of the column and runs through the neighbouring blocks to the
  boundary, the same propagation as before, now without the one-layer limit.
- **Why a sheet is always needed.** In a hex replaced in place, opposite edges
  carry the same count (N_BC = N_AM, N_CA = N_MB), so the tri-center counts
  give n_C = (N_BC + N_CA − N_AB) / 2 = 0: the new count cannot fit the old
  ones. The chain through the block must be split by a sheet (now) or stopped
  at a T-junction (step 2, with the integer count solve, after T-mesh
  quantization in Kobbelt's and Campen's groups). The block model shares
  faces whole today, so T-junctions need support in the topology and the
  grid builder first.
