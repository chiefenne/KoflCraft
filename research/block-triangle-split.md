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

## T-junctions and the count solve (step 2, 2026-11-27)

Step 2 of the 3D insertion plan, part 1 of three pull requests (maintainer
card "Three PRs"). Requirement BS-21.

**What a T-junction can and cannot do.** In a mesh whose cells meet node to
node, every dual sheet (a layer of hexahedra) ends on the boundary or closes
on itself. So the n_C new cells of a triangle split cannot stop at a
T-junction: they run along the neighbours' chain as extra count. What a
T-junction stops is the block cut. With the sheet split of step 1 every
block on the row's way is cut in two; with a T-junction the neighbour keeps
its block, a corner of the new pieces sits on its edge, and its grid spreads
the extra cells by its own law. The maintainer chose this reading on a card
("Stop block cuts, cells still run"), against closing the row around the
triangle's corner (a pillow, local but with irregular vertices).

**Model.** The topology lists `junctions`: an edge a–b and the vertices it
passes through. The check requires every piece to be an edge of a block and
the counts to add up (whole = sum of pieces). A block face of one block, in
no patch, whose region the faces of other blocks tile in strips (cuts from
one side to the opposite one) is a coarse face; the strips are its
sub-faces. The grid builder grids the coarse faces first (the whole edges
are the pieces' nodes joined, the interior by TFI), stores the sub-faces and
cuts as parts of that grid, and the blocks then read them like any shared
face. The hexahedra therefore conform without welding. This is the T-mesh of
layout quantisation (Lyon, Campen and Kobbelt, CGF 40(2), 2021; Brückler,
Bommes and Campen, ACM TOG 41(4), 2022), limited to strip covers.

**Count solve.** A sheet is a class of block axes with one count (the chains
of BS-05). Without T-junctions each sheet is free. Each junction adds the sum
rule count(whole) − Σ count(pieces) = 0 over sheets, so the tied sheets lie on
the integer kernel of that matrix. `block::counts::fit` takes an integer
basis of the kernel by column Hermite reduction (size-reduced pairwise), the
least-squares point of the relative deviations (count − target) / target on
it, rounds it and improves it by single and paired basis moves, with every
count at least 1. The one-integer-per-sheet form follows the advice of the
3D blocking survey (integer-sheet-pump, Brückler et al., SGP 2024). No code
was read.

**Triangle split with `stop: t_junction`.** The column is split into its
pieces as before; the row edges B–C and C–A of every layer get their new
vertex. Each of those edges that a neighbour still has whole becomes a
junction; no neighbour is cut. The count solve then starts from the counts
before the split (the pieces' own counts for them) and makes every sum add
up.

| Case | Blocks | Cells | Count solve |
|---|---|---|---|
| Two-layer wedge (171° triangle + quad, 2 layers), `cells` 3, sheet split | 10 | 150 | – |
| Same, `stop: t_junction` | 8 | 126 | 3 junctions, 7 sheets (3 tied); the quad's chain goes 4 → 6, n_C 3 → 2 (largest change 50 % on a 4-cell chain) |

The neighbour's chain must carry n_B + n_C, which the solve shares out
between n_C and the neighbour's count. On counts this small the relative
change is large; it shrinks with the counts.

Limits of part 1: strip covers only (no 2 × 2 split of a face, no vertex
inside a coarse face); elliptic smoothing, the layout optimiser and the
edits that add, remove or split blocks or set chain counts refuse
T-junctions; the 3D trigger (part 2) and the tetrahedron insertion (part 3)
are next.

## The 3D trigger (step 2, part 2)

`block.optimize_layout` keeps the 2D path (MB1g optimiser plus the MB1h
trigger) for topologies extruded along z. A topology that is not, or that has
T-junctions, gets the 3D trigger (`block_insert::insert_prism_splits`) and no
layout optimisation, because the MB1g optimiser is a 2D model; the report says
so (`mode: 3d`, `layout_optimised: false`).

- **Measure.** `block::quality3d::measure` works on any block grid: for every
  inner face of the hexahedra the angle between its normal and the line
  through the two cell centroids and the ratio of the two centroid distances,
  for every corner the scaled Jacobian, whole and per block.
- **Candidates.** The flat block edges (`flat_edges`, 150° or more at both
  ends), flattest first, four per round, a column tried once.
- **Trial.** Split with each row end asked for (default both: `t_junction`
  and `sheet`), build the TFI grid, measure.
- **Keep rule.** The blocks the split touched (changed or new) get a worst
  angle lower by at least 2°; their worst scaled Jacobian is at least 0.4
  (the threshold HexIt uses, from the 3D blocking survey) and no worse than
  before; the whole grid's worst angle rises by at most 0.05°, no more cells
  pass 70°, nothing folds, the growth stays within max(2, start) × 1.05, and
  the split adds at most 10 % of the grid's cells. Of the trials that pass,
  the larger gain wins, then the fewer cells.

| Case | Before | After | Kept |
|---|---|---|---|
| Two-layer wedge turned so the flat edge runs along x (not extruded along z), cap lifted | 48.4°, SJ 0.16, 72 cells | 30.0°, SJ 0.52, 126 cells | `t_junction` (sheet: 150 cells) |
| Same, default cap 10 % | 48.4° | unchanged | none: "adds 54 cells, more than 10 % of 72" for `t_junction`, 78 for `sheet` |

The cap is the plan's default. On a fixture this small any split passes it;
on a real topology a split adds a few per cent.

## The tet split (step 2, part 3)

A hexahedral block shaped like a tetrahedron is split into the four
hexahedra of the tetrahedron (`block_tet::split_tet`, tool `block.split_tet`).

- **Detection.** Over the 24 right-handed frames of the block, the frame
  whose new slots 2, 5 and 7 are flat across axes 2, 1 and 0 (the flat
  corners of the triangles k = 0, j = 0 and i = 0 at P0 = v0) and whose v6
  has a scaled Jacobian of at most 0.5 in magnitude (the three faces there
  nearly coplanar: they are the fourth triangle's three quadrilaterals). The
  angles are read once in the block's own frame and mapped through each
  rotation. 150° for the detection, 120° by hand.
- **Counts.** The block's counts (a, b, c) along P0P1, P0P2, P0P3 are the
  corner counts n1, n2, n3 of the tri-center map of the tetrahedron, which
  the face P1P2P3 already carries; n0 is new, the mean of a, b and c by
  default, and the edges at P0 get n0 more cells, the new vertex n0 nodes
  from the far end (so for equal counts it lies at the middle).
- **Centres.** A triangle's centre at the mean of its three edge points (the
  two new vertices and the flat one), not of its corners: the mean of three
  points is inside the triangle they span, so the three quadrilaterals are
  convex for any n0. With the corner centroid an explicit n0 = 3 on counts
  (4, 5, 6) left a corner at 180° and a scaled Jacobian of 0.06; with the
  edge-point mean it is 0.35, and 0.49 for the default n0 = 5. The triangle
  split (`split_triangle`) now places its centre the same way; on the 2D
  wedge the kept split is unchanged and the new blocks no longer have a flat
  corner left to try.
- **Neighbours (maintainer card, "Conformal, prism columns").** A T-junction
  face cannot hold the three quadrilaterals of a triangle (its centre has
  valence 3 inside the face), so a block on a triangle at P0 is split with
  the tet: a triangle split along its edge away from the tet, with every
  block stacked on it. The columns of all three triangles and the tet's own
  edges at P0 share one row of n0 cells (`split_prisms`: several columns, one
  row, blocks to leave out, extra row edges), which runs on through the
  blocks on those edges as a sheet or stops there with T-junctions and the
  count solve. The columns reuse the tet's new edge vertices; the tet reuses
  the columns' centres at its faces.

| Case | Before | After |
|---|---|---|
| Tet (counts 4, 5, 6), boundary faces, default n0 = 5 | SJ about 0, flat corners | 4 blocks, SJ 0.49, worst angle 50.8°, no folds |
| Tet on a slab (the slab under the triangle P0P1P2, counts 4, 5, 6), trigger with the cap lifted | 81.6°, SJ 0.04, 33 cells above 70°, 240 cells | 50.8°, SJ 0.49, none above 70°, 880 cells; slab split into three |
| Same, default cap 10 % | | left: adds more than 10 % |
| Box on the edge P0P1 only | | `sheet` cuts it in two; `t_junction` keeps it whole with one junction and a count solve |

The fixture's tet is a single block of a tiny grid, so every split adds many
cells; the cap is meant for real topologies.
