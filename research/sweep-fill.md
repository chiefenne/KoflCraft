# Sweep fill: swept O-grids along the flow potential

Date: 2026-11-27. Requirements: BS-18 (automatic block decomposition),
BS-03 (templates), BS-17 and MD-25 (blocks next to the polyhedral mesher).
The join of the region behind the cuts as polyhedra (step 2) is in its own
section below.
Tool: `regions.fill`. Follows `regions.propose`
([field-region-proposal.md](field-region-proposal.md)). Decided by the
maintainer on a card (2026-11-27): the stations of a sweep are level sets of
the flow potential, and each section is meshed with the template's section.

## Question

How is a sweep region of the proposal filled with hexahedra, so that the
cells follow the duct through bends, end in calm planar sections where
another region begins, and keep the wall layers the polyhedral mesher will
continue?

## Method

All fields come from `regions.propose`: the tet scaffold, the potential p
(0 at inlets, 1 at outlets, no flux at walls) and the region labels.

1. **Stations on level sets of p.** A harmonic potential has level sets
   that never cross, meet the wall at right angles and, in a torus bend,
   are exactly the meridian planes of the elbow template. This is the
   harmonic volume parameterisation of generalised cylinders (Martin,
   Cohen & Kirby 2008, *Volumetric parameterization and trivariate
   B-spline fitting using harmonic functions*, SPM; Li, Guo & Qin 2007,
   *Harmonic volumetric mapping for solid modeling applications*, SPM),
   used here like the skeleton-driven sweep of tubular shapes (Livesu,
   Muntoni, Puppo & Scateni 2016, *Skeleton-driven adaptive hexahedral
   meshing of tubular shapes*, CGF 35(7)), with level sets in place of
   planes normal to the skeleton. Each station's wall loop is traced by
   marching triangles on the scaffold wall and projected onto the STL.
2. **Calm cuts.** A sweep ends at its port or at a cut before a box,
   junction or remainder. 96 level sets are sampled; a cut is the first
   one from that end whose wall loop is flat (largest distance from its
   best plane at most 0.001 of the diameter) and square to the centreline
   (tilt at most 0.5°). If none is calm, the calmest of that half is used
   and the report says so. Every node of a cut or port section is put into
   its plane, wall nodes by alternating projection onto plane and STL, so
   the cut is an exact planar interface for the join.
3. **Station spacing** at equal centreline length (through the loop
   centroids), with the axial size from the first section's wall spacing.
4. **Mitre corners.** Sharp wall edges (normals more than 30° apart)
   running more than 60° from the flow direction are grouped into chains;
   the station nearest to a chain's mean level is pinned to it, the loop is
   snapped onto the edges within 0.75 h, and when the snapped loop is flat
   (0.02) the station is the plane through the edges: the mitre plane.
5. **Section template**: the pipe O-grid ("butterfly"), an H-core of
   n x n cells (n = cells around / 4) and four ring blocks. With four sharp
   section corners (a rectangular duct) the ring corners sit on them;
   otherwise at 45° + 90° m of a rotation-minimising frame carried by
   double reflection (Wang, Jüttler, Zheng & Liu 2008, *Computation of
   rotation minimizing frames*, ACM TOG 27(1)). The first rows are wall
   offsets at the layer heights (first height, growth; mitred at sharp
   corners), the remaining rows grow geometrically to the core; the
   automatic row count keeps that growth at most 1.25.
6. **Interior nodes onto the level set.** On a curved station each
   interior node moves along the station normal to where p equals the
   station level: an outward search for a sign change, then regula falsi
   (Illinois), at most 0.45 of the way to the same node of the neighbouring
   stations. Nodes without a crossing in reach stay where the template put
   them and are counted.
7. **Hexahedra** join neighbouring sections node by node; the export goes
   through the validated FoehnCraft hand-off (closed cells, S·d > 0 on
   every interior face, positive volumes), with patches for walls, ports
   and each cut (`cut_<region>_start|end`, role interface).

## Results

Release build. Default settings: 32 cells around, 5 layer rows, growth 1.2,
first height 0.5 % of the first section's diameter, core at half the way to
the ring corners.

| Case | Hexahedra | Stations | Section | Min scaled Jacobian | Max non-orthogonality | Volume vs region | Ends | Time (scaffold) |
|---|---|---|---|---|---|---|---|---|
| Elbow (square duct, mitre bend) | 17 472 | 40, one the mitre plane | 4 corners, 12 rows, core 8 x 8 | 0.26 | 66.1° | -0.003 % | port, port | 2.7 s (1.9 s) |
| Manifold, inlet pipe with bend | 41 984 | 83 | round, 14 rows, core 8 x 8 | 0.70 | 23.5° (both) | -0.34 % | port, calm cut (flatness 2e-4) | 12.4 s (10.4 s) for both |
| Manifold, outlet pipe | 27 648 | 55 | round | 0.71 | | -0.31 % | calm cut, port | |

No inverted face in either mesh; the Manifold box (75.4 % of the volume)
is left for the polyhedral mesher. The round-pipe volume deficit is the
faceting of 32 cells around (a 32-gon has 0.64 % less area than its
circle; the scaffold is faceted too).

What was tried and dropped:

- **Newton steps along grad p** to put interior nodes on the level set:
  at the outer corner of the mitre grad p vanishes (stagnation), the steps
  jump, and the elbow had cells of negative volume. Replaced by the
  bounded search along the station normal.
- **Station spacing by the largest wall spacing** (Hausdorff distance
  between neighbouring loops) to keep the outer-wall cells of a bend short:
  around a mitre it crowds the stations at the inner corner, where the level
  sets bunch, into slivers (one inverted face on the elbow). Equal
  centreline spacing is kept.
- **A looser calm test** (flatness 0.01 only) put the Manifold cuts 3°
  askew inside the bend's exit; the tilt limit fixed it.

Benchmark against the classical sweep (ALG-16), planes normal to the
centreline: in a torus bend both give the meridian planes, which the
measured station tilt confirms. In the mitre the classical planes of the
two legs meet only at the mitre plane; the classical structured mitre
(sections blended from a square section to the mitre plane) gives
trapezoidal cells, while the level sets fan around the inner corner
(non-orthogonality 66° there against 24° in the Manifold bend). Not yet
measured numerically: it is the first follow-up below.

## Step 2: the region behind the cuts as polyhedra, joined at the cuts

Decided by the maintainer on two cards (2026-11-27): next is the join of the
box to the pipes, and at a cut the pipe adapts to the polyhedral mesher (the
mesher stays as it is; the last section of the pipe takes the positions the
mesher gives the cap).

### Method

1. **Cap of a cut (`CutCap`).** The polyhedral mesher (median dual of the tet
   scaffold with prism layers, [m4](m4-voronoi-core.md)) shrinks a planar
   port cap by one similarity about the mean c of its rim, x -> c + (x - c) s
   with s = 1 - H / r (H the layer stack, r the smallest distance of c from a
   rim edge), and puts the layer nodes of a rim point w at height h at
   c + (w - c)(1 - h / r). The cap given to the mesher is the preimage of the
   cut section under that map: the rim is the wall nodes of the cut, the
   inner points are the section's nodes from the layer top inwards, scaled
   by 1 / s. The layer rows of the cut section are moved to the mesher's
   layer nodes. Then both sides meet in the same points without changing
   the mesher.
2. **Quads split by the Delaunay diagonal** (Lawson's angle test). A split
   about the quad centre gives the thin rows next to the layer top flat
   triangles (an apex of about 164°) whose circumcircle holds the next row:
   no tetrahedralization has such a triangle as a face, so the boundary
   recovery of the scaffold (Shewchuk 2002) split them without end (the
   Manifold reached 750,000 splits and 6 GB in pass 15). With the diagonal
   split the recovery needs no split at all.
3. **The surface behind the cuts.** At each cut the STL facets that cross
   the cut plane, joined over the edges that cross it, are the curves where
   the plane meets the STL. The band is the one curve closest to the cut's
   axis that stays within the rim's reach plus one facet. The band and
   everything the band separates on the sweep's side (a flood over facet
   edges) are dropped. Taking every crossing facet within reach instead also
   took the Manifold's other pipe, which crosses the outlet cut's plane 72
   units away, and flooded the whole STL. The hole is zipped to the cap's rim
   by a strip of triangles, greedily by the shorter diagonal (Fuchs, Kedem &
   Uselton 1977, *Optimal surface reconstruction from planar contours*, CACM
   20(10)). The result is closed and outward: the kept STL, the strips (wall)
   and the caps (port patches the remesh keeps).
4. **Polyhedral mesh** of that surface with the standard pipeline: wall
   remesh at the cut's wall spacing (tolerance 0.01, minimum 1/8 of it), a
   Delaunay scaffold, the median dual, and prism layers with the sweeps' layer
   heights. When the sweeps' first heights differ (they follow each
   section's diameter), all sweeps are re-run with the smallest one, so one
   stack fits every cut.
5. **Wall offset at a box corner.** The mesher offsets the scaffold by the
   layer stack and lets the layers yield where a tetrahedron would turn over.
   Only the tetrahedron's own wall vertices yielded. On the Manifold box
   (H = 1.86, five times the stack of the polyhedral acceptance run) one
   tetrahedron at a corner of three walls has three vertices on two walls
   and its fourth inside, near the third wall. That node moves with the
   third wall, so the tetrahedron stayed turned over in every round. Now a
   tetrahedron still turned over a round later also lets the wall vertices
   next to its interior nodes yield. The case was fixed in three rounds.
   Meshes where every bad tetrahedron is fixed in one round do not change.
6. **The join.** The cap points of the polyhedral mesh are welded to the
   cut's nodes (grid hash, 1e-8 of the diameter). Each polyhedral cap face (a
   median-dual piece of a cap triangle, or the end of a layer prism) lies in
   one hex face of the cut, found by its centroid (point in polygon in the
   cut plane). The hex face is replaced by these pieces, reversed, and checked
   to have the same area. The points the pieces put on hex-face edges are
   spliced into every other face along that edge (swept faces and wall faces),
   so all cells stay closed. The FoehnCraft hand-off then merges faces two
   cells share more than once and validates the whole mesh.

### Results (Manifold, defaults)

| | Value |
|---|---|
| Cells | 132,520: 69,632 hexahedra (pipes), 62,888 polyhedra (box) |
| Validity | 0 inverted faces, closure 1.9e-14, smallest cell 0.026 |
| Volume | -0.17 % of the STL (pipe faceting); polyhedra -0.22 % of the region behind the cuts |
| At the cuts | 9,088 cap faces in the hex faces, 1,760 points spliced into edges, weld 6e-14, area match 3e-12 |
| Non-orthogonality (cell centroids) | at the cuts 48.2° (none over 70°), hexahedra 23.5°, polyhedra 87.7° (393 faces over 70°, the box mesh itself) |
| FoehnCraft read-back | valid, 0 negative volumes, 0 inverted faces (by hand, see below) |
| Time | 30 s (field scaffold 9.6 s, polyhedral box 15.2 s) |

The elbow has no cut and is unchanged. `join: false` gives the step 1 mesh
(the cuts as interface patches).

Pictures: `koflcraft-images/sweep/manifold-join-1-plan.png` (plan view
cut through both pipe axes), `-2-cut.png` (close-up at the outlet cut) and
`-3-interface.png` (the cut seen from the box: each O-grid face split into
the cap faces it holds).

What was tried and dropped:

- **Cap quads split about their centres:** the boundary recovery never
  converged (item 2).
- **The band as every crossing facet within reach:** it took the other pipe
  (item 3).

## Blended stations at mitre corners

Date: 2026-11-27. The level sets of the potential crowd at the inner corner
of a mitre (the potential is singular at a reentrant edge), so the stations
fan around it: on the elbow the cells next to the inner corner are short at
the inner wall and long at the outer wall, with 66.1° non-orthogonality.

### Method

The classical structured mitre (as in the swept and extruded block meshes
of mitred ducts): on each side of a mitre station m, the stations from a
start station a up to m become planes. Their normals turn from a's normal
to the mitre plane's normal by spherical linear interpolation (Shoemake
1985, *Animating rotation with quaternion curves*, SIGGRAPH), their centres
run on the chord from a's centroid to m's, both in equal steps. Each plane's
wall loop is traced on the scaffold wall (the same marching triangles as for
the level sets, on the signed distance to the plane), the loop nearest the
centre is kept and projected onto the STL in the plane.

- **Run length.** Over a run of z stations the planes turn by θ/z each. A
  wall point at distance r across the turn axis then moves by
  spacing ± r θ / z between stations; the inner wall keeps at least half
  the outer wall's spacing when z ≥ 3 r θ / spacing. On the elbow (r = 0.5,
  θ = 45°, spacing 0.125) that is 10 stations. The run stops at the sweep's
  ends and halfway to the next mitre's run.
- **Fallback.** When a plane meets the wall in no closed loop near its
  centre, or two neighbouring stations would come closer than 0.1 of the
  spacing anywhere on the wall, that side keeps the level sets and the
  report counts it (`unblended_mitre_sides`).
- `blended_stations` in the report counts the planes.

### Results (elbow, defaults)

| | Level sets | Blended planes |
|---|---|---|
| Max non-orthogonality (KoflCraft) | 66.1° | 59.1° |
| Max non-orthogonality (FoehnCraft) | 81.9° | 59.1° |
| Min scaled Jacobian | 0.26 | 0.37 |
| Volume vs region | -0.003 % | -0.002 % |
| Blended stations | 0 | 17 (both sides) |

Run length (factor k in z ≥ k r θ / spacing): k = 0.5 gives runs too short
to blend (66.1°), k = 1.5 blends one side and the other's planes would
cross (66.1°), k = 3 gives 59.1°, k = 6 and more (the whole leg) 59.6°.
k = 3 is kept. The Manifold has no mitre and does not change.

The 59° that is left sits where the square O-grid section's diagonal ring
edges meet the mitre: the section template alone gives 45.0° in the
straight legs (the corner cells of the butterfly are rhombi), and the
mitre's tilt adds to it. Picture:
`koflcraft-images/sweep/elbow-mitre-blend.png` (the floor wall, before and
after; `render_floor.py` beside it).

## H-grid for rectangular sections

Decided by the maintainer on a card (2026-11-28): rectangular sections
become a tensor grid clustered to all four walls; round sections keep the
O-grid.

### Method

The section keeps its topology (a ring of four blocks around an n × n
core), so the sweep, the join and the cut cap work unchanged:

1. **Ring.** The ring holds only the wall layers and one row more. Its
   columns leave the wall nodes along the wall normal (along the mitred
   bisector at a corner, so the corner node sits at the layer height from
   both walls).
2. **Wall nodes.** Each side has the corner, then nodes at top + q, where
   top is the ring's height and q is the core's geometric law. q starts with
   the next step of the layers' growth and grows to the middle of the side,
   at most 1.25 per cell. The law is mirrored, so opposite sides match. The
   number of cells per side is the fewest that keeps that growth on the
   longest side of any station (at least an eighth of `cells_around`).
3. **Core.** The core's boundary is the ring's top row; inside, it is the
   transfinite interpolation of its four sides. In a rectangle that is
   exactly the tensor grid, so every core line runs straight from a wall
   node to the opposite one.
4. **Cut.** The polyhedral mesher shrinks a cap by a similarity about the
   mean of its rim. Scaled out again, the straight core leans against
   the rim by up to the layer stack at the corners. With the wall nodes
   left in place, the cap was far from Delaunay there (48 edges per cap,
   up to 83° over) and the boundary recovery never converged. At a cut of
   an H-grid the wall nodes now slide along their sides to where the
   scaled-out core lines meet the wall. The cap is then the tensor grid
   and its rim, which is Delaunay. The layer rows follow the mesher's
   layer nodes as before. The stations before the cut take a share of the
   slide that falls off linearly, so the wall faces shear by at most 10°.

### Results (elbow, defaults)

| | O-grid | H-grid |
|---|---|---|
| Cells | 17,472 | 39,468 |
| Cells around a section | 32 | 88 (22 per side) |
| Max non-orthogonality, inner faces | 59.1° | 49.5° (at the mitre) |
| Max in the straight legs | 45.0° | 26.6° |
| Min scaled Jacobian | 0.37 | 0.49 |
| Volume vs region | -0.002 % | -0.002 % |
| FoehnCraft (with boundary faces) | 59.1° | 59.1°, valid, 0 negative volumes |

The 26.6° in the straight legs (atan ½) sits only in the corner cells of
the wall layers, where the layer rows of two walls meet on the diagonal.
The other faces in the legs are at about 0° (384 of 22,924 faces over 5°).
FoehnCraft's maximum is a boundary face: the wall cells at the mitre's
outer corner, measured from the cell centroid to the face centroid, are
at 59.1° in both grids. The cell count doubles because the core grows
from the walls at most 1.25 per cell in both directions.

Square T and fork (`ducts::tee`, `ducts::fork`: square stubs on a
junction cube): the cut caps are now Delaunay and the scaffold builds.
With two wall layers the T joined, with three to five the polyhedral
mesher's layer front collided next to the cuts. The cap's similarity moves
a rim vertex x by (c - x) / r per unit height. That is one along the wall
normal, and for a rim that is not a circle about c also up to one along
the wall (at a square's corners). The wall vertex next to the rim moved
along the normal only, so the layer prisms between them sheared by that
sideways shift over one wall edge. Thin scaffold tetrahedra there turned
over, the layers yielded and a cap quad folded.

**Fading the cap's shift into the walls** (maintainer's card,
2026-11-28). Each wall vertex reached from the rim over wall edges within
four layer stacks of the cap plane, on walls across the plane, now takes
the tangential part of the similarity at its place. That part is scaled
by a smoothstep from one at the plane to zero four stacks away. The part
along the wall normals is taken out, so the layer heights do not change.
For a round cap about its centre the tangential part is zero, so round
pipes do not change.

| Default 5 layers | Turned-over tetrahedra (round 0) | Result |
|---|---|---|
| T, no fade | 28 | layer front collision |
| T, fade | 11, yielded in 4 rounds | valid, 80,389 cells, 41.4° at the cuts |
| Fork, fade | 0 | valid, 106,713 cells |

## More than two openings behind the cuts

The polyhedral mesher had one cap similarity per port patch id, and it had
two ids (inlet and outlet). With three or more openings the fill refused
the join. Now:

1. **One similarity per cap.** `front::advance` splits each cap id into
   its connected caps (union-find over the shared edges of its facets) and
   gives each cap its own centre, normal and in-plane similarity
   (x -> c + (x - c)(1 - H/r)) with its own radius. Two caps of one id that
   touch at a vertex are refused. One similarity for the whole id would
   shrink every cap towards the mean centre of all of them.
2. **Shared ids.** The cuts take the free ids first; past two openings the
   rest share the first free id.
3. **The join tells the cuts apart by their planes.** A boundary face of
   the polyhedral mesh belongs to a cut when its id matches and every node
   lies in the cut's plane (100 times the weld tolerance); among coplanar
   cuts it goes to the nearest centre.
4. **No ring of STL nodes next to the cut plane.** A cut is a level set,
   so it can pass within a few thousandths of a ring of STL nodes, or cross
   it when it is slightly tilted. The band of plane-crossing facets then
   zigzags through that ring, and the strip to the rim holds needles (179°
   on the three-pipe box). The scaffold's boundary recovery split those
   without end and ran out of memory. Now STL nodes within a quarter of the
   rim spacing in front of the plane count as the sweep's side, so the
   strip reaches back to the ring before.

### Results (three-pipe box, defaults)

The fixture is `ducts::pipe_box`: the box [-1, 1]^3 with round pipes
(radius 0.5, length 4, 64 sides) on its -x face (inlet), +x and +y faces
(outlets). Every pipe is swept from its port to a calm cut; behind the
cuts are the box and the junction, with three openings.

| | Value |
|---|---|
| Cells | 105,638: 53,760 hexahedra (pipes), 51,878 polyhedra (box) |
| Validity | 0 inverted faces, closure 6.8e-15, min scaled Jacobian of the hexahedra 0.71 |
| At the cuts | 13,632 cap faces joined, weld 4.5e-16, area match 5.8e-14, non-orthogonality 47.3° (none over 70°) |
| Hexahedra | non-orthogonality 23.7° |
| Polyhedra | 87.8° (786 faces over 70°, the box mesh itself, as on the Manifold) |
| Not meshed | 0.19 % (pipe faceting) |
| FoehnCraft read-back | valid, 0 negative volumes, volume equal to the boundary volume; at most 87.7° on its measure (the joined Manifold: 88.8°) |
| Time | 27 s |

Square sections (the T and the fork of `ducts.rs`) failed: their O-grid
cut caps were not Delaunay (about 110 edges per cap, up to 95° over, near
the square's corners, where the O-grid's first rows sit off the rim
points), so the scaffold's boundary recovery did not converge. A 2 %
jitter of the inner nodes did not help. The H-grid (above) fixes the caps.

## Tubes cut at both ends, several shells behind the cuts

Once the proposal labels parallel tubes as sweeps (research
`field-region-proposal.md`, 2026-11-28), a tube between two headers is
swept from cut to cut, and the region behind the cuts is two closed
shells, one per header. Two changes:

1. **Bands first** (`remainder.rs`): the ring of plane-crossing facets of
   every cut is marked before any sweep side is flooded. Before, the band
   of each cut was found and flooded in turn: the flood from a tube's
   first cut ran through the tube, past its other cut (not marked yet) and
   on through the whole STL, which left 0 holes ("the cuts leave 0 holes
   in the STL, expected 14").
2. **One mesh per shell** (`fill.rs`): the closed surface behind the cuts
   is split into its connected shells; each is meshed by the median dual
   with the caps it holds (its own cap ids and rim spacing) and joined in
   turn. The mesher takes one connected shell.

| Case | Result |
|---|---|
| Three tubes between two headers (`ducts::header_tubes`, square) | 3 tubes swept cut to cut, 2 shells joined, 0 inverted faces, the cuts interior |
| Heat exchanger tube side (Gmsh fluid from the open STEP, 82,044 facets) | 7 tubes swept (679,936 hexahedra), 2 bonnets (691,642 polyhedra), 1,371,578 cells, 0 inverted faces, 50.1° at the cuts, 89.3° inside the bonnets, 230 s |
| Five round tubes between two box headers (Gmsh) | the tubes are swept; the polyhedral layers fail at one vertex of a header ("no common inward wall direction") |

A square tube flush with a header wall (the first `header_tubes` layout,
its outer tubes flush with the header's end walls) collides in the
polyhedral layer front at the cut corner on the flush side; the square
offset tee with its branch flush with the floor does not finish its
polyhedral build. Both are open (the same gap as the square tee with a
round branch).

## Limits and next steps

- Behind the cuts the join takes any number of openings, round or
  square. A rectangle much longer than wide fails at the cut: the
  mesher's single similarity per cap scales by the short side, so the
  scaled core leaves the rim on the long one. Symmetry planes behind the
  cuts are refused.
- The acceptance runner reads the elbow and the joined Manifold back with
  FoehnCraft (#87). The Manifold is valid: 132,520 cells (133,056 since
  the proposal's box test by section area), at most 88.8° on
  FoehnCraft's non-orthogonality measure (the plain polyhedral Manifold:
  90.0°).
- The polyhedral box keeps the mesher's own quality: 87.7° non-orthogonality
  inside the box, none at the cuts.
- At a mitre the stations next to the mitre plane are blended planes, and
  a rectangular section is an H-grid (see above). The elbow's inner faces
  reach 49.5° at the mitre; its wall cells at the outer corner reach 59.1°.
  An H-grid has twice the cells of the O-grid.
- Section rows of the ring blocks bow slightly inwards (each column is
  graded on its own); a 2D smoothing of the section is open.
- A sweep needs one closed wall loop per level set and the same section
  type throughout (four sharp corners or none). Junctions, obstacle
  passages and boxes have no template yet; 3D tetrahedron and prism blocks
  are planned in the thread on 3D tetrahedron insertion.

## Provenance

Clean room: built from the papers cited above (and Lawson's angle test and
Shewchuk 2002, *Constrained Delaunay tetrahedralizations and provably good
boundary recovery*, IMR, for the cap split) and documented formats; no
source code of any meshing or geometry project was read. Fixtures: the
elbow and the Manifold from `test-fixtures/` (provenance beside them); the
confidential cooling jacket was not used.
