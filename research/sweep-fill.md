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

## Limits and next steps

- The join needs at most two openings in the region behind the cuts (the
  polyhedral mesher's two port patch ids). STL ports behind the cuts take
  them first, and symmetry planes behind the cuts are refused.
- `regions.fill` is not in the acceptance runner yet, so CI does not read a
  joined mesh back with FoehnCraft. By hand (`foehn mesh` on the Manifold
  VTU, FoehnCraft main of 2026-11-27) it is valid: 132,520 cells (68,664
  hexahedra, 63,856 polyhedra, among them the hexahedra with split cut
  faces), 0 negative volumes, 0 inverted faces, closure 3.4e-16. On
  FoehnCraft's non-orthogonality measure it is at most 88.8°, with 97,304
  faces over 70°. The plain polyhedral Manifold of the acceptance run
  measures 90.0° and 271,252 faces on the same measure.
- The polyhedral box keeps the mesher's own quality: 87.7° non-orthogonality
  inside the box, none at the cuts.
- At a mitre the stations next to the mitre plane are blended planes (see
  above); the square O-grid section keeps 45° at its diagonals in the
  straight legs, which adds to the mitre's tilt (59° on the elbow).
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
