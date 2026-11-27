# Sweep fill: swept O-grids along the flow potential

Date: 2026-11-27. Requirements: BS-18 (automatic block decomposition),
BS-03 (templates), BS-17 and MD-25 (blocks next to the polyhedral mesher).
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

## Limits and next steps

- Only the sweeps are meshed. Next: the remainder (the Manifold box)
  through the polyhedral mesher with the conformal join at the cuts
  (MD-25): the remainder takes the cut's quads as its port mesh, the hex
  faces split into the median-dual pieces, and the near-wall rows of the
  sweep are the layer heights, so the poly layers end on the same points.
- At a mitre the stations fan around the inner corner; blending the
  stations between a calm section and the mitre plane (the classical mitre
  mesh) would give better cells there.
- Section rows of the ring blocks bow slightly inwards (each column is
  graded on its own); a 2D smoothing of the section is open.
- A sweep needs one closed wall loop per level set and the same section
  type throughout (four sharp corners or none). Junctions, obstacle
  passages and boxes have no template yet; 3D tetrahedron and prism blocks
  are planned in the thread on 3D tetrahedron insertion.

## Provenance

Clean room: built from the papers cited above and documented formats; no
source code of any meshing or geometry project was read. Fixtures: the
elbow and the Manifold from `test-fixtures/` (provenance beside them); the
confidential cooling jacket was not used.
