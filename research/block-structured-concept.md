# Block-structured mesher variant: concept

**Status:** concept, 2026-11-24. Requested by the maintainer
("make a concept for a block-structured mesher as a mesher variant").
The open questions are put to the maintainer one decision card at a time;
each answer is recorded in section 9 and turned into a requirement of the
group BS (section 8). Nothing here is built yet.

## 1. Why a second mesher variant

The polyhedral variant (OBJ-20: clipped Voronoi core, poly-prism layers)
meshes any closed geometry automatically. A block-structured mesher gives
something it cannot: hexahedral cells aligned with the flow, strict
control of wall-normal spacing and orthogonality, low cell counts for a
given accuracy in ducts, pipes, channels, airfoils and blade passages, and
meshes whose topology an engineer can reason about.

## 2. Fit with the existing structure

The variant is a second mesher behind the same engine, not a second
program (OBJ-05, OBJ-06):

- **One Rust engine.** Block topology, edge distributions and in-block
  grid generation live in `mesh-algorithms` (module `block`); the block
  data model in `mesh-core`; orchestration and evidence in `mesh-engine`.
  No new crate unless the module outgrows it (as with boundary layers).
- **Every action a tool.** `block.*` tools in the registry (AG-12), the
  same ids from CLI and MCP, each step journalled and replayable (AG-13),
  each tool declaring its effects for the approval table (AG-15).
- **Headless first.** A block topology is a versioned JSON file
  (`koflcraft.block-topology`, `format`/`format_version`, SW-20), written
  and changed only through tools, like the patch map (BC-02) and the
  domain map (MD-15).
- **Same geometry and boundary inputs.** STL and STEP (section 4b), patch
  maps (BC-01..06), domain maps (MD-15), grading sources (GRD-02/03).
- **Same output path.** The result is a `mesh-core` volume mesh (hexahedra
  as general polyhedra with 6 quad faces), validated by the same checks
  (AGENTS.md: closed cells, two-sided interior faces, positive volumes,
  conservation) and written by the existing VTU / MSH 4.1 / polyMesh
  writers for FoehnCraft (FC-01..06). Block structure (i, j, k ranges per
  block) is kept in the evidence so a structured writer can follow.

## 3. Data model (proposal)

```text
BlockTopology (koflcraft.block-topology v1)
  vertices   : id, position, optional geometry snap (point, curve, surface)
  edges      : vertex pair, optional curve (line, arc, spline, CAD edge,
               polyline on STL), node distribution (count + law)
  faces      : 4 edges, optional surface (plane, CAD face, STL patch),
               boundary patch name or interface
  blocks     : 8 vertices (hex), local axes, cell counts (ni, nj, nk),
               domain name (MD-06), optional collapse (prism/O-grid core)
  connections: derived; block faces shared by two blocks, with the
               orientation map between their local (i, j) frames
```

Edge node distributions: uniform, geometric (first spacing + ratio),
two-sided tanh (Vinokur), and "first cell height from the layer settings"
(BL-01) so the structured wall spacing uses the same inputs as prism
layers. Edge counts are propagated through opposite edges of each block
(the structured consistency constraint) and checked before any grid is
built.

## 4. Pipeline

1. **Topology:** define or generate blocks (section 5).
2. **Consistency check:** every block a valid hexahedron, shared faces
   matched, edge counts consistent along every chain of opposite edges,
   every boundary face assigned to a patch, every block to a domain.
3. **Edge grids:** distribute nodes on each edge along its curve.
4. **Face grids:** transfinite interpolation (Coons patch) of each block
   face from its four edge grids, projected onto the face's surface;
   optional smoothing on the surface.
5. **Volume grids:** 3D transfinite interpolation per block, then
   smoothing (section 6).
6. **Assembly:** merge nodes of shared faces exactly (they are built
   once per shared face, so no tolerance welding), build `mesh-core`
   cells, faces and patches.
7. **Validation and quality:** the standard checks plus structured
   metrics: Jacobian (min and per-block), skewness, orthogonality at
   walls, aspect ratio, growth ratio along grid lines.
8. **Export** through the existing writers, evidence with per-block
   statistics.

## 5. Ways to obtain the block topology

From least to most automatic; card 2 decides which come first.

- **Explicit:** an agent or a person adds vertices, edges and blocks
  through `block.*` tools (as `blockMesh` dictionaries do, but as tool
  calls with checks after every step).
- **Templates:** parametric topologies for standard cases: box, channel,
  pipe O-grid (butterfly / H-core), elbow, T-junction, annulus, C-grid and
  O-grid around an airfoil, blade passage (H/O/C).
  A template is a tool that writes a topology the explicit tools can then
  edit.
- **Sweep and extrude:** a 2D block topology swept along a path or
  rotated about an axis (ducts, axisymmetric parts).
- **Automatic decomposition (research):** graphics literature, first
  place to look under OBJ-19:
  - Kobbelt's group and its alumni: *Integer-Grid Maps* and
    *Mixed-Integer Quadrangulation* (Bommes, Zimmer, Kobbelt 2009;
    Bommes, Campen et al. 2013) for quad layouts on surfaces; *HexEx*
    (Lyon, Bommes, Kobbelt 2016) for extracting hexes from a volumetric
    parametrisation; *The 3D Motorcycle Complex* (Brückler, Gupta, Campen
    2022) for a block decomposition from a seamless parametrisation, the
    closest match to "find the blocks automatically".
  - Lévy's group: hex-dominant meshing from a frame field (Sokolov, Ray,
    Untereiner, Lévy 2016); frame fields and polycube maps.
  - Crane's group: direction fields (Knöppel, Crane, Pinkall, Schröder
    2013), *Boundary First Flattening* (Sawhney, Crane 2017) and the
    rectangular surface parameterization already listed for BL work,
    for surface block faces.
  - Outside the key groups: polycube maps (Gregson, Sheffer, Zhang 2011;
    Livesu et al. 2013), the classical alternative for ALG-16.
  These need a later screening note (LIC-12..15 apply).

## 6. Grid generation inside a block

Classical baseline (ALG-16 reference): transfinite interpolation followed
by elliptic smoothing (Winslow / Thompson–Thames–Mastin with control
functions for wall clustering and orthogonality, Sorenson's method in
2D). Graphics candidates to benchmark against it: harmonic and
conformal maps (Crane's group), and the volumetric parametrisation behind
integer-grid maps, which gives smooth grids across block faces without
per-block seams. Validity first: no negative Jacobian anywhere, then
orthogonality and smoothness.

## 6a. Taken from PyAero v3.0

The maintainer asked (2026-11-24) to check PyAero's `v3.0` branch for
modern parts and take them into use. That branch (last commit 2026-11-07)
holds a structured grid framework for 2D airfoil meshes, designed in its
`docs/.../specs/2026-08-22-structured-grid-framework-design.md`. The ideas
below are adopted as design input; KoflCraft writes them fresh in Rust
for 3D (no line-by-line port).

- **Orthogonal choices.** Topology, outer-boundary shape, volume
  algorithm and post-smoother are independent, first-class choices. A
  topology builder produces *frames* (the boundary polylines of a block
  plus metadata such as periodic seams and wake cuts); a volume algorithm
  fills a frame without knowing which topology made it. In KoflCraft:
  templates produce block topologies, `block.grid` takes the algorithm,
  `block.smooth` the smoother, each a separate tool argument.
- **Wall-distribution invariant.** Grid lines leave the wall exactly at
  the prepared wall nodes, never resampled; a test checks node-for-node
  equality for every algorithm.
- **Three volume algorithms:** TFI (standard Boolean sum and a Hermite
  variant with boundary derivative terms for tangent control), elliptic
  (Winslow / Thompson–Thames–Mastin with Steger–Sorenson control functions
  for wall spacing and angle), and **hyperbolic marching** (Steger–Chan:
  cell-volume specification plus orthogonality, implicit smoothing against
  front crossing, finished by an algebraic hand-off onto the prescribed
  outer distribution). Hyperbolic marching is added to the KoflCraft
  baseline next to card 4's TFI + elliptic.
- **Perpendicular near-wall block.** N layers built along exact surface
  normals with first height h1 and growth g; the offset is checked against
  the local radius of curvature and the block is thinned there with a
  warning instead of folding; normals are blended toward the bisector at
  sharp corners. Its outer rim is the inner boundary the volume algorithm
  fills. Same inputs as BL-01.
- **Honest controls.** A setting an algorithm cannot honour (an
  orthogonal-angle request to standard TFI) is refused or disabled, never
  silently ignored (AG-08).
- **Corner and trailing-edge handling.** Sharp versus blunt trailing edges
  are detected, not switched; a blunt edge changes the topology (an extra
  wake-strip block) instead of distorting the main block.
- **Post-smoothers:** constrained Laplacian, elliptic and angle-based,
  with a quality monitor (skewness, orthogonality, aspect ratio) before
  and after.

## 7. Open questions (one card each)

1. Main purpose and scope.
2. How topologies come in first (explicit, templates, sweeps, automatic).
3. Block connection: point-matched only, or also non-matching faces.
4. In-block grid method and smoothing.
5. Output: unstructured hexes only, or also structured multi-block files
   (CGNS structured, Plot3D).
6. Coupling with the polyhedral variant (per-domain mesher choice,
   structured layers inside a polyhedral mesh).
7. Use of the maintainer's PyAero and pyaeroflow repositories.
8. Milestone and priority.
9. Where this concept lives (research note plus a requirements section).

## 8. Proposed requirements

Recorded as section 4c of [the requirements](../docs/requirements.md)
(BS-01 to BS-19), with OBJ-21 and milestone MB1 (own track after M5).

## 9. Decisions

| Card | Question | Answer | Date |
|---|---|---|---|
| 1 | Main purpose | Classic CFD topologies first: hex multi-block for ducts, pipes (O-grid), airfoils (C/O/H), channels; the user or an agent gives the block topology; automatic hex decomposition stays a later P item | 2026-11-24 |
| 2 | Topology input first | Templates plus explicit tools: parametric templates (channel, pipe O-grid, elbow, T-junction, airfoil C/O grid) write a topology that `block.*` tools can then edit | 2026-11-24 |
| 3 | Block connection | Point-matched only, first: shared block faces have identical nodes; counts propagate through edge chains; hanging nodes and non-matching faces are later L items | 2026-11-24 |
| 4 | In-block grid method | TFI plus elliptic smoothing (Winslow / Thompson–Thames–Mastin, wall orthogonality control) as baseline; harmonic and conformal maps of the key groups benchmarked against it per ALG-16 | 2026-11-24 |
| 5 | Output formats | Unstructured hexahedra through the existing VTU / MSH 4.1 / polyMesh writers; per-block i, j, k ranges kept in the evidence; structured CGNS and Plot3D later (P). FoehnCraft refuses structured CGNS (its D188) | 2026-11-24 |
| 6 | Coupling with the polyhedral variant | Standalone now: first a pure multi-block mesh of the whole domain; mesher choice per domain (MD) and structured zones joined to a polyhedral core are later L items | 2026-11-24 |
| 7 | PyAero and pyaeroflow | First "don't use them"; refined the same minute: check PyAero branch `v3.0` for modern parts and take them into use. Done in section 6a (design ideas, no code port); pyaeroflow and PyAero master stay unused | 2026-11-24 |
| 8 | Milestone | Own track after M5: new milestone MB1 (templates, explicit tools, TFI + elliptic, export) starts once M5 boundary layers settle; polyhedral work keeps priority | 2026-11-24 |
| 9 | Where it lives | This research note plus requirements v1.0: section 4c (BS-01 to BS-19), OBJ-21, milestone MB1 | 2026-11-24 |

## 10. Build steps

| Step | Scope | Date |
|---|---|---|
| Start | The maintainer asked to start MB1 at once ("can we start with the blockmesher?"), although card 8 had placed it after M5: MB1 runs as a parallel track next to M5 | 2026-11-24 |
| MB1a | Card: "Channel + pipe O-grid end to end". Topology model (`mesh-core::block`) and file `koflcraft.block-topology` 1.0 (`mesh-engine::block`); check of indices, handedness, curves, laws, count chains, faces and patches; uniform, geometric and two-sided tanh (Vinokur) laws with achieved spacings reported; lines and circular arcs; 2D TFI of faces and 3D TFI of blocks with parameters blended from the edge distributions; every edge and face gridded once and read through the orientation map of each block's frame (no welding); corner Jacobians, aspect ratios, conservation against the boundary-enclosed volume; VTU, MSH 4.1, FoehnCraft VTU and evidence with block ranges; templates `channel` and `pipe_ogrid`; tools `block.template`, `block.check`, `block.grid` | 2026-11-24 |

MB1a results (`block.grid`, release build): pipe R = 0.05, L = 0.5, 8 core cells per quarter, 10 radial, 40 axial, wall spacing 5e-4: 15,360 hexahedra, minimum scaled Jacobian 0.707 (the core corners), volume 0.641 % below πR²L, which is the chord error of 32 wall nodes per circle; the sum of cell volumes equals the boundary-enclosed volume to rounding. Channel 1 × 0.1 × 0.05 with 50 × 20 × 12 cells and wall spacing 5e-4: 12,000 hexahedra, exact volume, scaled Jacobian 1.

Open after MB1a, each to be put on a card: elliptic smoothing (Winslow / Thompson–Thames–Mastin with Steger–Sorenson control), the orthogonal wall block, the explicit `block.*` edit tools (BS-04), elbow and airfoil C-grid templates, polyMesh in the bundle, the FoehnCraft round trip in the acceptance runner, and the ALG-16 benchmark against harmonic maps.

### MB1b: elliptic smoothing (2026-11-24)

Cards: "Elliptic smoothing + wall orthogonality" as the next step, and "Smooth across interfaces too" for the faces between blocks.

Built (`mesh-algorithms::block::smooth`, `block.grid` with `algorithm: {"kind": "elliptic"}`): every node off the boundary moves. Nodes inside a block or a shared face solve the Thompson–Thames–Mastin system in the index space of their home block; across a shared face the stencil reads the neighbour block's node one step off the same face node, so any pair of local frames joins into one structured neighbourhood. Nodes on inner block edges (the O-grid core corners, valence 3) move to the average of their neighbours in the plane normal to the edge. Distribution controls are frozen from the TFI grid, so a stretched orthogonal grid is a fixed point. At walls the system evaluated with the target slope (inward normal times the TFI first spacing; Steger and Sorenson 1979, 3DGRAPE 1989) gives the tangential controls. The off-wall control is the TFI distribution control plus an integral feedback on log(h / h_target).

Two findings, measured:
- **The Sorenson system alone does not pin the first spacing.** Solving all three controls from the wall equation with a one-sided cubic for r_ζζ converged with the first spacing 20 % (face interior) to 40 % (block edges) too large. Every spacing has a matching off-wall control that is self-consistent. Fix: the off-wall control comes from spacing feedback, and only the tangential ones from the system. Result: spacing kept to 1e-6 on the face interior.
- **Control relaxation 0.2 oscillates on coarse grids** (3 radial cells: last change 0.95 of the shortest edge after 3,000 sweeps); 0.05 converges in 618. The default is now 0.05.

Pipe R = 0.05, L = 0.5, 8/10/40 cells, wall spacing 5e-4 (release build, 8.8 s):

| | TFI | Elliptic |
|---|---|---|
| Wall angle of first grid line, max / mean | 3.12° / 2.26° | 0.37° / 0.18° |
| First spacing change (wall face interior) | — | 5e-7 |
| Minimum scaled Jacobian | 0.707 | 0.695 |
| FoehnCraft non-orthogonality, max / mean | 23.4° / 2.12° | 19.4° / 1.34° |
| FoehnCraft skewness, max / mean | 0.328 / 0.025 | 0.275 / 0.018 |
| Volume, cells and faces | equal | equal |

The minimum scaled Jacobian is set by the O-grid core corners on the inlet and outlet. Those are boundary nodes and do not move.

### MB1c: elbow and airfoil C-grid templates (2026-11-24)

Cards: "Elbow + airfoil C-grid templates" as the next step, "Project faces onto analytic surfaces" for the elbow wall, and "NACA 4-digit + coordinate file" for the airfoil shape.

Built:
- **Elbow** (`templates::elbow`): the pipe O-grid swept through stations; the bend is split into segments of at most 30°, so every longitudinal block edge is an exact circle about the bend axis; optional straight legs. Axial cells are shared out by centreline length.
- **Face projection** (BS-07, the blockMesh `projectFaces` idea): the topology gains `surfaces` (cylinder, torus) with the boundary faces projected onto them. Face interior nodes are projected after the face TFI and before the block TFI, so the block interiors see the exact wall. The check refuses a surface on a shared face, a face projected twice, and a face whose corners are off its surface.
- **Spline edges**: a natural cubic spline in the chord-length parameter through given points, evaluated by arc length from a table of 64 samples per segment.
- **Handedness along tangents**: the right-handed check uses each curved edge's tangent at the corner instead of its chord, so a block may start at an airfoil's leading edge, where the chords of the surface edge and the normal edge are collinear.
- **Airfoil C-grid** (`templates::airfoil_cgrid`), extruded along +z: a front and a rear block on each side (split at the data point nearest x/c = 0.3, because the upper and lower surface between the same two vertices would be one edge) and a wake block above and below the wake cut; a half circle about the trailing edge and straight lines to the wake end. The surface edges cluster at the leading and trailing edges (two-sided tanh), the C's arcs carry the same relative distribution, the normal edges grow geometrically from the wall and the wake cut. Sections: NACA 4-digit with the closed trailing edge (160 cosine-spaced points per side, analytic area), or a Selig or Lednicer file, moved to the leading edge at the origin, turned and scaled to unit chord, with a blunt trailing edge closed by a linear thickness correction.

Results:

| Case | Cells | Wall | Volume against the reference |
|---|---|---|---|
| Elbow R = 0.5, bend 1.5, 90°, legs 1, TFI without projection | 1,920 / 15,360 / 122,880 | 5e-3 off the torus at every level | stalls with the wall error |
| Same with projection | same | 4e-16 | 2.7 %, 0.67 %, 0.17 % short (second order) |
| NACA 0012 C-grid, farfield and wake 5 chords | 864 / 3,456 / 13,824 | 2.6e-7 chords off the section | 0.33 %, 0.095 %, 0.024 % short |
| NACA 2412, same | same | 5.3e-7 | same (the C's chords dominate) |

FoehnCraft read-back by hand: elbow R = 0.05 (15,360 cells) valid, TFI non-orthogonality max / mean 23.7° / 2.17°; NACA 2412 with farfield 10 chords and wall spacing 5e-4 chords (11,520 cells) valid, non-orthogonality max / mean 60.4° / 3.0°, skewness max 0.27.

Two findings:
- **Elliptic smoothing makes the elbow's worst cell worse.** Wall angle max 3.17° → 1.84°, but the FoehnCraft max non-orthogonality rises from 23.7° to 38.7°. The worst cells are not located yet; candidates are the joins between the straight legs and the bend, where the frozen TFI controls jump.
- **An extruded section with one span cell is not smoothed at all**: every node lies on a symmetry plane, and the smoother keeps all boundary nodes fixed. The airfoil's worst TFI cells (60° near the leading edge) stay. Sliding nodes on planar symmetry patches (mirror ghost points) would fix it.

### MB1d: smoother base, symmetry sliding and the elbow (2026-11-24)

Cards: "Smoother: symmetry sliding + elbow" as the next step, then "Thomas–Middlecoff controls" to try first for curved ducts. After the measurement below, a third card asks for the smoother base; this step was built on its recommended option (full residual freeze plus wall controls).

Where the elbow's worst cells were: the O-grid core corners at the joins between bend segments. Interior block vertices were held fixed, and the Laplace part of the system straightens the swept lines, so the whole core drifted about 5 % of the radius toward the bend centre and sheared the cells next to the fixed vertices.

Built:
- **Interior block vertices move**, normal to the grid line through them, and nodes on inner block edges move in the plane normal to the edge, both with the TFI offset frozen.
- **Full residual freeze.** The residual of every node at the TFI grid is stored and subtracted at every update. The TFI grid becomes a fixed point, and only the wall controls move nodes.
- **Symmetry sliding.** A patch with role `symmetry` that is planar (to 1e-9 of the domain extent) lets its nodes slide in the plane. The stencil reaches across the plane to the mirror images of the interior nodes, the node is projected back onto the plane after each update, and the wall control continues to the rims of wall faces on such a plane. A symmetry patch that is not planar keeps its nodes fixed.
- **Line Gauss–Seidel.** Each grid line inside a block is solved at once with the Thomas algorithm, and the line axis turns every sweep. Convergence: the largest move below 1e-4 of the shortest TFI edge at the node, and the first spacing off every wall within 1e-4 in log. The default iteration cap went from 3,000 to 10,000.

FoehnCraft max non-orthogonality (release build; pipe R = 0.05 with 8/10/40 cells and wall spacing 5e-4, elbow R = 0.05 at the same resolution, 15,360 cells each):

| Variant | Pipe | Elbow |
|---|---|---|
| TFI | 23.4° | 23.7° |
| MB1b smoother | 19.3° | 34.5° to 38.7° |
| Thomas–Middlecoff (residual kept at fixed nodes, Laplace-interpolated inside) | 23.1° | 27.3° to 33° |
| Full residual freeze + walls (built) | 23.3° (mean 1.99°) | 23.45° (mean 2.07°) |

| Case (built smoother) | Sweeps | Time | Wall angle max / mean, TFI → smoothed | FoehnCraft non-orthogonality max / mean | Skewness max |
|---|---|---|---|---|---|
| Pipe | 129 | 1.0 s | 3.12° / 2.26° → 0.44° / 0.29° | 23.3° / 1.99° | 0.328 |
| Elbow | 196 | 1.8 s | 3.17° / 2.25° → 0.63° / 0.37° | 23.45° / 2.07° | 0.327 |
| NACA 2412 C-grid, farfield 10, 80/48/40 cells, wall spacing 5e-4, one span cell | 2,478 | 35 s | 60.4° / 14.8° → 34.5° / 4.0° | 55.9° / 2.73° (TFI 60.4° / 2.99°) | 0.155 (TFI 0.267) |

No inverted faces in any case. The elbow's minimum scaled Jacobian is 0.704 against 0.706 for TFI.

What the freeze costs: interior skew that TFI puts away from the walls stays. On the pipe the MB1b smoother reached 19.3° because it was free to move the core; the frozen one keeps 23.3°. On the airfoil the worst cell sits in the lower front block below the leading edge, in the frozen TFI interior. The airfoil needs about 2,500 sweeps because the wall controls near the trailing edge and the wake cut settle slowly. Next candidates, each for a card: a 2D cap-smoothing pass for extruded sections, a line solver with relaxed freeze away from walls, and an O-grid around the airfoil.

Relaxed freeze, measured after MB1d and not kept: the frozen residual weighted by exp(−a (d − 1)), d the hop distance from the nearest wall node. FoehnCraft max non-orthogonality for a = 0 / 0.1 / 0.25 / 1: pipe 23.3° / 20.7° / 19.3° / 19.4°, elbow 23.45° / 26.0° / 31.0° / 34.8°; the NACA 2412 C-grid did not converge for any a > 0 (last move 5 % of the shortest edge after 10,000 sweeps). The elbow core drifts toward the bend centre again as soon as the interior is free.

### MB1e: edit tools (2026-11-24)

Cards: "block.* edit tools (BS-04)" as the next step; scope "Both in one step" (editing existing blocks was built first while the card was open, then the topology tools).

Built (`mesh-engine::block_edit`, tools in `mesh-tools::block`):
- `block.list`: vertices with points, blocks with vertices and counts, edges with curve and law, patches with role and face count, surfaces, the reference and the check with its chains.
- `block.set_count`: the count of a whole chain of opposite edges, found from the named edge by growing over shared block edges, so a count edit can never break BS-05.
- `block.set_distribution`, `block.set_edge_curve`: one edge, named by its two vertices in either order. A law or spline given from the other end is stored reversed, so it means the same spacing at the same vertex.
- `block.move_vertex` (`to` or `by`) and `block.set_patch` (rename, role).
- `block.add_vertex`, `block.remove_vertex` (only a vertex nothing uses; higher indices move down).
- `block.add_block`: faces that meet an existing block's face become shared and leave their patches (a patch left empty is removed); the other faces go to a named patch.
- `block.remove_block`: faces it shared open into a named patch; its own boundary faces leave their patches and surfaces; unused edge entries are dropped.
- `block.split_block`: split after cell `at` along an axis. The split runs through the whole chain, oriented edge by edge, so faces stay point-matched. New vertices sit on the TFI grid nodes. Each half of a split edge keeps a line as a line, gets an arc through one of its nodes, or keeps the spline's own points on its side; its law is fitted to the node spacings measured along the new curve (tanh with Vinokur's slopes corrected until the end cells match, geometric for two cells or a geometric original). The edges of the cut face are fitted to the TFI line the same way.

Split results: pipe O-grid along its axis and channel across its clustered direction keep every grid node to 3e-14 (a piece of a Vinokur distribution is again one). NACA 2412 C-grid, front block along the wall: wall nodes within 1.5e-8 chords; with chord-based spacings it was 1.2e-3, so the arc-length fit matters. Interior nodes follow the TFI of the new blocks, which is not the TFI of the whole: up to 2.2e-3 chords on the airfoil; on the elbow the projected wall face interiors slide along the torus by 4.8e-5 at radius 0.5.

Every edit checks the input, applies the change and checks the result; a refused edit writes nothing. Edits that change the shape (curve, vertex) remove the template's `reference`, since its analytic volume no longer applies. A chain of template, list, edits, a refused move and a grid replays from the journal with every result digest the same.

Next, each for a card: projection surface tools, merging blocks, collapsed (prism) blocks, and the FoehnCraft round trip for block meshes in the acceptance runner.

### FoehnCraft check in the acceptance runner (2026-11-24)

Card: "FoehnCraft check in acceptance". `cargo m1 --require-foehn` now builds pipe (8/10/40, wall spacing 5e-4), elbow (R = 0.05, bend 0.15, legs 0.1) and an extruded NACA 2412 C-grid (40/20/20, farfield 5) with TFI and with elliptic smoothing, and a pipe split through its chain with `block.split_block`. FoehnCraft reads every `foehncraft.vtu` (same cells and faces, volume and patch areas to 1e-12, no inverted face) and `mesh.msh` (volume and areas) back. The run fails if smoothing raises the max non-orthogonality FoehnCraft measures above the TFI grid's. First run with FoehnCraft 1.0.1: pipe 23.39° → 23.30°, elbow 23.77° → 23.58°, airfoil 60.54° → 47.00°; 56 files read and matched in all.

### MB1f: wall nodes slide along curved walls (2026-11-24)

The maintainer asked for this with multi-element airfoils in mind. Block vertices on an element's wall pin the grid lines there, so the grid kinks at them. With `wall_sliding: true` on the elliptic algorithm, the wall nodes of curved, extruded walls slide along the section, block vertices included.

- **Contours.** Wall face lines are taken along the curved direction of ruled wall faces, where one direction is bent and the other is straight. They are joined across blocks when the line turns by less than 30° and exactly two lines meet. The chain ends and the sharper corners (a trailing edge) stay fixed. Each contour is a dense polyline (32 samples per cell) of a spline through the TFI wall nodes.
- **Rule.** After each sweep, a node moves halfway to the foot point on the contour of its first off-wall node. That is a Neumann orthogonality condition. The Steger–Sorenson controls stay on.
- **Window (contours).** A node stays within half a TFI cell of its start, toward either neighbour, and 10 % of a cell away from its current neighbours.
  - With the foot point alone, the whole wall distribution drifted toward the convex leading edge, by 13 % of chord on the NACA 2412, and the grid folded.
  - Switching off the tangential wall controls on sliding nodes did not stop the drift. It also lost the orthogonality the controls give (mean wall angle 14.8°).
  - Blending the foot point with the midpoint of the neighbours pulled the clustering apart, and the largest angle rose to 21°.
  - A window of 0.25 gains less; 0.9 does not converge in 10,000 sweeps.
- **Measured.** On the NACA 2412 C-grid with 48 surface cells, the largest wall angle goes from 15.5° to 11.0° and the min scaled Jacobian from 0.58 to 0.62. On the acceptance C-grid, FoehnCraft's max non-orthogonality goes from 47.0° to 44.7° (TFI: 60.5°). The wall stays within 5e-6 chords of the section.
- **Projection surfaces (maintainer's card: keep working until they work).** The pipe cylinder and the elbow torus needed a different rule, measured in three steps:
  1. **Node by node.** The contour rule applied to every surface node, with a window in both directions, made the pipe wall angle go from 0.42° to 9.5° and folded the elbow. Halving the window gave 4.8° on the pipe, and the elbow still folded. On face interiors the Sorenson controls put the off-wall node over the wall node, so the foot point is wherever the node is. The node has no preferred place and drifts with every transient.
  2. **Block edges slide, face nodes follow.** Only the wall nodes on block edges slide to the foot point. Face nodes follow by transfinite interpolation of the rim displacement, at their TFI arc-length fractions, and are projected back onto the surface. The pipe holds (0.42° to 0.63°), but the elbow still folds. The follow nodes carry their off-wall nodes along, so the foot point of a rim node follows the rim node itself: a loop gain of one.
  3. **Anchor (kept).** A rim node now moves toward the mean of the foot point and its TFI position (`SURFACE_ANCHOR` 1). That makes the loop gain one half. Pipe and elbow then converge in 127 and 182 sweeps (fixed walls: 127 and 155), with the same wall angles (pipe 0.42°, elbow 0.59° to 0.58°) and Jacobians, and nodes on the surface to 1e-12. Anchor weights of 0.5 and 3 give the same.
  
  The same anchor on the section contours cuts the NACA gain (15.5° to 15.3° at weight 1, 14.8° at 0.25, against 11.0° without), so contours keep the window alone.

### MB1g: wall rings and layout optimisation (2026-11-25)

The multi-element template's first H-grid carried every wall spacing straight to the box and sheared the gap blocks: 76° to 89° of non-orthogonality, up to 240 cells over 70°, and two of nine test layouts refused. The maintainer chose "Ring + layout optimizer (paper)" on a card, general for every block topology. Each element now sits in a wall ring of six blocks, with a cove where a flap's nose sits under the main element, and the new `block::layout` optimiser moves the vertices, edge shapes and end spacings of the blocks off the rings toward orthogonal cells by a pattern search on the TFI grid (after Sauer & Morsbach 2023). All nine layouts then mesh with 47° to 68° and no cell over 70°. Design note, candidates and benchmark: [`block-layout-optimisation.md`](block-layout-optimisation.md).

### MB1h: triangle split (tri-center) and its trigger (2026-11-25)

The maintainer asked for the tri-center split as a general tool of the block mesher, with an automatic criterion for where it goes in, and chose "Also tetrahedron to 4 hexes" on the scope card. `block.split_triangle` turns a block with a flat corner into three blocks around a centre of valence 3. The new cell count runs as a row through the neighbouring blocks to the boundary. `block.optimize_layout` runs the MB1g optimiser on any extruded topology, together with the trigger. The trigger keeps a split only when the optimised grid measurably improves without new size jumps. `block::tri_center` builds a prism as three hexahedra and a tetrahedron as four, as building blocks for 3D templates. On a 171° triangle the worst angle falls from 62° to 26°. The nine airfoil layouts of MB1g have no triangle block left and stay unchanged. Design note and measurements: [`block-triangle-split.md`](block-triangle-split.md).
