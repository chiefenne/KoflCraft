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
