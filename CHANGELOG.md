# Changelog

## Unreleased

Nothing yet.

## Summary

| Version | Milestone | Highlights |
|---|---|---|
| 0.3.0 (2026-11-28) | Block-dominant meshing, first steps | Region proposal (`regions.propose`: heat-method wall distance, flow potential, skeleton, regions) and sweep fill (`regions.fill`: swept O-grids of hexahedra joined conformally to polyhedra, MD-25), with the elbow and the joined Manifold read back by FoehnCraft in the acceptance run; block T-junctions with a count solve, triangle and tet splits in 3D and a measured 3D trigger; ready-to-run programs for macOS, Linux and Windows on each release; README hero pictures. |
| 0.2.0 (2026-11-27) | Galleries and install guide | Example gallery (eight block-mesher and eleven polyhedral entries), recipes and an application gallery (heat exchanger shell side, seven-tube furnace, three-element airfoil), every entry with its exact commands, checked numbers and a MurmelCraft picture; install guide for macOS, Linux and Windows. No change to the mesher, its file formats or its tools. |
| 0.1.0 (2026-11-26) | M1 done, M2 to M5 and MB1 as prototypes | First release. Clipped Voronoi polyhedra in a box (M1, verified), convex STL meshing with size fields, geometry-derived sizes and boundary patch tools (M2), nonconvex ducts on a Rust scaffold with prism layers to the ports, named wall patches, flush openings, zero-thickness shells, polygon-prism layers and the solver rules of production codes (M3 to M5 prototypes), a block-structured mesher with templates (pipe, elbow, airfoil C-grid, multi-element airfoil), edit tools, elliptic smoothing and a layout optimiser (MB1), the blocked VTKHDF hand-off read by FoehnCraft and MurmelCraft, and one tool registry behind the CLI and the MCP server. **Every mesh checked before it is written; FoehnCraft reads every acceptance mesh back in CI; the seven-tube furnace (1.9 million cells) passes the solver rules.** |

## 0.3.0 (2026-11-28) – block-dominant meshing: regions, sweep fill and 3D block splits

**In one paragraph.** KoflCraft 0.3.0 takes the first steps towards
block-dominant meshing of real ducts. From a closed STL and its patch map,
`regions.propose` finds the duct's skeleton and marks each region as sweep,
junction, obstacle, box or remainder, and `regions.fill` meshes the sweep
regions with swept O-grids of hexahedra and joins them conformally to a
polyhedral mesh of the rest (MD-25): the elbow is hexahedra from port to
port, the Manifold one valid mesh of hexahedra in the pipes and polyhedra in
the box, and FoehnCraft reads both back in the acceptance run. The block
mesher gains T-junctions with an integer count solve, triangle and tet
splits that work in 3D, and a 3D trigger that keeps a split only when the
measured grid gets better. Releases now carry ready-to-run programs for
macOS, Linux and Windows.

### Highlights

- **Region proposal** (`regions.propose`). Wall distance by the heat method,
  the inlet-to-outlet potential, the skeleton with its end and branch
  points, regions marked sweep, junction, obstacle, box or remainder, and a
  second pass along each box's own long axis.
- **Sweep fill and join** (`regions.fill`). Swept O-grids along the flow
  potential with blended stations at mitre corners (elbow worst
  non-orthogonality 66.1° to 59.1°), the box behind the cuts meshed as
  polyhedra and joined node for node (Manifold: 132,520 cells, no inverted
  face, volume within 0.2 % of the STL), and more than two openings behind
  the cuts (a box with three round pipes: 105,638 cells, valid in
  FoehnCraft).
- **Block mesher in 3D.** T-junctions (BS-21) with an exact integer count
  solve, `block.split_triangle` on any topology, `block.split_tet` into four
  hexahedra, and a 3D trigger in `block.optimize_layout` that measures each
  candidate split.
- **Release downloads.** `kofl` and `mesh-mcp` for macOS (universal),
  Linux x86-64 and Windows x86-64 with `SHA256SUMS.txt`, attached to the
  release.

### Known limits

- Square pipe sections (the T-junction duct) still fail at the cut of the
  sweep join; round pipes work.
- `regions.propose` only proposes; `regions.fill` meshes the sweep regions
  and the box behind them, while junction and obstacle regions are left to
  the polyhedral mesher.
- The tube side and the solids of the heat exchanger still need the
  multi-domain mesher (MD-06 to MD-18), and there is no STEP reader yet.
- The Mac programs are signed ad hoc only, so macOS asks once before the
  first start.

### All changes

#### Region proposal and sweep fill

- **Region proposal for block-dominant meshing** (`regions.propose`). From a
  closed STL and its patch map, KoflCraft builds a tet scaffold, solves the
  wall distance (heat method) and the inlet-to-outlet potential, finds the
  skeleton with its end and branch points, and marks each region as sweep,
  junction, obstacle, box or remainder (left to the polyhedral mesher). A
  second pass tries each box and remainder region again along its own long
  axis. Writes a report, the fields (VTU) and a picture (SVG); nothing is
  meshed yet.
- **Sweep regions meshed with hexahedra** (`regions.fill`). The sweep
  regions of the region proposal are filled with swept O-grids: the
  stations are level sets of the inlet-to-outlet potential (the meridian
  planes in a bend, the mitre plane at a mitre corner), each section is the
  pipe template's O-grid with the first rows at the layer heights, and a
  sweep ends at its port or at a calm, exactly planar cut in front of the
  next region. The elbow is swept from port to port (17,472 hexahedra,
  volume within 0.003 %); in the Manifold both pipes are swept and the box
  is left for the polyhedral mesher. Writes the validated FoehnCraft VTU,
  a report and the mesh evidence.
- **Box joined to the swept pipes** (`regions.fill`, MD-25). The region
  behind the cuts of the sweeps is now meshed as polyhedra (median dual
  with the sweeps' wall layers) and joined conformally to the hexahedra:
  each cut gives the polyhedral mesher a planar cap, the pipe's last
  section takes the positions the mesher gives that cap, and every hex
  face at a cut is split into the polyhedral faces it holds. The Manifold
  is one valid mesh of 132,520 cells (69,632 hexahedra in the pipes,
  62,888 polyhedra in the box), with no inverted face and its volume
  within 0.2 % of the STL. `join: false` meshes only the sweeps, as
  before.
- **Polyhedral wall layers at box corners.** Where the wall offset kept a
  scaffold tetrahedron turned over at a corner of three walls, the layers
  of the third wall now yield too, so thick layer stacks no longer stop
  the mesh there.
- **Sweep fill in the acceptance run** (`regions.fill`, FC-06). The full
  CI run now meshes the elbow (hexahedra port to port) and the joined
  Manifold (swept pipes and the polyhedral box) and has FoehnCraft read
  both back: same cells, faces, patches and volume, no inverted face, no
  cut left as a patch. Each case's non-orthogonality on both measures is
  saved in `fill/quality.json`.
- **`regions.fill` names what it wrote.** The joined Manifold was told
  that its box was not meshed, and the elbow that it held polyhedra; the
  next step now says which of the three meshes the file is.
- **Gmsh file of the polyhedral Manifold readable again.** Four core
  cells had a corner as their kernel point, so their tetrahedral fan in
  `mesh.msh` held tetrahedra of zero volume and FoehnCraft refused the
  file. Those cells now get a kernel point strictly inside.
- **Blended stations at mitre corners** (`regions.fill`). Next to a mitre
  plane the sweep's stations are now planes that turn step by step from
  the leg's section to the mitre plane, the classical structured mitre,
  instead of level sets that crowd at the inner corner. On the elbow the
  worst non-orthogonality drops from 66.1° to 59.1° (FoehnCraft's measure
  from 81.9° to 59.1°) and the smallest scaled Jacobian rises from 0.26
  to 0.37. The report counts the blended stations.
- **More than two openings behind the cuts** (`regions.fill` with join).
  The polyhedral mesher now gives each connected cap of a port id its own
  layer similarity, so several cuts can share one id, and the join tells
  them apart by their planes. A box with three round pipes joins into one
  closed, valid mesh of 105,638 cells, and FoehnCraft reads it as valid.
  A cut that lies next to a ring of STL nodes no longer leaves a strip of
  needles that the scaffold could not recover. Square sections (the T)
  still fail at the cut.

#### Block mesher

- **Triangle split in 3D.** `block.split_triangle` now works on any block
  topology, not only extruded ones. The flat corner is a block edge (new
  input `along`, or the flattest edge at `vertex`): every block stacked
  along that edge becomes three blocks around a centre line, and the new
  cell count runs as a sheet through the neighbouring blocks to the
  boundary, so all faces stay point-matched. Blocks in any local frame and
  any direction in space are handled.
- **T-junctions in block topologies.** A block corner may now sit inside an
  edge of a bigger block (`junctions` in the topology file, BS-21). The big
  block's face is gridded first and the smaller blocks take its nodes, so the
  hexahedra still meet node for node. `block.check` lists the junction edges
  and the faces they split.
- **Count solve.** One integer per chain of opposite edges, with the sum rule
  at every T-junction, fitted to target counts by an exact integer lattice
  search (`block::counts::fit`).
- **Triangle split that stops at the neighbours.** `block.split_triangle`
  takes `"stop": "t_junction"`: the neighbouring blocks stay whole with a
  T-junction on their edges and carry the new cells in their counts, instead
  of being cut in two all the way to the boundary. On the two-layer wedge:
  8 blocks and 126 cells against 10 blocks and 150 cells for the sheet split.
- **3D trigger of the prism split.** `block.optimize_layout` on a topology
  that is not extruded along z (or has T-junctions) now finds prism-shaped
  blocks by geometry, tries the split with the T-junction stop and with the
  sheet, measures the 3D grid (non-orthogonality, growth, scaled Jacobian)
  and keeps a split only when its blocks get at least 2° better with a worst
  scaled Jacobian of 0.4 or more, while the whole grid gets no worse and
  grows by at most 10 %. The report (`mode: 3d`) names every split kept and
  every candidate left, with the reason.
- **Tet split of a block.** `block.split_tet` splits a block shaped like a
  tetrahedron into the four hexahedra of the tetrahedron. A block on one of
  the three split triangles is split into three with it (a prism column), so
  the cells still meet node for node, and the new cells run on as one row,
  cut as a sheet or stopped with T-junctions. The 3D trigger of
  `block.optimize_layout` tries tet-shaped blocks too and keeps a split only
  when it measurably helps (on the test tet: worst scaled Jacobian 0.04 to
  0.49, worst angle 82° to 51°).
- The triangle split now puts its centre at the mean of the three points on
  the triangle's sides, which keeps its three quadrilaterals convex.

#### Requirements

- **Requirements v1.2.** MD-25 (block next to poly) now gives the
  polyhedral scaffold each interface quad split along its Delaunay
  diagonal instead of about its centre, as the sweep join builds it: the
  centre split stalls the scaffold's boundary recovery at thin wall rows.

#### Releases and documents

- **Release downloads**: the new workflow Executables builds `kofl` and
  `mesh-mcp` for macOS (one universal build for Apple silicon and Intel),
  Linux x86-64 and Windows x86-64 and attaches them with `SHA256SUMS.txt`
  to the release; the Release workflow starts it after each release, and it
  can be run by hand for an existing one (docs/release.md). The README's
  "Download a release" says which file to take and how to start the
  unsigned Mac programs the first time.
- **README hero pictures.** Six square pictures on white, like FoehnCraft's
  hero row, each linked to its gallery or application entry: prism layers
  around a furnace tube, the NACA 2412 C-grid, the Manifold bend in
  non-orthogonality, the exploded elbow blocks, the heat exchanger shell side
  and the three-element airfoil blocks. They replace the Manifold slideshow,
  which stays on the tour page. `tools/readme-hero/make_heroes.py` renders
  them again with MurmelCraft from the gallery runs' meshes.

## 0.2.0 (2026-11-27) – example gallery, application gallery and install guide

**In one paragraph.** KoflCraft 0.2.0 shows what the mesher does and how to
run it. The new [example gallery](docs/gallery.md) has eight block-mesher
and eleven polyhedral entries, the [recipes](docs/recipes.md) the steps
behind them and the [application gallery](docs/applications.md) three
engineering cases (the shell side of a shell-and-tube heat exchanger from
its STEP file, the seven-tube furnace of FoehnCraft's gallery and a
three-element airfoil). Every entry carries the exact commands that made it,
its checked numbers and a picture made headless by MurmelCraft, and a script
in `tools/gallery/` makes every entry again. The README has an install
section and [docs/install.md](docs/install.md) a step-by-step guide for
macOS, Linux and Windows. The mesher itself, its file formats and its tools
are unchanged from 0.1.0.

### Highlights

- **Example gallery and recipes.** Block meshes (channel, pipe O-grid and
  split, elbow, NACA C-grid, multi-element airfoil, triangle splits) shown
  with their blocks coloured and exploded; polyhedral meshes from a box to
  the Manifold, the furnace and a pipe with a baffle, a plate and a probe.
- **Application gallery.** Heat exchanger shell side (177,994 cells),
  seven-tube furnace (222,715 cells) and three-element airfoil, each
  stopping at the checked mesh.
- **Install guide.** Prerequisites and commands per system; Gmsh is
  optional.

### Corrections

- The 0.1.0 notes give the seven-tube furnace as 1.9 million cells. Since
  the polygon-prism layers of BL-11 (#69, part of 0.1.0) it builds with
  222,715 cells; the README now says so.

### Known limits

- Unchanged from 0.1.0. The tube side and the solids of the heat exchanger
  need the multi-domain mesher (MD-06 to MD-18), which is not built; the
  STEP file is read through Gmsh in the gallery script, since KoflCraft has
  no STEP reader yet.
- No prebuilt binaries yet: build from source.

### All changes

#### Documents

- **Install section in the README.** How to build KoflCraft from source
  (prerequisites, clone, build, the three programs in `target/release/`,
  `cargo install`, update), with a placeholder for release downloads. The
  minimum Rust is corrected from 1.85 to 1.88 (`sha2` 0.11 needs 1.88).
- **Install guide for macOS, Linux and Windows.** `docs/install.md` lists
  the tools with the install commands per system (Xcode tools; apt and dnf;
  winget and Chocolatey), makes clear that Gmsh is optional (only the
  `gmsh -check` of `cargo m1 --require-gmsh`; KoflCraft writes `.msh` files
  itself) with the Gmsh install per system, and how to point the
  acceptance run at Foehncraft. The README's install section gives the
  short version per system.

#### Galleries and recipes

- **Example gallery, recipes and application gallery.** `docs/gallery.md`
  with a page per topic (first: eight block-mesher entries, each with its
  `kofl call` commands, checked numbers and a MurmelCraft picture of the
  blocks coloured and exploded), `docs/recipes.md` with the matching
  recipes, `docs/applications.md` (entries to come), linked from a bar at
  the top of the README. `tools/gallery/block_gallery.py` makes every
  entry and picture again; the README slideshow stills moved to
  `docs/gallery/readme-tour.md`.
- **Polyhedral gallery.** `docs/gallery/polyhedral.md` with eleven
  entries (a box, a layered channel, analytic sweeps, a convex STL with a
  size field, an elbow in three modes, the Manifold by hand and from a
  bare STL, the furnace with one and seven tubes, yielding layers and a
  hex core, and a pipe with a baffle, a plate and a probe made tool by
  tool), each with its commands, checked numbers and a MurmelCraft
  picture; `docs/recipes/polyhedral.md` with six recipes;
  `tools/gallery/poly_gallery.py` makes them again.
- **Application gallery.** `docs/applications.md` with three entries, each
  with its commands, checks and a MurmelCraft picture: the shell side of the
  shell-and-tube heat exchanger from its STEP file (the fluid prepared with
  Gmsh as a separate process, then ports, setup and the polyhedral duct in
  KoflCraft: 177,994 cells), the seven-tube furnace of FoehnCraft's gallery
  259 (222,715 cells) and the three-element airfoil. Every entry stops at
  the mesh. `tools/gallery/app_gallery.py` makes them again.
- **README.** The status block points to the galleries, the furnace cell
  count is corrected, the build commands use `kofl`, Leif Kobbelt is named
  with the graphics principle and the install guide is in the
  documentation table.

## 0.1.0 (2026-11-26) – first release: polyhedral and block-structured meshes for FoehnCraft

**In one paragraph.** KoflCraft 0.1.0 is the first tagged release of the
mesher. Milestone M1 is done and verified: clipped Voronoi polyhedra in a
box, validated cell by cell, exported to VTU, Gmsh MSH 4.1, OpenFOAM and
FoehnCraft's own hand-off, and driven from the same engine by the command
line and an MCP server. Beyond M1 the release carries working prototypes of
the later milestones, each with its grade and its gaps on the
[status page](docs/status.md): convex STL meshing with size control (M2),
nonconvex ducts with prism layers carried out to the ports (M3 to M5) and a
second, block-structured mesher (MB1). No result is presented as more than
it is: 9 features are graded verified, 24 partial, 1 unverified.

### Highlights

- **Polyhedral core.** Voronoi cells clipped by their nearest sites with a
  security-radius completeness test (about 25 times faster than all pairs at
  16,000 sites, the same cells), exact side predicates and symbolic vertices
  so neighbouring cells always agree, size fields from datum objects and
  coordinate regions with density-weighted CVT seeding, and sizes derived
  from wall curvature, gaps between walls and feature edges.
- **Nonconvex ducts.** A Rust volume scaffold (wall remesh to a chord
  tolerance, Delaunay tetrahedralization with exact predicates and boundary
  recovery) replaces Gmsh on the Manifold duct; a Voronoi or median-dual
  core with prism layers that follow the core cells and end on the port
  planes; named wall patches and flush openings, walls without layers and
  walls whose layers yield, zero-thickness baffles and plates, polygon-prism
  layers (BL-11), grading on the duct path and an opt-in hex core.
  Validity follows the rules of production finite-volume codes (closed
  shell, positive volume and face pyramids, non-orthogonality below 90°).
- **Boundary patches from a bare STL.** Ten `patch.*` tools, port detection,
  patch proposals (planes, cylinders, spheres, openings, box sides) and a
  FoehnCraft boundary-condition skeleton with the type of every role.
- **Block-structured mesher (MB1).** Topology file with checks, edge laws,
  transfinite interpolation, Thompson–Thames–Mastin elliptic smoothing with
  wall orthogonality and sliding wall nodes, templates for the pipe, the
  elbow, the airfoil C-grid and the multi-element airfoil, edit tools,
  triangle splits and a layout optimiser that cuts the worst
  non-orthogonality of the test layouts from 76°–89° to well below.
- **Hand-off.** `foehncraft.vtkhdf` (blocked VTKHDF agreed with FoehnCraft
  and MurmelCraft) is the default, `foehncraft.vtu` on request, with
  `mesh-evidence.json` and BLAKE3 hashes of every file, streamed to disk as
  written. `cargo m1 --require-foehn` has FoehnCraft read every acceptance
  mesh back in CI.
- **Agent-ready.** One tool registry (dotted ids, JSON Schemas, diagnostic
  codes) behind `kofl` and the MCP server with read and write roots,
  journal and replay, exact estimates before a build, approvals for large
  builds and a generated status page of verification grades.

### Known limits

- Nonconvex ducts are prototypes: the scaffold keeps the ports as they come,
  some CAD-like parts still leave a few cells that fail the final checks
  (the Voronoi core hardening is paused, see `research/voronoi-symbolic.md`),
  and the flow accuracy of the duct meshes is not yet tested.
- Several regions in one mesh (multi-domain, MD-01 to MD-27), STEP input and
  sliding interfaces are decided but not built.
- No prebuilt binaries: build from source with `cargo build --release`.
- Every gap of every feature is listed on the [status page](docs/status.md).

### All changes

### Polyhedral core: Voronoi cells, sizes and exact geometry

- Voronoi cells are clipped by their nearest sites only (ALG-04): a k-d tree gives the neighbours in order of distance, and a security-radius test (Ray, Sokolov, Lefebvre & Lévy 2018) proves each cell complete before it stops (ALG-13). Shared simple vertices are recomputed from the sites and planes that define them, so neighbouring cells agree to the last bit. The all-pairs clipper stays as the reference (`Strategy::AllPairs`) and a test compares both. About 25× faster at 16,000 sites; 256,000 cells build in about a minute (`research/alg-04-neighbour-clipping.md`, `examples/clipping_benchmark.rs`). Algorithm ids are now `..._v2`. Size fields may ask for 200,000 seeds (was 20,000) and the MCP site limit is 20,000 (was 256); inline site lists over MCP are still bounded by the 256 KiB request size.
- Size control for convex STL meshing (ALG-05, GRD-01 to GRD-03, GRD-05, GRD-08): an STL config can give a `size_field` instead of `sites`. Datum objects (point, segment, plane, box, slab, sphere, cylinder, cone) and coordinate regions (box, slab, cylinder in named frames) set sizes; the field is a closed-form, gradation-limited minimum. Seeds follow it by density-weighted CVT, chosen over DistMesh springs in `research/m2-size-control.md`. Results and artifacts carry the achieved sizes per source; `mesh.estimate` counts the seeds exactly; new diagnostic MESH007 `invalid_size_field`.
- Geometry-derived size sources (GRD-04): `curvature` (one cell spans at most `angle` of the wall normal's turn; normal-cycle curvature tensor, Cohen-Steiner & Morvan 2003), `gap` (at least `cells_across` cells between opposite walls; medial-ball thickness, exact for convex domains) and `feature_edges` (a size along edges sharper than `feature_angle`). They keep the closed-form gradation limit as weighted surface primitives in a pruned hierarchy. The `sizing` evidence adds per derived source the method, the pieces held at `min_size` and, for gaps, the cells counted across along the wall normals of the built mesh; shortfalls are listed as limitations. Benchmark against cotangent and normal-deviation curvature and the ray-cast thickness in `research/grd-04-geometry-sizing.md` (`examples/geometry_sizing_benchmark.rs`).
- Polyhedral duct core (`volume_clipped_voronoi_cvt_v1`, now the default): a Voronoi diagram restricted to the scaffold tetrahedra, relaxed by ten density-weighted Lloyd iterations (Yan, Wang, Lévy & Liu). Elbow 1,355 cells instead of 7,005, Manifold 92,058 instead of 417,457. Both stay valid in `foehn mesh`. Max non-orthogonality is 85.9° and 87.7° (was 87.0° and 89.6°). Manifold faces above 70° drop from 3,038 to 605. Core sphericity p05 rises to 0.78 (elbow) and 0.81 (Manifold). At reentrant wall edges, sites are split into mirrored pairs and the diagram is rebuilt. A repaired cell keeps one face per neighbour, slivers are merged, and a layout FoehnCraft would reject is rebuilt with one more Lloyd iteration (`retries`). `"core": "median_dual"` in `volume_scaffold` keeps the old core. `lloyd_iterations` sets the relaxation, at most 200. The hand-off now merges the boundary faces one cell has on one patch, and the Gmsh export splits polygons with more than four nodes. `mesh.estimate` counts the Voronoi core, and `cells_basis` notes the cells added at reentrant edges. Note: `research/m4-voronoi-core.md`.
- Voronoi core, consistent cuts (ALG-16): in the final diagram every vertex of a cell piece is named by the scaffold nodes and sites that define it, and its coordinates and its side of a further bisector follow from the name, so two neighbouring scaffold tetrahedra can no longer disagree about a thin piece. This was the cause of "piece of site S in tetrahedron T has no continuation" on thin CAD-like parts. All points are welded once, globally and independent of order. A film left by the weld joins its largest neighbour, a needle wall polygon joins the cell's other polygon on that wall triangle, and a pinched ring is cut instead of dropped. The report's `voronoi_core` gains `welded_vertices`, `merged_films`, `merged_wall_rings` and `turned_wall_rings`. Manifold stays valid with 95,938 polyhedra (95,666 before). Jacket-like parts now pass the core but still stop in the checks of the wall columns. [Note](../research/voronoi-symbolic.md)
- Voronoi core, exact sides (ALG-16): the side of a named vertex of a cell piece relative to a further bisector is now the exact sign of a quotient of two determinants of squared distances (floating-point filter, then expansion arithmetic), after the side predicates of Lévy's PCK. Both cells of a face and every tetrahedron take the same decision, with no tolerance. Wall seeds on the ball of one wall vertex share a stand-in distance from that vertex, so their designed ties are exact zeros. The security radius uses a safe upper bound of each vertex's distance instead of its rounded coordinates. Elbow and Manifold are unchanged in their counts. On a jacket-like part the open cells are gone (3 pairs before); the run still stops in the last checks (two faces with a fan triangle seen from behind by a rounding amount, three wall columns with an inverted sub-wedge, one volume mismatch). [Note](../research/voronoi-symbolic.md)
- Voronoi core hardening paused: the research note records the 56 cells of a jacket-like part that still fail the final checks (18 by far more than rounding), the cause traced for the wall columns (exact names, rounded coordinates, true edges shorter than the rounding) and the open options. Earlier counts of 2, 3 and 6 cells were too low. No code change. [Note](../research/voronoi-symbolic.md)

### Nonconvex ducts: Rust scaffold, wall layers and solver rules

- Rust volume scaffold (OPEN-07, M4): `volume_scaffold.build` in a duct partition remeshes the STL wall to a chord tolerance (Botsch–Kobbelt loop with the Dunyach curvature target, features and ports kept) and tetrahedralizes it in Rust (Bowyer–Watson with exact predicates, conforming boundary recovery, Delaunay refinement, sliver removal by edge removal). The Manifold duct no longer needs Gmsh; its polyhedral wall lies within 0.050 of the STL (Gmsh scaffold: 0.543). The report names the scaffold source and the measured wall deviation.
- Scaffold boundary recovery converges on regular CAD walls: the wall remesher ends with Delaunay flips between nearly coplanar facets (within the chord tolerance), and recovery splits edges at an input vertex on concentric shells. On the synthetic cooling jacket the recovery diverged and now converges in 12 passes; Manifold needs no recovery points (229,343 cells, smallest 5.6e-3). Tilted port caps from single-precision STL are accepted (port planarity to 1e-6 of the extent). Layer front collisions name their location.
- Gmsh export of polyhedral ducts: when two faces of a cell share more than one edge, one of them is fanned about its centre (in both cells that hold it), so no two subdivision elements share two faces; FoehnCraft refused the Manifold MSH otherwise.
- Thin gaps and thin walls (GRD-04, wall remesh, layers): the wall remesher takes a local-thickness size source for surfaces of any shape (`volume_scaffold.build.thickness`: shrinking-ball medial thickness on the fluid and the solid side, one edge across by default), runs until no edge is far below its target, and keeps a clear zone around sharp edges. `mesh.setup_duct` turns the thickness source on and never asks for a chord tolerance below the STL's own facet error. Prism layers are squeezed locally, with the same layer count, where they would fill more than 30% of a gap, where fronts collide, or where the wall offset would turn a scaffold tetrahedron over. `KOFL_TRACE=1` prints progress lines (counts and lengths, no coordinates). New fixture `cooling-jacket-thin.stl` (0.5 mm gap and 0.5 mm rib). The Voronoi core does not yet complete on the jackets.
- Polyhedral duct core (`volume_median_dual_v2`): median-dual cells of wall vertices on reentrant wall edges are split along scaffold faces into one part per wall plane, so no face between a wall prism and the core has S·d ≤ 0. FoehnCraft's `foehn mesh` now accepts the polyhedral elbow (10 inverted faces before) and Manifold (18 before). `cargo m1 --require-foehn` also reads the Manifold ducts back and fails on any export FoehnCraft calls invalid. `mesh.estimate` labels polyhedral duct counts in `cells_basis`, since the split adds a few cells after the estimate. Note: `research/m4-reentrant-dual-split.md`.
- Polyhedral ducts: the wall layers now follow the Voronoi core cells (wall columns, `research/m4-wall-columns.md`). Each core cell's wall polygon in a wall triangle is extruded through the layers at fixed barycentric coordinates, and the polygons of one core cell over coplanar wall triangles merge into one column. Over curved walls no layer cell spans two triangles, because FoehnCraft fans a bent face about its node average. FoehnCraft reads both fixtures as valid with exact volumes and no inverted faces. Max non-orthogonality is 80.1° on the elbow (was 85.9°, faces above 70° 44 → 2) and 82.8° on Manifold (was 87.7°, 605 → 67). Manifold has 162,435 cells instead of 92,058, and sliver layer cells remain on its curved walls (smallest cell about 5e-11). `"wall_join": "split_cap"` in `volume_scaffold` keeps the earlier join. The hand-off merges a cell's boundary faces on one patch only when the result is flat. `mesh.estimate` counts two layer cells per wall triangle and layer for columns and labels `cells_basis` approximate. The Gmsh export writes each wall column as its native prisms. A Clippy lint of the current stable toolchain in `mesh-acceptance` is fixed.
- Polyhedral ducts: the columns wall join now builds a valid Manifold mesh on the Rust scaffold (0 inverted faces in FoehnCraft, exact volume). A merged wall column whose footprint is nonconvex is split back into its units when FoehnCraft's test (S·d) between it and a neighbour fails or comes within 5° of failing on any layer. Before, every layout was rejected. Core faces merge only when every vertex of every merged ring lies in the plane, not only the outer loop, and only to the tolerance of the shell check. A tent of rings with a flat rim was flattened, and the core cell lost volume (MESH006 on the third layout). The hand-off's merge of boundary faces has the same all-vertex test. Manifold keeps `"wall_join": "split_cap"` for now, because columns there still have sliver cells (volume ratio about 1e13).
- Polyhedral ducts: the wall columns no longer leave sliver cells on curved walls. Wall triangles at edges where columns end (curved walls and reentrant edges, not sharp convex corners) each get a fixed Voronoi site on the balls of their three corners (conforming wall seeds, after VoroCrust; `research/m4-wall-columns.md`), so the core cell over such a triangle meets the wall in exactly that triangle. The other sites keep 1.5 ball radii away, wall-vertex sites without a ball stay on their march line, and a column below 5% of its wall triangle joins a neighbour in its plane. Two cells that met in two polygons round a reentrant wall edge are split. FoehnCraft on Manifold with columns: 239,347 cells (before 505,937), 0 inverted faces, exact volume, maximum non-orthogonality 74.5° (87.9°), 3 faces above 70° (412), skewness 1.44 (2.59), volume ratio 1,882 (1e13), smallest cell 1.1e-3 (9e-13), valid at the first layout in 130 s. Elbow: 1,004 cells, 56.6° (80.1°), volume ratio 77 (727). The Manifold fixture now uses the columns join (split cap: 89.1°, 1,931 faces above 70°, volume ratio 285); `"wall_join": "split_cap"` stays available. The core cells at seeded walls are lower: core sphericity on Manifold is 0.74 at the median (0.85), on the elbow 0.62 at the minimum (0.73).
- Mesh check: a long thin planar polyhedron face far from the origin is no longer rejected as nonplanar. The face normal is summed about the face's first vertex, and offsets are measured from the least-squares plane. The report gains `wall_seeds`, `seed_flips`, `unseeded_wall_triangles`, `torn_contacts`, `absorbed_columns` and `collapsed_pieces` under `voronoi_core`.
- Wall columns, nonconvex wall polygons: a wall polygon of a core cell is triangulated by ear clipping instead of a fan from its first corner, with straight corners put back by an edge split, so a polygon with a notch or a speck no longer gives inverted layer sub-wedges. The patch area check uses the vector area of a boundary face, which is the true area of a planar face that is not convex. On a jacket-like part all layer cells now pass (23 sub-wedges were inverted); the run still stops in the last check on a few core cells (fan triangles seen from behind by a rounding amount, three missing faces between cells). The checker's planarity limit for polyhedron faces goes from 1e-10 to 1e-7 of the face size (maintainer's choice): the core welds points up to 1e-10 of the model size apart, and a face within the limit is fanned from its vertex mean as FoehnCraft does. Elbow and Manifold counts are unchanged. [Note](../research/voronoi-symbolic.md)
- **Solver rules decide duct validity** (decision 2026-11-25): a duct mesh is valid when every cell has a closed shell, positive volume, positive face pyramids on its centroid and faces at less than 90° non-orthogonality, as in production finite-volume codes. Flat faces are no longer required (`maximum_warping` reports the largest warping angle). Patch areas and volume against the STL, face pairing and winding are still checked.
- **Median dual is the default polyhedral core**, as the dual of a tet mesh in production polyhedral meshers. A dual cell that breaks the solver rules is split into valid parts: its convex tet-corner pieces, regrown while the rules hold and each face shared with a neighbouring cell stays one polygon. Elbow: 7,011 cells (5 split), 0 invalid. Standard cooling jacket at size factor 2: 1,194,032 cells, 0 invalid (it failed with 58 cells on the Voronoi core).
- The Voronoi core (`"core": "voronoi"`) is on hold: its wall-join cells can break the solver rules (3 on the elbow with columns, 1 with split caps), and the build rejects them.

The solver rules now judge the mesh as it is exported: the faces between two
cells merged into one polygon and the cell centroid taken exactly from the
fan-triangulated shell, as the hand-off and FoehnCraft take them. Cells are
split against that same view. The elbow median dual splits one cell
(238 polyhedra); the Voronoi columns join passes there again.

Cells that still fail on the exported mesh (against a neighbour) are split
again into halves by volume until the mesh passes. The cooling jacket
(1.19 million cells) now builds with 0 failing cells and the export finds 0
faces turned away. When layer cells fail, the message says that the wall
layers collide or fold.
- Solver-rule count beside the strict duct checks: every duct mesh is also counted against the validity rules of production finite-volume codes (closed shell, positive volume, positive face pyramids, faces at less than 90° non-orthogonality, warped faces allowed, geometry from the vertex-mean triangle fan). `Measurements::solver_rules` holds the counts and `KOFL_TRACE=1` prints them. Report only; the strict checks still decide.
- **Named wall patches and flush openings** in the polyhedral duct: a patch map may name any number of wall patches. Each one keeps its name and role in the report, the Gmsh physical groups and the FoehnCraft hand-off. An inlet or outlet that sits flush with its wall grows wall layers instead of contracting in its plane. New fixture: the box furnace from FoehnCraft gallery 259.
- **Wall patches without layers** (`no_layers` in the partition): the layers of the walls around run straight into a planar patch that grows none, so on the furnace the tube layers end on the side walls as rings around each tube. The patch's own facets slide in its plane to make room. The furnace fixture names its side walls `sides`.
- **Walls whose layers yield** (`layers_yield` in the partition): a planar wall keeps its layers, and the walls that cross it run their own layers straight through them down to it, with no mitre. On the furnace, `"layers_yield": ["sides"]` keeps the side walls' layers and lets the tube layers win at the junction: a ring of (tube layer x side layer) hexahedra around each tube.
- **Rim vertices of a tube keep its polygon**: every vertex of a source feature line is now kept as a corner of the wall remesh, so the rim of a tube (or any polygonal opening) follows the STL's segments and no long chord cuts across the hole. On the furnace the layer rings at the side walls sat up to 4.5 mm inside the tube circle; they now stay on the polygon.
- **Polygon-prism wall layers** (BL-11): with the median dual core, the six sub-wedges per wall prism and layer are merged into one polygon prism per wall vertex sector and layer, as in production dual meshers; the sub-wedges stay the cell's exact prism decomposition (the Gmsh file still writes them as prisms). Columns on curved walls whose merged cell would fail the solver rules start as corner prisms, and any cell that fails as exported is split back. The standard cooling jacket has 186,716 cells instead of 1,187,199 (all pass the solver rules) and builds in 53 s instead of 57 s; Manifold has about 138,000 instead of 1.08 million. The report gains `layer_polyhedra`, the normal thickness ranges now include prism polyhedra (they were empty for those layers), and the cell estimate counts the polygon prisms. The Gmsh writer, the hand-off, the re-check after splits and the face pairing of the final check are faster.
- Zero-thickness shells in the scaffold (SHL-01 to SHL-03): `shells: [{name, polygon}]` in a duct partition puts each baffle or plate into the tet scaffold of the median-dual core. The wall is cut along the shell before the remesh, the shell is meshed at the local size and the scaffold is split into two sides, so each shell comes out as the wall patches `<name>` and `<name>_back` with one owner cell per face; free edges inside the fluid stay closed and the wall layers run into the shell. The build report lists each shell with its areas, faces and free-edge length. The Voronoi core keeps the exact cut of its convex cells. A shell over a port cap or touching another shell is refused with its name for now.
- Shell tools and shell sizes (SHL-01): `shell.add`, `shell.list` and `shell.remove` edit the shells of a duct partition by name (a repeated name, a nonconvex polygon and a shell touching another are refused; with the STL, `shell.list` shows the default sizes). Each shell has its own cell `size` (default: an eighth of its narrow width, capped by the build's `max_size`) and `edge_size` on its free edges (default: half the size), both in the build's size field, so the wall line and the tetrahedra around the shell are refined too. The shell is graded from its wall-line edges and splits encroached boundary edges (Ruppert), wall-line edges once with their wall triangles: on the pipe test the plate goes from 12 to 390 triangles and the baffle's smallest angle from 4.5 to 15.4 degrees. The report gives each shell's sizes, where they came from and the achieved edge lengths.
- Grading on the duct path (GRD-02, GRD-03, GRD-08): `volume_scaffold.build.grading` refines the wall and the core with datum objects and coordinate regions, the report gives each source's achieved cell sizes, and `grading.add`, `grading.list` and `grading.remove` edit the sources of a partition by name.
- Hex core with layers (prototype, opt-in): `"core": "hex_core"` in `volume_scaffold` fills the bulk of a duct with the cubes of a 2:1 balanced octree and the gap to the wall layers with the median dual, joined face to face (hanging-node faces at size jumps). On the furnace at an interior size of 0.0625 it builds in 26 s instead of 42 s with 7% fewer faces but 3% more cells, so it stays out of the region menu for now ([note](research/hex-core.md)).
- **Faster duct build with the median dual core**: the repeated split of cells that fail the solver rules on the exported mesh no longer checks the whole mesh after every round. The face table is kept up to date and only the cells within two face neighbours of a split are checked again, with the same verdicts as a full check. The full check itself runs on all cores with no allocation per face. The standard cooling jacket at size factor 2 builds in 64 s instead of 271 s (4 cores) and its files are byte-identical to before (1,187,199 cells, 0 failing).
- **Faster duct export and final check**: the volume, Gmsh, patch and report files are written side by side while the FoehnCraft hand-off mesh is built, and the element quality measures of the final check run on all cores. The standard cooling jacket at size factor 2 builds in 56 s instead of 64 s; every mesh file is byte-identical. `KOFL_TRACE=1` now also prints the wall-clock seconds of each stage.

### Boundary patches

- Boundary patch tools (BC-01, BC-02): ten registry tools `patch.init`, `patch.list`, `patch.select_plane`, `patch.select_box`, `patch.select_normal`, `patch.grow` (feature-angle region growing), `patch.split`, `patch.merge`, `patch.rename` and `patch.set_role` define the boundary of a closed STL (convex or not) step by step. Each reads a `koflcraft.patch-map` and writes a new one; every read checks that the map belongs to the STL and gives every facet exactly one patch (BC-06), and `patch.list` reports area, mean normal, bounds, regions and neighbours per patch. `generate-stl` and `inspect-stl` take `--patches MAP` (tools: `patches`): the map's names, which are FoehnCraft boundary names, name the faces in MSH 4.1, the FoehnCraft VTU and OpenFOAM, and its roles are saved with the artifact (BC-04). New codes PATCH001 to PATCH004, new MCP limit `max_patch_facets` (200,000). `cargo m1 --require-foehn` reads a patched box in all three formats with `foehn mesh`. Method and checks: `research/m2-patch-tools.md`.
- Bare STL to duct mesh (BC-03, ports): `patch.detect_ports` finds the planar end caps of a closed STL's pipes and writes a patch map with `walls`, `inlet` and `outlet`, with the ranked candidates as evidence; `mesh.setup_duct` writes a polyhedral partition sized from the port radius (`size_factor`); `generate-duct` and `mesh.generate_duct` read a `koflcraft.patch-map` in place of a port file, so the ports need not be contiguous facet ranges. New synthetic fixture `test-fixtures/cooling-jacket`.
- Patch proposals and boundary-condition skeleton (BC-03, BC-05): `patch.propose` splits a bare closed STL at feature edges into planes, cylinders, spheres and freeform faces, marks openings and bounding-box sides and writes a patch map to rename and give roles; `bc.skeleton` writes a FoehnCraft case with the type of every role and TODO comments in place of physical values.

### Block-structured mesher (MB1)

- Concept and requirements for a block-structured mesher variant: `research/block-structured-concept.md`, requirements v1.0 section 4c (BS-01 to BS-19), OBJ-21 and milestone MB1.
- Block-structured mesher, first step MB1a (OBJ-21, BS-01..BS-14 in part):
  the block topology file `koflcraft.block-topology` 1.0, a check of
  indices, handedness, edge curves and laws, count chains, faces and
  patches; uniform, geometric and two-sided tanh edge laws; lines and arcs;
  transfinite interpolation of faces and blocks; point-matched hexahedra
  built without welding; validation with corner Jacobians and conservation;
  VTU, MSH 4.1, FoehnCraft VTU and evidence with block ranges; templates
  `channel` and `pipe_ogrid`; tools `block.template`, `block.check`,
  `block.grid` and the diagnostic BLOCK001.
- MB1 runs as a parallel track next to M5 (requirements section 4c).
- Block mesher MB1b: `block.grid` takes `algorithm: {"kind": "elliptic"}` for Thompson–Thames–Mastin smoothing of every node off the boundary, across block faces, with wall orthogonality and first-spacing control at wall patches (BS-09); wall orthogonality without a wall patch is refused (BS-12). The report gives the wall angles before and after.
- Block mesher MB1c: `block.template` writes an `elbow` (pipe O-grid swept through a bend with optional straight legs) and an extruded `airfoil_cgrid` (NACA 4-digit or a Selig/Lednicer coordinate file, blunt trailing edges closed). The topology file gains spline edges and `surfaces` (cylinder, torus) that boundary faces are projected onto (BS-07); the handedness check follows curve tangents.
- Block mesher MB1d: the elliptic smoother freezes the TFI residual so only the wall controls move nodes, moves interior block vertices, lets nodes on planar `symmetry` patches slide in their plane (an extruded airfoil with one span cell is now smoothed), and solves grid lines with line Gauss–Seidel. The elbow no longer gets a worse cell than TFI (FoehnCraft max non-orthogonality 23.45° against 23.7°). The default iteration cap is 10,000 and the `block.grid` report gains `moved_vertex_nodes`, `sliding_patches` and `moved_sliding_nodes`.
- Block mesher MB1e: edit tools for a block topology (BS-04). `block.list` shows vertices, blocks, edges, patches and chains. `block.set_count`, `block.set_distribution`, `block.set_edge_curve`, `block.move_vertex`, `block.set_patch`, `block.add_vertex`, `block.remove_vertex`, `block.add_block`, `block.remove_block` and `block.split_block` each write a new checked topology or refuse the edit with nothing written. A split carries through the chain of opposite edges and keeps the grid nodes on block edges; the steps replay from the journal.
- Acceptance (FC-06 for block meshes): `cargo m1 --require-foehn` builds the pipe, elbow and airfoil C-grid block templates with TFI and elliptic smoothing, plus a pipe split by `block.split_block`, has FoehnCraft read each `foehncraft.vtu` and `mesh.msh` back against KoflCraft's evidence, and fails if elliptic smoothing raises FoehnCraft's max non-orthogonality above the TFI grid's.
- Block mesher (MB1f): `block.grid` with `{"kind": "elliptic", "wall_sliding": true}` lets wall nodes slide along curved walls toward orthogonal grid lines, block vertices included. Sharp corners such as trailing edges stay fixed.
  - **Extruded sections:** nodes slide along the section contour. On the NACA 2412 C-grid, FoehnCraft's max non-orthogonality drops from 47.0° to 44.7°.
  - **Cylinder and torus walls:** the block-edge nodes slide on the surface and the face nodes follow. The pipe and the elbow stay at their quality.
- Block mesher: `block.template` with `multi_element_airfoil` writes an extruded H-grid around a main element with an optional slat and flap. Each element is a NACA 4-digit section or a Selig/Lednicer file, placed by chord, deflection, gap and overlap; trailing edges are blunt (0.5 % of the chord by default).
  - The layout needs a negative overlap (the slat's trailing edge ahead of the main element's nose, the flap's nose behind the main trailing edge). Other layouts are refused with a message that names the gap to change; a slot topology for overlapping elements is not built yet.
  - Elliptic smoothing on this template uses a faster control decay and up to 40,000 sweeps, and wall nodes whose spacing control sits at its bound (the blunt trailing-edge faces) no longer keep the smoother from converging.
- Multi-element airfoil template (MB1g): each element sits in a wall ring of six blocks, with a cove where a flap's nose sits under the main element, and a new layout optimiser (`block::layout`, after Sauer & Morsbach 2023) moves the corners, edge shapes and end spacings of the other blocks toward orthogonal cells. On nine test layouts the largest non-orthogonality falls from 76°–89° (up to 240 cells over 70°, two layouts refused) to 47°–68° with no cell over 70°; elliptic smoothing converges in 1,310 sweeps instead of 40,000. `block.template` reports the optimisation under `layout`. Design note: `research/block-layout-optimisation.md`.
- Triangle split for block topologies (MB1h): `block.split_triangle` turns a block with a flat corner into three blocks around a centre (the tri-center map). The new cells run as a row through the neighbouring blocks to the boundary. `block.optimize_layout` optimises any extruded topology and tries every corner of 150° or more as a split. It keeps a split only when the worst angle of the split blocks drops by 2° or more and the whole grid gets no worse in angle, cells over 70°, folds or cell growth. On a triangle beside a quadrilateral, the worst angle falls from 62° to 26°. The multi-element airfoil runs the same trigger, but its layout has no triangle blocks, so nothing changes there. New templates `triangle` (three hexahedra) and `tetrahedron` (four hexahedra) are building blocks for 3D templates. Design note: `research/block-triangle-split.md`.

### Hand-off to FoehnCraft and file formats

- Gmsh export writes MSH 4.1 ASCII by default, the version FoehnCraft reads (FC-01, OPEN-11). `gmsh-2.2` keeps the old version. Box, STL, channel, sweep and duct meshes share one element model, so both versions carry the same elements. The M1 acceptance runner reads 4.1 independently.
- Every mesh type writes `foehncraft.vtu`, a single VTU with polyhedral cells, boundary faces, region numbers and patch and zone names that FoehnCraft reads (FC-02), and `mesh-evidence.json` with the summary FoehnCraft must reproduce and the SHA-256 of each file (FC-05). Box and STL meshes export `foehncraft` and `openfoam` (`polyMesh`, FC-03). Patch names become TOML bare keys with a role (FC-04, `docs/foehncraft.md`). Cell pairs that share several faces get one merged face, as FoehnCraft requires.
- `cargo m1 --require-foehn` reads every acceptance mesh exported for FoehnCraft with `foehn mesh` and compares counts, volume, patch areas, closure and inverted faces with KoflCraft's summary (FC-06). CI builds FoehnCraft and runs it.
- Streamed output (CORE-06, FC-05): every file of the sweep, channel and duct bundles and of `export` (VTU, Gmsh, FoehnCraft VTU, OpenFOAM) now goes to disk as it is written, hashed with SHA-256 on a second thread on the way, and `mesh-evidence.json` and the `mesh.export` result take their digests from that stream. Files are byte for byte unchanged. A 2.5 million cell sweep peaks at 3.5 GB instead of 6.5 GB and runs in 29 s instead of 35 s. A failed file still removes the whole bundle. [Pitfall](../docs/pitfalls.md#whole-files-in-memory)
- Blocked VTKHDF hand-off (FC-02, FC-05): `--handoff vtkhdf` on `generate-sweep`, `generate-channel` and `generate-duct` (tool input `"handoff": "vtkhdf"`) writes `foehncraft.vtkhdf` in place of `foehncraft.vtu`, and `export ARTIFACT foehncraft-vtkhdf OUT` writes it for box and STL meshes. The file is layout 1 of the blocked VTKHDF agreed with MurmelCraft and FoehnCraft: VTKHDF 2.4, cells in Hilbert order cut into partitions of 65,536 with one ghost layer, global ids, block bounds, and one PolyData block per patch with its owner cells and role, so a viewer can read only the blocks it needs. KoflCraft writes the HDF5 itself (`mesh_io::hdf5`, from the HDF5 File Format Specification), streamed and hashed like every other file. The default stays `foehncraft.vtu`. [Format](../docs/foehncraft.md#blocked-vtkhdf)
- The blocked binary `foehncraft.vtkhdf` is now the default FoehnCraft file of the sweep, channel and duct bundles, as the maintainer chose once FoehnCraft could read it (FoehnCraft PR #252); `--handoff vtu` (tool input `"handoff": "vtu"`) writes `foehncraft.vtu` as before. The FoehnCraft round trip of the full CI run now reads the VTKHDF bundles and exports too. The block mesher still writes `foehncraft.vtu`. [Format](../docs/foehncraft.md#blocked-vtkhdf)
- Every file KoflCraft writes is hashed with BLAKE3 instead of SHA-256 (the maintainer's choice, 2026-11-24): `mesh-evidence.json` is format 1.1 with `file_hash: "blake3"` (1.0 meant SHA-256 and stays checkable in MurmelCraft), and the `files` of tool results carry the key `blake3` instead of `sha256`. Provenance digests of small inputs (`source_sha256`, `config_sha256`, artifact and journal digests) stay SHA-256.

### Agents: tool registry, MCP, journal and verification grades

- New crate `mesh-tools`: the tool registry (AG-12) with 13 dotted tool ids, input and output JSON Schemas, action kinds and the result envelope `koflcraft.tool-result`, and the registry of diagnostic codes (AG-18, `diag.codes`). New binary `kofl` with `kofl commands [--markdown]`, `kofl call ID --input FILE`, `kofl diag codes`, `kofl info`; `docs/tools.md` is generated from the registry. `mesh-cli` runs every command through the registry; error objects gain the registered `diagnostic` code. `generate-duct` now prints `{schema_version, output, report}` like the MCP tool did.
- Journal and replay (AG-13): `kofl --journal FILE` or `KOFL_JOURNAL` records every tool call in a `koflcraft.journal` file; `kofl run --script FILE [--keep-going]` (tool `journal.replay`) runs it again and reports each call whose status, exit code or result digest differs (TOOL005). `mesh.export` results now name the format and the SHA-256 of every written file.
- The MCP server runs through the tool registry (AG-12, AG-14): `kofl mcp --read-root DIR --write-root DIR [--journal FILE]` (and `mesh-mcp` with the same arguments) refuses to start without both roots, reads only below read roots and writes only below write roots. `tools/list` names every tool by its dotted id with annotations and an envelope output schema; `tools/call` returns the `koflcraft.tool-result` envelope as structured content and still accepts the old tool names. `mesh-mcp --output-dir DIR` is kept as short for one root of both kinds. Calls can be journaled with `via: mcp`.
- `cargo xtask layering` (SW-19) checks the crate layers of requirements section 6 and that the front ends reach the engine only through `mesh-tools`. CI gains a `checks` job: format, layering, clippy with `-D warnings` and the test suite. Requirements section 1a and the M1b status are updated.
- Verification grades (AG-17): `docs/features.toml` holds each feature's grade (verified, partial, unverified, open), how it was checked, its tests and gaps; `cargo xtask status` generates `docs/status.md`, and `cargo xtask status --check` (in CI and as a test) fails when the page is out of date, a cited test does not exist, or a requirement id is unknown. Requirements section 1a now points to the page. New dependency `toml` in `xtask` only.
- `mesh.report` version 2 (AG-16): over MCP it returns a summary (counts and extremes of the cell and face metrics) with one page of rows and a `next_cursor`; pages hold at most 500 rows (`page_size`). The command line keeps the full report. The acceptance runner follows the cursors and compares the reassembled rows with the CLI report.
- `mesh.estimate` and `kofl --check COMMAND` (AG-15): the exact cell count of every build tool's input before it runs, with rough memory and output bytes and the files it would write. The engine gains `channel::cell_count`, `sweep::cell_count` (2D section only) and `duct::cell_count` (from the surface, partition or scaffold); a test compares each estimate with the cells then built.
- MCP approvals (AG-15): before a build runs over MCP, the server estimates it. Builds within `max_build_cells` and `max_output_bytes` run; larger ones (or ones over an input limit) ask the user through MCP elicitation and run once on approval, with a limitation in the result. New codes TOOL009 (approval needed but the client cannot ask), TOOL010 (declined) and TOOL011 (refused by policy). `kofl mcp --policy FILE` reads a TOML policy (allow, ask or refuse per action kind, and the limits); the file must lie outside the write roots.

### Documents, colour maps and CI

- Multi-region meshing strategy decided on decision cards: `research/multi-region-strategy.md`, requirements v1.1 MD-20 to MD-27 (one scaffold for all polyhedral regions, region mesher menu and defaults, thin swept prisms, conformal block–poly joining, proposal table, hex-core track); MD-05 sliding after M5b.
- README slideshow: three pictures of the Manifold mesh (patches, cell types under a clip plane, non-orthogonality) made by MurmelCraft's headless `export.image`, an animated WebP near the top of the README with alt text, the stills with captions in [docs/gallery.md](../docs/gallery.md), and the recipe in [tools/readme-tour](../tools/readme-tour/README.md).
- Colour maps shared with FoehnCraft and MurmelCraft: `docs/colour-maps.md` sets Rainbow Desaturated (ParaView's eight control points) as the map of every picture, Inferno as the second, and holds the shared table that `tests/colour_maps.rs` checks.
- `docs/colour-maps.md` names Rainbow Uniform (Colin Ware) as the map of signed values centred on zero, as in MurmelCraft and FoehnCraft.
- **Faster CI**: the tests run with cargo-nextest, which runs the test binaries side by side instead of one after the other. Locally the fast suite went from 97 s to 48 s.

