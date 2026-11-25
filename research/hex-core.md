# Hex core with layers (poly-hexcore) for large fluid regions

Design note of ALG-16 for a fifth entry in the region mesher menu of the multi-region strategy: **"Hex core with layers"**, a Cartesian octree of cubes in the bulk of a large fluid region, a polyhedral transition and poly-prism layers at the walls, in the manner of Fluent's Mosaic poly-hexcore. Question: on the project fixtures, does it beat the median-dual polyhedral core ("Polyhedral with layers") on validity, cell count, maximum non-orthogonality and build time? It joins the menu only if it does. Requested by the maintainer on a card of the multi-region strategy thread (2026-11-25).

No source code was read. Everything below is built from the papers, from product documentation at the level of a user and from first principles (LIC-12 to LIC-15).

## Candidates

| | Method | Leads |
|---|---|---|
| Classical reference | **Mosaic poly-hexcore** (Fluent): octree hexes in the bulk, poly-prism layers on the walls, the gap between them filled by polyhedra that connect to the hexes face to face | Zore, Sasanapuri, Parkhi & Varghese (2019), *Ansys Mosaic Poly-Hexcore Mesh for High-Lift Aircraft Configuration*, 21st AeSI CFD Symposium; Fluent *hexcore* (Cartesian core with a tet and pyramid transition) as its predecessor |
| Classical, other route | **Trimmed (cut-cell) Cartesian** mesh with prism layers (STAR-CCM+ trimmer), **octree plus snapping** (snappyHexMesh) | Product documentation only; snappyHexMesh is GPL, so its code is not read (LIC-03) |
| Graphics, bulk (OBJ-19) | **Adaptive octrees with 2:1 balance**, linear (Morton-keyed) storage | Sundar, Sampath & Biros (2008), *Bottom-up construction and 2:1 balance refinement of linear octrees in parallel*; Morton keys as in our scaffold insertion order |
| Graphics, boundary (OBJ-19, ALG-08) | **Dual contouring and feature-sensitive extraction** on the octree: cut or snap the boundary cubes at the wall | Kobbelt, Botsch, Schwanecke & Seidel (2001), *Feature sensitive surface extraction from volume data*; Ju, Losasso, Schaefer & Warren (2002), *Dual contouring of Hermite data*; Schaefer & Warren (2004), *Dual marching cubes*; Trettner & Kobbelt (2020), *Fast and robust QEF minimization using probabilistic quadrics*; exact cut cells through Nehring-Wirxel, Trettner & Kobbelt (2021), *Fast exact Booleans for iterated CSG using octree-embedded BSPs* |
| Graphics, hex-dominant (Lévy) | Hex-dominant meshes with polyhedral fill of the gaps | Sokolov, Ray, Untereiner & Lévy (2016), *Hexahedral-dominant meshing* (ACM TOG); Gao, Shen, Panozzo et al. (2019), *Feature preserving octree-based hexahedral meshing* |
| Graphics, queries (Crane, nv-tlabs/wosx) | Inside test and wall distance for the octree: fast winding numbers, closest-point queries; walk-on-spheres as a grid-free distance and size field | Barill, Dickson, Schmidt, Levin & Jacobson (2018), *Fast winding numbers*; Sawhney & Crane (2020), *Monte Carlo geometry processing* |

**Why not cut cells first.** A cut-cell or dual-contouring boundary puts arbitrary polyhedra (cut cubes, small slivers) right at the wall, where the prism layers must start from a clean surface. STAR's trimmer handles this by extruding the layers first and trimming against the layer top, which needs exact polyhedral clipping of every boundary cube (the octree-embedded BSPs above). KoflCraft already has a layer-top surface and a robust polyhedral core between it and an inner surface (the scaffold and its median dual). The Mosaic route reuses both: only the bulk changes. Cut cells stay the lead for a later step (ALG-08) if the transition costs too many cells.

## Construction (prototype `volume_hex_core_v1`)

1. **Size field and octree.** The scaffold's own spacing (wall size plus `rate` times wall distance, capped by `max_size` and the grading field) times `size_factor` gives a target size t at each point. Cubes of a lattice of edge `size_factor · max_size` covering the domain are split into eight while their edge is larger than √2 · t at their centre, so the cube edge is the power of two nearest to t on a log scale (an earlier rule, split while larger than t, made the cubes 1.5 times finer than the dual cells and tripled their count). The octree is then **2:1 balanced** across faces, edges and corners, so every hanging node is the midpoint of one edge or one face.
2. **Hex region H.** A leaf joins H when it is inside the fluid and its centre is farther from the wall than the layer stack plus `margin` times its own edge (the transition then has room for about one tetrahedron of the local size). Leaves that would make the boundary of H non-manifold (two cubes meeting only along an edge or at a corner) are dropped until the boundary is a manifold surface.
3. **H shell into the scaffold.** Each boundary face of H (a square, plus the hanging midpoints on its edges) becomes two triangles, or a fan from its centre when it has a hanging midpoint. The shell, facing out of H, goes into the scaffold Delaunay with the wall surface and is recovered like the wall. Tetrahedra inside H are dropped before refinement, so the scaffold fills only the gap between the wall and H. Refinement keeps its points clear of the shell as it keeps them clear of the wall.
4. **Layers and transition.** The wall offset and layers are unchanged; the H shell is a fixed inner boundary without layers. The median dual fills the gap. Each dual cell ends on the shell with small triangles (vertex, edge midpoint, face centroid of a shell triangle).
5. **Join.** Every cube of H becomes one cell. A face to a cube of the same size is one square. A face to a cube of half the size is that small cube's square on the big cube's side, so a big cube has four faces there (a hanging-node face, accepted by FoehnCraft's general polyhedra, its D23). A face on the shell is replaced by the dual triangles that lie on it, shared with the dual cells, exact because they lie in the cube's plane (the coupling chosen on card 4 of the strategy). Finally every hanging point on an edge of a cube face (hanging nodes, shell vertices on cube edges) is inserted into that face, so all cells close edge to edge. Lattice coordinates are computed by one formula, so points on a cube edge match by exact coordinates.
6. **Checks.** The solver rules (one face per cell pair, closed cells, S·d > 0, positive volume and pyramids) run on the whole mesh, as for the median dual; volume and patch areas must match the wall surface without the shell.

Code: [`hexcore.rs`](../crates/mesh-algorithms/src/hexcore.rs) (octree, balance, shell), [`duct/hex_join.rs`](../crates/mesh-algorithms/src/duct/hex_join.rs) (join), the hole in [`scaffold.rs`](../crates/mesh-algorithms/src/scaffold.rs) and the pinned shell in [`duct/volume_dual.rs`](../crates/mesh-algorithms/src/duct/volume_dual.rs). The join runs before the "split again until clean" loop, so dual cells next to the cubes that fail the solver rules are split like any other (on the furnace about 120 at first, 0 after six rounds). Partition setting: `"core": "hex_core"` in `volume_scaffold`, optional `"hex_core": {"size_factor": 1.0, "margin": 1.0}`; built scaffold only.

## Benchmark

Furnace fixture (one tube, 13,660 facets, 3 layers, the fixture's `partition.json`) with the interior `max_size` lowered to 0.0625 so that the bulk carries weight; release build, 4 cores, Linux x86_64. Both meshes read back by FoehnCraft `foehn mesh` 1.0.1 (main of 2026-11-25). Bulk cell size = median cbrt(volume) of the core cells in the box x 1..3, y 0.5..1.5, z 0.5..1.5, away from every wall and the tube.

| Core | Bulk size | Cells | Core cells | Faces | Inverted | Max / mean non-orth. | Faces > 70° | Skewness > 0.9 | Build |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Median dual v2 | 0.0683 | 408,237 | 72,699 | 1,459,107 | 0 | 87.7° / 10.6° | 657 | 1,024 | 41.7 s |
| Hex core, size factor 1.09, margin 1 | 0.0681 | 421,737 | 86,199 (44,575 cubes) | 1,363,815 | 0 | 88.6° / 11.0° | 702 | 1,029 | 26.1 s |
| Hex core, size factor 1.09, margin 0.5 | 0.0681 | 422,049 | 86,511 (46,573 cubes) | 1,356,594 | 0 | 89.9° / 11.0° | 878 | 1,028 | 23.2 s |
| Hex core, size factor 1, margin 1 | 0.0625 | 435,994 | 100,456 (59,444 cubes) | 1,405,815 | 0 | 86.1° / 10.7° | 768 | 1,028 | 24.4 s |
| Hex core, size factor 1.3, margin 0.5 | 0.0813 | 399,023 | 63,485 (25,673 cubes) | — | 0 | 88.7° / 11.0° | 825 | — | 23.9 s |

At the fixture's own `max_size` 0.25 the median dual has 3,917 core cells (bulk 0.269) and the hex core 5,193 (380 cubes of 0.25 to 0.0625, bulk 0.25), 4.2 s against 3.5 s; FoehnCraft reads both as valid.

**Findings.**
- **Validity:** every hex-core mesh passes the solver rules and reads back valid in FoehnCraft (no inverted faces, no negative volume, closure to 1e-16).
- **Cells:** at equal bulk size the hex core has 3% more cells. A cube holds as much volume as a median-dual cell of the same size, so the bulk gains no cells; the poly transition between the layers and the cubes adds about 41,600 dual cells, most of them on the shell, which is triangulated at the cube size.
- **Faces:** 7% fewer. A cube has 6 faces (more where a neighbour is smaller); a bulk dual cell of the median dual has many more neighbours (not counted here).
- **Non-orthogonality:** the worst faces are at the wall layers in both cores, so the maximum does not change; the faces between two cubes of one size are at 0°, and about 25° across a 2:1 size jump (geometry, not measured). The mean is 0.4° higher because of the transition cells.
- **Build time:** 38% faster at 0.0625, because the scaffold does not refine the bulk; the octree itself takes 0.06 s to 0.3 s.

**Outcome (maintainer, card of 2026-11-25):** the hex core wins on build time and faces but not on cells, so it is merged as an opt-in expert setting (`"core": "hex_core"`) and does **not** join the region menu yet. Next lead to make it compete: a cheaper transition, either a coarser shell triangulation (cube faces as single polygons, fewer shell vertices) or a cut-cell boundary on the layer top (ALG-08, STAR trimmer style) in place of the dual gap.

**Limits of the prototype.** One fluid region; median dual only (no Voronoi core, no `layers_yield`); built scaffold only; the octree does not follow the grading field's anisotropy; the inside test is ray parity (a ray along +x, nudged off the lattice); cubes at a non-manifold edge or corner of the shell are dropped (43 to 80 on the furnace) rather than added.
