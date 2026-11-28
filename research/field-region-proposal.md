# Field-guided region proposal for block-dominant meshing

Date: 2026-11-27. Requirement: BS-18 (automatic block decomposition), with
BS-17 (blocks next to the polyhedral mesher). Tool: `regions.propose`.
Decided by the maintainer on a card (2026-11-27): block-dominant meshing,
blocks where templates fit and the rest polyhedral with imprint; this step
builds the fields and the region proposal in Rust, without meshing.

## Question

Which parts of a CFD domain can a block template fill (a swept duct, a
junction, an O-grid about a body, a box chamber), and which parts stay for
the polyhedral mesher? The answer has to come from the geometry and the
ports alone, with no user topology, and must be cheap enough to run before
every mesh.

## Method

All fields are linear finite elements on the KoflCraft tetrahedral
scaffold (`scaffold::build` after a uniform wall remesh at spacing h,
default about 150 000 tetrahedra): cotangent stiffness K and lumped mass M.
Solver: conjugate gradients with a Jacobi preconditioner and Dirichlet
rows eliminated (`fields::sparse`).

1. **Wall distance d** by the heat method (Crane, Weischedel & Wardetzky
   2013, *Geodesics in Heat: A New Approach to Computing Distance Based on
   Heat Flow*, ACM TOG 32(5)): (M + tK) u = 0 with u = 1 on wall nodes,
   X = -grad u / |grad u| per tetrahedron, K d = div X with d = 0 on walls.
   t = (m h)^2 with m = 3 and h the mean edge length. Inlets, outlets and
   symmetry planes are not walls.
2. **Flow potential p**: K p = 0, p = 0 on inlets, 1 on outlets, no flux on
   walls and symmetry planes: the streamwise coordinate of potential flow.
3. **Clash** c = |cos(grad d, grad p)| per tetrahedron. A block template
   needs the wall-distance shells (O-grid rings, wall layers) to cross the
   streamwise slices at right angles; c >= 0.5 is a clash. Where |grad d|
   < 0.3 (on the medial ridge, the core of an O-grid) c is 0; where p has no
   gradient (a dead pocket) c is 1.
4. **Skeleton**: p is cut into slabs about 2h thick (count = L / 2h, with
   the streamwise length L = sum V|grad p| / sum V|grad p|^2), and each slab
   into its face-connected pieces. Pieces below 5 % of an average slab are
   merged into their neighbour. The pieces with their face adjacency are
   the Reeb graph of p (Reeb 1946; on meshes Ni, Garland & Hart 2004, *Fair
   Morse Functions for Extracting the Topological Structure of a Surface
   Mesh*). Each piece is drawn at the centroid of its ridge faces, the
   interior faces where grad d flips (d grows towards the face from both
   sides), i.e. on the medial ridge of the heat distance. End points: pieces
   with nothing above or below; branch points: pieces with two or more above
   or below; loops: the cycle rank of the graph (bridges by Tarjan 1974).
5. **Regions**: neighbouring pieces of one category form a region.
   Junction (a branch piece); box (a piece with at least 4 times the median
   piece's section area: a chamber where the potential is almost flat);
   sweep (clash on less than 30 % of the piece); remainder (a clashing run).
   The section area of a piece is A = sum V |grad p| / dp over its slab
   step dp, exact for any section, because |grad p| integrates to dp along
   every streamline. (Up to 2026-11-28 the test compared piece volumes; a
   slab of a wide duct is also shorter along the flow, by 1/|grad p|, so
   its volume grows with A squared and the wide pipe of a reducer with an
   area ratio of 2.8 read as a box at 7.7 times the volume.) A run of box
   pieces at least four equivalent diameters long is a wide duct, its own
   sweep region. A clean run is a sweep when its level sets, one per slab,
   meet the wall in one closed loop, and an obstacle when a second closed
   loop shows a body inside (the furnace's tube row); a run on a loop of
   the graph with one loop per level set is a parallel passage, a sweep
   (up to 2026-11-28 every run on a graph loop was an obstacle, so the
   tubes between two headers were never swept). A connected cluster of clashing
   tetrahedra inside a sweep, larger than 0.5 % of the volume, becomes a
   remainder of its own (clash zone): a slab would average a local clash
   away.
6. **Second pass** (asked by the maintainer the same day): every box and
   remainder region is taken as its own mesh (faces to other regions are
   cuts: no flux, not walls) and classified again with the same wall
   distance and its **Fiedler vector** (the lowest non-constant mode of the
   Laplacian with no flux everywhere; Fiedler 1973; Lévy 2006,
   *Laplace-Beltrami Eigenfunctions: Towards an Algorithm that Understands
   Geometry*) in place of the potential, by inverse iteration with the
   constant deflated. Reason: re-solving the potential with the cut faces as
   new inlets and outlets gives back the old potential, because the cuts lie
   on its level sets (the harmonic function restricted to the region already
   satisfies the new boundary conditions), so it cannot find anything new;
   the Fiedler vector follows the region's own long axis instead of the flow.

## Results

Release build, one core. Python reference (voxels, round 2):
`/mnt/project-files/koflcraft-fields/` in the project files.

| Case | Tetrahedra | Time | Blockable (c < 0.5) | Python | Skeleton | Regions (share of volume) | Second pass |
|---|---|---|---|---|---|---|---|
| Elbow STL | 93 283 | 3 s | 99.5 % | 99.8 % | 2 ends, 0 branches | sweep 100 % | — |
| T (square duct, `test-fixtures/ducts.rs`) | 20 613 | <1 s | 97.5 % | 99.0 % (round T) | 3 ends, 1 branch at the junction | sweep 90 %, junction 10 % | — |
| Y (fork) | 29 335 | <1 s | 96.4 % | 99.0 % (round Y) | 3 ends, 1 branch at the split | sweep 89 %, junction 11 % | — |
| Manifold | 110 103 | 17 s | 91.8 % | 89.7 % | 2 ends, 0 branches | sweep 28 %, box 72 % (chamber) | 58 % of the box becomes sweep (42 % of the domain) |
| Seven-tube furnace (`make_furnace.py … 7`) | 153 023 | 12 s | 87.6 % | 85.6 % | 1 end, 1 branch, 1 loop (at the outlet) | sweep 86 %, junction 1 %, remainder 13 % (burner, behind the first tubes, outlet) | 9 % of the domain inside the remainder becomes sweep or box |

Benchmark of the distance (ALG-16, from the Python round 2 on the same
cases): the exact (plain) distance gives 2093 skeleton end points on the
furnace, the heat distance with m = 3 gives 62; on the T the plain
distance loses an arm. The heat method is kept. The Rust skeleton uses the
Reeb graph of the potential instead of thinning the ridge voxels, so its
end and branch counts are those of the flow topology (1 to 3 ends here)
and the ridge only places the nodes.

### Parallel tubes, reducers, the furnace chamber (2026-11-28)

| Case | Before | After |
|---|---|---|
| Heat exchanger tube side (7 tubes, two bonnets) | obstacle 25 %, sweep 0 % | 7 sweeps (25 %), 2 bonnets box, 2 junctions |
| Five tubes between two headers (Gmsh) | obstacle 25 %, sweep 15 % | sweep 40 % (tubes and both pipes) |
| Three tubes between headers (`ducts::header_tubes`) | 3 obstacles | 3 sweeps |
| Round reducer r 0.5 to 0.3 (`ducts::reducer`) | sweep 31 %, box 69 % | sweep 100 % |
| Square reducer 1 to 0.6 (Gmsh) | sweep 31 %, box 69 % | sweep 100 % |
| Furnace, 1 and 7 tubes | sweep 89 % and 86 % | obstacle 89 % and 86 % |
| Manifold | box 72 %, second pass 42 % of the domain | box 66 %, second pass 36 % |
| Elbow, T, fork | unchanged | unchanged |

## Limits and next steps

- Nothing is meshed. The next step is the cut: sweep regions bounded by
  potential level sets outside the wall layers, the templates of BS-03, and
  the remainder through the polyhedral mesher with imprint (MD-25).
- "Box" means a piece much wider than the median piece; when the chamber
  is the whole domain (the furnace) it is an obstacle (its level sets have
  the tubes inside), not a box.
- A reducer with an area ratio of 4 or more reads as a box unless its wide
  run is four diameters long.
- Second-pass labels on small clash zones are noisy (few tetrahedra per
  slab).
- No cross-section shape is chosen yet (O-grid for a round section, H for a
  rectangle); the region centroid and slab range are in the report.
- Jacobi CG is enough up to about 10^6 nodes (here 200 to 500 iterations
  per solve).

## Provenance

Clean room: built from the papers cited above; no source code of any
meshing or geometry project was read. The T and Y fixtures are unions of
unit cubes written in code (`test-fixtures/ducts.rs`); the confidential
cooling jacket was not used.
