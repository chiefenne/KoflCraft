# Multi-region meshing strategy

Status: decided by the maintainer on decision cards, 2026-11-25. This note
records how the regions of a model are meshed, how they connect, and which
mesh type is used where. It sharpens MD-01 to MD-19, SHL-01 to SHL-05 and
BS-17 in [`docs/requirements.md`](../docs/requirements.md) as
MD-20 to MD-27 and updates MD-05 and BS-17. Nothing here is built yet.

## 1. Starting point (2026-11-25)

- The polyhedral variant is the vertex-centred median dual of a tetrahedral
  scaffold with poly-prism layers (decision of 2026-11-24, as STAR-CCM+ and
  Fluent). The scaffold (`mesh-algorithms/src/scaffold.rs`) is built from one
  closed, tagged surface: Delaunay insertion, then boundary recovery by
  **Steiner points on the surface triangles**, then a flood from the bounding
  box that drops the outside, then Delaunay refinement.
- The block-structured variant (MB1) meshes whole models from a block
  topology and is standalone (BS-17 open).
- No multi-domain code exists. The heat exchanger of MD-18 is the reference.
- FoehnCraft (checked 2026-11-25 in its repository): conjugate heat transfer
  runs on one conformal mesh with cell zones grouped into regions (its D112,
  7d1 done). Imprinting of non-matching faces exists only for sliding rotors on
  surfaces of revolution (`imprint.rs`, its D182 variant B and D186).
  General static imprinting (its D58, stage 7d3) and baffles with shell
  conduction (its D189, 7t2) are planned and were not found in the code;
  one-dimensional thin walls on wall boundaries (7t1) are there. The
  multi-file assembly of its D113 was not found in the code.

## 2. The pipeline

1. **Regions.** Bodies become regions with a kind and a material tag through
   the `domain.*` tools (MD-06, MD-07, MD-15). Contacts between bodies are
   found within the contact tolerance and imprinted on both surfaces
   (MD-09, MD-10), so that every interface exists once.
2. **One watertight surface mesh.** Walls, ports, shells and interfaces are
   remeshed together. Each interface is one triangulation that both of its
   regions use, so the two sides are conformal by construction.
3. **A mesher per region**, proposed automatically and changeable per region
   (section 3).
4. **Order of construction.** Thin regions (swept prisms) and block regions
   are built first; their faces on interfaces are frozen and become
   constrained faces of the next step. Then **all polyhedral regions are built
   together in one tetrahedral scaffold**: walls, shells and interfaces are
   inner constrained surfaces, the flood that today separates inside from
   outside labels each tetrahedron with its region (from a material point or a
   body), and the median dual is cut into one cell per vertex and region.
   Prism layers grow from walls, from both sides of shells, and from the fluid
   side of fluid–solid interfaces before the scaffold of the core is built;
   the other side conforms to the layer base faces (MD-13).
5. **Checks and export.** The solver rules and MD-14 per region, then one file
   with one zone per region (FC-12).

**Why one scaffold (card 1).** Boundary recovery adds points on the surface.
With one scaffold per region, the two sides of an interface would get
different points and stop matching, and every added point would have to be
synchronised to the neighbour. With one scaffold the points are shared, and
the cut along the interface is exact. This is how Fluent's watertight
workflow builds a multi-region volume mesh in one step.

## 3. Region menu and defaults

The user sees these names; each is a mesher setting of one region.

| Mesher | Used for | Default when |
|---|---|---|
| Polyhedral with layers | fluids | every fluid region |
| Polyhedral | bulky solids | every solid region that is not thin |
| Thin, swept prisms | tube walls, plates, fins | a solid thinner than about 4 cells of its target size with two opposite faces |
| Thin wall, not meshed | very thin parts whose inner temperature field is not needed | never by default; offered as an alternative in the proposal |
| Block-structured | pipe, channel, annulus, elbow, airfoil | only when the user or the agent picks a template |

Decided (card 2, rules as proposed): no automatic block suggestion and no
default layers on the solid side; solid-side layers stay an option (MD-13).

**Thin solids (card 3): swept prisms, at least 3 cells across.** The
interface triangulation on one face of the thin part is swept to the opposite
face through the thickness; the cells are poly-prisms whose top and bottom
are the dual polygons of that triangulation, so the fluids on both sides
conform to them without a cut. This is the thin mesher of STAR-CCM+ and the
stacked layers of other CFD meshers. The thickness comes from the two-sided
measure of `sizing/thickness.rs` (Ma, Bae & Choi 2012). A part that cannot be
swept (no matching opposite face, branching) falls back to polyhedral with
gap sizing (MD-11), and the proposal says why.

## 4. How regions connect

| Connection | Where | Written as |
|---|---|---|
| Conformal (default) | every interface in one model | internal faces between zones, or a pair of interface patches with matching faces (MD-02, FC-12) |
| Conformal block–poly (card 4) | block region next to a polyhedral one | the hex face on the interface is split into the pieces of its polyhedral neighbours; the hex becomes a polyhedron with a split face |
| Non-conformal pair | rotor and stator (section 5); later independent parts | two patches paired in the evidence; FoehnCraft imprints them |
| Thin shell | zero-thickness baffle in a fluid | two boundary patches, `<shell>` and `<shell>_back` (SHL-02) |

**Block next to poly (card 4).** The block face is frozen first. Its quads are
given to the scaffold as constrained faces split into triangles about the quad
centre; the polyhedral cells on that side get the median-dual pieces of those
triangles, and each hex face is split into the same pieces. The pieces lie in
the plane of the quad (in the bilinear surface for a warped quad, with the
split made from the same points on both sides), so the two sides match
exactly and no tolerance weld is needed. FoehnCraft reads general polyhedra
and polygons with any number of nodes (its D23, D57).

## 5. Sliding interfaces (card 5)

After the heat exchanger (M5b): a rotor region and a stator region whose
interface is a surface of revolution about the rotor axis (cylinder, disc,
cone and chains of them), written as the two boundaries FoehnCraft's sliding
mesh reads (its D182 variant B, D186). Each side is meshed in the shared
scaffold with the interface surface as a split, so the two sides start
matched at angle zero. Fully general non-conformal interfaces (MD-04) wait for
FoehnCraft's general imprinting (its D58, 7d3).

## 6. Setup for the engineer (card 6)

One proposal table, accepted with one click. `domain.propose_meshing`
(name to be fixed when it is built) returns one row per region: name, kind,
material tag, suggested mesher, layers (count and first height), its
interfaces with their connection, and the reason for the suggestion
("thin: 0.6 mm wall, 2 cells at the target size"). The user or the agent
accepts all rows or changes single rows; each change writes a new version of
the domain map, journalled and replayable (MD-15, AG-12, AG-13). The
`domain.*`, `patch.*` and `block.*` tools stay available for experts. The
same table is the MCP result for an agent and a view in a later GUI.

## 7. Hex core for large fluid volumes (card 7)

Decided: start now, in its own track next to the heat-exchanger work. A
**poly-hexcore** region mesher (the default bulk type of Fluent's Mosaic
meshing since 2019): an octree of Cartesian hexahedra fills the bulk of a
large fluid region, polyhedra fill the transition to the poly-prism layers
at the walls, and each hex next to a smaller cell is a polyhedron with split
faces (hanging nodes are ordinary polygons for FoehnCraft, its D57). It joins
the region menu as "Hex core with layers" only after it beats the pure
polyhedral core on the project fixtures for validity, cell count, quality and
build time (ALG-16); the graphics-first leads (OBJ-19) are adaptive octrees
and dual contouring (ALG-08) and the Mosaic transition as the classical
method. Its own research note comes first.

## 8. Order of work

1. Multi-region scaffold: interfaces as inner surfaces, region flood, dual
   cut per region, MD-14 checks, zones in the export. First fixture: two
   boxes in contact, then a tube in a box (fluid, solid, fluid).
2. Proposal table and domain map (`domain.*`).
3. Thin swept prisms, on the tube walls of the heat exchanger.
4. Layers on the fluid side of interfaces.
5. Heat exchanger from STEP (MD-18), the M5b acceptance.
6. Block–poly coupling (BS-17), then sliding (section 5).

## 9. References

- Fluent watertight geometry workflow and its multi-region volume mesh;
  STAR-CCM+ thin mesher, prism layer mesher and polyhedral mesher: see
  the research notes on robust polyhedral meshing of 2026-11-24 (project
  files, `koflcraft/research/robust-polyhedral-meshing/commercial_meshers.md`)
  and the dual-mesh decision of 2026-11-24.
- Ma, Bae & Choi (2012), *3D medial axis point approximation using nearest
  neighbors and the normal field* (thickness measure, `sizing/thickness.rs`).
- Farrell, Piggott, Pain, Gorman & Wilson (2009), *Conservative interpolation
  between unstructured meshes via supermesh construction* (FoehnCraft's
  imprinting, D58).
- Murphy, Mount & Gable (2001); Cohen-Steiner, Colin de Verdière & Yvinec
  (2002): conforming boundary recovery used by the scaffold.
