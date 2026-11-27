# KoflCraft — 3D Polyhedral CFD and Multi-Physics Mesh Generator
## Software Requirements Specification — v1.2
**Date:** 2026-11-27
**Date:** 2026-11-23

**Purpose:** A Rust-native mesh generator for CFD and multi-physics simulation that can grow into a high-performance, agent-controllable polyhedral meshing system. It is the mesher of the Craft family and supplies meshes to the companion solver [FoehnCraft](https://github.com/chiefenne/Foehncraft).  
**Document status:** Decisions consolidated 2026-11-17, with the M1 export scope and inlet/outlet layer treatment updated 2026-11-18. **v0.7 (2026-11-22):** reworked to follow the development method of the [ArtCraft Craft Anything apps](https://getartcraft.com/apps/craft-anything) in the same way FoehnCraft does (its decisions D196–D203 and its [ArtCraft assessment](https://github.com/chiefenne/Foehncraft/blob/main/docs/assessments/artcraft-development-method.md)), and to define the interface to FoehnCraft. New sections: 1a (implementation status), 7 (rewritten), 7a (repository rules), 16 (FoehnCraft interface). OPEN-07 to OPEN-12 were decided by the maintainer the same day, each as proposed. **v0.8 (2026-11-22):** Leif Kobbelt (RWTH Aachen) added as third key person of the graphics principle (OBJ-19); online screening of the polyhedral mesher direction ([`research/graphics-screening-2026-11.md`](../research/graphics-screening-2026-11.md)) recorded in section 9a; new section 4a with requirements the maintainer named for mesh grading, thin shells, multi-domain meshing and agent-defined boundary faces, and BL-11/BL-12 for poly-prism layers; clean-room handling of the key groups' publications and software (LIC-12 to LIC-15) with a licence register of the screened software (section 14a). **v0.9 (2026-11-23):** the maintainer asked to check and complete the requirements for multi-domain meshing (heat exchangers with several fluid passes and several solid materials) and for STEP input. Multi-domain: MD-01 to MD-05 sharpened and MD-06 to MD-19 added (domain definition, kinds and material tags, separate fluid passes, interface detection, naming and imprinting, thin solids, interface sizing and layers, per-domain validation, export, `domain.*` tools, scale, periodic slices, the heat-exchanger fixture `SHELL_TUBE_HX.step` supplied by the maintainer, port caps), with FC-12 for the multi-domain hand-off. STEP: new section 4b (STP-01 to STP-13: schemas, units, assemblies, CAD faces and edges as patch and feature candidates, watertight tessellation, healing, reader safety, provenance, library route, material properties), STEP software added to the licence register, OPEN-13 (STEP library route, decided by the maintainer the same day as proposed: a pure-Rust reader crate with project-owned tessellation) and OPEN-14 (STEP timing, decided the same day as proposed: a first STEP slice moves into M5b if the heat-exchanger example is a STEP file). Items marked **O** remain proposals. **D** marks a commitment regardless of implementation status. **v1.0 (2026-11-24):** the maintainer asked for a concept of a **block-structured mesher** as a second mesher variant and decided its requirements on decision cards the same day: new section 4c (BS-01 to BS-19), OBJ-21, milestone MB1; concept and literature in [`research/block-structured-concept.md`](../research/block-structured-concept.md). **v1.1 (2026-11-25):** multi-region strategy decided by the maintainer on seven decision cards: one scaffold for all polyhedral regions, shared interface surface, order of construction, region mesher menu with defaults, thin swept prisms, block next to poly, proposal table, hex-core track (MD-20 to MD-27); MD-05 sliding moved to after M5b; BS-17 joining method; note [`research/multi-region-strategy.md`](../research/multi-region-strategy.md). **v1.2 (2026-11-27):** MD-25 splits the interface quads along their Delaunay diagonal instead of about the quad centre (maintainer, decision card after the sweep join, PR #84).

### Requirement-status convention

- **D — Decided:** committed to the scope; not a claim that implementation is complete.
- **L — Planned for later:** decided, but excluded from milestone M1.
- **P — Possible later:** identified option, not committed.
- **X — Not planned:** explicitly outside current scope.
- **O — Open:** pending a decision; the present proposal is stated alongside it.

**Status is distinct from milestone and implementation progress.** All requirements below have an ID, a status, and, where relevant, a milestone.

## 1. Project objectives and decisions

| ID | Status | Requirement | Timing / notes |
|---|---|---|---|
| OBJ-01 | D | Build a standalone 3D, cell-centered finite-volume CFD mesh generator. | Project goal |
| OBJ-02 | D | Support general polyhedral cells with variable numbers of polygonal faces. | Mesh core from M1 |
| OBJ-03 | D | Implement project-owned application and numerical code entirely in Rust. | From M1; no Python or C++ application/kernel code |
| OBJ-04 | D | Begin with a constrained, testable prototype rather than a production mesher. | M1 |
| OBJ-05 | D | Follow the development method of the ArtCraft Craft Anything apps, as FoehnCraft does: one Rust engine with thin front ends, every behaviour a registered tool, headless first, honest status, rules and lessons in the repository, clean room. Work in a *Study → Plan → Build → Check → Approve* loop with reproducible evidence. | Reworked 2026-11-22; details in sections 7 and 7a. The method is adopted, not the products' behaviour |
| OBJ-06 | D | Expose the same underlying operations via Rust library, CLI, and MCP server. | Minimal set from M1 |
| OBJ-07 | D | Assess suitable computer-graphics and computational-geometry research before committing to major meshing algorithms. | Ongoing; sharpened into the guiding principle OBJ-19 on 2026-11-22 |
| OBJ-08 | L | Support demanding industrial mesh scales of approximately 1–10 million cells. | Performance milestone, not M1 |
| OBJ-09 | L | Generate prismatic boundary-layer meshes with conformal polyhedral core transitions. | After basic bulk meshing |
| OBJ-10 | L | Provide GPU-accelerated execution for algorithms shown to benefit from it. | CPU baseline first; optional GPU backend |
| OBJ-11 | X | Depend on commercial CAD/meshing tools or commercial runtime licenses. | Not planned |
| OBJ-13 | D | Release project-authored software under the MIT license. | Root `LICENSE`; repository policy |
| OBJ-14 | D | Implement algorithms independently from published mathematics and algorithms; do not copy or translate third-party source code. | Clean-room engineering process |
| OBJ-15 | D | Keep a traceable, versioned research-to-requirement-to-implementation record. | Publications and design rationale, no source-code reuse |
| OBJ-16 | L | Investigate wall-normal refinement using Voronoi seed distributions *without* a conventional prism stack. | Alternative near-wall meshing, not a replacement for committed prism-layer support |
| OBJ-19 | D | **Guiding principle (maintainer, 2026-11-22): the algorithms of KoflCraft come from the computer graphics and geometry processing community.** The key people are **Keenan Crane (Carnegie Mellon University)**, **Bruno Lévy (INRIA)** and **Leif Kobbelt (RWTH Aachen)** and their groups (Kobbelt added 2026-11-22). Their research and the software they publish are the first place to look for every geometry or meshing topic, and each design note asks where their methods could improve or speed up the meshing task. The principle is a starting point, not an exclusive rule: a method from elsewhere (another graphics group or a classical CFD mesher) is chosen where it is better suited, after the measured comparison of ALG-16 that records why. | Section 9a lists the topics, the candidates and the recommended direction. Publications and software are used under the clean-room rules LIC-12 to LIC-15 |
| OBJ-12 | P | Offer an interactive 3D GUI. | Headless-first architecture. If built, it dispatches the same tool ids as the CLI and MCP (AG-12) and follows FoehnCraft's GUI rules (its L7: translation keys, keyboard use, accessibility tree, theme tokens) |
| OBJ-17 | D | Be the mesher of the Craft family: every KoflCraft volume mesh is directly usable by FoehnCraft without conversion by a third tool. | Section 16. FoehnCraft defers its own general mesher (its L8) and reads meshes from files |
| OBJ-20 | D | The polyhedral mesher follows the direction recorded in section 9a: a Voronoi and power-diagram core with per-cell clipping and a completeness test, seed-based grading, mirrored seed pairs for boundaries, thin shells and region interfaces, exact plane-based predicates, and extruded poly-prism layers. Changing direction needs a new screening note and a maintainer decision. | Screening 2026-11-22; it is a direction, each step still passes ALG-16 |
| OBJ-21 | D | A **block-structured variant** is the second mesher of KoflCraft (maintainer, 2026-11-24): multi-block hexahedral meshes of classic CFD topologies (channels, ducts, pipes, elbows, airfoils), with the block topology given by templates or by an agent through tools. It shares the engine, registry, validation and FoehnCraft exports with the polyhedral variant (OBJ-20). | Section 4c; concept in [`research/block-structured-concept.md`](../research/block-structured-concept.md); milestone MB1 |
| OBJ-18 | L | Support multi-physics meshes: several named cell zones (fluid, solid, porous) in one conformal mesh, with named region interfaces, for conjugate heat transfer and the region model of FoehnCraft (its D112, D113). The reference case is a heat exchanger with several fluid passes and several solid materials (maintainer, 2026-11-23). | After single-region boundary-layer meshing (M5). Detailed in section 4a, MD-01 to MD-19 |

## 1a. Implementation status (2026-11-22)

Grades follow AG-17: **verified** (checked against an exact or independent reference by a named test), **partial** (works on the named fixtures, known gaps), **unverified** (built, not yet checked), **open** (not built). How each was checked is stated. A partial result is never presented as success.

**How this was checked (2026-11-22, after nearest-site clipping, ALG-04):** `cargo test --release --locked --workspace --all-targets` passed on Linux x86_64 (Rust 1.97): 77 passed, 0 failed, 1 ignored (the dormant box-only OpenFOAM test, which passes when run on its own). `cargo m1 --require-foehn` passed in the container with `foehn` 0.7.0 built from FoehnCraft's main branch (29 exports read back); Gmsh was not installed there, so the Gmsh part ran in CI only. CI (`.github/workflows/m1.yml`) runs `cargo fmt --check`, `cargo xtask layering`, clippy with `-D warnings`, the test suite, and `cargo m1 --require-gmsh --require-foehn`.

The grade of every feature, how it was checked and its gaps are in the generated [status page](status.md), from [features.toml](features.toml) (AG-17). This section used to hold that table by hand; `cargo xtask status --check` now fails in CI when the page is out of date or cites a test that does not exist.

## 2. Milestone M1 — minimum working generator

**Goal:** Generate a valid, inspectable general-polyhedral mesh in a convex rectangular 3D domain using one Rust codebase, and operate it from CLI or MCP.

| ID | Status | Requirement | Measurable result |
|---|---|---|---|
| M1-01 | D | Define an axis-aligned 3D box as the fluid domain. | `xmin..zmax` from configuration |
| M1-02 | D | Construct a general polyhedral mesh from a small set of interior sites. | At least 16, 64, and 256 site tests |
| M1-03 | D | Use clipped Voronoi cells: start with the bounding convex box and clip each cell by pairwise perpendicular-bisector half-spaces. | **Decided:** first algorithm; use exhaustive site pairs only for small M1 tests |
| M1-04 | D | Represent variable-length faces and cell connectivity without one heap-allocated Rust object per entity in the canonical mesh. | Packed connectivity arrays |
| M1-05 | D | Compute face normals/areas, cell centroids, and positive volumes. | Analytical box tests and conservation checks |
| M1-06 | D | Validate closed cells, shared internal faces, outward orientation, no gaps, no overlaps, and region-volume conservation. | Zero topology defects on M1 tests |
| M1-07 | D | Export Gmsh MSH with a conforming tetrahedral subdivision and six named boundary groups. | Stable Gmsh 4.x compatibility and independently checked per-cell volume and topology |
| M1-08 | D | Export a visualization-friendly polyhedral mesh, preferably VTK legacy or VTU with true polyhedron cell topology. | Openable in independent viewer |
| M1-09 | D | Provide a Rust CLI and minimal MCP interface for create/check/export/report. | Same core engine, equivalent outputs |
| M1-10 | D | Produce machine-readable generation evidence, diagnostics, and benchmark metadata. | JSON reports, fixture IDs, seed and version |
| M1-11 | D | Repeatable generation for fixed inputs and deterministic site order. | Identical topology digest in repeat tests on same build/platform |
| M1-12 | X | STL/STEP input, arbitrary concave domains, boundary layers, 10M cells, or GPU speedups are M1 acceptance requirements. | Explicitly excluded from M1 |

### M1 exclusions are milestone exclusions, not project cancellations

- **L:** STL import; STEP import; complex surfaces; boundary-layer extrusion; GPU processing; million-cell scaling.
- **P:** GUI, direct non-Voronoi polyhedral meshing, non-convex cells.

### M1 acceptance criteria

1. Valid meshes of a cube using 16, 64, and 256 distinct interior sites, including nonuniform distributions.
2. Total cell volume agrees with the prescribed box volume within a documented floating-point tolerance.
3. Every interior face has exactly two opposite cell incidences; every boundary face has exactly one.
4. Zero negative-volume, zero unclosed or overlapping cells, zero unclassified external faces.
5. Basic finite-volume metrics are computed and reported for every face/cell.
6. VTU and Gmsh exports are independently read and validated. The CI acceptance command also runs Gmsh's mesh consistency check.
7. The equivalent operation can be invoked from CLI and from MCP without divergent engine logic.
8. The build and test suite run without a GPU and without a commercial dependency.

## 3. Algorithms and geometry

| ID | Status | Requirement | Current proposal |
|---|---|---|---|
| ALG-01 | D | Initial bulk-polyhedron construction strategy. | **Clipped Voronoi on convex box**, with robust face reconstruction and shared-face reconciliation. Retain alternative backend interface. |
| ALG-02 | D | Ensure geometrical validity and consistent global topology at every committed mesh stage. | Validate before export and after any modifying stage. |
| ALG-03 | D | Use double-precision CPU coordinates and explicit scale-aware tolerances. | `f64` in mesh core. |
| ALG-04 | L | Efficient neighbor search instead of all-pairs clipping. | Spatial index/candidate culling; essential before large meshes. Built 2026-11-22: k-d tree and security radius ([note](../research/alg-04-neighbour-clipping.md)) |
| ALG-05 | L | Variable mesh density controlled by a spatial size field. | Seed density and relaxation. Built 2026-11-22 for convex STL domains: a closed-form gradation-limited field and density-weighted CVT seeding ([note](../research/m2-size-control.md)) |
| ALG-06 | P | Weighted Voronoi/power diagrams. | Evaluate for size control and boundary alignment. |
| ALG-07 | P | Delaunay tetrahedralization followed by polyhedral agglomeration. | Alternative backend if preferable for complex boundaries. |
| ALG-08 | P | Octree/cut-cell / dual contouring bulk meshing. | Alternative for highly parallel construction. |
| ALG-09 | L | Reliable geometric predicates and ambiguity resolution. | Research exact/adaptive predicates; must preserve topology. |
| ALG-10 | P | Non-convex general polyhedra. | Allow only if geometry/FVM constraints justify them. |
| ALG-11 | L | CAD/STL geometry surface queries, patch identifiers, feature edges, and inside/outside classification. | Added progressively. |
| ALG-12 | D | Separate construction of candidate convex cells from global face reconciliation and finite-volume mesh certification. | Locally valid polyhedra alone do not establish conforming global topology. |
| ALG-13 | D | Accept no missing cell fragments or non-conservative topology arising from heuristic neighbor cutoffs. | Any accelerated candidate search must be verified against an exact/reliable geometric completeness criterion or CPU reference. The nearest-site clipper stops only on the security radius and is compared with the all-pairs reference (2026-11-22) |
| ALG-14 | P | Conforming mirror-seed/VoroCrust-style construction rather than surface clipping. | Research alternative for future geometry, **not** the M1 algorithm. |
| ALG-15 | P | Voronoi seeds with wall-normal anisotropic spacing instead of explicit prism stacks. | Compare in CFD validation for WMLES and wall-resolved requirements. |
| ALG-16 | D | **Graphics-first selection rule.** Each algorithm choice has a design note in `research/` that names the computer graphics / geometry processing candidate (OBJ-19), the classical alternative, and a benchmark on the project fixtures that compares them for validity first (closed cells, conservation, no inverted cells), then robustness on the hard fixtures, then speed and memory. The note records the outcome even when the classical method wins. | Same rule as FoehnCraft's 2026-11-21 decision on geometry topics (robust mesh Booleans, fast winding numbers, dual contouring, narrow-band distance grids, heat and signed heat methods). Extends OBJ-15 |

**Important algorithm caveat:** The clipped-Voronoi M1 algorithm assumes a convex box. Simply clipping Voronoi cells against arbitrary non-convex domains or triangle soups does not guarantee valid CFD cells; extension needs a robust boundary-conforming strategy.

## 4. Geometry and boundary layers

| ID | Status | Requirement | Timing / proposal |
|---|---|---|---|
| GEO-01 | D | Support a programmatic analytic box domain. | M1 |
| GEO-02 | L | Import ASCII and binary STL, check closure/orientation, preserve named surface patches. | M2 |
| GEO-03 | L | Import STEP/BREP geometry while retaining CAD face/topology identification. | Detailed in section 4b (STP-01 to STP-13); timing is OPEN-14 |
| GEO-04 | L | Handle CAD tessellation with configurable geometric tolerances. | Detailed in STP-05 and STP-06 |
| GEO-05 | P | Repair damaged, non-manifold or self-intersecting inputs automatically. | Require evidence of correct geometry |
| GEO-06 | L | Assign an explicit boundary role to named patches, including wall, inlet, and outlet. A cap on a closed input shell can carry an inlet or outlet role. | Patch roles are supplied independently of triangle winding and plane direction |
| BL-01 | L | Extrude wall-adjacent prism layers with user-defined first thickness, number of layers, and growth ratio. | Initial planar-wall prototype first |
| BL-02 | L | Join the last prism face layer to a conformal general-polyhedral core. | After planar-wall demonstration |
| BL-03 | L | Support curved walls and independently configured boundary patches. | Later |
| BL-04 | L | Detect layer self-intersections, pinch-off, opposite-wall collisions, and inverted cells. | Before complex geometry support |
| BL-05 | L | Provide strict vs adaptive behavior, local thickness reduction/termination, and achieved-coverage reporting. | Production BL milestone |
| BL-06 | P | Infer first-layer thickness from y+ and flow assumptions. | Convenience only |
| BL-07 | X | Guarantee full boundary-layer coverage on every arbitrary, degenerate input without compromises. | Must report infeasible cases |
| BL-08 | P | Offer near-wall anisotropic Voronoi refinement without explicit prism cells. | Evaluate especially for wall-modeled LES; keep true prism layers as committed future capability. |
| BL-09 | O | Define separate near-wall mesh quality and accuracy criteria for wall-resolved RANS/LES versus WMLES. | Proposal: begin with channel flow, measure wall-normal spacing, surface flux, and y+ targets at solver level. |
| BL-10 | L | At an inlet or outlet, grow prism layers from adjoining wall patches up to the opening plane. Keep the wall-rim layer nodes on that plane. Label the terminal prism faces and the remaining core opening with the port patch. | Each terminal face has one incident cell and the boundary faces tile the original opening. Exercise on a straight channel in M3 and on the Manifold fixture in M5. |
| BL-11 | L | Layer cells on polyhedral surfaces are **poly-prisms**: each polygonal wall face (any number of edges) is extruded along per-vertex marching directions into a prism with that polygon as top and bottom, so a layer stack is aligned with the wall and the flow along it. Triangle and quad prisms are special cases. The outer face of the last layer is the inner boundary the polyhedral core conforms to (BL-02). | Maintainer 2026-11-22: "mostly prismatic if derived from polys". FoehnCraft reads general polyhedra (its D23), so no splitting is needed for the solver |
| BL-12 | L | Marching directions, layer thickness limits in gaps and at concave corners come from smooth fields computed on the surface or in the volume (candidates: heat method, signed heat method, harmonic fields by Walk on Spheres), not from face normals alone, and the method used is recorded in the evidence. | Section 9a; LIC-08 note on the boundary-layer patent in the screening |

## 4a. Grading, thin shells, multi-domain meshes and boundary definition

Requirements named by the maintainer on 2026-11-22. Candidate algorithms for each group are in section 9a and the [screening note](../research/graphics-screening-2026-11.md).

### Mesh grading

| ID | Status | Requirement | Notes |
|---|---|---|---|
| GRD-01 | L | One scalar **size field** h(x) drives the mesh density of the core, combined from several sources by the minimum, with a global minimum and maximum size. | Implements ALG-05. The seeds follow it (density ∝ h⁻³) |
| GRD-02 | L | **Datum objects** in the domain set a target size within a distance and a growth rate away from it: point, line segment, plane, box, sphere, cylinder, cone, and an STL surface or closed body. | Datum objects are not part of the fluid boundary. Duct path built 2026-11-25: `volume_scaffold.build.grading` and `grading.add/list/remove` ([doc](nonconvex-duct.md#local-sizes-datum-objects-and-coordinate-regions)) |
| GRD-03 | L | **Coordinate regions** set a size inside a region given by coordinate bounds (axis-aligned box, slab, cylinder in a named coordinate system). | —. Duct path built 2026-11-25: `volume_scaffold.build.grading` and `grading.add/list/remove` ([doc](nonconvex-duct.md#local-sizes-datum-objects-and-coordinate-regions)) |
| GRD-04 | L | Geometry-derived sizes: surface curvature, proximity of opposite walls (gap resolution with a minimum cell count across) and feature edges. | — |
| GRD-05 | L | **Gradation limit**: neighbouring cells differ in size by at most a configured ratio; the size field is smoothed to satisfy it before seeding. | Candidates: gradient limiting, diffusion of the field (also by Walk on Spheres) |
| GRD-06 | P | Anisotropic size (a metric tensor instead of a scalar), for wakes and near-wall Voronoi refinement. | Candidate Lp-CVT; links to BL-08, ALG-15 |
| GRD-07 | O | Further grading methods ("other tbd" of the maintainer). **Proposal:** solution-based refinement from a FoehnCraft result file, decided after GRD-01 to GRD-05 work. | — |
| GRD-08 | D | Every grading source is a typed, named input of the mesh tool (AG-12), reported in the evidence with the achieved sizes per source. | Agents can add, list and remove sources by name |

### Thin shells (zero-thickness walls)

| ID | Status | Requirement | Notes |
|---|---|---|---|
| SHL-01 | L | A surface inside the fluid (a baffle, plate or fin of zero thickness) is imprinted into the volume mesh: no cell crosses it and the cell faces on it tile it. | — |
| SHL-02 | L | The faces on a shell are written as **two boundary faces**, one per side, with one owner cell each, so the solver treats the shell as a wall; the two sides are named patches (`<shell>` and `<shell>_back` unless named otherwise) and are paired in the evidence. | Coordinate the pairing with FoehnCraft (FC-11) |
| SHL-03 | L | A shell may have free edges inside the fluid and may touch the outer boundary or another shell; cells at free edges stay closed and valid. | — |
| SHL-04 | L | Prism layers grow on both sides of a shell (BL-11), terminating at free edges with the same rules as at ports (BL-10). | — |
| SHL-05 | P | A shell may be marked porous or as a "jump" interface instead of a wall. | Needs FoehnCraft support |

### Multi-domain meshes

The reference case (maintainer, 2026-11-23) is a **heat exchanger**: several fluid passes (for example a hot and a cold stream) that never touch each other, separated by solid parts of several materials (tubes, fins, plates, casing), meshed for conjugate heat transfer in FoehnCraft. FoehnCraft groups cell zones into named regions with a kind (its D112), assigns materials to regions by name (its D100), couples regions through interfaces (its D117) and can assemble one mesh from several files (its D113).

| ID | Status | Requirement | Notes |
|---|---|---|---|
| MD-01 | L | Several domains (fluid, solid, porous) meshed in one run, each a named cell zone. There is no fixed limit on the number of domains or materials. | Extends OBJ-18 |
| MD-02 | L | **Conformal** interfaces: the two domains share every interface face; the face is written once per side as a named interface patch or once as an internal face, as the export format requires. | M5b. Candidate: mirrored seed pairs across the interface (section 9a), clipping as fallback |
| MD-03 | L | Each domain has its own size field and layer settings; layers can grow on either side of an interface. Grading sources (GRD-02, GRD-03) can be limited to named domains. | — |
| MD-04 | L | **Non-conformal** interfaces: domains meshed independently, the interface written as two patches with the exact intersection of their faces (or the overlap data the solver needs), conservative to round-off. | Later than MD-02. FoehnCraft D113 already imprints non-matching parts of a multi-file mesh; which side computes the intersection is decided with FoehnCraft (FC-11) |
| MD-05 | L | Moving or sliding non-conformal interfaces (rotor and stator): a rotor region and a stator region whose interface is a surface of revolution about the rotor axis (cylinder, disc, cone and chains of them), written as the two boundaries FoehnCraft's sliding mesh reads (its D182 variant B, D186). | Maintainer 2026-11-25 (strategy card 5): after the heat exchanger (M5b). Fully general non-conformal interfaces stay MD-04 |
| MD-06 | L | **Domain definition.** A domain is a closed body of the input (an STL solid or closed shell, a STEP solid, STP-03) or the region enclosed by several surfaces, picked by a material point or by body name. A fluid domain can be the complement of solid bodies inside an enclosure (the fluid around tubes and fins inside a casing). Bodies may be nested and may touch; every point of the meshed volume belongs to exactly one domain, and overlapping bodies are a diagnostic, not a silent choice. | Material point as in classical CFD meshers; inside/outside by winding numbers (FC-09) |
| MD-07 | L | **Kind and material tag.** Each domain carries a kind (`fluid`, `solid`, `porous`, the FoehnCraft region kinds of its D112) and a material tag (a name such as `water`, `oil`, `copper`, `aluminium`). Both are written with the mesh (FC-12). KoflCraft never sets material properties; FoehnCraft assigns them to the region by name (its D100). | As BC-05: no second source of truth. Maintainer 2026-11-23: a material named in a STEP file is accepted as the default tag of its bodies (STP-13), and the domain tools (MD-15) set or override tags where the file has none or a wrong one (Fusion's STEP exporter does not export physical materials, so the heat-exchanger tags are set this way). No reader for proprietary CAD formats such as Fusion `.f3d` |
| MD-08 | L | **Separate fluid passes.** Fluid domains that do not share a face stay separate domains and zones even when their material is the same, each with its own inlets and outlets. Two fluid domains touch only where an interface between them is declared; an undeclared fluid–fluid contact (a leak through a missing tube wall, a gap below the contact tolerance) is reported as a diagnostic with its location and area. | Typical CAD error in heat-exchanger models |
| MD-09 | L | **Interface detection and naming.** Interfaces between touching bodies are found geometrically within a contact tolerance (gap and overlap, relative to the local size and stated in the evidence). Each side gets a patch with role `interface`, named `<a>_to_<b>` and `<b>_to_<a>` unless named otherwise, paired in the evidence with both domain names. Interface kinds: fluid–solid, solid–solid (contact between materials; FoehnCraft applies any contact resistance), fluid–porous, and fluid–fluid only where declared. | BC-04 role `interface`; names follow FC-04 |
| MD-10 | L | **Imprinting of touching bodies.** Where two bodies touch with surfaces that do not match (a tube through a hole in a fin, a plate on a casing, CAD contact faces), the common area is imprinted on both surfaces before meshing, so the conformal interface of MD-02 exists; the parts that do not touch stay walls of their own domain. Gaps and overlaps above the contact tolerance are reported, never closed silently. | Exact plane-based Booleans and predicates (section 9a, Kobbelt's group, Lévy 2025) |
| MD-11 | L | **Thin solid parts.** Thin solids (tube walls, plates, fins) get a configured minimum number of cells across their thickness, by the gap sizing of GRD-04 applied per domain. | Heat exchangers have many walls well under 1 mm |
| MD-12 | P | A solid thinner than a threshold may be represented as a thermal shell (SHL-01, SHL-05) with its thickness and material tag recorded, instead of being meshed. | Needs FoehnCraft support for a conducting shell |
| MD-13 | L | **Sizes and layers at interfaces.** The gradation limit (GRD-05) holds across interfaces, so the two sides of a conformal interface differ in size by at most the configured ratio. Prism layers (BL-11) grow on the fluid side of fluid–solid interfaces by default and optionally on the solid side; the base faces of a layer stack are the interface faces, and the cells of the other domain conform to them. | Conjugate heat transfer needs the wall-normal resolution on the fluid side |
| MD-14 | D | **Validation per domain.** Every cell belongs to exactly one domain; the cell volumes of each domain sum to the volume of its body; every face between two domains is an interface face; each conformal interface face has one cell on each side with opposite orientation; the areas of the two sides of an interface agree to round-off. Part of mesh validation (CORE-03). | Commitment for when multi-domain meshing is built |
| MD-15 | L | **Domain tools.** Domains are defined through registered tools (`domain.*`): list bodies with volume, bounding box and touching neighbours; propose domains; name a domain and set its kind and material tag; merge and split; list interfaces with their areas. Each step writes a new versioned `koflcraft.domain-map` file, journalled and replayable, as the patch tools do (BC-02, AG-12, AG-13). | Agent-driven definition, as for boundary patches |
| MD-16 | L | **Scale.** A heat exchanger may have tens of domains and hundreds of interface patches (one per tube and fin contact). Reports summarise per domain and per interface kind and page the detail (AG-16); zone and patch names stay unique and valid FoehnCraft names (FC-04). | — |
| MD-17 | P | **Periodic slices.** A repeated section (one tube row of a heat exchanger) can be meshed with periodic patch pairs whose faces match under the given translation or rotation. | Needs FoehnCraft periodic support |
| MD-18 | L | **Heat-exchanger acceptance fixture.** The maintainer's shell-and-tube heat exchanger in [`test-fixtures/heat-exchanger/`](../test-fixtures/heat-exchanger/README.md) (version 2, 2026-11-23, AP242, millimetres) is the M5b acceptance fixture, in two variants: `SHELL_TUBE_HX.step` with open nozzles (tests MD-19) and `SHELL_TUBE_HX_closed_fluid_domains.step` with four named cap faces (`water_inlet`, `water_outlet`, `air_inlet`, `air_outlet`). Both hold 13 named solid bodies in three parts (`TUBES`: 7 tubes; `HOUSING`: `housing`; `BAFFLES`: 2 tube sheets and 3 baffles) and **no fluid bodies**, so the fluid domains `tube_side` and `shell_side` come from material points or from the regions the caps close (MD-06). Acceptance: both fluids and every solid domain are their own zones, the fluids never share a face (MD-08), every contact listed in the fixture README is imprinted and conformal, named on both sides and passes MD-14, the four ports carry the cap names in both variants, the analytic tube volumes and cap areas in the README hold, and `foehn mesh` loads the result with one region per domain and the same counts (FC-06). Both files name `Kunststoff` (density 1290) as the material of every part, so that is the expected material read from the file (STP-13); the intended tags (maintainer, 2026-11-23) are set through the domain tools: `copper` for the tubes and baffles, `steel` for the housing and tube sheets. Tube sheets and baffles share one part, so solid domains are defined per body. | Replaces the synthetic tube-in-tube stand-in. Also the first STEP fixture (OPEN-14) |
| MD-19 | L | **Port caps.** A fluid opening that the CAD leaves open (a nozzle ending in an annular wall face) is closed by a planar cap spanned by its rim edges, so the fluid domain is closed; the cap becomes an inlet or outlet patch (GEO-06, BC-04). Caps are proposed from planar rim loops and accepted by name through the patch or domain tools. Caps that the CAD already provides as named surface bodies (as in the closed heat-exchanger file) are used as they are, oriented by their domain and not by their normal in the file. | The open heat-exchanger file needs it at all four nozzles |
| MD-20 | D | **One scaffold for all polyhedral regions.** All regions meshed with the polyhedral mesher are built together in one tetrahedral scaffold: walls, shells and interfaces are inner constrained surfaces, the flood that separates inside from outside labels every tetrahedron with its region (material point or body), and the median dual is cut into one cell per vertex and region. Points added by boundary recovery on an interface are therefore shared by both sides. | Maintainer 2026-11-25 (strategy card 1). Note: [`research/multi-region-strategy.md`](../research/multi-region-strategy.md) |
| MD-21 | D | **Shared interface surface.** The whole model is remeshed into one watertight surface; each interface is one triangulation used by both of its regions, so conformal interfaces (MD-02) hold by construction and need no tolerance weld. | Strategy note, section 2 |
| MD-22 | D | **Order of construction.** Thin-swept and block regions are built first and their interface faces frozen as constrained faces of the scaffold; prism layers (walls, both sides of shells, the fluid side of fluid–solid interfaces) come next; the polyhedral regions follow (MD-20). | Strategy note, section 2 |
| MD-23 | L | **Region mesher menu.** Each region has one mesher: *polyhedral with layers*, *polyhedral*, *thin, swept prisms*, *thin wall, not meshed* (MD-12) or *block-structured* (OBJ-21). **Defaults:** fluids polyhedral with layers; bulky solids polyhedral without layers; thin solids swept prisms; block only when the user or the agent picks a template; *thin wall* is never a default. | Maintainer 2026-11-25 (strategy card 2: rules as proposed) |
| MD-24 | L | **Thin swept prisms.** A solid thinner than about 4 cells of its target size with two opposite faces is swept from one face to the other with at least 3 cells across; the cells are poly-prisms whose ends are the dual polygons of the interface triangulation, so the regions on both sides conform without a cut. The thickness comes from the two-sided measure of GRD-04; a part that cannot be swept falls back to polyhedral with gap sizing (MD-11), and the proposal says why. | Maintainer 2026-11-25 (strategy card 3). As the thin mesher of STAR-CCM+ |
| MD-25 | L | **Block next to poly.** The interface faces of a block region are frozen and given to the scaffold, each quad split along its Delaunay diagonal (Lawson's angle test); each hex face on the interface is split into the median-dual pieces of its polyhedral neighbours (the hex becomes a polyhedron with split faces), so the mesh stays conformal without imprinting in the solver. | Maintainer 2026-11-25 (strategy card 4). Details BS-17. Quad split changed from the quad centre to the Delaunay diagonal by the maintainer on a card, 2026-11-27: a split about the centre gives the thin wall rows flat triangles that no tetrahedralization has as faces, and the scaffold's boundary recovery then splits without end ([sweep-fill.md](../research/sweep-fill.md), step 2) |
| MD-26 | L | **Proposal table.** A tool (`domain.propose_meshing`, name fixed when built) returns one row per region: name, kind, material tag, suggested mesher, layers, its interfaces and their connection, and the reason for the suggestion. The user or agent accepts all rows in one step or changes single rows; each change writes a new domain map (MD-15). The `domain.*`, `patch.*` and `block.*` tools stay available for experts. | Maintainer 2026-11-25 (strategy card 6). UX first: engineers, simple and guided |
| MD-27 | L | **Hex-core region mesher.** An octree of Cartesian hexahedra in the bulk of a large fluid region, polyhedra in the transition to the poly-prism layers, hanging-node faces as polygons (FoehnCraft D23, D57). It joins the menu of MD-23 as *hex core with layers* only after it beats the polyhedral core on the project fixtures (ALG-16). | Maintainer 2026-11-25 (strategy card 7): start now as its own track; Fluent Mosaic poly-hexcore as classical method, octrees and dual contouring (ALG-08) as graphics-first leads |

### Boundary faces for FoehnCraft boundary conditions

| ID | Status | Requirement | Notes |
|---|---|---|---|
| BC-01 | L | Boundary patches can be defined after import: by STL solid names, by feature-angle region growing, by plane or box selection, by normal direction, and by merging or splitting patches. | Extends GEO-06. Built 2026-11-22 as `koflcraft.patch-map` and `patch.*` ([note](../research/m2-patch-tools.md)) |
| BC-02 | L | Each selection operation is a registered tool (`patch.*`), so an agent can define the patches through MCP and every step is journalled and replayable (AG-12, AG-13). | Agent-driven definition, maintainer 2026-11-22. Built 2026-11-22: ten tools, each writes a new map |
| BC-03 | L | KoflCraft can **propose** a patch segmentation (feature edges, planar and cylindrical regions) with geometric evidence (area, normal, bounding box, neighbours) so an agent or a person names and accepts patches instead of picking triangles. | Candidate: variational shape approximation. Ports first (2026-11-23): `patch.detect_ports` proposes planar pipe-end caps with evidence and writes inlet, outlet and walls. Built 2026-11-25: `patch.propose` (feature-edge regions classified as plane, cylinder, sphere or freeform from the Gauss map, openings, box sides; [note](../research/bc-03-patch-proposals.md)) |
| BC-04 | D | Each patch carries a role (wall, inlet, outlet, symmetry, interface, shell) that maps to a FoehnCraft boundary-condition kind, and the mapping is written with the mesh (FC-04). | Patch maps carry explicit roles into STL meshes and their exports (2026-11-22); `shell` waits for SHL work |
| BC-05 | L | KoflCraft can write a FoehnCraft boundary-condition skeleton (patch names, roles, placeholder values) that FoehnCraft's own tools fill in; KoflCraft never sets physical values. | Avoids a second source of truth. Built 2026-11-25: `bc.skeleton` writes a FoehnCraft case (format 0.7) with types from roles and TODO comments for values |
| BC-06 | D | The surface patch assignment is checked: every surface face belongs to exactly one patch, and the check is part of mesh validation (CORE-03). | Every read of a patch map checks it (2026-11-22) |

## 4b. STEP input

Requested by the maintainer on 2026-11-23. STEP (ISO 10303-21 clear-text files) is the usual exchange format of CAD models such as heat exchangers, and its B-rep keeps the faces, edges and bodies that STL loses. These requirements detail GEO-03 and GEO-04. The library route is OPEN-13 and the timing OPEN-14. The first fixture is the maintainer's heat exchanger (MD-18): AP242 from Autodesk Fusion (version 1 was AP214), with named surface models for the caps in its closed variant, millimetres, 13 solids in an assembly of three parts, planes, cylinders, cones, tori and B-spline surfaces, lines, circles and B-spline curves.

| ID | Status | Requirement | Notes |
|---|---|---|---|
| STP-01 | L | Read STEP files of the application protocols AP203, AP214 and AP242: manifold solid B-reps, closed shells, and open shells (candidates for thin shells, SHL-01). Other content (wireframe, PMI, drawings) is skipped with a listed diagnostic, never silently. The minimum entity set is what the heat-exchanger fixture uses: planes, cylinders, cones, tori and B-spline surfaces with knots, lines, circles and B-spline curves, complex (multi-type) entity instances. | ISO 10303-21 file structure; the entity schemas come from the published EXPRESS schemas |
| STP-02 | L | **Units.** The length unit of the file (SI unit with prefix, or a conversion-based unit such as inch) is read and the geometry converted to the target unit of the configuration (proposal: metres). Source unit, target unit and factor are reported; a model whose size looks inconsistent with its declared unit is flagged. | STL has no units, which is one reason STEP is wanted |
| STP-03 | L | **Assemblies and bodies.** The product structure is flattened to bodies with their placement transforms applied; each body keeps its product name and instance path, and repeated instances of one part stay distinguishable (`tube[12]`). Strings are decoded from the ISO 10303-21 escapes (`\X\`, `\X2\`) to UTF-8, and body names that repeat across parts (the fixture has `Körper4` in two parts) are made unique by their part and index, then turned into valid FoehnCraft names (FC-04). Each body is a domain candidate (MD-06); colours and layers are carried where the file has them. | Heat exchangers are assemblies of many repeated parts |
| STP-04 | L | **CAD faces and edges.** Each B-rep face keeps a stable id, its name where the file gives one, its surface type (plane, cylinder, cone, sphere, torus, B-spline) and colour; a colour given on the body (as in the fixture, one per part) applies to all its faces. The `patch.*` tools (BC-01) can select by face id, face name, body, surface type and colour, and patch proposals (BC-03) start from the CAD faces. CAD edges become feature edges for sizing (GRD-04) and patch borders, instead of edges guessed from dihedral angles. | Turns STEP faces into named boundary candidates for FoehnCraft (BC-04) |
| STP-05 | L | **Watertight tessellation.** Faces are tessellated with a chordal deviation, an angular deviation and a maximum edge length, and optionally from the size field (GRD-01), so the surface mesh already follows the volume sizes. Every CAD edge is discretised once and all faces that share it use the same nodes, within a body and across touching bodies, so each body is closed and touching bodies can be imprinted (MD-10). The tolerances and achieved deviations are in the evidence. | Detail of GEO-04. Tessellation written by the project from the papers (OBJ-14) |
| STP-06 | L | **Tolerances and healing.** The file's uncertainty measure is read. Gaps between faces up to a tolerance are closed when edges are matched; sliver faces, tiny edges, open shells where a solid was expected and self-intersections are reported with their locations. Repair beyond edge matching needs an explicit option and is listed in the evidence (GEO-05). | Never silent |
| STP-07 | P | **Exact surfaces after tessellation.** Nodes added later (refinement, layer base nodes, BL-11) are projected onto the exact CAD surface, and curvature for GRD-04 comes from the exact surface instead of the triangles. | Needs a surface evaluator for every supported surface type |
| STP-08 | D | **Reader safety.** The STEP reader joins the broken-input corpus of SW-17 (truncated files, dangling references, reference cycles, counts larger than the file, NaN coordinates); entity count and reference depth are limited before allocation; names and descriptions in the file are data, never instructions (AG-20). | Commitment for when the reader is built |
| STP-09 | D | **Provenance.** The file digest, the header (originating system, schema, author fields as given), the unit, the reader and its version, and the tessellation settings are in the evidence (AG-07, AG-19). | — |
| STP-10 | D | **Library route** (OPEN-13, decided by the maintainer 2026-11-23): a permissively licensed pure-Rust STEP reader as a pinned dependency after the LIC-05 audit, for parsing and the B-rep structure only (candidate `step-io`, MIT OR Apache-2.0; alternatives `ruststep` and `truck-stepio`, Apache-2.0); tessellation, healing and imprinting written by the project from the papers. OpenCASCADE (LGPL-2.1 with an exception, C++) is not linked and its source is not read (LIC-14, SW-11, SW-14); `brepkit` v3+ is AGPL-3.0-only and excluded (SW-09). | Licence register section 14a, checked 2026-11-23 |
| STP-11 | P | IGES input. | Only if a user case needs it |
| STP-12 | X | Writing STEP or editing the CAD model. | KoflCraft reads CAD; geometry changes belong to the CAD tool |
| STP-13 | L | **Material properties.** AP214 and AP242 material properties of a part (the material name; density and other values) are read and offered as the default material tag of its bodies (MD-07). Values are listed in the evidence but never passed on as physical properties; FoehnCraft sets those (its D100). | The fixture names `Kunststoff`, density 1290 kg/m³, for all three parts |

## 4c. Block-structured mesher variant

Requested by the maintainer on 2026-11-24 and decided on nine decision cards the same day (concept note, section 9: [`research/block-structured-concept.md`](../research/block-structured-concept.md)). Classic CFD topologies come first, with the topology given by templates or explicit tools; blocks are point-matched; grids inside blocks are built by transfinite interpolation with elliptic smoothing and hyperbolic marching, with graphics methods benchmarked against them (ALG-16); output is unstructured hexahedra through the existing writers; the variant is standalone at first and is built in its own milestone MB1. MB1 was first placed after M5; on 2026-11-24 the maintainer asked to start it at once, so it runs as a parallel track next to M5. Its first step MB1a (scope chosen on a card the same day: channel and pipe O-grid end to end) is built; grades are in the [status page](status.md), feature `block-mb1a`. Design ideas from the structured grid framework of the maintainer's PyAero `v3.0` branch are taken in as design input (concept section 6a); the code is written fresh in Rust for 3D.

| ID | Status | Requirement | Notes |
|---|---|---|---|
| BS-01 | L | The block-structured variant meshes a whole domain as point-matched multi-block hexahedra from a **block topology**, as a second mesher behind the same engine, tool registry, validation and writers as the polyhedral variant. | OBJ-21. Code in `mesh-algorithms` (module `block`), data model in `mesh-core`; no new crate unless the module outgrows it |
| BS-02 | L | A block topology is a versioned JSON file `koflcraft.block-topology` (`format`, `format_version`, SW-20): vertices, edges (a curve and a node distribution each), block faces (a surface and a boundary patch or an interface each), and hexahedral blocks (eight vertices, cell counts ni, nj, nk, domain name). It is written and changed only through tools. | Like `koflcraft.patch-map` (BC-02) and `koflcraft.domain-map` (MD-15) |
| BS-03 | L | **Templates.** Parametric tools write a complete topology for standard cases: box and channel, pipe O-grid with H-core, elbow, T-junction, annulus, and C-grid and O-grid around an airfoil section (extruded or swept). A template's result is editable by the explicit tools (BS-04). | Card 2 |
| BS-04 | L | **Explicit tools** `block.*`: add, move and remove vertices and blocks, set edge curves, counts and distributions, associate with geometry, split a block, list, check. Each step writes a new versioned topology, is journalled and replayable (AG-12, AG-13), and declares its effects (AG-15). | Card 2. Agent-driven definition, as for patches and domains |
| BS-05 | D | **Point-matched connections.** Two blocks that share a face share every node of it. Cell counts propagate along each chain of opposite edges; an inconsistent chain is a diagnostic that names the chain and the blocks, before any grid is built. Shared faces are built once, so no tolerance welding is needed. | Card 3. Hanging nodes and non-matching faces: BS-16 |
| BS-06 | L | **Edge distributions:** uniform, geometric (first spacing and ratio), two-sided hyperbolic tangent (Vinokur), and a first-cell height taken from the layer settings of BL-01; the achieved first spacing and growth ratio are reported per edge. | Grading sources (GRD-02, GRD-03) may later set counts; decided when MB1 starts |
| BS-07 | L | **Geometry association.** Vertices, edges and block faces can be tied to points, curves and surfaces of the input (STL patches, STEP faces and edges, analytic curves); boundary grids are projected onto them and the deviation from the geometry is reported. | GEO-02, section 4b |
| BS-08 | D | **Wall-node invariant.** Grid lines leave a wall exactly at the prescribed wall nodes, never resampled; a test checks node-for-node equality for every volume algorithm. | From PyAero `v3.0` |
| BS-09 | L | **Volume algorithms**, a separate choice from the topology: transfinite interpolation (standard Boolean sum and a Hermite variant with boundary derivative terms), elliptic smoothing (Winslow / Thompson–Thames–Mastin with Steger–Sorenson control functions for wall spacing and orthogonality), and hyperbolic marching (Steger–Chan) for wall-bounded blocks. Harmonic and conformal maps and volumetric parametrisations of the key groups (OBJ-19) are benchmarked against them per ALG-16. | Card 4; hyperbolic marching from PyAero `v3.0` |
| BS-10 | L | **Orthogonal wall block.** N layers built along surface normals with first height and growth ratio (the inputs of BL-01); the offset is limited by the local radius of curvature, and the block is thinned there with a warning instead of folding; normals blend toward the bisector at sharp edges. | From PyAero `v3.0` |
| BS-11 | L | **Post-smoothers** (constrained Laplacian, elliptic, angle-based) as a separate choice, with quality metrics before and after in the evidence. | From PyAero `v3.0` |
| BS-12 | D | A setting that the chosen algorithm cannot honour (an orthogonality request to standard TFI) is refused with a diagnostic, never silently ignored. | AG-08 |
| BS-13 | D | **Validation.** The mesh checks of CORE-03 plus a positive Jacobian at every corner of every hexahedron; per block the evidence gives minimum Jacobian, skewness, wall orthogonality, aspect ratio and growth ratio. | — |
| BS-14 | D | **Export** as unstructured hexahedra through the VTU, MSH 4.1 and `polyMesh` writers (FC-01 to FC-03), with patches from block faces and one cell zone per domain; the i, j, k range of each block is kept in the evidence. | Card 5. FoehnCraft refuses structured CGNS (its D188) |
| BS-15 | P | Structured multi-block files (CGNS structured zones with 1-to-1 connectivity, Plot3D). | Card 5; for other solvers |
| BS-16 | L | Hanging nodes (2:1 across a block face, the coarse face written as a polygon) and non-matching block faces (written as paired interface patches as MD-04). | Card 3; after MB1 |
| BS-17 | L | Coupling with the polyhedral variant: a mesher choice per domain of a multi-domain mesh (MD-01, MD-23), and structured blocks joined to a polyhedral core conformally by splitting the hex faces on the interface (MD-25). | Card 6; after MB1. Joining method decided 2026-11-25 (strategy card 4) |
| BS-18 | P | **Automatic block decomposition** of general bodies. Leads, first place to look (OBJ-19): the 3D motorcycle complex (Brückler, Gupta, Campen 2022), integer-grid maps and HexEx (Kobbelt's group), frame-field hex-dominant meshing (Lévy's group), direction fields and boundary-first flattening (Crane's group); polycube maps as classical alternative. | Card 1; needs a screening note (LIC-12 to LIC-15) |
| BS-19 | L | **MB1 acceptance fixtures:** a channel and a pipe O-grid (volume against the analytic value), an elbow, and an extruded airfoil C-grid; each passes BS-13 and the FoehnCraft round trip (FC-06). | — |

## 5. Data model and quality

| ID | Status | Requirement | Details |
|---|---|---|---|
| CORE-01 | D | `f64` vertex coordinates and indexed, offset-based variable-length face/cell connectivity. | Rust `Vec`/slices; compact serialization. |
| CORE-02 | D | One canonical internal face, owner/neighbour relationship, named boundary patches. | Common mesh representation. |
| CORE-03 | D | Mesh validation independently callable without regeneration. | Shared Rust API/CLI/MCP. |
| CORE-04 | D | Positive volume, manifold connectivity, closed cell faces, consistent normals, total-volume conservation. | M1 hard checks. |
| CORE-05 | D | Compute non-orthogonality, skewness, face size, volume and adjacency metrics. | M1 initial metrics; expand later. |
| CORE-06 | L | Efficient large-mesh memory layout, parallel iteration, 64-bit-safe indexing/offsets. | Necessary for millions of cells. |
| CORE-07 | L | Local mesh optimization and quality-constrained topological changes. | After first mesher. |
| CORE-08 | X | Accept invalid negative-volume or overlapping cells to satisfy requested mesh count. | Always reject or flag failure. |

## 6. Rust software architecture

One Cargo workspace, with a domain engine reusable by all interfaces. The crates are in layers, lowest first; a crate may depend only on lower layers, and the front ends reach the engine only through the tool registry (SW-19, AG-12). Three same-layer uses are allowed and listed with their reason in `crates/xtask`: `mesh-quality` uses `mesh-geometry` (validation against the surface), `mesh-io` uses `mesh-algorithms` (the writers read the meshes the algorithms build), and `mesh-cli` uses `mesh-mcp` (`kofl mcp`). `exists` marks crates in the workspace on 2026-11-22.

```text
koflcraft/
  crates/
    L0 mesh-core/        # geometry and indexed cell/face data          (exists)
    L1 mesh-geometry/    # STL import and surface checks; later STEP     (exists)
    L1 mesh-quality/     # validation, geometric/FVM quality             (exists)
    L2 mesh-algorithms/  # Voronoi clipping, sweeps, duct layers, dual   (exists)
    L2 mesh-io/          # VTU / Gmsh / later OpenFOAM writers, readers  (exists)
    L3 mesh-engine/      # operations, plans, orchestration, evidence    (exists)
    L4 mesh-tools/       # tool registry, journal, approvals (AG-12..16) (exists)
    L5 mesh-cli/         # command line, through mesh-tools only         (exists)
    L5 mesh-mcp/         # MCP server, the same tool ids                 (exists)
    -- mesh-acceptance/  # M1 acceptance runner (test tool)              (exists)
    -- xtask/            # repository tasks: `cargo xtask layering` (SW-19) (exists)
    later: mesh-gpu/     # optional GPU kernels (SW-07)
  examples/  test-fixtures/  research/  docs/
  AGENTS.md  docs/pitfalls.md                                           (new, PROC-01, PROC-02)
```

Boundary-layer code lives in `mesh-algorithms` (`sweep`, `channel`, `duct`); a separate `mesh-boundary` crate is no longer planned unless the module outgrows it.

| ID | Status | Requirement | Implementation note |
|---|---|---|---|
| SW-01 | D | Rust-only project implementation for CPU meshing, geometry, I/O, and automation. | `cargo build`, `cargo test` |
| SW-02 | D | Stable, typed Rust library API independent of CLI/MCP. | Meshing operations implemented only once. |
| SW-03 | D | Configuration and reports are schema-versioned and machine-readable. | `serde` JSON initially. |
| SW-04 | D | No implicit shell execution by agent tools and no uncontrolled writes. | Explicit allowed paths and approval mechanism. |
| SW-05 | D | Test units, invariants, property-based/fuzz tests for geometry. | CI with reproducible seeds. |
| SW-06 | L | Parallel CPU execution. | `rayon` or equivalent after correctness baseline. |
| SW-07 | L | GPU backend with CPU fallback. | Feature- and capability-gated. |
| SW-08 | D | GPU kernels must be authored in Rust; no project-authored WGSL, CUDA C/C++, GLSL, or HLSL kernels. | Investigate `rust-gpu` Rust-to-SPIR-V; defer GPU acceleration if toolchain cannot meet the requirement. |
| SW-09 | D | Evaluate pure-Rust `brepkit` as a STEP research candidate only, at a much later stage. | **Not an approved dependency:** current brepkit v3+ is AGPL-3.0-only or commercial (`brepkit-io` 4.11.41 is AGPL-3.0-only on crates.io, 2026-11-23). Its earlier MIT/Apache releases require version/license audit; no dependency now. The STEP route is OPEN-13 (STP-10) |
| SW-10 | P | Interactive viewer and web frontend. | Not needed for M1. |
| SW-11 | X | A Python orchestration layer or C++ computational kernel. | Not planned. |
| SW-12 | D | Original source code only in the project repository, with a documented independent design process. | No copying, porting, translation, or derived adaptation of third-party implementation code. |
| SW-13 | O | Third-party *infrastructure* crates and exact license allowlist. | Proposal: independently developed geometry/meshing code; separately version-pinned MIT-licensed crates allowed only after provenance/license review. |
| SW-14 | D | No AGPL, GPL, or commercial runtime dependencies in deliverables. | Includes future CAD/GPU toolchain assessment. Copyleft programs (for example Gmsh) may run as separate processes in tests and fixture preparation only (LIC-09). |
| SW-15 | D | `unsafe_code = "forbid"` workspace-wide. If a later GPU or SIMD module needs `unsafe`, it is named in the workspace configuration and every block carries a `// SAFETY:` comment. | As FoehnCraft D202. No `unsafe` exists today, so the forbid costs nothing now |
| SW-16 | D | No panics in non-test code of the engine, the readers and the tools: `unwrap`, `expect` and `panic!` are denied crate by crate (new code first, old code when touched), so an agent always gets a diagnostic code and never a dead process. | As FoehnCraft D202. 47 `unwrap`/`panic!` sites in library sources on 2026-11-22 |
| SW-17 | D | Every reader (STL, patch map, partition, scaffold MSH, configuration JSON, later STEP) has a corpus of broken inputs that must each end in a diagnostic: truncated files, declared counts larger than the file, NaN coordinates, wrong winding, non-manifold edges. | As FoehnCraft D202. Allocation from a declared count is checked against the file size first ([pitfalls](pitfalls.md)) |
| SW-18 | D | Tests follow the keep rule: a test stays only if its invariant, a plausible change that breaks it, and why no other test catches it can be named. Assert relations (closure, volume partition, matching faces), not chosen constants or timings. A bug fix adds a test that fails without it; a crash fix adds its input to the corpus of SW-17. Never delete, weaken or skip a test to get green. | As ArtCraft `testing.md` and FoehnCraft D202. Geometric reference values (analytic volumes, port areas) are relations and stay |
| SW-19 | D | `cargo xtask layering` checks the crate layers of section 6 and that the front ends call the engine only through `mesh-tools`; it runs in CI. Repository tasks are Rust `xtask` subcommands. | As FoehnCraft D203. Python stays allowed for user-side scripts and fixture preparation (OPEN-08) |
| SW-20 | D | Machine output follows FoehnCraft's conventions: `--json` gives one JSON object on standard output, human text and errors go to standard error, no interactive prompts; every machine format carries `format` and `format_version` with only additive changes within a major version; dates come from `SOURCE_DATE_EPOCH` when set, so identical inputs give byte-identical files. | As FoehnCraft D79, D80, D84. Today: `schema_version` integers and JSON on stdout; extend, do not break |

## 7. Agent-first operation (ArtCraft method, aligned with FoehnCraft)

KoflCraft adopts the development method of the ArtCraft Craft Anything apps as FoehnCraft did on 2026-11-21: the practices are listed as A1–A21 in FoehnCraft's [ArtCraft assessment](https://github.com/chiefenne/Foehncraft/blob/main/docs/assessments/artcraft-development-method.md) and were decided there as D196–D203. KoflCraft takes the same decisions so that one person or agent can drive both programs the same way. It adopts the method, not the behaviour of any product. The working loop stays **Study → Plan → Build → Check → Approve**.

The table keeps AG-01 to AG-11 and adds AG-12 to AG-21. The FoehnCraft column names the matching FoehnCraft decision.

| ID | Status | Requirement | Example / FoehnCraft |
|---|---|---|---|
| AG-01 | D | Study operation records imported geometry/configuration and known limitations. | `inspect_domain` |
| AG-02 | D | Build operation runs mesh generation and creates a versioned artifact. | `generate_mesh` |
| AG-03 | D | Check operation runs hard validity checks and returns structured diagnostics. | `check_mesh` |
| AG-04 | D | Agent can request statistics and quality reports. | `mesh_report` |
| AG-05 | D | Export must be explicit and confined to allowed project/output paths. | `export_mesh` |
| AG-06 | D | Every tool has documented JSON input/output schema, stable errors, and a bounded effect. | Shared Rust engine |
| AG-07 | D | Evidence records include input digest, config digest, algorithm identifier, random seed, runtime, versions, metrics and validation result. | Report JSON |
| AG-08 | D | No silent success when validation fails. | Failed build cannot be tagged solver-ready |
| AG-09 | L | Agent-driven experimental parameter sweeps and mesh-quality comparisons. | After reproducible M1 |
| AG-10 | P | Full autonomous multi-step optimization and agent-selected algorithm switching. | Human approval for consequential operations |
| AG-11 | X | Unrestricted code execution through the MCP server by default. | Security boundary |
| AG-12 | D | **One registry of tools.** Every engine operation a front end can call is a registered tool with a stable dotted id, a contract version, a typed input and output with JSON Schema, and an action kind (`read`, `write-files`, `start-build`). `kofl commands` prints the catalogue as JSON, `kofl commands --markdown` writes `docs/tools.md`, `kofl call <id> --input file.json` runs any tool and prints a result envelope (`koflcraft.tool-result`). The CLI and the MCP server dispatch the same ids; nothing exists in only one front end. | FoehnCraft D196. Ids such as `domain.inspect`, `mesh.generate`, `mesh.generate_stl`, `mesh.generate_sweep`, `mesh.generate_channel`, `mesh.generate_duct`, `mesh.check`, `mesh.report`, `mesh.export`, `system.info`, `system.commands`, `diag.codes`. Replaces today's separate CLI and MCP dispatch (`generate-duct` vs `generate_duct`) |
| AG-13 | D | **Journal and replay.** Every top-level tool call of a session, by a person or an agent, can be appended with its typed input, status, exit code and duration to a journal (`--journal FILE` or `KOFL_JOURNAL`); `kofl run --script FILE` (tool `journal.replay`) replays it without a window and reports each call whose result digest differs from its record. | FoehnCraft D197. Mesh digests (AG-07) make replay checks exact |
| AG-14 | D | **Mandatory roots for agents.** `kofl mcp` refuses to start without `--read-root` and `--write-root`, and refuses every path outside them; every write creates a new entry, never overwrites. The plain CLI may run without roots. KoflCraft starts no processes and opens no network connections. | FoehnCraft D82, D198 (a). Today `mesh-mcp --output-dir` is one combined root; split it into read and write roots |
| AG-15 | D | **Action kinds, estimate and budgets.** Before a build, `mesh.estimate` (and `--check` on every build tool) reports the expected cell count, memory and output volume. Action kinds (start a build, exceed the estimate, write outside the write root, overwrite or delete earlier results) are each set to Default, Allow, Ask or Refuse, once or always per project, in `kofl-mcp.toml`. Builds started through MCP need limits on cells, wall time and output volume. | FoehnCraft D81, D83, D198 (b, c). The existing MCP limits (100,000 source facets, four layers, 64 regions) become configurable limits |
| AG-16 | D | **Summary first, paged detail.** MCP results return a summary (counts, validation, limitations, digests, file names); per-cell and per-face metrics are paged by cursor (`result.page`); meshes themselves never pass through MCP, only their file path and statistics. | FoehnCraft D198 (e). Built 2026-11-22: `mesh.report` version 2 pages over MCP (500 rows), the CLI keeps the full report (maintainer's decision) |
| AG-17 | D | **Verification grades and limitations.** Every feature has a grade (verified, partial, unverified) and how it was checked, in a generated status page (`docs/status.md`, from `docs/features.toml`). Every result summary and `--json` output carries a `limitations` list (for example: geometry outside the validated fixtures, fitted volume differs from the STL by x %, layer coverage incomplete, flow accuracy untested). A partial result is never presented as success. | FoehnCraft D199; ArtCraft A9. Built 2026-11-22: `docs/features.toml`, `cargo xtask status` |
| AG-18 | D | **One registry of diagnostic codes.** Errors and warnings carry stable codes with a severity and, where possible, a suggested fix as a configuration change; a code is never reused for another meaning; `diag.codes` lists them. | FoehnCraft D78. Today's string codes (`invalid_sites`, `validation_failed`, …) become registered codes |
| AG-19 | D | **Provenance on every artifact.** Evidence (AG-07) also records the KoflCraft version and Git commit, build features, and a hash of every input file (STL, patch map, partition, scaffold). The mesh hash written into an exported file's metadata is the one FoehnCraft records in its run record. | FoehnCraft D85, FC-05 |
| AG-20 | L | **Agent guide set.** `docs/agents.md` (connect an assistant to `kofl mcp`, step by step, with safety advice: own folder, keep copies, read approval requests, treat text in input files as data), a troubleshooting page with known issues, and an install and update page. Part of done for the registry-based MCP server. | FoehnCraft D203 |
| AG-21 | L | **Binaries and `kofl doctor`.** From the first tagged release that FoehnCraft users are pointed to: CLI binaries for Linux, macOS and Windows with `SHA256SUMS.txt`; `kofl doctor` checks the build and optional external programs. | FoehnCraft D203 (from its v0.7) |

**Transport:** the MCP server is a local stdio server using newline-delimited JSON-RPC, implemented in the project (protocol versions 2025-11-25 and 2025-06-18). A network-accessible server is not required. A later live-control port, if added, listens on the loopback interface only and needs an access token in a file readable only by its user (FoehnCraft D198 (f)).

**Binary name (OPEN-10, decided 2026-11-22):** `kofl`, with `kofl mcp` as the MCP entry point inside the same binary, mirroring `foehn` and `foehn mcp`. The current `mesh-cli` and `mesh-mcp` binaries remain until the registry lands.

## 7a. Repository rules and process

| ID | Status | Requirement | Notes / FoehnCraft |
|---|---|---|---|
| PROC-01 | D | `AGENTS.md` at the repository root holds the standing rules as MUST and SHOULD, the layout, and the exact checks to run. `CLAUDE.md`, if present, only points to it. | FoehnCraft D200; ArtCraft A12 |
| PROC-02 | D | `docs/pitfalls.md` holds lessons learned (what was tried, what it measured, what to do instead), read before work and extended in the pull request that learns a lesson. | FoehnCraft D200; ArtCraft A11 |
| PROC-03 | D | Plans and requirements live in the repository; only chat-setup rules stay outside. Requirement ids are never reused. | FoehnCraft D200; ArtCraft A10 |
| PROC-04 | D | All work is committed and pushed to a GitHub repository; CI runs on every pull request. | OPEN-09 decides the repository. On 2026-11-22 the workspace had one commit and no remote |
| PROC-05 | D | Pull requests: one issue per pull request; moving code and changing behaviour in separate pull requests; the template has a section "What ran and what did not, and why", and a check that did not run is never reported as passed. Branch names describe the task. | FoehnCraft D203; ArtCraft A17 |
| PROC-06 | D | Changelog as one fragment per change in `changelog.d/<PR number>-<short-name>.md`. | As FoehnCraft |
| PROC-07 | D | An agent asks before steps a new commit cannot undo (publishing, releasing, deleting a branch or tag). Merging policy follows the maintainer's FoehnCraft rule: the agent may squash-merge its own pull request when CI is green. | FoehnCraft F12 records the same difference from ArtCraft A18 |
| PROC-08 | D | Code rules: lints are fixed, not silenced; `cargo xtask` tasks run in CI; no secrets in commits. | ArtCraft A16 |

## 8. GPU and performance strategy

**Principle:** GPU work is an optimization path, not a precondition for a valid mesh.

| ID | Status | Requirement | Proposed implementation |
|---|---|---|---|
| GPU-01 | D | Preserve a correct CPU reference algorithm for every accelerated operation. | Common input/output contract |
| GPU-02 | L | Evaluate parallel site clipping, candidate-site queries, distance fields and quality calculations on GPU. | Independent kernels/batches |
| GPU-03 | L | Benchmark end-to-end time including host/device transfer and memory allocation. | CPU vs GPU on identical fixtures |
| GPU-04 | L | Detect device availability, fall back to CPU, and report path used. | Portable behavior |
| GPU-05 | O | Define quantitative GPU adoption threshold. | **Proposal:** at least 2× end-to-end speedup on a representative workload, no validity regression, acceptable memory overhead; adoption itself is **D: benchmarks required**. |
| GPU-06 | O | Select Rust-authored GPU compilation/runtime strategy. | **Proposal:** prototype `rust-gpu` + `wgpu`/Vulkan where supported; establish device availability, `f64` and geometric predicate reliability. |
| GPU-07 | P | Device-specific optimized backends beyond a portable GPU interface. | Only if benchmarks warrant |
| GPU-08 | X | Assume that all topological operations are efficiently parallelizable on GPU. | CPU retains irregular updates and robust fallbacks |
| GPU-09 | D | Introduce a GPU path only if representative measured end-to-end performance improves and validity is preserved. | Speedup must include transfers, allocations and topology reconciliation. |
| GPU-10 | D | All project-authored GPU kernels are Rust source. | Device intermediates/binaries emitted by toolchains are not project-authored shader source. |

**Numerics:** Use `f64` and robust predicates on CPU. GPU device/shader `f64` availability must be probed, not assumed. Candidate generation and coarse geometric screening may use reduced precision *only* if final topology decisions are recomputed robustly on CPU or equivalently guaranteed.

## 9. Research watchlist and decision evidence

Research is an explicit development workstream, not a fixed algorithm selection. Papers below are relevant leads, not proof of full CFD meshing capability. Under OBJ-19 the computer graphics and geometry processing literature, above all the work of Keenan Crane (CMU), Bruno Lévy (INRIA) and Leif Kobbelt (RWTH Aachen), is where algorithm candidates are looked for first; section 9a maps their work to project topics.

| Status | Paper / group | Why relevant | Project interpretation |
|---|---|---|---|
| D | Yan, Wang, Lévy & Liu (2010), *Efficient Computation of 3D Clipped Voronoi Diagram* | Nonconvex restricted Voronoi on a tetrahedralized domain | User-supplied `main.pdf`; future geometric domain restriction |
| D | Liu, Ma, Guo & Yan (2022), *Parallel Computation of 3D Clipped Voronoi Diagrams* | GPU cell/tet clipping and neighbor pairing; reports error modes | User-supplied `VolumeVoronoiGPU_RG.pdf`; do not accept heuristic volume losses |
| D | VoroCrust (2019/2020), *Voronoi Meshing Without Clipping* | Conforming cells via controlled surface sampling and seed placement | Alternative to clipped boundary cells, **not** M1 |
| D | Cadence Stitch / Numerical (2026) | Industrial clipped Voronoi and seed-controlled near-wall WMLES | Investigate non-prismatic wall resolution; vendor claims require independent validation |
| D | Bernardo Taveira et al. (2026), *Scalable GPU Construction of 3D Voronoi and Power Diagrams* | Independently clipped convex cells with accelerated candidate search | **High priority:** assess for M1 clipper scalability and later GPU |
| D | Julius Nehring-Wirxel et al., RWTH Aachen (2025), *Exact and Efficient Mesh-Kernel Generation* | Robust half-space geometry and exact plane predicates | **High priority:** reliability techniques, not complete mesh generator |
| D | Bruno Lévy (2025), *Exact Predicates, Exact Constructions and Combinatorics for Mesh CSG* | Robust mesh intersections and complex boundaries | **High priority:** future STL/boundary conformity |
| D | Wuheng Gao & Renjie Chen (2025), *Parallel 3D Delaunay Triangulation on the GPU* | GPU construction of Delaunay dual topology | Assess alternate later core backend |
| D | Nathan Tihon, Jonathan Lambrechts, Jean-François Remacle (2026), *Robust and Efficient 3D α-Shapes on GPU* | Exact predicates on massively parallel geometry operations | General techniques, indirect relevance |
| D | Nicole Feng & Keenan Crane (2024), *A Heat Method for Generalized Signed Distance* | Robust distance/inside-outside on imperfect surfaces | Possible future geometry processing |
| D | Hossein Baktash, Mark Gillespie, Keenan Crane (2025/2026 listing), *Subgrid Marching Tetrahedra* | Locally parallel surface reconstruction that handles thin features | Possible future STL reconstruction, not direct polyhedral volume mesh |
| D | *Anisotropic boundary layer surface mesh generation by improved advancing-front/local-reconnection approach* (2026) | Advancing layer and local reconnection ideas | Relevant to BL research, **surface** layering specifically |
| P | RWTH Aachen field-aligned and hexahedral meshing work | Structure/quality insights | Lower priority for unrestricted polyhedra |
| D | Ray, Sokolov, Lefebvre & Lévy (2018), *Meshless Voronoi on the GPU*; Basselin et al. (2021), *Restricted Power Diagrams on the GPU* | Per-cell clipping by k nearest seeds with a security-radius completeness test | **High priority:** replaces all-pairs clipping (ALG-04) and gives the completeness proof of ALG-13 |
| D | Lévy & Liu (2010), *Lp CVT*; Lévy (2022), *Partial optimal transport for a constant-volume Lagrangian mesh* | Anisotropic and volume-prescribed Voronoi / power diagrams | Grading (GRD-01, GRD-06) |
| D | Trettner, Nehring-Wirxel & Kobbelt (2022), *EMBER*; Nehring-Wirxel, Trettner & Kobbelt (2021), *Fast exact Booleans … octree-embedded BSPs* | Exact plane-based mesh Booleans | Shell imprinting (SHL-01), port termination, robust clipping |
| D | Sawhney & Crane (2020), *Monte Carlo geometry processing*; Sawhney et al. (2023), *Walk on Stars*; `nv-tlabs/wosx` (Apache-2.0) | Grid-free point-wise PDE solves on the raw boundary, GPU | Size-field smoothing and layer direction fields (GRD-05, BL-12); fields only, never topology |
| D | Botsch & Kobbelt (2004), isotropic remeshing; Botsch, Kobbelt, Pauly, Alliez & Lévy (2010), *Polygon Mesh Processing* | Surface preparation under layers | M4, M5 |
| P | Diazzi & Attene (2021), *Convex polyhedral meshing for robust solid modeling* | Exact convex decomposition of a nonconvex domain | Alternative to the Gmsh scaffold (OPEN-07); its code is LGPL and is not read |
| P | Bruno Lévy and collaborators' Geogram Voronoi/power diagram papers | Mature conceptual foundations | Reimplement selected ideas in Rust, subject to licenses |

## 9a. Computer graphics principle: key groups and topics

Decided by the maintainer on 2026-11-22 (OBJ-19), with Leif Kobbelt (RWTH Aachen) added the same day. The table maps topics KoflCraft needs to published work of the three key groups. It lists leads for the design notes of ALG-16 and commits to no algorithm; the 2026-11-22 [screening note](../research/graphics-screening-2026-11.md) gives the candidates per requirement, the expected speed-ups and the sources.

### Recommended direction (screening 2026-11-22, OBJ-20)

1. **Core:** clipped Voronoi and power diagrams, each cell clipped by its k nearest seeds and accepted only when a security radius proves completeness (Lévy's group, *Meshless Voronoi on the GPU*). This is the M1 clipper made scalable and parallel, not a change of course.
2. **Grading through seeds:** size field (GRD) → seed density → density-weighted CVT; semi-discrete optimal transport (Lévy) where cell volumes must match the size field; Lp-CVT for anisotropy.
3. **Mirrored seed pairs** across boundaries, thin shells and region interfaces put Voronoi faces exactly on those surfaces (SHL, MD-02); clipping stays the fallback. Patent check required (LIC-08).
4. **Exact plane-based predicates and Booleans** (Kobbelt's group; Lévy 2025) for clipping, shell imprinting, port termination and layer collisions.
5. **Poly-prism layers** extruded from the polygonal wall faces (BL-11), directions and limits from heat-method and harmonic fields (BL-12).
6. **Walk on Spheres** (Crane's group, `wosx`) for fields before a mesh exists: size-field smoothing, direction fields, inside/outside. Never for topology.
7. The Delaunay-tetrahedra dual (ALG-07) stays the classical alternative to benchmark against.

| Topic | Candidate from the key groups | Where it enters | Status |
|---|---|---|---|
| Clipped and restricted Voronoi, power diagrams, efficient neighbour search | Lévy and collaborators: Yan, Wang, Lévy & Liu (2010); the Voronoi and power diagram papers behind Lévy's Geogram library | M1 clipper (done, all-pairs, kept as reference), ALG-04 neighbour search (done 2026-11-22), ALG-06, nonconvex domains in M4 | Studied (2010 paper); the Geogram line is **P** |
| Robust geometric predicates and exact constructions | Lévy (2025), *Exact Predicates, Exact Constructions and Combinatorics for Mesh CSG* | ALG-09, OPEN-06, STL clipping and boundary conformity (M4) | **High priority** to study |
| Mesh Booleans and intersections on imperfect surfaces | Lévy (2025), as above | Port termination and layer collision handling (M5), imprinting for FoehnCraft's immersed-boundary plan (FC-08) | Lead |
| Signed distance, inside/outside, wall distance on imperfect surfaces | Feng & Crane (2024), *A Heat Method for Generalized Signed Distance*; Crane, Weischedel & Wardetzky (2013), the heat method for distance, as FoehnCraft cites it | Layer growth directions and local thickness limits in gaps (BL-04, BL-05), FC-08 | Lead; FoehnCraft's D108/L13 plan names the same method |
| Surface reconstruction that keeps thin features | Baktash, Gillespie & Crane, *Subgrid Marching Tetrahedra* | Preparing damaged STL input (GEO-05) | **P** |
| Surface parameterization and remeshing under layers | Corman & Crane (2025), rectangular surface parameterization, as FoehnCraft cites it; Botsch & Kobbelt (2004), isotropic remeshing | Surface mesh quality under prism layers (BL-03, BL-11, M5) | Lead |
| Scalable, parallel Voronoi cells with a completeness proof | Ray, Sokolov, Lefebvre & Lévy (2018); Basselin et al. (2021) | ALG-04, ALG-13, GPU-02 | **High priority** |
| Grading with prescribed cell sizes, anisotropy | Lévy (2015, 2022) semi-discrete and partial optimal transport; Lévy & Liu (2010) Lp-CVT | GRD-01, GRD-06 | Density-weighted CVT chosen for GRD-01 after the 2026-11-22 benchmark against DistMesh; optimal transport and Lp-CVT remain leads |
| Exact Booleans and plane-based geometry | Kobbelt's group: Nehring-Wirxel, Trettner & Kobbelt (2021); Trettner, Nehring-Wirxel & Kobbelt (2022), *EMBER*; Nehring-Wirxel et al. (2025) | ALG-09, SHL-01, MD-04, port termination | **High priority** with Lévy (2025) |
| Mesh-free field solves on the raw surface | Sawhney & Crane (2020); Sawhney et al. (2023); `nv-tlabs/wosx` | GRD-05, BL-12 | Lead |
| Patch proposals for agent-defined boundaries | Feature-angle region growing; variational shape approximation (Cohen-Steiner, Alliez & Desbrun 2004, outside the key groups) | BC-03 | Lead |

**Fit with the existing policy.**
- **Watchlist:** before the 2026-11-22 screening it contained five works of the first two groups (rows for Yan et al. 2010, Lévy 2025, Feng & Crane 2024, Baktash et al., and the Geogram line). OBJ-19 turns them from leads into the first place to look. Non-graphics entries (Cadence, Discreetize, boundary-layer surface mesh paper) stay as leads for CFD-specific needs the graphics literature does not cover, mainly layers and near-wall resolution.
- **Clean room and MIT (section 14):** algorithms are built from the papers and their mathematics (LIC-02, LIC-03). How publications and software of the key groups are used is fixed in LIC-12 to LIC-15; the licences of the software screened so far are in the register of section 14a. An author's reference implementation is not a licence to copy it.
- **Limits of the principle:** a graphics paper often assumes a clean surface and aims at visual quality. KoflCraft needs closed, conforming, conserving cells (ALG-12, ALG-13), so a graphics method is accepted after the validity benchmark of ALG-16, never on the paper's claims alone. The patent assessment of LIC-08 still applies.

### Sources and starting URLs

1. Artcraft Craft Anything: https://getartcraft.com/apps/craft-anything
2. Artcraft apps, agent-ready architecture: https://getartcraft.com/
3. Taveira et al., 2026: https://arxiv.org/abs/2605.06408
4. RWTH, Nehring-Wirxel et al., 2025: https://www.graphics.rwth-aachen.de/publication/03355/
5. Lévy, 2025: https://doi.org/10.1145/3744642
6. Gao & Chen, 2025: https://www.sciencedirect.com/science/article/pii/S0010448525000946
7. Tihon et al., 2026: https://epubs.siam.org/doi/abs/10.1137/1.9781611979138.7
8. Crane research: https://www.cs.cmu.edu/~kmcrane/
9. RWTH publications: https://www.graphics.rwth-aachen.de/publications/
10. Bruno Lévy: https://brunolevy.github.io/
11. Boundary-layer surface paper (2026): https://www.sciencedirect.com/science/article/pii/S0045793026000745
12. `rmcp` Rust MCP SDK: https://docs.rs/rmcp/latest/rmcp/
13. Pure-Rust STEP candidate: https://docs.rs/brepkit-io/latest/brepkit_io/
14. Alternative `truck-stepio` status (STEP input incomplete per docs): https://docs.rs/truck-stepio/latest/truck_stepio/
15. `wgpu` GPU abstraction: https://wgpu.rs/doc/wgpu/
16. Rust-authored GPU kernels: https://rust-gpu.github.io/
17. Yan et al. 2010, DOI: https://doi.org/10.1007/978-3-642-13411-1_18
18. Liu et al. 2022, DOI: https://doi.org/10.1109/TVCG.2020.3012288
19. VoroCrust: https://doi.org/10.1145/3337680
20. Discreetize: https://discreetize.ai/technology/
21. Numerical 2026: https://www.numerical.it/en/clipped-voronoi-meshing-industrial-les-fidelity-charles/
22. Cadence primary account: https://community.cadence.com/cadence_blogs_8/b/cfd/posts/clipped-voronoi-meshing-industrial-les
23. brepkit licensing: https://github.com/andymai/brepkit#license
24. Screening note 2026-11-22: ../research/graphics-screening-2026-11.md
25. Ray et al. 2018, *Meshless Voronoi on the GPU*: https://hal.archives-ouvertes.fr/hal-01927559
26. Basselin et al. 2021, *Restricted Power Diagrams on the GPU*: https://diglib7.eg.org/bitstream/handle/10.1111/cgf142610/v40i2pp001-012.pdf
27. Lévy & Liu 2010, *Lp CVT*: https://hal.archives-ouvertes.fr/inria-00600251
28. Lévy 2022, partial optimal transport: https://arxiv.org/abs/2106.03936
29. Trettner, Nehring-Wirxel & Kobbelt 2022, *EMBER*: https://www.graphics.rwth-aachen.de/publication/03339
30. Nehring-Wirxel, Trettner & Kobbelt 2021, octree-embedded BSPs: https://arxiv.org/abs/2103.02486
31. Walk on Stars: https://www.cs.cmu.edu/~kmcrane/Projects/WalkOnStars/
32. WoSX: https://github.com/nv-tlabs/wosx
33. Monte Carlo geometry processing resources: https://rohan-sawhney.github.io/mcgp-resources/
34. Leif Kobbelt, RWTH Visual Computing Institute: https://www.graphics.rwth-aachen.de/
35. VoroCrust patent US 10,776,540 B1: https://patents.google.com/patent/US10776540B1/en
36. `step-io` STEP reader (MIT OR Apache-2.0): https://github.com/elgar328/step-io
37. `ruststep` and `truck` (Apache-2.0): https://github.com/ricosjp/ruststep, https://github.com/ricosjp/truck
38. OpenCASCADE licence (LGPL-2.1 with the Open CASCADE exception): https://github.com/Open-Cascade-SAS/OCCT/blob/master/OCCT_LGPL_EXCEPTION.txt
39. FoehnCraft regions, materials and multi-file meshes (its D100, D112, D113, D117): https://github.com/chiefenne/Foehncraft/blob/main/docs/requirements.md

## 10. Planned delivery milestones

| Milestone | Status | Scope | Exit evidence | Progress (2026-11-22) |
|---|---|---|---|---|
| **M1** | D | Rust polyhedral core; box domain; small clipped-Voronoi mesh; QA; VTU/Gmsh; CLI/MCP | 16/64/256-site fixtures pass invariants | Done, verified |
| **M1b** | D | Craft alignment: tool registry, journal, roots, grades and limitations, diagnostic codes, `AGENTS.md`, pitfalls, layering check, repository and CI (sections 7, 7a) | `kofl commands` lists every tool; CLI and MCP give identical result envelopes; `cargo xtask layering` green | Mostly done 2026-11-22 (PRs #5 to #11): registry, journal and replay, roots, diagnostic codes, layering check, FoehnCraft exports and round trip. Then AG-17 status page (#12), AG-16 paged report (#13), AG-15 estimates (#14) and MCP approvals (#15). Open: wall-time limits and persistent approvals |
| **M2** | L | STL input; boundary patch labeling and the first `patch.*` tools (BC-01, BC-02); surface validation; basic mesh size control (GRD-01 to GRD-05); FoehnCraft hand-off FC-01 to FC-03 | Box/convex-polyhedron STL cases; every fixture loads in `foehn mesh` (FC-06) | STL import and convex meshing partial; size control partial (GRD-01 to GRD-05; density-weighted CVT, [note](../research/m2-size-control.md); curvature, gap and feature-edge sizes, [note](../research/grd-04-geometry-sizing.md)); hand-off partial (M1b); patch tools (BC-01, BC-02) verified 2026-11-22 ([note](../research/m2-patch-tools.md)): patch maps, ten `patch.*` tools, names and roles in every export, loaded by `foehn mesh` |
| **M3** | L | Simple planar-wall prism boundary layers, conformal core | [Straight channel benchmark and port-area checks](../research/open-port-layers.md#run-the-straight-channel-benchmark) | Prototype partial (channel, swept block and pipe) |
| **M4** | L | Robust STL/triangulated geometry, feature handling, clipped Voronoi in nonconvex domains | Feature and narrow-gap benchmarks, including the [Manifold pipe-and-box surface](../test-fixtures/manifold/README.md) | Prototype partial via partitioned dual; the Manifold scaffold is built in Rust (remesh plus Delaunay, [`research/m4-rust-scaffold.md`](../research/m4-rust-scaffold.md)), the elbow still uses its Gmsh scaffold (OPEN-07) |
| **M5** | L | Robust prism layer termination, collision handling and complex geometry; first FoehnCraft flow run on a KoflCraft mesh | Narrow channels, corners, curved walls, the Manifold port termination checks, and a FoehnCraft pipe-flow case on the Manifold mesh | Port termination prototype partial |
| **M5b** | L | Multi-region meshes for conjugate heat transfer (OBJ-18, FC-04, FC-12, MD-01 to MD-03, MD-06 to MD-11, MD-13 to MD-16, MD-18 to MD-24, MD-26, the first STEP slice of OPEN-14) and thin shells (SHL-01 to SHL-04) | The maintainer's shell-and-tube heat exchanger read from STEP (MD-18, STP-01 to STP-05, STP-08, STP-09, STP-13): both fluids and every solid domain their own zone, conformal named interfaces, capped ports, loaded by FoehnCraft with one region per domain | Open; fixture added 2026-11-23 |
| **MB1** | L | Block-structured variant, parallel track next to M5 since 2026-11-24 (first placed after M5) (section 4c, OBJ-21): topology file and `block.*` tools, templates, point-matched blocks, edge distributions, TFI, elliptic and hyperbolic volume algorithms, orthogonal wall block, smoothers, validation and export (BS-01 to BS-14, BS-19) | The fixtures of BS-19 pass validation and the FoehnCraft round trip | Concept and requirements 2026-11-24 |
| **M6** | L | Large meshes, CPU optimization, parallelization, GPU studies | 1–10 million cells, audited time and memory |
| **M7** | L | STEP import (section 4b, STP-01 to STP-09) with a licence-audited Rust reader (OPEN-13); its first slice is in M5b, since the heat-exchanger example is STEP (OPEN-14) | STEP vs STL geometry benchmarks; the same case from STEP and from STL gives the same patches and volumes |
| **M8** | P | GUI, adaptive meshing, advanced solver formats, general non-conformal interfaces (MD-04); sliding rotor–stator interfaces (MD-05) come after M5b | Separate future decision |

The milestone ordering after M1 may be revised based on benchmark evidence. **STEP is deferred to M7**, except the first slice that OPEN-14 moves into M5b when the heat-exchanger example is a STEP file; evaluating a library is not an agreement to depend upon it.

## 11. Resolved decisions and remaining open requirements

### 11.1 Decisions adopted 2026-11-17

| ID | Status | Decision | Consequence |
|---|---|---|---|
| DEC-01 | **D** | Use clipped Voronoi for the first volume generator. | Box-domain half-space clipping in M1. |
| DEC-02 | **D** | Require Rust-authored GPU kernels. | Evaluate Rust-to-GPU toolchains; no alternative shader-language source authored by project. |
| DEC-03 | **D** | Evaluate `brepkit`, but implement STEP only much later. | License and quality research, no early dependency. |
| DEC-04 | **D** | Use GPU acceleration only after measured benefit. | CPU is correctness reference; track transfer, memory and end-to-end throughput. |
| DEC-05 | **D** | Project software under MIT; independent implementation from publications only. | Original source, documentation of references, third-party license and patent due diligence. |

### 11.2 Open engineering details (not changes to the decisions)

| ID | Status | Current proposal |
|---|---|---|
| OPEN-01 | O | Use 2× end-to-end acceleration as a starting benchmark target; waive only with explicit evidence. |
| OPEN-02 | O | Evaluate `rust-gpu` for Rust-authored shader kernels; choose precise device backend after a technical spike. |
| OPEN-03 | O | Accept audited MIT-licensed Rust infrastructure crates, while writing geometry/meshing algorithms independently; prohibit foreign meshing algorithm code. |
| OPEN-04 | O | Keep conventional prismatic boundary layers as **L** and investigate seed-controlled Voronoi near-wall resolution as **P**; compare for intended turbulence models. |
| OPEN-05 | O | Name the project when the initial architecture and license files are settled. |
| OPEN-06 | O | Decide independently implemented robust predicates (`f64` filters + adaptive/exact fallbacks) before general STL clipping. |
| OPEN-07 | **D** (maintainer, 2026-11-22) | The nonconvex duct path needs a Gmsh volume scaffold. Decided as proposed: keep Gmsh as an external fixture tool (allowed by LIC-09) for now, and plan a Rust tetrahedralization or scaffold-free construction for M4 so a user needs no external mesher. Update 2026-11-23: the Manifold scaffold is built in Rust (`volume_scaffold.build`: wall remesh within a chord tolerance, Delaunay tetrahedralization with boundary recovery, refinement and sliver removal); the elbow keeps its Gmsh file until ports are remeshed. |
| OPEN-08 | **D** (maintainer, 2026-11-22) | Python in fixture preparation (`prepare_scaffold.py`, `check_sections.py`) versus OBJ-03. Decided as proposed, as FoehnCraft D203: OBJ-03 covers application and numerical code; Python stays allowed for fixture preparation and ParaView checks, and moves to Rust `xtask` only when next changed. |
| OPEN-09 | **D** (maintainer, 2026-11-22) | Repository: `chiefenne/KoflCraft` on GitHub, public, MIT, with the existing CI workflow, created by the maintainer; then commit the current workspace in reviewable steps. |
| OPEN-10 | **D** (maintainer, 2026-11-22) | Binary name: `kofl` with `kofl mcp`, mirroring `foehn` (section 7). |
| OPEN-11 | **D** (maintainer, 2026-11-22) | Gmsh format, as proposed: write MSH 4.1 as the default (FoehnCraft D33 reads only 4.1) and keep 2.2 as an option while the acceptance runner uses it. |
| OPEN-12 | **D** (maintainer, 2026-11-22) | Primary polyhedral hand-off to FoehnCraft: **single-file VTU first.** Rationale as proposed: VTU with cells and boundary faces in one file (FC-02) as the first format, because it keeps true polyhedra and is already written; OpenFOAM `polyMesh` (FC-03) second, since a writer exists for the box. Gmsh MSH stays the tetrahedral-subdivision format. |
| OPEN-05 | — | Closed 2026-11-22: the project is named KoflCraft (Cargo package `koflcraft`). |
| OPEN-13 | **D** (maintainer, 2026-11-23) | STEP library route (STP-10). Decided as proposed: pure Rust. Parse STEP with a permissively licensed reader crate (first candidate `step-io`, MIT OR Apache-2.0, experimental, created 2026-04; alternatives `ruststep` and `truck-stepio`, Apache-2.0, last released 2024-09), pinned after the LIC-05 audit and wrapped behind a project-owned B-rep interface so it can be replaced; write tessellation, healing and imprinting ourselves from the papers. Alternatives: a project-owned reader written from ISO 10303-21 and the EXPRESS schemas (no dependency, more work); OpenCASCADE as an external converter process (works soonest, but C++ and LGPL at the user's side, which SW-11 and SW-14 exclude today); no STEP, users convert to STL with named solids. |
| OPEN-14 | **D** (maintainer, 2026-11-23) | STEP timing. The heat-exchanger example arrived as STEP the same day, so the first slice is part of M5b. Decided as proposed: keep the full STEP import in M7, but move a first slice (STP-01 to STP-05, STP-08, STP-09: read bodies, units, faces and edges, tessellate watertight) into M5b if the heat-exchanger example (MD-18) arrives as STEP, because it is the natural source of its bodies and material names. |

## 12. New external research: Discreetize and industrial LES

### 12.1 Discreetize — conforming Voronoi *without* clipping

**Reference:** https://discreetize.ai/technology/

**Observed published claims:** Their feature-mirroring approach places seeds across boundaries to form conforming Voronoi faces without cutting cells. The site also describes automated feature detection/defeaturing, handling of narrow regions, anisotropic resolution, seed-based random mesh ensembles and local mesh edits.

**Project assessment:** A valuable **P** alternative for *later* nonconvex and feature-rich domain meshing. It is **not** proof that our simple clipped-Voronoi implementation will inherit these features. Discreetize's proprietary details are outside our source basis. Review the public academic paper *VoroCrust: Voronoi Meshing Without Clipping* (2019/2020), DOI 10.1145/3337680, for independently implementable published principles. The 2020 work specifies boundary-conforming seed placement, protection of sharp features and an input triangulated boundary; it is not a ready-made prismatic CFD boundary-layer generator.

### 12.2 Cadence Fidelity Stitch / CharLES — clipped Voronoi and no-prism near-wall meshing

**References:**
- Numerical, 20 Oct 2026: https://www.numerical.it/en/clipped-voronoi-meshing-industrial-les-fidelity-charles/
- Primary vendor article: https://community.cadence.com/cadence_blogs_8/b/cfd/posts/clipped-voronoi-meshing-industrial-les

**Observed published claims:** Stitch generates clipped Voronoi meshes from seed-point clouds and adjusts wall-normal resolution with point spacing instead of a conventional prism stack; the vendor describes applications to wall-modeled LES. The reported example of about 30 million cells in four minutes used 96 CPU cores: it is not evidence of GPU-accelerated **meshing** or the performance of our future implementation. Published solver GPU throughput is a separate metric.

**Project assessment:** **P** path for non-prismatic near-wall meshes, especially wall-modeled LES; must be assessed separately from explicit prismatic BL requirements. For wall-resolved RANS/LES, first-cell distance, extreme anisotropy, wall-normal continuity and numerical diffusion require dedicated validation. No third-party implementation details or source code are adopted.

## 13. Uploaded scientific papers — technical interpretation

### 13.1 Yan, Wang, Lévy & Liu (2010): *Efficient Computation of 3D Clipped Voronoi Diagram*

**User-supplied file:** `main.pdf`  
**DOI:** 10.1007/978-3-642-13411-1_18

The paper treats the nonconvex-volume problem by representing the domain as a tetrahedral decomposition, computing Voronoi regions and finding intersections with local tetrahedra; it also describes boundary-cell detection and centroidal Voronoi optimization applied to **tetrahedral** mesh generation.

**Implications:** Relevant mathematical principles for future nonconvex domain restriction. It assumes a prior tetrahedralization of the volume, which may be expensive and is not part of M1. It does not, by itself, specify a complete conformal FVM polyMesh exporter with layer control.

### 13.2 Liu, Ma, Guo & Yan (2022): *Parallel Computation of 3D Clipped Voronoi Diagrams*

**User-supplied file:** `VolumeVoronoiGPU_RG.pdf`  
**DOI:** 10.1109/TVCG.2020.3012288 (journal publication 2022)

The GPU algorithm pre-tetrahedralizes the domain, pairs candidate cells and tetrahedra, and clips Voronoi cells using four tetrahedron planes. The parallel independent intersections and volume/centroid reductions are valuable candidates for Rust-authored GPU experimentation.

**Critical limitation:** The paper reports nonzero total-volume errors and describes missing or erroneously included cell contributions when nearest-neighbor heuristics fail, particularly with nonuniform sites or sliver tetrahedra. Figure 9 shows aggregate errors up to 0.20% in the illustrated examples, and Figure 11 illustrates severe cell-specific error. This is **not acceptable evidence of a conservative conformal CFD mesh**. The M1 and later CFD mesh paths must preserve strict topological and volume invariants. Such heuristic approximations may be used for acceleration proposals, not final unchecked mesh construction.

### 13.3 Distinguish three algorithm families

| Family | Boundary approach | Initial project use | Research state |
|---|---|---|---|
| Convex-domain clipped Voronoi | Clip Voronoi cells against box/domain half-spaces | **D**, M1 | Implement and validate |
| Nonconvex clipped Voronoi | Intersect with a tetrahedralized or otherwise decomposed domain | **L**, post-M1 | Study 2010 and 2022 papers |
| Conforming Voronoi without clipping | Deliberately place mirrored/surface-conforming seeds | **P**, alternative | Study VoroCrust / Discreetize claims |

**Technical invariant:** a collection of local clipped convex polyhedra is not automatically a unique, consistent, globally conformal mesh with matching internal faces. Explicit face reconciliation, geometric completeness, shared vertex handling and mesh conservation checks are mandatory.

## 14. MIT and clean-room policy

| ID | Status | Requirement |
|---|---|---|
| LIC-01 | D | The repository's project-authored source code and documentation are published under MIT. |
| LIC-02 | D | Meshing/geometry algorithms are implemented independently using mathematical principles and published papers. |
| LIC-03 | D | Do not copy, paste, translate, port or adapt third-party implementation source code, including online examples, repositories, or source snippets embedded in papers. |
| LIC-04 | D | Maintain research citations, new original design notes, and commit/test traceability. |
| LIC-05 | D | Perform dependency, license and provenance checks before introducing crates or distributing releases. |
| LIC-06 | D | Do not introduce AGPL, GPL, closed-source, or paid-runtime dependencies. |
| LIC-07 | O | Define whether the dependency allowlist is strictly `MIT`/`MIT OR Apache-2.0` (with MIT selected), or allows other MIT-compatible permissive infrastructure crates. **Proposal:** require explicit MIT option and reviewed transitive dependencies. |
| LIC-08 | L | Conduct a separate freedom-to-operate/patent assessment of material algorithm choices before commercial distribution. MIT licensing alone cannot guarantee absence of patent constraints. |
| LIC-09 | D | Do not read copyleft (GPL, LGPL, AGPL) meshing or geometry source code at all. Permissively licensed source may be read and is never copied. Copyleft programs, such as Gmsh, may run as separate processes for tests and fixture preparation; their output is data with recorded provenance. | Stricter clean room as FoehnCraft D201 and ArtCraft `clean-room.md`; replaces "licence permitting" readings |
| LIC-10 | D | Every fixture and asset records its provenance (author, tool and version, licence) beside it. Geometry or numbers produced by proprietary meshers need a decision of their own before use as a reference. | ArtCraft A13; FoehnCraft D201 |
| LIC-11 | D | Commercial meshers and codes may be named where appropriate (a file format, a published reference, a comparison), without breaking LIC-03 or their licences; trademarks used as format names are marked. | FoehnCraft owner decision 2026-11-22 |
| LIC-12 | D | The graphics-community principle (OBJ-19) works through papers and mathematics. **Publications** of the key groups (papers, preprints, theses, talks, project pages, supplementary derivations) may be read freely, cited in `research/` and used as the source of an independent implementation; pseudo-code is re-expressed from the mathematics, and code listings inside a paper count as source code (LIC-03). **Software** of the key groups is third-party source like any other and follows LIC-03, LIC-09, LIC-13 and LIC-14. | Maintainer principle 2026-11-22, refined the same day; section 9a |
| LIC-13 | D | **Papers first, code last.** A design note starts from the papers. Permissively licensed source of a screened project is opened only when the note records why the papers are not enough (for example an unstated tolerance), what was read, and that nothing was copied. READMEs, licence files, documentation and published benchmark numbers may always be read. Published binaries may be run as separate processes to produce comparison data, with provenance (LIC-10). | Screening 2026-11-22 read READMEs and licences only |
| LIC-14 | D | **Weak copyleft counts as copyleft for reading.** Source under MPL-2.0, EPL, CDDL or LGPL (including LGPL-dual-licensed projects where the permissive option is not explicit for the file) is not read. A repository whose top-level licence is permissive but which bundles copyleft parts (for example Geogram's optional TetGen and Triangle, libigl's `copyleft` directory) may be read only in its permissive files, and only after its licence file has been checked. | Stricter than LIC-09; register in section 14a |
| LIC-15 | D | The licence of every screened or used third-party project is recorded in the register of section 14a with the date, the licence and how it was checked (licence file read, GitHub detection, crates.io record), before anyone opens its source or adds it as a dependency. Patents found during screening are recorded there too and feed LIC-08. | Same day as the screening |

## 14a. Licence register of screened software

Checked 2026-11-22. "File" means the repository's `LICENSE` file was read; "GitHub" means the GitHub licence detection; "crates.io" means the published crate record. No source code of any project below has been read. "Read source" is what LIC-09, LIC-13 and LIC-14 allow, not a plan to read it.

| Project | Group | Licence | Checked by | Read source? | Use for KoflCraft |
|---|---|---|---|---|---|
| `nv-tlabs/wosx` (Walk on Spheres Extensions) | Sawhney (NVIDIA), Crane line | Apache-2.0 (patent grant); deps FCPW MIT, nanoflann BSD, oneTBB Apache-2.0, pcg32 Apache-2.0 | File, README | Allowed after LIC-13 note | Field solves (GRD-05, BL-12) from the papers |
| `rohan-sawhney/zombie` (origin of WoSX) | Sawhney, Crane | MIT | File | Allowed after LIC-13 note | As WoSX |
| `rohan-sawhney/fcpw` (closest point, ray queries) | Sawhney | MIT | File | Allowed after LIC-13 note | Geometry queries (FC-09) |
| `nmwsharp/geometry-central` | Sharp, Crane group | MIT | File | Allowed after LIC-13 note | Heat method, intrinsic triangulations |
| `nzfeng/signed-heat-3d` | Feng, Crane | MIT | GitHub | Allowed after LIC-13 note | Signed heat method (BL-12) |
| `nmwsharp/nonmanifold-laplacian` | Sharp, Crane | MIT | GitHub | Allowed after LIC-13 note | Laplacian on dirty meshes |
| `nmwsharp/polyscope` (viewer) | Sharp | MIT | GitHub | Allowed | Viewing only, never a dependency of the engine |
| `rohan-sawhney/mcgp-resources` | Sawhney, Crane | CC-BY-4.0 (tutorial) | GitHub | Text may be read and cited | Background for WoS |
| `BrunoLevy/geogram` | Lévy | BSD-3-Clause; bundles LM6 (LGPL), optional TetGen (AGPL) and Triangle (non-commercial) | File, `doc/devkit/license.dox` | Permissive files only (LIC-14) | Voronoi, power diagrams, OT, exact predicates from the papers |
| `BrunoLevy/GraphiteThree` | Lévy | GPL-2.0 | GitHub | **No** | — |
| `geogram_predicates` crate (Rust port of Geogram PCK) | third party | LGPL-3.0 OR MIT | crates.io | Not read: a port of another code base, and LGPL-dual (LIC-14) | Possible benchmark oracle only |
| `meshless_voronoi` crate (`yuyttenhove/meshless_voro`) | third party, implements Ray et al. 2018 | MIT OR Apache-2.0 OR LGPL-3.0 | crates.io | Not read (LIC-14); its binary output may serve as comparison data | Benchmark oracle for ALG-04 |
| OpenMesh | Kobbelt group | BSD-3-Clause (project licence page) | Search result, page not reachable from CI proxy | Allowed after LIC-13 note, once the file is read | Surface data structures |
| OpenVolumeMesh | Kobbelt group | LGPL-3.0 with exceptions | File | **No** | — |
| `sandialabs/vorocrust-meshing` | Sandia (not a key group) | BSD-3-Clause; method patented (US 10,776,540 B1) | File | Allowed after LIC-13 note; patent check first (LIC-08) | Mirrored seeds (SHL, MD-02) |
| `MarcoAttene/VolumeMesher` | Diazzi & Attene | LGPL-3.0 | GitHub | **No** | Paper only |
| `wildmeshing/fTetWild` | Hu, Panozzo et al. | MPL-2.0 | GitHub | **No** (LIC-14) | Paper only |
| `libigl/libigl` | Jacobson, Panozzo et al. | MPL-2.0 with a GPL `copyleft` directory (GitHub reports GPL-3.0) | GitHub | **No** (LIC-14) | Paper only |
| `robust` crate | georust | MIT OR Apache-2.0 | crates.io | Allowed | Candidate dependency for Shewchuk predicates under LIC-05, LIC-07 |
| `wgpu` | gfx-rs | MIT OR Apache-2.0 (crate 30.0.2; GitHub detects Apache-2.0) | crates.io, GitHub | Allowed | GPU path (GPU-06) |
| `step-io` (`elgar328/step-io`) | third party | MIT OR Apache-2.0 (crate 0.3.1, 2026-11-18; README marks it experimental) | crates.io, `LICENSE-MIT`, README (2026-11-23) | Allowed after LIC-13 note | First candidate STEP reader (STP-10, OPEN-13), as a dependency |
| `ruststep` (`ricosjp/ruststep`) | RICOS | Apache-2.0 (crate 0.4.0, 2024-09) | crates.io (2026-11-23) | Allowed after LIC-13 note | Alternative STEP reader |
| `truck-stepio`, `truck-modeling` (`ricosjp/truck`) | RICOS | Apache-2.0 (crates 0.3.0 and 0.6.0, 2024-09) | crates.io, README (2026-11-23) | Allowed after LIC-13 note | Alternative STEP reader; STEP input reported incomplete (source 14) |
| `iso-10303` crate | third party | MIT (0.5.0, 2021) | crates.io (2026-11-23) | Allowed | Reader-code generator, unmaintained; not a candidate |
| OpenCASCADE (OCCT) | Open Cascade SAS | LGPL-2.1 with the Open CASCADE exception | `LICENSE_LGPL_21.txt`, `OCCT_LGPL_EXCEPTION.txt` (2026-11-23) | **No** (LIC-14) | At most an external converter in tests (LIC-09); not linked |
| `opencascade` crate (`bschwind/opencascade-rs`) | third party | LGPL-2.1 (0.3.0) | crates.io (2026-11-23) | **No** (LIC-14) | — |
| `brepkit-io` | `andymai/brepkit` | AGPL-3.0-only (4.11.41, 2026-11-23) | crates.io (2026-11-23) | **No** | Excluded (SW-09, SW-14) |

**Patents seen during screening:** US 10,776,540 B1 (conforming Voronoi mesh, VoroCrust line); US 12,094,059 (anisotropic polyhedra boundary layer adaptation, not read). Both go to the LIC-08 assessment before the related work is distributed commercially.

**Dependency caveat:** As of 2026-11-17, `brepkit` v3+ is AGPL-3.0-only or commercial; versions through 2.129.x are reported MIT OR Apache-2.0. This only establishes possible *license* eligibility for a pinned historical version; it does not demonstrate suitability, maintenance, numerical robustness or permission to copy its implementation. The project will not depend on it at M1.

## 15. Immediate M1 engineering sequence (completed; kept for traceability)

Next sequence after v0.7: (1) repository and first commits (OPEN-09, PROC-04); (2) `AGENTS.md` and `docs/pitfalls.md` (PROC-01, PROC-02); (3) `mesh-tools` registry with `kofl commands` and `kofl call`, CLI and MCP moved onto it (AG-12), with `cargo xtask layering` (SW-19); (4) roots, diagnostic codes, limitations and grades (AG-14, AG-17, AG-18); (5) MSH 4.1 and single-file VTU export with the FoehnCraft round-trip test (FC-01, FC-02, FC-06); then continue M2.

Original M1 sequence:

1. Define `Vec3`, `Plane`, convex polyhedron and immutable mesh connectivity contracts in Rust.
2. Independently implement convex polyhedron/half-space intersection using published geometry principles.
3. Generate clipped Voronoi cells from deterministic seeds inside an axis-aligned box.
4. Reconcile shared faces/vertices, construct owner/neighbour topology and classify six box boundaries.
5. Compute volumes, centroids, face areas, surface normals, nonorthogonality and skewness.
6. Validate conservation, closedness and topology with 16/64/256-site fixtures and difficult near-degenerate seeds.
7. Export true polyhedra to VTU and a conforming tetrahedral subdivision to Gmsh. Verify both with independent readers.
8. Provide `mesh-cli` plus safe `mesh-mcp` study/build/check/report/export tools from the shared engine.
9. Produce evidence JSON containing build fingerprint, cited design version, seed set, result digest, timing and all validator results.
10. Begin CPU profiling; defer GPU kernels and STEP until measured need and their later milestones.

**M1 exit gate:** No invalid cells, no missing or duplicated shared faces, full volume partition within numerical tolerance, deterministic tested artifacts, and validated VTU and Gmsh output. The project is a *small but correct* mesh generator, not a marketing-grade demonstration.

## 16. Interface to FoehnCraft

FoehnCraft is the companion general-purpose CFD and multi-physics solver (MIT, Rust). It solves on a cell-centred, face-based mesh of arbitrary polyhedra (its D23), reads Gmsh MSH 4.1, SU2, CGNS, OpenFOAM `polyMesh`, Fluent `.msh` and VTK/VTU behind one reader interface (its D33, D187, D188, D192), and defers a built-in general mesher (its L8). KoflCraft fills that gap. The two programs stay separate repositories and binaries; they meet only through files and tool calls, never through shared source.

| ID | Status | Requirement | Notes |
|---|---|---|---|
| FC-01 | D | Export Gmsh MSH 4.1 (ASCII, optionally binary) with physical groups as named boundary-face sets and named volume groups, so that `foehn mesh file.msh` reads it unchanged. | FoehnCraft D33 reads 4.1 only; today KoflCraft writes 2.2 (OPEN-11) |
| FC-02 | D | Export one VTU file per volume mesh holding the polyhedral cells (type 42) and prisms together with the boundary faces as surface cells, a region array (`boundary`/`zone`) and the string arrays `boundary_names` and `zone_names`, as the FoehnCraft VTK reader expects. The current multi-file output (`mesh.vtm` plus patch files) stays for viewing. | FoehnCraft D192; keeps true polyhedra, unlike MSH (OPEN-12) |
| FC-03 | L | Export OpenFOAM `polyMesh` (`points`, `faces`, `owner`, `neighbour`, `boundary`) for every mesh type, not only the box. | FoehnCraft D188 reads it; a box writer exists |
| FC-04 | D | Patch and zone naming: boundary names follow FoehnCraft's rules for boundary names, with the patch role (wall, inlet, outlet, symmetry, interface) recorded in the evidence and, where the format allows, in the file; cell zones carry region names for FoehnCraft regions (its D112). Region interfaces of multi-region meshes (OBJ-18) are conformal and named on both sides. | GEO-06 roles map to FoehnCraft boundary-condition kinds; the mapping table is part of the interface doc |
| FC-05 | D | The mesh digest of AG-07 and the input digests of AG-19 are written into the exported file's metadata where the format allows, and into a side-car `koflcraft.mesh-evidence` JSON otherwise, so a FoehnCraft run record can name the exact KoflCraft build and inputs of its mesh. | FoehnCraft D85 hashes the mesh file; this links that hash to the generator |
| FC-06 | D | Round-trip acceptance: every acceptance and benchmark fixture exported for FoehnCraft is loaded by `foehn mesh` (run as a separate process when it is installed, as Gmsh is today) and must report the same cell, face and boundary counts, every cell closed (FoehnCraft's Σ n·dS check) and the same total volume and patch areas to relative 1e-12. CI runs it when a FoehnCraft binary is available and reports it as not run otherwise. | Adds FoehnCraft to the independent readers of M1-07 and M1-08 |
| FC-07 | D | Tool-id and envelope conventions match FoehnCraft (`system.*`, `diag.codes`, `result.page`, `journal.replay`; `format`/`format_version`; diagnostic-code style), so one agent can drive both through MCP and FoehnCraft's study layer (its D62) can call `mesh.*` tools for parametric mesh variants. | AG-12 to AG-18 |
| FC-08 | L | Meshes for FoehnCraft's immersed-boundary plan: a thin body-fitted layer mesh grown from a surface, whose outer surface is handed to FoehnCraft for imprinting (its L13 plan). The open problems FoehnCraft lists for that mesher (column insertion at convex edges, merging at concave corners and in gaps, surface quality under the layers) are KoflCraft BL work items. | Coordinate with FoehnCraft's [L13 plan](https://github.com/chiefenne/Foehncraft/blob/main/docs/requirements/l13-immersed-boundary-plan.md) |
| FC-09 | P | A shared geometry-query layer (bounding volume hierarchy, closest point, winding numbers) published once as an MIT crate and used by both projects. | FoehnCraft D108 builds one; decide only when both need the same queries |
| FC-11 | L | Thin shells (SHL-02) and non-conformal interfaces (MD-04) are handed over in a form FoehnCraft reads: paired boundary patches for the two sides of a shell, and interface patches with the pairing recorded in the evidence and, where the format allows, in the file. | FoehnCraft's support for paired patches is to be confirmed with FoehnCraft before SHL work starts |
| FC-12 | L | **Multi-domain hand-off.** A multi-domain mesh (MD-01) is written either as one file with every domain a named cell zone (VTU `zone_names`, MSH named volume groups, `polyMesh` `cellZones`) or as one file per domain whose interface nodes match, for FoehnCraft's multi-part assembly (its D113) and OpenFOAM's multi-region layout (`constant/<region>/polyMesh`). Kind and material tag of each domain (MD-07) and the interface pairing (MD-09) are in the evidence and, where the format allows, in the file. | FoehnCraft regions and materials: its D100, D112, D113 |
| FC-10 | X | KoflCraft depends on FoehnCraft source crates, or FoehnCraft links KoflCraft into its solver binary. | Files and tool calls only |

