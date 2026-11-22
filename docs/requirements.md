# KoflCraft — 3D Polyhedral CFD and Multi-Physics Mesh Generator
## Software Requirements Specification — v0.7

**Date:** 2026-11-22

**Purpose:** A Rust-native mesh generator for CFD and multi-physics simulation that can grow into a high-performance, agent-controllable polyhedral meshing system. It is the mesher of the Craft family and supplies meshes to the companion solver [FoehnCraft](https://github.com/chiefenne/Foehncraft).  
**Document status:** Decisions consolidated 2026-11-17, with the M1 export scope and inlet/outlet layer treatment updated 2026-11-18. **v0.7 (2026-11-22):** reworked to follow the development method of the [ArtCraft Craft Anything apps](https://getartcraft.com/apps/craft-anything) in the same way FoehnCraft does (its decisions D196–D203 and its [ArtCraft assessment](https://github.com/chiefenne/Foehncraft/blob/main/docs/assessments/artcraft-development-method.md)), and to define the interface to FoehnCraft. New sections: 1a (implementation status), 7 (rewritten), 7a (repository rules), 16 (FoehnCraft interface). OPEN-07 to OPEN-12 were decided by the maintainer the same day, each as proposed. Items marked **O** remain proposals. **D** marks a commitment regardless of implementation status.

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
| OBJ-19 | D | **Guiding principle (maintainer, 2026-11-22): the algorithms of KoflCraft come from the computer graphics and geometry processing community.** The key people are **Keenan Crane (Carnegie Mellon University)** and **Bruno Lévy (INRIA)** and their groups. For every geometry or meshing topic, their published methods and the papers that build on them are the first candidates; a classical CFD-mesher method replaces them only where a measured benchmark shows the graphics method falls short, and the decision records why (ALG-16). | Section 9a lists the topics and the candidates. Fits the clean-room policy (section 14) because the sources are published papers, not code |
| OBJ-12 | P | Offer an interactive 3D GUI. | Headless-first architecture. If built, it dispatches the same tool ids as the CLI and MCP (AG-12) and follows FoehnCraft's GUI rules (its L7: translation keys, keyboard use, accessibility tree, theme tokens) |
| OBJ-17 | D | Be the mesher of the Craft family: every KoflCraft volume mesh is directly usable by FoehnCraft without conversion by a third tool. | Section 16. FoehnCraft defers its own general mesher (its L8) and reads meshes from files |
| OBJ-18 | L | Support multi-physics meshes: several named cell zones (fluid, solid, porous) in one conformal mesh, with named region interfaces, for conjugate heat transfer and the region model of FoehnCraft (its D112, D113). | After single-region boundary-layer meshing (M5) |

## 1a. Implementation status (2026-11-22)

Grades follow AG-17: **verified** (checked against an exact or independent reference by a named test), **partial** (works on the named fixtures, known gaps), **unverified** (built, not yet checked), **open** (not built). How each was checked is stated. A partial result is never presented as success.

**How this was checked (2026-11-22):** `cargo test --release --locked --workspace --all-targets` passed on macOS (Rust 1.99.1) and on Linux x86_64 (Rust 1.97): 31 passed, 0 failed, 1 ignored. The ignored test (`openfoam_files_parse_as_a_conforming_partition`, a dormant OpenFOAM `polyMesh` writer for the box) also passes when run on its own. `cargo m1 --require-gmsh` passed on macOS with Gmsh 4.15.2 on all six 16/64/256-site fixtures. No FoehnCraft run has used a KoflCraft mesh yet.

| Area | Grade | Evidence | Gaps |
|---|---|---|---|
| M1 box clipped Voronoi, validation, VTU and Gmsh export, CLI and MCP | verified | M1 acceptance (`cargo m1`): 16/64/256-site lattice and clustered fixtures, independent VTU and MSH readers, equal topology digests from CLI and MCP | All-pairs clipping only (ALG-04 open) |
| M2 STL import and convex STL meshing | partial | `tests/stl_workflow.rs`: unit and skewed boxes, malformed, inward, nonconvex and non-manifold rejection | Size control (ALG-05) open |
| Swept block and pipe with prism layers | partial | `mesh-cli` workflow tests, `docs/analytic-development.md` (exact block volume, pipe 0.0564 % below the cylinder) | Analytic shapes only |
| Straight channel with port-terminated layers (M3, BL-10) | partial | `tests/channel_gmsh.rs`, `research/open-port-layers.md` | Planar walls only |
| Nonconvex duct, bulk, layered and polyhedral (M4/M5 prototype) | partial | `nonconvex_elbow_and_solidcraft_…` workflow test: volume, patch areas to 1e-9, every layer reaches both ports | Needs a Gmsh scaffold and a Python preparation script (OPEN-07); fitted volume 0.183 % below the STL; flow accuracy untested |
| FoehnCraft hand-off (section 16) | open | — | MSH 2.2 written, FoehnCraft reads 4.1 only; no round-trip test |
| OpenFOAM `polyMesh` export (box only) | unverified in CI | Ignored test passes when run alone | Box patches only; not in acceptance |
| Tool registry, journal, approvals, grades (section 7) | open | MCP already confines writes to `--output-dir` | CLI and MCP dispatch separately, with different command names |
| Repository rules (section 7a) | open | — | Most of the code is not committed; no remote repository |

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
| ALG-04 | L | Efficient neighbor search instead of all-pairs clipping. | Spatial index/candidate culling; essential before large meshes. |
| ALG-05 | L | Variable mesh density controlled by a spatial size field. | Seed density and relaxation. |
| ALG-06 | P | Weighted Voronoi/power diagrams. | Evaluate for size control and boundary alignment. |
| ALG-07 | P | Delaunay tetrahedralization followed by polyhedral agglomeration. | Alternative backend if preferable for complex boundaries. |
| ALG-08 | P | Octree/cut-cell / dual contouring bulk meshing. | Alternative for highly parallel construction. |
| ALG-09 | L | Reliable geometric predicates and ambiguity resolution. | Research exact/adaptive predicates; must preserve topology. |
| ALG-10 | P | Non-convex general polyhedra. | Allow only if geometry/FVM constraints justify them. |
| ALG-11 | L | CAD/STL geometry surface queries, patch identifiers, feature edges, and inside/outside classification. | Added progressively. |
| ALG-12 | D | Separate construction of candidate convex cells from global face reconciliation and finite-volume mesh certification. | Locally valid polyhedra alone do not establish conforming global topology. |
| ALG-13 | D | Accept no missing cell fragments or non-conservative topology arising from heuristic neighbor cutoffs. | Any accelerated candidate search must be verified against an exact/reliable geometric completeness criterion or CPU reference. |
| ALG-14 | P | Conforming mirror-seed/VoroCrust-style construction rather than surface clipping. | Research alternative for future geometry, **not** the M1 algorithm. |
| ALG-15 | P | Voronoi seeds with wall-normal anisotropic spacing instead of explicit prism stacks. | Compare in CFD validation for WMLES and wall-resolved requirements. |
| ALG-16 | D | **Graphics-first selection rule.** Each algorithm choice has a design note in `research/` that names the computer graphics / geometry processing candidate (OBJ-19), the classical alternative, and a benchmark on the project fixtures that compares them for validity first (closed cells, conservation, no inverted cells), then robustness on the hard fixtures, then speed and memory. The note records the outcome even when the classical method wins. | Same rule as FoehnCraft's 2026-11-21 decision on geometry topics (robust mesh Booleans, fast winding numbers, dual contouring, narrow-band distance grids, heat and signed heat methods). Extends OBJ-15 |

**Important algorithm caveat:** The clipped-Voronoi M1 algorithm assumes a convex box. Simply clipping Voronoi cells against arbitrary non-convex domains or triangle soups does not guarantee valid CFD cells; extension needs a robust boundary-conforming strategy.

## 4. Geometry and boundary layers

| ID | Status | Requirement | Timing / proposal |
|---|---|---|---|
| GEO-01 | D | Support a programmatic analytic box domain. | M1 |
| GEO-02 | L | Import ASCII and binary STL, check closure/orientation, preserve named surface patches. | M2 |
| GEO-03 | L | Import STEP/BREP geometry while retaining CAD face/topology identification. | **Much later than M1 and boundary-layer prototypes**; no fixed early milestone |
| GEO-04 | L | Handle CAD tessellation with configurable geometric tolerances. | After STEP import |
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
| BL-10 | L | At an inlet or outlet, grow prism layers from adjoining wall patches up to the opening plane. Keep the wall-rim layer nodes on that plane. Label the terminal prism faces and the remaining core opening with the port patch. | Each terminal face has one incident cell and the boundary faces tile the original opening. Exercise on a straight channel in M3 and on the SolidCraft fixture in M5. |

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

One Cargo workspace, with a domain engine reusable by all interfaces. The crates are in layers, lowest first; a crate may depend only on lower layers, and the front ends reach the engine only through the tool registry (SW-19, AG-12). `exists` marks crates in the workspace on 2026-11-22.

```text
koflcraft/
  crates/
    L0 mesh-core/        # geometry and indexed cell/face data          (exists)
    L1 mesh-geometry/    # STL import and surface checks; later STEP     (exists)
    L1 mesh-quality/     # validation, geometric/FVM quality             (exists)
    L2 mesh-algorithms/  # Voronoi clipping, sweeps, duct layers, dual   (exists)
    L2 mesh-io/          # VTU / Gmsh / later OpenFOAM writers, readers  (exists)
    L3 mesh-engine/      # operations, plans, orchestration, evidence    (exists)
    L4 mesh-tools/       # tool registry, journal, approvals (AG-12..16) (new)
    L5 mesh-cli/         # command line, through mesh-tools only         (exists)
    L5 mesh-mcp/         # MCP server, the same tool ids                 (exists)
    -- mesh-acceptance/  # M1 acceptance runner (test tool)              (exists)
    later: mesh-gpu/     # optional GPU kernels (SW-07)
  xtask/                 # repository tasks: layering, docs (SW-19)     (new)
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
| SW-09 | D | Evaluate pure-Rust `brepkit` as a STEP research candidate only, at a much later stage. | **Not an approved dependency:** current brepkit v3+ is AGPL-3.0-only or commercial. Its earlier MIT/Apache releases require version/license audit; no dependency now. |
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
| AG-16 | D | **Summary first, paged detail.** MCP results return a summary (counts, validation, limitations, digests, file names); per-cell and per-face metrics are paged by cursor (`result.page`); meshes themselves never pass through MCP, only their file path and statistics. | FoehnCraft D198 (e). Today `mesh_report` returns every cell and face metric in one object |
| AG-17 | D | **Verification grades and limitations.** Every feature has a grade (verified, partial, unverified) and how it was checked, in a generated status page (`docs/status.md`, from `docs/features.toml`). Every result summary and `--json` output carries a `limitations` list (for example: geometry outside the validated fixtures, fitted volume differs from the STL by x %, layer coverage incomplete, flow accuracy untested). A partial result is never presented as success. | FoehnCraft D199; ArtCraft A9. Section 1a is the first manual version |
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

Research is an explicit development workstream, not a fixed algorithm selection. Papers below are relevant leads, not proof of full CFD meshing capability. Under OBJ-19 the computer graphics and geometry processing literature, above all the work of Keenan Crane (CMU) and Bruno Lévy (INRIA), is where algorithm candidates are looked for first; section 9a maps their work to project topics.

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
| P | Bruno Lévy and collaborators' Geogram Voronoi/power diagram papers | Mature conceptual foundations | Reimplement selected ideas in Rust, subject to licenses |

## 9a. Computer graphics principle: key groups and topics

Decided by the maintainer on 2026-11-22 (OBJ-19). The table maps topics KoflCraft needs to published work of the two key groups. It lists leads for the design notes of ALG-16 and commits to no algorithm. Entries cite only works already in the watchlist above or named in FoehnCraft's requirements; further papers are added to `research/` when a design note reads them.

| Topic | Candidate from the key groups | Where it enters | Status |
|---|---|---|---|
| Clipped and restricted Voronoi, power diagrams, efficient neighbour search | Lévy and collaborators: Yan, Wang, Lévy & Liu (2010); the Voronoi and power diagram papers behind Lévy's Geogram library | M1 clipper (done, all-pairs), ALG-04 neighbour search, ALG-06, nonconvex domains in M4 | Studied (2010 paper); the Geogram line is **P** |
| Robust geometric predicates and exact constructions | Lévy (2025), *Exact Predicates, Exact Constructions and Combinatorics for Mesh CSG* | ALG-09, OPEN-06, STL clipping and boundary conformity (M4) | **High priority** to study |
| Mesh Booleans and intersections on imperfect surfaces | Lévy (2025), as above | Port termination and layer collision handling (M5), imprinting for FoehnCraft's immersed-boundary plan (FC-08) | Lead |
| Signed distance, inside/outside, wall distance on imperfect surfaces | Feng & Crane (2024), *A Heat Method for Generalized Signed Distance*; Crane, Weischedel & Wardetzky (2013), the heat method for distance, as FoehnCraft cites it | Layer growth directions and local thickness limits in gaps (BL-04, BL-05), FC-08 | Lead; FoehnCraft's D108/L13 plan names the same method |
| Surface reconstruction that keeps thin features | Baktash, Gillespie & Crane, *Subgrid Marching Tetrahedra* | Preparing damaged STL input (GEO-05) | **P** |
| Surface parameterization and remeshing under layers | Corman & Crane (2025), rectangular surface parameterization, as FoehnCraft cites it | Surface mesh quality under prism layers (BL-03, M5) | Lead |

**Fit with the existing policy.**
- **Watchlist:** it already contains five works of the two groups (rows for Yan et al. 2010, Lévy 2025, Feng & Crane 2024, Baktash et al., and the Geogram line). OBJ-19 turns them from leads into the first place to look. Non-graphics entries (Cadence, Discreetize, boundary-layer surface mesh paper) stay as leads for CFD-specific needs the graphics literature does not cover, mainly layers and near-wall resolution.
- **Clean room and MIT (section 14):** nothing changes. Algorithms are built from the papers and their mathematics (LIC-02, LIC-03). The groups also publish software; LIC-09 applies to it like any other source: do not read copyleft code, permissively licensed code may be read and is never copied, and each library's licence is checked at the LIC-05 audit before anyone opens it. An author's reference implementation is not a licence to copy it.
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

## 10. Planned delivery milestones

| Milestone | Status | Scope | Exit evidence | Progress (2026-11-22) |
|---|---|---|---|---|
| **M1** | D | Rust polyhedral core; box domain; small clipped-Voronoi mesh; QA; VTU/Gmsh; CLI/MCP | 16/64/256-site fixtures pass invariants | Done, verified |
| **M1b** | D | Craft alignment: tool registry, journal, roots, grades and limitations, diagnostic codes, `AGENTS.md`, pitfalls, layering check, repository and CI (sections 7, 7a) | `kofl commands` lists every tool; CLI and MCP give identical result envelopes; `cargo xtask layering` green | New, open. Do before more meshing features, as FoehnCraft did before its MCP server |
| **M2** | L | STL input; boundary patch labeling; surface validation; basic mesh size control; FoehnCraft hand-off FC-01 to FC-03 | Box/convex-polyhedron STL cases; every fixture loads in `foehn mesh` (FC-06) | STL import and convex meshing partial; size control and hand-off open |
| **M3** | L | Simple planar-wall prism boundary layers, conformal core | [Straight channel benchmark and port-area checks](../research/open-port-layers.md#run-the-straight-channel-benchmark) | Prototype partial (channel, swept block and pipe) |
| **M4** | L | Robust STL/triangulated geometry, feature handling, clipped Voronoi in nonconvex domains | Feature and narrow-gap benchmarks, including the [SolidCraft pipe-and-box surface](../test-fixtures/solidcraft-vs-foehncraft/README.md) | Prototype partial via partitioned dual on a Gmsh scaffold (OPEN-07) |
| **M5** | L | Robust prism layer termination, collision handling and complex geometry; first FoehnCraft flow run on a KoflCraft mesh | Narrow channels, corners, curved walls, the SolidCraft port termination checks, and a FoehnCraft pipe-flow case on the SolidCraft mesh | Port termination prototype partial |
| **M5b** | L | Multi-region meshes for conjugate heat transfer (OBJ-18, FC-04) | Fluid and solid zones with a conformal interface load in FoehnCraft as two regions | Open |
| **M6** | L | Large meshes, CPU optimization, parallelization, GPU studies | 1–10 million cells, audited time and memory |
| **M7** | L | Investigate STEP import with a license-audited Rust CAD solution, after core and boundary-layer maturity | STEP vs STL geometry benchmarks |
| **M8** | P | GUI, adaptive meshing, advanced solver formats | Separate future decision |

The milestone ordering after M1 may be revised based on benchmark evidence. **STEP is explicitly deferred to M7 or later**; evaluating a library is not an agreement to depend upon it.

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
| OPEN-07 | **D** (maintainer, 2026-11-22) | The nonconvex duct path needs a Gmsh volume scaffold. Decided as proposed: keep Gmsh as an external fixture tool (allowed by LIC-09) for now, and plan a Rust tetrahedralization or scaffold-free construction for M4 so a user needs no external mesher. |
| OPEN-08 | **D** (maintainer, 2026-11-22) | Python in fixture preparation (`prepare_scaffold.py`, `check_sections.py`) versus OBJ-03. Decided as proposed, as FoehnCraft D203: OBJ-03 covers application and numerical code; Python stays allowed for fixture preparation and ParaView checks, and moves to Rust `xtask` only when next changed. |
| OPEN-09 | **D** (maintainer, 2026-11-22) | Repository: `chiefenne/KoflCraft` on GitHub, public, MIT, with the existing CI workflow, created by the maintainer; then commit the current workspace in reviewable steps. |
| OPEN-10 | **D** (maintainer, 2026-11-22) | Binary name: `kofl` with `kofl mcp`, mirroring `foehn` (section 7). |
| OPEN-11 | **D** (maintainer, 2026-11-22) | Gmsh format, as proposed: write MSH 4.1 as the default (FoehnCraft D33 reads only 4.1) and keep 2.2 as an option while the acceptance runner uses it. |
| OPEN-12 | **D** (maintainer, 2026-11-22) | Primary polyhedral hand-off to FoehnCraft: **single-file VTU first.** Rationale as proposed: VTU with cells and boundary faces in one file (FC-02) as the first format, because it keeps true polyhedra and is already written; OpenFOAM `polyMesh` (FC-03) second, since a writer exists for the box. Gmsh MSH stays the tetrahedral-subdivision format. |
| OPEN-05 | — | Closed 2026-11-22: the project is named KoflCraft (Cargo package `koflcraft`). |

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
| LIC-12 | D | The graphics-community principle (OBJ-19) works through papers and mathematics only. Software published by the key groups is treated like any third-party source under LIC-03 and LIC-09. | Maintainer principle 2026-11-22; section 9a |

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
| FC-10 | X | KoflCraft depends on FoehnCraft source crates, or FoehnCraft links KoflCraft into its solver binary. | Files and tool calls only |

