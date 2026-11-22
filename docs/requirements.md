# 3D Polyhedral CFD Mesh Generator
## Software Requirements Specification — v0.6

**Date:** 2026-11-18

**Purpose:** A small, Rust-native CFD mesh generator that can grow into a high-performance, agent-controllable polyhedral meshing system.  
**Document status:** Decisions consolidated 2026-11-17, with the M1 export scope and inlet/outlet layer treatment updated 2026-11-18. Items marked **O** remain proposals. **D** marks a commitment regardless of implementation status.

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
| OBJ-05 | D | Use an Artcraft-inspired *Study → Plan → Build → Check* workflow, with reproducible evidence. | From M1 |
| OBJ-06 | D | Expose the same underlying operations via Rust library, CLI, and MCP server. | Minimal set from M1 |
| OBJ-07 | D | Assess suitable computer-graphics and computational-geometry research before committing to major meshing algorithms. | Ongoing |
| OBJ-08 | L | Support demanding industrial mesh scales of approximately 1–10 million cells. | Performance milestone, not M1 |
| OBJ-09 | L | Generate prismatic boundary-layer meshes with conformal polyhedral core transitions. | After basic bulk meshing |
| OBJ-10 | L | Provide GPU-accelerated execution for algorithms shown to benefit from it. | CPU baseline first; optional GPU backend |
| OBJ-11 | X | Depend on commercial CAD/meshing tools or commercial runtime licenses. | Not planned |
| OBJ-13 | D | Release project-authored software under the MIT license. | Root `LICENSE`; repository policy |
| OBJ-14 | D | Implement algorithms independently from published mathematics and algorithms; do not copy or translate third-party source code. | Clean-room engineering process |
| OBJ-15 | D | Keep a traceable, versioned research-to-requirement-to-implementation record. | Publications and design rationale, no source-code reuse |
| OBJ-16 | L | Investigate wall-normal refinement using Voronoi seed distributions *without* a conventional prism stack. | Alternative near-wall meshing, not a replacement for committed prism-layer support |
| OBJ-12 | P | Offer an interactive 3D GUI. | Headless-first architecture |

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

One Cargo workspace, with a domain engine reusable by all interfaces:

```text
mesh-workspace/
  crates/
    mesh-core/       # geometry and indexed cell/face data
    mesh-algorithms/ # Voronoi clipping, candidate queries; backend trait
    mesh-quality/    # validation, geometric/FVM quality
    mesh-io/         # VTK / Gmsh / native formats
    mesh-engine/     # commands, plans, orchestration, evidence
    mesh-cli/        # command line
    mesh-mcp/        # MCP server, same engine
    mesh-geometry/   # later STL/STEP support
    mesh-boundary/   # later boundary layers
    mesh-gpu/        # later optional GPU kernels
  examples/
  test-fixtures/
  research/
  docs/
```

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
| SW-14 | D | No AGPL, GPL, or commercial runtime dependencies in deliverables. | Includes future CAD/GPU toolchain assessment. |

## 7. Agentic/MCP operation and Artcraft-inspired workflow

The architectural analogy is **Study → Plan → Build → Check → Approve**. It is a design influence, not a claim that Artcraft itself is a CFD mesh generator.

| ID | Status | Requirement | Example |
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

**Transport proposal:** use the maintained Rust MCP SDK `rmcp`, with a local stdio server initially. A network-accessible server is not required for M1.

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

Research is an explicit development workstream, not a fixed algorithm selection. Papers below are relevant leads, not proof of full CFD meshing capability.

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

| Milestone | Status | Scope | Exit evidence |
|---|---|---|---|
| **M1** | D | Rust polyhedral core; box domain; small clipped-Voronoi mesh; QA; VTU/Gmsh; CLI/MCP | 16/64/256-site fixtures pass invariants |
| **M2** | L | STL input; boundary patch labeling; surface validation; basic mesh size control | Box/convex-polyhedron STL cases |
| **M3** | L | Simple planar-wall prism boundary layers, conformal core | [Straight channel benchmark and port-area checks](../research/open-port-layers.md#run-the-straight-channel-benchmark) |
| **M4** | L | Robust STL/triangulated geometry, feature handling, clipped Voronoi in nonconvex domains | Feature and narrow-gap benchmarks, including the [SolidCraft pipe-and-box surface](../test-fixtures/solidcraft-vs-foehncraft/README.md) |
| **M5** | L | Robust prism layer termination, collision handling and complex geometry | Narrow channels, corners, curved walls, and the SolidCraft port termination checks |
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

**Dependency caveat:** As of 2026-11-17, `brepkit` v3+ is AGPL-3.0-only or commercial; versions through 2.129.x are reported MIT OR Apache-2.0. This only establishes possible *license* eligibility for a pinned historical version; it does not demonstrate suitability, maintenance, numerical robustness or permission to copy its implementation. The project will not depend on it at M1.

## 15. Immediate M1 engineering sequence

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
