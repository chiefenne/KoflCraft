# Screening: computer graphics algorithms for the polyhedral mesher

**Date:** 2026-11-22 · **Status:** screening note under OBJ-19 and ALG-16. It names candidates and a direction; it selects no algorithm. Each candidate still needs its own design note and benchmark (ALG-16) before it is built.

**How this was done.** Web search of the publications of the three key groups (Keenan Crane, CMU; Bruno Lévy, INRIA; Leif Kobbelt, RWTH Aachen) and of related work, plus the README of `nv-tlabs/wosx`. Only papers, abstracts, project pages, READMEs and licence files were read. **No third-party source code was opened** (LIC-09, LIC-13). Licences were taken from each repository's `LICENSE` file where it could be fetched, otherwise from the GitHub licence detection or the crates.io record; the method is named per row in the software register of requirements section 14a.

## 1. Recommended direction

Keep the **Voronoi / power-diagram family** as the core of the polyhedral mesher, and build it the way Lévy's group computes Voronoi cells on the GPU: each cell is clipped independently from its k nearest seeds and is accepted only when a **security radius** proves no further seed can cut it. This is what M1 already does in its all-pairs form, so it is an upgrade path, not a change of course.

Around that core:

1. **Seeds, not cells, carry the design intent.** Grading (datum objects, coordinate regions, curvature, proximity) becomes a size field; the size field becomes a seed density; seeds are relaxed by a density-weighted centroidal Voronoi tessellation (CVT). Where cell volumes must follow the size field exactly, Lévy's semi-discrete **optimal transport** gives a **power diagram with prescribed cell volumes**. Anisotropy near walls can come from **Lp-CVT** (Lévy & Liu 2010).
2. **Boundaries, thin shells and region interfaces by mirrored seed pairs.** A seed pair mirrored across a surface puts the shared Voronoi face on that surface. The same mechanism gives boundary-conforming cells (VoroCrust-style), zero-thickness baffles (duplicate the shared face as two boundary faces) and conformal multi-region interfaces. Clipping against the surface remains the fallback where the seed spacing cannot resolve a feature.
3. **Exact geometry from Kobbelt's and Lévy's robustness work.** Plane-based exact predicates (Nehring-Wirxel, Trettner & Kobbelt; Lévy 2025) for clipping, mesh Booleans, port termination and layer collisions.
4. **Prismatic layers as extruded polygons.** The wall-adjacent faces of the polyhedral surface mesh are extruded into poly-prisms (n-gon prisms); the last layer is the inner boundary that the Voronoi core is clipped against. Marching directions and local thickness limits come from smooth fields on the surface and in the volume (heat method, signed heat method).
5. **Grid-free Monte Carlo (Walk on Spheres, `wosx`) for fields before a mesh exists.** Size-field smoothing, harmonic marching-direction fields and inside/outside queries can be evaluated point by point on the raw surface, in parallel and on the GPU, without a background mesh. It is a field tool, never a topology tool, because its answers are noisy.
6. **Surfaces first.** Isotropic remeshing (Botsch & Kobbelt 2004), intrinsic Delaunay ideas (Crane's group) and patch segmentation give the surface on which layers grow and on which an agent names boundary patches.

Where the graphics literature does not reach, classical CFD methods stay candidates: the Delaunay-tetrahedra dual (the polyhedral approach of the large commercial meshers) as the alternative backend of ALG-07, and advancing-layer heuristics for layer collisions.

## 2. Candidates by requirement

| Requirement | Key-group candidate | Other candidate | Speed or quality gain expected | Risk |
|---|---|---|---|---|
| Polyhedral core, scaling (ALG-01, ALG-04, ALG-13) | Ray, Sokolov, Lefebvre & Lévy (2018), *Meshless Voronoi on the GPU*: per-cell clipping by k nearest neighbours with a security-radius completeness test; Basselin et al. (2021), *Restricted Power Diagrams on the GPU* | Taveira et al. (2026), *Scalable GPU construction of 3D Voronoi and power diagrams* (already watchlist) | Replaces all-pairs O(n²) with near O(n) and makes the cell loop embarrassingly parallel (CPU threads now, GPU later); the security radius is the completeness proof ALG-13 asks for | Degenerate seeds (co-spherical lattices) need exact predicates for shared-face reconciliation |
| Nonconvex domains (M4, OPEN-07) | Yan, Wang, Lévy & Liu (2010), clipped Voronoi on a tetrahedralised domain; restricted Voronoi diagrams | Diazzi & Attene (2021), *Convex polyhedral meshing for robust solid modeling*; mirrored-seed conformity (VoroCrust) | Removes the Gmsh scaffold once a Rust tetrahedral or convex decomposition exists | VoroCrust conforming construction is patented (LIC-08, below) |
| Mesh grading (GRD-01 to GRD-07) | Density-weighted CVT; Lévy (2015), semi-discrete optimal transport in 3D; Lévy (2022), *Partial optimal transport for a constant-volume Lagrangian mesh with free boundaries* (power diagrams clipped by a domain); Lévy & Liu (2010), *Lp CVT* for anisotropy | Persson (2006), mesh size functions with gradient limiting | OT gives cell volumes that match the size field directly instead of by trial; Lp-CVT gives direction-aligned cells | OT needs a Newton solve over all seeds; cost to be measured |
| Size-field evaluation and smoothing (GRD-05) | Sawhney & Crane (2020), *Monte Carlo geometry processing*; Sawhney et al. (2023), *Walk on Stars*; `wosx` | Background octree with gradient limiting | Point-wise, mesh-free, GPU-parallel evaluation on the raw STL | Monte Carlo noise; needs a fixed RNG seed for M1-11 determinism |
| Prismatic layers (BL-01 to BL-05, BL-11, BL-12) | Crane, Weischedel & Wardetzky (2013), heat method; Feng & Crane (2024), signed heat method, for marching directions, gap detection and local thickness limits; harmonic direction fields by WoS | Advancing-layer methods with local reconnection (watchlist, 2026) | Smooth, collision-aware marching directions from one linear solve or point-wise Monte Carlo | Layer collision at concave corners remains combinatorial work |
| Thin shells / baffles (SHL-01 to SHL-05) | Mirrored seed pairs across the shell; exact mesh Booleans to imprint the shell into the core (Trettner, Nehring-Wirxel & Kobbelt 2022, *EMBER*; Lévy 2025) | Commercial baffle meshing (clip and duplicate faces) | One mechanism for boundaries, shells and interfaces | Shell edges (free boundaries inside the fluid) need a local rule; layers on both sides of a shell |
| Multi-domain (MD-01 to MD-05) | Conformal: one Voronoi diagram, cells tagged by region, interfaces by mirrored seeds. Non-conformal: independent diagrams, interface faces with exact polygon intersection (Kobbelt/Lévy robust Booleans) | Arbitrary mesh interface (AMI) style overlap weights | Conformal case comes nearly for free from the core | Non-conformal interface weights must be conservative to round-off |
| Robust predicates (ALG-09, OPEN-06) | Nehring-Wirxel et al. (2025), *Exact and efficient mesh-kernel generation*; Nehring-Wirxel, Trettner & Kobbelt (2021), *Fast exact Booleans for iterated CSG using octree-embedded BSPs*; Lévy (2025), *Exact predicates, exact constructions and combinatorics for mesh CSG* | Shewchuk adaptive predicates (`robust` crate, MIT OR Apache-2.0) | Plane-based representation keeps clipped vertices exact without growing coordinates | Bit-length budget of plane-based predicates must fit `i128`/expansions |
| Surface preparation | Botsch & Kobbelt (2004), isotropic remeshing; Botsch, Kobbelt, Pauly, Alliez & Lévy (2010), *Polygon Mesh Processing* (book, the shared reference of two key groups); Feng, Gillespie & Crane (2023), winding numbers on discrete surfaces; Baktash, Gillespie & Crane, *Subgrid Marching Tetrahedra* | Jacobson, Kavan & Sorkine (2013), generalized winding numbers | Robust inside/outside on dirty STL | — |
| Boundary faces for FoehnCraft by agents (BC-01 to BC-06) | Feature-angle region growing and variational shape approximation (Cohen-Steiner, Alliez & Desbrun 2004) to propose patches | — | The agent names and merges proposals instead of picking triangles | Proposals must be reproducible and journalled (AG-13) |

## 3. Where these algorithms could speed up the meshing task

- **All-pairs clipping to k-nearest clipping** is the largest single gain: the current M1 clipper visits every seed pair; a k-d tree or grid plus the security radius visits about 30 to 60 neighbours per cell.
- **Per-cell independence** makes the core parallel without locks; global face reconciliation stays a separate pass (ALG-12).
- **GPU** later for the cell loop (Ray et al. 2018; Taveira et al. 2026) and for WoS field evaluation (`wosx` has a GPU path), under GPU-01 and GPU-09.
- **Optimal transport instead of iterative seed tuning** for grading: fewer full rebuilds of the diagram.
- **Exact plane-based predicates** avoid the slow path of arbitrary-precision coordinates.

## 4. `nv-tlabs/wosx` in detail

Read: README and `LICENSE` only (2026-11-22). WoSX is a header-only C++17 library with GPU support and Python bindings for Walk on Spheres-style Monte Carlo PDE solvers, forked from Rohan Sawhney's *Zombie* (MIT). It queries the boundary representation directly and estimates the solution only where it is needed. **Licence: Apache-2.0** (with an explicit patent grant). Dependencies named in the README: FCPW (MIT), nanoflann (BSD), oneTBB (Apache-2.0), pcg32 (Apache-2.0); demos use polyscope (MIT), stb (MIT), nlohmann/json (MIT); bindings use nanobind (BSD). No copyleft dependency.

Uses for KoflCraft: size-field diffusion from datum objects (Laplace/Poisson with Dirichlet data on the datum surfaces), harmonic fields for layer marching directions, smooth wall-distance-like fields, and possibly a potential-flow initial guess for FoehnCraft (outside KoflCraft's scope). Not a use: any decision about topology, cell validity or conservation. Clean room: we implement WoS from the papers (Muller 1956; Sawhney & Crane 2020; Sawhney et al. 2023), not from the WoSX source. Apache-2.0 would allow reading it under LIC-09, but LIC-13 asks us to start from the papers and open the code only when a design note records why.

## 5. Patents found during screening (for LIC-08)

- US 10,776,540 B1, *Constructing a conforming Voronoi mesh for an arbitrarily-shaped enclosed geometric domain* (DOE/Sandia record on OSTI). It covers the VoroCrust line of mirrored-seed conformity. The mirrored-seed items above need a freedom-to-operate check before commercial distribution.
- US 12,094,059, *Anisotropic polyhedra boundary layer adaptation method*: found by search, not read; to be assessed before BL-12 work.

## 6. Sources

Works named in the tables without a link (Lévy 2015 on semi-discrete optimal transport in 3D, Botsch & Kobbelt 2004, the *Polygon Mesh Processing* book, Crane et al. 2013, Persson 2006, Cohen-Steiner et al. 2004, Muller 1956, Jacobson et al. 2013) are standard references cited from the literature; they were not re-fetched for this note, and the design note that uses one adds its DOI.


- Ray, Sokolov, Lefebvre, Lévy (2018), *Meshless Voronoi on the GPU*, ACM TOG: https://hal.archives-ouvertes.fr/hal-01927559
- Basselin et al. (2021), *Restricted Power Diagrams on the GPU*, CGF: https://diglib7.eg.org/bitstream/handle/10.1111/cgf142610/v40i2pp001-012.pdf
- Lévy & Liu (2010), *Lp Centroidal Voronoi Tessellation and its applications*: https://hal.archives-ouvertes.fr/inria-00600251
- Lévy (2022), *Partial optimal transport for a constant-volume Lagrangian mesh with free boundaries*, JCP: https://arxiv.org/abs/2106.03936
- Trettner, Nehring-Wirxel, Kobbelt (2022), *EMBER: exact mesh booleans via efficient and robust local arrangements*: https://www.graphics.rwth-aachen.de/publication/03339
- Nehring-Wirxel, Trettner, Kobbelt (2021), *Fast exact Booleans for iterated CSG using octree-embedded BSPs*: https://arxiv.org/abs/2103.02486
- Nehring-Wirxel et al. (2025), *Exact and efficient mesh-kernel generation*: https://www.graphics.rwth-aachen.de/publication/03355/
- Diazzi & Attene (2021), *Convex polyhedral meshing for robust solid modeling*: https://ar5iv.labs.arxiv.org/html/2109.14434
- Sawhney et al. (2023), *Walk on Stars*: https://www.cs.cmu.edu/~kmcrane/Projects/WalkOnStars/
- Monte Carlo geometry processing resources (CC-BY-4.0): https://rohan-sawhney.github.io/mcgp-resources/
- WoSX: https://github.com/nv-tlabs/wosx
- VoroCrust seed placement and patent: https://www.osti.gov/biblio/1735266-constructing-conforming-voronoi-mesh-arbitrarily-shaped-enclosed-geometric-domain, https://patents.google.com/patent/US10776540B1/en
- Keenan Crane: https://www.cs.cmu.edu/~kmcrane/ · Bruno Lévy: https://brunolevy.github.io/ · RWTH Visual Computing: https://www.graphics.rwth-aachen.de/publications/
