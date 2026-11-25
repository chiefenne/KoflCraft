# BC-03 patch proposals and BC-05 boundary-condition skeleton

Date: 2026-11-25. Tools `patch.propose` and `bc.skeleton`; engine modules
`mesh_engine::proposals` and `mesh_engine::skeleton`.

## Goal

An engineer (or an agent) with a bare closed STL should name boundary
faces, not pick triangles. `patch.propose` splits the surface into the faces
a person would name, says what each one is and writes a patch map; the only
decisions left are the roles of the openings and of the planes on the
bounding box. `bc.skeleton` then writes the FoehnCraft case skeleton.

Guided flow: `patch.propose` → `patch.set_role` / `patch.rename` /
`patch.merge` (as the result's `next` line says) → `bc.skeleton`.

## Method (clean room, papers only)

1. **Regions** are the facets connected across smooth edges (dihedral turn
   at most `feature_angle_deg`, default 30°), the sharp-edge step of
   hierarchical face clustering (Garland, Willmott, Heckbert, *Hierarchical
   face clustering on polygonal surfaces*, I3D 2001).
2. **Primitive per region** from the area-weighted Gauss map
   N = Σ a nnᵀ (Attene, Falcidieno, Spagnuolo, *Hierarchical mesh
   segmentation based on fitting primitives*, The Visual Computer 2006;
   normal-cone tests of Lukács, Martin, Marshall, *Faithful least-squares
   fitting of spheres, cylinders, cones and tori*, ECCV 1998):
   - plane: every normal within `plane_angle_deg` (1°) of the mean normal;
   - cylinder: the smallest eigenvalue of N/A (mean sin² between normals
     and the axis) is at most sin²(max(1°, 0.57°)); the axis is its
     eigenvector; the axis point is the least-squares meeting point of the
     normal lines in the plane across the axis; the radius is the mean
     corner distance, accepted when the corner distances spread by at most
     `fit_tolerance` (2 %);
   - sphere: the least-squares meeting point of all normal lines, the same
     corner test;
   - otherwise freeform.
   Corners, not facet centroids, give the radius: on a CAD strip cylinder
   every corner lies on the surface, every centroid inside it.
3. **Openings** are planar regions the port detection (`ports.rs`) scores at
   least 0.9 (round pipe ends: circularity × pipe-end rim), or planar caps
   whose whole rim turns away at 60° to 120° into exactly one neighbouring
   region (a pipe of any section). **Box sides** are the other planes lying
   on a side of the bounding box (`xmin` … `zmax`).
4. **Map**: openings `opening_1…` (by area), box sides by side, everything
   else in `walls` (role wall), or one patch per region (`plane_1`,
   `cylinder_1`, `sphere_1`, `wall_1`) with `merge_walls: false`. Openings
   and box sides get no role.

## Why not variational shape approximation (VSA) here

VSA (Cohen-Steiner, Alliez, Desbrun 2004) with L²¹ (normal) proxies is the
candidate the requirement names. It partitions by planar proxies, so a
cylinder becomes k flat bands and the result depends on k and the seeds.
For boundary conditions the faces end where the engineer's faces end, at
creases; within a smooth region everything is usually one wall. The
feature-edge split is exact and deterministic, and the Gauss-map fit gives
the evidence (axis, radius) VSA would not. VSA stays the candidate for
splitting smooth regions (a plane blended into a cylinder) if a fixture
needs it.

## Results (release build, container)

| Surface | facets | regions | time | proposal |
|---|---:|---:|---:|---|
| unit box | 12 | 6 | <0.1 ms | `xmin` … `zmax`, all without role |
| CAD round pipe, 64 sides | 256 | 3 | <1 ms | 2 openings (score > 0.99), cylinder r = 0.5 exactly |
| icosphere, 3 levels | 1,280 | 1 | <1 ms | sphere with the analytic centre and radius (1e-9) |
| cooling jacket | 24,732 | 48 | 16 ms | the two 96-gon pipe ends are the openings; the eight square block faces (score 0.57) are not; cylinders r = 60, 64, 66 |
| Manifold | 90,894 | 9 | 70 ms | openings = the hand-recorded port caps; 4 outer box faces as box sides; rest `walls` |

The classical alternative, feature-angle splitting alone (`patch.split`),
gives the same regions on these surfaces but names them `<patch>_1…`
without saying which is an opening, a box side or a cylinder.

## Choices

- Box sides get no role. Decided by the maintainer on a card on 2026-11-25
  ("no role, ask the engineer"): a wrong silent wall on a symmetry plane or
  a far field is worse than one more question.
- Openings get no role; the inlet/outlet guess stays in
  `patch.detect_ports`.
- Skeleton format: FoehnCraft case format 0.7, `[mesh] file` defaults to
  `foehncraft.vtu`; inlets `velocity-inlet` with a TODO for
  `normal_velocity`, outlets `static-pressure-outlet` with a TODO for
  `pressure`; `[fluid]` only as a TODO comment. A patch without a role gets
  an empty `[boundary.<name>]` table, so FoehnCraft refuses the case until
  it has a type.

## Open

- Smooth regions mixing primitives (VSA or primitive-fitting merges).
- Cones and tori.
- Reading the skeleton with `foehn case check` in the acceptance job.
