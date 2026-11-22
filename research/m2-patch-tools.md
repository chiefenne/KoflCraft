# M2 boundary patch tools (BC-01, BC-02)

Design note of ALG-16 for defining boundary patches on an imported STL after import, so that an agent can name the boundary faces of a FoehnCraft case. The outcome is built in [`mesh_geometry::patches`](../crates/mesh-geometry/src/patches.rs) (the geometry), [`mesh_engine::patches`](../crates/mesh-engine/src/patches.rs) (the patch map and its operations) and [`mesh-tools` `patch.*`](../crates/mesh-tools/src/patch.rs) (the tools); usage is in [m2-stl.md](../docs/m2-stl.md#boundary-patches).

## What BC-01 asks for, and why no algorithm choice is open

BC-01 names the operations: STL solid names, feature-angle region growing, plane or box selection, normal direction, merging and splitting. Each is either an exact predicate on the facets or a flood fill over the shell's edge adjacency:

| Operation | Rule |
|---|---|
| Solid names (`patch.init`) | One patch per STL solid, the name made a FoehnCraft boundary name (`handoff_name`, FC-04), the role from the name |
| Plane (`patch.select_plane`) | All three corners within `tolerance` of the plane (default 1e-6 of the longest extent, above the f32 rounding of binary STL) |
| Box (`patch.select_box`) | Centroid in the closed box, or all corners (`contain: all_corners`) |
| Normal (`patch.select_normal`) | Angle between the unit outward facet normal and the direction at most `max_angle_deg` |
| Region growing (`patch.grow`) | Flood fill from a seed facet across every edge whose dihedral angle (the angle between the two facet normals) is at most `feature_angle_deg`; the seed is a facet index or the facet nearest to a point by exact point–triangle distance (Ericson 2005, section 5.1.5) |
| Split (`patch.split`) | The same flood fill restricted to the patch, every region a new patch; without an angle, edge connectivity only |
| Merge, rename, set role | Relabelling |

Feature-angle growing is the classical rule that CAD and CFD preprocessors use for feature edges. The graphics literature has no replacement for an *explicit* selection: the operations are defined by the requirement, they are exact, and an agent must be able to predict their result. So ALG-16's comparison against a graphics candidate does not apply to BC-01; the check is against an independent reference instead (below).

Where the graphics literature does compete is **BC-03, proposing a segmentation**: variational shape approximation (Cohen-Steiner, Alliez & Desbrun 2004) and the feature-aware chart growing of Lévy, Petitjean, Ray & Maillot (2002), both covered with feature detection in Botsch, Kobbelt, Pauly, Alliez & Lévy, *Polygon Mesh Processing* (2010). They are the candidates to benchmark against plain feature-angle splitting when BC-03 is built. No source code was read for this note (LIC-13).

## Checks

- **Independent reference.** The SolidCraft surface (90,894 facets, sharp rims, rounded pipe bends) has its two port caps recorded by hand in [its patch map](../test-fixtures/solidcraft-vs-foehncraft/patches.json). Region growing at 30° from a point on each port returns exactly the recorded facet ranges (510 and 512 facets), the inlet's area matches the recorded area to 1e-12 and its mean normal is exactly −Y; the plane `y = 0` selects the same inlet set (`region_growing_and_planes_recover_the_solidcraft_ports`).
- **Exact evidence.** On the unit box the summary reports areas 1, axis normals with zero spread, one region and four neighbours per side (`solid_names_start_the_map_with_exact_evidence`). Splitting the whole skewed box at 10° returns its six sides, two facets each.
- **Every facet in one patch (BC-06).** Each read of a map checks the STL digest, the facet count, one valid index per facet, unique valid names and no empty patch; an edited map fails (`invalid_maps_and_empty_operations_are_diagnosed`).
- **Names in every format.** A mesh built with a map carries its names in MSH 4.1 physical groups, the FoehnCraft VTU and OpenFOAM `boundary`, and its roles in the VTU, the OpenFOAM types and the evidence (`patch_map_names_and_roles_reach_every_export`). `cargo m1 --require-foehn` reads such a mesh in all three formats with `foehn mesh` and compares counts, volume and patch areas.
- **Cost.** On the SolidCraft surface, release build, Linux x86_64, 2026-11-22: `patch.init` 0.21 s, `patch.grow` 0.18 s and `patch.split` at 30° (eight regions: the outlet cap, the pipe and box walls, the flat box sides) 0.29 s per call, each including the STL import and the summary.

## Choices made without a decision card

- One map per call: every tool writes a new file and never edits its input, as every KoflCraft output is new; the journal of the calls is the history.
- A selection moves the facets into the target patch, whatever patch they were in, unless `within` restricts it; patches left empty are dropped and reported.
- A rename keeps the role. A new patch without a `role` gets the role its name suggests, as the STL solids do.
- An operation that selects or changes nothing is an error (PATCH004), so an agent never mistakes a no-op for a step.

## Open

- The nonconvex duct path (`mesh.generate_duct`) still reads its older patch map (a default patch and two port ranges). Moving it to `koflcraft.patch-map` is a separate change.
- BC-03 (proposals), BC-05 (boundary-condition skeleton) and the `shell` role of BC-04 are not built.
