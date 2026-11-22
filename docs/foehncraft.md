# Meshes for FoehnCraft

How KoflCraft hands a volume mesh to the companion solver
[FoehnCraft](https://github.com/chiefenne/Foehncraft) (requirements section
16). The two programs meet only through files.

## Files

| File | Where | What FoehnCraft gets |
|---|---|---|
| `foehncraft.vtu` | every bundle of `generate-channel`, `generate-sweep`, `generate-duct`; `export ARTIFACT foehncraft OUT.vtu` for box and STL meshes | The primary hand-off (FC-02, OPEN-12): every cell a VTK polyhedron (type 42), every boundary face a polygon (type 7), the cell array `region` (zone index for cells, patch index for faces), and the field-data string arrays `zone_names`, `boundary_names` and `boundary_roles` (`id=name`). True polyhedra, no subdivision. |
| `*.msh` | every bundle (`mesh.msh`, `channel.msh`); `export ARTIFACT gmsh OUT.msh` | Gmsh MSH 4.1 ASCII (FC-01): a tetrahedral subdivision of the polyhedra (prisms and pyramids where the mesh has them), one physical surface group per patch and the volume group `fluid`. `gmsh-2.2` writes the older version. |
| `export ARTIFACT openfoam DIR` | box and STL meshes | An OpenFOAM `polyMesh` directory (`points`, `faces`, `owner`, `neighbour`, `boundary`, `cellZones`) of the same face-based mesh (FC-03, so far for box and STL meshes only), with `mesh-evidence.json`. FoehnCraft reads the directory as it stands. |
| `mesh-evidence.json` | every bundle and OpenFOAM export | `koflcraft.mesh-evidence` 1.0: KoflCraft version, generator, the summary below, the patches with their roles and source names, the SHA-256 of every other file, and the limitations (FC-05). |

`foehncraft.vtu` of a box or STL export carries the artifact's
`koflcraft_mesh_digest`, `koflcraft_topology_digest`,
`koflcraft_input_digest` and `koflcraft_algorithm` as further field-data
strings, so a FoehnCraft run record that hashes the mesh file can be traced
to the KoflCraft build and inputs. MSH files carry no metadata; their hash
is in the bundle's `mesh-evidence.json`.

## One face between two cells

FoehnCraft, like most finite-volume solvers, refuses two cells that share
more than one face. KoflCraft's polyhedral duct dual has such pairs (the
pieces of a bent interface). The hand-off merges the pieces into one
polygon: edges shared by two pieces cancel, the rest must form one loop,
otherwise the export fails. The area vector of the merged face is the sum of
the pieces, so volumes and fluxes are unchanged. `merged_faces` in the
summary counts the pieces merged away.

## Patch names and roles (FC-04)

FoehnCraft names boundaries in case files as `[boundary.<name>]`, so a
hand-off name is a TOML bare key: ASCII letters, digits, `_` and `-` are
kept and every other character becomes `_` (`west inlet` → `west_inlet`);
a repeated name gets a suffix (`_2`). The original name is kept as
`source_name` in the evidence. The role comes from the name and suggests
the FoehnCraft boundary-condition type a case starts from:

| Name contains | Role | FoehnCraft `type` | OpenFOAM `type` |
|---|---|---|---|
| `inlet` | inlet | `velocity-inlet` | `patch` |
| `outlet` | outlet | `static-pressure-outlet` | `patch` |
| `wall` | wall | `wall` | `wall` |
| `symmetry` | symmetry | `symmetry` | `symmetry` |
| `interface` | interface | (region interface, OBJ-18) | `patch` |
| anything else (`xmin` … `zmax` of the box) | unspecified | none: the case decides | `patch` |

The role is a suggestion recorded in the files, not a boundary condition;
the case file decides. Cell zones are named (`fluid` today) for FoehnCraft
regions (its D112).

## What FoehnCraft must reproduce (FC-06)

The summary in `mesh-evidence.json` has the counts of cells, faces,
interior and boundary faces, the total volume, the faces and area of every
patch, the largest relative closure |Σ n·dS| of a cell, and
`inverted_faces`, the interior faces with S·d ≤ 0 (area vector against the
line from owner to neighbour centroid). KoflCraft computes the geometry with
the same decomposition as FoehnCraft (tetrahedra from a cell point, the face
node average and an edge), so `foehn mesh` must report the same counts and
the same volume and areas to round-off, and calls a mesh invalid exactly
when `inverted_faces` is not zero.

Checked by hand on 2026-11-22 with FoehnCraft main: box (4 sites, OpenFOAM),
straight channel, analytic block and pipe sweeps, elbow duct bulk and
layered all load valid with equal counts and volume and areas within
1.7e-12. The polyhedral elbow loads with equal counts but 10 faces with
S·d ≤ 0, which FoehnCraft (and KoflCraft's prediction) flag as invalid; that
is a quality gap of the polyhedral dual prototype (M4), not of the files.
