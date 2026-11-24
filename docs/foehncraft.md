# Meshes for FoehnCraft

How KoflCraft hands a volume mesh to the companion solver
[FoehnCraft](https://github.com/chiefenne/Foehncraft) (requirements section
16). The two programs meet only through files.

## Files

| File | Where | What FoehnCraft gets |
|---|---|---|
| `foehncraft.vtu` | bundles of `generate-channel`, `generate-sweep`, `generate-duct` with `--handoff vtu` (tool input `"handoff": "vtu"`), and every bundle of the block mesher; `export ARTIFACT foehncraft OUT.vtu` for box and STL meshes | The ASCII hand-off (FC-02, OPEN-12): every cell a VTK polyhedron (type 42), every boundary face a polygon (type 7), the cell array `region` (zone index for cells, patch index for faces), and the field-data string arrays `zone_names`, `boundary_names` and `boundary_roles` (`id=name`). True polyhedra, no subdivision. |
| `foehncraft.vtkhdf` | every bundle of `generate-channel`, `generate-sweep`, `generate-duct` (the default since PR #43, the maintainer's choice once FoehnCraft read it); `export ARTIFACT foehncraft-vtkhdf OUT.vtkhdf` | The primary hand-off: the same mesh, names and roles as a blocked binary VTKHDF file (see [Blocked VTKHDF](#blocked-vtkhdf)), which FoehnCraft and MurmelCraft read block by block and ParaView 5.13 or later opens as a partitioned data set. |
| `*.msh` | every bundle (`mesh.msh`, `channel.msh`); `export ARTIFACT gmsh OUT.msh` | Gmsh MSH 4.1 ASCII (FC-01): a tetrahedral subdivision of the polyhedra (prisms and pyramids where the mesh has them), one physical surface group per patch and the volume group `fluid`. `gmsh-2.2` writes the older version. |
| `export ARTIFACT openfoam DIR` | box and STL meshes | An OpenFOAM `polyMesh` directory (`points`, `faces`, `owner`, `neighbour`, `boundary`, `cellZones`) of the same face-based mesh (FC-03, so far for box and STL meshes only), with `mesh-evidence.json`. FoehnCraft reads the directory as it stands. |
| `mesh-evidence.json` | every bundle and OpenFOAM export | `koflcraft.mesh-evidence` 1.0: KoflCraft version, generator, the summary below, the patches with their roles and source names, the SHA-256 of every other file, and the limitations (FC-05). |

`foehncraft.vtu` of a box or STL export carries the artifact's
`koflcraft_mesh_digest`, `koflcraft_topology_digest`,
`koflcraft_input_digest` and `koflcraft_algorithm` as further field-data
strings, so a FoehnCraft run record that hashes the mesh file can be traced
to the KoflCraft build and inputs. `foehncraft.vtkhdf` carries the same
strings as attributes of its `/VTKHDF` group. MSH files carry no metadata;
their hash is in the bundle's `mesh-evidence.json`.

## Writing large meshes

Every file of a bundle and of an export is streamed to disk as it is
written: a second thread hashes each block of 1 MiB with SHA-256 and writes
it while the next block is formatted, and `mesh-evidence.json` (and the
`sha256` in the `mesh.export` result) takes the digests from that stream,
so no file is ever held in memory. A failed file removes the whole bundle
or export. The files are byte for byte the ones the in-memory writers made
before (checked on the fixtures and on a 1.28 million cell sweep). On a
2.5 million cell sweep block (3.4 GB of files) the peak memory fell from
6.5 GB to 3.5 GB and the run from 35 to 29 s; what remains is the meshes
themselves (see [pitfalls](pitfalls.md#whole-files-in-memory)).

## Blocked VTKHDF

`foehncraft.vtkhdf` is layout 1 of the blocked VTKHDF file agreed with
MurmelCraft and FoehnCraft (MurmelCraft's `docs/blocked-vtkhdf.md`), so a
viewer can open a mesh larger than its memory and read only the blocks
under a region. It is a VTKHDF 2.4 file (VTK documentation, "VTK File
Formats", VTKHDF section) in HDF5 written by KoflCraft's own writer
(`mesh_io::hdf5`, from the HDF5 File Format Specification 3.0).

- `/VTKHDF` is a `PartitionedDataSetCollection` with the attributes
  `Version` = [2, 4], `murmelcraft_blocks` = 1, `block_order` = `hilbert`,
  `koflcraft_version` and, for an export, the four digests above. Its
  `Assembly` group links every block.
- Block `internal` (`Index` 0, `block_cells` = 65,536) holds the cells,
  one partition per block. The cells are sorted along the Hilbert curve
  of their centroids (the mean of each cell's points; 21 bits per axis
  over the bounding box of all points, Skilling's transpose algorithm, ties
  by KoflCraft's cell number), and each run of 65,536 is a block.
- Each partition lists its owned cells, then one layer of ghost cells:
  every cell of another block that shares a point with an owned cell.
  Points are local to the partition and repeated across partitions.
- Every cell is a polyhedron (type 42) in the 2.4 arrays (`FaceOffsets`,
  `FaceConnectivity`, `PolyhedronOffsets`, `PolyhedronToFaces`); its faces
  are listed per cell and wound outward from it, so an interior face
  appears once in each of its two cells, and its `Connectivity` holds its
  distinct points in the order its faces first use them. Every integer is 64-bit except
  `Types` and `vtkGhostType` (8-bit) and `zone` (32-bit).
- `CellData`: `zone` (the `region` of the VTU), `GlobalCellIds` (the
  cell's position among the owned cells of the file, also for a ghost),
  `OriginalCellIds` (KoflCraft's cell number, the order of
  `foehncraft.vtu`), `vtkGhostType` (1 for a ghost). `PointData`:
  `GlobalPointIds` (the point index of `foehncraft.vtu`) and `vtkGhostType`
  (1 where a lower block owns the point). `FieldData`: `BlockOwnedCells`,
  `BlockBounds` (x, y, z minimum then maximum of each block's owned cells)
  and `zone_names` (`id=name`).
- Each patch with faces is a `PolyData` block named after the patch
  (`Index` 1, 2, … in patch order), with the attribute `role` (`wall`,
  `inlet`, `outlet`, `symmetry`, `interface`; none when unspecified), its
  faces as `Polygons` wound out of the domain, `PointData/GlobalPointIds`
  and `CellData/OwnerCell` (the owner's `GlobalCellIds`). A patch may not
  be named `internal` or `Assembly`.
- HDF5: superblock 2, version 2 object headers, compact links with their
  creation order (VTK lists the blocks in that order), every dataset
  contiguous and uncompressed, fixed-length strings, little-endian. The
  file is written in one pass (all headers first, then the data of each
  dataset in turn) and hashed as it goes.

On the 2.5 million cell sweep block the file is 1.7 GB (the VTU 1.3 GB:
every id takes eight bytes, and the ghosts and repeated points add about
10 %). Writing it takes about 7 s (4 s to order the cells and find the
ghosts, 3 s to write); the whole run takes 34 s instead of 31.5 s with the
VTU, at the same peak memory.

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

A patch map from the `patch.*` tools ([m2-stl.md](m2-stl.md#boundary-patches),
BC-01, BC-04) sets the role of each patch explicitly instead; its names are
already hand-off names, so every format carries them unchanged.

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

`cargo m1 --require-foehn` checks this for every acceptance fixture (35 files) and records each in `evidence.json` under `foehncraft`; CI runs it. First checked on 2026-11-22 with FoehnCraft main: box (4 sites, OpenFOAM),
straight channel, analytic block and pipe sweeps, elbow duct bulk and
layered all load valid with equal counts and volume and areas within
1.7e-12. The polyhedral elbow then loaded with equal counts but 10 faces with
S·d ≤ 0 (Manifold: 18), which FoehnCraft flagged as invalid. Since the
split of dual cells at reentrant wall edges (2026-11-22,
[note](../research/m4-reentrant-dual-split.md)) both load valid, and the
round trip fails on any export FoehnCraft calls invalid.
