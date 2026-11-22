# SolidCraft vs FoehnCraft

The user-supplied [binary STL](SolidCraft_vs_FoehnCraft.stl) is stored byte for byte. It contains 90,894 triangles and 45,449 distinct vertices. It has one closed nonconvex shell with pipe roundings and sharp wall edges. Its bounds span `x = −120 … 380`, `y = 0 … 400`, and `z = −40 … 80`, in model units.

The STL already includes the two triangulated circular end caps. [The patch map](patches.json) records the source SHA-256 and assigns an explicit role to each cap. The ranges use zero-based facet indices with exclusive upper bounds, in STL record order.

| Patch | Plane and outward normal | Center | Cap triangles |
|---|---|---|---:|
| Inlet | `y = 0`, `−Y` | `(0, 0, 0)` | 510 |
| Outlet | `x = −120`, `−X` | `(−120, 360, 40)` | 512 |
| Walls | Remaining 89,872 triangles | | |

Both caps have diameter 50 and a rim of 240 edges. The inlet/outlet assignment is a fixture convention. The geometry supplies two outward directions and the patch map supplies the flow roles.

`cargo test --locked -p mesh-geometry --test pipe_ports` verifies the original bytes, closed shell topology, cap classification and winding, area, circular rims, and a mix of sharp and smoothly varying wall dihedrals. The fixture also serves as a negative case for the convex volume generator.

The partition defines the six regions for the tetrahedral benchmarks and identifies the Gmsh scaffold for the distributed polyhedral core. The scaffold was prepared separately from this STL with a target bulk spacing of six model units and local refinement at the port rims. [Partition](partition.json) · [Preparation script](prepare_scaffold.py) · [Generation commands and measured port coverage](../../docs/nonconvex-duct.md)
