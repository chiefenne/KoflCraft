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

The partition defines the six regions for the tetrahedral benchmarks and how the polyhedral core gets its volume scaffold. Since 2026-11-23 the scaffold is built in Rust from this STL (`volume_scaffold.build`): the wall is remeshed onto the STL to a chord tolerance of 0.05 model units, with spacing 0.65 to 6 and 0.85 at the port rims, and the volume is tetrahedralized and refined. `max_wall_deviation` (0.06) fails the run if the wall leaves the STL by more. No Gmsh file or preparation script is needed; the earlier Gmsh scaffold followed the STL only to 0.54 units ([note](../../research/m4-rust-scaffold.md)). [Partition](partition.json) · [Generation commands and measured port coverage](../../docs/nonconvex-duct.md)
