- **Box joined to the swept pipes** (`regions.fill`, MD-25). The region
  behind the cuts of the sweeps is now meshed as polyhedra (median dual
  with the sweeps' wall layers) and joined conformally to the hexahedra:
  each cut gives the polyhedral mesher a planar cap, the pipe's last
  section takes the positions the mesher gives that cap, and every hex
  face at a cut is split into the polyhedral faces it holds. The Manifold
  is one valid mesh of 132,520 cells (69,632 hexahedra in the pipes,
  62,888 polyhedra in the box), with no inverted face and its volume
  within 0.2 % of the STL. `join: false` meshes only the sweeps, as
  before.
- **Polyhedral wall layers at box corners.** Where the wall offset kept a
  scaffold tetrahedron turned over at a corner of three walls, the layers
  of the third wall now yield too, so thick layer stacks no longer stop
  the mesh there.
