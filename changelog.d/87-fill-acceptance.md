- **Sweep fill in the acceptance run** (`regions.fill`, FC-06). The full
  CI run now meshes the elbow (hexahedra port to port) and the joined
  Manifold (swept pipes and the polyhedral box) and has FoehnCraft read
  both back: same cells, faces, patches and volume, no inverted face, no
  cut left as a patch. Each case's non-orthogonality on both measures is
  saved in `fill/quality.json`.
- **`regions.fill` names what it wrote.** The joined Manifold was told
  that its box was not meshed, and the elbow that it held polyhedra; the
  next step now says which of the three meshes the file is.
- **Gmsh file of the polyhedral Manifold readable again.** Four core
  cells had a corner as their kernel point, so their tetrahedral fan in
  `mesh.msh` held tetrahedra of zero volume and FoehnCraft refused the
  file. Those cells now get a kernel point strictly inside.
