- **Wall layers retract locally where their top layer folds** (polyhedral
  mesher, `regions.fill`). A layer cell that still fails the solver rules
  when nothing is left to split (the top of the stack in a concave
  corner, as where a tube passes through a tube plate) no longer fails
  the mesh: the layers under it are squeezed at its wall vertices and the
  mesh is built again, up to three times, keeping the full layer height
  everywhere else. The shell side of the heat exchanger now meshes with
  its five layers (before: folded layer cells at the tube plates).
