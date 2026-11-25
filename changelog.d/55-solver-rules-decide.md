- **Solver rules decide duct validity** (decision 2026-11-25): a duct mesh is valid when every cell has a closed shell, positive volume, positive face pyramids on its centroid and faces at less than 90° non-orthogonality, as in production finite-volume codes. Flat faces are no longer required (`maximum_warping` reports the largest warping angle). Patch areas and volume against the STL, face pairing and winding are still checked.
- **Median dual is the default polyhedral core**, as the dual of a tet mesh in production polyhedral meshers. A dual cell that breaks the solver rules is split into valid parts: its convex tet-corner pieces, regrown while the rules hold and each face shared with a neighbouring cell stays one polygon. Elbow: 7,011 cells (5 split), 0 invalid. Standard cooling jacket at size factor 2: 1,194,032 cells, 0 invalid (it failed with 58 cells on the Voronoi core).
- The Voronoi core (`"core": "voronoi"`) is on hold: its wall-join cells can break the solver rules (3 on the elbow with columns, 1 with split caps), and the build rejects them.

The solver rules now judge the mesh as it is exported: the faces between two
cells merged into one polygon and the cell centroid taken exactly from the
fan-triangulated shell, as the hand-off and FoehnCraft take them. Cells are
split against that same view. The elbow median dual splits one cell
(238 polyhedra); the Voronoi columns join passes there again.
