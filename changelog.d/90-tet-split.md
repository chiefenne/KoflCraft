- **Tet split of a block.** `block.split_tet` splits a block shaped like a
  tetrahedron into the four hexahedra of the tetrahedron. A block on one of
  the three split triangles is split into three with it (a prism column), so
  the cells still meet node for node, and the new cells run on as one row,
  cut as a sheet or stopped with T-junctions. The 3D trigger of
  `block.optimize_layout` tries tet-shaped blocks too and keeps a split only
  when it measurably helps (on the test tet: worst scaled Jacobian 0.04 to
  0.49, worst angle 82° to 51°).
- The triangle split now puts its centre at the mean of the three points on
  the triangle's sides, which keeps its three quadrilaterals convex.
