- **Non-conformal tet split** (`block.split_tet`, BS-16, step 3b). When a
  block on one of the tet's triangles cannot be split into a prism column,
  the split now leaves it whole and makes the shared face a non-conformal
  interface (`join` `auto`, the default; `conformal_only` refuses as
  before). KoflCraft overlays the two grids on that face and writes one
  conformal mesh whose cells at the interface are polyhedra; FoehnCraft
  reads it as it is. The tet takes the coarsest count that keeps its cells
  at a scaled Jacobian of 0.4 or more. `block.optimize_layout` uses the
  same fallback.
