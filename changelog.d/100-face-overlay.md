- **Face overlay for non-conformal block joins** (BS-16, step 3a). A new
  kernel overlays the two grids that the two sides of a block face put on
  it, in the face's own parameters, with exact predicates and a 2 % snap
  against sliver faces. Each overlap polygon becomes one face between its
  two cells, so the flux through the interface is conserved to round-off.
  The block mesher uses it in the next step.
