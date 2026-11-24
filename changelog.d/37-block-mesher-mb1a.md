### Added

- Block-structured mesher, first step MB1a (OBJ-21, BS-01..BS-14 in part):
  the block topology file `koflcraft.block-topology` 1.0, a check of
  indices, handedness, edge curves and laws, count chains, faces and
  patches; uniform, geometric and two-sided tanh edge laws; lines and arcs;
  transfinite interpolation of faces and blocks; point-matched hexahedra
  built without welding; validation with corner Jacobians and conservation;
  VTU, MSH 4.1, FoehnCraft VTU and evidence with block ranges; templates
  `channel` and `pipe_ogrid`; tools `block.template`, `block.check`,
  `block.grid` and the diagnostic BLOCK001.
- MB1 runs as a parallel track next to M5 (requirements section 4c).
