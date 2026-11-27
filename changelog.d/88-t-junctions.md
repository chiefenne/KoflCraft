- **T-junctions in block topologies.** A block corner may now sit inside an
  edge of a bigger block (`junctions` in the topology file, BS-21). The big
  block's face is gridded first and the smaller blocks take its nodes, so the
  hexahedra still meet node for node. `block.check` lists the junction edges
  and the faces they split.
- **Count solve.** One integer per chain of opposite edges, with the sum rule
  at every T-junction, fitted to target counts by an exact integer lattice
  search (`block::counts::fit`).
- **Triangle split that stops at the neighbours.** `block.split_triangle`
  takes `"stop": "t_junction"`: the neighbouring blocks stay whole with a
  T-junction on their edges and carry the new cells in their counts, instead
  of being cut in two all the way to the boundary. On the two-layer wedge:
  8 blocks and 126 cells against 10 blocks and 150 cells for the sheet split.
