- **Triangle split in 3D.** `block.split_triangle` now works on any block
  topology, not only extruded ones. The flat corner is a block edge (new
  input `along`, or the flattest edge at `vertex`): every block stacked
  along that edge becomes three blocks around a centre line, and the new
  cell count runs as a sheet through the neighbouring blocks to the
  boundary, so all faces stay point-matched. Blocks in any local frame and
  any direction in space are handled.
