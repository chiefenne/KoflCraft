- **Sweep regions meshed with hexahedra** (`regions.fill`). The sweep
  regions of the region proposal are filled with swept O-grids: the
  stations are level sets of the inlet-to-outlet potential (the meridian
  planes in a bend, the mitre plane at a mitre corner), each section is the
  pipe template's O-grid with the first rows at the layer heights, and a
  sweep ends at its port or at a calm, exactly planar cut in front of the
  next region. The elbow is swept from port to port (17,472 hexahedra,
  volume within 0.003 %); in the Manifold both pipes are swept and the box
  is left for the polyhedral mesher. Writes the validated FoehnCraft VTU,
  a report and the mesh evidence.
