- **Region proposal for block-dominant meshing** (`regions.propose`). From a
  closed STL and its patch map, KoflCraft builds a tet scaffold, solves the
  wall distance (heat method) and the inlet-to-outlet potential, finds the
  skeleton with its end and branch points, and marks each region as sweep,
  junction, obstacle, box or remainder (left to the polyhedral mesher). A
  second pass tries each box and remainder region again along its own long
  axis. Writes a report, the fields (VTU) and a picture (SVG); nothing is
  meshed yet.
