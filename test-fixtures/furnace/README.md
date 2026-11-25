# Box furnace

The geometry of FoehnCraft's radiation gallery furnace (example 259, `examples/meshes/furnace.py` there), fluid volume only, in m: a chamber 4 long (x), 2 deep (y) and 2.5 high (z), horizontal tubes of 0.16 diameter through it under the roof, a burner disc of 0.5 diameter in the front wall and a rectangular opening in the back wall under the roof. The burner and the opening are flush with their walls. The fixture has one tube (13,660 facets); the gallery furnace has seven.

`make_furnace.py` builds the STL and its patch map (`furnace.patches.json`: `walls`, `front`, `tubes`, `burner` as inlet, `outlet`) with Gmsh: `python3 -m pip install gmsh numpy && python3 make_furnace.py furnace.stl 1`. `partition.json` is a coarse built scaffold with three layers for the fast test. [Flush openings and named walls](../../docs/nonconvex-duct.md#flush-openings-and-named-walls)
