# README slideshow

Part of the [Example gallery](../gallery.md). The slides of the README's slideshow as still images. Every picture is
KoflCraft's Manifold mesh opened in [MurmelCraft](https://github.com/chiefenne/MurmelCraft)
and written by its headless `export.image` (2560 x 1600, scaled to
1440 x 900), set up by [tools/readme-tour](../../tools/readme-tour/README.md).
The mesh comes from the fixture in this repository; the command is in
[make-tour.sh](../../tools/readme-tour/make-tour.sh).

## The Manifold duct

![The Manifold duct in the dark theme, coloured by patch: the inlet pipe with its bend, the box and the outlet pipe in the wall colour, the inlet cap in orange](../../assets/readme/tour/01-manifold-patches.jpg)

The user-supplied Manifold STL (a pipe with a bend, a box, a rounded
junction and an outlet pipe) meshed by `generate-duct … polyhedral`:
229,615 cells with the inlet, outlet and wall patches.

## Cell types under a clip plane

![The Manifold cut in half along its height and coloured by cell type: polyhedra in teal fill the core of the pipes and the box, prisms in orange line the walls of the inlet pipe](../../assets/readme/tour/02-manifold-cell-types.jpg)

Voronoi polyhedra (89,467) in the core and prism layers (140,105) at the
walls, cut by a crinkle clip that keeps whole cells.

## Cell quality

![The Manifold cut flat along its height in the light theme, coloured by non-orthogonality on a pale yellow to dark red colour bar from 0 to above 60 degrees; most cells are pale, the darker ones lie along the pipe walls](../../assets/readme/tour/03-manifold-quality.jpg)

Face non-orthogonality of every cell on the cut, the check a finite-volume
solver cares about most.
