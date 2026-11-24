# Gallery

The slides of the README's slideshow as still images. Every picture is a
KoflCraft mesh opened in [MurmelCraft](https://github.com/chiefenne/MurmelCraft)
and written by its headless `export.image` (2560 x 1600, scaled to
1440 x 900), set up by [tools/readme-tour](../tools/readme-tour/README.md).
The meshes come from the fixtures in this repository; the commands are in
[make-tour.sh](../tools/readme-tour/make-tour.sh).

## The Manifold duct

![The Manifold duct in the dark theme, coloured by patch: the inlet pipe with its bend, the box and the outlet pipe in the wall colour, the inlet cap in orange](../assets/readme/tour/01-manifold-patches.jpg)

The user-supplied Manifold STL (a pipe with a bend, a box, a rounded
junction and an outlet pipe) meshed by `generate-duct … polyhedral`:
229,615 cells with the inlet, outlet and wall patches.

## Cell types under a clip plane

![The Manifold cut in half along its height and coloured by cell type: polyhedra in teal fill the core of the pipes and the box, prisms in orange line the walls of the inlet pipe](../assets/readme/tour/02-manifold-cell-types.jpg)

Voronoi polyhedra (89,467) in the core and prism layers (140,105) at the
walls, cut by a crinkle clip that keeps whole cells.

## Prism layers at the wall

![A close view of a flat cut across the inlet pipe: large polyhedral cells in the core, a band of thin columns towards the wall and the stacked orange prism layers along the wall](../assets/readme/tour/03-manifold-prism-layers.jpg)

A flat cut across the inlet pipe near the wall: the prism layers grow off
the wall and join the polyhedral core conformally.

## Cell quality

![The Manifold cut flat along its height in the light theme, coloured by non-orthogonality on a pale yellow to dark red colour bar from 0 to above 60 degrees; most cells are pale, the darker ones lie along the pipe walls](../assets/readme/tour/04-manifold-quality.jpg)

Face non-orthogonality of every cell on the cut, the check a finite-volume
solver cares about most.

## A swept block

![A close view of a corner of the swept block: the top face shows hexagonal polyhedra in teal, the side walls show the grey hexahedral wall layers](../assets/readme/tour/05-block-sweep.jpg)

`generate-sweep test-fixtures/analytic/block.json`: a polyhedral section
swept through the height, with hexahedral layers at the walls (161,135
cells).

## A swept pipe with an O-grid

![The swept pipe seen from its end: grey hexahedral layers in a ring at the wall and a polyhedral core in teal](../assets/readme/tour/06-pipe-o-grid.jpg)

`generate-sweep test-fixtures/analytic/pipe.json`: hexahedral O-grid layers
round a polyhedral core (191,280 cells).

## A straight channel

![The straight channel in the light theme, coloured by cell type and cut lengthwise: orange prism layers on all four walls round a teal polyhedral core](../assets/readme/tour/07-channel.jpg)

`generate-channel test-fixtures/straight-channel.json`: three prism layers
on each wall round the core (776 cells).

## An elbow

![The L-shaped elbow cut flat at mid-height, coloured by cell type: teal polyhedra in the core and thin columns along the walls](../assets/readme/tour/08-elbow.jpg)

The elbow fixture meshed by `generate-duct … polyhedral` (1,004 cells), cut
flat at mid-height.
