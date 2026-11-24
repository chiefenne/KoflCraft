# Thin gaps and thin walls

Design note for the thin-part work of 2026-11-23/24. Question: why does the scaffold's boundary recovery blow up on a cooling jacket with 0.5 mm gaps and ribs (splits grow about 3.5 times per pass until memory runs out), and what makes such parts mesh reliably? The fixture is [`cooling-jacket-thin.stl`](../test-fixtures/cooling-jacket/README.md): one rib leaves a 0.5 mm gap under the outer wall, one rib is 0.5 mm wide.

## Local thickness for any shape

The gap source of [GRD-04](grd-04-geometry-sizing.md) is exact for convex domains only. [`sizing::thickness`](../crates/mesh-algorithms/src/sizing/thickness.rs) finds the medial ball at a wall point by the shrinking-ball method (Ma, Bae & Choi 2012) with exact closest points on the triangles: the ball on the fluid side gives the gap, the ball on the solid side the thickness of a thin wall. Only opposite walls stop the ball (normals within 45° of opposite), the same rule as in GRD-04, so sharp edges do not collapse the thickness. On the fixture it reports 0.4997 and 0.5000.

## What made recovery blow up

Measured on the thin jacket at `size_factor` 2 (missing = wall edges absent from the Delaunay tetrahedralization of the wall vertices before any split):

| Thickness sizing | Chord tolerance | Remesh rounds | Missing edges | Recovery |
|---|---|---|---|---|
| off | 0.016 | 10 | 724 | blows up |
| on | 0.016 | 10 | 174 | blows up |
| off | 0.04 | 25 | 237 | blows up after 55 passes |
| on | 0.016 | 25 | 55 | blows up after 30 passes |
| on | 0.04 | 25 | 0 | not needed |

Three causes, each needed:

1. **Wall edges longer than the local thickness.** The thickness source sizes the wall at one edge across.
2. **A remesh that had not converged.** With a strongly graded target, ten rounds left 19,000 edges below a quarter of their target; the loop now goes on until none remain (25 to 27 rounds here, 40% fewer triangles).
3. **A tolerance below the STL's own facet error.** The jacket's 96-sided cylinders are 0.035 off the round wall. With a tolerance of 0.016 the deviation pass bisected edges towards the facet lines and left ladders of needles. `mesh.setup_duct` now floors the tolerance at 1.25 times the facet error.

In addition the remesher keeps the diametral ball of every sharp edge free of the apexes of its two facets (the segment rule of Ruppert 1995 on the surface): a vertex 0.04 beside a 2 long sharp edge kept that edge out of every tetrahedralization.

**Still open:** any leftover missing edge ends in a blowup, sooner or later. The guarantee is protecting balls with a weighted Delaunay triangulation (Cohen-Steiner, de Verdière & Yvinec 2002; Cheng, Dey & Shewchuk 2012), chosen by the maintainer and not built yet.

## Layers in thin places

The maintainer chose local squeezing with the same layer count everywhere (the height scaling of Garimella & Shephard 2000 and Dyedov et al. 2009) over dropping layers. Each wall vertex gets a share of the layer height: at most 30% of the gap in front of it, times 0.7 for the vertices of front facets that still intersect, times 0.7 for the wall vertices of scaffold tetrahedra that the wall offset turns over, with the height changing by at most 0.3 per unit length along the wall. Port rims keep the full height. On the thin jacket 41 of 71,432 wall vertices are squeezed, the least to 34%.

## Not solved

The Voronoi core does not complete on either jacket ("a wall ring edge is used twice", "piece has no continuation"). Two core failures were fixed on the way: a convex piece thinner than the kernel search margin (its centroid stands in), and two scaffold points one rounding step apart from cospherical CAD rings (points within 1e-12 of the box diagonal of a vertex now count as that vertex).
