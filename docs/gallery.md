# Example gallery

What KoflCraft meshes, one picture per example, with the commands that
make it again. Command sequences for everyday tasks are in
[Recipes](recipes.md); real engineering geometries are in the
[Application gallery](applications.md).

All commands run from the folder of the example's inputs in
`examples/gallery/`, with `kofl` built by `cargo build --release -p mesh-cli`
(`target/release/kofl`). The pictures are the meshes opened in
[MurmelCraft](https://github.com/chiefenne/MurmelCraft) and written by its
headless `export.image`; the scripts that make them are in
[tools/gallery](../tools/gallery/README.md).

The gallery is split by topic; each page has its own table of contents.
An example's number is the pull request that added the feature it shows
("gallery 71" always means the same entry); its recipe carries the same
number.

| Page | What it covers | Examples |
|---|---|---|
| [Block-structured meshes](gallery/block.md) | multi-block hexahedral meshes from templates (channel, pipe O-grid, elbow, airfoil C-grid, three-element airfoil, triangle and tetrahedron), elliptic smoothing with wall sliding, the edit tools, the tri-center split of flat corners; blocks coloured and exploded | 37, 37b, 44, 44b, 48, 67, 71, 71b |
| [Polyhedral and layered meshes](gallery/polyhedral.md) | Voronoi and median-dual meshes: a box, a channel, analytic sweeps, a convex STL with a size field, an elbow in three modes, the Manifold duct (by hand and from a bare STL), FoehnCraft's furnace with one and seven tubes, yielding layers, a hex core, and a pipe with a baffle, a plate and a probe made tool by tool | M1, M1b, M1c, M2, M3, 26, 30, 57, 57b, 60, 70 |
| [README slideshow](gallery/readme-tour.md) | the Manifold duct (polyhedral core, prism layers) as still images with captions | — |

## All examples by number

| # | Example | Page |
|---|---|---|
| M1 | [Four Voronoi cells in a box](gallery/polyhedral.md#box) | [Polyhedral and layered meshes](gallery/polyhedral.md) |
| M1b | [A layered channel](gallery/polyhedral.md#channel) | [Polyhedral and layered meshes](gallery/polyhedral.md) |
| M1c | [A graded pipe and block](gallery/polyhedral.md#sweep) | [Polyhedral and layered meshes](gallery/polyhedral.md) |
| M2 | [A convex STL with a size field](gallery/polyhedral.md#convex-stl) | [Polyhedral and layered meshes](gallery/polyhedral.md) |
| M3 | [One elbow, three modes](gallery/polyhedral.md#elbow) | [Polyhedral and layered meshes](gallery/polyhedral.md) |
| 26 | [The Manifold duct](gallery/polyhedral.md#manifold) | [Polyhedral and layered meshes](gallery/polyhedral.md) |
| 30 | [From a bare STL](gallery/polyhedral.md#bare-stl) | [Polyhedral and layered meshes](gallery/polyhedral.md) |
| 37 | [A channel with two-sided wall clustering](gallery/block.md#channel) | [Block-structured meshes](gallery/block.md) |
| 37b | [A pipe as an O-grid](gallery/block.md#pipe-ogrid) | [Block-structured meshes](gallery/block.md) |
| 44 | [A 90° elbow](gallery/block.md#elbow) | [Block-structured meshes](gallery/block.md) |
| 44b | [A C-grid around a NACA 2412](gallery/block.md#naca-cgrid) | [Block-structured meshes](gallery/block.md) |
| 48 | [Splitting the core of a pipe](gallery/block.md#pipe-split) | [Block-structured meshes](gallery/block.md) |
| 57 | [A furnace with flush openings](gallery/polyhedral.md#furnace) | [Polyhedral and layered meshes](gallery/polyhedral.md) |
| 57b | [The seven-tube furnace](gallery/polyhedral.md#furnace-7) | [Polyhedral and layered meshes](gallery/polyhedral.md) |
| 60 | [Yielding layers and a hex core](gallery/polyhedral.md#furnace-variants) | [Polyhedral and layered meshes](gallery/polyhedral.md) |
| 67 | [A three-element airfoil with an optimised layout](gallery/block.md#multi-element) | [Block-structured meshes](gallery/block.md) |
| 70 | [Baffle, plate and probe in a pipe](gallery/polyhedral.md#pipe-shells) | [Polyhedral and layered meshes](gallery/polyhedral.md) |
| 71 | [A flat corner split into three blocks](gallery/block.md#triangle-split) | [Block-structured meshes](gallery/block.md) |
| 71b | [A triangle and a tetrahedron as hexahedra](gallery/block.md#tri-center) | [Block-structured meshes](gallery/block.md) |
