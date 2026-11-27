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
| [README slideshow](gallery/readme-tour.md) | the Manifold duct (polyhedral core, prism layers) as still images with captions | — |

## All examples by number

| # | Example | Page |
|---|---|---|
| 37 | [A channel with two-sided wall clustering](gallery/block.md#channel) | [Block-structured meshes](gallery/block.md) |
| 37b | [A pipe as an O-grid](gallery/block.md#pipe-ogrid) | [Block-structured meshes](gallery/block.md) |
| 44 | [A 90° elbow](gallery/block.md#elbow) | [Block-structured meshes](gallery/block.md) |
| 44b | [A C-grid around a NACA 2412](gallery/block.md#naca-cgrid) | [Block-structured meshes](gallery/block.md) |
| 48 | [Splitting the core of a pipe](gallery/block.md#pipe-split) | [Block-structured meshes](gallery/block.md) |
| 67 | [A three-element airfoil with an optimised layout](gallery/block.md#multi-element) | [Block-structured meshes](gallery/block.md) |
| 71 | [A flat corner split into three blocks](gallery/block.md#triangle-split) | [Block-structured meshes](gallery/block.md) |
| 71b | [A triangle and a tetrahedron as hexahedra](gallery/block.md#tri-center) | [Block-structured meshes](gallery/block.md) |
