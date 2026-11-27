# Block-structured meshes

Part of the [Recipes](../recipes.md): short command sequences for the
block mesher. The inputs, templates and tools are described in
[block-mesher.md](../block-mesher.md); pictures of complete cases are on
the gallery page [Block-structured meshes](../gallery/block.md).
`kofl` is `target/release/kofl`; every tool reads one JSON file, and
outputs must not exist yet (IO003: delete or rename the old one).

| # | Recipe |
|---|---|
| 37 | [Mesh a template](#block-template) |
| 44b | [Smooth the grid and let the wall nodes slide](#block-smooth) |
| 48 | [Edit a topology: split, recount, move](#block-edit) |
| 71 | [Split a flat corner into three blocks](#triangle-split) |
| NNN | [Look at the blocks: colour by block and explode in MurmelCraft](#explode-blocks) |

---

<a id="block-template"></a>
## 37. Mesh a template

Write the topology from a template, then build the grid:

```json
{"config": {"template": "pipe_ogrid", "radius": 0.05, "length": 0.3,
            "core_cells": 12, "radial_cells": 22, "axial_cells": 30,
            "wall_spacing": 0.00015},
 "output": "pipe.json"}
```

```
kofl call block.template --input template.json
kofl call block.grid --input grid.json      # {"topology": "pipe.json", "output": "pipe-mesh"}
```

`pipe-mesh/` holds `foehncraft.vtu` (for FoehnCraft and MurmelCraft),
`mesh.msh`, `volume.vtu` (with the cell fields `BlockId`, `I`, `J`, `K`),
`report.json` and `mesh-evidence.json`. Templates: `channel`, `pipe_ogrid`,
`elbow`, `airfoil_cgrid`, `multi_element_airfoil`, `triangle`,
`tetrahedron` ([block-mesher.md](../block-mesher.md#templates)).

What to read in `report.json`: `summary.inverted_faces` (must be 0),
`min_scaled_jacobian` (1 is a box, near 0 a collapsed cell),
`max_aspect_ratio`, and per block `blocks[]`.

<a id="block-smooth"></a>
## 44b. Smooth the grid and let the wall nodes slide

```json
{"topology": "naca.json", "output": "naca-mesh",
 "algorithm": {"kind": "elliptic", "wall_sliding": true}}
```

`elliptic` smooths every node off the boundary with orthogonality and
first-spacing control at the walls; `wall_sliding` moves wall nodes along
curved walls. `report.json` `smoothing.wall_angle_deg` gives the largest
angle between grid line and wall normal before (`max_tfi`) and after
(`max`).

<a id="block-edit"></a>
## 48. Edit a topology: split, recount, move

Each edit reads a topology and writes a new, checked one:

```json
{"topology": "pipe.json", "block": 0, "axis": 2, "at": 10, "output": "pipe-split.json"}
```

```
kofl call block.split_block --input split.json
kofl call block.list --input list.json       # {"topology": "pipe-split.json"}
```

Other edits: `block.set_count`, `block.set_distribution`,
`block.set_edge_curve`, `block.move_vertex`, `block.add_block`,
`block.remove_block`, `block.set_patch` ([block-mesher.md](../block-mesher.md#editing-a-topology)).

<a id="triangle-split"></a>
## 71. Split a flat corner into three blocks

Let the trigger decide (it splits where a corner reaches 150° and the
optimised grid gets measurably better):

```json
{"topology": "wedge.json", "output": "wedge-split.json"}
```

```
kofl call block.optimize_layout --input optimize.json
```

`result.layout.triangle_splits` lists what was `kept` (block, vertex,
corner angle, worst angle before and after) and what was `skipped`, with
the reason. `"triangle_splits": false` runs the optimiser alone.

Split by hand (any corner of 120° or more):

```json
{"topology": "wedge.json", "output": "wedge-hand.json", "block": 0, "vertex": 1}
```

```
kofl call block.split_triangle --input split-by-hand.json
```

Limits: one layer of blocks (2-D extruded topologies). The tetrahedron
and prism templates are standalone shapes; a tetrahedron is not yet
injected into a larger 3-D topology.

<a id="explode-blocks"></a>
## NNN. Look at the blocks: colour by block and explode in MurmelCraft

Open the grid's `volume.vtu` (it carries `BlockId`) in
[MurmelCraft](https://github.com/chiefenne/MurmelCraft) and choose
**Colour by: Category** and **Explode** in the toolbar; the number beside
Explode sets how far the blocks move. 2-D meshes stay in their plane.

Headless, as the gallery does it (`murmel run --script shot.json`):

```json
[{"tool": "mesh.open", "input": {"path": "pipe-mesh/volume.vtu"}},
 {"tool": "display.color_by", "input": {"by": "category", "field": "BlockId"}},
 {"tool": "display.explode", "input": {"factor": 0.35}},
 {"tool": "view.standard", "input": {"view": "front_right"}},
 {"tool": "view.fit", "input": {}},
 {"tool": "export.image", "input": {"path": "pipe-exploded.png", "width": 1600, "height": 1000, "scale": 2, "theme": "light"}}]
```

Needs MurmelCraft with `display.explode` (MurmelCraft PR #124). To make
every gallery picture again:
`MURMEL=/path/to/murmel python3 tools/gallery/block_gallery.py target/release/kofl target/gallery`.
