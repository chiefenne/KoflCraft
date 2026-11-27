# Block-structured meshes

KoflCraft's second mesher (OBJ-21, milestone MB1, requirements section 4c)
builds multi-block hexahedral meshes of classic CFD shapes: channels,
pipes, elbows and airfoils. A **block topology** (`koflcraft.block-topology`
1.0, a JSON file) holds the vertices, the hexahedral blocks with their cell
counts, the edge curves and spacing laws and the boundary patches.
`block.grid` turns it into a validated mesh. Everything runs through the
tool registry, so the same calls work from `kofl call`, from a journal
script and from an agent over MCP ([tools.md](tools.md)).

Pictures and commands of every template: [Example gallery](gallery/block.md); short recipes: [recipes/block.md](recipes/block.md).

The concept and the literature are in
[research/block-structured-concept.md](../research/block-structured-concept.md);
the grade and the gaps of every part are on the [status page](status.md).

## From a template to a mesh

Build `kofl` once (`cargo build --release -p mesh-cli`, binary
`target/release/kofl`). Every tool reads one JSON input file; outputs must
not exist yet.

1. Write a topology from a template. `template-airfoil.json`:

   ```json
   {"config": {"template": "airfoil_cgrid", "naca": "2412", "span": 0.1,
               "farfield": 5, "wake_length": 5, "surface_cells": 40,
               "normal_cells": 20, "wake_cells": 20, "wall_spacing": 2e-3},
    "output": "naca2412.json"}
   ```

   ```
   kofl call block.template --input template-airfoil.json
   ```

   The topology is checked before it is written (6 blocks, 2,400 cells).

2. Check or list it at any time (`{"topology": "naca2412.json"}`):

   ```
   kofl call block.check --input check.json
   kofl call block.list --input check.json
   ```

3. Build the grid. `grid.json`:

   ```json
   {"topology": "naca2412.json", "output": "naca2412-mesh",
    "algorithm": {"kind": "elliptic", "wall_sliding": true}}
   ```

   ```
   kofl call block.grid --input grid.json
   ```

   The directory `naca2412-mesh` holds `foehncraft.vtu` (read by
   FoehnCraft's `foehn mesh` and by MurmelCraft), `mesh.msh` (Gmsh MSH
   4.1), `volume.vtu`, `report.json` (per block: cells, Jacobians, aspect
   ratio; wall angles before and after smoothing) and `mesh-evidence.json`
   with the BLAKE3 hash of every file. A grid with a non-positive corner
   Jacobian, a gap or a volume error is refused, never written.

## Templates

| `template` | Shape | Main inputs |
|---|---|---|
| `channel` | one block, optional two-sided wall clustering | `length`, `width`, `height`, `cells` (three counts), `wall_spacing` |
| `pipe_ogrid` | H-core and four O-grid blocks | `radius`, `length`, `core_cells`, `radial_cells`, `axial_cells`, `wall_spacing` |
| `elbow` | the pipe O-grid swept through a bend, optional straight legs; walls projected onto cylinders and the torus | as `pipe_ogrid`, plus `bend_radius`, `angle`, `inlet_length`, `outlet_length` |
| `airfoil_cgrid` | extruded C-grid of six blocks | `naca` (4 digits) or `coordinates` (Selig or Lednicer file), `span`, `farfield`, `wake_length`, `surface_cells`, `normal_cells`, `wake_cells`, `wall_spacing` |
| `multi_element_airfoil` | extruded H-grid around a main element with an optional slat and flap; layout optimised, triangle blocks split where it helps | `main`, `slat`, `flap` (each `naca` or `coordinates`, with `chord`, `deflection`, `gap` and a negative `overlap`), and the C-grid inputs |
| `triangle`, `tetrahedron` | three or four hexahedra around the centroid | `corners`, `cells` per corner; the triangle also `span`, `span_cells` |

A three-element airfoil (40 blocks, 3,432 cells, elliptic grid in about 7 s):

```json
{"config": {"template": "multi_element_airfoil",
            "main": {"naca": "0012"},
            "slat": {"naca": "4410", "chord": 0.16, "deflection": 25, "gap": 0.05, "overlap": -0.04},
            "flap": {"naca": "4412", "chord": 0.3, "deflection": 30, "gap": 0.03, "overlap": -0.02},
            "span": 0.1, "farfield": 5, "wake_length": 5, "surface_cells": 24,
            "normal_cells": 12, "wake_cells": 12, "wall_spacing": 0.002},
 "output": "three.json"}
```

## Grid algorithms

- `{"kind": "tfi"}` (default): transfinite interpolation of faces and
  blocks from the edge curves and laws.
- `{"kind": "elliptic"}`: TFI, then Thompson–Thames–Mastin smoothing of
  every node off the boundary, across block faces, with orthogonality and
  first-spacing control at wall patches (`wall_orthogonality: false`
  switches that off). Nodes on planar `symmetry` patches slide in their
  plane.
- `"wall_sliding": true` lets wall nodes slide along curved walls
  (section contours, cylinders, the torus) toward orthogonal grid lines;
  corners sharper than 30° stay fixed.

The acceptance run (`cargo m1 --require-foehn`) builds the pipe, the elbow
and the airfoil with all three and fails if smoothing raises FoehnCraft's
largest non-orthogonality above the TFI grid's.

## Editing a topology

Each edit tool reads a topology and writes a new, checked one, or refuses
the edit and names the blocks, edge or face concerned:

| Tool | Edit |
|---|---|
| `block.set_count` | cell count of the chain of opposite edges through an edge |
| `block.set_distribution` | spacing law of an edge: `uniform`, `geometric`, `tanh` |
| `block.set_edge_curve` | `line`, `arc` through a point, `spline` through points |
| `block.move_vertex` | move a vertex `to` a point or `by` an offset |
| `block.add_vertex`, `block.remove_vertex` | vertices |
| `block.add_block`, `block.remove_block` | blocks, with the patches of their faces |
| `block.split_block` | split a block (and its chain) along an axis after `at` cells |
| `block.split_triangle` | turn a block with a flat corner into three around a centre |
| `block.optimize_layout` | move inner vertices and bend inner edges toward orthogonal cells |
| `block.set_patch` | rename a patch or set its role |

For example, split the core block of a pipe halfway along its axis:

```json
{"topology": "pipe.json", "block": 0, "axis": 2, "at": 20, "output": "pipe-split.json"}
```

```
kofl call block.split_block --input split.json
```

## T-junctions

A block corner may sit inside an edge of a bigger block: a T-junction
(BS-21). The topology lists such edges under `junctions`, each with the
vertices it passes through in order:

```json
"junctions": [{"edge": [1, 2], "through": [8]}]
```

The big block keeps its edge 1–2 whole; on the other side, 1–8 and 8–2 are
edges of two smaller blocks. The big edge's cell count is the sum of the
pieces' counts, and the face of the big block is gridded first, so the
smaller blocks' faces take its nodes: the cells still meet node for node,
and the mesh is plain hexahedra for every solver. `block.check` lists the
junction edges and the coarse faces they make (`t_junctions`).

`block.split_triangle` with `"stop": "t_junction"` uses them: the new row of
cells stops at the neighbouring blocks instead of cutting every block on its
way to the boundary. The neighbours keep their blocks and carry the new cells
in their counts; a count solve picks one integer per chain, close to the
counts before the split, so that every sum adds up. The result reports what
the solve did (`counts`). In a conforming mesh the new cells still run along
the chain to the boundary; what stops at the T-junction are the block cuts.

Faces must be covered in strips (the cuts run from one side of the big face
to the other). Elliptic smoothing and the edits that add, remove or split
blocks or set chain counts refuse a topology with T-junctions for now; use
the TFI grid.

## Limits

- Hexahedral blocks only; one layer of blocks for the 2D-extruded tools
  (`split_triangle`, `optimize_layout`).
- The block mesher writes `foehncraft.vtu`, not yet the blocked
  `foehncraft.vtkhdf` of the polyhedral bundles.
- Joining a block region to a polyhedral one (MD-23) is decided but not
  built.
