# Polyhedral and layered meshes

Part of the [Recipes](../recipes.md): short command sequences for the
duct mesher (`generate-duct`) and the tools around it. The modes and
partition keys are described in [nonconvex-duct.md](../nonconvex-duct.md);
pictures of complete cases are on the gallery page
[Polyhedral and layered meshes](../gallery/polyhedral.md). `kofl` is
`target/release/kofl`; every tool reads one JSON file, and outputs must
not exist yet (IO003: delete or rename the old one).

| # | Recipe |
|---|---|
| 30 | [Mesh a closed STL with two pipe ends](#bare-stl) |
| 57 | [Name walls and keep flush openings](#named-walls) |
| 60 | [Choose how layers meet a flat wall](#layer-ends) |
| 62 | [Propose patches and write a FoehnCraft skeleton](#propose) |
| 64 | [Refine where you say](#grading) |
| 70 | [Add a baffle or a plate](#shells) |

---

<a id="bare-stl"></a>
## 30. Mesh a closed STL with two pipe ends

```
echo '{"surface": "part.stl", "output": "ports.json"}' > detect.json
kofl call patch.detect_ports --input detect.json
echo '{"surface": "part.stl", "patches": "ports.json", "output": "partition.json"}' > setup.json
kofl call mesh.setup_duct --input setup.json
kofl generate-duct part.stl ports.json partition.json polyhedral mesh
```

Look at `result.detection.candidates` first: the two best caps become
inlet and outlet; pass `"inlet"` and `"outlet"` (candidate ids) to swap or
pick others. `"size_factor": 2` in `setup.json` makes the mesh twice as
coarse for a first look; `layers`, `first_thickness` and `growth` set the
layers. `kofl --check generate-duct ...` gives the estimate (cells,
memory, files) without building. Gallery [30](../gallery/polyhedral.md#bare-stl).

<a id="named-walls"></a>
## 57. Name walls and keep flush openings

A patch map may hold any number of wall patches, one inlet and one
outlet. Make one from the STL's solids with `patch.init`, then move facets
with `patch.select_plane`, `patch.select_box`, `patch.select_normal` or
`patch.grow`, and set roles with `patch.set_role`. An inlet or outlet
flush with its wall (a burner disc) grows layers like the wall. Every name
reaches FoehnCraft as its own boundary. Gallery
[57](../gallery/polyhedral.md#furnace).

<a id="layer-ends"></a>
## 60. Choose how layers meet a flat wall

In the partition, for a planar wall patch:

```
"no_layers": ["sides"]       // no layers there; the other walls' layers run into it
"layers_yield": ["sides"]    // its layers stay; crossing layers run through them
```

A patch can't be in both. Without either, the layers of two walls mitre
at their corner. Gallery [60](../gallery/polyhedral.md#furnace-variants).

<a id="propose"></a>
## 62. Propose patches and write a FoehnCraft skeleton

```
kofl call patch.propose --input propose.json    # {"surface": "part.stl", "output": "proposed.json"}
kofl call patch.rename --input rename.json      # {"surface": ..., "patches": "proposed.json", "from": "opening_1", "to": "inlet", "output": ...}
kofl call patch.set_role --input role.json      # {"surface": ..., "patches": ..., "name": "inlet", "role": "inlet", "output": ...}
kofl call bc.skeleton --input skeleton.json     # {"surface": ..., "patches": "patches.json", "output": "case.toml"}
```

`patch.propose` splits the STL at feature edges, says what each face is
(plane, cylinder, ...) and lists in `needs_role` the patches still
without a role. `case.toml` has one `[boundary.<patch>]` per patch with
its FoehnCraft type and a TODO for each value. Gallery
[70](../gallery/polyhedral.md#pipe-shells).

<a id="grading"></a>
## 64. Refine where you say

```json
{"partition": "partition.json", "output": "graded.json",
 "source": {"kind": "datum", "name": "probe", "size": 0.25, "distance": 1,
            "shape": {"type": "sphere", "centre": [0, 0, 32], "radius": 1}}}
```

```
kofl call grading.add --input probe.json     # grading.list, grading.remove by name
```

A `datum` (point, segment, plane, box, slab, sphere, cylinder, cone) holds
`size` within `distance` and grows away at the build's rate; a `region`
(box, slab, cylinder) holds `size` inside. The build report's `grading`
object gives the achieved cell sizes per source.

<a id="shells"></a>
## 70. Add a baffle or a plate

```json
{"partition": "partition.json", "output": "partition-1.json",
 "shell": {"name": "baffle", "polygon": [[-9, -9, 24.3], [0.7, -9, 24.3], [0.7, 9, 24.3], [-9, 9, 24.3]]}}
```

```
kofl call shell.add --input baffle.json      # shell.list, shell.remove by name
```

A planar convex polygon; the part inside the fluid is meshed, and it may
reach past the wall. Each shell becomes two wall patches, `<name>` and
`<name>_back`. Optional `size` and `edge_size` set its cell size and the
size along its free edges. Gallery [70](../gallery/polyhedral.md#pipe-shells).
