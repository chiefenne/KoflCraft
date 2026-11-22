# Block and straight-pipe development meshes

These are the primary cases for core spacing, wall-layer growth and open-patch inspection. Coordinates and lengths are in millimetres. Rust generates both meshes directly from their analytic dimensions.

| Case | Geometry and patches | Core spacing | Wall layers, starting at the wall |
|---|---|---:|---|
| Block | 100 × 100 × 10 mm. The 100 × 100 mm face at z = 0 is the inlet. The other five faces are walls. | 1.00 mm | Five, starting at 0.38 mm, growing by 1.20 |
| Pipe | Diameter 50 mm, length 150 mm, centred on the x-axis. The inlet is at x = 0 and the outlet at x = 150 mm. | 1.25 mm | Five, starting at 0.48 mm, growing by 1.20 |

## Generate and inspect

From the repository root, choose a new output directory for each run.

```sh
cargo build --locked --release -p mesh-cli
mkdir -p target/development-meshes
target/release/mesh-cli generate-sweep test-fixtures/analytic/block.json target/development-meshes/block
target/release/mesh-cli generate-sweep test-fixtures/analytic/pipe.json target/development-meshes/pipe
```

Open `inlet.vtu` in ParaView, select **Surface With Edges**, and color by **LayerIndex**. Zero identifies the core. Layers 1 through 5 count inward from the wall. The inlet is a separately tagged section needed to close the volume mesh. Open `walls.vtu` alone for an unobstructed view into the opening. Open `mesh.vtm` for the volume and named surface patches. In the block, a longitudinal slice also shows the layers growing off the blind end.

For polygonal sections, load `volume.vtu`, apply **Slice**, and clear **Triangulate the Slice**. To verify the exported data and save both sections automatically, run the following with ParaView's Python interpreter.

```sh
pvpython test-fixtures/analytic/check_sections.py target/development-meshes/block target/development-meshes/pipe
```

The script measures cell volumes and adjacent-cell size changes from VTU connectivity, integrates the VTK polyhedron face streams, and checks patch ownership, opening areas and exposed layers. It writes `cross-section.vtu` and `longitudinal-section.vtu`. The Gmsh export is `mesh.msh`. It preserves the fitted geometry and patches using wedges and hexahedra. Its `ParentCellId` field maps each volume element to its zero-based VTU cell.

## Mesh construction and sizing

The core consists of convex polygonal prisms. The generator clips a two-dimensional Voronoi tessellation to the inset cross-section and sweeps it along the opening normal. The block uses staggered site rows. The pipe uses circumferential rings with decreasing site counts towards the axis. Each wall prism shares its complete interface with its neighbour. At the blind end the axial intervals follow the same wall-layer thickness schedule. At an intersection of two wall stacks the smaller layer index identifies the nearest wall layer.

`core_spacing_mm` sets the nominal site spacing and maximum axial core step. Change it to refine the same physical geometry. To refine the complete mesh approximately proportionally, scale `first_layer_mm` by the same factor. The block's `inlet` can be set to `x_min`, `y_min` or `z_min`. For MCP, pass the same configuration to `generate_sweep` with a new output name.

The volume VTU contains the following per-cell fields.

| Field | Meaning |
|---|---|
| `CellSizeMm` | Equivalent edge length, the cube root of measured cell volume |
| `MaxNeighbourSizeRatio` | Largest ratio of the larger to the smaller equivalent edge length across any internal face of the cell |
| `MaxInternalNonOrthogonalityDeg` | Largest angle between an internal face normal and the line joining the two volume centroids |
| `MaxInternalSkewness` | Largest distance from an internal face centre to the centroid-line intersection, divided by the distance between volume centroids |

Both fixtures limit neighbouring equivalent sizes to a ratio of 1.15. Generation fails before export if the measured maximum exceeds the configured limit. This limit measures volume change. The separately configured factor of 1.20 controls the physical wall-normal thickness growth.

| Verified result | Block | Pipe |
|---|---:|---:|
| Core cells / all volume cells | 82,840 / 161,135 | 126,480 / 191,280 |
| Median equivalent core size | 0.917 mm | 1.197 mm |
| Largest neighbour size ratio | 1.110 | 1.055 |
| Axial intervals | 8 core + 5 blind-end layers | 120 core |

The block fills its exact analytic volume. The pipe uses 108 planar wall facets around the circumference, with vertices on the specified radius. Its fluid volume is 0.0564% below the ideal cylinder. These measurements establish mesh geometry and grading. Solver convergence requires a separate study.
