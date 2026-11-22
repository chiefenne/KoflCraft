# Wall layers at inlet and outlet patches

[GEO-06 and BL-10](../docs/requirements.md) separate a face’s geometric role in closing the volume from its boundary condition. The inlet and outlet caps close the STL shell for meshing. They describe openings in the flow domain and have their own named boundary roles.

Wall prism layers start from the wetted walls and reach the port rim. At a planar opening, constrain the rim displacements to the port plane. Each layer then ends flush at the opening, with its end faces exposed to the flow boundary. The interior part of the opening is filled by the core cell boundary faces. Both groups of faces carry the port patch so the inlet or outlet condition applies across the entire opening. Dyedov et al. describe projecting displacements on inlet/outlet rims onto the plane of the constrained opening. [Dyedov et al., Section 5.3](https://www.ams.sunysb.edu/~jiao/papers/boundary_layer.pdf)

The straight-channel benchmark checks the following conditions.

- Open patches constrain the rim and seed no wall extrusion of their own.
- The terminal prism faces and core opening form a conforming planar tiling of the original cap. Every face is outward-wound with exactly one volume-cell incidence.
- All requested layers reach both openings. The channel generator rejects a requested thickness stack if opposing walls would leave a zero-width or inverted core.
- The sum of exposed prism-end and core-face areas matches the original port area within the mesh conservation tolerance. Volume cells occupy the prescribed domain.

## Run the straight-channel benchmark

From the repository root, choose a new output directory.

```sh
cargo run --locked -p mesh-cli -- generate-channel test-fixtures/straight-channel.json target/my-channel
```

Open `target/my-channel/inlet.vtu` in ParaView, select **Surface With Edges**, color by **LayerIndex**, and view along the X axis. Zero is the core. Layers 1, 2 and 3 run from the wall inward. Open `channel.vtm` for the volume and all three named boundary blocks, or `channel.msh` in Gmsh for the `inlet`, `outlet`, `walls` and `fluid` physical groups.

Edit a copy of the [channel configuration](../test-fixtures/straight-channel.json) to vary length (X), width (Y), height (Z), axial subdivision count, first-layer thickness, layer count and growth ratio. Each successive thickness equals the previous thickness times the growth ratio. Their sum must leave a positive core between opposing walls. The CLI returns measured areas for every layer on both ports, core areas, total volume and the smallest cell volume. The MCP tool `generate_channel` takes the same inline `config` and an `output` directory name. Its per-job limits are eight layers and 32 axial subdivisions.

## Mesh and reference values

On each flat wall the generator extrudes triangular panels in the wall-normal direction. Where two stacks meet in a right-angle corner, a graded grid of axial wedges joins their quadrilateral sides. Port faces are quadrilaterals at the ends of wall-normal stacks, triangles at the ends of corner wedges, and a quadrilateral over the core. One polyhedral core cell per axial subdivision shares eight triangular faces with the adjoining stacks. All interfaces use common vertex and face IDs.

The reference channel measures 4 × 2 × 1. Its three layers have thicknesses 0.04, 0.06 and 0.09. Each opening has the following measured areas in square model units.

| Region | Faces per port | Area per port |
|---|---:|---:|
| First layer | 44 | 0.2336 |
| Second layer | 28 | 0.3264 |
| Third layer | 12 | 0.4356 |
| Core | 1 | 1.0044 |
| Whole port | 85 | 2.0000 |

The [generator](../crates/mesh-algorithms/src/channel.rs) pairs opposing face loops. The engine validates the resulting mesh against the prescribed box and checks the measured area of each port ring against the difference between successive inset rectangles. For this fixture the 768 prisms and eight core cells have total volume 8. ParaView 6.2 read back the same cell types and volume. Every face in the three boundary blocks matched exactly one volume-cell face.

In VTU, `OwnerCellId` on each boundary face refers to the zero-based cell index in `volume.vtu`. `PatchId` values 1, 2 and 3 identify the inlet, outlet and walls in both formats. The Gmsh writer preserves each prism. It subdivides each core polyhedron into eight tetrahedra and two pyramids sharing its centre, retaining the original port quads and prism interfaces. The [Gmsh regression](../tests/channel_gmsh.rs) checks the signed volume of each element, per-parent volume, face incidence, patch areas and port winding after serialization. The exports follow the [VTK XML specification](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html) and the [Gmsh MSH specification](https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2-_0028Legacy_0029).

## SolidCraft fixture

The [SolidCraft surface](../test-fixtures/solidcraft-vs-foehncraft/README.md) has two circular ports with nonparallel normals and a connected wall surface with sharp edges and roundings. The [nonconvex duct benchmark](../docs/nonconvex-duct.md) now exercises this shell as a bulk mesh and with three prism layers, including exposed port faces and a collision rejection case.
