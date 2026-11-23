# Nonconvex STL volumes and exposed prism layers

The Manifold and rectangular elbow fixtures support a distributed polyhedral core with wall prisms exposed at both ports. The elbow includes a Gmsh volume scaffold referenced with its SHA-256 by its partition file. The Manifold fixture (formerly Manifold, `manifold`) builds its scaffold in Rust from the STL (`volume_scaffold.build`, see [Scaffold built in Rust](#scaffold-built-in-rust)). Rust builds the polyhedral core as a Voronoi diagram restricted to the scaffold and relaxed by Lloyd iterations. One cell per scaffold vertex is the start. Setting `"core": "median_dual"` in `volume_scaffold` gives the earlier median-dual core instead. The `bulk` and `layered` modes use the original STL and a partition into star-shaped regions.

## Generate and inspect

From the repository root, build the CLI and choose a new output directory for each run.

```sh
cargo build --locked --release -p mesh-cli
case_dir=test-fixtures/manifold
target/release/mesh-cli generate-duct "$case_dir/manifold.stl" "$case_dir/patches.json" "$case_dir/partition.json" bulk target/my-manifold-bulk
target/release/mesh-cli generate-duct "$case_dir/manifold.stl" "$case_dir/patches.json" "$case_dir/partition.json" layered target/my-manifold-layered
target/release/mesh-cli generate-duct "$case_dir/manifold.stl" "$case_dir/patches.json" "$case_dir/partition.json" polyhedral target/my-manifold-polyhedral

case_dir=test-fixtures/elbow
target/release/mesh-cli generate-duct "$case_dir/elbow.stl" "$case_dir/patches.json" "$case_dir/partition.json" layered target/my-elbow
```

Open `mesh.vtm` in ParaView for the volume and three named boundary blocks. Open `inlet.vtu` or `outlet.vtu` separately, choose **Surface With Edges**, and color by **LayerIndex**. Zero marks a core face. Positive indices count the prism layers inward from the wall. `OwnerCellId` gives the adjacent cell in `volume.vtu`. For an interior slice, open `volume.vtu` and apply **Slice**. The elbow midplane is at `z = 0.5`.

In `polyhedral` mode, `volume.vtu` contains VTK polyhedra (type 42) and six-node prisms (type 13). To see the core faces in ParaView, slice the volume, clear **Triangulate the Slice**, and choose **Surface With Edges**. A horizontal slice at `z = 3.7` cuts through the Manifold interior. The wall surface shows the ends of the wall prisms. A second horizontal slice at `z = 40` passes through the outlet pipe.

The Gmsh file `mesh.msh` preserves wall prisms and subdivides general polyhedra into native pyramids and tetrahedra. In polyhedral mode, `ParentCellId` maps each volume element back to its zero-based cell ID in `volume.vtu`. Physical groups 1, 2, 3 and 4 are inlet, outlet, walls and fluid. Each volume element's elementary entity is its one-based partition region. The JSON report records source and configuration digests, cell counts, measured volume, patch coverage and per-layer wall-normal thickness ranges.

For MCP, place the STL, both JSON files and, if the partition names one, the `volume-scaffold.msh` in the server's authorized directory. Call `generate_duct` with their filenames (`surface`, `patches`, `partition`), `mode` (`bulk`, `layered` or `polyhedral`), and a new output directory name (`output`). The same Rust engine serves the CLI and MCP. The MCP job limit is 100,000 source facets, four layers and 64 regions.

## From a bare STL

A closed STL with two pipe ends needs no hand-written port file or partition. `patch.detect_ports` finds the planar end caps and writes a patch map with `walls`, `inlet` and `outlet`; `mesh.setup_duct` writes a polyhedral partition sized from the smaller port; `generate-duct` reads the patch map in place of a port file, so the cap facets may lie anywhere in the STL.

```sh
cargo build --locked --release -p mesh-cli
case_dir=test-fixtures/cooling-jacket
mkdir -p target/jacket && cp "$case_dir/cooling-jacket.stl" target/jacket/
cd target/jacket
echo '{"surface": "cooling-jacket.stl", "output": "ports.json"}' > detect.json
../release/kofl call patch.detect_ports --input detect.json
echo '{"surface": "cooling-jacket.stl", "patches": "ports.json", "output": "partition.json", "size_factor": 2}' > setup.json
../release/kofl call mesh.setup_duct --input setup.json
../release/kofl generate-duct cooling-jacket.stl ports.json partition.json polyhedral mesh
```

**Detection.** The surface is split into coplanar regions by a flood fill from each seed facet to neighbours whose normal is within `planar_angle_deg` (default 1°) of the seed's and whose corners lie on the seed's plane (to 1e-5 of the longest extent); comparing with the seed keeps a finely tessellated cylinder from creeping into one region. A region whose rim is one connected loop (a disc) is a candidate. Its score is the isoperimetric quotient 4πA/P² (1 for a circle, 0.785 for a square) times the share of rim edges where the wall turns away from the cap at 60° to 120°, as at the flat end of a pipe. The result lists the best `candidates` (default 8) with facets, area, centre, normal, equivalent radius, circularity, rim angles and score, and `margin`, the lead of the weaker picked cap over the third candidate. The two best are inlet and outlet; which is which is a guess, so check the centres and pass `inlet` and `outlet` (candidate ids) to swap them or pick others. Any `patch.*` tool can edit the map afterwards; `generate-duct` takes a map with exactly one patch of role inlet, one of role outlet and every other patch a wall, and names them `inlet`, `outlet` and `walls`.

**Setup.** Every length is a fixed fraction of the smaller port's equivalent radius r, the fractions of the Manifold fixture (r = 25): chord tolerance r/500, wall sizes r/38.5 to r/4.17, rim size r/29.4, spacing rate 0.4, first layer r/312.5 with growth 1.4 and three layers, a wall deviation limit of 1.5 times the tolerance and a geometry error limit of 0.2%. `size_factor` multiplies the lengths (2 is twice as coarse); `layers`, `first_thickness` and `growth` override the layers. The partition has no regions, so it serves the polyhedral mode only. Edit it by hand for finer control.

On the synthetic [cooling jacket](../test-fixtures/cooling-jacket/README.md) (24,732 facets) both caps score 1.000 and the next candidate, a block face, 0.57. On the Manifold surface the two best candidates are exactly the port caps of its hand-made port file. Run `cargo test --release --locked -p mesh-engine --test ports`.

## Fixture results

| Case | Core tetrahedra | Core polyhedra | Wall prisms | Fluid volume | Inlet / outlet area |
|---|---:|---:|---:|---:|---:|
| 90° elbow, bulk | 52 | 0 | 0 | 5 | 1 / 1 |
| 90° elbow, layered | 52 | 0 | 120 | 5 | 1 / 1 |
| 90° elbow, polyhedral | 0 | 1,004 | 0 | 5 | 1 / 1 |
| Manifold, bulk | 93,326 | 0 | 0 | 5,290,738.125459 | 1,963.240067 / 1,963.240102 |
| Manifold, layered | 93,326 | 0 | 269,616 | 5,290,738.125459 | 1,963.240067 / 1,963.240102 |
| Manifold, polyhedral | 0 | 108,025 | 131,322 | 5,289,262.968278 | 1,963.240067 / 1,963.240102 |

The elbow consists of five unit cubes arranged in an L. Its two external openings are perpendicular. The Manifold partition uses six regions spanning the inlet pipe, bend, straight connector, box, rounded outlet junction and outlet pipe. Inlet and outlet roles follow the supplied fixture convention. [Elbow partition](../test-fixtures/elbow/partition.json) · [Manifold partition](../test-fixtures/manifold/partition.json) · [Manifold patch map](../test-fixtures/manifold/patches.json)

Manifold uses layer increments of 0.08, 0.112 and 0.1568 model units. Its polyhedral mesh retains the original two 240-edge port outlines and both port areas to roundoff. Its wall is remeshed onto the STL to a chord tolerance of 0.05 model units (56,606 wall triangles): the largest sampled distance between the wall and the STL is 0.050 in both directions, reported as `scaffold.wall_deviation` and bounded by `max_wall_deviation`. The earlier Gmsh scaffold was 0.54 off ([note](../research/m4-rust-scaffold.md)). The fitted fluid volume is 0.028% smaller than the STL volume (Gmsh scaffold: 0.183%). Wall-normal advances in the first layer range from 0.0800 to 0.0836 units. Manifold joins its layers with wall columns, the default. On its curved walls and at its reentrant edges the core has one conforming seed per wall triangle (46,571 seeds for 56,606 wall triangles), so the core cell over such a triangle meets the wall in exactly that triangle, and 43,774 triangles keep their native prisms. Flat walls and sharp convex corners have no seeds; columns below 5% of their wall triangle join a neighbour in their plane ([note](../research/m4-wall-columns.md)). FoehnCraft reads the mesh as valid: 239,347 cells, no inverted face, maximum non-orthogonality 74.5° with 3 faces above 70°, skewness 1.44, largest volume ratio 1,882, smallest cell 1.1e-3. With `"wall_join": "split_cap"` it has 224,028 cells, 89.1°, 1,931 faces above 70° and volume ratio 285. A build takes about 2.5 minutes.

In polyhedral mode the polyhedra count includes the wall-layer columns under the core cells (795 layer cells on the elbow, 34,596 on Manifold). The core itself has 209 cells on the elbow and 73,429 on Manifold. The previous radial core spanned up to 255.6 units in a single cell, measured by the diagonal of its axis-aligned bounding box. The distributed core stays below 25 units. On Manifold the Voronoi core's median sphericity is 0.744, its fifth percentile is 0.644 and its minimum is 0.353 (median dual: 0.833, 0.695 and 0.344); two thirds of the core cells are the low cells of the wall seeds. Sphericity equals the surface area of a sphere of equal volume divided by the actual cell area. FoehnCraft measures a maximum non-orthogonality of 56.6° on the elbow and 74.5° on Manifold, with 0 and 3 faces above 70° (split cap: 85.9° and 89.1°, 44 and 1,931; median dual: 87.0° and 89.6°). [Design note and benchmark](../research/m4-voronoi-core.md) · [Wall columns](../research/m4-wall-columns.md) Flow accuracy remains untested.

The original `layered` mode retains every source wall triangle. It exposes 240 quadrilaterals per layer at each port, and measured normal advances stay within 0.05% of the requested thicknesses. The six radial regions produce 93,326 core tetrahedra, with 91,868 below a mean ratio of 0.1.

## Scaffold built in Rust

With `volume_scaffold.build` the scaffold is made from the STL and needs no Gmsh. The wall facets (not the port caps, which are kept as they are) are remeshed onto the STL by split, collapse, flip, tangential smoothing and projection (Botsch & Kobbelt 2004), with target lengths from the curvature so that a facet stays within `tolerance` of the STL (Dunyach et al. 2013), between `min_size` and `max_size`, `rim_size` at the port rims and graded by `rate`. Edges sharper than `feature_angle` (default 30°) stay as feature lines. Facets still farther than the tolerance are split until they are not. The volume is then a Delaunay tetrahedralization with exact predicates, the wall recovered by midpoint splits on it, refined at circumcentres to the spacing grown from the wall, and cleaned of sliver tetrahedra by edge removal.

```json
"volume_scaffold": {"build": {"tolerance": 0.05, "min_size": 0.65, "max_size": 6, "rate": 0.4, "rim_size": 0.85},
                    "max_relative_geometry_error": 0.002, "max_wall_deviation": 0.06}
```

The report's `scaffold` object names the source (`rust_remesh_delaunay_v1` or `file`), the build statistics and the measured `wall_deviation` both ways; `max_wall_deviation` fails the run when either exceeds it, for a file scaffold too. [Design note and benchmark](../research/m4-rust-scaffold.md) · [Remesher](../crates/mesh-algorithms/src/remesh.rs) · [Tetrahedralization](../crates/mesh-algorithms/src/scaffold.rs)

## Scaffold from a file

Supply a Gmsh MSH 2.2 ASCII file containing only linear tetrahedra (physical volume 4) and outward-oriented boundary triangles (physical surfaces 1 inlet, 2 outlet and 3 walls). The size and distribution of the scaffold vertices control polyhedral spacing. Reference it in the partition JSON with `volume_scaffold.file`, `volume_scaffold.sha256` and `volume_scaffold.max_relative_geometry_error`. The CLI resolves the filename relative to the partition. The geometry error limit applies separately to total fluid volume and each patch area. Port planes and orientations are checked against the source STL. Existing configurations without a scaffold retain the regional fan construction under the algorithm name `stl_partitioned_dual_v1`. Gmsh optimization can change node positions and counts between preparations; the checked-in elbow scaffold pins the input for reproducible output. [Elbow sizing](../test-fixtures/elbow/prepare-scaffold.geo)

## Construction and checks

- The importer verifies the source hash against both maps, reads the closed shell and checks that the facet ranges partition it exactly once. The backend detects intersections using a bounding-volume tree and triangle contact tests. Adjoining source regions must share one simple planar loop whose center sees all rim edges. Both regions use the same triangulated interface disk. The specified region point must have strictly positive clearance from every core triangle. A fan from that point produces the core tetrahedra. [Polyhedron kernel computation](https://arxiv.org/abs/2202.06625)
- With a volume scaffold, Rust offsets the wall and extends its displacement through interior points by inverse-squared edge-length weighting. Every deformed tetrahedron must retain positive orientation. Edge midpoints, face centroids and tetrahedron centroids define a quadrilateral between each pair of adjacent vertex cells inside each tetrahedron. The quadrilateral lies in the plane where the two barycentric coordinates are equal. The checker compares every cell's signed shell volume with one quarter of the sum of its incident tetrahedron volumes. At a reentrant wall edge a wall vertex's dual cell is L-shaped and its centroid falls behind the wall planes, so the faces to the wall prisms there would have S·d ≤ 0; such a cell is split along scaffold faces through the vertex into one part per wall plane (10 cells on the elbow, 3 on Manifold). [Note](../research/m4-reentrant-dual-split.md) Each coarse wall prism splits into six subprisms whose end triangles match the dual shell exactly. [Volume dual construction](../crates/mesh-algorithms/src/duct/volume_dual.rs)
- Voronoi core (`volume_clipped_voronoi_cvt_v1`, default): each scaffold tetrahedron is clipped by the bisector planes of its nearest sites, with the security-radius stop of ALG-04. The sites move ten times to the density-weighted centroids of their cells (Lloyd). The density follows the scaffold spacing. A site's pieces are joined into star-shaped cells. A cell that wraps a reentrant wall edge has its site replaced by a pair mirrored across the plane that halves the wall angle, and the diagram is rebuilt. The core cells cut the innermost wall surface into one polygon per core cell and wall triangle. Each polygon is extruded through the layers at fixed barycentric coordinates, so it splits the coarse wall wedges exactly. The polygons of one core cell over triangles that are coplanar on every layer surface merge into one column, so on a flat wall each core cell sits on one layer cell per layer. On curved walls each polygon stays its own column, because FoehnCraft fans a bent face about its node average. So that these polygons are whole wall triangles and not slivers, the wall triangles at edges where columns end (curved walls, reentrant edges; not sharp convex corners) each get a fixed site on the balls of their three corners (conforming wall seeds, after VoroCrust), the other sites keep out of the balls, and a column below 5% of its wall triangle joins a neighbour in its plane. Every layer cell is checked through its sub-wedges' Jacobians and volumes. `"wall_join": "split_cap"` keeps the earlier join, where the innermost wedge shared by several core cells became one polyhedron with a split wall face. [Wall columns](../crates/mesh-algorithms/src/duct/columns.rs) · [Note](../research/m4-wall-columns.md) The median-dual bullet above applies with `"core": "median_dual"`. [Voronoi core](../crates/mesh-algorithms/src/duct/voronoi_core.rs) · [Note](../research/m4-voronoi-core.md)
- Wall nodes advance along directions derived from angle-weighted incident triangle normals. At sharp edges, the algorithm offsets the distinct wall planes together. Each planar port contracts about its rim center in its own plane. This preserves its original triangulation. The strips between adjacent port rims carry the same physical patch as the core opening.
- Every layer front is checked for self-intersection. For each prism the Jacobian determinant is affine across a triangle at fixed height and quadratic along each vertical edge. Checking the three edge minima over the full height detects interior inversions. Integrating these quadratics gives the prism volume. Strict mode rejects collisions, invalid kernel clearance and nonpositive elements. [Dyedov et al., prismatic boundary layers and Jacobian validity](https://pmc.ncbi.nlm.nih.gov/articles/PMC2745959/)
- The checker independently pairs opposite internal face cycles, matches the owner and winding of every boundary face, and verifies volumes and all three patch areas against the meshed shell to relative tolerance `1e-9`. For a remeshed scaffold, the report gives both the original and fitted volume and patch areas. The fixture regression also checks achieved normal thickness and that every layer reaches both ports. Run `cargo test --release --locked -p mesh-cli nonconvex_ -- --include-ignored` (the Manifold case is ignored in the fast run). [Mesher and checks](../crates/mesh-algorithms/src/duct/mod.rs)
