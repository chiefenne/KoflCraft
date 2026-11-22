# Nonconvex STL volumes and exposed prism layers

The SolidCraft and rectangular elbow fixtures support a distributed polyhedral core with wall prisms exposed at both ports. Each fixture includes a Gmsh volume scaffold and a partition file with its SHA-256. Rust builds one median-dual polyhedron per scaffold vertex. The `bulk` and `layered` modes use the original STL and a partition into star-shaped regions.

## Generate and inspect

From the repository root, build the CLI and choose a new output directory for each run.

```sh
cargo build --locked --release -p mesh-cli
case_dir=test-fixtures/solidcraft-vs-foehncraft
target/release/mesh-cli generate-duct "$case_dir/SolidCraft_vs_FoehnCraft.stl" "$case_dir/patches.json" "$case_dir/partition.json" bulk target/my-solidcraft-bulk
target/release/mesh-cli generate-duct "$case_dir/SolidCraft_vs_FoehnCraft.stl" "$case_dir/patches.json" "$case_dir/partition.json" layered target/my-solidcraft-layered
target/release/mesh-cli generate-duct "$case_dir/SolidCraft_vs_FoehnCraft.stl" "$case_dir/patches.json" "$case_dir/partition.json" polyhedral target/my-solidcraft-polyhedral

case_dir=test-fixtures/elbow
target/release/mesh-cli generate-duct "$case_dir/elbow.stl" "$case_dir/patches.json" "$case_dir/partition.json" layered target/my-elbow
```

Open `mesh.vtm` in ParaView for the volume and three named boundary blocks. Open `inlet.vtu` or `outlet.vtu` separately, choose **Surface With Edges**, and color by **LayerIndex**. Zero marks a core face. Positive indices count the prism layers inward from the wall. `OwnerCellId` gives the adjacent cell in `volume.vtu`. For an interior slice, open `volume.vtu` and apply **Slice**. The elbow midplane is at `z = 0.5`.

In `polyhedral` mode, `volume.vtu` contains VTK polyhedra (type 42) and six-node prisms (type 13). To see the core faces in ParaView, slice the volume, clear **Triangulate the Slice**, and choose **Surface With Edges**. A horizontal slice at `z = 3.7` cuts through the SolidCraft interior. The wall surface shows the ends of the wall prisms. A second horizontal slice at `z = 40` passes through the outlet pipe.

The Gmsh file `mesh.msh` preserves wall prisms and subdivides general polyhedra into native pyramids and tetrahedra. In polyhedral mode, `ParentCellId` maps each volume element back to its zero-based cell ID in `volume.vtu`. Physical groups 1, 2, 3 and 4 are inlet, outlet, walls and fluid. Each volume element's elementary entity is its one-based partition region. The JSON report records source and configuration digests, cell counts, measured volume, patch coverage and per-layer wall-normal thickness ranges.

For MCP, place the STL, both JSON files and the configured `volume-scaffold.msh` in the server's authorized directory. Call `generate_duct` with their filenames (`surface`, `patches`, `partition`), `mode` (`bulk`, `layered` or `polyhedral`), and a new output directory name (`output`). The same Rust engine serves the CLI and MCP. The MCP job limit is 100,000 source facets, four layers and 64 regions.

## Fixture results

| Case | Core tetrahedra | Core polyhedra | Wall prisms | Fluid volume | Inlet / outlet area |
|---|---:|---:|---:|---:|---:|
| 90° elbow, bulk | 52 | 0 | 0 | 5 | 1 / 1 |
| 90° elbow, layered | 52 | 0 | 120 | 5 | 1 / 1 |
| 90° elbow, polyhedral | 0 | 227 | 6,768 | 5 | 1 / 1 |
| SolidCraft, bulk | 93,326 | 0 | 0 | 5,290,738.125459 | 1,963.240067 / 1,963.240102 |
| SolidCraft, layered | 93,326 | 0 | 269,616 | 5,290,738.125459 | 1,963.240067 / 1,963.240102 |
| SolidCraft, polyhedral | 0 | 26,962 | 390,492 | 5,281,069.900754 | 1,963.240067 / 1,963.240102 |

The elbow consists of five unit cubes arranged in an L. Its two external openings are perpendicular. The SolidCraft partition uses six regions spanning the inlet pipe, bend, straight connector, box, rounded outlet junction and outlet pipe. Inlet and outlet roles follow the supplied fixture convention. [Elbow partition](../test-fixtures/elbow/partition.json) · [SolidCraft partition](../test-fixtures/solidcraft-vs-foehncraft/partition.json) · [SolidCraft patch map](../test-fixtures/solidcraft-vs-foehncraft/patches.json)

SolidCraft uses layer increments of 0.08, 0.112 and 0.1568 model units. Its polyhedral mesh retains the original two 240-edge port outlines and both port areas to roundoff. The wall is remeshed along the curved pipes, fillets and chamber. The fitted fluid volume is 0.183% smaller than the STL volume. In a two-way check of wall vertices and triangle centroids, the largest sampled point-to-surface distance was 0.516 model units. Wall-normal advances in the first layer range from 0.07966 to 0.08976 units. All 130,164 innermost prism faces match a core face by vertex ID and opposite winding.

The previous radial core spanned up to 255.6 units in a single cell, measured by the diagonal of its axis-aligned bounding box. The distributed core's maximum is 24.2 units. Its median sphericity is 0.833, its fifth percentile is 0.695 and its minimum is 0.344. Sphericity equals the surface area of a sphere of equal volume divided by the actual cell area. The cells contain 20 to 144 planar faces because each shared surface between neighbouring cells consists of several quadrilaterals. Flow accuracy remains untested.

The original `layered` mode retains every source wall triangle. It exposes 240 quadrilaterals per layer at each port, and measured normal advances stay within 0.05% of the requested thicknesses. The six radial regions produce 93,326 core tetrahedra, with 91,868 below a mean ratio of 0.1.

## Scaffold for a new shape

Supply a Gmsh MSH 2.2 ASCII file containing only linear tetrahedra (physical volume 4) and outward-oriented boundary triangles (physical surfaces 1 inlet, 2 outlet and 3 walls). The size and distribution of the scaffold vertices control polyhedral spacing. Reference it in the partition JSON with `volume_scaffold.file`, `volume_scaffold.sha256` and `volume_scaffold.max_relative_geometry_error`. The CLI resolves the filename relative to the partition. The geometry error limit applies separately to total fluid volume and each patch area. Port planes and orientations are checked against the source STL. Existing configurations without a scaffold retain the regional fan construction under the algorithm name `stl_partitioned_dual_v1`.

The supplied scaffold is sufficient to run the Rust CLI. To prepare another SolidCraft resolution, edit the size field in the Gmsh geometry file and run `python3 test-fixtures/solidcraft-vs-foehncraft/prepare_scaffold.py /path/to/new-scaffold.msh`. Gmsh optimization can change node positions and counts between preparations. The checked-in scaffold pins the input for reproducible Rust output. Use the new file and its printed digest in the partition. The preparation script groups smooth half-pipes and planar chamber faces into surface charts and preserves the original port rim segments. [SolidCraft fixture preparation](../test-fixtures/solidcraft-vs-foehncraft/prepare_scaffold.py) · [SolidCraft sizing](../test-fixtures/solidcraft-vs-foehncraft/prepare-scaffold.geo) · [Elbow sizing](../test-fixtures/elbow/prepare-scaffold.geo)

## Construction and checks

- The importer verifies the source hash against both maps, reads the closed shell and checks that the facet ranges partition it exactly once. The backend detects intersections using a bounding-volume tree and triangle contact tests. Adjoining source regions must share one simple planar loop whose center sees all rim edges. Both regions use the same triangulated interface disk. The specified region point must have strictly positive clearance from every core triangle. A fan from that point produces the core tetrahedra. [Polyhedron kernel computation](https://arxiv.org/abs/2202.06625)
- With a volume scaffold, Rust offsets the wall and extends its displacement through interior points by inverse-squared edge-length weighting. Every deformed tetrahedron must retain positive orientation. Edge midpoints, face centroids and tetrahedron centroids define a quadrilateral between each pair of adjacent vertex cells inside each tetrahedron. The quadrilateral lies in the plane where the two barycentric coordinates are equal. The checker compares every cell's signed shell volume with one quarter of the sum of its incident tetrahedron volumes. Each coarse wall prism splits into six subprisms whose end triangles match the dual shell exactly. [Volume dual construction](../crates/mesh-algorithms/src/duct/volume_dual.rs)
- Wall nodes advance along directions derived from angle-weighted incident triangle normals. At sharp edges, the algorithm offsets the distinct wall planes together. Each planar port contracts about its rim center in its own plane. This preserves its original triangulation. The strips between adjacent port rims carry the same physical patch as the core opening.
- Every layer front is checked for self-intersection. For each prism the Jacobian determinant is affine across a triangle at fixed height and quadratic along each vertical edge. Checking the three edge minima over the full height detects interior inversions. Integrating these quadratics gives the prism volume. Strict mode rejects collisions, invalid kernel clearance and nonpositive elements. [Dyedov et al., prismatic boundary layers and Jacobian validity](https://pmc.ncbi.nlm.nih.gov/articles/PMC2745959/)
- The checker independently pairs opposite internal face cycles, matches the owner and winding of every boundary face, and verifies volumes and all three patch areas against the meshed shell to relative tolerance `1e-9`. For a remeshed scaffold, the report gives both the original and fitted volume and patch areas. The fixture regression also checks achieved normal thickness and that every layer reaches both ports. Run `cargo test --release --locked -p mesh-cli nonconvex_elbow_and_solidcraft`. [Mesher and checks](../crates/mesh-algorithms/src/duct/mod.rs)
