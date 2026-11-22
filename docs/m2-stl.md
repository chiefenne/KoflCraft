# M2 convex STL meshing

The STL workflow accepts a single connected, closed convex shell in ASCII or binary STL. It uses the facet vertex order to determine winding and welds vertices with exactly identical coordinates. The mesher uses the original boundary triangles as clipping faces, so patch divisions on the same plane are preserved.

## Generate a mesh

```sh
cargo build --locked -p mesh-cli -p mesh-mcp
mkdir -p mesh-output
target/debug/mesh-cli inspect-stl test-fixtures/unit-box.stl test-fixtures/stl-sites.json
target/debug/mesh-cli generate-stl test-fixtures/unit-box.stl test-fixtures/stl-sites.json mesh-output/stl-box.mesh.json
target/debug/mesh-cli check mesh-output/stl-box.mesh.json
target/debug/mesh-cli report mesh-output/stl-box.mesh.json
target/debug/mesh-cli export mesh-output/stl-box.mesh.json vtu mesh-output/stl-box.vtu
target/debug/mesh-cli export mesh-output/stl-box.mesh.json gmsh mesh-output/stl-box.msh
```

Use `test-fixtures/skewed-box.stl` with a new artifact name to exercise angled boundary faces. VTU retains the original polyhedral cells. Gmsh exports a conforming tetrahedral subdivision and preserves the named boundary patches as physical groups.

## Sites and patch names

The [STL configuration](../test-fixtures/stl-sites.json) has `schema_version: 1`, an ordered `sites` array, and optional `fixture_id` and `seed`. All sites must lie inside the actual surface. The inspection report gives surface bounds, enclosed volume, triangle and vertex counts, and area by named patch.

ASCII `solid` labels supply patch names. Binary facets and unnamed ASCII solids default to `surface`. To supply or override labels, set `facet_patches` to an array with one name per triangle in file order. Adjacent triangles can share a name or belong to different surface patches. Use nonempty names without quotes or control characters.

STL artifacts use schema version 2 and embed the indexed shell and patch table. The artifact's SHA-256 input and config digests cover this surface. Checks and exports revalidate the stored shell and mesh directly, so they can run after the source STL has moved. The box workflow uses schema version 1.

## MCP

Launch the server with the same `--output-dir` argument as the [M1 workflow](m1-usage.md). Place the STL in that directory before calling the tools.

| Tool | Arguments |
|---|---|
| `inspect_stl_domain` | `surface` filename and inline STL `config` |
| `generate_stl_mesh` | `surface` filename, inline STL `config`, and new `artifact` filename |

Reuse `check_mesh`, `mesh_report`, and `export_mesh` with the generated artifact. The MCP server limits this prototype to 2,048 STL facets and 256 sites per job. The CLI uses the same engine and accepts the STL path as an explicit argument.

## Verification

`cargo test --locked --test stl_workflow` compares each unit-box STL cell with the analytic reference, checks boundary patch areas, verifies binary STL with explicit labels, and exercises the skewed fixture plus malformed surfaces. `cargo m1 --require-gmsh` includes these tests in the workspace regression suite.

The Rust `import_closed_stl` entry point also checks shell connectivity, winding, and positive signed volume for nonconvex test surfaces. The [SolidCraft pipe-and-box case](../test-fixtures/solidcraft-vs-foehncraft/README.md) uses it to verify its inlet, outlet, and wall labels before convexity enters the volume meshing path.

The [partitioned nonconvex workflow](nonconvex-duct.md) accepts the elbow and SolidCraft fixtures and exports bulk or layered volume meshes.
