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

## Size control

Instead of `sites`, a config can give a `size_field` that places them ([example](../test-fixtures/stl-size-field.json)); give one or the other. The size h is the edge of the equivalent cube, cbrt(cell volume).

```sh
target/debug/kofl generate-stl test-fixtures/unit-box.stl test-fixtures/stl-size-field.json mesh-output/stl-graded.mesh.json
```

| Field | Meaning |
|---|---|
| `min_size`, `max_size` | Bounds of h; `max_size` is also the size where no source applies |
| `gradation` | Largest neighbour size ratio; the field's slope is capped at `gradation - 1` |
| `relaxation_iterations` | Lloyd (CVT) iterations, default 20 |
| `frames` | Named coordinate systems: `name`, `origin`, `x_axis`, `y_axis` (`z = x × y`); `global` is predefined |
| `sources` | Datum objects and coordinate regions, each with a unique `name` |

A **datum** (`"kind": "datum"`, GRD-02) has a `shape`, a `size`, an optional `distance` within which the size holds (default 0), an optional `growth` ratio beyond it (default and cap: the gradation) and an optional `frame`. Shapes: `point` (`point`), `segment` (`start`, `end`), `plane` (`point`, `normal`; two-sided), `box` (`min`, `max`), `slab` (`axis` x/y/z, `min`, `max`), `sphere` (`centre`, `radius`), `cylinder` (`start`, `end`, `radius`) and `cone` (`start`, `end`, `start_radius`, `end_radius`). Bodies are solid. A **region** (`"kind": "region"`, GRD-03) is a `box`, `slab` or `cylinder` in a frame with a `size` inside it.

The field is the minimum over all sources, clamped to `[min_size, max_size]`. The seed count is the volume integral of h⁻³; `kofl --check` and `mesh.estimate` report it exactly. Seeds are drawn from the field with the config's `seed` (default 0) and relaxed into a density-weighted centroidal Voronoi tessellation; the choice and its benchmark are in [research/m2-size-control.md](../research/m2-size-control.md). The artifact stores the placed sites, so checks and exports work as with explicit sites.

The result and the artifact carry a `sizing` block: the field as given, the seed count, the iterations, the size ratio cbrt(V)/h and the neighbour size ratio (quantiles), the interior faces whose neighbour ratio exceeds the gradation, and per source the target, the achieved sizes and the number of cells it governs (GRD-08). The gradation bounds the field, not each pair of cells: about one face in ten exceeds it on the fixture, and the result lists that as a limitation.

Limits: datum objects from STL surfaces, and curvature, gap and feature-edge sizes (GRD-04), are open. Clipping is all-pairs (ALG-04), so each relaxation iteration costs O(N²); a field may ask for at most 20,000 seeds, and the MCP site limit applies to the seed count.

## MCP

Launch the server with read and write roots as in the [M1 workflow](m1-usage.md). Place the STL below a read root before calling the tools.

| Tool | Arguments |
|---|---|
| `inspect_stl_domain` | `surface` filename and inline STL `config`; with a size field it reports the seed count as `site_count` |
| `generate_stl_mesh` | `surface` filename, inline STL `config`, and new `artifact` filename |

Reuse `check_mesh`, `mesh_report`, and `export_mesh` with the generated artifact. The MCP server limits this prototype to 2,048 STL facets and 256 sites per job. The CLI uses the same engine and accepts the STL path as an explicit argument.

## Verification

`cargo test --release --locked -p mesh-engine --test size_field` builds graded meshes on both fixtures and checks the achieved sizes, determinism, the estimate and the evidence digest. `cargo test --locked --test stl_workflow` compares each unit-box STL cell with the analytic reference, checks boundary patch areas, verifies binary STL with explicit labels, and exercises the skewed fixture plus malformed surfaces. `cargo m1 --require-gmsh` includes these tests in the workspace regression suite.

The Rust `import_closed_stl` entry point also checks shell connectivity, winding, and positive signed volume for nonconvex test surfaces. The [SolidCraft pipe-and-box case](../test-fixtures/solidcraft-vs-foehncraft/README.md) uses it to verify its inlet, outlet, and wall labels before convexity enters the volume meshing path.

The [partitioned nonconvex workflow](nonconvex-duct.md) accepts the elbow and SolidCraft fixtures and exports bulk or layered volume meshes.
