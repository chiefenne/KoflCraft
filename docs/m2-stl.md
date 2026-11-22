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

Limits: datum objects from STL surfaces, and curvature, gap and feature-edge sizes (GRD-04), are open. Each cell is clipped by its nearest sites only (ALG-04, [note](../research/alg-04-neighbour-clipping.md)), so a relaxation iteration costs about O(N log N); a field may ask for at most 200,000 seeds, and the MCP site limit applies to the seed count.

## Boundary patches

The `patch.*` tools (BC-01, BC-02) name the boundary faces after import, so an agent or a person can define the patches of a FoehnCraft case without editing the STL. They work on any closed shell, convex or not. Each call reads the STL and a patch map, applies one operation and writes a **new** patch map (`koflcraft.patch-map` 1.0): the STL's SHA-256, the facet count, the patches with `name`, `role` and an optional `source_name`, and `facet_patch`, the patch index of every facet in file order. Every read checks that the map belongs to the STL and gives every facet exactly one patch with a unique name (BC-06). Names are FoehnCraft boundary names (1 to 64 ASCII letters, digits, `_`, `-`), so MSH 4.1, the FoehnCraft VTU and OpenFOAM carry the same names.

| Tool | What it does |
|---|---|
| `patch.init` | A map from the STL solids: one patch per solid (`west inlet` becomes `west_inlet`, binary STL gives `surface`), roles from the names |
| `patch.list` | The check and the evidence of each patch: role and FoehnCraft type, facets, area, mean normal and its largest deviation, bounds, regions, neighbouring patches |
| `patch.select_plane` | Facets with all corners on a plane (`point`, `normal`, optional `tolerance`) |
| `patch.select_box` | Facets in a box (`min`, `max`; by centroid, or `contain: all_corners`) |
| `patch.select_normal` | Facets whose outward normal is within `max_angle_deg` of `direction` |
| `patch.grow` | Feature-angle region growing from `seed_facet` or the facet nearest to `seed_point`, across edges with a dihedral angle up to `feature_angle_deg` |
| `patch.split` | A patch into its regions (separated by edges sharper than `feature_angle_deg`, or by connectivity), named `names` or `<patch>_1`, … |
| `patch.merge` | The patches in `merge` into `name` |
| `patch.rename`, `patch.set_role` | Rename (the new name must be free); set the role: wall, inlet, outlet, symmetry, interface or unspecified ([mapping](foehncraft.md#patch-names-and-roles-fc-04)) |

A selection moves the facets it picks into the patch `name` (created when new, with `role` or the role its name suggests), from any patch or only from the patches in `within`. Patches it empties are dropped; the result lists `change.created`, `change.removed` and `change.facets_changed` with the new summary. An operation that selects or changes nothing writes no file and fails with PATCH004. A rename keeps the role.

```sh
cat > init.json <<'JSON'
{"surface": "test-fixtures/unit-box.stl", "output": "mesh-output/p0.json"}
JSON
cat > inlet.json <<'JSON'
{"surface": "test-fixtures/unit-box.stl", "patches": "mesh-output/p0.json", "output": "mesh-output/p1.json",
 "name": "inlet", "role": "inlet", "point": [0, 0, 0], "normal": [1, 0, 0]}
JSON
cat > walls.json <<'JSON'
{"surface": "test-fixtures/unit-box.stl", "patches": "mesh-output/p1.json", "output": "mesh-output/p2.json",
 "merge": ["south_wall", "north_wall", "floor", "ceiling"], "name": "walls", "role": "wall"}
JSON
echo '{"surface": "test-fixtures/unit-box.stl", "patches": "mesh-output/p2.json"}' > list.json
target/debug/kofl --journal mesh-output/patches.jsonl call patch.init --input init.json
target/debug/kofl --journal mesh-output/patches.jsonl call patch.select_plane --input inlet.json
target/debug/kofl --journal mesh-output/patches.jsonl call patch.merge --input walls.json
target/debug/kofl --journal mesh-output/patches.jsonl call patch.list --input list.json
target/debug/kofl generate-stl test-fixtures/unit-box.stl test-fixtures/stl-sites.json mesh-output/patched.mesh.json --patches mesh-output/p2.json
target/debug/kofl export mesh-output/patched.mesh.json openfoam mesh-output/patched-foam
```

`kofl run --script mesh-output/patches.jsonl` replays the steps. With `--patches` (tools: `patches`), `generate-stl` and `inspect-stl` take the names from the map instead of the solids, and the map's roles are saved with the artifact and used by every export; `facet_patches` and a map cannot be combined. Over MCP the patch tools accept up to 200,000 facets (`max_patch_facets`). The method and its checks are in [research/m2-patch-tools.md](../research/m2-patch-tools.md).

## MCP

Launch the server with read and write roots as in the [M1 workflow](m1-usage.md). Place the STL below a read root before calling the tools.

| Tool | Arguments |
|---|---|
| `inspect_stl_domain` | `surface` filename and inline STL `config`; with a size field it reports the seed count as `site_count` |
| `generate_stl_mesh` | `surface` filename, inline STL `config`, new `artifact` filename, optional `patches` map |
| `patch.*` | `surface`, `patches` and new `output` filenames and the operation's fields ([tools](tools.md)) |

Reuse `check_mesh`, `mesh_report`, and `export_mesh` with the generated artifact. To chain patch calls over MCP, make the write root a read root too, so a call can read the map the previous one wrote. The MCP server limits this prototype to 2,048 STL facets and 20,000 sites per job. The CLI uses the same engine and accepts the STL path as an explicit argument.

## Verification

`cargo test --release --locked -p mesh-engine --test patches` checks the patch tools against the SolidCraft ports recorded by hand, the exact evidence on the boxes, the map checks and the names and roles in every export. `cargo test --release --locked -p mesh-engine --test size_field` builds graded meshes on both fixtures and checks the achieved sizes, determinism, the estimate and the evidence digest. `cargo test --locked --test stl_workflow` compares each unit-box STL cell with the analytic reference, checks boundary patch areas, verifies binary STL with explicit labels, and exercises the skewed fixture plus malformed surfaces. `cargo m1 --require-gmsh` includes these tests in the workspace regression suite.

The Rust `import_closed_stl` entry point also checks shell connectivity, winding, and positive signed volume for nonconvex test surfaces. The [SolidCraft pipe-and-box case](../test-fixtures/solidcraft-vs-foehncraft/README.md) uses it to verify its inlet, outlet, and wall labels before convexity enters the volume meshing path.

The [partitioned nonconvex workflow](nonconvex-duct.md) accepts the elbow and SolidCraft fixtures and exports bulk or layered volume meshes.
