# M1 box meshing

Build both interfaces with `cargo build --locked -p mesh-cli -p mesh-mcp`; this builds `kofl` too, the command line of the tool registry (`kofl commands` lists every tool, `kofl call ID --input FILE` runs one and prints its result envelope, see [tools.md](tools.md)). `mesh-cli` accepts the same commands as `kofl`. Run the test suite with `cargo test --locked --workspace --all-targets`.

For the primary block and pipe development cases, see [Analytic development meshes](analytic-development.md).

For triangulated domains, see [M2 convex STL meshing](m2-stl.md).

Run `cargo m1 --require-gmsh` for the release acceptance gate. The runner checks 16-, 64-, and 256-site lattice and clustered fixtures through CLI and MCP, reads both export formats independently, and runs `gmsh -check` on each MSH file. It saves configurations, generated meshes, reports, and `evidence.json` under a fresh directory in `target/m1-evidence/`. Each fixture includes one CLI and one MCP measurement of mesh generation and validation time in microseconds. With `--require-foehn` (FoehnCraft as `$FOEHN` or `foehn` on PATH) it also runs the FoehnCraft round trip of FC-06: every box fixture as `foehncraft.vtu`, OpenFOAM `polyMesh` and MSH 4.1, and the straight channel and the bulk, layered and polyhedral elbow ducts as `foehncraft.vtu` and MSH 4.1, are read by `foehn mesh` and must match KoflCraft's summary (see [foehncraft.md](foehncraft.md)); without FoehnCraft the evidence says `not run`. CI builds FoehnCraft from its main branch and runs `cargo m1 --require-gmsh --require-foehn`. `cargo m1` runs the in-repo checks and uses the Gmsh executable when it is on PATH. [Gmsh consistency check](https://gmsh.info/doc/texinfo/gmsh.html#Command_002dline-options)

## CLI

```sh
mkdir -p mesh-output
target/debug/mesh-cli inspect test-fixtures/box-4.json
target/debug/mesh-cli generate test-fixtures/box-4.json mesh-output/box.mesh.json
target/debug/mesh-cli check mesh-output/box.mesh.json
target/debug/mesh-cli report mesh-output/box.mesh.json > mesh-output/report.json
target/debug/mesh-cli export mesh-output/box.mesh.json vtu mesh-output/box.vtu
target/debug/mesh-cli export mesh-output/box.mesh.json gmsh mesh-output/box.msh
```

Open `mesh-output/box.vtu` in ParaView or `mesh-output/box.msh` in Gmsh. The VTU file preserves the polyhedral cells. The MSH 4.1 ASCII file (the format FoehnCraft reads) uses a conforming tetrahedral subdivision with six named surface groups and a `fluid` volume group. A tetrahedron’s elementary entity tag is the original zero-based cell ID plus one. `export ARTIFACT.json gmsh-2.2 OUT.msh` writes the older MSH 2.2 instead. Each generate or export destination must be new, with an existing parent directory. [Gmsh MSH format](https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format)

## Journal and replay

`kofl --journal FILE COMMAND` (or the environment variable `KOFL_JOURNAL`) appends every tool call to a journal (`koflcraft.journal`, JSON lines: a header with the KoflCraft version, then one row per call with the tool id and version, the input, status, exit code, duration and the SHA-256 of the result without its timings). A journal that cannot be written adds the warning TOOL003; the call itself stands.

```sh
target/debug/kofl --journal run.jsonl generate test-fixtures/box-4.json box.mesh.json
target/debug/kofl --journal run.jsonl export box.mesh.json foehncraft box.vtu
mkdir again && cd again
../target/debug/kofl run --script ../run.jsonl [--keep-going]
```

`kofl run --script` (tool `journal.replay`) calls every recorded tool again and prints the `koflcraft.replay` report: for each call the recorded and the replayed status, exit code and result digest. It stops at the first call that ends differently (TOOL005, exit 1) unless `--keep-going` is given. Calls that failed with TOOL001 or TOOL002 and earlier replays are skipped. Configs are recorded inline, but files a call reads (STL, artifacts) and writes are paths: relative paths resolve from the directory the replay runs in, and outputs must not exist yet, so replay in a fresh directory that has the same input files. Export results carry the SHA-256 of every written file, so a replay also compares the file content.

## JSON contract

The analytic box workflow uses `schema_version: 1` and zero-based mesh indices. [The example configuration](../test-fixtures/box-4.json) defines an axis-aligned box with `min` and `max` triples and an ordered list of interior sites. Optional `fixture_id` and `seed` describe how the supplied sites were obtained. The engine uses the supplied coordinates directly.

The saved artifact contains `config`, `build`, and `mesh`. The mesh stores triples in `vertices` and the `topology` object stores `face_offsets`, `face_vertices`, `adjacency`, and `cell_count`. Each adjacency entry is either `{"kind":"interior","owner":0,"neighbour":1}` or `{"kind":"boundary","owner":0,"patch":"xmin"}`. Checks reconstruct cell incidence from these arrays and validate the saved geometry.

Build evidence records the algorithm revision, engine version, OS and architecture, debug-assertion status, and elapsed mesh construction and validation time in microseconds. SHA-256 digests cover the compact JSON serialization of the domain and ordered sites, the full config, the packed topology, and the complete packed mesh respectively. The topology digest excludes point coordinates and build timing. Its equality is checked for repeated builds on the same build and platform.

| CLI command | Tool id (MCP name before the registry) | Successful result |
|---|---|---|
| `inspect` | `domain.inspect` (`inspect_domain`) | Box bounds, site count, fixture and seed, algorithm and limitations |
| `generate` | `mesh.generate` (`generate_mesh`) | Artifact name and build evidence |
| `check` | `mesh.check` (`check_mesh`) | `valid`, mesh counts, domain and mesh volumes, relative volume error, tolerances, topology digest |
| `report` | `mesh.report` (`mesh_report`) | Build evidence, validation result, fixture and seed, and indexed arrays of cell and face metrics |
| `export` | `mesh.export` (`export_mesh`) | Output name, format and the SHA-256 of every written file |

The CLI writes successful JSON to stdout and versioned error objects to stderr with a nonzero exit status. The error object has `code` and `message` fields. Input parsing and file errors can also include a path. Engine codes include `unsupported_schema`, `invalid_domain`, `invalid_sites`, `generation_failed`, `invalid_connectivity`, `digest_mismatch`, `validation_failed`, `unsupported_format`, `export_failed`, `io_error`, and `serialization_failed`.

## MCP

Configure the MCP client to launch `target/debug/kofl mcp --read-root /abs/inputs --write-root /abs/mesh-output` (`mesh-mcp` takes the same arguments). The server refuses to start without at least one read root and one write root (AG-14); both flags may repeat, and `mesh-mcp --output-dir DIR` is short for `--read-root DIR --write-root DIR`. Files are read only below a read root and created only below a write root (symbolic links are resolved first); relative paths resolve from the first write root, and every write creates a new entry. `--journal FILE` records every call with `via: mcp` (see Journal and replay). The stdio server supports protocol versions `2025-11-25` and `2025-06-18`. [MCP stdio transport](https://modelcontextprotocol.io/specification/2025-11-25/basic/transports), [MCP lifecycle](https://modelcontextprotocol.io/specification/2025-11-25/basic/lifecycle), [MCP tools](https://modelcontextprotocol.io/specification/2025-11-25/server/tools)

`tools/list` lists every tool of the registry by its dotted id ([tools.md](tools.md)) with its input schema, an output schema of the result envelope around the tool's output, annotations (`readOnlyHint` for read tools; no tool overwrites, so `destructiveHint` is false) and `_meta` with the action kind, contract version and the input fields that name written files. `tools/call` also accepts the names from before the registry (`generate_mesh`, `export_mesh`, …). Every call returns the `koflcraft.tool-result` envelope as structured content and as text; `isError` is true when its status is `error`, and the registered diagnostic code (AG-18) says why. An unknown tool is a JSON-RPC error `-32602`.

`mesh.report` returns a summary over MCP (AG-16): build evidence, validation, counts and extremes of the cell and face metrics (`summary`), and one page of rows (`page`: `rows` is `cells` or `faces`, then `offset`, `count`, `items` and `next_cursor`). Pass `next_cursor` back as `cursor` to get the next page, until it is null; the cells come first, then the faces. `page_size` sets the rows per page, at most 500. A cursor carries the start of the mesh digest and is refused (TOOL002) for another or a changed artifact. `kofl report` on the command line prints the full report; `kofl call mesh.report` pages only when given `cursor` or `page_size`.

The server bounds each request to 256 KiB, file reads to 8 MiB, and reference workloads (256 sites, 2048 STL facets, 8 channel layers and 32 axial cells, 500 000 sweep cells, duct facets, layers and regions); a larger input fails with TOOL008. Path errors are TOOL006.
