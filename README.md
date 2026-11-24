<p align="center">
  <img src="assets/logo/koflcraft.png" alt="KoflCraft logo: a stag" width="180" height="180">
</p>

<h1 align="center">KoflCraft</h1>

<p align="center">
  An open-source polyhedral mesh generator for CFD and multi-physics simulations, written in Rust. The companion mesher of <a href="https://github.com/chiefenne/Foehncraft">Foehncraft</a>.
</p>

<p align="center">
  <a href="docs/gallery.md"><img src="assets/readme/tour.webp" width="880" alt="A slideshow of KoflCraft's Manifold mesh, three slides that change every four seconds, pictured with MurmelCraft: the Manifold duct coloured by patch; the Manifold cut in half and coloured by cell type, polyhedra in the core and prisms along the walls; the non-orthogonality of the cells on a flat cut, in the light theme"></a>
</p>

<sub>KoflCraft's Manifold mesh pictured with <a href="https://github.com/chiefenne/MurmelCraft">MurmelCraft</a>:
229,615 cells, Voronoi polyhedra in the core with prism layers at the
walls. The slides are also [still images with captions](docs/gallery.md);
how they were made is in [tools/readme-tour](tools/readme-tour/README.md).</sub>

---

> **Status: milestone M1 done (clipped Voronoi polyhedra in a box, validated,
> exported to VTU and Gmsh, driven from a CLI and an MCP server), M2 under way.**
> STL import, convex STL meshing, prism boundary layers on swept and channel
> geometries, and a first nonconvex duct mesher with layers carried out to the
> ports exist as prototypes. The milestone plan and every decision are in
> [docs/requirements.md](docs/requirements.md).

[![Rust 1.85+](https://img.shields.io/badge/rust-1.85%2B-orange.svg)](Cargo.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## What it is for

KoflCraft makes the volume meshes that a cell-centred finite-volume solver
needs: general polyhedral cells, prism layers at walls, named boundary
patches and, later, several cell zones for conjugate heat transfer. Its first
user is [Foehncraft](https://github.com/chiefenne/Foehncraft), which reads
the meshes directly.

The guiding targets are:

- **Correct before fast.** Every mesh is checked before it is written:
  closed cells, matching interior faces, positive volumes, and conservation
  of the total volume and of every patch area. A failed check is never
  reported as success.
- **Polyhedral and layered.** Clipped Voronoi cells in the core, prism layers
  at walls that end cleanly at inlets and outlets, and conformal joins
  between them.
- **Agent-ready.** One Rust engine behind a command line and an MCP server,
  JSON in and out, every result with its evidence: input digests, algorithm,
  version, timings and validation.
- **Graphics-grade geometry.** Algorithms come first from the computer
  graphics and geometry processing community, above all the work of Keenan
  Crane (CMU) and Bruno Lévy (INRIA), and are benchmarked before they are
  adopted.
- **Clean room.** Algorithms implemented from published papers and documented
  formats, MIT licensed, no copyleft code.

## Build and run

Needs Rust 1.85 or newer. Gmsh on the `PATH` is optional; the acceptance run
uses it for an extra check.

```
cargo test --release --workspace
cargo build --release -p mesh-cli -p mesh-mcp
target/release/mesh-cli generate test-fixtures/box-4.json box.mesh.json
target/release/mesh-cli check box.mesh.json
target/release/mesh-cli export box.mesh.json vtu box.vtu
cargo m1 --require-gmsh
```

## Documentation

| Kind | Where |
|---|---|
| Requirements, decisions and status | [docs/requirements.md](docs/requirements.md) |
| Box meshing, CLI and MCP | [docs/m1-usage.md](docs/m1-usage.md) |
| STL input | [docs/m2-stl.md](docs/m2-stl.md) |
| Block and pipe meshes with wall layers | [docs/analytic-development.md](docs/analytic-development.md) |
| Nonconvex ducts with exposed prism layers | [docs/nonconvex-duct.md](docs/nonconvex-duct.md) |
| Colour maps shared with FoehnCraft and MurmelCraft | [docs/colour-maps.md](docs/colour-maps.md) |
| Rules for agents and contributors | [AGENTS.md](AGENTS.md), [docs/pitfalls.md](docs/pitfalls.md) |

## How it is built

KoflCraft is developed in the way of the
[ArtCraft Craft Anything apps](https://getartcraft.com/apps/craft-anything),
like Foehncraft: Rust, one engine with thin front ends, built step by step
with a coding agent in a study, plan, build, check loop, with every decision
written down in the requirements.

## License

[MIT](LICENSE)
