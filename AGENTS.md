# Rules for agents and contributors

This file holds the standing rules of KoflCraft (PROC-01 in
[docs/requirements.md](docs/requirements.md)). It is written for coding
agents first and applies to every contributor. **MUST** rules are never
broken; **SHOULD** rules are followed unless the pull request says why not.
Lessons from earlier work are in [docs/pitfalls.md](docs/pitfalls.md).
The rules follow the development method of the ArtCraft Craft Anything apps
in the same form as the companion solver
[FoehnCraft](https://github.com/chiefenne/Foehncraft/blob/main/AGENTS.md).

## Scope and approval

- **MUST** start code changes only after the maintainer's go for the task.
  Documenting decisions already taken needs no extra go.
- **MUST** ask before any step a new commit cannot undo: publishing,
  releasing, deleting a branch or tag.
- **MUST** keep every function reachable without a window: scriptable, with
  JSON inputs and outputs, no interactive-only paths.
- **MUST** add every engine operation a front end calls as a tool of the
  registry (AG-12) once `mesh-tools` exists, and call it from the CLI and
  the MCP server only through the registry. Until then, a new operation
  is added to both `mesh-cli` and `mesh-mcp` with the same behaviour and a
  test that compares their outputs.
- **MUST** give a tool that writes files or starts a build its effects (the
  files it writes, the expected cell count), so the MCP approval table
  (AG-15) can judge the call before it runs.

## Correctness

- **MUST NOT** present a partial or failed result as success, in the CLI,
  MCP output or documentation (AG-08, AG-17). Say what was checked and how.
- **MUST** validate every mesh before export: closed cells, each interior
  face shared by exactly two cells with opposite orientation, each boundary
  face by one, positive volumes, volume and patch-area conservation.
- **MUST NOT** accept heuristic neighbour cutoffs that can drop cell
  fragments (ALG-13).
- **MUST NOT** add `unsafe` (SW-15). **SHOULD NOT** add `unwrap`, `expect`
  or `panic!` outside tests; return a diagnostic instead (SW-16).
- **MUST** check counts a file declares against the file size before
  allocating (SW-17).

## Clean room and licences

- **MUST NOT** read copyleft (GPL, LGPL, AGPL) meshing or geometry source.
- **MAY** read permissively licensed source; **MUST NOT** copy, port or
  translate it (LIC-03, LIC-09).
- **MUST** build algorithms from papers and documented formats and cite them
  in `research/` or the relevant doc.
- **MUST** look for algorithms first in the computer graphics and geometry
  processing literature, above all the work of Keenan Crane (CMU), Bruno
  Lévy (INRIA) and Leif Kobbelt (RWTH Aachen) (OBJ-19), and benchmark the
  candidate against the classical method before choosing, with the outcome in
  a `research/` note (ALG-16). A better method from elsewhere wins if the
  benchmark shows it.
- **MUST** start from the papers. Open permissively licensed source of a
  screened project only when the design note says why the papers are not
  enough and records what was read (LIC-13).
- **MUST NOT** read weak-copyleft source either (MPL, EPL, CDDL, LGPL, or
  LGPL-dual without an explicit permissive option) (LIC-14).
- **MUST** check and record a project's licence in the register of
  requirements section 14a before opening its source (LIC-15).
- **MAY** run Gmsh and other copyleft programs as separate processes for
  tests and fixture preparation only; never link them.
- **MUST** record the provenance of every fixture beside it (LIC-10).
- **MUST NOT** add AGPL, GPL or commercial runtime dependencies; review the
  licence of every new crate (LIC-05).

## Tests

- **MUST** follow the keep rule (SW-18): name the invariant, a change that
  breaks it, and why no other test catches it.
- **MUST** add a failing test with every bug fix, and the input of every
  crash to the reader's broken-input corpus.
- **MUST NOT** delete, weaken, ignore or skip a test to get green.

## Pull requests

- **MUST** name branches by the task (`docs-requirements-craft-alignment`).
- **MUST** keep one issue per pull request; moving code and changing
  behaviour go in separate pull requests.
- **MUST** fill in "What ran and what did not, and why"; a check that did
  not run is never reported as passed (PROC-05).
- **MUST** write the changelog entry as a fragment
  `changelog.d/<PR number>-<short-name>.md` (PROC-06).
- **SHOULD** squash-merge when CI is green and delete the branch.

## This repository

Layout: see section 6 of [docs/requirements.md](docs/requirements.md).
Usage: [docs/m1-usage.md](docs/m1-usage.md),
[docs/m2-stl.md](docs/m2-stl.md),
[docs/foehncraft.md](docs/foehncraft.md),
[docs/analytic-development.md](docs/analytic-development.md),
[docs/nonconvex-duct.md](docs/nonconvex-duct.md).

Checks to run before a pull request:

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --release --locked --workspace --all-targets
cargo m1 --require-gmsh        # needs gmsh on PATH; without it: cargo m1
```

When `cargo xtask layering` exists (SW-19), run it too.
