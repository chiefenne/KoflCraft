# Recorded chats

The journals of the chats on [docs/chats.md](../../docs/chats.md), as
`kofl mcp --journal` wrote them (one line per tool call, `via: mcp`), with
the cloud folder of the recording replaced by the folder of this
directory. Next to them are the files the engineer or the assistant wrote
outside the tools:

| Folder | Chat | Also here |
|---|---|---|
| `pipe-baffle/` | [1. A pipe with a baffle](../../docs/chats.md#pipe-baffle) | `case.toml`: the FoehnCraft skeleton from `bc.skeleton`, with the two baffle patches the assistant added |
| `naca-cgrid/` | [2. A NACA 2412 C-grid for y+ about 1](../../docs/chats.md#naca-yplus) | |
| `manifold-sweep/` | [3. The Manifold with hexahedra in the pipes](../../docs/chats.md#manifold-sweep) | |
| `furnace-7/` | [4. The seven-tube furnace](../../docs/chats.md#furnace-7) | `kofl-mcp.toml`: the policy that raises the read limit |

To replay a journal, put its inputs into an empty folder (the pipe STL
from [make_pipe.py](../gallery/polyhedral/pipe-shells/make_pipe.py), the
Manifold from [test-fixtures/manifold](../../test-fixtures/manifold/), the
furnace from the fixture's Gmsh script as in
[application A2](../../docs/applications.md#furnace-7)), copy the journal
there and run, from that folder:

```
echo '{"script": "journal.jsonl", "keep_going": true}' > replay.json
kofl call journal.replay --input replay.json
```

The replay runs every call again and reports each one whose status, exit
code or result digest differs from its record. A duct build is always
reported as different (`TOOL005`), because its result holds the build's
timings; its mesh files are the same to the bit, which the BLAKE3 hashes
in `mesh/mesh-evidence.json` show (only `report.json`, which holds the
timings, differs). Checked for `pipe-baffle/` on 2026-11-28: calls 2 to 9
the same, call 10 different only in `report.json`. The chats' file edits (the
`no_layers` partition of the furnace, the baffle lines of `case.toml`) are
not tool calls: make them by hand before the call that reads them.
