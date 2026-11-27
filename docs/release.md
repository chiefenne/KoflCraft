# Making a release

A release is the annotated Git tag `v<version>` on `main` plus a GitHub
release whose notes are that version's section of
[CHANGELOG.md](../CHANGELOG.md). Versions follow semantic versioning; before
1.0 a minor step (0.1 to 0.2) may change file formats and tool inputs.

1. In a pull request (branch `release-v<version>-changelog`):
   - set `version` in `[workspace.package]` of `Cargo.toml` and run
     `cargo build` so `Cargo.lock` follows;
   - move every fragment of `changelog.d/` into a new section
     `## <version> (YYYY-MM-DD) – <title>` of `CHANGELOG.md`: a paragraph
     "In one paragraph", the highlights, the known limits and all changes
     grouped by topic; delete the fragments (keep `changelog.d/README.md`);
   - add the version's row at the top of the `## Summary` table;
   - update the status block of the README.
2. Merge it when CI is green.
3. Run the **Release** workflow on `main` (Actions ▸ Release ▸ Run workflow,
   or through the GitHub API) with `version` = the number without the `v`,
   first with *dry run* checked (it prints the notes and the commit it would
   tag), then without. It refuses a version whose changelog section is
   missing or that differs from `Cargo.toml`.

Tags cannot be pushed from the cloud sessions; the workflow creates them.
Locally, after the release: `git checkout main && git pull --tags`.

## Executables

After it creates the release, the Release workflow starts **Executables**
(`.github/workflows/executables.yml`), which builds `kofl` and `mesh-mcp`
and attaches them to the release:
`koflcraft-<version>-macos-universal.tar.gz` (Apple silicon and Intel),
`koflcraft-<version>-linux-x86_64.tar.gz`,
`koflcraft-<version>-windows-x86_64.zip` and `SHA256SUMS.txt`.

To attach them to a release that exists already, run Executables by hand
with `tag` = `v<version>` and *upload* checked (unchecked, it only builds,
and the files are an artifact of the run). A pull request that changes
`executables.yml` builds everything without uploading. The Mac builds are
signed ad hoc only, not with an Apple Developer ID, so macOS asks once
before the first start (README, "Download a release").
