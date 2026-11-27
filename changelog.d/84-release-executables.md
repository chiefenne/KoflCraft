- **Release downloads**: the new workflow Executables builds `kofl` and
  `mesh-mcp` for macOS (one universal build for Apple silicon and Intel),
  Linux x86-64 and Windows x86-64 and attaches them with `SHA256SUMS.txt`
  to the release; the Release workflow starts it after each release, and it
  can be run by hand for an existing one (docs/release.md). The README's
  "Download a release" says which file to take and how to start the
  unsigned Mac programs the first time.
