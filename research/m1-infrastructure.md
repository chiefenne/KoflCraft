# M1 infrastructure dependencies

This review applies to the versions selected in [Cargo.lock](../Cargo.lock). Registry metadata and upstream repository locations for all 59 registry packages were checked with `cargo metadata`. License declarations were reviewed alongside the license files available in the downloaded packages. Project-owned clipping, topology, validation, and export algorithms remain in the workspace.

| Direct dependency | Version and role | Upstream |
|---|---|---|
| `serde` | 1.0.229 for typed configuration and evidence | [serde-rs/serde](https://github.com/serde-rs/serde) |
| `serde_json` | 1.0.154 for JSON and binary64 round trips | [serde-rs/json](https://github.com/serde-rs/json) |
| `sha2` | 0.11.1 for SHA-256 artifact and evidence digests | [RustCrypto/hashes](https://github.com/RustCrypto/hashes) |
| `cap-std` | 4.0.3 for directory-scoped filesystem access | [Bytecode Alliance cap-std](https://github.com/bytecodealliance/cap-std) |

The lockfile dependencies offer MIT except the Windows-only `winx` platform helper, which uses Apache-2.0 with the LLVM exception. The `unicode-ident` package also applies Unicode-3.0 to its Unicode-derived tables. These permissive terms are accepted for this dependency set. The project-wide allowlist in LIC-07 remains open. [winx upstream repository](https://github.com/sunfishcode/winx), [unicode-ident licensing](https://docs.rs/crate/unicode-ident/1.0.27)

The local MCP transport uses the published newline-delimited JSON-RPC protocol over stdio. It calls the shared Rust engine and keeps filesystem operations inside a directory capability selected at launch. [MCP transports](https://modelcontextprotocol.io/specification/2025-11-25/basic/transports), [cap-std directory API](https://docs.rs/cap-std/latest/cap_std/fs/struct.Dir.html)
