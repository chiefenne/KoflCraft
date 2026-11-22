//! `mesh-mcp`: the MCP server of KoflCraft (same as `kofl mcp`).

fn main() -> std::process::ExitCode {
    mesh_mcp::main(std::env::args_os().skip(1).collect())
}
