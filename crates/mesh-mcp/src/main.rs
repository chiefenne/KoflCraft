mod protocol;
mod tools;

use std::io::{self, BufRead, Read, Write};
use std::path::Path;
use std::process::ExitCode;

use mesh_engine::MeshDirectory;

use crate::protocol::{Session, rpc_error};

const MAX_REQUEST_BYTES: usize = 256 * 1024;
const USAGE: &str = "Usage: mesh-mcp --output-dir EXISTING_DIRECTORY";

fn main() -> ExitCode {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        eprintln!("{USAGE}");
        return ExitCode::SUCCESS;
    }
    let [flag, directory] = args.as_slice() else {
        eprintln!("{USAGE}");
        return ExitCode::FAILURE;
    };
    if flag != "--output-dir" {
        eprintln!("{USAGE}");
        return ExitCode::FAILURE;
    }
    match MeshDirectory::open(Path::new(directory))
        .and_then(|dir| serve(dir, io::stdin().lock(), io::stdout().lock()))
    {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("mesh-mcp: {error}");
            ExitCode::FAILURE
        }
    }
}

fn serve(
    directory: MeshDirectory,
    mut input: impl BufRead,
    mut output: impl Write,
) -> io::Result<()> {
    let mut session = Session::new(directory);
    let mut frame = Vec::new();
    loop {
        frame.clear();
        let bytes = (&mut input)
            .take(MAX_REQUEST_BYTES as u64 + 1)
            .read_until(b'\n', &mut frame)?;
        if bytes == 0 {
            return Ok(());
        }
        let response = if bytes > MAX_REQUEST_BYTES {
            if frame.last() != Some(&b'\n') {
                input.skip_until(b'\n')?;
            }
            Some(rpc_error(
                serde_json::Value::Null,
                -32600,
                "request exceeds the size limit",
            ))
        } else {
            session.handle(&frame)
        };
        if let Some(response) = response {
            serde_json::to_writer(&mut output, &response)?;
            output.write_all(b"\n")?;
            output.flush()?;
        }
    }
}
