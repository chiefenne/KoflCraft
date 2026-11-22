//! The MCP server of KoflCraft (`kofl mcp`, `mesh-mcp`). It lists the tools
//! of the registry (`mesh-tools`, AG-12) by their dotted ids and calls them
//! in a confined context: reads only below the read roots, writes only
//! below the write roots (AG-14), with the MCP limits (AG-15). Every
//! result is the registry's envelope (`koflcraft.tool-result`).

mod policy;
mod protocol;

use std::ffi::OsString;
use std::io::{self, BufRead, Read, Write};
use std::path::PathBuf;
use std::process::ExitCode;

use mesh_tools::{Context, Journal};
use serde_json::Value;

pub use crate::policy::{Policy, Rule};
use crate::protocol::{Peer, Session, rpc_error};

const MAX_REQUEST_BYTES: usize = 256 * 1024;
pub const USAGE: &str = "Usage: kofl mcp --read-root DIR [--read-root DIR …] --write-root DIR [--write-root DIR …] [--policy FILE] [--journal FILE]
       mesh-mcp --output-dir DIR   (the same as --read-root DIR --write-root DIR)

Relative paths in tool inputs resolve from the first write root. The
policy file (TOML, outside the write roots) sets `allow`, `ask` or
`refuse` for builds and changes the limits; without it, builds within
the limits run and larger ones ask through the client (AG-15).";

/// The server's options.
#[derive(Debug, Default, PartialEq)]
pub struct Options {
    pub read_roots: Vec<PathBuf>,
    pub write_roots: Vec<PathBuf>,
    pub journal: Option<PathBuf>,
    pub policy: Option<PathBuf>,
}

/// Parse the arguments; the server refuses to start without both a read
/// and a write root (AG-14).
pub fn parse(args: &[OsString]) -> Result<Options, String> {
    let mut o = Options::default();
    let mut it = args.iter();
    while let Some(flag) = it.next() {
        let value = || -> Result<PathBuf, String> {
            Err(format!("{} needs a value", flag.to_string_lossy()))
        };
        let mut next = || it.next().map(PathBuf::from).map_or_else(value, Ok);
        match flag.to_str() {
            Some("--read-root") => o.read_roots.push(next()?),
            Some("--write-root") => o.write_roots.push(next()?),
            Some("--output-dir") => {
                let d = next()?;
                o.read_roots.push(d.clone());
                o.write_roots.push(d);
            }
            Some("--journal") => o.journal = Some(next()?),
            Some("--policy") => o.policy = Some(next()?),
            _ => return Err(format!("unknown argument {}", flag.to_string_lossy())),
        }
    }
    if o.read_roots.is_empty() || o.write_roots.is_empty() {
        return Err("the MCP server needs at least one --read-root and one --write-root".into());
    }
    Ok(o)
}

/// Run the server on standard input and output.
pub fn main(args: Vec<OsString>) -> ExitCode {
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        eprintln!("{USAGE}");
        return ExitCode::SUCCESS;
    }
    let options = match parse(&args) {
        Ok(o) => o,
        Err(e) => {
            eprintln!("mcp: {e}\n{USAGE}");
            return ExitCode::FAILURE;
        }
    };
    let policy = match &options.policy {
        None => Policy::default(),
        Some(path) => match Policy::load(path, &options.write_roots) {
            Ok(p) => p,
            Err(e) => {
                eprintln!("mcp: {e}");
                return ExitCode::FAILURE;
            }
        },
    };
    let ctx = match Context::confined(
        options.read_roots,
        options.write_roots,
        Some(policy.limits.clone()),
        "mcp",
    ) {
        Ok(c) => c,
        Err(d) => {
            eprintln!("mcp: {d}");
            return ExitCode::FAILURE;
        }
    };
    let journal = match Journal::from_env(options.journal.as_deref()) {
        None => None,
        Some(Ok(j)) => Some(j),
        Some(Err(e)) => {
            eprintln!("mcp: cannot open the journal: {e}");
            return ExitCode::FAILURE;
        }
    };
    match serve(
        Session::new(ctx, journal, policy),
        io::stdin().lock(),
        io::stdout().lock(),
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("mcp: {error}");
            ExitCode::FAILURE
        }
    }
}

/// Standard input and output as the transport, with the client as a peer
/// for requests the server sends while a call waits.
struct Transport<R, W> {
    input: R,
    output: W,
    frame: Vec<u8>,
}

impl<R: BufRead, W: Write> Transport<R, W> {
    /// The next frame, or `None` at the end of the input.
    fn read(&mut self) -> io::Result<Option<Result<Vec<u8>, ()>>> {
        self.frame.clear();
        let bytes = (&mut self.input)
            .take(MAX_REQUEST_BYTES as u64 + 1)
            .read_until(b'\n', &mut self.frame)?;
        if bytes == 0 {
            return Ok(None);
        }
        if bytes > MAX_REQUEST_BYTES {
            if self.frame.last() != Some(&b'\n') {
                self.input.skip_until(b'\n')?;
            }
            return Ok(Some(Err(())));
        }
        Ok(Some(Ok(std::mem::take(&mut self.frame))))
    }

    fn write(&mut self, value: &Value) -> io::Result<()> {
        serde_json::to_writer(&mut self.output, value)?;
        self.output.write_all(b"\n")?;
        self.output.flush()
    }
}

impl<R: BufRead, W: Write> Peer for Transport<R, W> {
    fn request(&mut self, request: Value) -> io::Result<Value> {
        self.write(&request)?;
        loop {
            let Some(frame) = self.read()? else {
                return Err(io::Error::other("the client closed the connection"));
            };
            let Ok(frame) = frame else {
                self.write(&rpc_error(
                    Value::Null,
                    -32600,
                    "request exceeds the size limit",
                ))?;
                continue;
            };
            let Ok(message) = serde_json::from_slice::<Value>(&frame) else {
                self.write(&rpc_error(Value::Null, -32700, "invalid JSON"))?;
                continue;
            };
            if message.get("method").is_none() && message["id"] == request["id"] {
                return Ok(message);
            }
            // another request while the approval is open: one call at a time
            if message.get("method").is_some()
                && let Some(id) = message.get("id")
            {
                self.write(&rpc_error(
                    id.clone(),
                    -32000,
                    "busy: waiting for the person to approve a build",
                ))?;
            }
        }
    }
}

fn serve(mut session: Session, input: impl BufRead, output: impl Write) -> io::Result<()> {
    let mut t = Transport {
        input,
        output,
        frame: Vec::new(),
    };
    loop {
        let Some(frame) = t.read()? else {
            return Ok(());
        };
        let response = match frame {
            Err(()) => Some(rpc_error(
                Value::Null,
                -32600,
                "request exceeds the size limit",
            )),
            Ok(frame) => session.handle(&frame, &mut t),
        };
        if let Some(response) = response {
            t.write(&response)?;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Invariant (AG-14): the server does not start without a read and a
    /// write root. Dropping the check would let `kofl mcp` run unconfined.
    #[test]
    fn roots_are_required() {
        let a = |v: &[&str]| v.iter().map(OsString::from).collect::<Vec<_>>();
        assert!(parse(&a(&[])).is_err());
        assert!(parse(&a(&["--read-root", "x"])).is_err());
        assert!(parse(&a(&["--write-root", "x"])).is_err());
        assert!(parse(&a(&["--read-root"])).is_err());
        let o = parse(&a(&[
            "--read-root",
            "r",
            "--write-root",
            "w",
            "--journal",
            "j",
        ]))
        .unwrap();
        assert_eq!(o.read_roots, vec![PathBuf::from("r")]);
        assert_eq!(o.journal, Some(PathBuf::from("j")));
        assert_eq!(
            parse(&a(&["--output-dir", "d"])).unwrap().write_roots.len(),
            1
        );
    }
}
