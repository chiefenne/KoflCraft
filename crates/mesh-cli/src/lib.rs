//! The command line of KoflCraft (`kofl`, and `mesh-cli` with the same
//! commands until it is retired). Every command builds the input of one
//! tool of the registry (`mesh-tools`, AG-12) and calls it; nothing here
//! reaches the engine directly.
//!
//! The classic commands (`generate`, `export`, …) print the tool's result
//! as before, and on failure an error object on standard error. `kofl call`
//! prints the whole result envelope (`koflcraft.tool-result`).

use std::ffi::OsString;
use std::io::{self, Read, Write};
use std::path::Path;
use std::process::ExitCode;

use mesh_tools::{Context, Diagnostic, Envelope, Journal, call_journaled};
use serde_json::{Value, json};

const USAGE: &str = "Usage: kofl [--journal FILE] [--check] COMMAND
  kofl commands [--markdown]
  kofl call TOOL_ID [--input FILE|-]
  kofl run --script JOURNAL [--keep-going]
  kofl mcp --read-root DIR --write-root DIR [--journal FILE]
  kofl diag codes
  kofl info
  kofl inspect CONFIG.json
  kofl generate CONFIG.json ARTIFACT.json
  kofl generate-duct SURFACE.stl PATCHES.json PARTITION.json bulk|layered|polyhedral NEW_OUTPUT_DIR
  kofl generate-sweep CONFIG.json NEW_OUTPUT_DIR
  kofl generate-channel CONFIG.json NEW_OUTPUT_DIR
  kofl inspect-stl SURFACE.stl CONFIG.json
  kofl generate-stl SURFACE.stl CONFIG.json ARTIFACT.json
  kofl check ARTIFACT.json
  kofl report ARTIFACT.json
  kofl export ARTIFACT.json vtu|gmsh|gmsh-2.2|foehncraft OUTPUT
  kofl export ARTIFACT.json openfoam NEW_OUTPUT_DIR

`--check` prints the estimate of a build command (`mesh.estimate`:
cells, memory, output, files) instead of running it.
`--journal FILE` (or `KOFL_JOURNAL`) appends every tool call to a
journal that `kofl run --script` replays. `mesh-cli` accepts the same
commands.";

/// What a command line asks for.
enum Request {
    /// A tool call whose result is printed as the classic commands did.
    Classic(&'static str, Value),
    /// `kofl call`: print the envelope.
    Call(String, Value),
    /// Text to print as it is.
    Text(String),
}

pub fn main(program: &str) -> ExitCode {
    let mut args: Vec<OsString> = std::env::args_os().skip(1).collect();
    let mut journal_flag = None;
    if args.len() >= 2 && args[0] == "--journal" {
        journal_flag = Some(std::path::PathBuf::from(args.remove(1)));
        args.remove(0);
    }
    let mut check = false;
    if args.first().is_some_and(|a| a == "--check") {
        check = true;
        args.remove(0);
    }
    if args.first().is_some_and(|a| a == "mcp") {
        let mut rest = args.split_off(1);
        if let Some(j) = journal_flag {
            rest.extend([OsString::from("--journal"), j.into_os_string()]);
        }
        return mesh_mcp::main(rest);
    }
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        println!("{}", USAGE.replace("kofl ", &format!("{program} ")));
        return ExitCode::SUCCESS;
    }
    let ctx = Context::cli();
    let journal = match Journal::from_env(journal_flag.as_deref()) {
        None => None,
        Some(Ok(j)) => Some(j),
        Some(Err(e)) => return fail(&Diagnostic::new("journal_failed", e)),
    };
    let journal = journal.as_ref();
    let request = parse(&args).map(|r| match (check, r) {
        (true, Request::Classic(id, input)) => {
            Request::Classic("mesh.estimate", json!({ "tool": id, "input": input }))
        }
        (true, Request::Call(id, input)) => Request::Call(
            "mesh.estimate".into(),
            json!({ "tool": id, "input": input }),
        ),
        (_, r) => r,
    });
    match request {
        Err(d) => fail(&d),
        Ok(Request::Text(text)) => {
            print!("{text}");
            ExitCode::SUCCESS
        }
        Ok(Request::Call(id, input)) => {
            let envelope = call_journaled(&id, input, &ctx, journal);
            if let Err(e) = print_json(&envelope) {
                return fail(&Diagnostic::new("io_error", e));
            }
            if envelope.ok() {
                ExitCode::SUCCESS
            } else {
                ExitCode::from(envelope.exit_code.clamp(1, 255) as u8)
            }
        }
        Ok(Request::Classic(id, input)) => classic(call_journaled(id, input, &ctx, journal)),
    }
}

fn classic(envelope: Envelope) -> ExitCode {
    if let Some(d) = envelope.diagnostics.first().filter(|_| !envelope.ok()) {
        return fail(d);
    }
    for d in &envelope.diagnostics {
        eprintln!("warning: {d}");
    }
    match print_json(&envelope.result) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => fail(&Diagnostic::new("io_error", e)),
    }
}

fn print_json(value: &impl serde::Serialize) -> io::Result<()> {
    let mut out = io::stdout().lock();
    serde_json::to_writer_pretty(&mut out, value)?;
    writeln!(out)
}

/// The classic error object: `code` is the snake-case name scripts read,
/// `diagnostic` the registered code (AG-18).
fn fail(d: &Diagnostic) -> ExitCode {
    let mut error = json!({ "code": d.name, "diagnostic": d.code, "message": d.message });
    if let Some(path) = &d.path {
        error["path"] = json!(path);
    }
    eprintln!("{}", json!({ "schema_version": 1, "error": error }));
    ExitCode::FAILURE
}

fn usage() -> Diagnostic {
    Diagnostic::new("usage", USAGE)
}

fn text(arg: &OsString) -> Result<&str, Diagnostic> {
    arg.to_str()
        .ok_or_else(|| Diagnostic::new("invalid_path", "arguments must be UTF-8"))
}

fn read_json(path: &OsString) -> Result<Value, Diagnostic> {
    let p = Path::new(path);
    let bytes =
        std::fs::read(p).map_err(|e| Diagnostic::new("io_error", e).with_path(p.display()))?;
    serde_json::from_slice(&bytes)
        .map_err(|e| Diagnostic::new("invalid_json", e).with_path(p.display()))
}

fn parse(args: &[OsString]) -> Result<Request, Diagnostic> {
    use Request::Classic as C;
    let Some(command) = args.first().map(text).transpose()? else {
        return Err(usage());
    };
    let rest = &args[1..];
    Ok(match (command, rest) {
        ("commands", []) => Request::Text(format!(
            "{}\n",
            serde_json::to_string_pretty(&mesh_tools::catalogue()).unwrap_or_default()
        )),
        ("commands", [flag]) if flag == "--markdown" => Request::Text(mesh_tools::markdown()),
        ("diag", [what]) if what == "codes" => Request::Call("diag.codes".into(), json!({})),
        ("info", []) => Request::Call("system.info".into(), json!({})),
        ("call", [id]) => Request::Call(text(id)?.into(), json!({})),
        ("call", [id, flag, file]) if flag == "--input" => {
            let input = if file == "-" {
                let mut s = String::new();
                io::stdin()
                    .read_to_string(&mut s)
                    .map_err(|e| Diagnostic::new("io_error", e))?;
                serde_json::from_str(&s).map_err(|e| Diagnostic::new("invalid_json", e))?
            } else {
                read_json(file)?
            };
            Request::Call(text(id)?.into(), input)
        }
        ("run", [flag, script]) if flag == "--script" => {
            Request::Call("journal.replay".into(), json!({ "script": text(script)? }))
        }
        ("run", [flag, script, keep]) if flag == "--script" && keep == "--keep-going" => {
            Request::Call(
                "journal.replay".into(),
                json!({ "script": text(script)?, "keep_going": true }),
            )
        }
        ("inspect", [config]) => C("domain.inspect", json!({ "config": read_json(config)? })),
        ("generate", [config, artifact]) => C(
            "mesh.generate",
            json!({ "config": read_json(config)?, "artifact": text(artifact)? }),
        ),
        ("generate-duct", [surface, patches, partition, mode, output]) => C(
            "mesh.generate_duct",
            json!({
                "surface": text(surface)?, "patches": text(patches)?, "partition": text(partition)?,
                "mode": text(mode)?, "output": text(output)?
            }),
        ),
        ("generate-sweep", [config, output]) => C(
            "mesh.generate_sweep",
            json!({ "config": read_json(config)?, "output": text(output)? }),
        ),
        ("generate-channel", [config, output]) => C(
            "mesh.generate_channel",
            json!({ "config": read_json(config)?, "output": text(output)? }),
        ),
        ("inspect-stl", [surface, config]) => C(
            "domain.inspect_stl",
            json!({ "surface": text(surface)?, "config": read_json(config)? }),
        ),
        ("generate-stl", [surface, config, artifact]) => C(
            "mesh.generate_stl",
            json!({ "surface": text(surface)?, "config": read_json(config)?, "artifact": text(artifact)? }),
        ),
        ("check", [artifact]) => C("mesh.check", json!({ "artifact": text(artifact)? })),
        ("report", [artifact]) => C("mesh.report", json!({ "artifact": text(artifact)? })),
        ("export", [artifact, format, output]) => C(
            "mesh.export",
            json!({ "artifact": text(artifact)?, "format": text(format)?, "output": text(output)? }),
        ),
        _ => return Err(usage()),
    })
}
