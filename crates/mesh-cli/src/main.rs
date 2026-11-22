use std::ffi::{OsStr, OsString};
use std::fmt;
use std::fs::{self, File};
use std::io::{self, BufReader, Write};
use std::path::Path;
use std::process::ExitCode;

use mesh_engine::channel::{ChannelConfig, generate_channel};
use mesh_engine::{
    EngineError, MeshArtifact, MeshConfig, MeshDirectory, SCHEMA_VERSION, StlMeshConfig,
    build_mesh, build_stl_mesh, check_mesh, export_mesh, inspect_domain, inspect_stl_domain,
    mesh_report,
};
use serde_json::{Value, json};

const USAGE: &str = "Usage:
  mesh-cli inspect CONFIG.json
  mesh-cli generate CONFIG.json ARTIFACT.json
  mesh-cli generate-duct SURFACE.stl PATCHES.json PARTITION.json bulk|layered|polyhedral NEW_OUTPUT_DIR
  mesh-cli generate-sweep CONFIG.json NEW_OUTPUT_DIR
  mesh-cli generate-channel CONFIG.json NEW_OUTPUT_DIR
  mesh-cli inspect-stl SURFACE.stl CONFIG.json
  mesh-cli generate-stl SURFACE.stl CONFIG.json ARTIFACT.json
  mesh-cli check ARTIFACT.json
  mesh-cli report ARTIFACT.json
  mesh-cli export ARTIFACT.json vtu OUTPUT.vtu
  mesh-cli export ARTIFACT.json gmsh|gmsh-2.2 OUTPUT.msh
  mesh-cli export ARTIFACT.json foehncraft OUTPUT.vtu
  mesh-cli export ARTIFACT.json openfoam NEW_OUTPUT_DIR";

fn main() -> ExitCode {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        println!("{USAGE}");
        return ExitCode::SUCCESS;
    }
    match dispatch(&args).and_then(|result| {
        let mut output = io::stdout().lock();
        serde_json::to_writer_pretty(&mut output, &result)
            .map_err(|error| Failure::message("io_error", error))?;
        writeln!(output).map_err(|error| Failure::message("io_error", error))
    }) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!(
                "{}",
                json!({ "schema_version": SCHEMA_VERSION, "error": error.0 })
            );
            ExitCode::FAILURE
        }
    }
}

fn dispatch(args: &[OsString]) -> Result<Value, Failure> {
    match args {
        [command, path] if command == "inspect" => {
            let config: MeshConfig = read_json(path.as_ref())?;
            Ok(json!(inspect_domain(&config)?))
        }
        [command, config, artifact] if command == "generate" => {
            let config: MeshConfig = read_json(config.as_ref())?;
            let artifact_path = Path::new(artifact);
            let artifact = build_mesh(config)?;
            let (dir, name) = target(artifact_path)?;
            artifact.write_new(&dir, name)?;
            Ok(json!({
                "schema_version": SCHEMA_VERSION,
                "artifact": artifact_path.display().to_string(),
                "build": artifact.build()
            }))
        }
        [command, surface, patches, plan, mode, output] if command == "generate-duct" => {
            let mode = match mode.to_str() {
                Some("bulk") => mesh_engine::duct::Mode::Bulk,
                Some("layered") => mesh_engine::duct::Mode::Layered,
                Some("polyhedral") => mesh_engine::duct::Mode::Polyhedral,
                _ => return Err(Failure::message("usage", USAGE)),
            };
            let (dir, name) = target(Path::new(output))?;
            let partition: mesh_engine::duct::Partition = read_json(plan.as_ref())?;
            let scaffold = partition
                .scaffold_file(mode)
                .map(|file| read_surface(&Path::new(plan).with_file_name(file)))
                .transpose()?;
            Ok(mesh_engine::duct::generate(
                &read_surface(Path::new(surface))?,
                &read_json(patches.as_ref())?,
                &partition,
                mode,
                scaffold.as_deref(),
                &dir,
                name,
            )?)
        }
        [command, config, output] if command == "generate-sweep" => {
            let config = read_json(config.as_ref())?;
            let path = Path::new(output);
            let (dir, name) = target(path)?;
            let report = mesh_engine::sweep::generate(&config, &dir, name, usize::MAX / 1024)?;
            Ok(
                json!({"schema_version":SCHEMA_VERSION,"output":path.display().to_string(),"report":report}),
            )
        }
        [command, config, output] if command == "generate-channel" => {
            let config: ChannelConfig = read_json(config.as_ref())?;
            let path = Path::new(output);
            let (dir, name) = target(path)?;
            let report = generate_channel(&config, &dir, name)?;
            Ok(
                json!({ "schema_version": SCHEMA_VERSION, "output": path.display().to_string(), "report": report }),
            )
        }
        [command, surface, config] if command == "inspect-stl" => {
            let config: StlMeshConfig = read_json(config.as_ref())?;
            Ok(json!(inspect_stl_domain(
                &read_surface(Path::new(surface))?,
                &config,
                None
            )?))
        }
        [command, surface, config, artifact] if command == "generate-stl" => {
            let config: StlMeshConfig = read_json(config.as_ref())?;
            let path = Path::new(artifact);
            let artifact = build_stl_mesh(&read_surface(Path::new(surface))?, config, None)?;
            let (dir, name) = target(path)?;
            artifact.write_new(&dir, name)?;
            Ok(json!({
                "schema_version": SCHEMA_VERSION,
                "artifact_schema_version": artifact.schema_version(),
                "artifact": path.display().to_string(),
                "build": artifact.build()
            }))
        }
        [command, path] if command == "check" => {
            let artifact: MeshArtifact = read_json(path.as_ref())?;
            Ok(json!(check_mesh(&artifact)?))
        }
        [command, path] if command == "report" => {
            let artifact: MeshArtifact = read_json(path.as_ref())?;
            Ok(json!(mesh_report(&artifact)?))
        }
        [command, artifact, format, output] if command == "export" => {
            let format = format.to_str().unwrap_or_default().parse()?;
            let artifact: MeshArtifact = read_json(artifact.as_ref())?;
            let path = Path::new(output);
            let data = export_mesh(&artifact, format)?;
            let (dir, name) = target(path)?;
            data.write_new(&dir, name)?;
            Ok(json!({ "schema_version": SCHEMA_VERSION, "output": path.display().to_string() }))
        }
        _ => Err(Failure::message("usage", USAGE)),
    }
}

fn read_surface(path: &Path) -> Result<Vec<u8>, Failure> {
    fs::read(path).map_err(|error| Failure::path("io_error", path, error))
}

fn read_json<T: serde::de::DeserializeOwned>(path: &Path) -> Result<T, Failure> {
    let file = File::open(path).map_err(|error| Failure::path("io_error", path, error))?;
    serde_json::from_reader(BufReader::new(file))
        .map_err(|error| Failure::path("invalid_json", path, error))
}

fn target(path: &Path) -> Result<(MeshDirectory, &OsStr), Failure> {
    let name = path
        .file_name()
        .ok_or_else(|| Failure::path("io_error", path, "expected an output name"))?;
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let dir =
        MeshDirectory::open(parent).map_err(|error| Failure::path("io_error", parent, error))?;
    Ok((dir, name))
}

struct Failure(Value);

impl Failure {
    fn message(code: &str, message: impl fmt::Display) -> Self {
        Self(json!({ "code": code, "message": message.to_string() }))
    }

    fn path(code: &str, path: &Path, message: impl fmt::Display) -> Self {
        Self(
            json!({ "code": code, "path": path.display().to_string(), "message": message.to_string() }),
        )
    }
}

impl From<EngineError> for Failure {
    fn from(error: EngineError) -> Self {
        Self(json!(error))
    }
}
