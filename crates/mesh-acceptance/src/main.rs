#[path = "../../../test-fixtures/sites.rs"]
mod fixtures;
mod foehn;
mod readers;

use std::error::Error;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode, Output, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

use mesh_core::geometry::Vec3;
use mesh_engine::{BoxBounds, MeshConfig, SCHEMA_VERSION};
use serde::Serialize;
use serde_json::{Value, json};

type RunResult<T> = Result<T, Box<dyn Error>>;

#[derive(Serialize)]
struct FixtureEvidence {
    fixture_id: String,
    seed: Option<u64>,
    cells: usize,
    topology_digest: Value,
    mesh_digest: Value,
    relative_volume_error: f64,
    vtu_polyhedra: usize,
    gmsh_tetrahedra: usize,
    cli_build: Value,
    mcp_build: Value,
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("M1 acceptance failed: {error}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> RunResult<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    const USAGE: &str = "usage: cargo m1 [--require-gmsh] [--require-foehn]";
    if args == ["--help"] {
        println!(
            "Usage: cargo m1 [--require-gmsh] [--require-foehn]\n\nFoehnCraft is `$FOEHN` or `foehn` on PATH."
        );
        return Ok(());
    }
    require(
        args.iter()
            .all(|a| a == "--require-gmsh" || a == "--require-foehn"),
        USAGE,
    )?;
    let gmsh = gmsh_version()?;
    require(
        !args.iter().any(|a| a == "--require-gmsh") || gmsh.is_some(),
        "gmsh must be on PATH",
    )?;
    let foehncraft = foehn::find()?;
    require(
        !args.iter().any(|a| a == "--require-foehn") || foehncraft.is_some(),
        "foehn must be on PATH or in $FOEHN",
    )?;

    let workspace = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()?;
    let target = workspace.join("target");
    let cargo = option_env!("CARGO").unwrap_or("cargo");
    let status = |args: &[&str]| {
        Command::new(cargo)
            .args(args)
            .arg("--target-dir")
            .arg(&target)
            .current_dir(&workspace)
            .status()
    };
    require(
        status(&[
            "test",
            "--locked",
            "--release",
            "--workspace",
            "--all-targets",
        ])?
        .success(),
        "workspace tests failed",
    )?;
    require(
        status(&[
            "build",
            "--locked",
            "--release",
            "-p",
            "mesh-cli",
            "-p",
            "mesh-mcp",
        ])?
        .success(),
        "could not build the release CLI and MCP binaries",
    )?;
    let cli = binary(&target, "mesh-cli");
    let mcp = binary(&target, "mesh-mcp");
    let time = SystemTime::now().duration_since(UNIX_EPOCH)?.as_millis();
    let root = target
        .join("m1-evidence")
        .join(format!("{time}-{}", std::process::id()));
    fs::create_dir_all(root.parent().unwrap())?;
    fs::create_dir(&root)?;
    println!("Writing acceptance artifacts to {}", root.display());

    let mut runs = Vec::new();
    for dimensions in fixtures::LATTICES {
        let lattice = fixtures::lattice(dimensions);
        let count = lattice.len();
        for (distribution, sites, seed) in [
            ("lattice", lattice, None),
            (
                "clustered",
                fixtures::clustered(count),
                Some(fixtures::SEED),
            ),
        ] {
            let fixture = format!("{distribution}-{count}");
            let directory = root.join(&fixture);
            fs::create_dir(&directory)?;
            let config = MeshConfig {
                schema_version: SCHEMA_VERSION,
                domain: BoxBounds {
                    min: [0.0; 3],
                    max: [1.0; 3],
                },
                sites: sites
                    .into_iter()
                    .map(|Vec3 { x, y, z }| [x, y, z])
                    .collect(),
                fixture_id: Some(fixture.clone()),
                seed,
            };
            let evidence = check_fixture(&cli, &mcp, &directory, config, gmsh.is_some())
                .map_err(|error| format!("{fixture}: {error}"))?;
            println!(
                "{fixture}: {} VTU polyhedra, {} Gmsh tetrahedra",
                evidence.vtu_polyhedra, evidence.gmsh_tetrahedra
            );
            runs.push(evidence);
        }
    }
    let foehncraft_evidence = match &foehncraft {
        Some((path, version)) => {
            let checks = round_trip(&cli, path, &root, &workspace)?;
            json!({ "status": "passed", "version": version, "tolerance": foehn::TOLERANCE, "checks": checks })
        }
        None => {
            println!("FoehnCraft round trip (FC-06) not run: no foehn on PATH or in $FOEHN");
            json!({ "status": "not run", "reason": "no foehn on PATH or in $FOEHN" })
        }
    };
    let summary = json!({ "schema_version": SCHEMA_VERSION, "profile": "release", "gmsh_version": gmsh, "fixtures": runs, "foehncraft": foehncraft_evidence });
    let evidence = root.join("evidence.json");
    write_json(&evidence, &summary)?;
    println!("M1 acceptance passed. Evidence at {}", evidence.display());
    Ok(())
}

fn check_fixture(
    cli: &Path,
    mcp: &Path,
    directory: &Path,
    config: MeshConfig,
    check_gmsh: bool,
) -> RunResult<FixtureEvidence> {
    write_json(&directory.join("config.json"), &config)?;
    let cli_inspect = run_cli(cli, directory, &["inspect", "config.json"])?;
    let cli_build = run_cli(
        cli,
        directory,
        &["generate", "config.json", "cli.mesh.json"],
    )?;
    run_cli(
        cli,
        directory,
        &["export", "cli.mesh.json", "vtu", "cli.vtu"],
    )?;
    run_cli(
        cli,
        directory,
        &["export", "cli.mesh.json", "gmsh", "cli.msh"],
    )?;
    let cli_report = run_cli(cli, directory, &["report", "cli.mesh.json"])?;
    write_json(&directory.join("report.json"), &cli_report)?;

    let mcp_results = run_mcp(mcp, directory, &config)?;
    let [mcp_inspect, mcp_build, mcp_check, _, _, mcp_report] = mcp_results.as_slice() else {
        return Err("unexpected MCP result count".into());
    };
    require(
        *mcp_inspect == cli_inspect,
        "CLI and MCP input inspections differ",
    )?;
    let cli_check = run_cli(cli, directory, &["check", "mcp.mesh.json"])?;
    require(
        cli_check == *mcp_check && cli_check == cli_report["validation"],
        "cross-interface mesh checks differ",
    )?;
    require(
        cli_check["valid"] == true && cli_check["cells"] == config.sites.len(),
        "invalid cell count",
    )?;

    let cli_build = &cli_build["build"];
    let mcp_build = &mcp_build["build"];
    for name in [
        "input_digest",
        "config_digest",
        "topology_digest",
        "mesh_digest",
    ] {
        require(
            cli_build[name]
                .as_str()
                .is_some_and(|digest| digest.len() == 64)
                && cli_build[name] == mcp_build[name],
            &format!("CLI and MCP {name} differ"),
        )?;
    }
    for name in ["validation", "cells", "faces", "fixture_id", "seed"] {
        require(
            cli_report[name] == mcp_report[name],
            &format!("CLI and MCP report {name} differ"),
        )?;
    }
    require(
        cli_build["debug_assertions"] == false && mcp_build["debug_assertions"] == false,
        "expected release build evidence",
    )?;
    require(
        cli_report["build"] == *cli_build && mcp_report["build"] == *mcp_build,
        "saved build evidence differs",
    )?;

    let relative_volume_error = cli_check["relative_volume_error"]
        .as_f64()
        .ok_or("missing volume error")?;
    require(
        relative_volume_error <= 1e-10,
        "box volume failed conservation",
    )?;
    let expected: Vec<_> = cli_report["cells"]
        .as_array()
        .ok_or("missing cell metrics")?
        .iter()
        .map(|cell| cell["volume"].as_f64().ok_or("missing cell volume"))
        .collect::<Result<_, _>>()?;
    let vtu = fs::read_to_string(directory.join("cli.vtu"))?;
    let msh = fs::read_to_string(directory.join("cli.msh"))?;
    require(
        vtu == fs::read_to_string(directory.join("mcp.vtu"))?,
        "CLI and MCP VTU files differ",
    )?;
    require(
        msh == fs::read_to_string(directory.join("mcp.msh"))?,
        "CLI and MCP Gmsh files differ",
    )?;
    let vtu_polyhedra = readers::vtu(&vtu, &expected)?;
    let gmsh_tetrahedra = readers::gmsh(&msh, &expected)?;
    if check_gmsh {
        let log = capture(
            Command::new("gmsh")
                .args(["cli.msh", "-check", "-nopopup"])
                .current_dir(directory)
                .output()?,
        )?;
        require(
            !log.lines()
                .any(|line| line.trim_start().starts_with("Error")),
            &log,
        )?;
        fs::write(directory.join("gmsh-check.log"), log)?;
    }

    Ok(FixtureEvidence {
        fixture_id: config.fixture_id.unwrap(),
        seed: config.seed,
        cells: config.sites.len(),
        topology_digest: cli_build["topology_digest"].clone(),
        mesh_digest: cli_build["mesh_digest"].clone(),
        relative_volume_error,
        vtu_polyhedra,
        gmsh_tetrahedra,
        cli_build: cli_build.clone(),
        mcp_build: mcp_build.clone(),
    })
}

/// FC-06: the box fixtures as `foehncraft.vtu`, `polyMesh` and MSH 4.1, and
/// the channel and duct bundles as `foehncraft.vtu` and MSH 4.1, each read
/// by FoehnCraft and compared with KoflCraft's summary.
fn round_trip(cli: &Path, foehn: &Path, root: &Path, workspace: &Path) -> RunResult<Vec<Value>> {
    let mut checks = Vec::new();
    let mut run = |dir: &Path, file: &str, evidence: &Path, kind: foehn::Kind| -> RunResult<()> {
        let expected: Value = serde_json::from_slice(&fs::read(evidence)?)?;
        let out = dir.join(format!("foehn-{}", file.replace('/', "-")));
        let check = foehn::check(foehn, &dir.join(file), &expected["summary"], kind, &out)?;
        let limitation = (check["foehncraft_valid"] == false).then(|| {
            format!(
                "FoehnCraft calls {} invalid: {} faces with S·d <= 0, as KoflCraft predicted",
                dir.join(file).display(),
                check["inverted_faces"]
            )
        });
        let mut check = check;
        if let Some(l) = limitation {
            println!("limitation: {l}");
            check["limitation"] = json!(l);
        }
        checks.push(check);
        Ok(())
    };
    for entry in fs::read_dir(root)? {
        let dir = entry?.path();
        if !dir.join("cli.mesh.json").is_file() {
            continue;
        }
        run_cli(
            cli,
            &dir,
            &["export", "cli.mesh.json", "foehncraft", "cli.foehn.vtu"],
        )?;
        run_cli(
            cli,
            &dir,
            &["export", "cli.mesh.json", "openfoam", "cli.foam"],
        )?;
        let evidence = dir.join("cli.foam/mesh-evidence.json");
        run(&dir, "cli.foehn.vtu", &evidence, foehn::Kind::Faces)?;
        run(&dir, "cli.foam", &evidence, foehn::Kind::Faces)?;
        run(&dir, "cli.msh", &evidence, foehn::Kind::Subdivision)?;
    }
    let bundles = root.join("foehncraft");
    fs::create_dir(&bundles)?;
    let fixtures = workspace.join("test-fixtures");
    let path = |p: &str| fixtures.join(p).display().to_string();
    run_cli(
        cli,
        &bundles,
        &[
            "generate-channel",
            &path("straight-channel.json"),
            "channel",
        ],
    )?;
    run(
        &bundles,
        "channel/foehncraft.vtu",
        &bundles.join("channel/mesh-evidence.json"),
        foehn::Kind::Faces,
    )?;
    run(
        &bundles,
        "channel/channel.msh",
        &bundles.join("channel/mesh-evidence.json"),
        foehn::Kind::Subdivision,
    )?;
    for mode in ["bulk", "layered", "polyhedral"] {
        let name = format!("elbow-{mode}");
        run_cli(
            cli,
            &bundles,
            &[
                "generate-duct",
                &path("elbow/elbow.stl"),
                &path("elbow/patches.json"),
                &path("elbow/partition.json"),
                mode,
                &name,
            ],
        )?;
        let evidence = bundles.join(&name).join("mesh-evidence.json");
        run(
            &bundles,
            &format!("{name}/foehncraft.vtu"),
            &evidence,
            foehn::Kind::Faces,
        )?;
        run(
            &bundles,
            &format!("{name}/mesh.msh"),
            &evidence,
            foehn::Kind::Subdivision,
        )?;
    }
    println!(
        "FoehnCraft round trip (FC-06): {} files read and matched",
        checks.len()
    );
    Ok(checks)
}

fn run_cli(cli: &Path, directory: &Path, args: &[&str]) -> RunResult<Value> {
    let output = Command::new(cli)
        .args(args)
        .current_dir(directory)
        .output()?;
    if !output.status.success() {
        return Err(format!(
            "mesh-cli {}: {}",
            args.join(" "),
            String::from_utf8_lossy(&output.stderr)
        )
        .into());
    }
    Ok(serde_json::from_slice(&output.stdout)?)
}

fn run_mcp(mcp: &Path, directory: &Path, config: &MeshConfig) -> RunResult<Vec<Value>> {
    let call = |id, name, arguments| {
        json!({
            "jsonrpc": "2.0", "id": id, "method": "tools/call", "params": { "name": name, "arguments": arguments }
        })
    };
    let requests = [
        json!({
            "jsonrpc": "2.0", "id": 0, "method": "initialize",
            "params": {
                "protocolVersion": "2025-11-25", "capabilities": {},
                "clientInfo": { "name": "mesh-acceptance", "version": "1" }
            }
        }),
        json!({ "jsonrpc": "2.0", "method": "notifications/initialized" }),
        call(1, "inspect_domain", json!({ "config": config })),
        call(
            2,
            "generate_mesh",
            json!({ "config": config, "artifact": "mcp.mesh.json" }),
        ),
        call(3, "check_mesh", json!({ "artifact": "cli.mesh.json" })),
        call(
            4,
            "export_mesh",
            json!({ "artifact": "mcp.mesh.json", "format": "vtu", "output": "mcp.vtu" }),
        ),
        call(
            5,
            "export_mesh",
            json!({ "artifact": "mcp.mesh.json", "format": "gmsh", "output": "mcp.msh" }),
        ),
        call(6, "mesh_report", json!({ "artifact": "mcp.mesh.json" })),
    ];
    let mut server = Command::new(mcp)
        .arg("--output-dir")
        .arg(directory)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;
    let mut input = server.stdin.take().unwrap();
    let writer = std::thread::spawn(move || -> io::Result<()> {
        for request in requests {
            writeln!(input, "{request}")?;
        }
        Ok(())
    });
    let output = server.wait_with_output()?;
    writer.join().map_err(|_| "MCP request writer panicked")??;
    require(
        output.status.success(),
        &String::from_utf8_lossy(&output.stderr),
    )?;
    let stdout = String::from_utf8(output.stdout)?;
    let mut responses = stdout.lines().map(serde_json::from_str::<Value>);
    let initialized = responses.next().ok_or("missing MCP initialization")??;
    require(
        initialized["id"] == 0 && initialized["result"]["protocolVersion"] == "2025-11-25",
        "MCP initialization failed",
    )?;
    responses
        .enumerate()
        .map(|(index, response)| {
            let response = response?;
            require(
                response["id"] == index + 1 && response["result"]["isError"] == false,
                &format!("MCP tool call failed: {response}"),
            )?;
            Ok(response["result"]["structuredContent"].clone())
        })
        .collect()
}

fn binary(target: &Path, name: &str) -> PathBuf {
    target
        .join("release")
        .join(format!("{name}{}", std::env::consts::EXE_SUFFIX))
}

fn gmsh_version() -> RunResult<Option<String>> {
    match Command::new("gmsh").arg("-version").output() {
        Ok(output) => Ok(Some(capture(output)?.trim().to_owned())),
        Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(None),
        Err(error) => Err(error.into()),
    }
}

fn capture(output: Output) -> RunResult<String> {
    let log = format!(
        "{}{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    require(output.status.success(), &log)?;
    Ok(log)
}

fn write_json(path: &Path, value: &impl Serialize) -> RunResult<()> {
    let mut data = serde_json::to_vec_pretty(value)?;
    data.push(b'\n');
    fs::write(path, data)?;
    Ok(())
}

fn require(condition: bool, message: &str) -> RunResult<()> {
    if condition {
        Ok(())
    } else {
        Err(message.to_owned().into())
    }
}
