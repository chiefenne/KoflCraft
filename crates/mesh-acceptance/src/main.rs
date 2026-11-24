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

/// FC-06: the box fixtures as `foehncraft.vtu`, `foehncraft.vtkhdf`,
/// `polyMesh` and MSH 4.1, and the channel and duct bundles as their
/// default `foehncraft.vtkhdf` and MSH 4.1 (the channel also with
/// `--handoff vtu`), each read by FoehnCraft and compared with KoflCraft's
/// summary.
fn round_trip(cli: &Path, foehn: &Path, root: &Path, workspace: &Path) -> RunResult<Vec<Value>> {
    let mut checks = Vec::new();
    let mut run = |dir: &Path, file: &str, evidence: &Path, kind: foehn::Kind| -> RunResult<()> {
        let expected: Value = serde_json::from_slice(&fs::read(evidence)?)?;
        let out = dir.join(format!("foehn-{}", file.replace('/', "-")));
        let check = foehn::check(foehn, &dir.join(file), &expected["summary"], kind, &out)?;
        // Every export must be valid for FoehnCraft: no face with S·d <= 0
        // (the elbow and Manifold polyhedral cores once had them).
        require(
            check["foehncraft_valid"] == true,
            &format!(
                "FoehnCraft calls {} invalid: {} faces with S·d <= 0",
                dir.join(file).display(),
                check["inverted_faces"]
            ),
        )?;
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
            &[
                "export",
                "cli.mesh.json",
                "foehncraft-vtkhdf",
                "cli.foehn.vtkhdf",
            ],
        )?;
        run_cli(
            cli,
            &dir,
            &["export", "cli.mesh.json", "openfoam", "cli.foam"],
        )?;
        let evidence = dir.join("cli.foam/mesh-evidence.json");
        run(&dir, "cli.foehn.vtu", &evidence, foehn::Kind::Faces)?;
        run(&dir, "cli.foehn.vtkhdf", &evidence, foehn::Kind::Faces)?;
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
        "channel/foehncraft.vtkhdf",
        &bundles.join("channel/mesh-evidence.json"),
        foehn::Kind::Faces,
    )?;
    run_cli(
        cli,
        &bundles,
        &[
            "generate-channel",
            &path("straight-channel.json"),
            "channel-vtu",
            "--handoff",
            "vtu",
        ],
    )?;
    run(
        &bundles,
        "channel-vtu/foehncraft.vtu",
        &bundles.join("channel-vtu/mesh-evidence.json"),
        foehn::Kind::Faces,
    )?;
    run(
        &bundles,
        "channel/channel.msh",
        &bundles.join("channel/mesh-evidence.json"),
        foehn::Kind::Subdivision,
    )?;
    for (case, stl, mode) in ["elbow", "manifold"]
        .into_iter()
        .zip(["elbow.stl", "manifold.stl"])
        .flat_map(|(c, s)| ["bulk", "layered", "polyhedral"].map(|m| (c, s, m)))
    {
        let name = format!("{case}-{mode}");
        run_cli(
            cli,
            &bundles,
            &[
                "generate-duct",
                &path(&format!("{case}/{stl}")),
                &path(&format!("{case}/patches.json")),
                &path(&format!("{case}/partition.json")),
                mode,
                &name,
            ],
        )?;
        let evidence = bundles.join(&name).join("mesh-evidence.json");
        run(
            &bundles,
            &format!("{name}/foehncraft.vtkhdf"),
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
    // MB1 (BS-13, BS-14): block-structured meshes from the templates, with
    // TFI and with elliptic smoothing, and one split with the edit tools;
    // FoehnCraft reads each back, and smoothing must not make the worst
    // non-orthogonality FoehnCraft measures larger than TFI's
    let blocks = bundles.join("block");
    fs::create_dir(&blocks)?;
    let block_cases = [
        (
            "pipe",
            json!({"template": "pipe_ogrid", "radius": 0.05, "length": 0.5,
                   "core_cells": 8, "radial_cells": 10, "axial_cells": 40, "wall_spacing": 5e-4}),
        ),
        (
            "elbow",
            json!({"template": "elbow", "radius": 0.05, "bend_radius": 0.15, "angle": 90,
                   "inlet_length": 0.1, "outlet_length": 0.1,
                   "core_cells": 8, "radial_cells": 10, "axial_cells": 40, "wall_spacing": 5e-4}),
        ),
        (
            "airfoil",
            json!({"template": "airfoil_cgrid", "naca": "2412", "span": 0.1, "farfield": 5,
                   "wake_length": 5, "surface_cells": 40, "normal_cells": 20, "wake_cells": 20,
                   "wall_spacing": 2e-3}),
        ),
    ];
    let non_orthogonality = |dir: &Path, file: &str| -> RunResult<f64> {
        let out = dir.join(format!("foehn-{}", file.replace('/', "-")));
        let report: Value = serde_json::from_slice(&fs::read(out.join("mesh_report.json"))?)?;
        report["quality"]["non_orthogonality"]["max"]
            .as_f64()
            .ok_or_else(|| "missing non_orthogonality in the FoehnCraft report".into())
    };
    let mut block_quality = Vec::new();
    for (case, config) in &block_cases {
        let topology = format!("{case}.json");
        let call = |tool: &str, file: &str, input: Value| -> RunResult<()> {
            write_json(&blocks.join(file), &input)?;
            let envelope = run_cli(cli, &blocks, &["call", tool, "--input", file])?;
            require(
                envelope["status"] == "finished",
                &format!("{file}: {envelope}"),
            )
        };
        call(
            "block.template",
            &format!("template-{case}.json"),
            json!({"config": config, "output": topology}),
        )?;
        let mut worst = Vec::new();
        // and once more with the wall nodes sliding on the wall
        let kinds = [
            ("tfi", json!({"kind": "tfi"})),
            ("elliptic", json!({"kind": "elliptic"})),
            ("sliding", json!({"kind": "elliptic", "wall_sliding": true})),
        ];
        for (kind, algorithm) in &kinds {
            let name = format!("{case}-{kind}");
            call(
                "block.grid",
                &format!("grid-{name}.json"),
                json!({"topology": topology, "output": name, "algorithm": algorithm}),
            )?;
            let evidence = blocks.join(&name).join("mesh-evidence.json");
            let vtu = format!("{name}/foehncraft.vtu");
            run(&blocks, &vtu, &evidence, foehn::Kind::Faces)?;
            run(
                &blocks,
                &format!("{name}/mesh.msh"),
                &evidence,
                foehn::Kind::Subdivision,
            )?;
            worst.push(non_orthogonality(&blocks, &vtu)?);
        }
        for (k, (kind, _)) in kinds.iter().enumerate().skip(1) {
            require(
                worst[k] <= worst[0] + 1e-9,
                &format!(
                    "{case}: {kind} smoothing raises the max non-orthogonality from {:.2}° (TFI) to {:.2}°",
                    worst[0], worst[k]
                ),
            )?;
        }
        println!(
            "block {case}: FoehnCraft max non-orthogonality TFI {:.2}°, elliptic {:.2}°{}",
            worst[0],
            worst[1],
            worst
                .get(2)
                .map(|w| format!(", sliding walls {w:.2}°"))
                .unwrap_or_default()
        );
        let mut q = json!({"case": case, "tfi_deg": worst[0], "elliptic_deg": worst[1]});
        if let Some(w) = worst.get(2) {
            q["sliding_deg"] = json!(w);
        }
        block_quality.push(q);
    }
    // BS-04: a pipe split through its chain by the edit tools still meshes
    write_json(
        &blocks.join("split-pipe.json"),
        &json!({"topology": "pipe.json", "output": "pipe-split.json", "block": 0, "axis": 2, "at": 15}),
    )?;
    let envelope = run_cli(
        cli,
        &blocks,
        &["call", "block.split_block", "--input", "split-pipe.json"],
    )?;
    require(
        envelope["status"] == "finished",
        &format!("block.split_block: {envelope}"),
    )?;
    write_json(
        &blocks.join("grid-pipe-split.json"),
        &json!({"topology": "pipe-split.json", "output": "pipe-split", "algorithm": {"kind": "elliptic"}}),
    )?;
    run_cli(
        cli,
        &blocks,
        &["call", "block.grid", "--input", "grid-pipe-split.json"],
    )?;
    run(
        &blocks,
        "pipe-split/foehncraft.vtu",
        &blocks.join("pipe-split/mesh-evidence.json"),
        foehn::Kind::Faces,
    )?;
    write_json(&blocks.join("quality.json"), &block_quality)?;
    // BC-01, BC-02: a boundary defined with the patch tools reaches
    // FoehnCraft under its names in every format
    let patched = bundles.join("patched-box");
    fs::create_dir(&patched)?;
    let stl = path("unit-box.stl");
    let steps = [
        ("patch.init", json!({"surface": stl, "output": "p0.json"})),
        (
            "patch.select_plane",
            json!({"surface": stl, "patches": "p0.json", "output": "p1.json",
            "name": "inlet", "role": "inlet", "point": [0, 0, 0], "normal": [1, 0, 0]}),
        ),
        (
            "patch.select_normal",
            json!({"surface": stl, "patches": "p1.json", "output": "p2.json",
            "name": "outlet", "role": "outlet", "direction": [1, 0, 0], "max_angle_deg": 1}),
        ),
        (
            "patch.merge",
            json!({"surface": stl, "patches": "p2.json", "output": "p3.json",
            "merge": ["south_wall", "north_wall", "floor", "ceiling"], "name": "walls", "role": "wall"}),
        ),
    ];
    for (k, (id, input)) in steps.iter().enumerate() {
        let file = format!("in{k}.json");
        write_json(&patched.join(&file), input)?;
        let envelope = run_cli(cli, &patched, &["call", id, "--input", &file])?;
        require(
            envelope["status"] == "finished",
            &format!("{id}: {envelope}"),
        )?;
    }
    run_cli(
        cli,
        &patched,
        &[
            "generate-stl",
            &stl,
            &path("stl-size-field.json"),
            "m.json",
            "--patches",
            "p3.json",
        ],
    )?;
    for (format, file) in [
        ("openfoam", "foam"),
        ("foehncraft", "m.vtu"),
        ("foehncraft-vtkhdf", "m.vtkhdf"),
        ("gmsh", "m.msh"),
    ] {
        run_cli(cli, &patched, &["export", "m.json", format, file])?;
    }
    let evidence = patched.join("foam/mesh-evidence.json");
    let written: Value = serde_json::from_slice(&fs::read(&evidence)?)?;
    let mut names: Vec<&str> = written["patches"]
        .as_array()
        .ok_or("missing patches in the evidence")?
        .iter()
        .filter_map(|p| p["name"].as_str())
        .collect();
    names.sort_unstable();
    require(
        names == ["inlet", "outlet", "walls"],
        &format!("patched box exports patches {names:?}"),
    )?;
    run(&bundles, "patched-box/m.vtu", &evidence, foehn::Kind::Faces)?;
    run(
        &bundles,
        "patched-box/m.vtkhdf",
        &evidence,
        foehn::Kind::Faces,
    )?;
    run(&bundles, "patched-box/foam", &evidence, foehn::Kind::Faces)?;
    run(
        &bundles,
        "patched-box/m.msh",
        &evidence,
        foehn::Kind::Subdivision,
    )?;
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
    use std::io::BufRead;
    let mut server = Command::new(mcp)
        .arg("--read-root")
        .arg(directory)
        .arg("--write-root")
        .arg(directory)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;
    let mut input = server.stdin.take().ok_or("no MCP stdin")?;
    let mut output = io::BufReader::new(server.stdout.take().ok_or("no MCP stdout")?);
    let mut next_id = 0;
    let mut send = |method: &str, params: Value, notify: bool| -> RunResult<Value> {
        let mut request = json!({ "jsonrpc": "2.0", "method": method, "params": params });
        if notify {
            writeln!(input, "{request}")?;
            input.flush()?;
            return Ok(Value::Null);
        }
        request["id"] = json!(next_id);
        next_id += 1;
        writeln!(input, "{request}")?;
        input.flush()?;
        let mut line = String::new();
        output.read_line(&mut line)?;
        let response: Value = serde_json::from_str(&line)?;
        require(response["id"] == request["id"], "MCP response out of order")?;
        Ok(response)
    };
    let initialized = send(
        "initialize",
        json!({
            "protocolVersion": "2025-11-25", "capabilities": {},
            "clientInfo": { "name": "mesh-acceptance", "version": "1" }
        }),
        false,
    )?;
    require(
        initialized["result"]["protocolVersion"] == "2025-11-25",
        "MCP initialization failed",
    )?;
    send("notifications/initialized", json!({}), true)?;
    let mut call = |name: &str, arguments: Value| -> RunResult<Value> {
        let response = send(
            "tools/call",
            json!({ "name": name, "arguments": arguments }),
            false,
        )?;
        require(
            response["result"]["isError"] == false,
            &format!("MCP tool call failed: {response}"),
        )?;
        let envelope = &response["result"]["structuredContent"];
        require(
            envelope["format"] == "koflcraft.tool-result" && envelope["status"] == "finished",
            &format!("MCP result is not a finished envelope: {envelope}"),
        )?;
        Ok(envelope["result"].clone())
    };
    let mut results = vec![
        call("domain.inspect", json!({ "config": config }))?,
        call(
            "mesh.generate",
            json!({ "config": config, "artifact": "mcp.mesh.json" }),
        )?,
        call("mesh.check", json!({ "artifact": "cli.mesh.json" }))?,
        call(
            "mesh.export",
            json!({ "artifact": "mcp.mesh.json", "format": "vtu", "output": "mcp.vtu" }),
        )?,
        call(
            "mesh.export",
            json!({ "artifact": "mcp.mesh.json", "format": "gmsh", "output": "mcp.msh" }),
        )?,
    ];
    // AG-16: the MCP report is a summary with pages; follow the cursors
    // and put the rows back together, so the comparison with the CLI's
    // full report covers the paging too.
    let mut report = call("mesh.report", json!({ "artifact": "mcp.mesh.json" }))?;
    let mut rows = json!({ "cells": [], "faces": [] });
    loop {
        let page = report["page"].take();
        let kind = page["rows"]
            .as_str()
            .ok_or("MCP report page without rows")?;
        let items = page["items"]
            .as_array()
            .ok_or("MCP report page without items")?;
        rows[kind]
            .as_array_mut()
            .ok_or("unknown row kind")?
            .extend(items.iter().cloned());
        let Some(cursor) = page["next_cursor"].as_str() else {
            break;
        };
        let next = call(
            "mesh.report",
            json!({ "artifact": "mcp.mesh.json", "cursor": cursor }),
        )?;
        report["page"] = next["page"].clone();
    }
    require(
        report["summary"]["cells"] == rows["cells"].as_array().map_or(0, Vec::len)
            && report["summary"]["faces"] == rows["faces"].as_array().map_or(0, Vec::len),
        "MCP report pages do not add up to the summary counts",
    )?;
    report["cells"] = rows["cells"].take();
    report["faces"] = rows["faces"].take();
    results.push(report);
    // closes the server's stdin; the closures that wrote to it are done
    drop(input);
    let output = server.wait_with_output()?;
    require(
        output.status.success(),
        &String::from_utf8_lossy(&output.stderr),
    )?;
    Ok(results)
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
