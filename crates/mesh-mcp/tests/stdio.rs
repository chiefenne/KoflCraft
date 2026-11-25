use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use mesh_engine::{ExportData, ExportFormat, GmshVersion, build_mesh, export_mesh};
use serde_json::{Value, json};

#[test]
fn stdio_tools_run_the_shared_engine_and_bound_writes() {
    let work = WorkDirectory::new();
    let config: Value =
        serde_json::from_slice(include_bytes!("../../../test-fixtures/box-4.json")).unwrap();
    let channel: Value = serde_json::from_slice(include_bytes!(
        "../../../test-fixtures/straight-channel.json"
    ))
    .unwrap();
    let elbow = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test-fixtures/elbow");
    for file in [
        "elbow.stl",
        "patches.json",
        "partition.json",
        "volume-scaffold.msh",
    ] {
        fs::copy(elbow.join(file), work.path().join(file)).unwrap();
    }
    let mut oversized = channel.clone();
    oversized["layers"] = json!(usize::MAX);
    let requests = [
        rpc(0, "tools/list", json!({})),
        rpc(
            1,
            "initialize",
            json!({
                "protocolVersion": "2025-11-25",
                "capabilities": {},
                "clientInfo": { "name": "test", "version": "1" }
            }),
        ),
        json!({ "jsonrpc": "2.0", "method": "notifications/initialized" }),
        rpc(2, "tools/list", json!({})),
        call(3, "domain.inspect", json!({ "config": config })),
        call(
            4,
            "mesh.generate",
            json!({ "config": config, "artifact": "mesh.json" }),
        ),
        call(5, "mesh.check", json!({ "artifact": "mesh.json" })),
        call(6, "mesh_report", json!({ "artifact": "mesh.json" })),
        call(
            7,
            "export_mesh",
            json!({ "artifact": "mesh.json", "format": "vtu", "output": "mesh.vtu" }),
        ),
        call(
            8,
            "export_mesh",
            json!({ "artifact": "mesh.json", "format": "gmsh", "output": "mesh.msh" }),
        ),
        call(
            9,
            "export_mesh",
            json!({ "artifact": "mesh.json", "format": "vtu", "output": "../escaped.vtu" }),
        ),
        call(
            10,
            "export_mesh",
            json!({ "artifact": "mesh.json", "format": "vtu", "output": "mesh.vtu" }),
        ),
        call(
            11,
            "generate_channel",
            json!({"config": channel, "output":"duct"}),
        ),
        call(
            12,
            "generate_channel",
            json!({"config": oversized, "output":"over"}),
        ),
        call(
            13,
            "generate_duct",
            json!({"surface":"elbow.stl","patches":"patches.json","partition":"partition.json","mode":"layered","output":"elbow"}),
        ),
        call(
            14,
            "generate_duct",
            json!({"surface":"elbow.stl","patches":"patches.json","partition":"partition.json","mode":"polyhedral","output":"elbow-poly"}),
        ),
        call(
            15,
            "generate_sweep",
            json!({"config":{"schema_version":1,"geometry":{"kind":"pipe","diameter_mm":12,"length_mm":6},"core_spacing_mm":1.25,"wall_layers":2,"first_layer_mm":0.55,"layer_growth":1.2,"max_neighbour_size_ratio":1.4},"output":"sweep"}),
        ),
    ];
    let responses = run(work.path(), &requests);
    assert_eq!(responses.len(), requests.len() - 1);
    for (index, response) in responses.iter().enumerate() {
        assert_eq!(response["id"], index);
    }
    assert_eq!(responses[0]["error"]["code"], -32002);
    assert_eq!(responses[1]["result"]["protocolVersion"], "2025-11-25");
    assert_eq!(
        responses[2]["result"]["tools"].as_array().unwrap().len(),
        mesh_tools::TOOLS.len()
    );
    let names: Vec<&str> = responses[2]["result"]["tools"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| t["name"].as_str().unwrap())
        .collect();
    assert!(names.contains(&"mesh.generate") && names.contains(&"journal.replay"));
    assert_eq!(
        responses[3]["result"]["structuredContent"]["result"]["site_count"],
        4
    );
    let build = &responses[4]["result"]["structuredContent"]["result"]["build"];
    assert_eq!(
        responses[5]["result"]["structuredContent"]["result"]["valid"],
        true
    );
    assert_eq!(
        responses[6]["result"]["structuredContent"]["result"]["build"],
        *build
    );
    // AG-16: over MCP the report is a summary and one page, never the
    // full per-cell and per-face arrays
    let report = &responses[6]["result"]["structuredContent"]["result"];
    assert!(report.get("cells").is_none() && report.get("faces").is_none());
    assert_eq!(report["summary"]["cells"], 4);
    assert_eq!(report["page"]["rows"], "cells");
    assert_eq!(report["page"]["count"], 4);
    assert!(report["page"]["next_cursor"].is_string());
    for response in &responses[7..=8] {
        assert_eq!(response["result"]["isError"], false);
    }
    for (response, expected) in responses[9..=10].iter().zip(["invalid_path", "io_error"]) {
        assert_eq!(response["result"]["isError"], true);
        let text = response["result"]["content"][0]["text"].as_str().unwrap();
        let envelope: Value = serde_json::from_str(text).unwrap();
        assert_eq!(envelope, response["result"]["structuredContent"]);
        assert_eq!(envelope["status"], "error");
        assert_eq!(envelope["diagnostics"][0]["name"], expected);
    }
    assert!(!work.path().parent().unwrap().join("escaped.vtu").exists());

    assert_eq!(
        responses[11]["result"]["structuredContent"]["result"]["report"]["prisms"],
        768
    );
    assert!(work.path().join("duct/inlet.vtu").is_file());
    assert_eq!(responses[12]["result"]["isError"], true);
    assert!(!work.path().join("over").exists());

    assert_eq!(
        responses[13]["result"]["structuredContent"]["result"]["report"]["prisms"],
        120
    );
    assert!(work.path().join("elbow/mesh.vtm").is_file());
    let poly = &responses[14]["result"]["structuredContent"]["result"]["report"];
    assert_eq!(poly["tetrahedra"], 0);
    // median dual core (default since 2026-11-25) and the polygon prisms
    // of its wall layers
    assert_eq!(poly["polyhedra"], 238 + 816);
    assert_eq!(poly["layer_polyhedra"], 816);
    assert!(
        responses[15]["result"]["structuredContent"]["result"]["report"]["core_cells"]
            .as_u64()
            .unwrap()
            > 0
    );
    let local = build_mesh(serde_json::from_value(config).unwrap()).unwrap();
    assert_eq!(build["topology_digest"], local.build().topology_digest);
    assert_eq!(build["mesh_digest"], local.build().mesh_digest);
    let ExportData::Vtu(expected) = export_mesh(&local, ExportFormat::Vtu).unwrap() else {
        panic!("expected VTU data");
    };
    assert_eq!(fs::read(work.path().join("mesh.vtu")).unwrap(), expected);
    let ExportData::Gmsh(expected) =
        export_mesh(&local, ExportFormat::Gmsh(GmshVersion::V4_1)).unwrap()
    else {
        panic!("expected Gmsh data");
    };
    assert_eq!(fs::read(work.path().join("mesh.msh")).unwrap(), expected);
}

#[cfg(unix)]
#[test]
fn external_symlinks_cannot_be_used_as_artifacts_or_outputs() {
    use std::os::unix::fs::symlink;

    let outside = WorkDirectory::new();
    let inside = WorkDirectory::new();
    let config: Value =
        serde_json::from_slice(include_bytes!("../../../test-fixtures/box-4.json")).unwrap();
    let path = outside.path().join("existing.json");
    fs::write(&path, b"keep").unwrap();
    symlink(&path, inside.path().join("link.json")).unwrap();
    let responses = run(
        inside.path(),
        &[
            rpc(
                1,
                "initialize",
                json!({
                    "protocolVersion": "2025-06-18",
                    "capabilities": {},
                    "clientInfo": { "name": "test", "version": "1" }
                }),
            ),
            json!({ "jsonrpc": "2.0", "method": "notifications/initialized" }),
            call(2, "check_mesh", json!({ "artifact": "link.json" })),
            call(
                3,
                "generate_mesh",
                json!({ "config": config, "artifact": "link.json" }),
            ),
        ],
    );
    assert_eq!(responses[0]["result"]["protocolVersion"], "2025-06-18");
    assert!(
        responses[1..]
            .iter()
            .all(|response| response["result"]["isError"] == true)
    );
    assert_eq!(fs::read(&path).unwrap(), b"keep");
}

#[test]
fn stl_tools_reuse_check_and_export_for_a_skewed_surface() {
    let work = WorkDirectory::new();
    let source = include_bytes!("../../../test-fixtures/skewed-box.stl");
    fs::write(work.path().join("skew.stl"), source).unwrap();
    let config: Value =
        serde_json::from_slice(include_bytes!("../../../test-fixtures/stl-sites.json")).unwrap();
    let responses = run(
        work.path(),
        &[
            rpc(
                0,
                "initialize",
                json!({
                    "protocolVersion": "2025-11-25",
                    "capabilities": {},
                    "clientInfo": { "name": "test", "version": "1" }
                }),
            ),
            json!({ "jsonrpc": "2.0", "method": "notifications/initialized" }),
            call(
                1,
                "inspect_stl_domain",
                json!({ "surface": "skew.stl", "config": config }),
            ),
            call(
                2,
                "generate_stl_mesh",
                json!({ "surface": "skew.stl", "config": config, "artifact": "stl.mesh.json" }),
            ),
            call(3, "check_mesh", json!({ "artifact": "stl.mesh.json" })),
            call(
                4,
                "export_mesh",
                json!({ "artifact": "stl.mesh.json", "format": "gmsh", "output": "stl.msh" }),
            ),
        ],
    );
    for response in &responses[1..] {
        assert_eq!(response["result"]["isError"], false, "{response}");
    }
    assert_eq!(
        responses[1]["result"]["structuredContent"]["result"]["patches"][0]["name"],
        "west inlet"
    );
    let local =
        mesh_engine::build_stl_mesh(source, serde_json::from_value(config).unwrap(), None, None)
            .unwrap();
    let actual = &responses[2]["result"]["structuredContent"]["result"];
    assert_eq!(actual["artifact_schema_version"], 2);
    assert_eq!(actual["build"]["mesh_digest"], local.build().mesh_digest);
    assert_eq!(
        responses[3]["result"]["structuredContent"]["result"]["valid"],
        true
    );
    let ExportData::Gmsh(expected) =
        export_mesh(&local, ExportFormat::Gmsh(GmshVersion::V4_1)).unwrap()
    else {
        panic!("expected Gmsh data");
    };
    assert_eq!(fs::read(work.path().join("stl.msh")).unwrap(), expected);
}

fn run(directory: &Path, requests: &[Value]) -> Vec<Value> {
    run_with(&["--output-dir".as_ref(), directory.as_os_str()], requests)
}

fn run_with(args: &[&std::ffi::OsStr], requests: &[Value]) -> Vec<Value> {
    let mut server = Command::new(env!("CARGO_BIN_EXE_mesh-mcp"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    {
        let mut input = server.stdin.take().unwrap();
        for request in requests {
            writeln!(input, "{request}").unwrap();
        }
    }
    let output = server.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect()
}

const INIT: fn() -> [Value; 2] = || {
    [
        rpc(
            0,
            "initialize",
            json!({
                "protocolVersion": "2025-11-25",
                "capabilities": {},
                "clientInfo": { "name": "test", "version": "1" }
            }),
        ),
        json!({ "jsonrpc": "2.0", "method": "notifications/initialized" }),
    ]
};

/// Invariant (AG-14): files are read only below a read root and written
/// only below a write root, and the server does not start without both.
/// A context that checked only the write side, or a server that fell back
/// to the current directory, breaks it.
#[test]
fn read_and_write_roots_confine_every_call() {
    let input = WorkDirectory::new();
    let output = WorkDirectory::new();
    let elsewhere = WorkDirectory::new();
    let source = include_bytes!("../../../test-fixtures/skewed-box.stl");
    fs::write(input.path().join("skew.stl"), source).unwrap();
    fs::write(elsewhere.path().join("skew.stl"), source).unwrap();
    let config: Value =
        serde_json::from_slice(include_bytes!("../../../test-fixtures/stl-sites.json")).unwrap();
    let inside = input.path().join("skew.stl");
    let outside = elsewhere.path().join("skew.stl");
    let mut requests = INIT().to_vec();
    requests.extend([
        call(
            1,
            "mesh.generate_stl",
            json!({ "surface": inside, "config": config, "artifact": "m.json" }),
        ),
        call(
            2,
            "domain.inspect_stl",
            json!({ "surface": outside, "config": config }),
        ),
        call(
            3,
            "mesh.export",
            json!({ "artifact": "m.json", "format": "vtu", "output": input.path().join("m.vtu") }),
        ),
        call(4, "mesh.check", json!({ "artifact": "m.json" })),
    ]);
    let responses = run_with(
        &[
            "--read-root".as_ref(),
            input.path().as_os_str(),
            "--read-root".as_ref(),
            output.path().as_os_str(),
            "--write-root".as_ref(),
            output.path().as_os_str(),
        ],
        &requests,
    );
    let envelope = |i: usize| &responses[i]["result"]["structuredContent"];
    assert_eq!(envelope(1)["status"], "finished", "{}", envelope(1));
    assert!(output.path().join("m.json").is_file());
    assert_eq!(envelope(2)["diagnostics"][0]["code"], "TOOL006");
    assert_eq!(envelope(3)["diagnostics"][0]["code"], "TOOL006");
    assert!(!input.path().join("m.vtu").exists());
    assert_eq!(envelope(4)["result"]["valid"], true);

    let refused = Command::new(env!("CARGO_BIN_EXE_mesh-mcp"))
        .arg("--write-root")
        .arg(output.path())
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(!refused.status.success());
}

fn call(id: u32, name: &str, args: Value) -> Value {
    rpc(id, "tools/call", json!({ "name": name, "arguments": args }))
}

fn rpc(id: u32, method: &str, params: Value) -> Value {
    json!({ "jsonrpc": "2.0", "id": id, "method": method, "params": params })
}

struct WorkDirectory(PathBuf);

impl WorkDirectory {
    fn new() -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "koflcraft-mcp-{}-{stamp}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&dir).unwrap();
        Self(dir)
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for WorkDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(self.path());
    }
}

fn init_with(elicitation: bool) -> [Value; 2] {
    let capabilities = if elicitation {
        json!({ "elicitation": {} })
    } else {
        json!({})
    };
    [
        rpc(
            0,
            "initialize",
            json!({
                "protocolVersion": "2025-11-25",
                "capabilities": capabilities,
                "clientInfo": { "name": "test", "version": "1" }
            }),
        ),
        json!({ "jsonrpc": "2.0", "method": "notifications/initialized" }),
    ]
}

fn approval(n: u32, approve: bool) -> Value {
    json!({ "jsonrpc": "2.0", "id": format!("koflcraft-approval-{n}"),
            "result": { "action": "accept", "content": { "approve": approve } } })
}

/// Invariant (AG-15): with the decided defaults a build within the limits
/// runs without a prompt; a build over them (by estimate, or by an input
/// limit) asks through MCP elicitation and runs only on an explicit
/// approval; without elicitation it fails with TOOL009, a decline with
/// TOOL010, a `refuse` rule with TOOL011, and nothing is written in those
/// cases. A policy file inside a write root is refused at start.
/// Dropping the estimate check, or treating a missing answer as approval,
/// breaks it; no other test drives the approval path.
#[test]
fn builds_over_the_limits_need_approval() {
    let work = WorkDirectory::new();
    let policies = WorkDirectory::new();
    let config: Value =
        serde_json::from_slice(include_bytes!("../../../test-fixtures/box-4.json")).unwrap();
    let policy = |name: &str, text: &str| {
        let path = policies.path().join(name);
        fs::write(&path, text).unwrap();
        path
    };
    let small = policy("small.toml", "[limits]\nmax_build_cells = 2\n");
    let sites = policy("sites.toml", "[limits]\nmax_sites = 2\n");
    let refuse = policy("refuse.toml", "[policy]\nstart_build = \"refuse\"\n");
    let generate = |id: u32, artifact: &str| {
        call(
            id,
            "mesh.generate",
            json!({ "config": config, "artifact": artifact }),
        )
    };
    let serve = |policy: Option<&Path>, requests: &[Value]| {
        let mut args: Vec<&std::ffi::OsStr> = vec![
            "--read-root".as_ref(),
            work.path().as_os_str(),
            "--write-root".as_ref(),
            work.path().as_os_str(),
        ];
        if let Some(p) = policy {
            args.extend(["--policy".as_ref(), p.as_os_str()]);
        }
        run_with(&args, requests)
    };
    let envelope = |r: &Value| r["result"]["structuredContent"].clone();

    // defaults: four cells are within the limits, no prompt
    let mut requests = init_with(true).to_vec();
    requests.push(generate(1, "a.json"));
    let r = serve(None, &requests);
    assert_eq!(r.len(), 2, "no elicitation request expected: {r:?}");
    assert_eq!(envelope(&r[1])["status"], "finished");

    // over max_build_cells, client cannot ask: TOOL009
    let mut requests = init_with(false).to_vec();
    requests.push(generate(1, "b.json"));
    let r = serve(Some(&small), &requests);
    assert_eq!(envelope(&r[1])["diagnostics"][0]["code"], "TOOL009");
    assert!(!work.path().join("b.json").exists());

    // over max_build_cells, the person approves: built, with a limitation
    let mut requests = init_with(true).to_vec();
    requests.extend([generate(1, "c.json"), approval(1, true)]);
    let r = serve(Some(&small), &requests);
    assert_eq!(r[1]["method"], "elicitation/create");
    assert!(
        r[1]["params"]["message"]
            .as_str()
            .unwrap()
            .contains("4 cells")
    );
    assert_eq!(envelope(&r[2])["status"], "finished", "{}", r[2]);
    assert!(work.path().join("c.json").is_file());
    assert!(
        envelope(&r[2])["limitations"]
            .to_string()
            .contains("with approval")
    );

    // the person declines: TOOL010, nothing written
    let mut requests = init_with(true).to_vec();
    requests.extend([generate(1, "d.json"), approval(1, false)]);
    let r = serve(Some(&small), &requests);
    assert_eq!(envelope(&r[2])["diagnostics"][0]["code"], "TOOL010");
    assert!(!work.path().join("d.json").exists());

    // an input limit (sites) counts as over the limits too
    let mut requests = init_with(false).to_vec();
    requests.push(generate(1, "e.json"));
    let r = serve(Some(&sites), &requests);
    assert_eq!(envelope(&r[1])["diagnostics"][0]["code"], "TOOL009");
    assert!(!work.path().join("e.json").exists());

    // a refuse rule
    let mut requests = init_with(true).to_vec();
    requests.push(generate(1, "f.json"));
    let r = serve(Some(&refuse), &requests);
    assert_eq!(envelope(&r[1])["diagnostics"][0]["code"], "TOOL011");
    assert!(!work.path().join("f.json").exists());

    // a policy file inside the write root is refused at start
    fs::write(work.path().join("kofl-mcp.toml"), "").unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_mesh-mcp"))
        .args([
            "--read-root".as_ref(),
            work.path().as_os_str(),
            "--write-root".as_ref(),
            work.path().as_os_str(),
            "--policy".as_ref(),
        ])
        .arg(work.path().join("kofl-mcp.toml"))
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(!out.status.success());
}
