//! Journal and replay (AG-13). Every top-level tool call can be appended to
//! a journal (`koflcraft.journal`, JSON lines): a header, then one row per
//! call with its input, status, exit code, duration and a digest of its
//! result without the timings. `journal.replay` runs the calls of a journal
//! again and reports each one that ends differently.

use std::fs::OpenOptions;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};

use crate::diag::Diagnostic;
use crate::{Context, Envelope, FORMAT_VERSION, call};

pub const JOURNAL_FORMAT: &str = "koflcraft.journal";
pub const REPLAY_FORMAT: &str = "koflcraft.replay";

/// Keys whose values depend on the run, not on the inputs.
const VOLATILE: &[&str] = &["elapsed_micros", "generation_seconds", "seconds"];

fn strip(value: &Value) -> Value {
    match value {
        Value::Object(map) => Value::Object(
            map.iter()
                .filter(|(k, _)| !VOLATILE.contains(&k.as_str()))
                .map(|(k, v)| (k.clone(), strip(v)))
                .collect(),
        ),
        Value::Array(a) => Value::Array(a.iter().map(strip).collect()),
        v => v.clone(),
    }
}

/// SHA-256 of a result with the timings removed (keys are sorted, so the
/// serialization is canonical).
pub fn result_digest(result: &Value) -> String {
    let bytes = serde_json::to_vec(&strip(result)).unwrap_or_default();
    hex(&Sha256::digest(&bytes))
}

pub(crate) fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// The date of a row: `SOURCE_DATE_EPOCH` when set (SW-20), else now, in
/// seconds since 1970.
fn date() -> u64 {
    std::env::var("SOURCE_DATE_EPOCH")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or_else(|| {
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |d| d.as_secs())
        })
}

/// One call row.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Row {
    #[serde(rename = "type")]
    pub kind: String,
    pub seq: u64,
    pub date: u64,
    pub via: String,
    pub tool: String,
    pub tool_version: u32,
    pub input: Value,
    pub status: String,
    pub exit_code: i32,
    pub seconds: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub result_digest: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub diagnostic: Option<Value>,
}

/// An open journal file.
pub struct Journal {
    path: PathBuf,
    seq: Mutex<u64>,
}

impl Journal {
    /// Open (or create) a journal; a new file gets the header row.
    pub fn open(path: &Path) -> std::io::Result<Self> {
        let existing = std::fs::read_to_string(path).unwrap_or_default();
        let seq = existing
            .lines()
            .filter_map(|l| serde_json::from_str::<Value>(l).ok())
            .filter_map(|v| v["seq"].as_u64())
            .max()
            .unwrap_or(0);
        if existing.trim().is_empty() {
            let mut f = OpenOptions::new().create(true).append(true).open(path)?;
            writeln!(
                f,
                "{}",
                json!({
                    "type": "header", "format": JOURNAL_FORMAT, "format_version": FORMAT_VERSION,
                    "koflcraft_version": env!("CARGO_PKG_VERSION"),
                    "commit": option_env!("KOFL_COMMIT").unwrap_or("unknown"),
                })
            )?;
        }
        Ok(Self {
            path: path.to_path_buf(),
            seq: Mutex::new(seq),
        })
    }

    /// The journal of `--journal FILE` or `KOFL_JOURNAL`.
    pub fn from_env(flag: Option<&Path>) -> Option<std::io::Result<Self>> {
        let path = flag
            .map(Path::to_path_buf)
            .or_else(|| std::env::var_os("KOFL_JOURNAL").map(PathBuf::from))?;
        Some(Self::open(&path))
    }

    fn append(
        &self,
        id: &str,
        input: Value,
        envelope: &Envelope,
        via: &str,
    ) -> std::io::Result<()> {
        let mut seq = self
            .seq
            .lock()
            .map_err(|_| std::io::Error::other("journal lock"))?;
        *seq += 1;
        let row = Row {
            kind: "call".into(),
            seq: *seq,
            date: date(),
            via: via.into(),
            tool: if envelope.tool.is_empty() {
                id.into()
            } else {
                envelope.tool.clone()
            },
            tool_version: envelope.tool_version,
            input,
            status: envelope.status.clone(),
            exit_code: envelope.exit_code,
            seconds: envelope.seconds,
            result_digest: envelope.ok().then(|| result_digest(&envelope.result)),
            diagnostic: envelope
                .diagnostics
                .first()
                .filter(|_| !envelope.ok())
                .and_then(|d| serde_json::to_value(d).ok()),
        };
        let mut f = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&self.path)?;
        writeln!(f, "{}", serde_json::to_string(&row)?)
    }
}

/// Call a tool and append the call to `journal`. A journal that cannot be
/// written adds the warning TOOL003; the call stands.
pub fn call_journaled(
    id: &str,
    input: Value,
    ctx: &Context,
    journal: Option<&Journal>,
) -> Envelope {
    let recorded = input.clone();
    let mut envelope = call(id, input, ctx);
    if let Some(j) = journal
        && let Err(e) = j.append(id, recorded, &envelope, ctx.via)
    {
        envelope
            .diagnostics
            .push(Diagnostic::new("journal_failed", e));
    }
    envelope
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ReplayInput {
    /// The journal to replay.
    pub script: String,
    /// Go on past a call that ends differently from its record.
    #[serde(default)]
    pub keep_going: bool,
}

/// The call rows of a journal file.
pub(crate) fn calls(ctx: &Context, script: &str) -> Result<Vec<Row>, Diagnostic> {
    let text = String::from_utf8(ctx.read(script)?)
        .map_err(|e| Diagnostic::new("invalid_script", e).with_path(script))?;
    let mut rows = Vec::new();
    for (n, line) in text
        .lines()
        .enumerate()
        .filter(|(_, l)| !l.trim().is_empty())
    {
        let v: Value = serde_json::from_str(line).map_err(|e| {
            Diagnostic::new("invalid_script", format!("line {}: {e}", n + 1)).with_path(script)
        })?;
        if v["type"] == "call" {
            let row: Row = serde_json::from_value(v).map_err(|e| {
                Diagnostic::new("invalid_script", format!("line {}: {e}", n + 1)).with_path(script)
            })?;
            rows.push(row);
        }
    }
    if rows.is_empty() {
        return Err(Diagnostic::new("invalid_script", "the script has no calls").with_path(script));
    }
    Ok(rows)
}

/// A replay skips earlier replays and calls that never reached a tool.
pub(crate) fn skipped(row: &Row) -> bool {
    row.tool == "journal.replay"
        || row
            .diagnostic
            .as_ref()
            .is_some_and(|d| matches!(d["code"].as_str(), Some("TOOL001" | "TOOL002")))
}

/// `journal.replay`.
pub(crate) fn replay(input: Value, ctx: &Context) -> Result<Value, Diagnostic> {
    let i: ReplayInput = crate::mesh::decode(input)?;
    let rows = calls(ctx, &i.script)?;
    let inner = Context {
        via: "replay",
        ..ctx.clone()
    };
    let mut calls = Vec::new();
    let (mut replayed, mut skipped, mut different) = (0, 0, 0);
    for row in rows {
        if row.tool == "journal.replay"
            || row
                .diagnostic
                .as_ref()
                .is_some_and(|d| matches!(d["code"].as_str(), Some("TOOL001" | "TOOL002")))
        {
            skipped += 1;
            calls.push(json!({"seq": row.seq, "tool": row.tool, "skipped": true}));
            continue;
        }
        let envelope = call(&row.tool, row.input.clone(), &inner);
        replayed += 1;
        let digest = envelope.ok().then(|| result_digest(&envelope.result));
        let same = envelope.status == row.status
            && envelope.exit_code == row.exit_code
            && (row.result_digest.is_none() || digest == row.result_digest);
        calls.push(json!({
            "seq": row.seq, "tool": row.tool,
            "recorded": {"status": row.status, "exit_code": row.exit_code, "result_digest": row.result_digest},
            "replayed": {"status": envelope.status, "exit_code": envelope.exit_code, "result_digest": digest,
                          "diagnostics": envelope.diagnostics},
            "same": same,
        }));
        if !same {
            different += 1;
            if !i.keep_going {
                break;
            }
        }
    }
    let report = json!({
        "format": REPLAY_FORMAT, "format_version": FORMAT_VERSION,
        "calls": calls, "replayed": replayed, "skipped": skipped, "different": different,
    });
    if different > 0 {
        let mut d = Diagnostic::new(
            "replay_differs",
            format!("{different} replayed calls ended differently from their record"),
        );
        d.path = Some(i.script);
        return Err(d.with_report(report));
    }
    Ok(report)
}
