use mesh_tools::{Action, Context, Diagnostic, Envelope, Journal, TOOLS, call_journaled};
use serde_json::{Value, json};

use crate::policy::{Policy, Rule};

/// The client side of the connection, for requests the server sends while
/// it handles a call (elicitation, AG-15).
pub trait Peer {
    /// Send a request and return the client's response to it.
    fn request(&mut self, request: Value) -> std::io::Result<Value>;
}

const PROTOCOL_VERSION: &str = "2025-11-25";

enum Phase {
    Initialize,
    InitializedNotification,
    Ready,
}

pub struct Session {
    ctx: Context,
    journal: Option<Journal>,
    policy: Policy,
    /// The client declared the elicitation capability.
    can_ask: bool,
    asked: u64,
    phase: Phase,
}

/// The registry's tools as MCP tools: name = dotted id, output schema =
/// the result envelope around the tool's own output schema.
fn listing() -> Vec<Value> {
    TOOLS
        .iter()
        .map(|t| {
            let read = t.action == Action::Read;
            json!({
                "name": t.id,
                "title": t.id,
                "description": t.summary,
                "inputSchema": t.input_schema(),
                "outputSchema": envelope_schema(t.output_schema()),
                "annotations": {
                    "readOnlyHint": read, "destructiveHint": false,
                    "idempotentHint": read, "openWorldHint": false
                },
                "_meta": {
                    "koflcraft/version": t.version,
                    "koflcraft/action": t.action,
                    "koflcraft/legacy_name": t.legacy_name,
                    "koflcraft/output_fields": t.outputs,
                }
            })
        })
        .collect()
}

fn envelope_schema(result: Value) -> Value {
    json!({
        "type": "object",
        "required": ["format", "format_version", "tool", "tool_version", "status",
                     "exit_code", "result", "diagnostics", "limitations"],
        "properties": {
            "format": {"const": mesh_tools::RESULT_FORMAT},
            "format_version": {"type": "string"},
            "tool": {"type": "string"},
            "tool_version": {"type": "integer"},
            "status": {"enum": ["finished", "error"]},
            "exit_code": {"type": "integer"},
            "result": {"anyOf": [result, {"type": "null"}]},
            "diagnostics": {"type": "array", "items": {"type": "object"}},
            "limitations": {"type": "array", "items": {"type": "string"}}
        }
    })
}

impl Session {
    pub fn new(ctx: Context, journal: Option<Journal>, policy: Policy) -> Self {
        let ctx = Context {
            limits: Some(policy.limits.clone()),
            ..ctx
        };
        Self {
            ctx,
            journal,
            policy,
            can_ask: false,
            asked: 0,
            phase: Phase::Initialize,
        }
    }

    /// Ask the person through the client (MCP elicitation) whether this
    /// build may start. `Ok(true)` only for an explicit accept with
    /// `approve: true`.
    fn ask(&mut self, peer: &mut dyn Peer, message: String) -> Result<bool, Diagnostic> {
        if !self.can_ask {
            return Err(Diagnostic::new("approval_required", message));
        }
        self.asked += 1;
        let response = peer
            .request(json!({
                "jsonrpc": "2.0",
                "id": format!("koflcraft-approval-{}", self.asked),
                "method": "elicitation/create",
                "params": {
                    "message": message,
                    "requestedSchema": {
                        "type": "object",
                        "properties": {
                            "approve": {"type": "boolean", "title": "Start this build",
                                         "description": "Run it once with the limits lifted"}
                        },
                        "required": ["approve"]
                    }
                }
            }))
            .map_err(|e| Diagnostic::new("approval_required", format!("{message} ({e})")))?;
        let result = &response["result"];
        Ok(result["action"] == "accept" && result["content"]["approve"] == true)
    }

    /// A tools/call under the policy: estimate a build, then allow, ask or
    /// refuse; a build that hits an input limit counts as over the limits.
    fn call(&mut self, name: &str, args: Value, peer: &mut dyn Peer) -> Envelope {
        let Some(tool) = mesh_tools::find(name) else {
            return call_journaled(name, args, &self.ctx, self.journal.as_ref());
        };
        if tool.action != Action::StartBuild {
            return call_journaled(name, args, &self.ctx, self.journal.as_ref());
        }
        let refuse = |d: Diagnostic| Envelope::error(tool.id, tool.version, d);
        // the estimate without the input limits, so it is not cut short
        let open = Context {
            limits: None,
            ..self.ctx.clone()
        };
        let estimate = mesh_tools::call(
            "mesh.estimate",
            json!({ "tool": tool.id, "input": args }),
            &open,
        );
        let limits = &self.policy.limits;
        let (over, what) = match estimate.result.as_object() {
            Some(e) => {
                let cells = e["cells"].as_u64().unwrap_or(0);
                let bytes = e["output_bytes"].as_u64().unwrap_or(0);
                (
                    cells > limits.max_build_cells as u64 || bytes > limits.max_output_bytes,
                    format!(
                        "{} would build {cells} cells and write about {} MB (limits: {} cells, {} MB)",
                        tool.id,
                        bytes / 1_000_000,
                        limits.max_build_cells,
                        limits.max_output_bytes / 1_000_000
                    ),
                )
            }
            // a failing estimate: the call itself reports the reason
            None => (false, format!("{} (no estimate)", tool.id)),
        };
        if !over {
            match self.policy.start_build {
                Rule::Refuse => {
                    return refuse(Diagnostic::new(
                        "policy_refused",
                        format!("start_build is refused: {what}"),
                    ));
                }
                Rule::Ask => match self.ask(peer, format!("Start this build? {what}")) {
                    Ok(true) => {}
                    Ok(false) => return refuse(Diagnostic::new("approval_declined", what)),
                    Err(d) => return refuse(d),
                },
                Rule::Allow | Rule::Default => {}
            }
            let envelope = call_journaled(name, args.clone(), &self.ctx, self.journal.as_ref());
            let hit_limit = envelope
                .diagnostics
                .first()
                .is_some_and(|d| d.code == "TOOL008");
            if !hit_limit {
                return envelope;
            }
            let reason = envelope.diagnostics[0].message.clone();
            return self.exceed(name, args, peer, format!("{what}; {reason}"));
        }
        self.exceed(name, args, peer, what)
    }

    fn exceed(&mut self, name: &str, args: Value, peer: &mut dyn Peer, what: String) -> Envelope {
        let tool = mesh_tools::find(name);
        let (id, version) = tool.map_or((name, 0), |t| (t.id, t.version));
        let refuse = |d: Diagnostic| Envelope::error(id, version, d);
        match self.policy.exceed_limits {
            Rule::Refuse => {
                return refuse(Diagnostic::new(
                    "limit_exceeded",
                    format!("over the MCP limits and exceed_limits is refuse: {what}"),
                ));
            }
            Rule::Ask => match self.ask(
                peer,
                format!("This build is over the MCP limits. Start it anyway? {what}"),
            ) {
                Ok(true) => {}
                Ok(false) => return refuse(Diagnostic::new("approval_declined", what)),
                Err(d) => return refuse(d),
            },
            Rule::Allow | Rule::Default => {}
        }
        // approved or allowed: once, with the limits lifted; the roots stay
        let open = Context {
            limits: None,
            ..self.ctx.clone()
        };
        let mut envelope = call_journaled(name, args, &open, self.journal.as_ref());
        envelope
            .limitations
            .push(format!("built over the MCP limits with approval: {what}"));
        envelope
    }

    pub fn handle(&mut self, frame: &[u8], peer: &mut dyn Peer) -> Option<Value> {
        let request: Value = match serde_json::from_slice(frame) {
            Ok(value) => value,
            Err(_) => return Some(rpc_error(Value::Null, -32700, "invalid JSON")),
        };
        let id = request.get("id").cloned();
        let method = request.get("method").and_then(Value::as_str);
        if request.get("jsonrpc").and_then(Value::as_str) != Some("2.0")
            || method.is_none()
            || id
                .as_ref()
                .is_some_and(|id| !(id.is_string() || id.is_i64() || id.is_u64()))
        {
            return Some(rpc_error(Value::Null, -32600, "invalid JSON-RPC request"));
        }
        let method = method.unwrap();
        let Some(id) = id else {
            if method == "notifications/initialized"
                && matches!(self.phase, Phase::InitializedNotification)
            {
                self.phase = Phase::Ready;
            }
            return None;
        };
        let params = request.get("params").unwrap_or(&Value::Null);
        let result = match method {
            "initialize" if matches!(self.phase, Phase::Initialize) => {
                let Some(version) = params.get("protocolVersion").and_then(Value::as_str) else {
                    return Some(rpc_error(id, -32602, "initialize requires protocolVersion"));
                };
                let client = &params["clientInfo"];
                if !params["capabilities"].is_object()
                    || !client["name"].is_string()
                    || !client["version"].is_string()
                {
                    return Some(rpc_error(
                        id,
                        -32602,
                        "initialize requires capabilities and clientInfo",
                    ));
                }
                let protocol_version = match version {
                    "2025-06-18" => version,
                    _ => PROTOCOL_VERSION,
                };
                self.can_ask = params["capabilities"]["elicitation"].is_object();
                self.phase = Phase::InitializedNotification;
                json!({
                    "protocolVersion": protocol_version,
                    "capabilities": { "tools": { "listChanged": false } },
                    "serverInfo": { "name": "koflcraft", "version": env!("CARGO_PKG_VERSION") }
                })
            }
            "initialize" => return Some(rpc_error(id, -32600, "server is already initialized")),
            "ping" => json!({}),
            _ if !matches!(self.phase, Phase::Ready) => {
                return Some(rpc_error(id, -32002, "server is not initialized"));
            }
            "tools/list" => json!({ "tools": listing() }),
            "tools/call" => {
                let Some(name) = params.get("name").and_then(Value::as_str) else {
                    return Some(rpc_error(id, -32602, "tools/call requires a tool name"));
                };
                let args = params
                    .get("arguments")
                    .cloned()
                    .unwrap_or_else(|| json!({}));
                if mesh_tools::find(name).is_none() {
                    return Some(rpc_error(id, -32602, &format!("unknown tool {name}")));
                }
                let envelope = self.call(name, args, peer);
                let body = serde_json::to_value(&envelope).unwrap_or(Value::Null);
                json!({
                    "content": [{ "type": "text", "text": body.to_string() }],
                    "structuredContent": body,
                    "isError": !envelope.ok()
                })
            }
            _ => return Some(rpc_error(id, -32601, "method not found")),
        };
        Some(json!({ "jsonrpc": "2.0", "id": id, "result": result }))
    }
}

pub fn rpc_error(id: Value, code: i32, message: &str) -> Value {
    json!({ "jsonrpc": "2.0", "id": id, "error": { "code": code, "message": message } })
}
