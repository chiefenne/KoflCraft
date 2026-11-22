use mesh_tools::{Action, Context, Journal, TOOLS, call_journaled};
use serde_json::{Value, json};

const PROTOCOL_VERSION: &str = "2025-11-25";

enum Phase {
    Initialize,
    InitializedNotification,
    Ready,
}

pub struct Session {
    ctx: Context,
    journal: Option<Journal>,
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
    pub fn new(ctx: Context, journal: Option<Journal>) -> Self {
        Self {
            ctx,
            journal,
            phase: Phase::Initialize,
        }
    }

    pub fn handle(&mut self, frame: &[u8]) -> Option<Value> {
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
                let envelope = call_journaled(name, args, &self.ctx, self.journal.as_ref());
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
