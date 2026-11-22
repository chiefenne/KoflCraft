//! The registry of diagnostic codes (AG-18). A code names one meaning
//! forever; it is never reused. `name` is the snake-case code KoflCraft
//! reported before the registry (kept in error objects for scripts that
//! read it).

use serde::Serialize;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Severity {
    Error,
    Warning,
    Note,
}

/// One registered code.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct Code {
    pub code: &'static str,
    pub name: &'static str,
    pub severity: Severity,
    pub meaning: &'static str,
    /// A suggested fix, where one can be given.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fix: Option<&'static str>,
}

const fn code(
    code: &'static str,
    name: &'static str,
    severity: Severity,
    meaning: &'static str,
    fix: Option<&'static str>,
) -> Code {
    Code {
        code,
        name,
        severity,
        meaning,
        fix,
    }
}

use Severity::{Error, Warning};

/// Every code, in the order of `diag.codes` and `docs/tools.md`.
pub const CODES: &[Code] = &[
    code(
        "TOOL001",
        "unknown_tool",
        Error,
        "No tool has this id.",
        Some("Run `kofl commands` for the ids."),
    ),
    code(
        "TOOL002",
        "invalid_arguments",
        Error,
        "The input does not match the tool's input schema.",
        Some("Compare the input with the schema in `kofl commands`."),
    ),
    code(
        "TOOL003",
        "journal_failed",
        Warning,
        "The journal could not be written: the call ran, but the journal misses it.",
        Some("Check that the journal path is writable."),
    ),
    code(
        "TOOL004",
        "invalid_script",
        Error,
        "A replay script is not JSON lines of journal calls, or has no calls.",
        None,
    ),
    code(
        "TOOL005",
        "replay_differs",
        Error,
        "A replayed call ended differently from its record (status, exit code or result digest).",
        None,
    ),
    code(
        "TOOL006",
        "invalid_path",
        Error,
        "A path is outside the read or write roots, or not a plain file name where one is needed.",
        Some("Use a path inside the roots given to `kofl mcp`."),
    ),
    code(
        "TOOL007",
        "usage",
        Error,
        "The command line is not one of the documented forms.",
        Some("Run `kofl --help`."),
    ),
    code(
        "TOOL008",
        "limit_exceeded",
        Error,
        "The input exceeds a limit of this front end (sites, facets, layers, cells, file size).",
        Some("Run the tool from the command line, which has no limits, or reduce the input."),
    ),
    code(
        "TOOL009",
        "approval_required",
        Error,
        "The build needs approval under the MCP policy (`ask`), but the MCP client cannot ask (it has no elicitation capability).",
        Some(
            "Use an MCP client with elicitation, set the rule to `allow` in the policy file, or run the build from the command line.",
        ),
    ),
    code(
        "TOOL010",
        "approval_declined",
        Error,
        "The person declined or cancelled the approval of this build.",
        None,
    ),
    code(
        "TOOL011",
        "policy_refused",
        Error,
        "The MCP policy refuses this action (`refuse`).",
        Some(
            "Change the rule in the policy file given with `--policy`, or run the build from the command line.",
        ),
    ),
    code(
        "IN001",
        "unsupported_schema",
        Error,
        "The `schema_version` of an input or artifact is not supported.",
        Some("Set `schema_version` to 1."),
    ),
    code(
        "IN002",
        "invalid_json",
        Error,
        "An input file is not valid JSON of the expected shape.",
        None,
    ),
    code(
        "MESH001",
        "invalid_domain",
        Error,
        "The domain (box, STL surface, patch map or partition) is invalid.",
        None,
    ),
    code(
        "MESH002",
        "invalid_sites",
        Error,
        "A site is outside the domain, not finite, or repeated.",
        Some("Keep every site strictly inside the domain and distinct."),
    ),
    code(
        "MESH003",
        "generation_failed",
        Error,
        "Mesh generation failed.",
        None,
    ),
    code(
        "MESH004",
        "invalid_connectivity",
        Error,
        "The mesh connectivity is inconsistent (artifact or hand-off mesh).",
        None,
    ),
    code(
        "MESH005",
        "digest_mismatch",
        Error,
        "An artifact's stored digest does not match its content: it was changed after it was built.",
        Some("Generate the artifact again."),
    ),
    code(
        "MESH006",
        "validation_failed",
        Error,
        "The mesh failed validation (closure, orientation, volume or area conservation, layer collision).",
        None,
    ),
    code(
        "MESH007",
        "invalid_size_field",
        Error,
        "The size field is invalid (sizes, gradation, frames, shapes) or asks for more seeds than allowed.",
        Some(
            "Keep min_size <= every source size <= max_size and gradation > 1; raise the sizes to need fewer seeds.",
        ),
    ),
    code(
        "IO001",
        "unsupported_format",
        Error,
        "The export format is not one of vtu, gmsh, gmsh-2.2, foehncraft, openfoam.",
        None,
    ),
    code(
        "IO002",
        "export_failed",
        Error,
        "A file could not be produced from a valid mesh.",
        None,
    ),
    code(
        "IO003",
        "io_error",
        Error,
        "A file could not be read or written (missing, exists already, permission).",
        Some("Outputs must be new; pick a name that does not exist."),
    ),
    code(
        "IO004",
        "serialization_failed",
        Error,
        "A result could not be serialized.",
        None,
    ),
];

/// The code of a snake-case name, or `MESH003` for an unknown engine code.
pub fn by_name(name: &str) -> &'static Code {
    CODES
        .iter()
        .find(|c| c.name == name)
        .or_else(|| CODES.iter().find(|c| c.code == "MESH003"))
        .unwrap_or(&CODES[0])
}

/// A diagnostic as it appears in a result envelope.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Diagnostic {
    pub code: &'static str,
    pub name: &'static str,
    pub severity: Severity,
    pub message: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fix: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
    /// A result that goes with the error (a replay report); `call` moves
    /// it into the envelope's `result`.
    #[serde(skip)]
    pub report: Option<Box<serde_json::Value>>,
}

impl Diagnostic {
    /// A diagnostic of the code with this snake-case `name`.
    pub fn new(name: &str, message: impl std::fmt::Display) -> Self {
        let c = by_name(name);
        Self {
            code: c.code,
            name: c.name,
            severity: c.severity,
            message: message.to_string(),
            fix: c.fix,
            path: None,
            report: None,
        }
    }

    pub fn with_path(mut self, path: impl std::fmt::Display) -> Self {
        self.path = Some(path.to_string());
        self
    }

    pub fn with_report(mut self, report: serde_json::Value) -> Self {
        self.report = Some(Box::new(report));
        self
    }
}

impl From<mesh_engine::EngineError> for Diagnostic {
    fn from(error: mesh_engine::EngineError) -> Self {
        let name = serde_json::to_value(error.code)
            .ok()
            .and_then(|v| v.as_str().map(str::to_string))
            .unwrap_or_else(|| "generation_failed".to_string());
        Self::new(&name, error.message)
    }
}

impl std::fmt::Display for Diagnostic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {}: {}", self.code, self.name, self.message)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Invariant: codes and names are unique, so no code has two meanings
    /// (AG-18). Adding a duplicate row breaks it.
    #[test]
    fn codes_and_names_are_unique() {
        let mut codes: Vec<_> = CODES.iter().map(|c| c.code).collect();
        let mut names: Vec<_> = CODES.iter().map(|c| c.name).collect();
        codes.sort_unstable();
        names.sort_unstable();
        codes.dedup();
        names.dedup();
        assert_eq!(codes.len(), CODES.len());
        assert_eq!(names.len(), CODES.len());
    }

    /// Invariant: every engine error code has a registered code. A new
    /// `ErrorCode` variant without a row would fall back to MESH003.
    #[test]
    fn every_engine_code_is_registered() {
        use mesh_engine::ErrorCode::*;
        for e in [
            UnsupportedSchema,
            InvalidDomain,
            InvalidSites,
            GenerationFailed,
            InvalidConnectivity,
            DigestMismatch,
            ValidationFailed,
            InvalidSizeField,
            UnsupportedFormat,
            ExportFailed,
            IoError,
            SerializationFailed,
        ] {
            let name = serde_json::to_value(e).unwrap();
            let name = name.as_str().unwrap();
            assert_eq!(by_name(name).name, name);
        }
    }
}
