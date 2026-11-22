//! Summary first, paged detail (AG-16). A paged `mesh.report` result keeps
//! the build evidence and validation, adds a summary of the per-cell and
//! per-face metrics, and carries one page of rows. The cells come first,
//! then the faces. A cursor names the row kind, the offset and the start
//! of the mesh digest, so a cursor from another or a changed artifact is
//! refused rather than returning rows of a different mesh.

use serde_json::{Map, Value, json};

use crate::diag::Diagnostic;

const KINDS: [&str; 2] = ["cells", "faces"];

fn bad_cursor(message: impl std::fmt::Display) -> Diagnostic {
    Diagnostic::new("invalid_arguments", format!("cursor: {message}"))
}

fn cursor(kind: &str, offset: usize, digest: &str) -> String {
    format!("{kind}:{offset}:{digest}")
}

fn parse(cursor: &str, digest: &str) -> Result<(usize, usize), Diagnostic> {
    let mut parts = cursor.split(':');
    let (Some(kind), Some(offset), Some(d), None) =
        (parts.next(), parts.next(), parts.next(), parts.next())
    else {
        return Err(bad_cursor("expected the `next_cursor` of a previous page"));
    };
    let kind = KINDS
        .iter()
        .position(|k| *k == kind)
        .ok_or_else(|| bad_cursor(format!("unknown row kind {kind}")))?;
    let offset = offset
        .parse()
        .map_err(|_| bad_cursor("the offset is not a number"))?;
    if d != digest {
        return Err(bad_cursor(
            "it belongs to another artifact, or the artifact changed since the cursor was issued",
        ));
    }
    Ok((kind, offset))
}

fn extreme(rows: &[Value], key: &str, max: bool) -> Value {
    let values = rows.iter().filter_map(|r| r[key].as_f64());
    let v = if max {
        values.fold(f64::NEG_INFINITY, f64::max)
    } else {
        values.fold(f64::INFINITY, f64::min)
    };
    if v.is_finite() { json!(v) } else { Value::Null }
}

/// The paged form of a full `mesh.report` result.
pub(crate) fn report_page(
    report: Value,
    at: Option<&str>,
    size: usize,
) -> Result<Value, Diagnostic> {
    let Value::Object(mut full) = report else {
        return Err(Diagnostic::new(
            "serialization_failed",
            "the report is not an object",
        ));
    };
    let digest: String = full["build"]["mesh_digest"]
        .as_str()
        .unwrap_or_default()
        .chars()
        .take(16)
        .collect();
    let mut rows: Vec<Vec<Value>> = KINDS
        .iter()
        .map(|k| match full.remove(*k) {
            Some(Value::Array(a)) => a,
            _ => Vec::new(),
        })
        .collect();
    let (mut kind, mut offset) = match at {
        Some(c) => parse(c, &digest)?,
        None => (0, 0),
    };
    if offset > rows[kind].len() {
        return Err(bad_cursor("the offset is past the last row"));
    }
    // an empty remainder of the cells moves on to the faces
    if kind == 0 && offset == rows[0].len() {
        kind = 1;
        offset = 0;
    }
    let summary = json!({
        "cells": rows[0].len(),
        "faces": rows[1].len(),
        "min_cell_volume": extreme(&rows[0], "volume", false),
        "max_cell_volume": extreme(&rows[0], "volume", true),
        "max_non_orthogonality_degrees": extreme(&rows[1], "non_orthogonality_degrees", true),
        "max_skewness": extreme(&rows[1], "skewness", true),
    });
    let end = (offset + size).min(rows[kind].len());
    let next = if end < rows[kind].len() {
        Some(cursor(KINDS[kind], end, &digest))
    } else if kind == 0 && !rows[1].is_empty() {
        Some(cursor(KINDS[1], 0, &digest))
    } else {
        None
    };
    let items: Vec<Value> = rows[kind].drain(offset..end).collect();
    let mut out = Map::new();
    for (k, v) in full {
        out.insert(k, v);
    }
    out.insert("summary".into(), summary);
    out.insert(
        "page".into(),
        json!({
            "rows": KINDS[kind], "offset": offset, "count": items.len(),
            "items": items, "next_cursor": next,
        }),
    );
    Ok(Value::Object(out))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn report(cells: usize, faces: usize) -> Value {
        json!({
            "schema_version": 1,
            "build": {"mesh_digest": "0123456789abcdef0123456789abcdef"},
            "cells": (0..cells).map(|i| json!({"volume": i as f64 + 1.0})).collect::<Vec<_>>(),
            "faces": (0..faces).map(|i| json!({"skewness": i as f64, "non_orthogonality_degrees": 0.0})).collect::<Vec<_>>(),
        })
    }

    /// Invariant (AG-16): following `next_cursor` from the first page
    /// returns every cell and then every face exactly once, in order, and
    /// ends with a null cursor. An off-by-one at a page or kind boundary
    /// breaks it; no other test walks the pages of every boundary case.
    #[test]
    fn pages_cover_every_row_once() {
        for (cells, faces, size) in [(5, 7, 2), (4, 4, 4), (3, 0, 5), (0, 3, 1), (6, 9, 3)] {
            let mut seen = [Vec::new(), Vec::new()];
            let mut at: Option<String> = None;
            for _ in 0..100 {
                let page = report_page(report(cells, faces), at.as_deref(), size).unwrap();
                assert_eq!(page["summary"]["cells"], cells);
                assert!(page.get("cells").is_none());
                let p = &page["page"];
                let k = KINDS.iter().position(|k| *k == p["rows"]).unwrap();
                assert!(p["count"].as_u64().unwrap() as usize <= size);
                seen[k].extend(p["items"].as_array().unwrap().iter().cloned());
                at = p["next_cursor"].as_str().map(str::to_string);
                if at.is_none() {
                    break;
                }
            }
            assert_eq!(
                seen[0],
                report(cells, faces)["cells"].as_array().unwrap().clone()
            );
            assert_eq!(
                seen[1],
                report(cells, faces)["faces"].as_array().unwrap().clone()
            );
        }
    }

    /// Invariant: a cursor of another artifact, or a malformed one, is
    /// refused (TOOL002) instead of paging a different mesh.
    #[test]
    fn foreign_and_malformed_cursors_are_refused() {
        for c in [
            "cells:0:ffffffffffffffff",
            "cells:x:0123456789abcdef",
            "edges:0:0123456789abcdef",
            "nope",
            "cells:99:0123456789abcdef",
        ] {
            let e = report_page(report(3, 3), Some(c), 2).unwrap_err();
            assert_eq!(e.code, "TOOL002", "{c}");
        }
    }
}
