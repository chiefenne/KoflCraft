//! Edits of a block topology file (`block.*` edit tools, BS-04). Each edit
//! takes a file and returns a new one with what changed; the new file is
//! checked before it is returned, so an edit that breaks the topology is
//! refused and nothing is written. Blocks and edges are named by their
//! vertex indices, as in the file.

use std::collections::BTreeSet;

use mesh_core::block::{block_edges, corner_slot};
use mesh_io::handoff::PatchRole;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::block::{CurveFile, DistributionFile, EdgeFile, TopologyFile, check};
use crate::{EngineError, ErrorCode};

/// One edit of a topology.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(tag = "op", rename_all = "snake_case", deny_unknown_fields)]
pub enum Edit {
    /// Set the cell count of the chain of opposite edges through `edge`
    /// (BS-05): every block edge in the chain gets `cells`.
    SetCount { edge: [usize; 2], cells: usize },
    /// Set the law of one edge, read from `edge[0]` to `edge[1]`.
    SetDistribution {
        edge: [usize; 2],
        distribution: DistributionFile,
    },
    /// Set the curve of one edge; spline points run from `edge[0]` to
    /// `edge[1]`.
    SetCurve { edge: [usize; 2], curve: CurveFile },
    /// Move a vertex to `to`, or by `by`.
    MoveVertex {
        vertex: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        to: Option<[f64; 3]>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        by: Option<[f64; 3]>,
    },
    /// Rename a patch and/or set its role.
    SetPatch {
        name: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        rename: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        role: Option<PatchRole>,
    },
}

/// What an edit changed.
#[derive(Clone, Debug, Default, Serialize, PartialEq)]
pub struct Change {
    pub blocks: Vec<usize>,
    pub edges: Vec<[usize; 2]>,
    pub vertices: Vec<usize>,
    pub patches: Vec<String>,
    /// The template's analytic volume no longer describes the shape after a
    /// geometric edit, so the file's `reference` is removed.
    pub reference_removed: bool,
}

fn refused(message: impl std::fmt::Display) -> EngineError {
    EngineError::new(ErrorCode::InvalidArguments, message)
}

fn key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

/// The edges of block `b` along each axis, as vertex keys.
fn axis_edges(file: &TopologyFile, b: usize) -> [Vec<(usize, usize)>; 3] {
    let v = &file.blocks[b].vertices;
    let mut out: [Vec<(usize, usize)>; 3] = Default::default();
    for (axis, c) in block_edges() {
        let mut e = c;
        e[axis] = 1;
        let (p, q) = (
            v[corner_slot(c[0], c[1], c[2])],
            v[corner_slot(e[0], e[1], e[2])],
        );
        out[axis].push(key(p, q));
    }
    out
}

/// The blocks and axes that carry `edge`; refused when it is no block edge.
fn carriers(file: &TopologyFile, edge: [usize; 2]) -> Result<Vec<(usize, usize)>, EngineError> {
    let k = key(edge[0], edge[1]);
    let found: Vec<(usize, usize)> = (0..file.blocks.len())
        .flat_map(|b| {
            let edges = axis_edges(file, b);
            (0..3)
                .filter(move |&a| edges[a].contains(&k))
                .map(move |a| (b, a))
                .collect::<Vec<_>>()
        })
        .collect();
    if found.is_empty() {
        return Err(refused(format!(
            "edge {}–{} is not an edge of any block",
            edge[0], edge[1]
        )));
    }
    Ok(found)
}

/// The same law read from the other end.
fn reversed_law(d: &DistributionFile) -> DistributionFile {
    match d {
        DistributionFile::Uniform => DistributionFile::Uniform,
        DistributionFile::Geometric { first, last } => DistributionFile::Geometric {
            first: *last,
            last: *first,
        },
        DistributionFile::Tanh { first, last } => DistributionFile::Tanh {
            first: *last,
            last: *first,
        },
    }
}

fn reversed_curve(c: &CurveFile) -> CurveFile {
    match c {
        CurveFile::Spline { through } => CurveFile::Spline {
            through: through.iter().rev().copied().collect(),
        },
        other => other.clone(),
    }
}

/// The file's entry for `edge`, created as a uniform line when absent,
/// and whether it runs the other way.
fn edge_entry(file: &mut TopologyFile, edge: [usize; 2]) -> (&mut EdgeFile, bool) {
    let k = key(edge[0], edge[1]);
    let i = match file
        .edges
        .iter()
        .position(|e| key(e.vertices[0], e.vertices[1]) == k)
    {
        Some(i) => i,
        None => {
            file.edges.push(EdgeFile {
                vertices: edge,
                curve: CurveFile::Line,
                distribution: DistributionFile::Uniform,
            });
            file.edges.len() - 1
        }
    };
    let e = &mut file.edges[i];
    let flipped = e.vertices != edge;
    (e, flipped)
}

/// Apply one edit; the result is checked (BS-05) or refused.
pub fn apply(file: &TopologyFile, edit: &Edit) -> Result<(TopologyFile, Change), EngineError> {
    // the input must be a valid topology, so a refusal names this edit
    check(file)?;
    let mut out = file.clone();
    let mut change = Change::default();
    match edit {
        Edit::SetCount { edge, cells } => {
            if *cells == 0 {
                return Err(refused("a cell count must be at least 1"));
            }
            let mut marked: BTreeSet<(usize, usize)> = carriers(file, *edge)?.into_iter().collect();
            let mut keys: BTreeSet<(usize, usize)> = BTreeSet::new();
            // grow the chain: a block axis joins when it shares an edge
            loop {
                for &(b, a) in &marked {
                    keys.extend(axis_edges(file, b)[a].iter().copied());
                }
                let before = marked.len();
                for b in 0..file.blocks.len() {
                    for (a, edges) in axis_edges(file, b).iter().enumerate() {
                        if edges.iter().any(|k| keys.contains(k)) {
                            marked.insert((b, a));
                        }
                    }
                }
                if marked.len() == before {
                    break;
                }
            }
            for &(b, a) in &marked {
                if out.blocks[b].cells[a] != *cells {
                    out.blocks[b].cells[a] = *cells;
                    change.blocks.push(b);
                }
            }
            change.blocks.dedup();
            change.edges = keys.into_iter().map(|(p, q)| [p, q]).collect();
        }
        Edit::SetDistribution { edge, distribution } => {
            carriers(file, *edge)?;
            let (e, flipped) = edge_entry(&mut out, *edge);
            e.distribution = if flipped {
                reversed_law(distribution)
            } else {
                distribution.clone()
            };
            change.edges.push(*edge);
        }
        Edit::SetCurve { edge, curve } => {
            carriers(file, *edge)?;
            let (e, flipped) = edge_entry(&mut out, *edge);
            e.curve = if flipped {
                reversed_curve(curve)
            } else {
                curve.clone()
            };
            change.edges.push(*edge);
            change.reference_removed = out.reference.take().is_some();
        }
        Edit::MoveVertex { vertex, to, by } => {
            let Some(p) = out.vertices.get_mut(*vertex) else {
                return Err(refused(format!(
                    "vertex {vertex} does not exist (the file has {})",
                    file.vertices.len()
                )));
            };
            *p = match (to, by) {
                (Some(t), None) => *t,
                (None, Some(d)) => [p[0] + d[0], p[1] + d[1], p[2] + d[2]],
                _ => return Err(refused("give exactly one of `to` and `by`")),
            };
            if p.iter().any(|c| !c.is_finite()) {
                return Err(refused("the new position is not finite"));
            }
            change.vertices.push(*vertex);
            change.reference_removed = out.reference.take().is_some();
        }
        Edit::SetPatch { name, rename, role } => {
            if rename.is_none() && role.is_none() {
                return Err(refused("give `rename`, `role` or both"));
            }
            let Some(i) = out.patches.iter().position(|p| &p.name == name) else {
                let known: Vec<&str> = file.patches.iter().map(|p| p.name.as_str()).collect();
                return Err(EngineError::new(
                    ErrorCode::UnknownPatch,
                    format!("no patch `{name}` (patches: {})", known.join(", ")),
                ));
            };
            if let Some(new) = rename {
                if new != name && out.patches.iter().any(|p| &p.name == new) {
                    return Err(refused(format!("a patch `{new}` exists already")));
                }
                out.patches[i].name = new.clone();
            }
            if let Some(r) = role {
                out.patches[i].role = Some(*r);
            }
            change.patches.push(out.patches[i].name.clone());
        }
    }
    check(&out)?;
    Ok((out, change))
}

/// The contents of a topology for `block.list`: vertices, blocks, edges
/// with their curves and laws, patches with roles, surfaces, the check.
pub fn list(file: &TopologyFile) -> Result<Value, EngineError> {
    let checked = check(file)?;
    let curve = |c: &CurveFile| match c {
        CurveFile::Line => "line",
        CurveFile::Arc { .. } => "arc",
        CurveFile::Spline { .. } => "spline",
    };
    Ok(json!({
        "vertices": file.vertices.iter().enumerate()
            .map(|(i, p)| json!({"index": i, "point": p})).collect::<Vec<_>>(),
        "blocks": file.blocks.iter().enumerate()
            .map(|(i, b)| json!({"index": i, "vertices": b.vertices, "cells": b.cells, "domain": b.domain}))
            .collect::<Vec<_>>(),
        "edges": file.edges.iter()
            .map(|e| json!({"vertices": e.vertices, "curve": curve(&e.curve), "distribution": e.distribution}))
            .collect::<Vec<_>>(),
        "patches": file.patches.iter()
            .map(|p| json!({"name": p.name, "role": p.role, "faces": p.faces.len()}))
            .collect::<Vec<_>>(),
        "surfaces": file.surfaces.iter()
            .map(|s| json!({"name": s.name, "faces": s.faces.len()})).collect::<Vec<_>>(),
        "reference": file.reference.as_ref().map(|r| json!({"template": r.template, "analytic_volume": r.analytic_volume})),
        "check": checked,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::block::{TemplateConfig, template};

    fn pipe() -> TopologyFile {
        template(&TemplateConfig::PipeOgrid {
            radius: 0.05,
            length: 0.2,
            core_fraction: None,
            core_cells: 4,
            radial_cells: 3,
            axial_cells: 6,
            wall_spacing: None,
        })
        .unwrap()
    }

    /// Invariant (BS-04, BS-05): a count set on one edge reaches every
    /// block of its chain, so the result checks, and the cell total follows.
    #[test]
    fn a_count_follows_its_chain() {
        let f = pipe();
        // the axial edge of the core block runs from vertex slot 0 to 4
        let v = f.blocks[0].vertices;
        let (g, c) = apply(
            &f,
            &Edit::SetCount {
                edge: [v[0], v[4]],
                cells: 9,
            },
        )
        .unwrap();
        assert_eq!(c.blocks, [0, 1, 2, 3, 4]);
        assert!(g.blocks.iter().all(|b| b.cells[2] == 9));
        let cells = check(&g).unwrap()["cells"].as_u64().unwrap();
        assert_eq!(cells, 9 * (16 + 4 * 12));
        assert_eq!(g.reference, f.reference, "a count keeps the shape");
    }

    /// Invariant: a law given from the other end of an edge is stored
    /// reversed, so it means the same spacing at the same vertex; a law
    /// on a vertex pair that is no block edge is refused.
    #[test]
    fn a_law_is_stored_in_the_edge_direction() {
        let f = pipe();
        let v = f.blocks[0].vertices;
        let (g, _) = apply(
            &f,
            &Edit::SetDistribution {
                edge: [v[4], v[0]],
                distribution: DistributionFile::Geometric {
                    first: Some(0.01),
                    last: None,
                },
            },
        )
        .unwrap();
        let e = g
            .edges
            .iter()
            .find(|e| key(e.vertices[0], e.vertices[1]) == key(v[0], v[4]))
            .unwrap();
        let expected = if e.vertices == [v[4], v[0]] {
            DistributionFile::Geometric {
                first: Some(0.01),
                last: None,
            }
        } else {
            DistributionFile::Geometric {
                first: None,
                last: Some(0.01),
            }
        };
        assert_eq!(e.distribution, expected);
        let err = apply(
            &f,
            &Edit::SetDistribution {
                edge: [v[0], v[6]],
                distribution: DistributionFile::Uniform,
            },
        )
        .unwrap_err();
        assert!(err.to_string().contains("not an edge"), "{err}");
    }

    /// Invariant: moving a vertex that a projected face needs off its
    /// surface is refused by the check, nothing else changes; a move that
    /// keeps the topology valid drops the template's analytic reference.
    #[test]
    fn moves_are_checked_and_drop_the_reference() {
        let f = template(&TemplateConfig::Elbow {
            radius: 0.5,
            bend_radius: 1.5,
            angle: 90.0,
            inlet_length: 1.0,
            outlet_length: 1.0,
            core_fraction: None,
            core_cells: 4,
            radial_cells: 3,
            axial_cells: 12,
            wall_spacing: None,
        })
        .unwrap();
        // a wall vertex is a corner of a projected face
        let wall = f.surfaces[0].faces[0][0];
        let err = apply(
            &f,
            &Edit::MoveVertex {
                vertex: wall,
                to: None,
                by: Some([0.0, 0.0, 0.01]),
            },
        )
        .unwrap_err();
        assert!(err.to_string().contains("surface"), "{err}");
        // a corner of the pipe's core block is on no surface: a small move
        // inside the core keeps the topology valid
        let p = pipe();
        let core = p.blocks[0].vertices[4];
        let (g, c) = apply(
            &p,
            &Edit::MoveVertex {
                vertex: core,
                to: None,
                by: Some([0.001, 0.0, 0.0]),
            },
        )
        .unwrap();
        assert_eq!(g.vertices[core][0], p.vertices[core][0] + 0.001);
        assert!(c.reference_removed && g.reference.is_none());
        assert!(
            apply(
                &p,
                &Edit::MoveVertex {
                    vertex: 10_000,
                    to: Some([0.0; 3]),
                    by: None
                }
            )
            .is_err()
        );
    }

    /// Invariant: a patch keeps its faces when renamed or given a role; a
    /// taken name and an unknown patch are refused.
    #[test]
    fn patches_are_renamed_and_given_roles() {
        let f = pipe();
        let (g, c) = apply(
            &f,
            &Edit::SetPatch {
                name: "walls".into(),
                rename: Some("pipe_wall".into()),
                role: Some(PatchRole::Wall),
            },
        )
        .unwrap();
        assert_eq!(c.patches, ["pipe_wall"]);
        let p = g.patches.iter().find(|p| p.name == "pipe_wall").unwrap();
        assert_eq!(p.role, Some(PatchRole::Wall));
        assert_eq!(
            p.faces,
            f.patches.iter().find(|p| p.name == "walls").unwrap().faces
        );
        let taken = Edit::SetPatch {
            name: "walls".into(),
            rename: Some("inlet".into()),
            role: None,
        };
        assert!(apply(&f, &taken).is_err());
        let unknown = Edit::SetPatch {
            name: "nope".into(),
            rename: None,
            role: Some(PatchRole::Wall),
        };
        assert!(apply(&f, &unknown).is_err());
    }
}
