//! Layout optimisation of an extruded block topology (BS-03, MB1g).
//!
//! The vertices and edge shapes of the blocks away from the walls move so
//! that the cells of the transfinite grid come closer to orthogonal in the
//! finite-volume sense: the angle between a face's normal and the line
//! through the centroids of its two cells (non-orthogonality). The idea is
//! the block-layout optimisation of Sauer and Morsbach ("An optimization
//! based multi-block-structured grid generation method", IJNME 2023,
//! doi:10.1002/nme.7308; `research/block-layout-optimisation.md`), which
//! moves the block vertices and the control points of the block edges so
//! that the inner cell angles approach 90°. Here the criterion is measured on the
//! grid itself, so a layout is judged by the cells a solver will see, and
//! it is minimised by a pattern search (Hooke and Jeeves, "Direct search
//! solution of numerical and statistical problems", J. ACM 8, 1961): the
//! worst angle is not smooth in the parameters, and one trial move only
//! re-grids the blocks next to it.
//!
//! On a topology extruded along z (two layers of vertices):
//! - a block with a face on a rigid patch (a wall) keeps its vertices and
//!   edges: that is the near-wall layer the template built, such as an
//!   O-ring;
//! - a vertex on the outer boundary slides along it where the boundary is
//!   straight; every other free vertex moves in the plane;
//! - an edge inside the domain bends by two Hermite modes on top of its
//!   template shape, which set the angles at which it leaves its ends.
//!
//! Every accepted move keeps every cell unfolded; folded cells of the
//! start layout carry a large cost, so the search also unfolds them.

use std::collections::{BTreeSet, HashMap};

use mesh_core::block::{BlockTopology, Curve, Distribution, EdgeEnd, EdgeSpec};
use mesh_core::geometry::Vec3;

use super::distribution::{Evaluator, fractions};
use super::{DEFAULT_EDGE, edge_key};

/// Angle (degrees) of the quadratic part of the face cost.
const MEAN_ANGLE: f64 = 45.0;
/// Angle (degrees) of the steep part of the face cost, which makes the
/// worst faces dominate.
const WORST_ANGLE: f64 = 40.0;
const WORST_POWER: i32 = 12;
/// Growth of the cells across a face that costs as much as a 45° face.
const GROWTH_SCALE: f64 = 2.5;
/// Cost of one folded cell corner.
const FOLD_COST: f64 = 1.0e4;
/// Corners sharper than this sine cost `SHARP_COST` times the squared
/// shortfall over it, so the search sees the way out of a fold.
const SHARP_SINE: f64 = 0.1;
const SHARP_COST: f64 = 100.0;
/// Farthest an end spacing may move from the template's, as a factor.
const MAX_LAW_FACTOR: f64 = 20.0;
/// First step of a vertex move, over its shortest edge, and of a bend.
const FIRST_STEP: f64 = 0.02;
/// The search on one parameter ends when its step has shrunk this much.
const STEP_RANGE: f64 = 100.0;
/// Largest bend, in chords of the edge.
const MAX_BEND: f64 = 0.5;
/// Shortest a boundary edge may get, over its start length.
const MIN_SIDE: f64 = 0.1;
/// Points of the reference polyline of a bent edge.
const SAMPLES: usize = 9;

/// What [`optimize`] may change and how long it searches.
#[derive(Clone, Debug)]
pub struct LayoutOptions {
    /// Patches whose blocks keep their vertices and edges.
    pub rigid_patches: Vec<usize>,
    /// Most sweeps over all parameters.
    pub max_sweeps: usize,
}

impl Default for LayoutOptions {
    fn default() -> Self {
        Self {
            rigid_patches: Vec::new(),
            max_sweeps: 40,
        }
    }
}

/// Non-orthogonality of the grid in the plane of the extrusion.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Orthogonality {
    /// Largest and mean angle (degrees) between a face's normal and the
    /// line through the centroids of its two cells, over inner faces.
    pub max_deg: f64,
    pub mean_deg: f64,
    /// Cells with a face above 70°.
    pub cells_over_70: usize,
    /// Largest ratio of the centroid distances from an inner face, ≥ 1.
    pub max_growth: f64,
    /// Cell corners where a cell is folded.
    pub folded_corners: usize,
}

/// What [`optimize`] did.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LayoutReport {
    /// Parameters searched: one per sliding vertex, two per free vertex
    /// and per bendable edge.
    pub parameters: usize,
    pub sweeps: usize,
    pub trials: usize,
    pub before: Orthogonality,
    pub after: Orthogonality,
}

type P = [f64; 2];

fn sub(a: P, b: P) -> P {
    [a[0] - b[0], a[1] - b[1]]
}

fn cross(a: P, b: P) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn norm(a: P) -> f64 {
    a[0].hypot(a[1])
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    /// A side of a rigid block.
    Rigid,
    /// On the outer boundary.
    Boundary,
    /// Inside the domain: bends.
    Flexible,
}

struct Edge {
    /// Topology vertex ids (layer 0) in the stored direction.
    a: usize,
    b: usize,
    curve: Curve,
    law: Distribution,
    /// The template's law.
    law0: Distribution,
    cells: usize,
    kind: Kind,
    /// Flexible edges: the template shape as a polyline with its ends and
    /// the chord fractions of its points, the unit normal and length of
    /// its chord.
    reference: Vec<P>,
    s: Vec<f64>,
    normal: P,
    chord: f64,
    straight: bool,
    bend: [f64; 2],
    nodes: Vec<P>,
    /// Unit tangents at a (into the edge) and at b (back into the edge).
    tangents: [P; 2],
}

struct Block {
    v: [usize; 4],
    n: [usize; 2],
    /// bottom v0→v1, right v1→v2, top v3→v2, left v0→v3: the edge and
    /// whether it runs that way.
    sides: [(usize, bool); 4],
    rigid: bool,
    sign: f64,
    nodes: Vec<P>,
    centroids: Vec<P>,
    folded: usize,
    /// Faces inside the block and folded corners.
    cost: f64,
}

/// The cells on both sides of one block edge shared by two blocks.
struct Interface {
    edge: usize,
    sides: [(usize, usize); 2],
    cost: f64,
}

#[derive(Clone, Copy, Debug)]
enum Param {
    /// Vertex `v` moves along the unit vector `dir`; `slide` holds the
    /// boundary edges it slides between.
    Move { v: usize, dir: P },
    /// Mode `mode` of edge `e`.
    Bend { e: usize, mode: usize },
    /// The logarithm of the first (`end` 0) or last spacing of edge `e`.
    Law { e: usize, end: usize },
}

struct Dof {
    param: Param,
    step: f64,
    min_step: f64,
    max_step: f64,
    sign: f64,
    edges: Vec<usize>,
    blocks: Vec<usize>,
    interfaces: Vec<usize>,
    /// Boundary edges a sliding vertex must keep long enough.
    sides: Vec<usize>,
}

struct Model {
    z: f64,
    pos: Vec<P>,
    reference: Vec<P>,
    edges: Vec<Edge>,
    blocks: Vec<Block>,
    interfaces: Vec<Interface>,
}

/// Cost of the face `p`–`q` between the cells with centroids `c0`, `c1`:
/// its non-orthogonality and the growth of the cells across it.
fn face_cost(p: P, q: P, c0: P, c1: P) -> f64 {
    let theta = face_angle(p, q, c0, c1);
    let a = theta / MEAN_ANGLE.to_radians();
    let w = theta / WORST_ANGLE.to_radians();
    let g = growth(p, q, c0, c1).ln() / GROWTH_SCALE.ln();
    a * a + w.powi(WORST_POWER) + g * g
}

/// Ratio of the distances of the two centroids from the face's middle.
fn growth(p: P, q: P, c0: P, c1: P) -> f64 {
    let m = [0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1])];
    let (d0, d1) = (norm(sub(c0, m)), norm(sub(c1, m)));
    let r = d0.max(d1) / d0.min(d1);
    if r.is_finite() { r } else { 1.0e6 }
}

/// `law` with the spacing at `end` (0 first, 1 last, as the edge runs)
/// times `factor`.
fn scale_law(law: Distribution, end: usize, factor: f64) -> Distribution {
    match law {
        Distribution::Tanh { first, last } => {
            if end == 0 {
                Distribution::Tanh {
                    first: first * factor,
                    last,
                }
            } else {
                Distribution::Tanh {
                    first,
                    last: last * factor,
                }
            }
        }
        Distribution::Geometric { spacing, at } => Distribution::Geometric {
            spacing: spacing * factor,
            at,
        },
        Distribution::Uniform => Distribution::Uniform,
    }
}

/// The spacing at `end` of a law, if it has one.
fn law_spacing(law: Distribution, end: usize) -> Option<f64> {
    match (law, end) {
        (Distribution::Tanh { first, .. }, 0) => Some(first),
        (Distribution::Tanh { last, .. }, _) => Some(last),
        (Distribution::Geometric { spacing, .. }, 0) => Some(spacing),
        _ => None,
    }
}

fn reverse_law(law: Distribution) -> Distribution {
    match law {
        Distribution::Tanh { first, last } => Distribution::Tanh {
            first: last,
            last: first,
        },
        Distribution::Geometric { spacing, at } => Distribution::Geometric {
            spacing,
            at: match at {
                EdgeEnd::Start => EdgeEnd::End,
                EdgeEnd::End => EdgeEnd::Start,
            },
        },
        Distribution::Uniform => Distribution::Uniform,
    }
}

/// Angle between the normal of the face `p`–`q` and the line `c0`–`c1`.
fn face_angle(p: P, q: P, c0: P, c1: P) -> f64 {
    let t = sub(q, p);
    let d = sub(c1, c0);
    let (lt, ld) = (norm(t), norm(d));
    if lt <= 0.0 || ld <= 0.0 {
        return std::f64::consts::FRAC_PI_2;
    }
    // |n · d| with n = (−t_y, t_x)
    let c = (t[0] * d[1] - t[1] * d[0]).abs() / (lt * ld);
    c.min(1.0).acos()
}

fn chord_fractions(p: &[P]) -> Vec<f64> {
    let mut s = vec![0.0; p.len()];
    for k in 1..p.len() {
        s[k] = s[k - 1] + norm(sub(p[k], p[k - 1]));
    }
    let total = s[p.len() - 1];
    if total > 0.0 {
        s.iter_mut().for_each(|v| *v /= total);
    }
    s
}

impl Block {
    fn cell(&self, side: usize, k: usize) -> usize {
        let [ni, nj] = self.n;
        match side {
            0 => k,
            1 => ni - 1 + ni * k,
            2 => k + ni * (nj - 1),
            _ => ni * k,
        }
    }

    /// Grid the block from its sides by 2D transfinite interpolation with
    /// the blended parameters of `grid::build`, then its centroids and
    /// cost.
    fn fill(&mut self, edges: &[Edge]) {
        let [ni, nj] = self.n;
        let side = |s: usize| -> Vec<P> {
            let (e, fwd) = self.sides[s];
            let mut v = edges[e].nodes.clone();
            if !fwd {
                v.reverse();
            }
            v
        };
        let (bottom, right, top, left) = (side(0), side(1), side(2), side(3));
        let (sb, st, sl, sr) = (
            chord_fractions(&bottom),
            chord_fractions(&top),
            chord_fractions(&left),
            chord_fractions(&right),
        );
        let w = ni + 1;
        let mut x = vec![[0.0; 2]; w * (nj + 1)];
        for i in 0..=ni {
            x[i] = bottom[i];
            x[i + w * nj] = top[i];
        }
        for j in 0..=nj {
            x[w * j] = left[j];
            x[ni + w * j] = right[j];
        }
        let (p00, p10, p11, p01) = (bottom[0], bottom[ni], top[ni], top[0]);
        for j in 1..nj {
            for i in 1..ni {
                let (db, dl) = (st[i] - sb[i], sr[j] - sl[j]);
                let u = (sb[i] + sl[j] * db) / (1.0 - dl * db);
                let v = sl[j] + u * dl;
                let mut p = [0.0; 2];
                for d in 0..2 {
                    p[d] = bottom[i][d] * (1.0 - v)
                        + top[i][d] * v
                        + left[j][d] * (1.0 - u)
                        + right[j][d] * u
                        - (p00[d] * ((1.0 - u) * (1.0 - v))
                            + p10[d] * (u * (1.0 - v))
                            + p11[d] * (u * v)
                            + p01[d] * ((1.0 - u) * v));
                }
                x[i + w * j] = p;
            }
        }
        let mut centroids = Vec::with_capacity(ni * nj);
        let mut folded = 0usize;
        let mut sharp = 0.0;
        for j in 0..nj {
            for i in 0..ni {
                let q = [
                    x[i + w * j],
                    x[i + 1 + w * j],
                    x[i + 1 + w * (j + 1)],
                    x[i + w * (j + 1)],
                ];
                let (mut a, mut cx, mut cy) = (0.0, 0.0, 0.0);
                for k in 0..4 {
                    let (p, r) = (q[k], q[(k + 1) % 4]);
                    let c = p[0] * r[1] - r[0] * p[1];
                    a += c;
                    cx += (p[0] + r[0]) * c;
                    cy += (p[1] + r[1]) * c;
                    let (u, w) = (sub(r, p), sub(q[(k + 3) % 4], p));
                    let sine = self.sign * cross(u, w) / (norm(u) * norm(w));
                    if sine.is_nan() || sine <= 0.0 {
                        folded += 1;
                    }
                    if sine.is_nan() || sine < SHARP_SINE {
                        let d = (SHARP_SINE - sine.max(-1.0)) / SHARP_SINE;
                        sharp += SHARP_COST * d * d;
                    }
                }
                centroids.push(if a.abs() > 0.0 {
                    [cx / (3.0 * a), cy / (3.0 * a)]
                } else {
                    [
                        0.25 * (q[0][0] + q[1][0] + q[2][0] + q[3][0]),
                        0.25 * (q[0][1] + q[1][1] + q[2][1] + q[3][1]),
                    ]
                });
            }
        }
        let mut cost = FOLD_COST * folded as f64 + sharp;
        for j in 0..nj {
            for i in 1..ni {
                cost += face_cost(
                    x[i + w * j],
                    x[i + w * (j + 1)],
                    centroids[i - 1 + ni * j],
                    centroids[i + ni * j],
                );
            }
        }
        for j in 1..nj {
            for i in 0..ni {
                cost += face_cost(
                    x[i + w * j],
                    x[i + 1 + w * j],
                    centroids[i + ni * (j - 1)],
                    centroids[i + ni * j],
                );
            }
        }
        self.nodes = x;
        self.centroids = centroids;
        self.folded = folded;
        self.cost = cost;
    }

    /// Angles (radians) of the faces inside the block, and the folded
    /// corners.
    fn angles(&self) -> (Vec<(usize, usize, f64, f64)>, usize) {
        let [ni, nj] = self.n;
        let w = ni + 1;
        let x = &self.nodes;
        let c = &self.centroids;
        let mut out = Vec::new();
        for j in 0..nj {
            for i in 1..ni {
                let (a, b) = (i - 1 + ni * j, i + ni * j);
                let (p, q) = (x[i + w * j], x[i + w * (j + 1)]);
                out.push((a, b, face_angle(p, q, c[a], c[b]), growth(p, q, c[a], c[b])));
            }
        }
        for j in 1..nj {
            for i in 0..ni {
                let (a, b) = (i + ni * (j - 1), i + ni * j);
                let (p, q) = (x[i + w * j], x[i + 1 + w * j]);
                out.push((a, b, face_angle(p, q, c[a], c[b]), growth(p, q, c[a], c[b])));
            }
        }
        (out, self.folded)
    }
}

impl Model {
    /// The 2D model of an extruded topology: blocks whose first four
    /// vertices lie on one plane z = const and whose last four repeat them
    /// on another.
    fn new(t: &BlockTopology, rigid_patches: &[usize]) -> Result<Self, String> {
        super::check(t).map_err(|issues| {
            issues
                .iter()
                .map(|i| i.message.clone())
                .collect::<Vec<_>>()
                .join("; ")
        })?;
        let z = t.vertices[t.blocks[0].vertices[0]].z;
        let mut layer0 = BTreeSet::new();
        for (b, block) in t.blocks.iter().enumerate() {
            for k in 0..4 {
                let (p, q) = (
                    t.vertices[block.vertices[k]],
                    t.vertices[block.vertices[k + 4]],
                );
                if p.z != z || q.z == z || p.x != q.x || p.y != q.y {
                    return Err(format!(
                        "block {b} is not extruded along z from the plane z = {z}"
                    ));
                }
                layer0.insert(block.vertices[k]);
            }
        }
        let pos: Vec<P> = t.vertices.iter().map(|p| [p.x, p.y]).collect();
        // the edges of the plane, with their specs
        let mut spec: HashMap<(usize, usize), usize> = HashMap::new();
        for (i, e) in t.edges.iter().enumerate() {
            spec.insert(edge_key(e.vertices[0], e.vertices[1]), i);
        }
        let mut edges: Vec<Edge> = Vec::new();
        let mut edge_of: HashMap<(usize, usize), usize> = HashMap::new();
        let mut blocks = Vec::with_capacity(t.blocks.len());
        for block in &t.blocks {
            let v = [
                block.vertices[0],
                block.vertices[1],
                block.vertices[2],
                block.vertices[3],
            ];
            let n = [block.cells[0], block.cells[1]];
            let ends = [
                (v[0], v[1], n[0]),
                (v[1], v[2], n[1]),
                (v[3], v[2], n[0]),
                (v[0], v[3], n[1]),
            ];
            let sides = ends.map(|(a, b, cells)| {
                let key = edge_key(a, b);
                let e = *edge_of.entry(key).or_insert_with(|| {
                    let (from, to, curve, law) = match spec.get(&key) {
                        Some(&i) => (
                            t.edges[i].vertices[0],
                            t.edges[i].vertices[1],
                            t.edges[i].curve.clone(),
                            t.edges[i].distribution,
                        ),
                        None => (key.0, key.1, DEFAULT_EDGE.0, DEFAULT_EDGE.1),
                    };
                    edges.push(Edge {
                        a: from,
                        b: to,
                        curve,
                        law,
                        law0: law,
                        cells,
                        kind: Kind::Flexible,
                        reference: Vec::new(),
                        s: Vec::new(),
                        normal: [0.0; 2],
                        chord: 0.0,
                        straight: false,
                        bend: [0.0; 2],
                        nodes: Vec::new(),
                        tangents: [[0.0; 2]; 2],
                    });
                    edges.len() - 1
                });
                (e, edges[e].a == a)
            });
            blocks.push(Block {
                v,
                n,
                sides,
                rigid: false,
                sign: 1.0,
                nodes: Vec::new(),
                centroids: Vec::new(),
                folded: 0,
                cost: 0.0,
            });
        }
        // rigid blocks: a side face on a rigid patch
        let mut rigid_keys = BTreeSet::new();
        for &p in rigid_patches {
            let Some(patch) = t.patches.get(p) else {
                return Err(format!("rigid patch {p} does not exist"));
            };
            for f in &patch.faces {
                let low: Vec<usize> = f.iter().copied().filter(|v| layer0.contains(v)).collect();
                if low.len() == 2 {
                    rigid_keys.insert(edge_key(low[0], low[1]));
                }
            }
        }
        let mut uses = vec![0usize; edges.len()];
        for b in blocks.iter_mut() {
            b.rigid = b
                .sides
                .iter()
                .any(|&(e, _)| rigid_keys.contains(&edge_key(edges[e].a, edges[e].b)));
            for &(e, _) in &b.sides {
                uses[e] += 1;
            }
        }
        for (e, edge) in edges.iter_mut().enumerate() {
            if uses[e] == 1 {
                edge.kind = Kind::Boundary;
            }
        }
        for b in &blocks {
            if b.rigid {
                for &(e, _) in &b.sides {
                    edges[e].kind = Kind::Rigid;
                }
            }
        }
        let mut model = Self {
            z,
            reference: pos.clone(),
            pos,
            edges,
            blocks,
            interfaces: Vec::new(),
        };
        for e in 0..model.edges.len() {
            if model.edges[e].kind == Kind::Flexible {
                model.set_reference(e)?;
            }
            (model.edges[e].nodes, model.edges[e].tangents) = model.edge_nodes(e)?;
        }
        // interfaces between blocks
        let mut owners: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
        for (b, block) in model.blocks.iter().enumerate() {
            for (s, &(e, _)) in block.sides.iter().enumerate() {
                owners.entry(e).or_default().push((b, s));
            }
        }
        let mut shared: Vec<(usize, Vec<(usize, usize)>)> = owners.into_iter().collect();
        shared.sort();
        for (e, o) in shared {
            if o.len() == 2 {
                model.interfaces.push(Interface {
                    edge: e,
                    sides: [o[0], o[1]],
                    cost: 0.0,
                });
            }
        }
        for b in 0..model.blocks.len() {
            let (edges, block) = (&model.edges, &mut model.blocks[b]);
            block.fill(edges);
            // the orientation the most cells have
            let w = block.n[0] + 1;
            let mut area = 0.0;
            for j in 0..block.n[1] {
                for i in 0..block.n[0] {
                    let x = &block.nodes;
                    let (p, q, r) = (x[i + w * j], x[i + 1 + w * j], x[i + w * (j + 1)]);
                    area += cross(sub(q, p), sub(r, p)).signum();
                }
            }
            block.sign = if area < 0.0 { -1.0 } else { 1.0 };
            block.fill(edges);
        }
        for i in 0..model.interfaces.len() {
            model.interfaces[i].cost = model.interface_cost(i);
        }
        Ok(model)
    }

    /// The template shape of flexible edge `e` as a polyline.
    fn set_reference(&mut self, e: usize) -> Result<(), String> {
        let z = self.z;
        let edge = &self.edges[e];
        let (a, b) = (self.pos[edge.a], self.pos[edge.b]);
        let v3 = |p: P| Vec3::new(p[0], p[1], z);
        let straight = edge.curve == Curve::Line;
        let reference: Vec<P> = match &edge.curve {
            Curve::Spline { through } if through.len() + 2 >= SAMPLES => {
                let mut r = vec![a];
                r.extend(through.iter().map(|q| [q.x, q.y]));
                r.push(b);
                r
            }
            curve => {
                let ev = Evaluator::new(v3(a), v3(b), curve)?;
                (0..SAMPLES)
                    .map(|k| {
                        if k == 0 {
                            a
                        } else if k == SAMPLES - 1 {
                            b
                        } else {
                            let q = ev.at(k as f64 / (SAMPLES - 1) as f64);
                            [q.x, q.y]
                        }
                    })
                    .collect()
            }
        };
        let c = sub(b, a);
        let chord = norm(c);
        let edge = &mut self.edges[e];
        edge.s = chord_fractions(&reference);
        edge.reference = reference;
        edge.normal = [-c[1] / chord, c[0] / chord];
        edge.chord = chord;
        edge.straight = straight;
        Ok(())
    }

    /// The corners of block `b` are right-handed with the edges' end
    /// tangents, as the topology check needs them.
    fn corners_turn(&self, b: usize) -> bool {
        let k = &self.blocks[b];
        // the tangent of side `s` leaving its start, or leaving its end
        // back along it
        let t = |s: usize, at_start: bool| {
            let (e, fwd) = k.sides[s];
            self.edges[e].tangents[usize::from(at_start != fwd)]
        };
        let neg = |p: P| [-p[0], -p[1]];
        let corners = [
            (t(0, true), t(3, true)),
            (neg(t(0, false)), t(1, true)),
            (neg(t(2, false)), neg(t(1, false))),
            (t(2, true), neg(t(3, false))),
        ];
        corners.iter().all(|&(i, j)| k.sign * cross(i, j) > 0.0)
    }

    fn moved(&self, v: usize) -> bool {
        self.pos[v] != self.reference[v]
    }

    /// Edge `e` changed from its template shape.
    fn changed(&self, e: usize) -> bool {
        let edge = &self.edges[e];
        edge.bend != [0.0; 2] || self.moved(edge.a) || self.moved(edge.b)
    }

    /// The bent polyline of flexible edge `e`.
    fn bent(&self, e: usize) -> Vec<P> {
        let edge = &self.edges[e];
        let da = sub(self.pos[edge.a], self.reference[edge.a]);
        let db = sub(self.pos[edge.b], self.reference[edge.b]);
        edge.reference
            .iter()
            .zip(&edge.s)
            .map(|(r, &s)| {
                // Hermite modes s(1 − s)² and s²(1 − s), scaled to peak 1
                let h0 = 6.75 * s * (1.0 - s) * (1.0 - s);
                let h1 = 6.75 * s * s * (1.0 - s);
                let w = edge.chord * (edge.bend[0] * h0 + edge.bend[1] * h1);
                [
                    r[0] + (1.0 - s) * da[0] + s * db[0] + w * edge.normal[0],
                    r[1] + (1.0 - s) * da[1] + s * db[1] + w * edge.normal[1],
                ]
            })
            .collect()
    }

    /// The curve of edge `e` as it now stands.
    fn curve(&self, e: usize) -> Curve {
        let edge = &self.edges[e];
        if edge.kind != Kind::Flexible || !self.changed(e) {
            return edge.curve.clone();
        }
        if edge.straight && edge.bend == [0.0; 2] {
            return Curve::Line;
        }
        let p = self.bent(e);
        Curve::Spline {
            through: p[1..p.len() - 1]
                .iter()
                .map(|q| Vec3::new(q[0], q[1], self.z))
                .collect(),
        }
    }

    fn edge_nodes(&self, e: usize) -> Result<(Vec<P>, [P; 2]), String> {
        let edge = &self.edges[e];
        let (a, b) = (self.pos[edge.a], self.pos[edge.b]);
        let v3 = |p: P| Vec3::new(p[0], p[1], self.z);
        let ev = Evaluator::new(v3(a), v3(b), &self.curve(e))?;
        let f = fractions(edge.law, edge.cells, ev.length())?;
        let n = edge.cells;
        let mut nodes = Vec::with_capacity(n + 1);
        nodes.push(a);
        for &s in &f[1..n] {
            let q = ev.at(s);
            nodes.push([q.x, q.y]);
        }
        nodes.push(b);
        let tangent = |start: bool| {
            let t = ev.end_tangent(start);
            [t.x, t.y]
        };
        Ok((nodes, [tangent(true), tangent(false)]))
    }

    fn interface_cost(&self, i: usize) -> f64 {
        let f = &self.interfaces[i];
        let edge = &self.edges[f.edge];
        let n = edge.cells;
        let ([b0, s0], [b1, s1]) = ([f.sides[0].0, f.sides[0].1], [f.sides[1].0, f.sides[1].1]);
        let (k0, k1) = (&self.blocks[b0], &self.blocks[b1]);
        let (fw0, fw1) = (k0.sides[s0].1, k1.sides[s1].1);
        (0..n)
            .map(|k| {
                let c0 = k0.centroids[k0.cell(s0, if fw0 { k } else { n - 1 - k })];
                let c1 = k1.centroids[k1.cell(s1, if fw1 { k } else { n - 1 - k })];
                face_cost(edge.nodes[k], edge.nodes[k + 1], c0, c1)
            })
            .sum()
    }

    fn orthogonality(&self) -> Orthogonality {
        let mut worst: Vec<Vec<f64>> = self
            .blocks
            .iter()
            .map(|b| vec![0.0; b.n[0] * b.n[1]])
            .collect();
        let (mut sum, mut count, mut max, mut folded) = (0.0, 0usize, 0.0f64, 0usize);
        let mut max_growth = 1.0f64;
        for (bi, b) in self.blocks.iter().enumerate() {
            let (faces, f) = b.angles();
            folded += f;
            for (c0, c1, t, g) in faces {
                max_growth = max_growth.max(g);
                sum += t;
                count += 1;
                max = max.max(t);
                worst[bi][c0] = worst[bi][c0].max(t);
                worst[bi][c1] = worst[bi][c1].max(t);
            }
        }
        for f in &self.interfaces {
            let edge = &self.edges[f.edge];
            let n = edge.cells;
            let ([b0, s0], [b1, s1]) = ([f.sides[0].0, f.sides[0].1], [f.sides[1].0, f.sides[1].1]);
            let (k0, k1) = (&self.blocks[b0], &self.blocks[b1]);
            for k in 0..n {
                let i0 = k0.cell(s0, if k0.sides[s0].1 { k } else { n - 1 - k });
                let i1 = k1.cell(s1, if k1.sides[s1].1 { k } else { n - 1 - k });
                let (p, q) = (edge.nodes[k], edge.nodes[k + 1]);
                let (c0, c1) = (k0.centroids[i0], k1.centroids[i1]);
                let t = face_angle(p, q, c0, c1);
                max_growth = max_growth.max(growth(p, q, c0, c1));
                sum += t;
                count += 1;
                max = max.max(t);
                worst[b0][i0] = worst[b0][i0].max(t);
                worst[b1][i1] = worst[b1][i1].max(t);
            }
        }
        let limit = 70f64.to_radians();
        Orthogonality {
            max_deg: max.to_degrees(),
            mean_deg: if count > 0 {
                (sum / count as f64).to_degrees()
            } else {
                0.0
            },
            cells_over_70: worst.iter().flatten().filter(|&&t| t > limit).count(),
            max_growth,
            folded_corners: folded,
        }
    }

    /// The parameters: sliding and free vertices, bendable edges.
    fn dofs(&self) -> Vec<Dof> {
        let mut incident: HashMap<usize, Vec<usize>> = HashMap::new();
        for (e, edge) in self.edges.iter().enumerate() {
            incident.entry(edge.a).or_default().push(e);
            incident.entry(edge.b).or_default().push(e);
        }
        let mut corner_of: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut side_of: HashMap<usize, Vec<usize>> = HashMap::new();
        for (b, block) in self.blocks.iter().enumerate() {
            if block.rigid {
                continue;
            }
            for &v in &block.v {
                corner_of.entry(v).or_default().push(b);
            }
            for &(e, _) in &block.sides {
                side_of.entry(e).or_default().push(b);
            }
        }
        let interfaces_of = |blocks: &[usize]| -> Vec<usize> {
            (0..self.interfaces.len())
                .filter(|&i| {
                    let f = &self.interfaces[i];
                    blocks.contains(&f.sides[0].0) || blocks.contains(&f.sides[1].0)
                })
                .collect()
        };
        let rigid_vertex: BTreeSet<usize> = self
            .blocks
            .iter()
            .filter(|b| b.rigid)
            .flat_map(|b| b.v)
            .collect();
        let mut dofs = Vec::new();
        let mut vertices: Vec<usize> = corner_of.keys().copied().collect();
        vertices.sort_unstable();
        for v in vertices {
            if rigid_vertex.contains(&v) {
                continue;
            }
            let es = &incident[&v];
            let boundary: Vec<usize> = es
                .iter()
                .copied()
                .filter(|&e| self.edges[e].kind == Kind::Boundary)
                .collect();
            let dirs: Vec<P> = match boundary.as_slice() {
                [] => vec![[1.0, 0.0], [0.0, 1.0]],
                &[e0, e1] => {
                    let straight = |e: usize| self.edges[e].curve == Curve::Line;
                    let dir = |e: usize| {
                        let edge = &self.edges[e];
                        let d = sub(self.pos[edge.b], self.pos[edge.a]);
                        let l = norm(d);
                        [d[0] / l, d[1] / l]
                    };
                    let (d0, d1) = (dir(e0), dir(e1));
                    if straight(e0) && straight(e1) && cross(d0, d1).abs() < 1e-9 {
                        vec![d0]
                    } else {
                        vec![]
                    }
                }
                _ => vec![],
            };
            let blocks = corner_of[&v].clone();
            let edges: Vec<usize> = es
                .iter()
                .copied()
                .filter(|&e| self.edges[e].kind != Kind::Rigid)
                .collect();
            let shortest = edges
                .iter()
                .map(|&e| {
                    let edge = &self.edges[e];
                    norm(sub(self.pos[edge.b], self.pos[edge.a]))
                })
                .fold(f64::INFINITY, f64::min);
            for dir in dirs {
                let step = FIRST_STEP * shortest;
                dofs.push(Dof {
                    param: Param::Move { v, dir },
                    step,
                    min_step: step / STEP_RANGE,
                    max_step: 0.25 * shortest,
                    sign: 1.0,
                    edges: edges.clone(),
                    blocks: blocks.clone(),
                    interfaces: interfaces_of(&blocks),
                    sides: boundary.clone(),
                });
            }
        }
        for (e, edge) in self.edges.iter().enumerate() {
            if edge.kind != Kind::Flexible {
                continue;
            }
            let blocks = side_of.get(&e).cloned().unwrap_or_default();
            for end in 0..2 {
                if law_spacing(edge.law, end).is_none() {
                    continue;
                }
                dofs.push(Dof {
                    param: Param::Law { e, end },
                    step: 0.1,
                    min_step: 0.1 / STEP_RANGE,
                    max_step: 0.5,
                    sign: 1.0,
                    edges: vec![e],
                    blocks: blocks.clone(),
                    interfaces: interfaces_of(&blocks),
                    sides: Vec::new(),
                });
            }
            for mode in 0..2 {
                dofs.push(Dof {
                    param: Param::Bend { e, mode },
                    step: FIRST_STEP,
                    min_step: FIRST_STEP / STEP_RANGE,
                    max_step: 0.1,
                    sign: 1.0,
                    edges: vec![e],
                    blocks: blocks.clone(),
                    interfaces: interfaces_of(&blocks),
                    sides: Vec::new(),
                });
            }
        }
        dofs
    }

    /// Change parameter `d` by `delta` and keep the change if the cost
    /// drops.
    fn trial(&mut self, d: &Dof, delta: f64) -> bool {
        let old: f64 = d.blocks.iter().map(|&b| self.blocks[b].cost).sum::<f64>()
            + d.interfaces
                .iter()
                .map(|&i| self.interfaces[i].cost)
                .sum::<f64>();
        let saved_pos;
        match d.param {
            Param::Move { v, dir } => {
                saved_pos = Some((v, self.pos[v]));
                self.pos[v] = [
                    self.pos[v][0] + delta * dir[0],
                    self.pos[v][1] + delta * dir[1],
                ];
                // a sliding vertex keeps its boundary edges long enough
                for &e in &d.sides {
                    let edge = &self.edges[e];
                    let now = sub(self.pos[edge.b], self.pos[edge.a]);
                    let then = sub(self.reference[edge.b], self.reference[edge.a]);
                    if now[0] * then[0] + now[1] * then[1] < MIN_SIDE * norm(then).powi(2) {
                        self.pos[v] = saved_pos.map_or(self.pos[v], |s| s.1);
                        return false;
                    }
                }
            }
            Param::Bend { e, mode } => {
                saved_pos = None;
                let b = self.edges[e].bend[mode] + delta;
                if b.abs() > MAX_BEND {
                    return false;
                }
                self.edges[e].bend[mode] = b;
            }
            Param::Law { e, end } => {
                saved_pos = None;
                let edge = &mut self.edges[e];
                let law = scale_law(edge.law, end, delta.exp());
                let (Some(now), Some(then)) = (law_spacing(law, end), law_spacing(edge.law0, end))
                else {
                    return false;
                };
                if (now / then).ln().abs() > MAX_LAW_FACTOR.ln() {
                    return false;
                }
                edge.law = law;
            }
        }
        let restore_param = |m: &mut Self| match d.param {
            Param::Move { v, .. } => {
                if let Some((_, p)) = saved_pos {
                    m.pos[v] = p;
                }
            }
            Param::Bend { e, mode } => m.edges[e].bend[mode] -= delta,
            Param::Law { e, end } => {
                let edge = &mut m.edges[e];
                edge.law = scale_law(edge.law, end, (-delta).exp());
            }
        };
        let mut new_nodes = Vec::with_capacity(d.edges.len());
        for &e in &d.edges {
            match self.edge_nodes(e) {
                Ok(n) => new_nodes.push(n),
                Err(_) => {
                    restore_param(self);
                    return false;
                }
            }
        }
        let old_nodes: Vec<(Vec<P>, [P; 2])> = d
            .edges
            .iter()
            .zip(new_nodes)
            .map(|(&e, (n, t))| {
                let edge = &mut self.edges[e];
                (
                    std::mem::replace(&mut edge.nodes, n),
                    std::mem::replace(&mut edge.tangents, t),
                )
            })
            .collect();
        // the topology check needs every block's corners right-handed
        if !d.blocks.iter().all(|&b| self.corners_turn(b)) {
            restore_param(self);
            for (&e, (n, t)) in d.edges.iter().zip(old_nodes) {
                self.edges[e].nodes = n;
                self.edges[e].tangents = t;
            }
            return false;
        }
        let old_blocks: Vec<(Vec<P>, Vec<P>, usize, f64)> = d
            .blocks
            .iter()
            .map(|&b| {
                let k = &self.blocks[b];
                (k.nodes.clone(), k.centroids.clone(), k.folded, k.cost)
            })
            .collect();
        for &b in &d.blocks {
            let (edges, block) = (&self.edges, &mut self.blocks[b]);
            block.fill(edges);
        }
        let old_interfaces: Vec<f64> = d
            .interfaces
            .iter()
            .map(|&i| self.interfaces[i].cost)
            .collect();
        for &i in &d.interfaces {
            self.interfaces[i].cost = self.interface_cost(i);
        }
        let new: f64 = d.blocks.iter().map(|&b| self.blocks[b].cost).sum::<f64>()
            + d.interfaces
                .iter()
                .map(|&i| self.interfaces[i].cost)
                .sum::<f64>();
        if new.is_finite() && new < old - 1e-12 * old.abs() {
            return true;
        }
        restore_param(self);
        for (&e, (n, t)) in d.edges.iter().zip(old_nodes) {
            self.edges[e].nodes = n;
            self.edges[e].tangents = t;
        }
        for (&b, (n, c, folded, cost)) in d.blocks.iter().zip(old_blocks) {
            let k = &mut self.blocks[b];
            k.nodes = n;
            k.centroids = c;
            k.folded = folded;
            k.cost = cost;
        }
        for (&i, c) in d.interfaces.iter().zip(old_interfaces) {
            self.interfaces[i].cost = c;
        }
        false
    }

    /// Write the moved vertices and bent edges into the topology, on both
    /// layers.
    fn write(&self, t: &mut BlockTopology) {
        let mut twin: HashMap<usize, usize> = HashMap::new();
        for b in &t.blocks {
            for k in 0..4 {
                twin.insert(b.vertices[k], b.vertices[k + 4]);
            }
        }
        for (&v, &w) in &twin {
            if self.moved(v) {
                for u in [v, w] {
                    t.vertices[u].x = self.pos[v][0];
                    t.vertices[u].y = self.pos[v][1];
                }
            }
        }
        let mut spec: HashMap<(usize, usize), usize> = HashMap::new();
        for (i, e) in t.edges.iter().enumerate() {
            spec.insert(edge_key(e.vertices[0], e.vertices[1]), i);
        }
        for (e, edge) in self.edges.iter().enumerate() {
            if edge.kind != Kind::Flexible || !(self.changed(e) || edge.law != edge.law0) {
                continue;
            }
            let curve = self.curve(e);
            for (a, b) in [(edge.a, edge.b), (twin[&edge.a], twin[&edge.b])] {
                let z = t.vertices[a].z;
                let curve = match &curve {
                    Curve::Spline { through } => Curve::Spline {
                        through: through.iter().map(|q| Vec3::new(q.x, q.y, z)).collect(),
                    },
                    other => other.clone(),
                };
                match spec.get(&edge_key(a, b)) {
                    Some(&i) => {
                        // keep the stored direction
                        let s = &mut t.edges[i];
                        (s.curve, s.distribution) = if s.vertices[0] == a {
                            (curve, edge.law)
                        } else {
                            (reverse(curve), reverse_law(edge.law))
                        };
                    }
                    None => t.edges.push(EdgeSpec {
                        vertices: [a, b],
                        curve,
                        distribution: edge.law,
                    }),
                }
            }
        }
    }
}

fn reverse(c: Curve) -> Curve {
    match c {
        Curve::Spline { mut through } => {
            through.reverse();
            Curve::Spline { through }
        }
        other => other,
    }
}

/// Non-orthogonality of the grid of an extruded topology.
pub fn orthogonality(t: &BlockTopology) -> Result<Orthogonality, String> {
    Ok(Model::new(t, &[])?.orthogonality())
}

/// Move the vertices and bend the edges of `t` away from the rigid
/// patches so that the grid comes closer to orthogonal; see the module
/// notes. A topology that is not extruded along z is refused.
pub fn optimize(t: &mut BlockTopology, o: &LayoutOptions) -> Result<LayoutReport, String> {
    let mut model = Model::new(t, &o.rigid_patches)?;
    let before = model.orthogonality();
    let mut dofs = model.dofs();
    let mut report = LayoutReport {
        parameters: dofs.len(),
        before,
        after: before,
        ..LayoutReport::default()
    };
    for _ in 0..o.max_sweeps {
        let mut active = false;
        for d in dofs.iter_mut() {
            if d.step < d.min_step {
                continue;
            }
            active = true;
            let s = d.step * d.sign;
            report.trials += 1;
            let mut ok = model.trial(d, s);
            if !ok {
                report.trials += 1;
                ok = model.trial(d, -s);
                if ok {
                    d.sign = -d.sign;
                }
            }
            d.step = if ok {
                (1.5 * d.step).min(d.max_step)
            } else {
                0.5 * d.step
            };
        }
        if !active {
            break;
        }
        report.sweeps += 1;
    }
    report.after = model.orthogonality();
    model.write(t);
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use mesh_core::block::{BlockSpec, PatchSpec};

    /// Two blocks side by side, extruded along z, whose shared edge leans
    /// over by 0.6 of the height; `wall` names the left side a patch of its
    /// own (index 0).
    fn leaning() -> BlockTopology {
        let plane = [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 1.0],
            [1.6, 1.0],
            [2.0, 1.0],
        ];
        let mut vertices: Vec<Vec3> = plane.iter().map(|p| Vec3::new(p[0], p[1], 0.0)).collect();
        vertices.extend(plane.iter().map(|p| Vec3::new(p[0], p[1], 0.1)));
        let block = |v: [usize; 4]| BlockSpec {
            vertices: [
                v[0],
                v[1],
                v[2],
                v[3],
                v[0] + 6,
                v[1] + 6,
                v[2] + 6,
                v[3] + 6,
            ],
            cells: [8, 8, 1],
            domain: "fluid".into(),
        };
        let patch = |name: &str, faces: Vec<[usize; 4]>| PatchSpec {
            name: name.into(),
            faces,
        };
        BlockTopology {
            vertices,
            edges: Vec::new(),
            blocks: vec![block([0, 1, 4, 3]), block([1, 2, 5, 4])],
            patches: vec![
                patch("left", vec![[0, 3, 9, 6]]),
                patch("right", vec![[2, 5, 11, 8]]),
                patch("bottom", vec![[0, 1, 7, 6], [1, 2, 8, 7]]),
                patch("top", vec![[3, 4, 10, 9], [4, 5, 11, 10]]),
                patch("low", vec![[0, 1, 4, 3], [1, 2, 5, 4]]),
                patch("high", vec![[6, 7, 10, 9], [7, 8, 11, 10]]),
            ],
            surfaces: Vec::new(),
        }
    }

    /// Invariant: the search straightens a leaning block edge by sliding
    /// its ends along the straight boundary, keeps the corners and the
    /// boundary lines, writes both layers alike, reports what the written
    /// topology grids to, and gives the same result every run. A change
    /// that let a vertex leave its boundary line, moved only one layer or
    /// reported the model instead of the file breaks it; the template
    /// tests see the optimiser only through the multi-element numbers.
    #[test]
    fn layout_straightens_a_leaning_edge_along_the_boundary() {
        let mut t = leaning();
        let r = optimize(&mut t, &LayoutOptions::default()).unwrap();
        assert_eq!(
            r.parameters, 4,
            "two sliding vertices and two bend modes of the shared edge"
        );
        assert!(r.before.max_deg > 25.0, "{:?}", r.before);
        assert!(r.after.max_deg < 1.0, "{:?}", r.after);
        assert_eq!(r.after.folded_corners, 0);
        let v = &t.vertices;
        assert!((v[1].x - v[4].x).abs() < 0.02, "{} {}", v[1].x, v[4].x);
        for (i, y) in [(0, 0.0), (1, 0.0), (2, 0.0), (3, 1.0), (4, 1.0), (5, 1.0)] {
            assert_eq!(v[i].y, y);
            assert_eq!((v[i + 6].x, v[i + 6].y), (v[i].x, v[i].y));
            assert_eq!(v[i + 6].z, 0.1);
        }
        for i in [0, 2, 3, 5] {
            assert_eq!(v[i].x, leaning().vertices[i].x, "corner {i} moved");
        }
        super::super::check(&t).unwrap();
        let again = orthogonality(&t).unwrap();
        assert!((again.max_deg - r.after.max_deg).abs() < 1e-9);
        let mut u = leaning();
        optimize(&mut u, &LayoutOptions::default()).unwrap();
        assert_eq!(t, u);
    }

    /// Invariant: blocks on a rigid patch keep their vertices and edges;
    /// with the leaning edge on such a block nothing is left to move and
    /// the topology comes back unchanged. A change that moved a vertex of
    /// a rigid block (a wall's near-wall layer) breaks it.
    #[test]
    fn layout_keeps_blocks_on_rigid_patches() {
        let mut t = leaning();
        let r = optimize(
            &mut t,
            &LayoutOptions {
                rigid_patches: vec![0],
                ..LayoutOptions::default()
            },
        )
        .unwrap();
        assert_eq!(r.parameters, 0);
        assert_eq!(t, leaning());
        assert_eq!(r.before, r.after);
    }
}
