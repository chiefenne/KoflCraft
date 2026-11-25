//! Join the cubes of a hex core to the median dual (`research/hex-core.md`).
//!
//! The dual cells end on the shell with small triangles in the planes of
//! the cube sides. Each such triangle becomes a face between its dual cell
//! and the cube behind it. A side between two cubes is one square (four
//! where the cube across is half the size), with every point that lies on
//! its edges (corners of smaller cubes, dual points on the shell) inserted,
//! so all cells close edge to edge.
use super::*;
use crate::hexcore::HexCore;
use std::collections::HashMap;

type Key = [u64; 3];

fn key(v: Vec3) -> Key {
    [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()]
}

fn coords(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

/// An axis line: the axis and the bits of the two other coordinates.
type Line = (usize, [u64; 2]);

/// Points by axis line: (coordinate along the axis, point id), sorted.
#[derive(Default)]
struct Lines(HashMap<Line, Vec<(f64, usize)>>);

impl Lines {
    fn add(&mut self, id: usize, v: Vec3) {
        let c = coords(v);
        for axis in 0..3 {
            let other = [c[(axis + 1) % 3].to_bits(), c[(axis + 2) % 3].to_bits()];
            self.0.entry((axis, other)).or_default().push((c[axis], id));
        }
    }

    fn sort(&mut self) {
        for v in self.0.values_mut() {
            v.sort_by(|a, b| a.0.total_cmp(&b.0));
        }
    }

    /// The points strictly between `a` and `b` on their axis line, in order
    /// from `a`.
    fn between(&self, a: Vec3, b: Vec3) -> Vec<usize> {
        let (ca, cb) = (coords(a), coords(b));
        let Some(axis) = (0..3).find(|&i| ca[i] != cb[i]) else {
            return Vec::new();
        };
        let other = [ca[(axis + 1) % 3].to_bits(), ca[(axis + 2) % 3].to_bits()];
        let Some(line) = self.0.get(&(axis, other)) else {
            return Vec::new();
        };
        let (lo, hi) = (ca[axis].min(cb[axis]), ca[axis].max(cb[axis]));
        let start = line.partition_point(|&(t, _)| t <= lo);
        let end = line.partition_point(|&(t, _)| t < hi);
        let mut out: Vec<usize> = line[start..end].iter().map(|&(_, i)| i).collect();
        if ca[axis] > cb[axis] {
            out.reverse();
        }
        out
    }
}

pub(super) fn join(
    core: &HexCore,
    p: &mut Vec<Vec3>,
    cells: &mut Vec<Cell>,
    boundary: &mut Vec<BoundaryFace>,
) -> Result<(), String> {
    let mut faces: Vec<Vec<Vec<usize>>> = vec![Vec::new(); core.cubes.len()];
    let mut at: HashMap<Key, usize> = HashMap::new();
    let mut kept = Vec::with_capacity(boundary.len());
    let probe = 0.25 * core.unit();
    for f in boundary.drain(..) {
        if !is_hex(f.patch) {
            kept.push(f);
            continue;
        }
        let ring: Vec<Vec3> = f.nodes.iter().map(|&n| p[n]).collect();
        let centre = ring.iter().fold(Vec3::ZERO, |s, &v| s + v) / ring.len() as f64;
        let mut normal = Vec3::ZERO;
        for k in 0..ring.len() {
            normal = normal + ring[k].cross(ring[(k + 1) % ring.len()]);
        }
        let length = normal.length();
        if length.is_nan() || length <= 0.0 {
            return Err("degenerate dual face on the hex core shell".into());
        }
        // the face points out of its dual cell, into the cube
        let cube = core
            .locate(centre + normal * (probe / length))
            .ok_or("a dual face on the hex core shell has no cube behind it")?;
        faces[cube].push(f.nodes.iter().rev().copied().collect());
        for &n in &f.nodes {
            at.entry(key(p[n])).or_insert(n);
        }
    }
    *boundary = kept;
    let mut lines = Lines::default();
    for (&k, &id) in &at {
        lines.add(
            id,
            Vec3::new(
                f64::from_bits(k[0]),
                f64::from_bits(k[1]),
                f64::from_bits(k[2]),
            ),
        );
    }
    // cube corners: the shell ones are scaffold points already
    let mut corner: HashMap<[i64; 3], usize> = HashMap::new();
    for c in 0..core.cubes.len() {
        let (lo, s) = core.cube_box(c);
        for k in 0..8 {
            let q = [
                lo[0] + s * (k & 1),
                lo[1] + s * ((k >> 1) & 1),
                lo[2] + s * ((k >> 2) & 1),
            ];
            if corner.contains_key(&q) {
                continue;
            }
            let v = core.position(q);
            let id = match at.get(&key(v)) {
                Some(&id) => id,
                None => {
                    p.push(v);
                    lines.add(p.len() - 1, v);
                    p.len() - 1
                }
            };
            corner.insert(q, id);
        }
    }
    lines.sort();
    for (c, cube_faces) in faces.iter_mut().enumerate() {
        for piece in core.cube_sides(c) {
            if piece.across.is_none() {
                continue;
            }
            let mut square = piece.square.corners().to_vec();
            if !piece.positive {
                square.reverse();
            }
            let mut ring = Vec::new();
            for k in 0..4 {
                let (a, b) = (square[k], square[(k + 1) % 4]);
                ring.push(corner[&a]);
                ring.extend(lines.between(core.position(a), core.position(b)));
            }
            cube_faces.push(ring);
        }
    }
    for (c, cube_faces) in faces.into_iter().enumerate() {
        let edge = core.cube_edge(c);
        cells.push(Cell {
            element: Element::Polyhedron(Box::new(Polyhedron::new(
                cube_faces,
                core.cube_centre(c),
                None,
                edge * edge * edge,
            ))),
            layer: 0,
            region: 0,
        });
    }
    Ok(())
}
