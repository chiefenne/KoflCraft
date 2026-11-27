//! T-junctions of a block topology (BS-21).
//!
//! A T-junction is a block vertex that sits inside an edge of another
//! block: the edge a–b of the block on one side passes through t, while on
//! the other side a–t and t–b are edges of two smaller blocks. A block face
//! whose region the faces of other blocks cover in strips is a *coarse*
//! face; the strips are its *sub-faces* and the edges between them its
//! *cuts*. The grid builder grids a coarse face first, by 2D transfinite
//! interpolation of its whole edges, and the sub-faces and cuts read its
//! nodes, so the hexahedra on both sides still meet node for node.
//!
//! The cell count of an edge with T-junctions is the sum of the counts of
//! its pieces. Every other rule of the topology stays: opposite edges of a
//! block carry one count, faces shared whole match.
//!
//! This is the T-mesh of quad and hex layout quantisation (Lyon, Campen and
//! Kobbelt, "Quad layouts via constrained T-mesh quantization", Computer
//! Graphics Forum 40(2), 2021; Brückler, Bommes and Campen, "Volume
//! parametrization quantization for hexahedral meshing", ACM TOG 41(4),
//! 2022), restricted to faces split in strips. No code was read.

use std::collections::{HashMap, HashSet};

use mesh_core::block::{BLOCK_FACES, BlockTopology, face_key};

use super::{Issue, edge_key, face_vertices, issue};

/// A coarse face and the sub-faces of other blocks that cover it.
#[derive(Clone, Debug, PartialEq)]
pub struct Cover {
    pub block: usize,
    /// Index into `BLOCK_FACES`.
    pub face: usize,
    /// The coarse face's corners in its local order (00, 10, 11, 01).
    pub corners: [usize; 4],
    /// True when the cuts run along the face's second free axis, so the
    /// strips follow each other along the first.
    pub along_first: bool,
    /// The sub-faces in the coarse face's local orientation (00, 10, 11,
    /// 01), in order along the strips.
    pub subs: Vec<[usize; 4]>,
    /// The cuts between neighbouring sub-faces, from the face's first
    /// side to the opposite one.
    pub cuts: Vec<(usize, usize)>,
}

/// The T-junctions of a topology, read and checked.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Junctions {
    /// Edge key → its vertices from `key.0` to `key.1`, ends included.
    pub through: HashMap<(usize, usize), Vec<usize>>,
    pub covers: Vec<Cover>,
    /// Face key → cover, for coarse faces and for sub-faces.
    pub coarse: HashMap<[usize; 4], usize>,
    pub sub: HashMap<[usize; 4], usize>,
    /// Cut edge key → cover.
    pub cuts: HashMap<(usize, usize), usize>,
}

impl Junctions {
    /// Read the junction list: valid vertices, each edge once, no vertex
    /// twice, no chain that comes back to itself.
    pub fn read(t: &BlockTopology, issues: &mut Vec<Issue>) -> Self {
        let nv = t.vertices.len();
        let mut j = Self::default();
        for (n, spec) in t.junctions.iter().enumerate() {
            let [a, b] = spec.edge;
            let mut chain = vec![a];
            chain.extend(&spec.through);
            chain.push(b);
            if chain.iter().any(|&v| v >= nv) {
                issues.push(issue(
                    "junction",
                    format!("junction {n} names a vertex that does not exist"),
                ));
                continue;
            }
            if spec.through.is_empty() {
                issues.push(issue(
                    "junction",
                    format!("junction {n} on edge {a}–{b} passes through no vertex"),
                ));
                continue;
            }
            if (0..chain.len()).any(|i| chain[..i].contains(&chain[i])) {
                issues.push(issue(
                    "junction",
                    format!("junction {n} on edge {a}–{b} repeats a vertex"),
                ));
                continue;
            }
            let key = edge_key(a, b);
            if a != key.0 {
                chain.reverse();
            }
            if j.through.insert(key, chain).is_some() {
                issues.push(issue(
                    "junction",
                    format!("edge {a}–{b} has two junction lists"),
                ));
            }
        }
        // nesting must end
        let keys: Vec<_> = j.through.keys().copied().collect();
        for key in keys {
            if j.expand_checked(key.0, key.1, 0).is_none() {
                issues.push(issue(
                    "junction",
                    format!(
                        "the junctions of edge {}–{} contain themselves",
                        key.0, key.1
                    ),
                ));
                j.through.remove(&key);
            }
        }
        j
    }

    /// The vertices of the edge `a` → `b` one level down: its junction
    /// list, or the two ends.
    pub fn chain(&self, a: usize, b: usize) -> Vec<usize> {
        let key = edge_key(a, b);
        match self.through.get(&key) {
            Some(c) => {
                let mut c = c.clone();
                if a != key.0 {
                    c.reverse();
                }
                c
            }
            None => vec![a, b],
        }
    }

    /// Every vertex on the edge `a` → `b`, nested junctions expanded.
    pub fn expand(&self, a: usize, b: usize) -> Vec<usize> {
        self.expand_checked(a, b, 0).unwrap_or_else(|| vec![a, b])
    }

    fn expand_checked(&self, a: usize, b: usize, depth: usize) -> Option<Vec<usize>> {
        if depth > 64 {
            return None;
        }
        let c = self.chain(a, b);
        if c.len() == 2 {
            return Some(c);
        }
        let mut out = vec![a];
        for w in c.windows(2) {
            let part = self.expand_checked(w[0], w[1], depth + 1)?;
            out.extend(&part[1..]);
        }
        Some(out)
    }

    /// Find the coarse faces: block faces of one block, in no patch, that
    /// the faces of other blocks cover in strips.
    pub fn find_covers(
        &mut self,
        t: &BlockTopology,
        valid_blocks: &[usize],
        faces: &HashMap<[usize; 4], Vec<usize>>,
        in_patch: &HashSet<[usize; 4]>,
        issues: &mut Vec<Issue>,
    ) {
        if self.through.is_empty() {
            return;
        }
        for &b in valid_blocks {
            for f in 0..6 {
                let c = face_vertices(&t.blocks[b].vertices, f);
                let key = face_key(c);
                if faces.get(&key).map(Vec::len) != Some(1) || in_patch.contains(&key) {
                    continue;
                }
                let sides = [(c[0], c[1]), (c[3], c[2]), (c[0], c[3]), (c[1], c[2])];
                if sides
                    .iter()
                    .all(|&(p, q)| !self.through.contains_key(&edge_key(p, q)))
                {
                    continue;
                }
                let found = [true, false]
                    .into_iter()
                    .find_map(|first| self.strips(c, first, b, faces).map(|s| (first, s)));
                let Some((along_first, subs)) = found else {
                    continue;
                };
                let n = self.covers.len();
                let mut cuts = Vec::new();
                for s in &subs[..subs.len() - 1] {
                    let cut = if along_first {
                        (s[1], s[2])
                    } else {
                        (s[3], s[2])
                    };
                    cuts.push(cut);
                    if self.through.contains_key(&edge_key(cut.0, cut.1)) {
                        issues.push(issue(
                            "junction",
                            format!(
                                "the cut {}–{} inside face {key:?} has junctions itself; split the face in strips only",
                                cut.0, cut.1
                            ),
                        ));
                    }
                    if self.cuts.insert(edge_key(cut.0, cut.1), n).is_some() {
                        issues.push(issue(
                            "junction",
                            format!("the cut {}–{} lies inside two faces", cut.0, cut.1),
                        ));
                    }
                }
                for s in &subs {
                    if self.sub.insert(face_key(*s), n).is_some() {
                        issues.push(issue(
                            "junction",
                            format!("face {:?} covers two coarse faces", face_key(*s)),
                        ));
                    }
                }
                self.coarse.insert(key, n);
                self.covers.push(Cover {
                    block: b,
                    face: f,
                    corners: c,
                    along_first,
                    subs,
                    cuts,
                });
            }
        }
        for (k, &n) in &self.sub {
            if self.coarse.contains_key(k) {
                issues.push(issue(
                    "junction",
                    format!(
                        "face {k:?} is coarse and a sub-face of block {}",
                        self.covers[n].block
                    ),
                ));
            }
        }
    }

    /// The strips of other blocks' faces that tile the face `c` of block
    /// `b`, or `None`. With `first`, the strips follow each other along
    /// the face's first free axis (between its sides c0–c1 and c3–c2).
    fn strips(
        &self,
        c: [usize; 4],
        first: bool,
        b: usize,
        faces: &HashMap<[usize; 4], Vec<usize>>,
    ) -> Option<Vec<[usize; 4]>> {
        let (lo, hi) = if first {
            (self.expand(c[0], c[1]), self.expand(c[3], c[2]))
        } else {
            (self.expand(c[0], c[3]), self.expand(c[1], c[2]))
        };
        if lo.len() == 2 && hi.len() == 2 {
            return None;
        }
        let other = |s: [usize; 4]| {
            faces
                .get(&face_key(s))
                .is_some_and(|v| v.len() == 1 && v[0] != b)
        };
        let (mut i, mut j) = (0, 0);
        let mut subs = Vec::new();
        while i + 1 < lo.len() || j + 1 < hi.len() {
            let mut next = None;
            'search: for i2 in i + 1..lo.len() {
                for j2 in j + 1..hi.len() {
                    let s = if first {
                        [lo[i], lo[i2], hi[j2], hi[j]]
                    } else {
                        [lo[i], hi[j], hi[j2], lo[i2]]
                    };
                    if other(s) {
                        next = Some((i2, j2, s));
                        break 'search;
                    }
                }
            }
            let (i2, j2, s) = next?;
            // a strip's own side must be the same run of vertices
            let side_ok = |p: usize, q: usize, run: &[usize]| self.expand(p, q) == run;
            if !side_ok(lo[i], lo[i2], &lo[i..=i2]) || !side_ok(hi[j], hi[j2], &hi[j..=j2]) {
                return None;
            }
            subs.push(s);
            (i, j) = (i2, j2);
        }
        (subs.len() >= 2).then_some(subs)
    }

    /// True when a topology has T-junctions.
    pub fn any(&self) -> bool {
        !self.through.is_empty()
    }

    /// For each block, which of its six faces is a coarse face or a
    /// sub-face.
    pub fn block_faces(&self, t: &BlockTopology) -> Vec<[bool; 6]> {
        t.blocks
            .iter()
            .map(|block| {
                std::array::from_fn(|f| {
                    let key = face_key(face_vertices(&block.vertices, f));
                    self.coarse.contains_key(&key) || self.sub.contains_key(&key)
                })
            })
            .collect()
    }

    /// The two free axes and the counts of a cover's coarse face.
    pub fn coarse_dims(&self, t: &BlockTopology, cover: &Cover) -> [usize; 2] {
        let (_, _, a, b) = BLOCK_FACES[cover.face];
        let n = t.blocks[cover.block].cells;
        [n[a], n[b]]
    }
}
