//! Region proposal from the wall distance d and a streamwise field p (the
//! flow potential, or the Fiedler vector of a region without ports).
//!
//! 1. Clash: c = |cos(grad d, grad p)| per tetrahedron. Where c is small
//!    the distance shells (O-grid, wall layers) and the streamwise slices
//!    cross at right angles, so a swept block template fits; c >= 0.5 is a
//!    clash.
//! 2. Slabs: p is cut into slabs about two edge lengths thick along the
//!    flow, and each slab into its connected pieces. The pieces with their
//!    face adjacency are the Reeb graph of p (Reeb 1946; on meshes Ni,
//!    Garland & Hart 2004): a chain of pieces is a duct, a piece with two
//!    pieces above or below it is a junction, and a cycle is a passage
//!    around a body. Tiny pieces (below a share of an average slab) are
//!    merged into their largest neighbour first, so wall noise does not
//!    make false branches.
//! 3. Skeleton: each piece is placed at the centroid of its ridge faces,
//!    the faces where grad d flips (d grows towards the face from both
//!    sides), i.e. the medial ridge of the heat distance; a piece without
//!    a ridge face is placed at its centroid.
//! 4. Regions: junction pieces next to each other form a junction; each
//!    piece much larger than the median piece (a chamber, where the
//!    potential is almost flat) is part of a box; the other pieces of a
//!    chain form runs of clean and clashing pieces. A clean run is a sweep,
//!    or an obstacle passage when it lies on a cycle; a clashing run is a
//!    remainder for the polyhedral mesher.

#![allow(clippy::needless_range_loop)] // parallel per-tet and per-piece arrays

use super::{Fem, TetMesh, submesh};
use mesh_core::geometry::Vec3;
use std::collections::HashMap;

/// The settings of a classification.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClassifyOptions {
    /// c at and above which a tetrahedron clashes.
    pub clash: f64,
    /// Share of clashing volume at and above which a piece clashes.
    pub clash_share: f64,
    /// Number of slabs; `None`: about two edge lengths per slab.
    pub slabs: Option<usize>,
    /// Pieces below this share of an average slab are merged.
    pub min_piece: f64,
    /// A cluster of clashing tetrahedra in a sweep at or above this share
    /// of the volume is a remainder of its own (clash zone).
    pub clash_zone: f64,
    /// A piece this many times larger than the median piece is part of a
    /// box (a chamber much wider than its ducts).
    pub box_growth: f64,
}

impl Default for ClassifyOptions {
    fn default() -> Self {
        Self {
            clash: 0.5,
            clash_share: 0.3,
            slabs: None,
            min_piece: 0.05,
            box_growth: 4.0,
            clash_zone: 0.005,
        }
    }
}

/// What a region is proposed for.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Kind {
    /// A swept block template along an unbranched duct.
    Sweep,
    /// A junction template (or the polyhedral mesher).
    Junction,
    /// A passage on a loop around a body (O-grid about it).
    Obstacle,
    /// A chamber much wider than its ducts (a box template).
    Box,
    /// Neither: the polyhedral mesher with imprint.
    Remainder,
}

impl Kind {
    pub const ALL: [Kind; 5] = [
        Kind::Sweep,
        Kind::Junction,
        Kind::Obstacle,
        Kind::Box,
        Kind::Remainder,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Kind::Sweep => "sweep",
            Kind::Junction => "junction",
            Kind::Obstacle => "obstacle",
            Kind::Box => "box",
            Kind::Remainder => "remainder",
        }
    }
}

/// A connected piece of one slab.
#[derive(Clone, Debug, Default)]
pub struct Piece {
    pub slab: usize,
    pub volume: f64,
    pub clash_volume: f64,
    pub centroid: Vec3,
    /// Centroid of the ridge faces in the piece.
    pub ridge: Option<Vec3>,
    pub down: Vec<usize>,
    pub up: Vec<usize>,
}

impl Piece {
    pub fn is_junction(&self) -> bool {
        self.down.len() >= 2 || self.up.len() >= 2
    }

    pub fn is_end(&self) -> bool {
        !self.is_junction() && (self.down.is_empty() || self.up.is_empty())
    }

    /// Where the skeleton passes through the piece.
    pub fn position(&self) -> Vec3 {
        self.ridge.unwrap_or(self.centroid)
    }
}

/// A proposed region.
#[derive(Clone, Debug)]
pub struct Region {
    pub kind: Kind,
    pub pieces: Vec<usize>,
    pub volume: f64,
    pub clash_volume: f64,
    pub centroid: Vec3,
    /// First and last slab.
    pub slabs: (usize, usize),
}

/// The result of [`classify`].
#[derive(Clone, Debug, Default)]
pub struct Classification {
    /// c per tetrahedron.
    pub clash: Vec<f64>,
    pub tet_piece: Vec<usize>,
    pub tet_region: Vec<usize>,
    pub pieces: Vec<Piece>,
    pub edges: Vec<(usize, usize)>,
    pub regions: Vec<Region>,
    pub slabs: usize,
    /// Streamwise length estimate, sum V |grad p| / sum V |grad p|^2.
    pub length: f64,
    pub ridge_faces: usize,
    pub volume: f64,
    /// Volume with c below the clash bound.
    pub blockable_volume: f64,
    /// Independent loops of the piece graph.
    pub loops: usize,
}

impl Classification {
    /// Volume of each kind, in [`Kind::ALL`] order.
    pub fn kind_volumes(&self) -> [f64; 5] {
        let mut v = [0.0; 5];
        for r in &self.regions {
            v[Kind::ALL.iter().position(|&k| k == r.kind).unwrap_or(4)] += r.volume;
        }
        v
    }

    /// Junction regions (each one branch point of the skeleton).
    pub fn branch_points(&self) -> usize {
        self.regions
            .iter()
            .filter(|r| r.kind == Kind::Junction)
            .count()
    }

    pub fn end_points(&self) -> usize {
        self.pieces.iter().filter(|p| p.is_end()).count()
    }
}

struct Union(Vec<usize>);

impl Union {
    fn new(n: usize) -> Self {
        Self((0..n).collect())
    }
    fn find(&mut self, mut a: usize) -> usize {
        while self.0[a] != a {
            self.0[a] = self.0[self.0[a]];
            a = self.0[a];
        }
        a
    }
    /// Join, keeping `keep`'s root.
    fn join(&mut self, keep: usize, other: usize) {
        let (a, b) = (self.find(keep), self.find(other));
        if a != b {
            self.0[b] = a;
        }
    }
}

/// Compact labels to 0..n.
fn compact(labels: &mut [usize]) -> usize {
    let mut map = HashMap::new();
    for l in labels.iter_mut() {
        let n = map.len();
        *l = *map.entry(*l).or_insert(n);
    }
    map.len()
}

fn face_area_normal(p: [Vec3; 3]) -> (f64, Vec3) {
    let n = (p[1] - p[0]).cross(p[2] - p[0]);
    let a = n.length();
    (a / 2.0, if a > 0.0 { n / a } else { Vec3::ZERO })
}

/// Bridges of an undirected simple graph (Tarjan 1974), iteratively.
fn bridges(n: usize, edges: &[(usize, usize)]) -> Vec<bool> {
    let mut adj = vec![Vec::new(); n];
    for (e, &(a, b)) in edges.iter().enumerate() {
        adj[a].push((b, e));
        adj[b].push((a, e));
    }
    let mut disc = vec![usize::MAX; n];
    let mut low = vec![0; n];
    let mut bridge = vec![false; edges.len()];
    let mut time = 0;
    for root in 0..n {
        if disc[root] != usize::MAX {
            continue;
        }
        // (node, parent edge, next adjacency index)
        let mut stack = vec![(root, usize::MAX, 0usize)];
        disc[root] = time;
        low[root] = time;
        time += 1;
        while let Some(&mut (v, pe, ref mut i)) = stack.last_mut() {
            if *i < adj[v].len() {
                let (w, e) = adj[v][*i];
                *i += 1;
                if e == pe {
                    continue;
                }
                if disc[w] == usize::MAX {
                    disc[w] = time;
                    low[w] = time;
                    time += 1;
                    stack.push((w, e, 0));
                } else {
                    low[v] = low[v].min(disc[w]);
                }
            } else {
                stack.pop();
                if let Some(&(u, _, _)) = stack.last() {
                    low[u] = low[u].min(low[v]);
                    if low[v] > disc[u] {
                        bridge[pe] = true;
                    }
                }
            }
        }
    }
    bridge
}

/// Classify the tetrahedra of `mesh` from node fields `d` (wall
/// distance) and `p` (streamwise, about 0 to 1).
pub fn classify(
    mesh: &TetMesh,
    fem: &Fem,
    adjacent: &[[Option<usize>; 4]],
    d: &[f64],
    p: &[f64],
    opts: ClassifyOptions,
) -> Classification {
    let nt = mesh.tets.len();
    let gd = fem.gradient(mesh, d);
    let gp = fem.gradient(mesh, p);
    let mut gp_len: Vec<f64> = gp.iter().map(|g| g.length()).collect();
    let (mut s1, mut s2) = (0.0, 0.0);
    for k in 0..nt {
        s1 += fem.volume[k] * gp_len[k];
        s2 += fem.volume[k] * gp_len[k] * gp_len[k];
    }
    let length = if s2 > 0.0 { s1 / s2 } else { 0.0 };
    let total: f64 = fem.volume.iter().sum();
    gp_len.sort_unstable_by(f64::total_cmp);
    let gp_ref = gp_len.get(nt / 2).copied().unwrap_or(0.0);
    let clash: Vec<f64> = (0..nt)
        .map(|k| {
            let (a, b) = (gd[k], gp[k]);
            let (la, lb) = (a.length(), b.length());
            if la < 0.3 {
                0.0 // on the medial ridge: the core of an O-grid
            } else if lb <= 1e-3 * gp_ref {
                1.0 // no streamwise direction (a dead pocket)
            } else {
                (a.dot(b) / (la * lb)).abs()
            }
        })
        .collect();
    let blockable_volume: f64 = (0..nt)
        .filter(|&k| clash[k] < opts.clash)
        .map(|k| fem.volume[k])
        .sum();
    let slabs = opts
        .slabs
        .unwrap_or_else(|| {
            if fem.h > 0.0 {
                (length / (2.0 * fem.h)).round() as usize
            } else {
                1
            }
        })
        .clamp(1, 512);
    let (lo, hi) = p
        .iter()
        .zip(&fem.used)
        .filter(|(_, u)| **u)
        .fold((f64::MAX, f64::MIN), |(lo, hi), (&v, _)| {
            (lo.min(v), hi.max(v))
        });
    let span = (hi - lo).max(f64::MIN_POSITIVE);
    let slab: Vec<usize> = mesh
        .tets
        .iter()
        .map(|t| {
            let pc = t.iter().map(|&v| p[v]).sum::<f64>() / 4.0;
            (((pc - lo) / span * slabs as f64).floor().max(0.0) as usize).min(slabs - 1)
        })
        .collect();
    // pieces: connected parts of each slab
    let mut u = Union::new(nt);
    for k in 0..nt {
        for o in adjacent[k].iter().flatten() {
            if slab[*o] == slab[k] {
                u.join(k, *o);
            }
        }
    }
    let mut label: Vec<usize> = (0..nt).map(|k| u.find(k)).collect();
    let mut count = compact(&mut label);
    let mut piece_slab = vec![0; count];
    let mut piece_volume = vec![0.0; count];
    for k in 0..nt {
        piece_slab[label[k]] = slab[k];
        piece_volume[label[k]] += fem.volume[k];
    }
    // merge tiny pieces into their neighbour with the most shared faces,
    // then pieces of one slab that touch
    let small = opts.min_piece * total / slabs as f64;
    for _ in 0..20 {
        let mut shared: HashMap<(usize, usize), usize> = HashMap::new();
        for k in 0..nt {
            for o in adjacent[k].iter().flatten() {
                let (a, b) = (label[k], label[*o]);
                if a != b {
                    *shared.entry((a, b)).or_default() += 1;
                }
            }
        }
        let mut best: Vec<Option<(usize, usize)>> = vec![None; count];
        for (&(a, b), &n) in &shared {
            if best[a].is_none_or(|(m, c)| n > m || (n == m && b < c)) {
                best[a] = Some((n, b));
            }
        }
        let mut pu = Union::new(count);
        let mut changed = false;
        let mut order: Vec<usize> = (0..count).collect();
        order.sort_unstable_by(|&a, &b| piece_volume[a].total_cmp(&piece_volume[b]));
        for &a in &order {
            if piece_volume[a] < small
                && let Some((_, b)) = best[a]
                && pu.find(a) != pu.find(b)
            {
                pu.join(b, a);
                changed = true;
            }
        }
        for &(a, b) in shared.keys() {
            if piece_slab[a] == piece_slab[b] && pu.find(a) != pu.find(b) {
                pu.join(a, b);
                changed = true;
            }
        }
        if !changed {
            break;
        }
        let roots: Vec<usize> = (0..count).map(|a| pu.find(a)).collect();
        for l in label.iter_mut() {
            *l = roots[*l];
        }
        count = compact(&mut label);
        // the slab of a merged piece is the slab holding most of its volume
        let mut votes: HashMap<(usize, usize), f64> = HashMap::new();
        let mut pv = vec![0.0; count];
        for k in 0..nt {
            *votes.entry((label[k], slab[k])).or_default() += fem.volume[k];
            pv[label[k]] += fem.volume[k];
        }
        let mut top = vec![(0.0, usize::MAX); count];
        for (&(l, s), &v) in &votes {
            if v > top[l].0 || (v == top[l].0 && s < top[l].1) {
                top[l] = (v, s);
            }
        }
        piece_slab = top.iter().map(|&(_, s)| s).collect();
        piece_volume = pv;
    }
    // piece data
    let mut pieces: Vec<Piece> = (0..count)
        .map(|a| Piece {
            slab: piece_slab[a],
            ..Default::default()
        })
        .collect();
    let mut ridge_sum = vec![(Vec3::ZERO, 0.0); count];
    let mut ridge_faces = 0;
    let mut links: HashMap<(usize, usize), usize> = HashMap::new();
    for k in 0..nt {
        let t = mesh.tets[k];
        let c = t.iter().fold(Vec3::ZERO, |s, &v| s + mesh.points[v]) / 4.0;
        let piece = &mut pieces[label[k]];
        piece.volume += fem.volume[k];
        piece.centroid = piece.centroid + c * fem.volume[k];
        if clash[k] >= opts.clash {
            piece.clash_volume += fem.volume[k];
        }
        for a in 0..4 {
            let Some(o) = adjacent[k][a] else { continue };
            if label[o] != label[k] {
                *links
                    .entry((label[k].min(label[o]), label[k].max(label[o])))
                    .or_default() += 1;
            }
            if o < k {
                continue;
            }
            let f = super::face_key(&t, a);
            let fp = f.map(|v| mesh.points[v]);
            let (area, mut n) = face_area_normal(fp);
            if n.dot(fp[0] - mesh.points[t[a]]) < 0.0 {
                n = Vec3::ZERO - n; // from k towards o
            }
            if gd[k].dot(n) > 0.0 && gd[o].dot(n) < 0.0 {
                ridge_faces += 1;
                let fc = (fp[0] + fp[1] + fp[2]) / 3.0;
                let s = &mut ridge_sum[label[k]];
                s.0 = s.0 + fc * area;
                s.1 += area;
            }
        }
    }
    for (a, piece) in pieces.iter_mut().enumerate() {
        if piece.volume > 0.0 {
            piece.centroid = piece.centroid / piece.volume;
        }
        let (s, w) = ridge_sum[a];
        if w > 0.0 {
            piece.ridge = Some(s / w);
        }
    }
    let mut edges: Vec<(usize, usize)> = links.keys().copied().collect();
    edges.sort_unstable();
    for &(a, b) in &edges {
        let (lo_p, hi_p) = if pieces[a].slab <= pieces[b].slab {
            (a, b)
        } else {
            (b, a)
        };
        pieces[lo_p].up.push(hi_p);
        pieces[hi_p].down.push(lo_p);
    }
    let bridge = bridges(count, &edges);
    let mut on_loop = vec![false; count];
    for (e, &(a, b)) in edges.iter().enumerate() {
        if !bridge[e] {
            on_loop[a] = true;
            on_loop[b] = true;
        }
    }
    let mut components = Union::new(count);
    for &(a, b) in &edges {
        components.join(a, b);
    }
    let parts = (0..count).filter(|&a| components.find(a) == a).count();
    let loops = edges.len() + parts - count;
    // regions
    let junction: Vec<bool> = pieces.iter().map(Piece::is_junction).collect();
    let clean: Vec<bool> = pieces
        .iter()
        .map(|p| p.volume > 0.0 && p.clash_volume / p.volume < opts.clash_share)
        .collect();
    let mut volumes: Vec<f64> = (0..count)
        .filter(|&a| !junction[a])
        .map(|a| pieces[a].volume)
        .collect();
    volumes.sort_unstable_by(f64::total_cmp);
    let reference = volumes.get(volumes.len() / 2).copied().unwrap_or(f64::MAX);
    // 0 junction, 1 box (much larger than the median piece), 2 clean, 3 clash
    let category: Vec<u8> = (0..count)
        .map(|a| {
            if junction[a] {
                0
            } else if pieces[a].volume >= opts.box_growth * reference {
                1
            } else if clean[a] {
                2
            } else {
                3
            }
        })
        .collect();
    let mut groups = Union::new(count);
    for &(a, b) in &edges {
        if category[a] == category[b] {
            groups.join(a, b);
        }
    }
    let mut by_group: HashMap<usize, Vec<usize>> = HashMap::new();
    for a in 0..count {
        by_group.entry(groups.find(a)).or_default().push(a);
    }
    let mut keys: Vec<usize> = by_group.keys().copied().collect();
    keys.sort_unstable();
    let mut regions = Vec::new();
    let mut piece_region = vec![0; count];
    for key in keys {
        let mut members = by_group.remove(&key).unwrap_or_default();
        members.sort_unstable_by_key(|&a| (pieces[a].slab, a));
        let volume: f64 = members.iter().map(|&a| pieces[a].volume).sum();
        let clash_volume: f64 = members.iter().map(|&a| pieces[a].clash_volume).sum();
        let centroid = members.iter().fold(Vec3::ZERO, |s, &a| {
            s + pieces[a].centroid * pieces[a].volume
        }) / volume.max(f64::MIN_POSITIVE);
        let first = members[0];
        let kind = match category[first] {
            0 => Kind::Junction,
            1 => Kind::Box,
            2 if members.iter().any(|&a| on_loop[a]) => Kind::Obstacle,
            2 => Kind::Sweep,
            _ => Kind::Remainder,
        };
        for &a in &members {
            piece_region[a] = regions.len();
        }
        regions.push(Region {
            kind,
            slabs: (
                pieces[first].slab,
                members.iter().map(|&a| pieces[a].slab).max().unwrap_or(0),
            ),
            pieces: members,
            volume,
            clash_volume,
            centroid,
        });
    }
    let mut tet_region: Vec<usize> = label.iter().map(|&l| piece_region[l]).collect();
    // clash zones: a connected cluster of clashing tetrahedra inside a
    // sweep or obstacle region, larger than the bound, is a remainder of
    // its own (a slab averages a local clash away, e.g. a tube row)
    let zone_min = opts.clash_zone * total;
    let mut zu = Union::new(nt);
    let clashing = |k: usize| clash[k] >= opts.clash;
    for k in 0..nt {
        if !clashing(k) {
            continue;
        }
        for &o in adjacent[k].iter().flatten() {
            if clashing(o) && tet_region[o] == tet_region[k] {
                zu.join(k, o);
            }
        }
    }
    let mut zones: HashMap<usize, Vec<usize>> = HashMap::new();
    for k in 0..nt {
        if clashing(k) && matches!(regions[tet_region[k]].kind, Kind::Sweep | Kind::Obstacle) {
            zones.entry(zu.find(k)).or_default().push(k);
        }
    }
    let mut zone_list: Vec<Vec<usize>> = zones
        .into_values()
        .filter(|z| z.iter().map(|&k| fem.volume[k]).sum::<f64>() >= zone_min)
        .collect();
    zone_list.sort_unstable_by_key(|z| z[0]);
    for z in zone_list {
        let id = regions.len();
        let (mut volume, mut centroid) = (0.0, Vec3::ZERO);
        let (mut s0, mut s1) = (usize::MAX, 0);
        for &k in &z {
            let v = fem.volume[k];
            let host = &mut regions[tet_region[k]];
            host.volume -= v;
            host.clash_volume -= v;
            volume += v;
            let t = mesh.tets[k];
            centroid = centroid + t.iter().fold(Vec3::ZERO, |a, &q| a + mesh.points[q]) * (v / 4.0);
            s0 = s0.min(slab[k]);
            s1 = s1.max(slab[k]);
            tet_region[k] = id;
        }
        regions.push(Region {
            kind: Kind::Remainder,
            pieces: Vec::new(),
            volume,
            clash_volume: volume,
            centroid: centroid / volume.max(f64::MIN_POSITIVE),
            slabs: (s0, s1),
        });
    }
    Classification {
        clash,
        tet_piece: label,
        tet_region,
        pieces,
        edges,
        regions,
        slabs,
        length,
        ridge_faces,
        volume: total,
        blockable_volume,
        loops,
    }
}

/// The second pass on one region: its own mesh, its Fiedler vector as the
/// streamwise field and the original wall distance.
pub struct SubPass {
    pub region: usize,
    pub nodes: usize,
    pub result: Classification,
    pub solve: super::FieldStats,
}

/// Classify each region of the given kinds again with its Fiedler vector
/// in place of the flow potential. (With its cut faces as inlets and
/// outlets the potential would be the same as before, because the cuts
/// lie on level sets of the potential, so nothing new could be found.)
pub fn second_pass(
    mesh: &TetMesh,
    adjacent: &[[Option<usize>; 4]],
    d: &[f64],
    first: &Classification,
    kinds: &[Kind],
    opts: ClassifyOptions,
) -> Result<Vec<SubPass>, String> {
    let mut out = Vec::new();
    for (r, region) in first.regions.iter().enumerate() {
        if !kinds.contains(&region.kind) {
            continue;
        }
        let keep: Vec<bool> = first.tet_region.iter().map(|&t| t == r).collect();
        if keep.iter().filter(|&&k| k).count() < 24 {
            continue;
        }
        let (sub, old) = submesh(mesh, &keep, adjacent);
        let fem = Fem::new(&sub)?;
        let (f, solve) = super::fiedler(&sub, &fem, 60);
        let sd: Vec<f64> = old.iter().map(|&v| d[v]).collect();
        let sub_adjacent = super::neighbours(&sub.tets);
        let result = classify(&sub, &fem, &sub_adjacent, &sd, &f, opts);
        out.push(SubPass {
            region: r,
            nodes: sub.points.len(),
            result,
            solve,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::super::tests::box_mesh;
    use super::super::{heat_distance, neighbours, potential};
    use super::*;

    /// Keep: a straight duct is one sweep with two skeleton ends and no
    /// branch; a wrong slab adjacency or piece merge makes false
    /// junctions, which no other test sees in a simple shape.
    #[test]
    fn straight_duct_is_one_sweep() {
        let m = box_mesh(6, 4);
        let fem = Fem::new(&m).unwrap();
        let adjacent = neighbours(&m.tets);
        let (d, _) = heat_distance(&m, &fem, 3.0).unwrap();
        let (p, _) = potential(&m, &fem).unwrap();
        let c = classify(&m, &fem, &adjacent, &d, &p, ClassifyOptions::default());
        assert_eq!(c.branch_points(), 0);
        assert_eq!(c.end_points(), 2);
        assert_eq!(c.loops, 0);
        let v = c.kind_volumes();
        assert!(v[0] / c.volume > 0.9, "sweep share {}", v[0] / c.volume);
    }

    /// Keep: bridges of a triangle with a tail are only the tail edge.
    #[test]
    fn bridges_of_a_cycle_with_a_tail() {
        let b = bridges(4, &[(0, 1), (1, 2), (2, 0), (2, 3)]);
        assert_eq!(b, vec![false, false, false, true]);
    }
}
