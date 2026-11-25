//! Octree hex core for large fluid regions ("hex core with layers", the
//! poly-hexcore route; design note `research/hex-core.md`).
//!
//! Cubes of a 2:1 balanced octree fill the bulk of the domain, away from
//! the wall by the layer stack plus a margin. Their outer boundary, the
//! *shell*, goes into the scaffold as an inner surface (tag [`HEX`]); the
//! median dual fills the gap between the wall layers and the shell, and
//! [`crate::duct`] joins the dual cells to the cubes face to face.
//!
//! Lattice points are always computed from integer coordinates in units of
//! the finest edge by one formula ([`HexCore::position`]), so points on a
//! common cube edge or face plane agree exactly in the fixed coordinates.

use crate::delaunay::{HashMap, HashSet};
use crate::remesh::ClosestSet;
use crate::sizing::SizeField;
use crate::sizing::surface::{Primitive, Weighted, WeightedSet};
use mesh_core::geometry::Vec3;

/// Scaffold boundary tag of the shell triangles (no patch, no layers).
pub const HEX: usize = 1 << 30;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HexCoreOptions {
    /// Cube edge over the local scaffold spacing.
    pub size_factor: f64,
    /// Gap between the layer stack and the cubes, in cube edges.
    pub margin: f64,
    /// Octree levels below the largest cube.
    pub max_depth: u32,
}

impl Default for HexCoreOptions {
    fn default() -> Self {
        Self {
            size_factor: 1.0,
            margin: 1.0,
            max_depth: 6,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct HexCoreStats {
    pub cubes: usize,
    /// Cubes per edge length, largest first: (edge, count).
    pub sizes: Vec<(f64, usize)>,
    pub removed_for_manifold: usize,
    pub shell_triangles: usize,
}

/// A square piece of a cube side: plane axis, plane coordinate and lower
/// corner in the two other axes (cyclic order), edge, all in lattice units.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Square {
    pub axis: usize,
    pub at: i64,
    pub low: [i64; 2],
    pub size: i64,
}

impl Square {
    /// The corners, counter-clockwise seen from +axis.
    pub fn corners(&self) -> [[i64; 3]; 4] {
        let (b, c) = ((self.axis + 1) % 3, (self.axis + 2) % 3);
        let [u, v] = self.low;
        let s = self.size;
        [(u, v), (u + s, v), (u + s, v + s), (u, v + s)].map(|(x, y)| {
            let mut p = [0; 3];
            p[self.axis] = self.at;
            p[b] = x;
            p[c] = y;
            p
        })
    }
}

/// One side piece of a cube: its square, whether the outward normal is
/// +axis, and the cube across it (`None`: the shell).
#[derive(Clone, Copy, Debug)]
pub struct SidePiece {
    pub square: Square,
    pub positive: bool,
    pub across: Option<usize>,
}

enum Found {
    Leaf(u32, [i64; 3], bool),
    Finer,
    Outside,
}

pub struct HexCore {
    origin: Vec3,
    unit: f64,
    depth: u32,
    top: [i64; 3],
    /// Leaves by (level, index at that level): whether a cube of the core.
    leaves: HashMap<(u32, [i64; 3]), bool>,
    /// The cubes of the core: (level, index).
    pub cubes: Vec<(u32, [i64; 3])>,
    ids: HashMap<(u32, [i64; 3]), usize>,
    pub shell_points: Vec<Vec3>,
    /// Shell triangles, normals pointing into the core.
    pub shell_triangles: Vec<[usize; 3]>,
    pub stats: HexCoreStats,
}

/// Ray parity inside test of a closed triangulated surface: rays along +x,
/// the query point nudged off the lattice in y and z.
struct Rays {
    triangles: Vec<[Vec3; 3]>,
    min: [f64; 2],
    cell: f64,
    n: [usize; 2],
    buckets: Vec<Vec<u32>>,
    nudge: f64,
}

impl Rays {
    fn new(vertices: &[Vec3], triangles: &[[usize; 3]], min: Vec3, max: Vec3) -> Self {
        let tris: Vec<[Vec3; 3]> = triangles.iter().map(|t| t.map(|v| vertices[v])).collect();
        let side = ((tris.len() as f64).sqrt().ceil() as usize).clamp(1, 512);
        let extent = (max.y - min.y).max(max.z - min.z).max(f64::MIN_POSITIVE);
        let cell = extent / side as f64 * 1.0001;
        let n = [
            (((max.y - min.y) / cell) as usize + 1).max(1),
            (((max.z - min.z) / cell) as usize + 1).max(1),
        ];
        let mut buckets = vec![Vec::new(); n[0] * n[1]];
        let at =
            |v: f64, lo: f64, k: usize| (((v - lo) / cell).floor().max(0.) as usize).min(k - 1);
        for (i, t) in tris.iter().enumerate() {
            let (y0, y1) = (
                t[0].y.min(t[1].y).min(t[2].y),
                t[0].y.max(t[1].y).max(t[2].y),
            );
            let (z0, z1) = (
                t[0].z.min(t[1].z).min(t[2].z),
                t[0].z.max(t[1].z).max(t[2].z),
            );
            for a in at(y0, min.y, n[0])..=at(y1, min.y, n[0]) {
                for b in at(z0, min.z, n[1])..=at(z1, min.z, n[1]) {
                    buckets[a * n[1] + b].push(i as u32);
                }
            }
        }
        Self {
            triangles: tris,
            min: [min.y, min.z],
            cell,
            n,
            buckets,
            nudge: (max - min).max_abs_component() * 1e-7,
        }
    }

    fn inside(&self, p: Vec3) -> bool {
        let (y, z) = (p.y + 1.234_567 * self.nudge, p.z + 2.236_067 * self.nudge);
        let a = ((y - self.min[0]) / self.cell).floor();
        let b = ((z - self.min[1]) / self.cell).floor();
        if a < 0. || b < 0. || a as usize >= self.n[0] || b as usize >= self.n[1] {
            return false;
        }
        let mut crossings = 0usize;
        for &i in &self.buckets[a as usize * self.n[1] + b as usize] {
            let [p0, p1, p2] = self.triangles[i as usize];
            let orient = |u: Vec3, v: Vec3| (v.y - u.y) * (z - u.z) - (v.z - u.z) * (y - u.y);
            let (s0, s1, s2) = (orient(p0, p1), orient(p1, p2), orient(p2, p0));
            if !((s0 > 0. && s1 > 0. && s2 > 0.) || (s0 < 0. && s1 < 0. && s2 < 0.)) {
                continue;
            }
            // x of the triangle plane at (y, z), by barycentric weights
            let total = s0 + s1 + s2;
            let x = (p2.x * s0 + p0.x * s1 + p1.x * s2) / total;
            if x > p.x {
                crossings += 1;
            }
        }
        crossings % 2 == 1
    }
}

const NEIGHBOURS: [[i64; 3]; 26] = {
    let mut out = [[0; 3]; 26];
    let mut k = 0;
    let mut i = 0;
    while i < 27 {
        let d = [
            (i / 9) as i64 - 1,
            ((i / 3) % 3) as i64 - 1,
            (i % 3) as i64 - 1,
        ];
        if !(d[0] == 0 && d[1] == 0 && d[2] == 0) {
            out[k] = d;
            k += 1;
        }
        i += 1;
    }
    out
};

impl HexCore {
    /// Edge of a cube of `level` in lattice units.
    fn units(&self, level: u32) -> i64 {
        1 << (self.depth - level)
    }

    pub fn edge(&self, level: u32) -> f64 {
        self.unit * self.units(level) as f64
    }

    /// The point at integer lattice coordinates.
    pub fn position(&self, p: [i64; 3]) -> Vec3 {
        Vec3::new(
            self.origin.x + p[0] as f64 * self.unit,
            self.origin.y + p[1] as f64 * self.unit,
            self.origin.z + p[2] as f64 * self.unit,
        )
    }

    pub fn unit(&self) -> f64 {
        self.unit
    }

    fn centre(&self, level: u32, idx: [i64; 3]) -> Vec3 {
        let s = self.units(level);
        let lo = self.position(idx.map(|v| v * s));
        let h = self.edge(level) / 2.;
        lo + Vec3::new(h, h, h)
    }

    /// The cube's lower corner and edge in lattice units.
    pub fn cube_box(&self, cube: usize) -> ([i64; 3], i64) {
        let (level, idx) = self.cubes[cube];
        let s = self.units(level);
        (idx.map(|v| v * s), s)
    }

    pub fn cube_centre(&self, cube: usize) -> Vec3 {
        let (level, idx) = self.cubes[cube];
        self.centre(level, idx)
    }

    pub fn cube_edge(&self, cube: usize) -> f64 {
        self.edge(self.cubes[cube].0)
    }

    fn find(&self, level: u32, idx: [i64; 3]) -> Found {
        let s = self.units(0);
        let top = idx.map(|v| (v * self.units(level)).div_euclid(s));
        if (0..3).any(|a| top[a] < 0 || top[a] >= self.top[a]) {
            return Found::Outside;
        }
        for l in (0..=level).rev() {
            let key = (l, idx.map(|v| v >> (level - l)));
            if let Some(&hex) = self.leaves.get(&key) {
                return Found::Leaf(l, key.1, hex);
            }
        }
        Found::Finer
    }

    /// The cube of the core that contains `p`, if any.
    pub fn locate(&self, p: Vec3) -> Option<usize> {
        let f = |v: f64, o: f64| ((v - o) / self.unit).floor() as i64;
        let idx = [
            f(p.x, self.origin.x),
            f(p.y, self.origin.y),
            f(p.z, self.origin.z),
        ];
        match self.find(self.depth, idx) {
            Found::Leaf(l, i, true) => self.ids.get(&(l, i)).copied(),
            _ => None,
        }
    }

    /// Whether `p` lies inside a cube of the core (on the shell counts as
    /// neither; callers query centroids).
    pub fn contains(&self, p: Vec3) -> bool {
        self.locate(p).is_some()
    }

    /// The pieces of one side of the cube `(level, idx)`: the outward
    /// side along `axis`, `positive` or not.
    fn side(&self, level: u32, idx: [i64; 3], axis: usize, positive: bool) -> Vec<SidePiece> {
        let s = self.units(level);
        let mut n = idx;
        n[axis] += if positive { 1 } else { -1 };
        let at = if positive {
            (idx[axis] + 1) * s
        } else {
            idx[axis] * s
        };
        let (b, c) = ((axis + 1) % 3, (axis + 2) % 3);
        let square = Square {
            axis,
            at,
            low: [idx[b] * s, idx[c] * s],
            size: s,
        };
        let mut out = Vec::new();
        self.pieces(level, n, square, positive, &mut out);
        out
    }

    /// Pieces of `square` against the region of cube `(level, idx)` across
    /// it; void pieces are merged as large as the leaves allow.
    fn pieces(
        &self,
        level: u32,
        idx: [i64; 3],
        square: Square,
        positive: bool,
        out: &mut Vec<SidePiece>,
    ) {
        match self.find(level, idx) {
            Found::Outside => out.push(SidePiece {
                square,
                positive,
                across: None,
            }),
            Found::Leaf(l, i, hex) => out.push(SidePiece {
                square,
                positive,
                across: if hex {
                    self.ids.get(&(l, i)).copied()
                } else {
                    None
                },
            }),
            Found::Finer => {
                let axis = square.axis;
                let (b, c) = ((axis + 1) % 3, (axis + 2) % 3);
                let half = square.size / 2;
                let mut sub = Vec::new();
                for (db, dc) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                    let mut child = idx.map(|v| v * 2);
                    // the children of the region next to the square
                    if !positive {
                        child[axis] += 1;
                    }
                    child[b] += db;
                    child[c] += dc;
                    let piece = Square {
                        axis,
                        at: square.at,
                        low: [square.low[0] + db * half, square.low[1] + dc * half],
                        size: half,
                    };
                    self.pieces(level + 1, child, piece, positive, &mut sub);
                }
                if sub.iter().all(|p| p.across.is_none()) {
                    out.push(SidePiece {
                        square,
                        positive,
                        across: None,
                    });
                } else {
                    out.extend(sub);
                }
            }
        }
    }

    /// All side pieces of a cube of the core.
    pub fn cube_sides(&self, cube: usize) -> Vec<SidePiece> {
        let (level, idx) = self.cubes[cube];
        let mut out = Vec::new();
        for axis in 0..3 {
            for positive in [false, true] {
                out.extend(self.side(level, idx, axis, positive));
            }
        }
        out
    }

    fn index_cubes(&mut self) {
        let mut cubes: Vec<(u32, [i64; 3])> = self
            .leaves
            .iter()
            .filter(|(_, hex)| **hex)
            .map(|(k, _)| *k)
            .collect();
        cubes.sort_unstable();
        self.ids = cubes.iter().enumerate().map(|(i, k)| (*k, i)).collect();
        self.cubes = cubes;
    }

    /// Lattice corners of every cube, sorted along each axis line.
    fn corner_lines(&self) -> HashMap<(usize, [i64; 2]), Vec<i64>> {
        let mut lines: HashMap<(usize, [i64; 2]), Vec<i64>> = HashMap::default();
        let mut seen: HashSet<[i64; 3]> = HashSet::default();
        for c in 0..self.cubes.len() {
            let (lo, s) = self.cube_box(c);
            for k in 0..8 {
                let p = [
                    lo[0] + s * (k & 1),
                    lo[1] + s * ((k >> 1) & 1),
                    lo[2] + s * ((k >> 2) & 1),
                ];
                if seen.insert(p) {
                    for axis in 0..3 {
                        let key = (axis, [p[(axis + 1) % 3], p[(axis + 2) % 3]]);
                        lines.entry(key).or_default().push(p[axis]);
                    }
                }
            }
        }
        for v in lines.values_mut() {
            v.sort_unstable();
        }
        lines
    }

    /// The polygon of a square with the corners on its edges inserted,
    /// counter-clockwise about the normal (+axis if `positive`).
    pub fn square_ring(
        square: &Square,
        positive: bool,
        lines: &HashMap<(usize, [i64; 2]), Vec<i64>>,
    ) -> Vec<[i64; 3]> {
        let corners = square.corners();
        let mut ring = Vec::new();
        for k in 0..4 {
            let (a, b) = (corners[k], corners[(k + 1) % 4]);
            ring.push(a);
            let axis = (0..3).find(|&i| a[i] != b[i]).unwrap_or(0);
            let key = (axis, [a[(axis + 1) % 3], a[(axis + 2) % 3]]);
            if let Some(line) = lines.get(&key) {
                let (lo, hi) = (a[axis].min(b[axis]), a[axis].max(b[axis]));
                let start = line.partition_point(|&v| v <= lo);
                let end = line.partition_point(|&v| v < hi);
                let mut between: Vec<i64> = line[start..end].to_vec();
                if a[axis] > b[axis] {
                    between.reverse();
                }
                for v in between {
                    let mut p = a;
                    p[axis] = v;
                    ring.push(p);
                }
            }
        }
        if !positive {
            ring.reverse();
        }
        ring
    }

    /// Shell polygons (normal into the core) with their owner cube.
    fn shell_polygons(&self) -> Vec<(Vec<[i64; 3]>, usize)> {
        let lines = self.corner_lines();
        let mut out = Vec::new();
        for c in 0..self.cubes.len() {
            for piece in self.cube_sides(c) {
                if piece.across.is_none() {
                    // outward of the cube is +axis when positive; the shell
                    // faces into the cube
                    out.push((Self::square_ring(&piece.square, !piece.positive, &lines), c));
                }
            }
        }
        out
    }

    /// Cubes at a non-manifold edge or vertex of the shell, one per defect.
    fn non_manifold(&self, polygons: &[(Vec<[i64; 3]>, usize)]) -> HashSet<usize> {
        let mut edges: HashMap<([i64; 3], [i64; 3]), Vec<usize>> = HashMap::default();
        let mut at_vertex: HashMap<[i64; 3], Vec<usize>> = HashMap::default();
        for (i, (ring, _)) in polygons.iter().enumerate() {
            for k in 0..ring.len() {
                let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
                edges.entry((a.min(b), a.max(b))).or_default().push(i);
                at_vertex.entry(a).or_default().push(i);
            }
        }
        let mut drop = HashSet::default();
        let pick = |list: &[usize]| {
            list.iter()
                .map(|&i| polygons[i].1)
                .max_by_key(|&c| (self.cubes[c].0, c))
        };
        for list in edges.values() {
            if list.len() != 2
                && let Some(c) = pick(list)
            {
                drop.insert(c);
            }
        }
        // vertex: the incident polygons must form one fan through edges
        for (v, list) in &at_vertex {
            if list.len() < 3 {
                continue;
            }
            let mut parent: HashMap<usize, usize> = list.iter().map(|&i| (i, i)).collect();
            fn root(parent: &mut HashMap<usize, usize>, mut i: usize) -> usize {
                while parent[&i] != i {
                    i = parent[&i];
                }
                i
            }
            for &i in list {
                let ring = &polygons[i].0;
                for k in 0..ring.len() {
                    let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
                    if a != *v && b != *v {
                        continue;
                    }
                    if let Some(others) = edges.get(&(a.min(b), a.max(b))) {
                        for &j in others {
                            if j != i && parent.contains_key(&j) {
                                let (x, y) = (root(&mut parent, i), root(&mut parent, j));
                                parent.insert(x, y);
                            }
                        }
                    }
                }
            }
            let mut roots: Vec<usize> = list.iter().map(|&i| root(&mut parent, i)).collect();
            roots.sort_unstable();
            roots.dedup();
            if roots.len() > 1 {
                // drop a cube of the smallest fan
                let mut best: Option<(usize, usize)> = None;
                for &r in &roots {
                    let fan: Vec<usize> = list
                        .iter()
                        .copied()
                        .filter(|&i| root(&mut parent, i) == r)
                        .collect();
                    if best.is_none_or(|(n, _)| fan.len() < n) {
                        best = pick(&fan).map(|c| (fan.len(), c));
                    }
                }
                if let Some((_, c)) = best {
                    drop.insert(c);
                }
            }
        }
        drop
    }
}

/// The scaffold spacing at a point: wall size plus `rate` times the
/// distance, capped by `max_size` and the grading field (as in
/// `scaffold::build_graded`).
fn spacing(vertices: &[Vec3], triangles: &[[usize; 3]], max_size: f64) -> WeightedSet {
    WeightedSet::new(
        triangles
            .iter()
            .map(|t| {
                let p = t.map(|v| vertices[v]);
                let l = ((p[0] - p[1]).length() + (p[1] - p[2]).length() + (p[2] - p[0]).length())
                    / 3.0;
                Weighted {
                    primitive: Primitive::Triangle(p),
                    size: l.min(max_size),
                    reach: 0.0,
                }
            })
            .collect(),
    )
}

/// Build the hex core of the closed wall surface (outward triangles).
/// `wall_offset` is the height of the layer stack.
#[allow(clippy::too_many_arguments)]
pub fn build(
    vertices: &[Vec3],
    triangles: &[[usize; 3]],
    max_size: f64,
    rate: f64,
    field: Option<&SizeField>,
    wall_offset: f64,
    opts: HexCoreOptions,
) -> Result<HexCore, String> {
    if vertices.is_empty() || triangles.is_empty() {
        return Err("hex core needs a wall".into());
    }
    if !(opts.size_factor > 0.0
        && opts.size_factor.is_finite()
        && opts.margin >= 0.0
        && opts.margin.is_finite()
        && (1..=12).contains(&opts.max_depth)
        && max_size > 0.0
        && max_size.is_finite()
        && wall_offset >= 0.0
        && wall_offset.is_finite())
    {
        return Err("invalid hex core options".into());
    }
    let mut min = vertices[0];
    let mut max = min;
    for v in vertices {
        min = Vec3::new(min.x.min(v.x), min.y.min(v.y), min.z.min(v.z));
        max = Vec3::new(max.x.max(v.x), max.y.max(v.y), max.z.max(v.z));
    }
    let top_edge = opts.size_factor * max_size;
    let depth = opts.max_depth;
    let unit = top_edge / (1u64 << depth) as f64;
    let extent = max - min;
    let top = [extent.x, extent.y, extent.z].map(|e| ((e / top_edge).ceil() as i64).max(1));
    if top.iter().product::<i64>() > 50_000_000 {
        return Err("hex core lattice too large".into());
    }
    let wall = ClosestSet::new(
        triangles
            .iter()
            .map(|t| Primitive::Triangle(t.map(|v| vertices[v])))
            .collect(),
    );
    let size = spacing(vertices, triangles, max_size);
    let rays = Rays::new(vertices, triangles, min, max);
    let target = |x: Vec3| {
        let h = size.envelope(x, rate, 0.0, max_size);
        opts.size_factor * field.map_or(h, |f| h.min(f.size(x)))
    };
    let mut core = HexCore {
        origin: min,
        unit,
        depth,
        top,
        leaves: HashMap::default(),
        cubes: Vec::new(),
        ids: HashMap::default(),
        shell_points: Vec::new(),
        shell_triangles: Vec::new(),
        stats: HexCoreStats::default(),
    };
    let distance = |p: Vec3| wall.closest(p).map_or(f64::INFINITY, |(_, d)| d);
    // whether the cube (level, idx) may join the core
    let qualifies = |core: &HexCore, level: u32, idx: [i64; 3]| {
        let c = core.centre(level, idx);
        let e = core.edge(level);
        let r = e * 3f64.sqrt() / 2.;
        let d = distance(c);
        (
            d >= wall_offset + opts.margin * e + r && rays.inside(c),
            d,
            r,
            c,
        )
    };
    // 1. refine from the top lattice
    let mut stack = Vec::new();
    for i in 0..top[0] {
        for j in 0..top[1] {
            for k in 0..top[2] {
                stack.push((0u32, [i, j, k]));
            }
        }
    }
    while let Some((level, idx)) = stack.pop() {
        let (hex, d, r, c) = qualifies(&core, level, idx);
        let e = core.edge(level);
        let t = target(c);
        let outside = d > r && !hex && !rays.inside(c);
        // a child (half the edge) can only qualify this far from the wall
        let child_may = d + r / 2. >= wall_offset + opts.margin * e / 2. + r / 2.;
        let split = level < depth && !outside && child_may && e > t * std::f64::consts::SQRT_2;
        if split {
            for k in 0..8 {
                stack.push((
                    level + 1,
                    [
                        idx[0] * 2 + (k & 1),
                        idx[1] * 2 + ((k >> 1) & 1),
                        idx[2] * 2 + ((k >> 2) & 1),
                    ],
                ));
            }
        } else {
            core.leaves.insert((level, idx), hex);
        }
    }
    // 2. 2:1 balance of the core cubes over faces, edges and corners
    let mut queue: Vec<(u32, [i64; 3])> = core
        .leaves
        .iter()
        .filter(|(_, h)| **h)
        .map(|(k, _)| *k)
        .collect();
    while let Some((level, idx)) = queue.pop() {
        if core.leaves.get(&(level, idx)) != Some(&true) {
            continue;
        }
        for d in NEIGHBOURS {
            loop {
                let n = [idx[0] + d[0], idx[1] + d[1], idx[2] + d[2]];
                let Found::Leaf(l, i, hex) = core.find(level, n) else {
                    break;
                };
                if l + 1 >= level {
                    break;
                }
                core.leaves.remove(&(l, i));
                for k in 0..8 {
                    let child = [
                        i[0] * 2 + (k & 1),
                        i[1] * 2 + ((k >> 1) & 1),
                        i[2] * 2 + ((k >> 2) & 1),
                    ];
                    let h = hex || qualifies(&core, l + 1, child).0;
                    core.leaves.insert((l + 1, child), h);
                    if h {
                        queue.push((l + 1, child));
                    }
                }
            }
        }
    }
    // 3. a manifold shell: drop cubes at non-manifold edges and vertices
    core.index_cubes();
    for _ in 0..64 {
        let polygons = core.shell_polygons();
        let drop = core.non_manifold(&polygons);
        if drop.is_empty() {
            break;
        }
        core.stats.removed_for_manifold += drop.len();
        for c in drop {
            let key = core.cubes[c];
            core.leaves.insert(key, false);
        }
        core.index_cubes();
    }
    if core.cubes.is_empty() {
        return Err("no cube fits the hex core: the domain is too small for the margin".into());
    }
    // 4. shell triangles: two per plain square, a fan from the centre of a
    // square with corners of smaller cubes on its edges
    let polygons = core.shell_polygons();
    if !core.non_manifold(&polygons).is_empty() {
        return Err("hex core shell is not manifold".into());
    }
    let mut point_ids: HashMap<[i64; 3], usize> = HashMap::default();
    let mut id = |core: &mut HexCore, p: [i64; 3]| {
        *point_ids.entry(p).or_insert_with(|| {
            core.shell_points.push(core.position(p));
            core.shell_points.len() - 1
        })
    };
    for (ring, _) in &polygons {
        let ids: Vec<usize> = ring.iter().map(|&p| id(&mut core, p)).collect();
        if ids.len() == 4 {
            core.shell_triangles.push([ids[0], ids[1], ids[2]]);
            core.shell_triangles.push([ids[0], ids[2], ids[3]]);
        } else {
            let mut lo = ring[0];
            let mut hi = ring[0];
            for p in ring {
                for a in 0..3 {
                    lo[a] = lo[a].min(p[a]);
                    hi[a] = hi[a].max(p[a]);
                }
            }
            let centre = id(&mut core, [0, 1, 2].map(|a| (lo[a] + hi[a]) / 2));
            for k in 0..ids.len() {
                core.shell_triangles
                    .push([centre, ids[k], ids[(k + 1) % ids.len()]]);
            }
        }
    }
    let mut sizes: Vec<(f64, usize)> = Vec::new();
    for l in 0..=depth {
        let n = core.cubes.iter().filter(|c| c.0 == l).count();
        if n > 0 {
            sizes.push((core.edge(l), n));
        }
    }
    core.stats.cubes = core.cubes.len();
    core.stats.sizes = sizes;
    core.stats.shell_triangles = core.shell_triangles.len();
    Ok(core)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A box [0, hi] as twelve outward triangles.
    pub(crate) fn box_surface(hi: Vec3) -> (Vec<Vec3>, Vec<[usize; 3]>) {
        let v: Vec<Vec3> = (0..8)
            .map(|k| {
                Vec3::new(
                    if k & 1 == 1 { hi.x } else { 0. },
                    if k & 2 == 2 { hi.y } else { 0. },
                    if k & 4 == 4 { hi.z } else { 0. },
                )
            })
            .collect();
        let quads = [
            [0, 2, 3, 1],
            [4, 5, 7, 6],
            [0, 1, 5, 4],
            [2, 6, 7, 3],
            [0, 4, 6, 2],
            [1, 3, 7, 5],
        ];
        let mut t = Vec::new();
        for q in quads {
            t.push([q[0], q[1], q[2]]);
            t.push([q[0], q[2], q[3]]);
        }
        (v, t)
    }

    #[test]
    fn box_core_is_a_closed_manifold_shell() {
        let hi = Vec3::new(8., 6., 5.);
        let (v, t) = box_surface(hi);
        let core = build(&v, &t, 0.5, 0.4, None, 0.1, HexCoreOptions::default()).unwrap();
        assert!(core.cubes.len() > 10, "{:?}", core.stats);
        // closed, consistently wound shell
        let mut edges: HashMap<(usize, usize), i32> = HashMap::default();
        for tri in &core.shell_triangles {
            for k in 0..3 {
                *edges.entry((tri[k], tri[(k + 1) % 3])).or_default() += 1;
            }
        }
        for (&(a, b), &n) in &edges {
            assert_eq!(n, 1);
            assert_eq!(edges.get(&(b, a)), Some(&1));
        }
        // the shell encloses the cubes' volume, normals into the core
        let vol: f64 = core
            .shell_triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| core.shell_points[i]);
                a.dot(b.cross(c)) / 6.
            })
            .sum();
        let cubes: f64 = (0..core.cubes.len())
            .map(|c| core.cube_edge(c).powi(3))
            .sum();
        assert!((vol + cubes).abs() < 1e-9 * cubes, "{vol} {cubes}");
        // every cube clear of the wall by layers plus margin
        for c in 0..core.cubes.len() {
            let p = core.cube_centre(c);
            let h = core.cube_edge(c) / 2.;
            let d = [p.x, hi.x - p.x, p.y, hi.y - p.y, p.z, hi.z - p.z]
                .into_iter()
                .fold(f64::INFINITY, f64::min);
            assert!(d - h >= 0.1 + core.cube_edge(c) - 1e-12);
            assert_eq!(core.locate(p), Some(c));
        }
    }
}
