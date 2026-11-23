//! Clipped Voronoi core of a nonconvex duct volume (M4, OBJ-20).
//!
//! The core domain is the wall-offset volume scaffold, a set of convex
//! tetrahedra. Following Yan, Wang, Lévy & Liu (2010/2013, clipped Voronoi
//! diagrams for mesh generation) every tetrahedron is split into its
//! intersections with the Voronoi cells of the sites: each piece is the
//! tetrahedron cut by bisector planes of its site's nearest sites, complete
//! once the next site lies beyond twice the piece's radius (the security
//! radius of Ray, Sokolov, Lefebvre & Lévy 2018, as in ALG-04). The sites
//! that reach a tetrahedron are found by a walk over bisector faces from the
//! site nearest its centroid, which covers every tetrahedron exactly once.
//! A cell is the union of its site's pieces, so cells follow the nonconvex
//! domain without any convexity assumption. Sites are moved to the
//! density-weighted centroids of their cells (Lloyd's iteration towards a
//! centroidal Voronoi tessellation, Du, Faber & Gunzburger 1999), so every
//! face between two core cells lies on the bisector of two sites that sit
//! close to the cell centroids: S·d > 0 by construction.
//!
//! Wall layers: a wall triangle covered by one core cell keeps its prism
//! stack. A triangle shared by several cells is cut into the cells' wall
//! polygons, and each polygon is extruded through the layers along the same
//! layer nodes (poly-prism pieces with planar caps and bilinear sides), so
//! every innermost layer face matches its core face exactly.
use super::convex::{Convex, Cut, Src, kernel_point};
use super::volume_dual::{Prepared, Scaffold, compact_points, map_scaffold, prepare};
use super::*;
use crate::neighbours::KdTree;
use std::collections::HashMap;

/// Options of the Voronoi core.
#[derive(Clone, Copy, Debug)]
pub struct VoronoiOptions {
    /// Lloyd iterations before the final diagram.
    pub lloyd_iterations: usize,
}

impl Default for VoronoiOptions {
    fn default() -> Self {
        Self {
            lloyd_iterations: DEFAULT_LLOYD_ITERATIONS,
        }
    }
}

pub const DEFAULT_LLOYD_ITERATIONS: usize = 10;

/// What the construction did, for the report.
#[derive(Clone, Debug, Default)]
pub struct VoronoiStats {
    pub sites: usize,
    pub core_cells: usize,
    /// Sites whose restricted cell was split into star-shaped parts.
    pub split_sites: usize,
    /// Sites whose restricted cell came out empty and were dropped.
    pub empty_sites: usize,
    pub pieces: usize,
    pub whole_wall_triangles: usize,
    pub split_wall_triangles: usize,
    pub layer_pieces: usize,
    pub lloyd_iterations: usize,
    /// Largest site move of the last Lloyd iteration over the local spacing.
    pub last_relative_shift: f64,
    pub merged_core_faces: usize,
    /// Cells split at reentrant wall edges.
    pub reentrant_splits: usize,
    /// Wall pieces handed to a neighbouring cell at reentrant edges.
    pub reassigned_pieces: usize,
    /// Sites split into a pair at reentrant wall edges.
    pub added_sites: usize,
    /// Wall pieces at a reentrant edge made cells of their own.
    pub detached_pieces: usize,
}

/// Rebuilds of the final diagram with sites added at reentrant edges.
const SITE_ROUNDS: usize = 3;
/// Distance of a split site pair from the halving plane of a reentrant
/// edge, relative to the size of the cell it replaces.
const REENTRANT_SITE_OFFSET: f64 = 0.3;
/// Volume share of its cell up to which a detached reentrant cell grows.
const DETACH_SHARE: f64 = 0.3;
/// A core cell whose centroid is closer than this (relative to its size)
/// to the plane of one of its wall faces, or behind it, is split.
const REENTRANT_CLEARANCE: f64 = 0.05;
/// Wall faces within 30° belong to one plane of a reentrant split.
const REENTRANT_PLANE_COS: f64 = 0.866_025_403_784_438_6;
const SECURITY_MARGIN: f64 = 1.0 + 1e-9;
const INITIAL_NEIGHBOURS: usize = 32;

struct Domain {
    tets: Vec<[Vec3; 4]>,
    /// Density at the four corners of each tetrahedron.
    rho: Vec<[f64; 4]>,
    tol: f64,
}

fn det(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    a.cross(b).dot(c)
}

fn bary(p: &[Vec3; 4], x: Vec3) -> [f64; 4] {
    let (e1, e2, e3, r) = (p[1] - p[0], p[2] - p[0], p[3] - p[0], x - p[0]);
    let v = det(e1, e2, e3);
    let b1 = det(r, e2, e3) / v;
    let b2 = det(e1, r, e3) / v;
    let b3 = det(e1, e2, r) / v;
    [1. - b1 - b2 - b3, b1, b2, b3]
}

struct Sites<'a> {
    points: &'a [Vec3],
    tree: KdTree,
    near: Vec<Vec<(f64, usize)>>,
}

impl<'a> Sites<'a> {
    fn new(points: &'a [Vec3]) -> Self {
        let tree = KdTree::new(points);
        let k = INITIAL_NEIGHBOURS.min(points.len().saturating_sub(1));
        let near = (0..points.len()).map(|i| tree.nearest(i, k)).collect();
        Self { points, tree, near }
    }

    /// The part of tetrahedron `t` in the Voronoi cell of site `s`.
    fn piece(&self, dom: &Domain, t: usize, s: usize) -> Result<Option<Convex>, String> {
        let mut c = Convex::tet(dom.tets[t]);
        let site = self.points[s];
        let radius = |c: &Convex| {
            c.pts
                .iter()
                .map(|&p| (p - site).length())
                .fold(0.0, f64::max)
        };
        let mut r = radius(&c);
        let others = self.points.len() - 1;
        let mut owned: Option<Vec<(f64, usize)>> = None;
        let mut i = 0;
        loop {
            let list = owned.as_deref().unwrap_or(&self.near[s]);
            if i == list.len() {
                if list.len() >= others {
                    break;
                }
                let k = (2 * list.len()).max(1).min(others);
                owned = Some(self.tree.nearest(s, k));
                continue;
            }
            let (d2, o) = list[i];
            i += 1;
            if d2.sqrt() > 2.0 * (r + dom.tol) * SECURITY_MARGIN {
                break;
            }
            let other = self.points[o];
            let n = other - site;
            let len = n.length();
            if len.is_nan() || len <= 0.0 {
                return Err(format!("sites {s} and {o} coincide"));
            }
            let n = n / len;
            match c
                .clip(n, n.dot((site + other) * 0.5), Src::Bis(o as u32), dom.tol)
                .map_err(|e| format!("site {s}, tetrahedron {t}: {e}"))?
            {
                Cut::Inside => {}
                Cut::Clipped => r = radius(&c),
                Cut::Empty => return Ok(None),
            }
        }
        Ok(Some(c))
    }

    /// Every non-empty piece of tetrahedron `t`, by site.
    fn split_tet(&self, dom: &Domain, t: usize) -> Result<Vec<(usize, Convex)>, String> {
        let p = &dom.tets[t];
        let centre = (p[0] + p[1] + p[2] + p[3]) * 0.25;
        let start = self.tree.closest(centre).ok_or("no sites")?;
        let mut seen = vec![start];
        let mut queue = vec![start];
        let mut out = Vec::new();
        while let Some(s) = queue.pop() {
            if let Some(c) = self.piece(dom, t, s)? {
                for (src, _) in &c.faces {
                    if let Src::Bis(o) = *src {
                        let o = o as usize;
                        if !seen.contains(&o) {
                            seen.push(o);
                            queue.push(o);
                        }
                    }
                }
                out.push((s, c));
            }
        }
        out.sort_by_key(|(s, _)| *s);
        Ok(out)
    }
}

/// Run `f` on every tetrahedron on all cores; results in tetrahedron order.
fn per_tet<T: Send>(
    count: usize,
    f: &(dyn Fn(usize) -> Result<T, String> + Sync),
) -> Result<Vec<T>, String> {
    let threads = std::thread::available_parallelism()
        .map_or(1, |n| n.get())
        .clamp(1, 64);
    let chunk = count.div_ceil(threads).max(1);
    let results: Vec<Result<Vec<T>, String>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..count)
            .step_by(chunk)
            .map(|start| scope.spawn(move || (start..(start + chunk).min(count)).map(f).collect()))
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().unwrap_or_else(|_| Err("worker panicked".into())))
            .collect()
    });
    let mut out = Vec::with_capacity(count);
    for r in results {
        out.extend(r?);
    }
    Ok(out)
}

/// Merge vertices closer than the tolerance into one id.
struct Weld {
    cell: f64,
    tol2: f64,
    grid: HashMap<(i64, i64, i64), Vec<usize>>,
}

impl Weld {
    fn key(&self, p: Vec3) -> (i64, i64, i64) {
        (
            (p.x / self.cell).floor() as i64,
            (p.y / self.cell).floor() as i64,
            (p.z / self.cell).floor() as i64,
        )
    }
    fn insert_new(&mut self, points: &mut Vec<Vec3>, p: Vec3) -> usize {
        let id = points.len();
        points.push(p);
        let k = self.key(p);
        self.grid.entry(k).or_default().push(id);
        id
    }
    fn register(&mut self, points: &[Vec3], id: usize) {
        let k = self.key(points[id]);
        self.grid.entry(k).or_default().push(id);
    }
    fn id(&mut self, points: &mut Vec<Vec3>, p: Vec3) -> usize {
        let (x, y, z) = self.key(p);
        let mut best: Option<(f64, usize)> = None;
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(ids) = self.grid.get(&(x + dx, y + dy, z + dz)) {
                        for &i in ids {
                            let d = (points[i] - p).dot(points[i] - p);
                            if d <= self.tol2 && best.is_none_or(|b| (d, i) < b) {
                                best = Some((d, i));
                            }
                        }
                    }
                }
            }
        }
        match best {
            Some((_, i)) => i,
            None => self.insert_new(points, p),
        }
    }
}

struct Piece {
    site: usize,
    tet: usize,
    faces: Vec<(Src, Vec<usize>)>,
    volume: f64,
    centroid: Vec3,
}

fn clean_ring(mut ring: Vec<usize>) -> Option<Vec<usize>> {
    ring.dedup();
    while ring.len() > 1 && ring.first() == ring.last() {
        ring.pop();
    }
    let mut seen = ring.clone();
    seen.sort_unstable();
    seen.dedup();
    (ring.len() >= 3 && seen.len() == ring.len()).then_some(ring)
}

/// Merge coplanar rings shared by the same two cells into one polygon when
/// they form one simple loop that is star-shaped from its vertex mean.
fn merge_rings(points: &[Vec3], rings: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    if rings.len() < 2 {
        return rings;
    }
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for ring in &rings {
        for k in 0..ring.len() {
            let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
            if let Some(i) = edges.iter().position(|&e| e == (b, a)) {
                edges.swap_remove(i);
            } else {
                edges.push((a, b));
            }
        }
    }
    let mut next = HashMap::new();
    for &(a, b) in &edges {
        if next.insert(a, b).is_some() {
            return rings;
        }
    }
    let Some(&start) = next.keys().min() else {
        return rings;
    };
    let mut ring = vec![start];
    let mut at = next[&start];
    while at != start {
        if ring.len() > next.len() {
            return rings;
        }
        ring.push(at);
        at = match next.get(&at) {
            Some(&n) => n,
            None => return rings,
        };
    }
    if ring.len() != next.len() {
        return rings;
    }
    let (n, c) = super::convex::ring_plane(points, &ring);
    let unit = n / n.length();
    let size = ring
        .iter()
        .map(|&v| (points[v] - c).length())
        .fold(0.0, f64::max);
    if ring
        .iter()
        .any(|&v| (points[v] - c).dot(unit).abs() > 1e-9 * size)
    {
        return rings;
    }
    for k in 0..ring.len() {
        let (a, b) = (points[ring[k]], points[ring[(k + 1) % ring.len()]]);
        if (a - c).cross(b - c).dot(n) <= 0.0 {
            return rings;
        }
    }
    vec![ring]
}

pub fn build_voronoi_core(
    scaffold: Scaffold,
    steps: &[f64],
    options: VoronoiOptions,
) -> Result<(DuctMesh, ClosedSurface, VoronoiStats), String> {
    let mut seeds = None;
    let mut added = 0;
    for round in 0..=SITE_ROUNDS {
        match build(
            &scaffold,
            steps,
            options,
            seeds.take(),
            round == SITE_ROUNDS,
        )? {
            Built::Done(done) => {
                let (mesh, surface, mut stats) = *done;
                stats.added_sites = added;
                return Ok((mesh, surface, stats));
            }
            Built::Retry {
                sites,
                extra,
                stats,
            } => {
                added += extra.len() / 2;
                let mut all = sites;
                all.extend(extra);
                seeds = Some((all, stats));
            }
        }
    }
    unreachable!("the last round never retries")
}

enum Built {
    Done(Box<(DuctMesh, ClosedSurface, VoronoiStats)>),
    /// Sites to add where a cell wraps a reentrant wall edge.
    Retry {
        sites: Vec<Vec3>,
        extra: Vec<Vec3>,
        stats: VoronoiStats,
    },
}

fn build(
    scaffold: &Scaffold,
    steps: &[f64],
    options: VoronoiOptions,
    seeds: Option<(Vec<Vec3>, VoronoiStats)>,
    last: bool,
) -> Result<Built, String> {
    let Prepared {
        surface,
        patches,
        mut front,
        surface_ids,
    } = prepare(scaffold, steps)?;
    let mapped = map_scaffold(scaffold, &surface_ids, &mut front)?;
    let layers = steps.len();
    let scale = surface.bounds().extent().max_abs_component();
    let mut stats = VoronoiStats::default();

    // Domain: the wall-offset scaffold tetrahedra and a density from the
    // scaffold spacing (CVT cell size h ∝ ρ^(-1/5) in 3D).
    let pos = |n: usize| front.points[mapped[n]];
    let mut spacing = vec![(0.0, 0usize); scaffold.points.len()];
    for t in &scaffold.tetrahedra {
        for i in 0..4 {
            for j in i + 1..4 {
                let l = (pos(t[i]) - pos(t[j])).length();
                for k in [t[i], t[j]] {
                    spacing[k].0 += l;
                    spacing[k].1 += 1;
                }
            }
        }
    }
    let h: Vec<f64> = spacing
        .iter()
        .map(|&(s, n)| if n > 0 { s / n as f64 } else { 0.0 })
        .collect();
    let h_min = h
        .iter()
        .copied()
        .filter(|&v| v > 0.0)
        .fold(f64::INFINITY, f64::min);
    let mut tets = Vec::with_capacity(scaffold.tetrahedra.len());
    let mut rho = Vec::with_capacity(scaffold.tetrahedra.len());
    for (index, t) in scaffold.tetrahedra.iter().enumerate() {
        let p = t.map(pos);
        let six = det(p[1] - p[0], p[2] - p[0], p[3] - p[0]);
        if !six.is_finite() || six <= 0. {
            return Err(format!(
                "inverted scaffold tetrahedron {index} after wall offset, determinant {six}"
            ));
        }
        tets.push(p);
        rho.push(t.map(|n| (h_min / h[n]).powi(5)));
    }
    let dom = Domain {
        tets,
        rho,
        tol: scale * 1e-12,
    };

    // Initial sites: the moved scaffold nodes that some tetrahedron uses.
    let mut used = vec![false; scaffold.points.len()];
    for t in &scaffold.tetrahedra {
        for &n in t {
            used[n] = true;
        }
    }
    let mut sites: Vec<Vec3> = (0..scaffold.points.len())
        .filter(|&n| used[n])
        .map(pos)
        .collect();
    let mut site_h: Vec<f64> = (0..scaffold.points.len())
        .filter(|&n| used[n])
        .map(|n| h[n])
        .collect();

    let lloyd = if let Some((fixed, previous)) = seeds {
        sites = fixed;
        stats = previous;
        0
    } else {
        options.lloyd_iterations
    };
    for _ in 0..lloyd {
        let s = Sites::new(&sites);
        let moments = per_tet(dom.tets.len(), &|t| {
            let p = &dom.tets[t];
            let r = dom.rho[t];
            let density = |x: Vec3| {
                let b = bary(p, x);
                b[0] * r[0] + b[1] * r[1] + b[2] * r[2] + b[3] * r[3]
            };
            Ok(s.split_tet(&dom, t)?
                .into_iter()
                .map(|(site, c)| {
                    let (_, mass, m) = c.moments(&density);
                    (site, mass, m)
                })
                .collect::<Vec<_>>())
        })?;
        let mut mass = vec![0.0; sites.len()];
        let mut moment = vec![Vec3::ZERO; sites.len()];
        for list in moments {
            for (site, m, x) in list {
                mass[site] += m;
                moment[site] = moment[site] + x;
            }
        }
        let mut shift = 0.0_f64;
        let mut next_sites = Vec::with_capacity(sites.len());
        let mut next_h = Vec::with_capacity(sites.len());
        for i in 0..sites.len() {
            if mass[i] > 0.0 {
                let c = moment[i] / mass[i];
                shift = shift.max((c - sites[i]).length() / site_h[i].max(f64::MIN_POSITIVE));
                next_sites.push(c);
                next_h.push(site_h[i]);
            } else {
                stats.empty_sites += 1;
            }
        }
        sites = next_sites;
        site_h = next_h;
        stats.lloyd_iterations += 1;
        stats.last_relative_shift = shift;
    }

    // Final diagram, welded into global vertex ids.
    let s = Sites::new(&sites);
    let mut weld = Weld {
        cell: scale * 1e-9,
        tol2: (scale * 1e-10).powi(2),
        grid: HashMap::new(),
    };
    for &m in &mapped {
        weld.register(&front.points, m);
    }
    let mut pieces: Vec<Piece> = Vec::new();
    let mut site_tet: HashMap<(usize, usize), usize> = HashMap::new();
    const CHUNK: usize = 20_000;
    for start in (0..dom.tets.len()).step_by(CHUNK) {
        let end = (start + CHUNK).min(dom.tets.len());
        let split = per_tet(end - start, &|i| s.split_tet(&dom, start + i))?;
        for (i, list) in split.into_iter().enumerate() {
            let t = start + i;
            for (site, c) in list {
                let ids: Vec<usize> = c
                    .pts
                    .iter()
                    .map(|&p| weld.id(&mut front.points, p))
                    .collect();
                let faces = c
                    .faces
                    .iter()
                    .filter_map(|(src, ring)| {
                        clean_ring(ring.iter().map(|&v| ids[v]).collect()).map(|r| (*src, r))
                    })
                    .collect();
                site_tet.insert((site, t), pieces.len());
                pieces.push(Piece {
                    site,
                    tet: t,
                    faces,
                    volume: c.volume(),
                    centroid: {
                        let (v, _, m) = c.moments(&|_| 1.0);
                        m / v
                    },
                });
            }
        }
    }
    stats.pieces = pieces.len();
    let points = &mut front.points;

    // Tetrahedron faces: neighbour across each face, or the boundary triangle.
    let sorted3 = |mut k: [usize; 3]| {
        k.sort_unstable();
        k
    };
    let local_face = |t: &[usize; 4], k: u8| {
        let k = k as usize;
        let mut f = [0; 3];
        let mut j = 0;
        for (i, &n) in t.iter().enumerate() {
            if i != k {
                f[j] = n;
                j += 1;
            }
        }
        sorted3(f)
    };
    let mut face_tets: HashMap<[usize; 3], Vec<usize>> = HashMap::new();
    for (i, t) in scaffold.tetrahedra.iter().enumerate() {
        for k in 0..4 {
            face_tets.entry(local_face(t, k)).or_default().push(i);
        }
    }
    let boundary_index: HashMap<[usize; 3], usize> = scaffold
        .boundary
        .iter()
        .enumerate()
        .map(|(i, &(t, _))| (sorted3(t), i))
        .collect();
    let across = |t: usize, k: u8| -> Option<usize> {
        face_tets[&local_face(&scaffold.tetrahedra[t], k)]
            .iter()
            .copied()
            .find(|&o| o != t)
    };

    // Parts: connected, star-shaped unions of a site's pieces.
    let mut by_site: Vec<Vec<usize>> = vec![Vec::new(); sites.len()];
    for (i, p) in pieces.iter().enumerate() {
        by_site[p.site].push(i);
    }
    let neighbours_of = |i: usize| -> Vec<usize> {
        let p = &pieces[i];
        p.faces
            .iter()
            .filter_map(|(src, _)| match *src {
                Src::Tet(k) => across(p.tet, k).and_then(|o| site_tet.get(&(p.site, o)).copied()),
                Src::Bis(o) => site_tet.get(&(o as usize, p.tet)).copied(),
                Src::Cut => None,
            })
            .collect()
    };
    let exposed = |group: &[usize]| -> Vec<Vec<usize>> {
        let mut out = Vec::new();
        for &i in group {
            let p = &pieces[i];
            for (src, ring) in &p.faces {
                if let Src::Tet(k) = *src
                    && let Some(o) = across(p.tet, k)
                    && let Some(q) = site_tet.get(&(p.site, o))
                    && group.contains(q)
                {
                    continue;
                }
                // a bisector towards a piece of the same part is interior
                if let Src::Bis(o) = *src
                    && let Some(q) = site_tet.get(&(o as usize, p.tet))
                    && group.contains(q)
                {
                    continue;
                }
                out.push(ring.clone());
            }
        }
        out
    };
    let mut part_of = vec![usize::MAX; pieces.len()];
    let mut parts: Vec<(Vec<usize>, Vec3)> = Vec::new();
    for list in &by_site {
        let mut left: Vec<usize> = list.clone();
        let before = parts.len();
        while let Some(&seed) = left.first() {
            // connected component
            let mut comp = vec![seed];
            let mut k = 0;
            while k < comp.len() {
                for n in neighbours_of(comp[k]) {
                    if left.contains(&n) && !comp.contains(&n) {
                        comp.push(n);
                    }
                }
                k += 1;
            }
            left.retain(|i| !comp.contains(i));
            if let Some(kernel) = kernel_point(points, &exposed(&comp)) {
                comp.sort_unstable();
                parts.push((comp, kernel));
                continue;
            }
            // Greedy: grow star-shaped parts from the largest pieces.
            let mut rest = comp;
            rest.sort_by(|&a, &b| {
                pieces[b]
                    .volume
                    .total_cmp(&pieces[a].volume)
                    .then(a.cmp(&b))
            });
            while !rest.is_empty() {
                let mut part = vec![rest.remove(0)];
                let mut kernel = kernel_point(points, &exposed(&part))
                    .ok_or("a convex Voronoi piece has no kernel")?;
                loop {
                    let mut grown = false;
                    for j in 0..rest.len() {
                        let q = rest[j];
                        if !part.iter().any(|&i| neighbours_of(i).contains(&q)) {
                            continue;
                        }
                        let mut trial = part.clone();
                        trial.push(q);
                        if let Some(kp) = kernel_point(points, &exposed(&trial)) {
                            part = trial;
                            kernel = kp;
                            rest.remove(j);
                            grown = true;
                            break;
                        }
                    }
                    if !grown {
                        break;
                    }
                }
                part.sort_unstable();
                parts.push((part, kernel));
            }
        }
        if parts.len() - before > 1 {
            stats.split_sites += 1;
        }
    }
    // A cell around a reentrant wall edge can have its centroid behind one
    // of its wall faces, so the face to the wall layer there would point
    // the wrong way (S·d ≤ 0). Split such a cell into one part per wall
    // plane, each piece joining the plane it lies farthest in front of (the
    // repair of the median dual, research/m4-reentrant-dual-split.md).
    let wall_rings = |group: &[usize]| -> Vec<(Vec3, Vec3, usize)> {
        let mut out = Vec::new();
        for &i in group {
            let p = &pieces[i];
            for (src, ring) in &p.faces {
                if let Src::Tet(k) = *src
                    && across(p.tet, k).is_none()
                    && let Some(&bi) =
                        boundary_index.get(&local_face(&scaffold.tetrahedra[p.tet], k))
                    && scaffold.boundary[bi].1 == 3
                {
                    let (n, c) = super::convex::ring_plane(points, ring);
                    out.push((n / n.length(), c, bi));
                }
            }
        }
        out
    };
    let centre = |group: &[usize]| {
        let v: f64 = group.iter().map(|&i| pieces[i].volume).sum();
        group.iter().fold(Vec3::ZERO, |s, &i| {
            s + pieces[i].centroid * pieces[i].volume
        }) / v
    };
    // centroid of the innermost layer prism over each wall triangle
    let wall_rank: HashMap<usize, usize> = (0..scaffold.boundary.len())
        .filter(|&bi| scaffold.boundary[bi].1 == 3)
        .enumerate()
        .map(|(r, bi)| (bi, r))
        .collect();
    let layer_centre = |bi: usize| -> Option<Vec3> {
        let r = *wall_rank.get(&bi)?;
        let n = (layers > 0).then(|| front.cells.get((layers - 1) * wall_rank.len() + r))??;
        Some(
            n.element
                .nodes()
                .iter()
                .fold(Vec3::ZERO, |s, &i| s + points[i])
                / 6.0,
        )
    };
    let behind_walls = |group: &[usize]| -> Vec<usize> {
        let c = centre(group);
        let size = group.iter().map(|&i| pieces[i].volume).sum::<f64>().cbrt();
        wall_rings(group)
            .iter()
            .filter(|(n, fc, bi)| {
                let outside = layer_centre(*bi).unwrap_or(*fc);
                n.dot(c - outside) > -REENTRANT_CLEARANCE * size
            })
            .map(|w| w.2)
            .collect()
    };
    let behind = |group: &[usize]| behind_walls(group).len();
    // First repair: hand the pieces with such a wall face to a neighbouring
    // cell whose centroid is in front of it (a small wall contact of a cell
    // that lies mostly beyond a reentrant edge).
    // Root repair: a cell that wraps a reentrant wall edge has its centroid
    // behind one of the two walls. Its site is replaced by two sites on
    // either side of the plane that halves the angle between the walls, so
    // the bisector between them runs along that plane, each new cell keeps
    // to one wall, and the diagram is rebuilt with all faces on bisectors.
    if !last {
        let mut drop = vec![false; sites.len()];
        let mut extra = Vec::new();
        for (group, _) in &parts {
            if behind(group) == 0 {
                continue;
            }
            let site = pieces[group[0]].site;
            if drop[site] {
                continue;
            }
            // the two wall planes with the largest contact
            let mut planes: Vec<(Vec3, Vec3, f64)> = Vec::new();
            for &i in group {
                for (src, ring) in &pieces[i].faces {
                    if let Src::Tet(k) = *src
                        && across(pieces[i].tet, k).is_none()
                        && let Some(&bi) =
                            boundary_index.get(&local_face(&scaffold.tetrahedra[pieces[i].tet], k))
                        && scaffold.boundary[bi].1 == 3
                    {
                        let (n, c) = super::convex::ring_plane(points, ring);
                        let area = n.length();
                        let unit = n / area;
                        match planes
                            .iter_mut()
                            .find(|p| p.0.dot(unit) > REENTRANT_PLANE_COS)
                        {
                            Some(p) => {
                                p.1 = p.1 + c * area;
                                p.2 += area;
                            }
                            None => planes.push((unit, c * area, area)),
                        }
                    }
                }
            }
            planes.sort_by(|x, y| y.2.total_cmp(&x.2));
            if planes.len() < 2 {
                continue;
            }
            let (na, pa) = (planes[0].0, planes[0].1 / planes[0].2);
            let (nb, pb) = (planes[1].0, planes[1].1 / planes[1].2);
            let g = na - nb;
            let len = g.length();
            if len.is_nan() || len <= 1e-6 {
                continue;
            }
            let m = g / len;
            let x = sites[site];
            let f = na.dot(x - pa) - nb.dot(x - pb);
            let on = x - m * (f / len);
            let size = group.iter().map(|&i| pieces[i].volume).sum::<f64>().cbrt();
            drop[site] = true;
            extra.push(on + m * (REENTRANT_SITE_OFFSET * size));
            extra.push(on - m * (REENTRANT_SITE_OFFSET * size));
        }
        if !extra.is_empty() {
            let sites = sites
                .iter()
                .zip(&drop)
                .filter(|(_, d)| !**d)
                .map(|(s, _)| *s)
                .collect();
            return Ok(Built::Retry {
                sites,
                extra,
                stats,
            });
        }
    }

    let mut owner_of = vec![usize::MAX; pieces.len()];
    for (id, (group, _)) in parts.iter().enumerate() {
        for &i in group {
            owner_of[i] = id;
        }
    }
    let wall_piece = |i: usize| !wall_rings(&[i]).is_empty();
    for id in 0..parts.len() {
        for _ in 0..4 {
            let before = behind(&parts[id].0);
            if before == 0 {
                break;
            }
            let offenders: Vec<usize> = parts[id]
                .0
                .iter()
                .copied()
                .filter(|&i| wall_piece(i))
                .collect();
            let mut done = false;
            'pieces: for i in offenders {
                let p = &pieces[i];
                let mut candidates: Vec<usize> = p
                    .faces
                    .iter()
                    .filter_map(|(src, _)| match *src {
                        Src::Bis(o) => site_tet.get(&(o as usize, p.tet)).map(|&q| owner_of[q]),
                        _ => None,
                    })
                    .collect();
                candidates.sort_unstable();
                candidates.dedup();
                for q in candidates {
                    if q == id {
                        continue;
                    }
                    let rest: Vec<usize> =
                        parts[id].0.iter().copied().filter(|&j| j != i).collect();
                    if rest.is_empty() {
                        continue;
                    }
                    let mut grown = parts[q].0.clone();
                    grown.push(i);
                    grown.sort_unstable();
                    if behind(&grown) > behind(&parts[q].0) || behind(&rest) >= before {
                        continue;
                    }
                    let (Some(kq), Some(kr)) = (
                        kernel_point(points, &exposed(&grown)),
                        kernel_point(points, &exposed(&rest)),
                    ) else {
                        continue;
                    };
                    // the rest must stay connected
                    let mut reach = vec![rest[0]];
                    let mut k = 0;
                    while k < reach.len() {
                        for n in neighbours_of(reach[k]) {
                            if rest.contains(&n) && !reach.contains(&n) {
                                reach.push(n);
                            }
                        }
                        k += 1;
                    }
                    if reach.len() != rest.len() {
                        continue;
                    }
                    parts[id] = (rest, kr);
                    parts[q] = (grown, kq);
                    owner_of[i] = q;
                    stats.reassigned_pieces += 1;
                    done = true;
                    break 'pieces;
                }
                // No neighbour can take it: the piece wedged into the
                // reentrant edge seeds a cell of its own, grown over
                // neighbouring pieces of the part while it stays
                // star-shaped and in front of its wall faces, so it does
                // not end up a sliver.
                let split_ok = |cluster: &[usize]| -> Option<(Vec<usize>, Vec3, Vec3)> {
                    let rest: Vec<usize> = parts[id]
                        .0
                        .iter()
                        .copied()
                        .filter(|j| !cluster.contains(j))
                        .collect();
                    if rest.is_empty() || behind(&rest) >= before || behind(cluster) > 0 {
                        return None;
                    }
                    let kc = kernel_point(points, &exposed(cluster))?;
                    let kr = kernel_point(points, &exposed(&rest))?;
                    let mut reach = vec![rest[0]];
                    let mut k = 0;
                    while k < reach.len() {
                        for n in neighbours_of(reach[k]) {
                            if rest.contains(&n) && !reach.contains(&n) {
                                reach.push(n);
                            }
                        }
                        k += 1;
                    }
                    (reach.len() == rest.len()).then_some((rest, kr, kc))
                };
                let Some(mut best) = split_ok(&[i]) else {
                    continue;
                };
                let mut cluster = vec![i];
                let total: f64 = parts[id].0.iter().map(|&j| pieces[j].volume).sum();
                loop {
                    let vol: f64 = cluster.iter().map(|&j| pieces[j].volume).sum();
                    if vol > DETACH_SHARE * total {
                        break;
                    }
                    let mut front: Vec<usize> = cluster
                        .iter()
                        .flat_map(|&j| neighbours_of(j))
                        .filter(|n| owner_of[*n] == id && !cluster.contains(n))
                        .collect();
                    front.sort_unstable();
                    front.dedup();
                    // largest neighbour first
                    front.sort_by(|&x, &y| pieces[y].volume.total_cmp(&pieces[x].volume));
                    let mut grown = false;
                    for n in front {
                        let mut next = cluster.clone();
                        next.push(n);
                        next.sort_unstable();
                        if let Some(r) = split_ok(&next) {
                            cluster = next;
                            best = r;
                            grown = true;
                            break;
                        }
                    }
                    if !grown {
                        break;
                    }
                }
                let (rest, kr, ki) = best;
                parts[id] = (rest, kr);
                for &j in &cluster {
                    owner_of[j] = parts.len();
                }
                parts.push((cluster, ki));
                stats.detached_pieces += 1;
                done = true;
                break 'pieces;
            }
            if !done {
                break;
            }
        }
    }
    let mut repaired = Vec::with_capacity(parts.len());
    for (group, kernel) in parts {
        let before = behind(&group);
        if before == 0 {
            repaired.push((group, kernel));
            continue;
        }
        let mut planes: Vec<Vec3> = Vec::new();
        for (n, _, _) in wall_rings(&group) {
            if !planes.iter().any(|m| m.dot(n) > REENTRANT_PLANE_COS) {
                planes.push(n);
            }
        }
        let c = centre(&group);
        let mut by_plane = vec![Vec::new(); planes.len()];
        for &i in &group {
            let r = pieces[i].centroid - c;
            let k = (0..planes.len())
                .min_by(|&a, &b| r.dot(planes[a]).total_cmp(&r.dot(planes[b])))
                .unwrap_or(0);
            by_plane[k].push(i);
        }
        let mut split = Vec::new();
        let mut ok = planes.len() > 1;
        for list in by_plane {
            let mut left = list;
            while ok && !left.is_empty() {
                let mut comp = vec![left[0]];
                let mut k = 0;
                while k < comp.len() {
                    for n in neighbours_of(comp[k]) {
                        if left.contains(&n) && !comp.contains(&n) {
                            comp.push(n);
                        }
                    }
                    k += 1;
                }
                left.retain(|i| !comp.contains(i));
                comp.sort_unstable();
                match kernel_point(points, &exposed(&comp)) {
                    Some(kp) => split.push((comp, kp)),
                    None => ok = false,
                }
            }
        }
        let after: usize = split.iter().map(|(g, _)| behind(g)).sum();
        if ok && split.len() > 1 && after < before {
            stats.reentrant_splits += 1;
            repaired.extend(split);
        } else {
            repaired.push((group, kernel));
        }
    }
    let parts = repaired;
    for (id, (group, _)) in parts.iter().enumerate() {
        for &i in group {
            part_of[i] = id;
        }
    }
    stats.sites = sites.len();
    stats.core_cells = parts.len();

    // Faces between parts, and the parts' boundary polygons.
    let mut pair_faces: BTreeMap<(usize, usize), Vec<Vec<usize>>> = BTreeMap::new();
    let mut records: BTreeMap<usize, Vec<(usize, Vec<usize>)>> = BTreeMap::new();
    for (i, p) in pieces.iter().enumerate() {
        let a = part_of[i];
        for (src, ring) in &p.faces {
            let other = match *src {
                Src::Bis(o) => site_tet.get(&(o as usize, p.tet)).copied().ok_or_else(|| {
                    format!(
                        "bisector face of site {} has no partner in tetrahedron {}",
                        p.site, p.tet
                    )
                })?,
                Src::Tet(k) => match across(p.tet, k) {
                    None => {
                        let key = local_face(&scaffold.tetrahedra[p.tet], k);
                        let bi = *boundary_index
                            .get(&key)
                            .ok_or("open scaffold face that is no boundary triangle")?;
                        records.entry(bi).or_default().push((a, ring.clone()));
                        continue;
                    }
                    Some(o) => site_tet.get(&(p.site, o)).copied().ok_or_else(|| {
                        format!(
                            "piece of site {} in tetrahedron {} has no continuation in {o}",
                            p.site, p.tet
                        )
                    })?,
                },
                Src::Cut => return Err("unexpected cut face".into()),
            };
            let b = part_of[other];
            if a < b {
                pair_faces.entry((a, b)).or_default().push(ring.clone());
            }
        }
    }
    let mut core_faces: Vec<Vec<Vec<usize>>> = vec![Vec::new(); parts.len()];
    for ((a, b), rings) in pair_faces {
        let before = rings.len();
        let merged = merge_rings(points, rings);
        stats.merged_core_faces += before - merged.len();
        for ring in merged {
            core_faces[b].push(ring.iter().rev().copied().collect());
            core_faces[a].push(ring);
        }
    }

    // Layers.
    let facets = surface.facets();
    let wall: Vec<usize> = (0..facets.len()).filter(|&i| patches[i] == 3).collect();
    let coarse = std::mem::take(&mut front.cells);
    let old_boundary = std::mem::take(&mut front.boundary);
    let mut stack = vec![vec![usize::MAX; surface.vertices().len()]; layers + 1];
    for (v, s0) in stack[0].iter_mut().enumerate() {
        *s0 = v;
    }
    let mut thickness = vec![[f64::INFINITY, 0.0_f64]; layers];
    let mut minimum_j = f64::INFINITY;
    for level in 0..layers {
        for (r, &bi) in wall.iter().enumerate() {
            let [a, b, c] = facets[bi].triangle.vertices;
            let Element::Wedge(n) = coarse[level * wall.len() + r].element else {
                return Err("expected a layer prism".into());
            };
            stack[level][a] = n[0];
            stack[level][c] = n[1];
            stack[level][b] = n[2];
            stack[level + 1][a] = n[3];
            stack[level + 1][c] = n[4];
            stack[level + 1][b] = n[5];
            let (jac, _, range) = verify::wedge(n.map(|i| points[i]));
            minimum_j = minimum_j.min(jac);
            thickness[level][0] = thickness[level][0].min(range[0]);
            thickness[level][1] = thickness[level][1].max(range[1]);
        }
    }
    let mut rim: HashMap<(usize, usize), usize> = HashMap::new();
    for (i, f) in facets.iter().enumerate() {
        if patches[i] == 3 {
            continue;
        }
        let [a, b, c] = f.triangle.vertices;
        for (u, v) in [(a, b), (b, c), (c, a)] {
            rim.insert((u.min(v), u.max(v)), patches[i]);
        }
    }
    let mut cells: Vec<Cell> = Vec::new();
    let mut boundary: Vec<BoundaryFace> = Vec::new();
    let mut core_wall: Vec<(usize, Vec<usize>)> = Vec::new();
    if layers > 0 {
        // Barycentric coordinates of the core wall vertices: the vertices on
        // a wall triangle edge, in order along it, are shared by the two
        // innermost prisms beside it.
        let mut edge_points: HashMap<(usize, usize), Vec<(f64, usize)>> = HashMap::new();
        let mut whole = vec![false; facets.len()];
        for &bi in &wall {
            let [a, b, c] = facets[bi].triangle.vertices;
            let corners = [a, b, c];
            let inner = corners.map(|v| stack[layers][v]);
            let tri = inner.map(|n| points[n]);
            let normal = (tri[1] - tri[0]).cross(tri[2] - tri[0]);
            let n2 = normal.dot(normal);
            let rec = records.get(&bi).map_or(&[][..], |v| v.as_slice());
            whole[bi] =
                rec.len() == 1 && rec[0].1.len() == 3 && rec[0].1.iter().all(|g| inner.contains(g));
            if whole[bi] {
                stats.whole_wall_triangles += 1;
                continue;
            }
            stats.split_wall_triangles += 1;
            for (_, ring) in rec {
                for &g in ring {
                    if inner.contains(&g) {
                        continue;
                    }
                    let x = points[g];
                    let la = (tri[1] - x).cross(tri[2] - x).dot(normal) / n2;
                    let lb = (tri[2] - x).cross(tri[0] - x).dot(normal) / n2;
                    let l = [la, lb, 1.0 - la - lb];
                    let Some(j) = (0..3).find(|&j| l[j].abs() < 1e-9) else {
                        continue;
                    };
                    let (u, v) = (corners[(j + 1) % 3], corners[(j + 2) % 3]);
                    let (lu, lv) = (l[(j + 1) % 3], l[(j + 2) % 3]);
                    let (key, t) = if u < v {
                        ((u, v), lv / (lu + lv))
                    } else {
                        ((v, u), lu / (lu + lv))
                    };
                    let list = edge_points.entry(key).or_default();
                    if !list.iter().any(|&(_, h)| h == g) {
                        list.push((t, g));
                    }
                }
            }
        }
        for list in edge_points.values_mut() {
            list.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        }
        // The innermost prism of a triangle shared by several core cells
        // keeps its shape; its inner end is split into the cells' wall
        // polygons, and its sides gain the vertices on their inner edges.
        let mut new_index = vec![usize::MAX; coarse.len()];
        let mut split_cell = vec![false; coarse.len()];
        for level in 0..layers {
            for (r, &bi) in wall.iter().enumerate() {
                let old = level * wall.len() + r;
                let Element::Wedge(n) = coarse[old].element else {
                    return Err("expected a layer prism".into());
                };
                new_index[old] = cells.len();
                if level + 1 < layers || whole[bi] {
                    cells.push(Cell {
                        element: Element::Wedge(n),
                        layer: level + 1,
                        region: 0,
                    });
                    continue;
                }
                split_cell[old] = true;
                let [a, b, c] = facets[bi].triangle.vertices;
                // n = [top a, top c, top b, bottom a, bottom c, bottom b]
                let surface_of = |i: usize| [a, c, b][i % 3];
                let mut faces = vec![vec![n[0], n[2], n[1]]];
                for [p0, p1, q1, q0] in [[0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]] {
                    // bottom edge runs from q1 to q0
                    let (x, y) = (surface_of(q1), surface_of(q0));
                    let mut ring = vec![n[p0], n[p1], n[q1]];
                    if let Some(list) = edge_points.get(&(x.min(y), x.max(y))) {
                        if x < y {
                            ring.extend(list.iter().map(|&(_, g)| g));
                        } else {
                            ring.extend(list.iter().rev().map(|&(_, g)| g));
                        }
                    }
                    ring.push(n[q0]);
                    if let Some(&patch) = rim.get(&(x.min(y), x.max(y))) {
                        boundary.push(BoundaryFace {
                            nodes: ring.clone(),
                            patch,
                            owner: cells.len(),
                        });
                    }
                    faces.push(ring);
                }
                for (_, ring) in records.get(&bi).into_iter().flatten() {
                    faces.push(ring.iter().rev().copied().collect());
                }
                let interior = n.iter().fold(Vec3::ZERO, |s, &i| s + points[i]) / 6.0;
                let expected = verify::wedge(n.map(|i| points[i])).1;
                let mut poly = Polyhedron::new(faces, interior, None, expected);
                poly.prism = Some(vec![n]);
                cells.push(Cell {
                    element: Element::Polyhedron(Box::new(poly)),
                    layer: level + 1,
                    region: 0,
                });
                stats.layer_pieces += 1;
            }
        }
        for f in old_boundary {
            // the port strips of split prisms were rebuilt above
            if split_cell[f.owner] && f.nodes.len() == 4 {
                continue;
            }
            boundary.push(BoundaryFace {
                nodes: f.nodes,
                patch: f.patch,
                owner: new_index[f.owner],
            });
        }
        for &bi in &wall {
            for (part, ring) in records.get(&bi).into_iter().flatten() {
                core_wall.push((*part, ring.clone()));
            }
        }
    } else {
        for &bi in &wall {
            for (part, ring) in records.get(&bi).into_iter().flatten() {
                core_wall.push((*part, ring.clone()));
            }
        }
    }
    let offset = cells.len();
    let interfaces = core_wall.len();
    for (part, ring) in core_wall {
        if layers == 0 {
            boundary.push(BoundaryFace {
                nodes: ring.clone(),
                patch: 3,
                owner: offset + part,
            });
        }
        core_faces[part].push(ring);
    }
    // port faces of each core cell
    let mut ports: BTreeMap<(usize, usize), Vec<Vec<usize>>> = BTreeMap::new();
    for (bi, rec) in &records {
        if patches[*bi] != 3 {
            for (part, ring) in rec {
                ports
                    .entry((*part, patches[*bi]))
                    .or_default()
                    .push(ring.clone());
            }
        }
    }
    for ((part, patch), rings) in ports {
        for ring in rings {
            boundary.push(BoundaryFace {
                nodes: ring.clone(),
                patch,
                owner: offset + part,
            });
            core_faces[part].push(ring);
        }
    }
    for (id, ((group, kernel), faces)) in parts.iter().zip(core_faces).enumerate() {
        let volume: f64 = group.iter().map(|&i| pieces[i].volume).sum();
        if faces.len() < 4 {
            return Err(format!("core cell {id} has {} faces", faces.len()));
        }
        cells.push(Cell {
            element: Element::Polyhedron(Box::new(Polyhedron::new(faces, *kernel, None, volume))),
            layer: 0,
            region: 0,
        });
    }
    compact_points(points, &mut cells, &mut boundary);
    let mut measurements = verify::check_with_interface_count(
        &surface, &patches, points, &cells, &boundary, layers, interfaces,
    )?;
    if layers > 0 {
        measurements.normal_thickness_ranges = thickness;
        measurements.minimum_wedge_jacobian = Some(minimum_j);
    }
    Ok(Built::Done(Box::new((
        DuctMesh {
            points: front.points,
            cells,
            boundary,
            measurements,
        },
        surface,
        stats,
    ))))
}
