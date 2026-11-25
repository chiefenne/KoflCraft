//! Zero-thickness shells inside the fluid (SHL-01 to SHL-03): a planar,
//! convex shell polygon (a baffle, a plate, a fin) is imprinted into a
//! finished duct mesh by cutting the cells it crosses.
//!
//! Method (the maintainer's choice of 2026-11-25: cut the finished cells,
//! in the plane-based style of exact mesh Booleans, Trettner,
//! Nehring-Wirxel and Kobbelt 2022):
//! 1. Vertices within a small share of their shortest edge from the plane
//!    are snapped onto it, so no cut leaves a sliver; port-cap vertices are
//!    never moved. Every other vertex has a clear side.
//! 2. A cell is cut when it has vertices on both sides and its cross-section
//!    meets the shell polygon. Each edge it has across the plane gets one
//!    intersection vertex, inserted into every face that has the edge, so
//!    neighbours stay conforming.
//! 3. Each face of a cut cell is split at its two on-plane vertices into an
//!    upper and a lower part (a face crossed more than once is refused); the
//!    on-plane segments form the cross-section loop.
//! 4. The loop is split along the shell outline into the part inside the
//!    shell (a wall) and the parts outside it (interior faces between the
//!    two halves, so a free edge inside the fluid leaves closed cells).
//! 5. On-plane vertices strictly inside the shell are doubled, one copy per
//!    side, so the two faces on the shell are distinct boundary faces: the
//!    side the normal points into is `<name>`, the other `<name>_back`.

use std::collections::HashMap;

use mesh_core::geometry::Vec3;

use super::polyhedral::Polyhedron;
use super::{BoundaryFace, Cell, DuctMesh, Element};

/// A planar shell: a convex polygon, counter-clockwise about its normal.
#[derive(Clone, Debug, PartialEq)]
pub struct PlanarShell {
    pub name: String,
    pub polygon: Vec<Vec3>,
}

/// What an imprint did.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ShellStats {
    pub cut_cells: usize,
    pub snapped_vertices: usize,
    pub doubled_vertices: usize,
    pub front_faces: usize,
    pub back_faces: usize,
    pub front_area: f64,
    pub back_area: f64,
    /// Area of the shell polygon.
    pub polygon_area: f64,
    /// Patch ids of the front and back sides.
    pub patches: [usize; 2],
    /// Cells split across the shell first, where it crossed them as a strip.
    pub through_cuts: usize,
}

/// Face rings on the upper and on the lower side (or inside and outside).
type Parts = (Vec<Vec<usize>>, Vec<Vec<usize>>);

/// Share of the shortest incident edge within which a vertex snaps. Off:
/// snapping moves vertices of faces that cross the plane, which makes faces
/// run along it where they cross; the exact cut takes the small pieces.
const SNAP: f64 = 0.0;

struct Frame {
    origin: Vec3,
    n: Vec3,
    u: Vec3,
    v: Vec3,
    /// The polygon in plane coordinates.
    poly: Vec<[f64; 2]>,
    scale: f64,
}

impl Frame {
    fn new(polygon: &[Vec3]) -> Result<Self, String> {
        if polygon.len() < 3 || polygon.iter().any(|p| !p.is_finite()) {
            return Err("a shell polygon needs at least three finite corners".into());
        }
        let origin = polygon[0];
        let mut normal = Vec3::ZERO;
        for i in 0..polygon.len() {
            let (a, b) = (
                polygon[i] - origin,
                polygon[(i + 1) % polygon.len()] - origin,
            );
            normal = normal + a.cross(b);
        }
        let area2 = normal.length();
        let scale = polygon
            .iter()
            .map(|p| (*p - origin).length())
            .fold(0.0, f64::max);
        if area2.is_nan() || area2 <= 1e-12 * scale * scale {
            return Err("the shell polygon has no area".into());
        }
        let n = normal / area2;
        let u = {
            let e = polygon[1] - polygon[0];
            e / e.length()
        };
        let v = n.cross(u);
        let poly: Vec<[f64; 2]> = polygon
            .iter()
            .map(|p| {
                let d = *p - origin;
                [d.dot(u), d.dot(v)]
            })
            .collect();
        for p in polygon {
            if ((*p - origin).dot(n)).abs() > 1e-9 * scale {
                return Err("the shell polygon is not planar".into());
            }
        }
        let k = poly.len();
        for i in 0..k {
            let (a, b, c) = (poly[i], poly[(i + 1) % k], poly[(i + 2) % k]);
            if cross2(sub2(b, a), sub2(c, b)) <= 0.0 {
                return Err("the shell polygon must be convex, its corners counter-clockwise about its normal, without repeated or collinear corners".into());
            }
        }
        Ok(Self {
            origin,
            n,
            u,
            v,
            poly,
            scale,
        })
    }

    fn distance(&self, p: Vec3) -> f64 {
        (p - self.origin).dot(self.n)
    }

    fn uv(&self, p: Vec3) -> [f64; 2] {
        let d = p - self.origin;
        [d.dot(self.u), d.dot(self.v)]
    }

    fn point(&self, q: [f64; 2]) -> Vec3 {
        self.origin + self.u * q[0] + self.v * q[1]
    }

    /// Signed distance of `q` inside edge `i` (positive inside).
    fn edge_side(&self, i: usize, q: [f64; 2]) -> f64 {
        let k = self.poly.len();
        let (a, b) = (self.poly[i], self.poly[(i + 1) % k]);
        let e = sub2(b, a);
        cross2(e, sub2(q, a)) / (e[0].hypot(e[1]))
    }

    /// Smallest edge distance: > 0 strictly inside, < 0 outside.
    fn inside(&self, q: [f64; 2]) -> f64 {
        (0..self.poly.len())
            .map(|i| self.edge_side(i, q))
            .fold(f64::INFINITY, f64::min)
    }

    fn area(&self) -> f64 {
        polygon_area2(&self.poly)
    }
}

fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn polygon_area2(p: &[[f64; 2]]) -> f64 {
    let mut s = 0.0;
    for i in 0..p.len() {
        s += cross2(p[i], p[(i + 1) % p.len()]);
    }
    0.5 * s
}

fn faces_of(cell: &Cell) -> Vec<Vec<usize>> {
    cell.element.faces()
}

fn ring_edges(r: &[usize]) -> impl Iterator<Item = (usize, usize)> + '_ {
    (0..r.len()).map(move |i| (r[i], r[(i + 1) % r.len()]))
}

fn key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

/// Insert the vertex of every split edge into a ring.
fn insert_splits(r: &[usize], splits: &HashMap<(usize, usize), Vec<usize>>) -> Vec<usize> {
    let mut out = Vec::with_capacity(r.len() + 2);
    for (a, b) in ring_edges(r) {
        out.push(a);
        if let Some(mid) = splits.get(&key(a, b)) {
            // the vertices are stored from the lower to the higher index
            if a < b {
                out.extend(mid.iter().copied());
            } else {
                out.extend(mid.iter().rev().copied());
            }
        }
    }
    out
}

/// A ring split by the plane into its parts on the upper and the lower side,
/// or put whole on one side.
enum Split {
    Upper(Vec<usize>),
    Lower(Vec<usize>),
    /// The ring lies in the plane.
    OnPlane,
    Both {
        upper: Vec<Vec<usize>>,
        lower: Vec<Vec<usize>>,
    },
}

/// Split a ring at its on-plane crossing vertices. A face crossed more than
/// once (not convex, or not flat) is split along chords that pair the
/// crossings in their order along the line where the plane meets the face,
/// as a simple polygon is split by a line.
fn split_ring(r: &[usize], side: &[i8], points: &[Vec3], n: Vec3) -> Result<Split, String> {
    let has = |s: i8| r.iter().any(|&v| side[v] == s);
    match (has(1), has(-1)) {
        (false, false) => return Ok(Split::OnPlane),
        (true, false) => return Ok(Split::Upper(r.to_vec())),
        (false, true) => return Ok(Split::Lower(r.to_vec())),
        _ => {}
    }
    let m = r.len();
    // a run of on-plane vertices where the side changes goes with the upper
    // part, but for its vertex next to the lower side, where the parts meet
    let mut vs: Vec<i8> = r.iter().map(|&v| side[v]).collect();
    if let Some(start) = (0..m).find(|&i| vs[i] != 0) {
        let mut i = 1;
        while i <= m {
            let at = (start + i) % m;
            if vs[at] != 0 {
                i += 1;
                continue;
            }
            let mut len = 0;
            while vs[(start + i + len) % m] == 0 {
                len += 1;
            }
            let before = vs[(start + i + m - 1) % m];
            let after = vs[(start + i + len) % m];
            if len >= 2 && before != after {
                for k in 0..len {
                    let pos = (start + i + k) % m;
                    let next_to_lower = (before == -1 && k == 0) || (after == -1 && k + 1 == len);
                    if !next_to_lower {
                        vs[pos] = 1;
                    }
                }
            }
            i += len;
        }
    }
    let s = |i: usize| vs[i % m];
    // crossings: on-plane vertices between an upper and a lower neighbour
    let cross: Vec<bool> = (0..m)
        .map(|i| {
            let (prev, next) = (s(i + m - 1), s(i + 1));
            s(i) == 0 && prev != 0 && next != 0 && prev != next
        })
        .collect();
    let ids: Vec<usize> = (0..m).filter(|&i| cross[i]).collect();
    if ids.len() < 2 || !ids.len().is_multiple_of(2) {
        return Err("a face crosses the shell plane between vertices".into());
    }
    let mut partner = vec![usize::MAX; m];
    if ids.len() == 2 {
        partner[ids[0]] = ids[1];
        partner[ids[1]] = ids[0];
    } else {
        // the line where the plane meets the face
        let mut normal = Vec3::ZERO;
        for (a, b) in ring_edges(r) {
            normal = normal + points[a].cross(points[b]);
        }
        let line = n.cross(normal);
        let mut order = ids.clone();
        order.sort_by(|&x, &y| points[r[x]].dot(line).total_cmp(&points[r[y]].dot(line)));
        for pair in order.chunks(2) {
            partner[pair[0]] = pair[1];
            partner[pair[1]] = pair[0];
        }
    }
    let mut parts = [Vec::new(), Vec::new()];
    for (k, target) in [(0, 1i8), (1, -1i8)] {
        let mut used = vec![false; m];
        for &i in &ids {
            if used[i] || s(i + 1) != target {
                continue;
            }
            let mut part = Vec::new();
            let mut j = i;
            loop {
                if s(j + 1) != target || part.len() > 2 * m {
                    return Err("the shell plane splits a face in a way that does not close; refine the mesh near the shell".into());
                }
                part.push(r[j]);
                used[j] = true;
                let mut q = (j + 1) % m;
                while !cross[q] {
                    part.push(r[q]);
                    q = (q + 1) % m;
                }
                part.push(r[q]);
                used[q] = true;
                let p = partner[q];
                if p == i {
                    break;
                }
                j = p;
            }
            parts[k].push(part);
        }
    }
    let [upper, lower] = parts;
    Ok(Split::Both { upper, lower })
}

/// Imprint `shell` into `mesh` and check the result: returns what was done.
/// Patches `<name>` (the side the normal points into) and `<name>_back` are
/// added as walls.
pub fn imprint_planar_shell(
    mesh: &mut DuctMesh,
    shell: &PlanarShell,
) -> Result<ShellStats, String> {
    let back_name = format!("{}_back", shell.name);
    if shell.name.is_empty() {
        return Err("a shell needs a name".into());
    }
    if mesh
        .patches
        .iter()
        .any(|(n, _)| *n == shell.name || *n == back_name)
    {
        return Err(format!("shell {}: a patch of that name exists", shell.name));
    }
    let f = Frame::new(&shell.polygon).map_err(|e| format!("shell {}: {e}", shell.name))?;
    let before = super::verify::check_edited(
        &mesh.points,
        &mesh.cells,
        &mesh.boundary,
        mesh.patches.len(),
    )?;
    let mut through_cuts = 0;
    for _ in 0..24 {
        match cut_by_plane(mesh, &f, Mode::Shell, &shell.name, back_name.clone())? {
            Outcome::Done(mut stats) => {
                check_shell(mesh, &f, &shell.name, &before, &stats)?;
                stats.through_cuts = through_cuts;
                return Ok(stats);
            }
            Outcome::Narrow(cells) => {
                // the shell crosses these cells as a strip: split each first
                // along the plane through one crossing outline edge, across
                // the shell, so every piece sees one outline crossing
                let mut by_edge: HashMap<usize, Vec<usize>> = HashMap::new();
                for (c, j) in cells {
                    by_edge.entry(j).or_default().push(c);
                }
                let mut edges: Vec<usize> = by_edge.keys().copied().collect();
                edges.sort_unstable();
                // cell ids shift with each cut: cut one edge per round
                let j = edges[0];
                let k = f.poly.len();
                let (a, b) = (shell.polygon[j], shell.polygon[(j + 1) % k]);
                let e = (b - a) / (b - a).length();
                let n = f.n.cross(e);
                let aux = Frame {
                    origin: a,
                    n,
                    u: e,
                    v: n.cross(e),
                    poly: Vec::new(),
                    scale: f.scale,
                };
                let mut cut_cells = by_edge.remove(&j).unwrap_or_default();
                cut_cells.sort_unstable();
                through_cuts += cut_cells.len();
                cut_by_plane(mesh, &aux, Mode::Through(cut_cells), "", String::new())?;
            }
        }
    }
    Err(format!(
        "shell {}: the shell is narrower than the cells it crosses; refine the mesh near it",
        shell.name
    ))
}

enum Mode {
    /// Imprint the shell polygon of the frame.
    Shell,
    /// Split these cells across the whole plane, with interior faces.
    Through(Vec<usize>),
}

enum Outcome {
    Done(ShellStats),
    /// Cells the shell crosses as a strip, each with a crossing outline edge.
    Narrow(Vec<(usize, usize)>),
}

/// The checks after an imprint: a conforming mesh that passes the solver
/// rules, the volume and the old patch areas kept, both sides equal.
fn check_shell(
    mesh: &mut DuctMesh,
    f: &Frame,
    name: &str,
    before: &super::verify::EditCheck,
    stats: &ShellStats,
) -> Result<(), String> {
    let after = super::verify::check_edited(
        &mesh.points,
        &mesh.cells,
        &mesh.boundary,
        mesh.patches.len(),
    )
    .map_err(|e| format!("shell {name}: {e}"))?;
    let close =
        |a: f64, b: f64| (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(f64::MIN_POSITIVE);
    if !close(after.volume, before.volume) {
        return Err(format!(
            "shell {name}: mesh volume {} differs from {} before the cut",
            after.volume, before.volume
        ));
    }
    for (i, (&a, &b)) in after.areas.iter().zip(&before.areas).enumerate() {
        if !close(a, b) && (a - b).abs() > 1e-12 * f.scale * f.scale {
            return Err(format!(
                "shell {name}: patch {} area {a} differs from {b} before the cut",
                i + 1
            ));
        }
    }
    if stats.front_faces == 0 {
        return Err(format!("shell {name} does not lie inside the mesh"));
    }
    if !close(stats.front_area, stats.back_area) {
        return Err(format!(
            "shell {name}: front area {} and back area {} differ",
            stats.front_area, stats.back_area
        ));
    }
    mesh.measurements.volume = after.volume;
    mesh.measurements.minimum_volume = after.minimum_volume;
    mesh.measurements.areas = after.areas;
    mesh.measurements.solver_rules = after.solver_rules;
    Ok(())
}

/// Cut the cells of `mesh` by the plane of `f` (see the module notes).
fn cut_by_plane(
    mesh: &mut DuctMesh,
    f: &Frame,
    mode: Mode,
    name: &str,
    back_name: String,
) -> Result<Outcome, String> {
    let shell_mode = matches!(mode, Mode::Shell);
    let mut stats = ShellStats {
        polygon_area: if shell_mode { f.area() } else { 0.0 },
        ..Default::default()
    };
    let cell_rings: Vec<Vec<Vec<usize>>> = mesh.cells.iter().map(faces_of).collect();
    // 1. sides; interior vertices close to the plane snap onto it
    let np = mesh.points.len();
    let mut moved: Vec<(usize, Vec3)> = Vec::new();
    let mut shortest = vec![f64::INFINITY; np];
    for rings in &cell_rings {
        for r in rings {
            for (a, b) in ring_edges(r) {
                let l = (mesh.points[a] - mesh.points[b]).length();
                shortest[a] = shortest[a].min(l);
                shortest[b] = shortest[b].min(l);
            }
        }
    }
    let mut on_boundary = vec![false; np];
    for b in &mesh.boundary {
        for &v in &b.nodes {
            on_boundary[v] = true;
        }
    }
    let mut side = vec![0i8; np];
    for v in 0..np {
        let d = f.distance(mesh.points[v]);
        let tol = SNAP * shortest[v];
        if d.abs() <= 1e-12 * f.scale {
            side[v] = 0;
        } else if shell_mode
            && d.abs() <= tol
            && !on_boundary[v]
            && f.inside(f.uv(mesh.points[v])) > -tol
        {
            moved.push((v, mesh.points[v]));
            mesh.points[v] = mesh.points[v] - f.n * d;
            stats.snapped_vertices += 1;
            side[v] = 0;
        } else {
            side[v] = if d > 0.0 { 1 } else { -1 };
        }
    }
    // 2. the cells to cut: both sides, and the cross-section meets the shell
    let crosses = |a: usize, b: usize, side: &[i8]| side[a] * side[b] == -1;
    let section_point = |a: usize, b: usize, p: &[Vec3]| {
        let (da, db) = (f.distance(p[a]), f.distance(p[b]));
        p[a] + (p[b] - p[a]) * (da / (da - db))
    };
    let tol = 1e-9 * f.scale;
    let mut cut = vec![false; mesh.cells.len()];
    for (c, rings) in cell_rings.iter().enumerate() {
        let mut up = false;
        let mut down = false;
        let mut section: Vec<[f64; 2]> = Vec::new();
        for r in rings {
            for (a, b) in ring_edges(r) {
                up |= side[a] == 1;
                down |= side[a] == -1;
                if side[a] == 0 {
                    section.push(f.uv(mesh.points[a]));
                }
                if crosses(a, b, &side) {
                    section.push(f.uv(section_point(a, b, &mesh.points)));
                }
            }
        }
        if !(up && down) {
            continue;
        }
        if let Mode::Through(cells) = &mode {
            cut[c] = cells.binary_search(&c).is_ok();
            continue;
        }
        let hull = convex_hull(&section);
        let meets = section.iter().any(|&q| f.inside(q) > tol)
            || (hull.len() >= 3
                && f.poly.iter().any(|&q| {
                    (0..hull.len()).all(|i| {
                        cross2(sub2(hull[(i + 1) % hull.len()], hull[i]), sub2(q, hull[i])) > tol
                    })
                }))
            || segments_cross(&hull, &f.poly, tol);
        cut[c] = meets;
    }
    stats.cut_cells = cut.iter().filter(|&&c| c).count();
    if stats.cut_cells == 0 && !shell_mode {
        return Err("internal: a through cut meets no cell".into());
    }
    // one intersection vertex on every crossing edge of a cut cell
    let mut splits: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (c, rings) in cell_rings.iter().enumerate() {
        if !cut[c] {
            continue;
        }
        for r in rings {
            for (a, b) in ring_edges(r) {
                if crosses(a, b, &side) && !splits.contains_key(&key(a, b)) {
                    let p = section_point(a, b, &mesh.points);
                    mesh.points.push(p);
                    side.push(0);
                    splits.insert(key(a, b), vec![mesh.points.len() - 1]);
                }
            }
        }
    }
    // 3. split the faces of the cut cells; chain their section loops
    let rings: Vec<Vec<Vec<usize>>> = cell_rings
        .iter()
        .map(|rs| rs.iter().map(|r| insert_splits(r, &splits)).collect())
        .collect();
    struct Piece {
        upper: Vec<Vec<usize>>,
        lower: Vec<Vec<usize>>,
        /// The piece each upper and each lower face goes to (a cell that is
        /// not convex can fall into several on a side).
        up_comp: Vec<usize>,
        lo_comp: Vec<usize>,
        pieces: [usize; 2],
        /// Directed edges of the upper and lower faces, with their piece.
        up_edge: HashMap<(usize, usize), usize>,
        lo_edge: HashMap<(usize, usize), usize>,
        /// Section loops, wound outward from the upper side (normal -n),
        /// with the upper and the lower piece they bound.
        loops: Vec<(Vec<usize>, usize, usize)>,
    }
    let mut pieces: HashMap<usize, Piece> = HashMap::new();
    // split faces, by their sorted nodes before the outline: (upper, lower)
    let mut face_splits: HashMap<Vec<usize>, Parts> = HashMap::new();
    for c in 0..rings.len() {
        if !cut[c] {
            continue;
        }
        let mut upper = Vec::new();
        let mut lower = Vec::new();
        let mut loop_edges = Vec::new();
        for r in &rings[c] {
            match split_ring(r, &side, &mesh.points, f.n)? {
                Split::Upper(r) => {
                    // an on-plane edge of an upper face bounds the section
                    for (a, b) in ring_edges(&r) {
                        if side[a] == 0 && side[b] == 0 {
                            loop_edges.push((b, a));
                        }
                    }
                    upper.push(r);
                }
                Split::Lower(r) => lower.push(r),
                Split::OnPlane => {
                    return Err("a cut cell has a face in the shell plane".into());
                }
                Split::Both { upper: u, lower: l } => {
                    // the chords and on-plane runs of the upper parts bound
                    // the section
                    for part in &u {
                        for (a, b) in ring_edges(part) {
                            if side[a] == 0 && side[b] == 0 {
                                loop_edges.push((b, a));
                            }
                        }
                    }
                    face_splits.insert(sorted(r), (u.clone(), l.clone()));
                    upper.extend(u);
                    lower.extend(l);
                }
            }
        }
        let (up_comp, nu) = components(&upper);
        let (lo_comp, nl) = components(&lower);
        let edge_map = |faces: &[Vec<usize>], comp: &[usize]| {
            let mut m = HashMap::new();
            for (r, &k) in faces.iter().zip(comp) {
                for e in ring_edges(r) {
                    m.insert(e, k);
                }
            }
            m
        };
        let up_edge = edge_map(&upper, &up_comp);
        let lo_edge = edge_map(&lower, &lo_comp);
        let mut loops = Vec::new();
        let uv = |v: usize| f.uv(mesh.points[v]);
        for lp in chain_loops(&loop_edges, &uv)? {
            // an upper face has the edge b -> a of the loop, a lower one a -> b
            let u = *ring_edges(&lp)
                .find_map(|(a, b)| up_edge.get(&(b, a)))
                .ok_or("internal: a section loop without its upper face")?;
            let l = *ring_edges(&lp)
                .find_map(|e| lo_edge.get(&e))
                .ok_or("internal: a section loop without its lower face")?;
            if loops.iter().any(|&(_, u0, l0)| (u0, l0) == (u, l)) {
                return Err("the shell cuts a cell into a ring; refine the mesh near it".into());
            }
            loops.push((lp, u, l));
        }
        pieces.insert(
            c,
            Piece {
                upper,
                lower,
                up_comp,
                lo_comp,
                pieces: [nu, nl],
                up_edge,
                lo_edge,
                loops,
            },
        );
    }
    // 4. split each section along the shell outline
    let mut outline = Outline {
        base: mesh.points.len(),
        ..Default::default()
    };
    let mut section_faces: HashMap<(usize, usize), Parts> = HashMap::new();
    let mut cut_ids: Vec<usize> = pieces.keys().copied().collect();
    cut_ids.sort_unstable();
    let mut narrow: Vec<(usize, usize)> = Vec::new();
    for &c in &cut_ids {
        for (li, (lp, _, _)) in pieces[&c].loops.iter().enumerate() {
            // the section seen from +n, counter-clockwise
            let ccw: Vec<usize> = lp.iter().rev().copied().collect();
            let q: Vec<[f64; 2]> = ccw.iter().map(|&v| f.uv(mesh.points[v])).collect();
            if polygon_area2(&q) <= 0.0 {
                return Err("a cell section has a hole or is wound the wrong way; refine the mesh near the shell".into());
            }
            if !shell_mode {
                section_faces.insert((c, li), (vec![ccw], Vec::new()));
                continue;
            }
            let (inside, outside) = split_section(f, &ccw, &mesh.points, &mut outline)?;
            if outside.len() > 1 {
                // an outline edge crossing this section: the integer part of
                // the place of a crossing point (corners have whole places)
                let j = outside[0]
                    .iter()
                    .filter(|&&v| v >= outline.base)
                    .map(|&v| outline.pos[v - outline.base])
                    .find(|p| p.fract() > 0.0)
                    .ok_or("internal: a strip crossing without an outline crossing")?;
                narrow.push((c, j as usize));
            }
            section_faces.insert((c, li), (inside, outside));
        }
    }
    // faces of uncut cells lying in the plane (snapped, or a plane through
    // mesh vertices) split along the outline too, as seen from the cell below
    let mut plane_faces: HashMap<Vec<usize>, Parts> = HashMap::new();
    if shell_mode {
        let mut uses: HashMap<Vec<usize>, usize> = HashMap::new();
        for (c, rs) in rings.iter().enumerate() {
            if cut[c] {
                continue;
            }
            for r in rs {
                if r.iter().all(|&v| side[v] == 0) {
                    *uses.entry(sorted(r)).or_default() += 1;
                }
            }
        }
        for (c, rs) in rings.iter().enumerate() {
            if cut[c] || rs.iter().flatten().map(|&v| side[v]).find(|&s| s != 0) != Some(-1) {
                continue;
            }
            for r in rs {
                let k = sorted(r);
                if uses.get(&k) != Some(&2) {
                    continue; // on the domain boundary, or not in the plane
                }
                let q: Vec<[f64; 2]> = r.iter().map(|&v| f.uv(mesh.points[v])).collect();
                if polygon_area2(&q) <= 0.0 {
                    return Err("a face in the shell plane is wound the wrong way".into());
                }
                let (inside, outside) = split_section(f, r, &mesh.points, &mut outline)?;
                if outside.len() > 1
                    && let Some(j) = outside[0]
                        .iter()
                        .filter(|&&v| v >= outline.base)
                        .map(|&v| outline.pos[v - outline.base])
                        .find(|p| p.fract() > 0.0)
                {
                    narrow.push((c, j as usize));
                }
                plane_faces.insert(k, (inside, outside));
            }
        }
    }
    if !narrow.is_empty() {
        // leave the mesh as it was
        mesh.points.truncate(np);
        for (v, p) in moved {
            mesh.points[v] = p;
        }
        return Ok(Outcome::Narrow(narrow));
    }
    for q in &outline.uv {
        mesh.points.push(f.point(*q));
        side.push(0);
    }
    let fix = |r: &[usize]| insert_splits(r, &outline.edge_points);
    let face_splits: HashMap<Vec<usize>, Parts> = face_splits
        .into_iter()
        .map(|(k, (u, l))| {
            let fixed = |p: Vec<Vec<usize>>| p.iter().map(|r| fix(r)).collect();
            (k, (fixed(u), fixed(l)))
        })
        .collect();
    // 5. rebuild the cells; a cut cell becomes its upper and lower pieces
    let mut first: Vec<usize> = Vec::with_capacity(mesh.cells.len());
    let mut next = 0;
    for c in 0..cut.len() {
        first.push(next);
        next += pieces.get(&c).map_or(1, |p| p.pieces[0] + p.pieces[1]);
    }
    let mut new_cells: Vec<Cell> = Vec::with_capacity(next);
    let mut cell_side: Vec<i8> = Vec::with_capacity(next);
    let old_cells = std::mem::take(&mut mesh.cells);
    for (c, cell) in old_cells.into_iter().enumerate() {
        if let Some(p) = pieces.get(&c) {
            let [nu, nl] = p.pieces;
            let mut up: Vec<Vec<Vec<usize>>> = vec![Vec::new(); nu];
            let mut lo: Vec<Vec<Vec<usize>>> = vec![Vec::new(); nl];
            for (r, &k) in p.upper.iter().zip(&p.up_comp) {
                up[k].push(fix(r));
            }
            for (r, &k) in p.lower.iter().zip(&p.lo_comp) {
                lo[k].push(fix(r));
            }
            for (li, &(_, u, l)) in p.loops.iter().enumerate() {
                let (inside, outside) = section_faces.remove(&(c, li)).unwrap_or_default();
                for s in inside.iter().chain(&outside) {
                    // s is counter-clockwise from +n: outward for the lower piece
                    lo[l].push(s.clone());
                    up[u].push(s.iter().rev().copied().collect());
                }
            }
            for faces in up {
                new_cells.push(poly_cell(faces, &mesh.points, cell.layer, cell.region));
                cell_side.push(1);
            }
            for faces in lo {
                new_cells.push(poly_cell(faces, &mesh.points, cell.layer, cell.region));
                cell_side.push(-1);
            }
        } else {
            // a neighbour of a cut cell takes the parts of the faces it split
            let mut rs: Vec<Vec<usize>> = Vec::new();
            for r in &rings[c] {
                let fixed = fix(r);
                if let Some((inside, outside)) = plane_faces.get(&sorted(r)) {
                    // counter-clockwise from +n: outward for the cell below
                    let below = r.iter().map(|&v| side[v]).all(|s| s == 0)
                        && rings[c]
                            .iter()
                            .flatten()
                            .map(|&v| side[v])
                            .find(|&s| s != 0)
                            == Some(-1);
                    for part in inside.iter().chain(outside) {
                        rs.push(if below {
                            part.clone()
                        } else {
                            part.iter().rev().copied().collect()
                        });
                    }
                } else if let Some((u, l)) = face_splits.get(&sorted(r)) {
                    let same = same_winding(&fixed, &u[0]);
                    for part in u.iter().chain(l) {
                        rs.push(if same {
                            part.clone()
                        } else {
                            part.iter().rev().copied().collect()
                        });
                    }
                } else {
                    rs.push(fixed);
                }
            }
            let s = rs
                .iter()
                .flatten()
                .map(|&v| side[v])
                .find(|&s| s != 0)
                .unwrap_or(0);
            if rs != cell.element.faces() {
                new_cells.push(poly_cell(rs, &mesh.points, cell.layer, cell.region));
            } else {
                new_cells.push(cell);
            }
            cell_side.push(s);
        }
    }
    // boundary faces of a cut cell go to the piece on their side, split
    // when they cross the plane
    let old_boundary = std::mem::take(&mut mesh.boundary);
    for b in old_boundary {
        for part in [b.nodes.clone()] {
            let nodes = insert_splits(&part, &splits);
            let mut push = |nodes: Vec<usize>, owner: usize| {
                mesh.boundary.push(BoundaryFace {
                    nodes: fix(&nodes),
                    patch: b.patch,
                    owner,
                });
            };
            let Some(p) = pieces.get(&b.owner) else {
                push(nodes, first[b.owner]);
                continue;
            };
            let owner_of = |r: &[usize], upper: bool| -> Result<usize, String> {
                let (map, offset) = if upper {
                    (&p.up_edge, 0)
                } else {
                    (&p.lo_edge, p.pieces[0])
                };
                ring_edges(r)
                    .find_map(|e| map.get(&e))
                    .map(|&k| first[b.owner] + offset + k)
                    .ok_or_else(|| {
                        "internal: a boundary face of a cut cell is not one of its faces"
                            .to_string()
                    })
            };
            match split_ring(&nodes, &side, &mesh.points, f.n)? {
                Split::Upper(r) => {
                    let o = owner_of(&r, true)?;
                    push(r, o);
                }
                Split::Lower(r) => {
                    let o = owner_of(&r, false)?;
                    push(r, o);
                }
                Split::OnPlane => {
                    return Err("a boundary face of a cut cell lies in the shell plane".into());
                }
                Split::Both { upper, lower } => {
                    for part in upper {
                        let o = owner_of(&part, true)?;
                        push(part, o);
                    }
                    for part in lower {
                        let o = owner_of(&part, false)?;
                        push(part, o);
                    }
                }
            }
        }
    }
    if !shell_mode {
        mesh.cells = new_cells;
        return Ok(Outcome::Done(stats));
    }
    // 6. double the on-plane vertices strictly inside the shell: the cells
    // below take the copies
    let mut back_copy: HashMap<usize, usize> = HashMap::new();
    let mut source_of: HashMap<usize, usize> = HashMap::new();
    for v in 0..mesh.points.len() {
        if side[v] == 0 && f.inside(f.uv(mesh.points[v])) > tol {
            mesh.points.push(mesh.points[v]);
            back_copy.insert(v, mesh.points.len() - 1);
            source_of.insert(mesh.points.len() - 1, v);
        }
    }
    stats.doubled_vertices = back_copy.len();
    for (c, cell) in new_cells.iter_mut().enumerate() {
        if cell_side[c] < 0 {
            relabel(cell, &back_copy);
        }
    }
    for b in &mut mesh.boundary {
        if cell_side[b.owner] < 0 {
            for v in &mut b.nodes {
                if let Some(&w) = back_copy.get(v) {
                    *v = w;
                }
            }
        }
    }
    let front = mesh.patches.len() + 1;
    mesh.patches.push((name.to_string(), "wall".into()));
    mesh.patches.push((back_name, "wall".into()));
    stats.patches = [front, front + 1];
    // every face in the plane and inside the shell is a face of the shell
    // (a boundary face of the domain in the plane stays what it is)
    let outer: std::collections::HashSet<Vec<usize>> =
        mesh.boundary.iter().map(|b| sorted(&b.nodes)).collect();
    for (c, cell) in new_cells.iter().enumerate() {
        if cell_side[c] == 0 {
            continue;
        }
        for r in cell.element.faces() {
            if outer.contains(&sorted(&r)) {
                continue;
            }
            let on_shell = r.iter().all(|&v| {
                let v0 = source_of.get(&v).copied().unwrap_or(v);
                v0 < side.len() && side[v0] == 0 && f.inside(f.uv(mesh.points[v])) >= -tol
            }) && {
                let centre = r.iter().fold(Vec3::ZERO, |s, &v| s + mesh.points[v]) / r.len() as f64;
                f.inside(f.uv(centre)) > tol
            };
            if !on_shell {
                continue;
            }
            let area = ring_area(&mesh.points, &r);
            let patch = if cell_side[c] > 0 {
                stats.front_faces += 1;
                stats.front_area += area;
                front
            } else {
                stats.back_faces += 1;
                stats.back_area += area;
                front + 1
            };
            mesh.boundary.push(BoundaryFace {
                nodes: r,
                patch,
                owner: c,
            });
        }
    }
    mesh.cells = new_cells;
    Ok(Outcome::Done(stats))
}

/// The connected groups of faces that share edges: the group of each face
/// and the number of groups.
fn components(faces: &[Vec<usize>]) -> (Vec<usize>, usize) {
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    fn root(p: &mut [usize], mut i: usize) -> usize {
        while p[i] != i {
            p[i] = p[p[i]];
            i = p[i];
        }
        i
    }
    let mut by_edge: HashMap<(usize, usize), usize> = HashMap::new();
    for (i, r) in faces.iter().enumerate() {
        for (a, b) in ring_edges(r) {
            if let Some(&j) = by_edge.get(&key(a, b)) {
                let (x, y) = (root(&mut parent, i), root(&mut parent, j));
                parent[x] = y;
            } else {
                by_edge.insert(key(a, b), i);
            }
        }
    }
    let mut label: HashMap<usize, usize> = HashMap::new();
    let mut comp = Vec::with_capacity(faces.len());
    for i in 0..faces.len() {
        let r = root(&mut parent, i);
        let n = label.len();
        comp.push(*label.entry(r).or_insert(n));
    }
    (comp, label.len())
}

/// Chain directed edges (wound clockwise seen from +n, the section on
/// their right) into closed loops. Where loops touch at a vertex, each
/// incoming edge goes on with the outgoing edge that keeps its section
/// tight, so touching loops stay apart.
fn chain_loops(
    edges: &[(usize, usize)],
    uv: &dyn Fn(usize) -> [f64; 2],
) -> Result<Vec<Vec<usize>>, String> {
    // reversed: counter-clockwise, the section on the left
    let mut out: HashMap<usize, Vec<usize>> = HashMap::new();
    for &(a, b) in edges {
        out.entry(b).or_default().push(a);
    }
    let mut starts: Vec<usize> = out.keys().copied().collect();
    starts.sort_unstable();
    let angle = |d: [f64; 2], o: [f64; 2]| {
        // clockwise angle from d to o, in (0, 2 pi]
        let t = (-cross2(d, o)).atan2(d[0] * o[0] + d[1] * o[1]);
        if t <= 0.0 {
            t + std::f64::consts::TAU
        } else {
            t
        }
    };
    let mut loops = Vec::new();
    for start in starts {
        while out.get(&start).is_some_and(|v| !v.is_empty()) {
            let mut lp = vec![start];
            let mut prev = start;
            let mut cur = out
                .get_mut(&start)
                .and_then(|v| v.pop())
                .ok_or("internal: empty edge list")?;
            while cur != start {
                if lp.len() > edges.len() {
                    return Err("the cut of a cell is not a set of loops".into());
                }
                lp.push(cur);
                let cands = out.get_mut(&cur).ok_or("the cut of a cell is open")?;
                let pick = match cands.len() {
                    0 => return Err("the cut of a cell is open".into()),
                    1 => 0,
                    _ => {
                        let (p, c) = (uv(prev), uv(cur));
                        let back = sub2(p, c);
                        let mut best = (f64::INFINITY, 0);
                        for (i, &o) in cands.iter().enumerate() {
                            let t = angle(back, sub2(uv(o), c));
                            if t < best.0 {
                                best = (t, i);
                            }
                        }
                        best.1
                    }
                };
                prev = cur;
                cur = cands.swap_remove(pick);
            }
            if lp.len() < 3 {
                return Err("the shell plane touches a cell along an edge only; refine the mesh near the shell".into());
            }
            // back to the clockwise winding of the input
            lp.reverse();
            loops.push(lp);
        }
    }
    Ok(loops)
}

fn sorted(r: &[usize]) -> Vec<usize> {
    let mut k = r.to_vec();
    k.sort_unstable();
    k.dedup();
    k
}

/// The points made where the shell outline meets the cell sections.
#[derive(Default)]
struct Outline {
    /// Id of the first new point.
    base: usize,
    /// New points, in plane coordinates.
    uv: Vec<[f64; 2]>,
    /// Their place along the outline (edge index + share of the edge).
    pos: Vec<f64>,
    /// New points on each section edge, from the lower to the higher id.
    edge_points: HashMap<(usize, usize), Vec<usize>>,
    /// Id of each polygon corner that lies inside a section.
    corners: HashMap<usize, usize>,
}

impl Outline {
    fn add(&mut self, q: [f64; 2], pos: f64) -> usize {
        self.uv.push(q);
        self.pos.push(pos);
        self.base + self.uv.len() - 1
    }
}

impl Frame {
    /// Place of a point on the outline (edge index + share), `None` when it
    /// is not within `tol` of it.
    fn perimeter(&self, q: [f64; 2], tol: f64) -> Option<f64> {
        let k = self.poly.len();
        let mut best: Option<(f64, f64)> = None;
        for j in 0..k {
            let (a, b) = (self.poly[j], self.poly[(j + 1) % k]);
            let e = sub2(b, a);
            let l = e[0].hypot(e[1]);
            let t = (sub2(q, a)[0] * e[0] + sub2(q, a)[1] * e[1]) / (l * l);
            if t * l < -tol || (t - 1.0) * l > tol {
                continue;
            }
            let d = self.edge_side(j, q).abs();
            if d <= tol && best.is_none_or(|(bd, _)| d < bd) {
                best = Some((d, (j as f64 + t.clamp(0.0, 1.0)).rem_euclid(k as f64)));
            }
        }
        best.map(|(_, p)| p)
    }
}

/// Whether `q` lies inside the simple polygon `s`, farther than `tol` from
/// its edges.
fn strictly_inside(s: &[[f64; 2]], q: [f64; 2], tol: f64) -> bool {
    let n = s.len();
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (s[i], s[(i + 1) % n]);
        let e = sub2(b, a);
        let l2 = e[0] * e[0] + e[1] * e[1];
        let t = ((q[0] - a[0]) * e[0] + (q[1] - a[1]) * e[1]) / l2;
        let c = [
            a[0] + e[0] * t.clamp(0.0, 1.0),
            a[1] + e[1] * t.clamp(0.0, 1.0),
        ];
        if (q[0] - c[0]).hypot(q[1] - c[1]) <= tol {
            return false;
        }
        if (a[1] > q[1]) != (b[1] > q[1]) && q[0] < a[0] + (q[1] - a[1]) / (b[1] - a[1]) * e[0] {
            inside = !inside;
        }
    }
    inside
}

/// Split a cell section (node ids, counter-clockwise from +n) into the part
/// inside the shell polygon and the part outside it, by a Weiler–Atherton
/// walk (Weiler and Atherton, *Hidden surface removal using polygon area
/// sorting*, SIGGRAPH 1977) over the section and the outline. Points where
/// the outline crosses a section edge are made once per edge, so the faces
/// that share the edge get the same ones.
fn split_section(
    f: &Frame,
    ring: &[usize],
    points: &[Vec3],
    o: &mut Outline,
) -> Result<Parts, String> {
    let tol = 1e-9 * f.scale;
    let k = f.poly.len();
    let uv = |v: usize, o: &Outline| {
        if v < o.base {
            f.uv(points[v])
        } else {
            o.uv[v - o.base]
        }
    };
    // the section with the outline crossings
    let n = ring.len();
    let mut aug: Vec<usize> = Vec::with_capacity(n + 4);
    for i in 0..n {
        let (a, b) = (ring[i], ring[(i + 1) % n]);
        aug.push(a);
        let e = key(a, b);
        if !o.edge_points.contains_key(&e) {
            let (pa, pb) = (uv(e.0, o), uv(e.1, o));
            let r = sub2(pb, pa);
            let lr = r[0].hypot(r[1]);
            let mut found: Vec<(f64, [f64; 2])> = Vec::new();
            for j in 0..k {
                let (c, d) = (f.poly[j], f.poly[(j + 1) % k]);
                let s = sub2(d, c);
                let ls = s[0].hypot(s[1]);
                let den = cross2(r, s);
                if den.abs() <= 1e-12 * lr * ls
                    || (f.edge_side(j, pa).abs() <= tol && f.edge_side(j, pb).abs() <= tol)
                {
                    continue; // parallel, or along the outline edge
                }
                let w = sub2(c, pa);
                let (ts, tp) = (cross2(w, s) / den, cross2(w, r) / den);
                if ts * lr <= tol
                    || (1.0 - ts) * lr <= tol
                    || tp * ls < -tol
                    || (tp - 1.0) * ls > tol
                {
                    continue;
                }
                if found.iter().any(|&(t0, _)| (t0 - ts).abs() * lr <= tol) {
                    continue; // a corner met by both of its edges
                }
                found.push((ts, [pa[0] + r[0] * ts, pa[1] + r[1] * ts]));
            }
            found.sort_by(|x, y| x.0.total_cmp(&y.0));
            let mut ids = Vec::with_capacity(found.len());
            for (_, q) in found {
                let pos = f
                    .perimeter(q, 2.0 * tol)
                    .ok_or("internal: an outline crossing is off the outline")?;
                ids.push(o.add(q, pos));
            }
            o.edge_points.insert(e, ids);
        }
        let mid = &o.edge_points[&e];
        if a < b {
            aug.extend(mid.iter().copied());
        } else {
            aug.extend(mid.iter().rev().copied());
        }
    }
    let m = aug.len();
    let q: Vec<[f64; 2]> = aug.iter().map(|&v| uv(v, o)).collect();
    // the outline points of the section and the corners inside it, in order
    // along the outline
    let mut along: Vec<(f64, usize)> = Vec::new();
    for (i, &v) in aug.iter().enumerate() {
        let pos = if v >= o.base {
            Some(o.pos[v - o.base])
        } else if f.inside(q[i]).abs() <= tol {
            f.perimeter(q[i], 2.0 * tol)
        } else {
            None
        };
        if let Some(p) = pos {
            along.push((p, v));
        }
    }
    let crossings = along.len();
    for (j, &corner) in f.poly.iter().enumerate() {
        if strictly_inside(&q, corner, tol) {
            let id = match o.corners.get(&j) {
                Some(&id) => id,
                None => {
                    let id = o.add(corner, j as f64);
                    o.corners.insert(j, id);
                    id
                }
            };
            along.push((j as f64, id));
        }
    }
    along.sort_by(|x, y| x.0.total_cmp(&y.0));
    along.dedup_by_key(|x| x.1);
    let corners_inside = along.len() > crossings;
    // which section edges lie inside the shell
    let mut edge_in = vec![false; m];
    for i in 0..m {
        let (a, b) = (q[i], q[(i + 1) % m]);
        let d = f.inside([(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0]);
        edge_in[i] = if d.abs() <= tol {
            // along the outline: inside when it runs the outline's way (both
            // interiors on its left), outside when it runs against it
            let mid = [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0];
            let pos = f
                .perimeter(mid, 2.0 * tol)
                .ok_or("internal: a section edge near the outline is off it")?;
            let j = (pos.floor() as usize).min(k - 1);
            let e = sub2(f.poly[(j + 1) % k], f.poly[j]);
            let r = sub2(b, a);
            e[0] * r[0] + e[1] * r[1] > 0.0
        } else {
            d > 0.0
        };
    }
    if crossings == 0 {
        return if edge_in.iter().all(|&x| x) {
            Ok((vec![aug], Vec::new()))
        } else if corners_inside {
            Err(
                "the shell lies inside one cell; refine the mesh near it or enlarge the shell"
                    .into(),
            )
        } else {
            Ok((Vec::new(), vec![aug]))
        };
    }
    let at_s: HashMap<usize, usize> = aug.iter().enumerate().map(|(i, &v)| (v, i)).collect();
    let at_p: HashMap<usize, usize> = along
        .iter()
        .enumerate()
        .map(|(i, &(_, v))| (v, i))
        .collect();
    let mut used_p = vec![false; along.len()];
    // walk the section edges of one kind; leave the section along the
    // outline, forward for the inside part and backward for the outside
    let mut walk = |inside: bool| -> Result<Vec<Vec<usize>>, String> {
        let mut used = vec![false; m];
        let mut faces = Vec::new();
        for s in 0..m {
            if edge_in[s] != inside || used[s] {
                continue;
            }
            let start = aug[s];
            let mut face = Vec::new();
            let mut v = start;
            loop {
                face.push(v);
                if face.len() > m + along.len() {
                    return Err("internal: the section split does not close".into());
                }
                let next = match at_s.get(&v) {
                    Some(&i) if edge_in[i] == inside => {
                        used[i] = true;
                        aug[(i + 1) % m]
                    }
                    _ => {
                        let &j = at_p
                            .get(&v)
                            .ok_or("internal: the section split left the outline")?;
                        used_p[j] = true;
                        let l = along.len();
                        along[if inside { (j + 1) % l } else { (j + l - 1) % l }].1
                    }
                };
                if next == start {
                    break;
                }
                v = next;
            }
            faces.push(face);
        }
        Ok(faces)
    };
    let inside_faces = walk(true)?;
    let outside_faces = walk(false)?;
    for (i, &(_, v)) in along.iter().enumerate() {
        if v >= o.base && !at_s.contains_key(&v) && !used_p[i] {
            return Err(
                "the shell lies inside one cell; refine the mesh near it or enlarge the shell"
                    .into(),
            );
        }
    }
    // the parts must tile the section
    let area = |r: &[usize]| polygon_area2(&r.iter().map(|&v| uv(v, o)).collect::<Vec<_>>());
    let whole = polygon_area2(&q);
    let mut sum = 0.0;
    for r in inside_faces.iter().chain(&outside_faces) {
        let a = area(r);
        if r.len() < 3 || a <= 0.0 {
            return Err("internal: a degenerate part of a cell section".into());
        }
        sum += a;
    }
    if (sum - whole).abs() > 1e-9 * whole.abs().max(tol * tol) {
        return Err("internal: the parts of a cell section do not tile it".into());
    }
    // a wall face needs a vertex strictly inside the shell, which the back
    // side takes a copy of; a face with every vertex on the outline (at a
    // corner, say) is fanned from its centre
    let mut walls = Vec::with_capacity(inside_faces.len());
    for r in inside_faces {
        if r.iter().any(|&v| f.inside(uv(v, o)) > tol) {
            walls.push(r);
            continue;
        }
        let pts: Vec<[f64; 2]> = r.iter().map(|&v| uv(v, o)).collect();
        let c = [
            pts.iter().map(|p| p[0]).sum::<f64>() / pts.len() as f64,
            pts.iter().map(|p| p[1]).sum::<f64>() / pts.len() as f64,
        ];
        let centre = o.add(c, -1.0);
        for (a, b) in ring_edges(&r) {
            walls.push(vec![centre, a, b]);
        }
    }
    Ok((walls, outside_faces))
}
fn relabel(cell: &mut Cell, map: &HashMap<usize, usize>) {
    let sub = |v: &mut usize| {
        if let Some(&w) = map.get(v) {
            *v = w;
        }
    };
    match &mut cell.element {
        Element::Tet(n) => n.iter_mut().for_each(sub),
        Element::Wedge(n) => n.iter_mut().for_each(sub),
        Element::Polyhedron(p) => {
            p.nodes.iter_mut().for_each(sub);
            p.nodes.sort_unstable();
            p.nodes.dedup();
            for r in &mut p.faces {
                r.iter_mut().for_each(sub);
            }
            if let Some(prism) = &mut p.prism {
                for w in prism {
                    w.iter_mut().for_each(sub);
                }
            }
        }
    }
}

fn same_winding(a: &[usize], b: &[usize]) -> bool {
    // b is a part of a after the splits: compare the direction of a shared edge
    for (x, y) in ring_edges(b) {
        for (p, q) in ring_edges(a) {
            if (p, q) == (x, y) {
                return true;
            }
            if (p, q) == (y, x) {
                return false;
            }
        }
    }
    true
}

fn ring_area(p: &[Vec3], r: &[usize]) -> f64 {
    let mut s = Vec3::ZERO;
    for (a, b) in ring_edges(r) {
        s = s + p[a].cross(p[b]);
    }
    0.5 * s.length()
}

fn poly_cell(faces: Vec<Vec<usize>>, points: &[Vec3], layer: usize, region: usize) -> Cell {
    let mut sum = Vec3::ZERO;
    let mut n = 0;
    for r in &faces {
        for &v in r {
            sum = sum + points[v];
            n += 1;
        }
    }
    let interior = sum / n.max(1) as f64;
    let volume = super::verify::shell_volume(points, &faces);
    Cell {
        element: Element::Polyhedron(Box::new(Polyhedron::new(faces, interior, None, volume))),
        layer,
        region,
    }
}

fn convex_hull(points: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut p = points.to_vec();
    p.sort_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
    p.dedup();
    if p.len() < 3 {
        return p;
    }
    let mut lower: Vec<[f64; 2]> = Vec::new();
    for &q in &p {
        while lower.len() >= 2
            && cross2(
                sub2(lower[lower.len() - 1], lower[lower.len() - 2]),
                sub2(q, lower[lower.len() - 2]),
            ) <= 0.0
        {
            lower.pop();
        }
        lower.push(q);
    }
    let mut upper: Vec<[f64; 2]> = Vec::new();
    for &q in p.iter().rev() {
        while upper.len() >= 2
            && cross2(
                sub2(upper[upper.len() - 1], upper[upper.len() - 2]),
                sub2(q, upper[upper.len() - 2]),
            ) <= 0.0
        {
            upper.pop();
        }
        upper.push(q);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

fn segments_cross(a: &[[f64; 2]], b: &[[f64; 2]], tol: f64) -> bool {
    let n = a.len();
    let m = b.len();
    if n < 2 || m < 2 {
        return false;
    }
    for i in 0..n {
        let (p, q) = (a[i], a[(i + 1) % n]);
        for j in 0..m {
            let (r, s) = (b[j], b[(j + 1) % m]);
            let d1 = cross2(sub2(q, p), sub2(r, p));
            let d2 = cross2(sub2(q, p), sub2(s, p));
            let d3 = cross2(sub2(s, r), sub2(p, r));
            let d4 = cross2(sub2(s, r), sub2(q, r));
            if d1 * d2 < -tol * tol && d3 * d4 < -tol * tol {
                return true;
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A unit cube of n³ hexahedral polyhedra, every boundary face on the
    /// wall patch.
    fn grid(n: usize) -> DuctMesh {
        let id = |i: usize, j: usize, k: usize| (k * (n + 1) + j) * (n + 1) + i;
        let h = 1.0 / n as f64;
        let mut points = Vec::new();
        for k in 0..=n {
            for j in 0..=n {
                for i in 0..=n {
                    points.push(Vec3::new(i as f64 * h, j as f64 * h, k as f64 * h));
                }
            }
        }
        let mut cells = Vec::new();
        let mut boundary = Vec::new();
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    let v = |a: usize, b: usize, c: usize| id(i + a, j + b, k + c);
                    // outward rings
                    let faces = vec![
                        (vec![v(0, 0, 0), v(0, 1, 0), v(1, 1, 0), v(1, 0, 0)], k == 0),
                        (
                            vec![v(0, 0, 1), v(1, 0, 1), v(1, 1, 1), v(0, 1, 1)],
                            k + 1 == n,
                        ),
                        (vec![v(0, 0, 0), v(1, 0, 0), v(1, 0, 1), v(0, 0, 1)], j == 0),
                        (
                            vec![v(0, 1, 0), v(0, 1, 1), v(1, 1, 1), v(1, 1, 0)],
                            j + 1 == n,
                        ),
                        (vec![v(0, 0, 0), v(0, 0, 1), v(0, 1, 1), v(0, 1, 0)], i == 0),
                        (
                            vec![v(1, 0, 0), v(1, 1, 0), v(1, 1, 1), v(1, 0, 1)],
                            i + 1 == n,
                        ),
                    ];
                    let owner = cells.len();
                    for (r, outer) in &faces {
                        if *outer {
                            boundary.push(BoundaryFace {
                                nodes: r.clone(),
                                patch: 3,
                                owner,
                            });
                        }
                    }
                    let rings: Vec<Vec<usize>> = faces.into_iter().map(|(r, _)| r).collect();
                    cells.push(poly_cell(rings, &points, 0, 0));
                }
            }
        }
        DuctMesh {
            points,
            cells,
            boundary,
            measurements: Measurements_default(),
            patches: super::super::classic_patches(),
        }
    }

    #[allow(non_snake_case)]
    fn Measurements_default() -> super::super::Measurements {
        super::super::Measurements::default()
    }

    fn square(z: f64, [x0, y0, x1, y1]: [f64; 4]) -> Vec<Vec3> {
        vec![
            Vec3::new(x0, y0, z),
            Vec3::new(x1, y0, z),
            Vec3::new(x1, y1, z),
            Vec3::new(x0, y1, z),
        ]
    }

    // Keep (SW-18): a baffle across the whole duct splits every cell it
    // crosses and leaves two faces per cut cell on its two sides. Breaks
    // when a split face is not shared by the neighbours (unpaired faces) or
    // a piece is wound inward; no other test cuts finished cells.
    #[test]
    fn a_spanning_baffle_cuts_every_cell_of_its_layer() {
        let mut mesh = grid(4);
        let shell = PlanarShell {
            name: "baffle".into(),
            polygon: square(0.43, [-1.0, -1.0, 2.0, 2.0]),
        };
        let stats = imprint_planar_shell(&mut mesh, &shell).unwrap();
        assert_eq!(stats.cut_cells, 16);
        assert_eq!((stats.front_faces, stats.back_faces), (16, 16));
        assert!((stats.front_area - 1.0).abs() < 1e-12, "{stats:?}");
        assert_eq!(mesh.cells.len(), 64 + 16);
        assert_eq!(mesh.patches.len(), 5);
        assert!((mesh.measurements.volume - 1.0).abs() < 1e-12);
    }

    // Keep (SW-18): a shell with free edges inside the fluid: the wall
    // faces cover exactly the polygon and the cells around the free edge
    // stay closed (the part of a section outside the polygon is an interior
    // face). Breaks when the section split loses or doubles an area.
    #[test]
    fn a_plate_with_free_edges_covers_its_polygon() {
        let mut mesh = grid(4);
        let shell = PlanarShell {
            name: "plate".into(),
            polygon: square(0.43, [0.13, 0.21, 0.71, 0.83]),
        };
        let stats = imprint_planar_shell(&mut mesh, &shell).unwrap();
        let area = 0.58 * 0.62;
        assert!((stats.polygon_area - area).abs() < 1e-12);
        assert!((stats.front_area - area).abs() < 1e-12, "{stats:?}");
        assert!((stats.back_area - area).abs() < 1e-12);
        assert_eq!(stats.cut_cells, 12);
        assert!((mesh.measurements.volume - 1.0).abs() < 1e-12);
    }

    // Keep (SW-18): a tilted triangle whose corners lie inside cells, and a
    // plane through grid vertices (snapped, no cut there). Breaks when a
    // corner inside a section or an on-plane vertex is mishandled.
    #[test]
    fn a_tilted_triangle_and_a_plane_through_vertices() {
        let mut mesh = grid(4);
        let shell = PlanarShell {
            name: "fin".into(),
            polygon: vec![
                Vec3::new(0.11, 0.17, 0.31),
                Vec3::new(0.87, 0.23, 0.52),
                Vec3::new(0.41, 0.79, 0.64),
            ],
        };
        let stats = imprint_planar_shell(&mut mesh, &shell).unwrap();
        assert!(
            (stats.front_area - stats.polygon_area).abs() < 1e-12 * stats.polygon_area.max(1.0),
            "{stats:?}"
        );
        let mut mesh = grid(4);
        let shell = PlanarShell {
            name: "mid".into(),
            polygon: square(0.5, [0.1, 0.1, 0.9, 0.9]),
        };
        let stats = imprint_planar_shell(&mut mesh, &shell).unwrap();
        assert!((stats.front_area - 0.64).abs() < 1e-12, "{stats:?}");
    }

    // Keep (SW-18): many tilted convex polygons at irregular places, each
    // inside the cube: the wall covers the polygon exactly. Breaks when a
    // corner, a strip crossing (through cuts) or a crossing near a vertex
    // goes wrong; the hand-picked cases above miss most of these.
    #[test]
    fn tilted_polygons_at_irregular_places_cover_themselves() {
        let mut seed = 0x2545_f491_4f6c_dd1du64;
        let mut rnd = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        for case in 0..40 {
            let mut mesh = grid(5);
            let centre = Vec3::new(0.3 + 0.4 * rnd(), 0.3 + 0.4 * rnd(), 0.3 + 0.4 * rnd());
            let n = Vec3::new(rnd() - 0.5, rnd() - 0.5, rnd() - 0.5 + 0.8);
            let n = n / n.length();
            let a = if n.x.abs() < 0.9 {
                Vec3::new(1.0, 0.0, 0.0)
            } else {
                Vec3::new(0.0, 1.0, 0.0)
            };
            let u = n.cross(a);
            let u = u / u.length();
            let v = n.cross(u);
            let sides = 3 + case % 5;
            let r = 0.12 + 0.15 * rnd();
            let phase = rnd();
            let polygon: Vec<Vec3> = (0..sides)
                .map(|i| {
                    let t = std::f64::consts::TAU * (i as f64 + phase) / sides as f64;
                    centre + u * (r * t.cos()) + v * (r * 0.7 * t.sin())
                })
                .collect();
            let shell = PlanarShell {
                name: "s".into(),
                polygon,
            };
            let stats = imprint_planar_shell(&mut mesh, &shell)
                .unwrap_or_else(|e| panic!("case {case}: {e}"));
            assert!(
                (stats.front_area / stats.polygon_area - 1.0).abs() < 1e-9,
                "case {case}: {stats:?}"
            );
        }
    }

    #[test]
    fn a_shell_on_the_domain_boundary_is_refused() {
        let mut mesh = grid(4);
        let shell = PlanarShell {
            name: "lid".into(),
            polygon: square(1.0, [0.1, 0.1, 0.9, 0.9]),
        };
        let e = imprint_planar_shell(&mut mesh, &shell).unwrap_err();
        assert!(e.contains("does not lie inside"), "{e}");
    }

    #[test]
    fn a_shell_inside_one_cell_is_refused() {
        let mut mesh = grid(4);
        let shell = PlanarShell {
            name: "tiny".into(),
            polygon: square(0.43, [0.3, 0.3, 0.4, 0.4]),
        };
        let e = imprint_planar_shell(&mut mesh, &shell).unwrap_err();
        assert!(e.contains("inside one cell"), "{e}");
    }
}
