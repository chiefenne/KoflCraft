//! The hand-off mesh of a block grid with non-conformal interfaces (BS-16,
//! step 3 of the 3D tetrahedron and prism insertion).
//!
//! An interface is a pair of patches that lie on one block face: the
//! *coarse* side is that whole face of one block, the *fine* side the
//! faces of the blocks on the other side, whose grid lines need not meet
//! the coarse ones. The two grids are overlaid in the coarse face's grid
//! parameters (its node (i, j) sits at (i / n_i, j / n_j)), with the fine
//! nodes placed there by inverting the coarse face's bilinear quads. Every
//! overlay face becomes one internal face between a coarse and a fine cell
//! ([`mesh_algorithms::block::overlay`]); a crossing point sits on the
//! coarse quad, so the coarse cells keep their shape. The points the
//! overlay puts inside an edge of either grid are put into every other face
//! of the mesh at that edge, so every cell closes. The cells along an
//! interface become polyhedra; all others stay hexahedra.

use std::collections::{HashMap, HashSet};

use mesh_algorithms::block::BlockGrid;
use mesh_algorithms::block::grid::HEX_FACES;
use mesh_algorithms::block::overlay::{Origin, P2, Tiling, overlay, ring_area2};
use mesh_core::block::BLOCK_FACES;
use mesh_core::geometry::Vec3;
use serde::Serialize;

use crate::block::domains;
use crate::handoff::{CellFaces, HandoffError, HandoffMesh, Patch};

/// A non-conformal interface: the indices of its two patches.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Interface {
    pub coarse: usize,
    pub fine: usize,
}

/// What the join of one interface did.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct JoinReport {
    pub coarse: String,
    pub fine: String,
    pub coarse_faces: usize,
    pub fine_faces: usize,
    /// Internal faces written for the interface.
    pub faces: usize,
    /// New points where grid lines of the two sides cross.
    pub crossings: usize,
    /// Largest snap in the face's parameters (fraction of the face).
    pub max_snap: f64,
    /// Largest |overlay area / cell face area − 1| in the parameters.
    pub max_area_error: f64,
    /// Largest distance of a fine node from the coarse face (the two sides
    /// of a curved face are different patches in 3D).
    pub max_gap: f64,
}

fn fail<T>(message: impl Into<String>) -> Result<T, HandoffError> {
    Err(HandoffError(message.into()))
}

/// The coarse face of an interface as a node grid.
struct CoarseFace {
    nu: usize,
    nv: usize,
    nodes: Vec<usize>,
}

impl CoarseFace {
    fn at(&self, iu: usize, iv: usize) -> usize {
        self.nodes[iv * (self.nu + 1) + iu]
    }

    /// The point at parameters (s, t), on the bilinear quad there.
    fn eval(&self, p: &[Vec3], st: P2) -> Vec3 {
        let x = (st[0] * self.nu as f64).clamp(0.0, self.nu as f64);
        let y = (st[1] * self.nv as f64).clamp(0.0, self.nv as f64);
        let iu = (x.floor() as usize).min(self.nu - 1);
        let iv = (y.floor() as usize).min(self.nv - 1);
        bilinear(p, self.quad(iu, iv), x - iu as f64, y - iv as f64)
    }

    fn quad(&self, iu: usize, iv: usize) -> [usize; 4] {
        [
            self.at(iu, iv),
            self.at(iu + 1, iv),
            self.at(iu + 1, iv + 1),
            self.at(iu, iv + 1),
        ]
    }

    /// The parameters of the point of the coarse face nearest to `x`, and
    /// its distance.
    fn invert(&self, p: &[Vec3], x: Vec3) -> (P2, f64) {
        let mut best = ([0.0, 0.0], f64::INFINITY);
        for iv in 0..self.nv {
            for iu in 0..self.nu {
                let q = self.quad(iu, iv);
                let (mut a, mut b) = (0.5, 0.5);
                for _ in 0..12 {
                    let r = x - bilinear(p, q, a, b);
                    let da = (p[q[1]] - p[q[0]]) * (1.0 - b) + (p[q[2]] - p[q[3]]) * b;
                    let db = (p[q[3]] - p[q[0]]) * (1.0 - a) + (p[q[2]] - p[q[1]]) * a;
                    let (m00, m01, m11) = (da.dot(da), da.dot(db), db.dot(db));
                    let (r0, r1) = (da.dot(r), db.dot(r));
                    let det = m00 * m11 - m01 * m01;
                    if det.abs() < f64::MIN_POSITIVE {
                        break;
                    }
                    a = (a + (m11 * r0 - m01 * r1) / det).clamp(0.0, 1.0);
                    b = (b + (m00 * r1 - m01 * r0) / det).clamp(0.0, 1.0);
                }
                let d = (x - bilinear(p, q, a, b)).length();
                if d < best.1 {
                    best = (
                        [
                            (iu as f64 + a) / self.nu as f64,
                            (iv as f64 + b) / self.nv as f64,
                        ],
                        d,
                    );
                }
            }
        }
        best
    }
}

fn bilinear(p: &[Vec3], q: [usize; 4], a: f64, b: f64) -> Vec3 {
    p[q[0]] * ((1.0 - a) * (1.0 - b))
        + p[q[1]] * (a * (1.0 - b))
        + p[q[2]] * (a * b)
        + p[q[3]] * ((1.0 - a) * b)
}

fn set_key(ring: &[usize]) -> Vec<usize> {
    let mut k = ring.to_vec();
    k.sort_unstable();
    k
}

/// Points to put inside mesh edges: undirected edge (lo, hi) → the points,
/// in order from lo.
#[derive(Default)]
struct Inserts(HashMap<(usize, usize), Vec<usize>>);

impl Inserts {
    fn add(&mut self, points: &[Vec3], u: usize, v: usize, inner: &[usize]) {
        if inner.is_empty() || u == v {
            return;
        }
        let (lo, hi) = (u.min(v), u.max(v));
        let mut list: Vec<usize> = if u == lo {
            inner.to_vec()
        } else {
            inner.iter().rev().copied().collect()
        };
        let entry = self.0.entry((lo, hi)).or_default();
        if !entry.is_empty() && *entry != list {
            list.extend(entry.iter().copied());
            let (p, d) = (points[lo], points[hi] - points[lo]);
            list.sort_by(|&x, &y| (points[x] - p).dot(d).total_cmp(&(points[y] - p).dot(d)));
            list.dedup();
        }
        *entry = list;
    }

    /// The ring with the points of its edges put in; true if it changed.
    fn expand(&self, ring: &mut Vec<usize>) -> bool {
        let n = ring.len();
        let mut out = Vec::with_capacity(n);
        let mut changed = false;
        let present: HashSet<usize> = ring.iter().copied().collect();
        for i in 0..n {
            let (u, v) = (ring[i], ring[(i + 1) % n]);
            out.push(u);
            if let Some(list) = self.0.get(&(u.min(v), u.max(v))) {
                let ordered: Vec<usize> = if u < v {
                    list.clone()
                } else {
                    list.iter().rev().copied().collect()
                };
                for w in ordered {
                    if !present.contains(&w) {
                        out.push(w);
                        changed = true;
                    }
                }
            }
        }
        *ring = out;
        changed
    }
}

/// Replace `ring`'s points by their aliases and drop repeats.
fn alias_ring(ring: &mut Vec<usize>, alias: &HashMap<usize, usize>) {
    for v in ring.iter_mut() {
        if let Some(&a) = alias.get(v) {
            *v = a;
        }
    }
    ring.dedup();
    while ring.len() > 1 && ring.first() == ring.last() {
        ring.pop();
    }
}

/// The FoehnCraft hand-off mesh of a block grid whose interfaces are
/// joined as the module notes say; the interface patches are not written.
/// `snap` is the overlay's relative snap distance.
pub fn handoff_joined(
    grid: &BlockGrid,
    patches: Vec<Patch>,
    interfaces: &[Interface],
    snap: f64,
) -> Result<(HandoffMesh, Vec<JoinReport>), HandoffError> {
    let (names, domain_of) = domains(grid);
    let mut points = grid.points.clone();
    let mut cells: Vec<CellFaces> = Vec::with_capacity(grid.hexes.len());
    for (b, r) in grid.blocks.iter().enumerate() {
        let count = r.cells.iter().product::<usize>();
        for h in &grid.hexes[r.first_cell..r.first_cell + count] {
            cells.push(CellFaces {
                faces: HEX_FACES.iter().map(|f| f.map(|s| h[s]).to_vec()).collect(),
                zone: domain_of[b],
            });
        }
    }
    let mut on_interface: HashSet<usize> = HashSet::new();
    for i in interfaces {
        if i.coarse >= patches.len() || i.fine >= patches.len() || i.coarse == i.fine {
            return fail("an interface names an unknown patch or one patch twice");
        }
        if !on_interface.insert(i.coarse) || !on_interface.insert(i.fine) {
            return fail("a patch is on two interfaces");
        }
    }
    let mut alias: HashMap<usize, usize> = HashMap::new();
    let mut inserts = Inserts::default();
    let mut reports = Vec::new();
    for i in interfaces {
        let name = |p: usize| patches[p].name.clone();
        // the coarse block face
        let faces: Vec<(usize, usize)> = grid
            .face_patches
            .iter()
            .enumerate()
            .flat_map(|(b, f)| {
                f.iter()
                    .enumerate()
                    .filter(|(_, p)| **p == Some(i.coarse))
                    .map(move |(g, _)| (b, g))
            })
            .collect();
        let [(o, g)] = faces[..] else {
            return fail(format!(
                "the coarse side {} of an interface must be one block face, not {}",
                name(i.coarse),
                faces.len()
            ));
        };
        let (fixed, side, u, v) = BLOCK_FACES[g];
        let n = grid.blocks[o].cells;
        let (nu, nv) = (n[u], n[v]);
        let mut nodes = Vec::with_capacity((nu + 1) * (nv + 1));
        let mut param: HashMap<usize, P2> = HashMap::new();
        for iv in 0..=nv {
            for iu in 0..=nu {
                let mut c = [0usize; 3];
                c[fixed] = side * n[fixed];
                c[u] = iu;
                c[v] = iv;
                let id = grid.block_corner(o, c[0], c[1], c[2]);
                nodes.push(id);
                param.insert(id, [iu as f64 / nu as f64, iv as f64 / nv as f64]);
            }
        }
        let coarse = CoarseFace { nu, nv, nodes };
        // the quads of each side with their cells
        let quads = |p: usize| -> Vec<([usize; 4], usize)> {
            grid.boundary
                .iter()
                .filter(|(_, patch, _)| *patch == p)
                .map(|&(q, _, c)| (q, c))
                .collect()
        };
        let (qa, qb) = (quads(i.coarse), quads(i.fine));
        if qb.is_empty() {
            return fail(format!("the fine side {} has no faces", name(i.fine)));
        }
        // tiling A in the coarse parameters
        let mut ga: Vec<usize> = Vec::new();
        let mut la: HashMap<usize, usize> = HashMap::new();
        let mut ta = Tiling::default();
        for (q, _) in &qa {
            ta.cells.push(
                q.iter()
                    .map(|&id| {
                        *la.entry(id).or_insert_with(|| {
                            ga.push(id);
                            ta.points.push(param[&id]);
                            ga.len() - 1
                        })
                    })
                    .collect(),
            );
        }
        // tiling B: fine nodes at the parameters of the nearest coarse point
        let mut gb: Vec<usize> = Vec::new();
        let mut lb: HashMap<usize, usize> = HashMap::new();
        let mut tb = Tiling::default();
        let mut max_gap: f64 = 0.0;
        for (q, _) in &qb {
            let mut ring = Vec::with_capacity(4);
            for &id in q {
                let l = match lb.get(&id) {
                    Some(&l) => l,
                    None => {
                        let st = match param.get(&id) {
                            Some(&st) => st,
                            None => {
                                let (st, d) = coarse.invert(&points, points[id]);
                                max_gap = max_gap.max(d);
                                st
                            }
                        };
                        gb.push(id);
                        tb.points.push(st);
                        lb.insert(id, gb.len() - 1);
                        gb.len() - 1
                    }
                };
                ring.push(l);
            }
            tb.cells.push(ring);
        }
        // fine nodes on the border lie on the face's border exactly
        let mut used: HashMap<(usize, usize), usize> = HashMap::new();
        for ring in &tb.cells {
            for k in 0..4 {
                let (x, y) = (ring[k], ring[(k + 1) % 4]);
                *used.entry((x.min(y), x.max(y))).or_default() += 1;
            }
        }
        for (&(x, y), &count) in &used {
            if count != 1 {
                continue;
            }
            for w in [x, y] {
                if param.contains_key(&gb[w]) {
                    continue;
                }
                let p = &mut tb.points[w];
                let d = [p[0], 1.0 - p[0], p[1], 1.0 - p[1]];
                let k = (0..4).min_by(|&a, &b| d[a].total_cmp(&d[b])).unwrap_or(0);
                match k {
                    0 => p[0] = 0.0,
                    1 => p[0] = 1.0,
                    2 => p[1] = 0.0,
                    _ => p[1] = 1.0,
                }
            }
        }
        // counterclockwise in the parameters
        let flip = |t: &mut Tiling| -> bool {
            let neg = t
                .cells
                .first()
                .is_some_and(|r| ring_area2(&t.points, r) < 0.0);
            if neg {
                for r in &mut t.cells {
                    r.reverse();
                }
            }
            neg
        };
        let (fa, fb) = (flip(&mut ta), flip(&mut tb));
        if fa == fb {
            return fail(format!(
                "the sides {} and {} of an interface face the same way",
                name(i.coarse),
                name(i.fine)
            ));
        }
        let o2 = overlay(&ta, &tb, snap).map_err(|e| {
            HandoffError(format!(
                "interface {} / {}: {e}",
                name(i.coarse),
                name(i.fine)
            ))
        })?;
        // the overlay's points in the mesh
        let mut crossings = 0;
        let global: Vec<usize> = o2
            .origin
            .iter()
            .zip(&o2.points)
            .map(|(org, &st)| match *org {
                Origin::A(k) => ga[k],
                Origin::B(k) => gb[k],
                Origin::Both(ka, kb) => {
                    if ga[ka] != gb[kb] {
                        alias.insert(gb[kb], ga[ka]);
                    }
                    ga[ka]
                }
                Origin::Crossing => {
                    crossings += 1;
                    points.push(coarse.eval(&points, st));
                    points.len() - 1
                }
            })
            .collect();
        // the cell faces on the interface become the overlay faces
        let mut by_a: Vec<Vec<Vec<usize>>> = vec![Vec::new(); qa.len()];
        let mut by_b: Vec<Vec<Vec<usize>>> = vec![Vec::new(); qb.len()];
        for f in &o2.faces {
            let ring: Vec<usize> = f.ring.iter().map(|&k| global[k]).collect();
            let rev: Vec<usize> = ring.iter().rev().copied().collect();
            by_a[f.a].push(if fa { rev.clone() } else { ring.clone() });
            by_b[f.b].push(if fb { rev } else { ring });
        }
        for (list, by) in [(&qa, by_a), (&qb, by_b)] {
            for ((q, c), polys) in list.iter().zip(by) {
                let key = set_key(q);
                let Some(at) = cells[*c].faces.iter().position(|r| set_key(r) == key) else {
                    return fail(format!("cell {c} has no face {q:?}"));
                };
                cells[*c].faces.remove(at);
                cells[*c].faces.extend(polys);
            }
        }
        for (edges, g) in [(&o2.a_edges, &ga), (&o2.b_edges, &gb)] {
            for (&(x, y), inner) in edges {
                let inner: Vec<usize> = inner.iter().map(|&k| global[k]).collect();
                inserts.add(&points, g[x], g[y], &inner);
            }
        }
        reports.push(JoinReport {
            coarse: name(i.coarse),
            fine: name(i.fine),
            coarse_faces: qa.len(),
            fine_faces: qb.len(),
            faces: o2.faces.len(),
            crossings,
            max_snap: o2.max_snap,
            max_area_error: o2.max_area_error,
            max_gap,
        });
    }
    // the other patches, renumbered
    let mut renumber = vec![usize::MAX; patches.len()];
    let mut kept = Vec::new();
    for (p, patch) in patches.into_iter().enumerate() {
        if !on_interface.contains(&p) {
            renumber[p] = kept.len();
            kept.push(patch);
        }
    }
    let mut boundary: Vec<(Vec<usize>, usize)> = grid
        .boundary
        .iter()
        .filter(|(_, p, _)| !on_interface.contains(p))
        .map(|&(q, p, _)| (q.to_vec(), renumber[p]))
        .collect();
    // aliases inside the insert lists, then into every ring
    let mut aliased = Inserts::default();
    for (&(x, y), list) in &inserts.0 {
        let a = |v: usize| alias.get(&v).copied().unwrap_or(v);
        let inner: Vec<usize> = list.iter().map(|&w| a(w)).collect();
        aliased.add(&points, a(x), a(y), &inner);
    }
    for c in &mut cells {
        for r in &mut c.faces {
            alias_ring(r, &alias);
        }
    }
    for (r, _) in &mut boundary {
        alias_ring(r, &alias);
    }
    for _ in 0..4 {
        let mut changed = false;
        for c in &mut cells {
            for r in &mut c.faces {
                changed |= aliased.expand(r);
            }
        }
        for (r, _) in &mut boundary {
            changed |= aliased.expand(r);
        }
        if !changed {
            break;
        }
    }
    let mesh = HandoffMesh::from_cells(points, cells, names, &boundary, kept)?;
    Ok((mesh, reports))
}

/// The cells of a joined mesh that are no longer plain hexahedra: those with
/// more than six faces or a face with more than four points.
pub fn polyhedral_cells(mesh: &HandoffMesh) -> Vec<bool> {
    let list = mesh.cell_face_list();
    (0..mesh.cell_count())
        .map(|c| list.len(c) != 6 || list.cell(c).any(|(f, _)| mesh.face(f).len() != 4))
        .collect()
}

/// MSH of a joined block mesh (card D of step 3): hexahedra stay
/// hexahedra; every polyhedron is written as pyramids and tetrahedra from
/// its centroid over its faces (a face with more than four points as a fan
/// from its own centroid), each piece carrying its polyhedron as parent,
/// as KoflCraft's polyhedral meshes are written. Boundary polygons are fans
/// from the same face centroids, so the pieces match them.
pub fn write_gmsh_joined(
    grid: &BlockGrid,
    mesh: &HandoffMesh,
    version: crate::msh::GmshVersion,
    out: &mut dyn std::io::Write,
) -> std::io::Result<()> {
    use crate::msh::Model;
    let poly = polyhedral_cells(mesh);
    let measures = mesh.cell_volumes();
    let mut points = mesh.points.clone();
    let mut centre: HashMap<usize, usize> = HashMap::new();
    let mut face_centre = |f: usize, points: &mut Vec<Vec3>| -> usize {
        *centre.entry(f).or_insert_with(|| {
            let r = mesh.face(f);
            let c = r.iter().fold(Vec3::ZERO, |s, &v| s + mesh.points[v]) / r.len() as f64;
            points.push(c);
            points.len() - 1
        })
    };
    // pieces as (kind, nodes, cell)
    let mut pieces: Vec<(u8, Vec<usize>, usize)> = Vec::new();
    let mut boundary: Vec<(u8, Vec<usize>, usize)> = Vec::new();
    let list = mesh.cell_face_list();
    for (p, _) in mesh.patches.iter().enumerate() {
        for f in mesh.patch_faces(p) {
            let r = mesh.face(f).to_vec();
            match r.len() {
                3 => boundary.push((2, r, p)),
                4 => boundary.push((3, r, p)),
                n => {
                    let fc = face_centre(f, &mut points);
                    for k in 0..n {
                        boundary.push((2, vec![fc, r[k], r[(k + 1) % n]], p));
                    }
                }
            }
        }
    }
    for c in 0..mesh.cell_count() {
        if !poly[c] {
            continue;
        }
        points.push(measures[c].1);
        let apex = points.len() - 1;
        for (f, outward) in list.cell(c) {
            let mut r = mesh.face(f).to_vec();
            if !outward {
                r.reverse();
            }
            // pieces need the face wound inward, towards the apex
            r.reverse();
            match r.len() {
                3 => pieces.push((4, vec![r[0], r[1], r[2], apex], c)),
                4 => pieces.push((7, vec![r[0], r[1], r[2], r[3], apex], c)),
                n => {
                    let fc = face_centre(f, &mut points);
                    for k in 0..n {
                        pieces.push((4, vec![fc, r[k], r[(k + 1) % n], apex], c));
                    }
                }
            }
        }
    }
    let mut model = Model::new(points);
    let np = mesh.patches.len();
    for (p, patch) in mesh.patches.iter().enumerate() {
        model.group(2, p + 1, &patch.name);
    }
    for (d, name) in mesh.zones.iter().enumerate() {
        model.group(3, np + 1 + d, name);
    }
    for (kind, nodes, p) in &boundary {
        model.element(*kind, p + 1, p + 1, nodes)?;
    }
    let any_poly = poly.iter().any(|&x| x);
    let mut next = 0;
    for (c, h) in grid.hexes.iter().enumerate() {
        let zone = mesh.cell_zone[c];
        let parent = any_poly.then_some(c);
        if poly[c] {
            while next < pieces.len() && pieces[next].2 == c {
                let (kind, nodes, _) = &pieces[next];
                model.element_of(*kind, np + 1 + zone, zone + 1, nodes, parent)?;
                next += 1;
            }
        } else {
            model.element_of(5, np + 1 + zone, zone + 1, h, parent)?;
        }
    }
    model.write(version, out)
}
