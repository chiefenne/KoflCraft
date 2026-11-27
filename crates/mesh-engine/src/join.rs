//! Join the polyhedral mesh of the region behind the cuts to the swept
//! hexahedra (MD-25, BS-17; strategy card 4 of 2026-11-25): every hex face
//! on a cut is replaced by the polyhedral cap faces that lie in it, and
//! every point of those faces on an edge of the hex face is put into the
//! other faces along that edge, so all cells close edge to edge. The cut's
//! cap was built from the cut section ([`CutCap`]), so the cap faces tile
//! each hex face exactly: no tolerance weld beyond the rounding of the
//! cap's similarity.
//!
//! [`CutCap`]: mesh_algorithms::fields::sweep::CutCap

use std::collections::HashMap;

use mesh_algorithms::duct::{DuctMesh, Element};
use mesh_core::geometry::Vec3;
use mesh_io::handoff::CellFaces;

/// One cut seen from the merged mesh.
pub(crate) struct Cut {
    /// The cap's patch id in the polyhedral mesh (1 or 2).
    pub patch: usize,
    pub normal: Vec3,
    /// Merged ids of the nodes of the cut section.
    pub nodes: Vec<usize>,
    /// The hex faces on the cut: (cell, face index in the cell).
    pub faces: Vec<(usize, usize)>,
}

/// What the join did.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Joined {
    pub cells: usize,
    pub pieces: usize,
    pub edge_points: usize,
    pub welded: usize,
    pub max_weld: f64,
    /// Largest |pieces area / hex face area - 1|.
    pub area_error: f64,
}

fn area_vector(points: &[Vec3], ring: &[usize]) -> Vec3 {
    let mut s = Vec3::ZERO;
    for k in 0..ring.len() {
        s = s + points[ring[k]].cross(points[ring[(k + 1) % ring.len()]]);
    }
    s * 0.5
}

fn cell_faces(element: &Element) -> Vec<Vec<usize>> {
    match element {
        Element::Tet(n) => [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]]
            .iter()
            .map(|f| f.iter().map(|&i| n[i]).collect())
            .collect(),
        Element::Wedge(n) => [
            vec![0, 2, 1],
            vec![3, 4, 5],
            vec![0, 1, 4, 3],
            vec![1, 2, 5, 4],
            vec![2, 0, 3, 5],
        ]
        .iter()
        .map(|f| f.iter().map(|&i| n[i]).collect())
        .collect(),
        Element::Polyhedron(p) => p.faces.clone(),
    }
}

fn key(ring: &[usize]) -> Vec<usize> {
    let mut k = ring.to_vec();
    k.sort_unstable();
    k
}

/// Two in-plane axes for `normal`.
fn axes(normal: Vec3) -> (Vec3, Vec3) {
    let a = if normal.x.abs() < 0.6 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = a - normal * a.dot(normal);
    let u = u / u.length();
    (u, normal.cross(u))
}

fn inside(poly: &[(f64, f64)], p: (f64, f64)) -> bool {
    let mut c = false;
    let n = poly.len();
    for i in 0..n {
        let (a, b) = (poly[i], poly[(i + 1) % n]);
        if (a.1 > p.1) != (b.1 > p.1) {
            let x = a.0 + (p.1 - a.1) / (b.1 - a.1) * (b.0 - a.0);
            if p.0 < x {
                c = !c;
            }
        }
    }
    c
}

/// Append the polyhedral mesh `duct` to the swept cells and join it at the
/// cuts. `patch_of` maps a polyhedral boundary patch id to a merged patch
/// (`None` for a cap of a cut, which becomes interior); `sweep_cells` is the
/// number of swept cells at the front of `cells`; `tol` the weld distance.
#[allow(clippy::too_many_arguments)]
pub(crate) fn join(
    points: &mut Vec<Vec3>,
    cells: &mut Vec<CellFaces>,
    boundary: &mut Vec<(Vec<usize>, usize)>,
    sweep_cells: usize,
    cuts: &[Cut],
    duct: &DuctMesh,
    patch_of: &dyn Fn(usize) -> Option<usize>,
    tol: f64,
) -> Result<Joined, String> {
    let mut out = Joined::default();
    // the polyhedral cells, wound outward
    let poly: Vec<CellFaces> = duct
        .cells
        .iter()
        .map(|c| {
            CellFaces {
                faces: cell_faces(&c.element),
                zone: 0,
            }
            .outward(&duct.points)
        })
        .collect();
    // weld the points of the cap faces onto the nodes of their cut
    let cell = 4.0 * tol;
    let quant = |x: Vec3| -> [i64; 3] {
        [
            (x.x / cell).floor() as i64,
            (x.y / cell).floor() as i64,
            (x.z / cell).floor() as i64,
        ]
    };
    let mut grid: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    for cut in cuts {
        for &v in &cut.nodes {
            grid.entry(quant(points[v])).or_default().push(v);
        }
    }
    let mut map: Vec<Option<usize>> = vec![None; duct.points.len()];
    for f in &duct.boundary {
        if !cuts.iter().any(|c| c.patch == f.patch) {
            continue;
        }
        for &v in &f.nodes {
            if map[v].is_some() {
                continue;
            }
            let x = duct.points[v];
            let q = quant(x);
            let mut best: Option<(f64, usize)> = None;
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(list) = grid.get(&[q[0] + dx, q[1] + dy, q[2] + dz]) {
                            for &w in list {
                                let d = (points[w] - x).length();
                                if d <= tol && best.is_none_or(|b| d < b.0) {
                                    best = Some((d, w));
                                }
                            }
                        }
                    }
                }
            }
            if let Some((d, w)) = best {
                map[v] = Some(w);
                out.welded += 1;
                out.max_weld = out.max_weld.max(d);
            }
        }
    }
    let mut ids = vec![0; duct.points.len()];
    for (v, m) in map.iter().enumerate() {
        ids[v] = match m {
            Some(w) => *w,
            None => {
                points.push(duct.points[v]);
                points.len() - 1
            }
        };
    }
    // the cap faces, wound out of their polyhedral cell, by cut
    let mut pieces: Vec<Vec<Vec<usize>>> = vec![Vec::new(); cuts.len()];
    for f in &duct.boundary {
        let Some(k) = cuts.iter().position(|c| c.patch == f.patch) else {
            if let Some(p) = patch_of(f.patch) {
                boundary.push((f.nodes.iter().map(|&v| ids[v]).collect(), p));
            }
            continue;
        };
        let want = key(&f.nodes);
        let ring = poly[f.owner]
            .faces
            .iter()
            .find(|r| key(r) == want)
            .ok_or("a cap face is not a face of its cell")?;
        pieces[k].push(ring.iter().map(|&v| ids[v]).collect());
    }
    let poly_start = cells.len();
    for c in poly {
        cells.push(CellFaces {
            faces: c
                .faces
                .iter()
                .map(|r| r.iter().map(|&v| ids[v]).collect())
                .collect(),
            zone: c.zone,
        });
    }
    out.cells = cells.len() - poly_start;
    // each piece into the hex face it lies in
    let mut replace: HashMap<(usize, usize), Vec<Vec<usize>>> = HashMap::new();
    let mut edge_points: HashMap<(usize, usize), Vec<(f64, usize)>> = HashMap::new();
    for (k, cut) in cuts.iter().enumerate() {
        let (u, w) = axes(cut.normal);
        let flat = |x: Vec3| (x.dot(u), x.dot(w));
        let polys: Vec<Vec<(f64, f64)>> = cut
            .faces
            .iter()
            .map(|&(c, fi)| {
                cells[c].faces[fi]
                    .iter()
                    .map(|&v| flat(points[v]))
                    .collect()
            })
            .collect();
        // a grid over the hex faces
        let (mut lo, mut hi) = (
            (f64::INFINITY, f64::INFINITY),
            (f64::NEG_INFINITY, f64::NEG_INFINITY),
        );
        for p in polys.iter().flatten() {
            lo = (lo.0.min(p.0), lo.1.min(p.1));
            hi = (hi.0.max(p.0), hi.1.max(p.1));
        }
        let side = ((polys.len() as f64).sqrt().ceil() as usize).max(1);
        let size = ((hi.0 - lo.0).max(hi.1 - lo.1) / side as f64).max(1e-300);
        let at = |p: (f64, f64)| {
            (
                (((p.0 - lo.0) / size) as usize).min(side - 1),
                (((p.1 - lo.1) / size) as usize).min(side - 1),
            )
        };
        let mut bins: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (q, poly) in polys.iter().enumerate() {
            let (mut a, mut b) = ((usize::MAX, usize::MAX), (0, 0));
            for &p in poly {
                let c = at(p);
                a = (a.0.min(c.0), a.1.min(c.1));
                b = (b.0.max(c.0), b.1.max(c.1));
            }
            for x in a.0..=b.0 {
                for y in a.1..=b.1 {
                    bins.entry((x, y)).or_default().push(q);
                }
            }
        }
        let mut area = vec![0.0; polys.len()];
        for ring in &pieces[k] {
            let c = ring.iter().fold(Vec3::ZERO, |s, &v| s + points[v]) / ring.len() as f64;
            let p = flat(c);
            let q = bins
                .get(&at(p))
                .and_then(|list| list.iter().copied().find(|&q| inside(&polys[q], p)))
                .ok_or("a cap face lies in no hex face of its cut")?;
            area[q] += area_vector(points, ring).dot(cut.normal);
            let (cell, fi) = cut.faces[q];
            let mut face = ring.clone();
            face.reverse();
            replace.entry((cell, fi)).or_default().push(face);
            out.pieces += 1;
            // points on the edges of the hex face
            let hex = &cells[cell].faces[fi];
            for i in 0..hex.len() {
                let (a, b) = (hex[i], hex[(i + 1) % hex.len()]);
                let (pa, pb) = (points[a], points[b]);
                let e = pb - pa;
                for &v in ring {
                    if v == a || v == b {
                        continue;
                    }
                    let t = (points[v] - pa).dot(e) / e.dot(e);
                    if t > 0.0 && t < 1.0 && (pa + e * t - points[v]).length() <= tol {
                        let (lo, hi, t) = if a < b { (a, b, t) } else { (b, a, 1.0 - t) };
                        let list = edge_points.entry((lo, hi)).or_default();
                        if !list.iter().any(|&(_, x)| x == v) {
                            list.push((t, v));
                        }
                    }
                }
            }
        }
        for (q, &(c, fi)) in cut.faces.iter().enumerate() {
            let own = area_vector(points, &cells[c].faces[fi])
                .dot(cut.normal)
                .abs();
            if own > 0.0 {
                out.area_error = out.area_error.max((area[q].abs() / own - 1.0).abs());
            }
        }
    }
    for list in edge_points.values_mut() {
        list.sort_by(|a, b| a.0.total_cmp(&b.0));
        out.edge_points += list.len();
    }
    // the hex faces on the cuts become their pieces
    // per cell: (face index, the pieces in its place)
    type Pieces = Vec<(usize, Vec<Vec<usize>>)>;
    let mut by_cell: HashMap<usize, Pieces> = HashMap::new();
    for ((c, fi), list) in replace {
        by_cell.entry(c).or_default().push((fi, list));
    }
    for (c, mut list) in by_cell {
        list.sort_by_key(|p| std::cmp::Reverse(p.0));
        for (fi, faces) in list {
            cells[c].faces.remove(fi);
            cells[c].faces.extend(faces);
        }
    }
    // put the edge points into every swept face along those edges
    let splice = |ring: &mut Vec<usize>| {
        let mut outr = Vec::with_capacity(ring.len() + 2);
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            outr.push(a);
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            if let Some(list) = edge_points.get(&(lo, hi)) {
                if a == lo {
                    outr.extend(list.iter().map(|&(_, v)| v));
                } else {
                    outr.extend(list.iter().rev().map(|&(_, v)| v));
                }
            }
        }
        *ring = outr;
    };
    if !edge_points.is_empty() {
        for c in cells.iter_mut().take(sweep_cells) {
            for ring in c.faces.iter_mut() {
                splice(ring);
            }
        }
        for (ring, _) in boundary.iter_mut() {
            splice(ring);
        }
    }
    Ok(out)
}
