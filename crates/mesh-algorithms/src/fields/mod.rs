//! Scalar fields on a tetrahedral scaffold that guide block-dominant
//! meshing (which parts of a domain a block template can fill).
//!
//! * Wall distance by the heat method (Crane, Weischedel & Wardetzky 2013,
//!   *Geodesics in Heat*, ACM TOG 32(5)): one short heat flow from the
//!   walls, its normalised gradient, one Poisson solve. The time step
//!   t = (m h)^2 with mean edge length h smooths the medial ridges, so the
//!   skeleton has far fewer spurious branches than the exact distance
//!   (prototype round 2: 2093 to 62 ends on the seven-tube furnace).
//! * Flow potential: the harmonic function 0 on inlets and 1 on outlets
//!   with no flux through walls, symmetry planes and cuts, the classical
//!   streamwise coordinate of potential flow and of harmonic skeletons
//!   (Ni, Garland & Hart 2004, *Fair Morse Functions for Extracting the
//!   Topological Structure of a Surface Mesh*).
//! * Fiedler vector: the lowest non-constant mode of the Laplacian with
//!   no flux everywhere (Fiedler 1973; Lévy 2006, *Laplace-Beltrami
//!   Eigenfunctions: Towards an Algorithm that Understands Geometry*), the
//!   natural long axis of a region without ports.
//!
//! All three use linear finite elements: the cotangent stiffness matrix
//! and the lumped (barycentric) mass of each node.

#![allow(clippy::needless_range_loop)] // parallel per-node and per-tet arrays

pub mod regions;
pub mod sparse;

use mesh_core::geometry::Vec3;
use sparse::{Csr, SolveStats};
use std::collections::HashMap;

/// What a boundary face is for the fields.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Role {
    /// Distance 0; no flux for the potential.
    Wall,
    /// Potential 0.
    Inlet,
    /// Potential 1.
    Outlet,
    /// No flux, not a wall.
    Symmetry,
    /// The face between two regions of a classification (no flux, not a
    /// wall).
    Cut,
}

/// A tetrahedral mesh with its boundary faces.
#[derive(Clone, Debug, Default)]
pub struct TetMesh {
    pub points: Vec<Vec3>,
    pub tets: Vec<[usize; 4]>,
    pub faces: Vec<([usize; 3], Role)>,
}

/// Linear finite elements on a [`TetMesh`].
#[derive(Clone, Debug)]
pub struct Fem {
    /// Gradients of the four hat functions of each tetrahedron.
    pub grads: Vec<[Vec3; 4]>,
    pub volume: Vec<f64>,
    pub stiffness: Csr,
    /// Lumped mass of each node (a quarter of each tetrahedron's volume).
    pub mass: Vec<f64>,
    /// Whether a node belongs to a tetrahedron.
    pub used: Vec<bool>,
    /// Mean edge length.
    pub h: f64,
}

/// Hat-function gradients and volume of a tetrahedron; `None` when flat.
pub fn tet_gradients(p: [Vec3; 4]) -> Option<([Vec3; 4], f64)> {
    let (e1, e2, e3) = (p[1] - p[0], p[2] - p[0], p[3] - p[0]);
    let det = e1.dot(e2.cross(e3));
    if det == 0.0 || !det.is_finite() {
        return None;
    }
    let g1 = e2.cross(e3) / det;
    let g2 = e3.cross(e1) / det;
    let g3 = e1.cross(e2) / det;
    let g0 = Vec3::ZERO - g1 - g2 - g3;
    Some(([g0, g1, g2, g3], det.abs() / 6.0))
}

impl Fem {
    pub fn new(mesh: &TetMesh) -> Result<Self, String> {
        let n = mesh.points.len();
        let mut grads = Vec::with_capacity(mesh.tets.len());
        let mut volume = Vec::with_capacity(mesh.tets.len());
        let mut mass = vec![0.0; n];
        let mut used = vec![false; n];
        let mut t = Vec::with_capacity(mesh.tets.len() * 16);
        let (mut edge_sum, mut edges) = (0.0, 0usize);
        for (k, tet) in mesh.tets.iter().enumerate() {
            if tet.iter().any(|&v| v >= n) {
                return Err(format!("tetrahedron {k} names a missing node"));
            }
            let p = tet.map(|v| mesh.points[v]);
            let (g, v) = tet_gradients(p).ok_or_else(|| format!("tetrahedron {k} is flat"))?;
            for a in 0..4 {
                mass[tet[a]] += v / 4.0;
                used[tet[a]] = true;
                for b in 0..4 {
                    t.push((tet[a], tet[b], v * g[a].dot(g[b])));
                }
                for b in a + 1..4 {
                    edge_sum += (p[a] - p[b]).length();
                    edges += 1;
                }
            }
            grads.push(g);
            volume.push(v);
        }
        // unused nodes get a unit diagonal so the matrix stays regular
        for i in 0..n {
            if !used[i] {
                t.push((i, i, 1.0));
            }
        }
        Ok(Self {
            grads,
            volume,
            stiffness: Csr::from_triplets(n, t),
            mass,
            used,
            h: if edges > 0 {
                edge_sum / edges as f64
            } else {
                0.0
            },
        })
    }

    /// The gradient of a node field in each tetrahedron.
    pub fn gradient(&self, mesh: &TetMesh, f: &[f64]) -> Vec<Vec3> {
        mesh.tets
            .iter()
            .zip(&self.grads)
            .map(|(t, g)| (0..4).fold(Vec3::ZERO, |s, a| s + g[a] * f[t[a]]))
            .collect()
    }
}

/// Nodes on boundary faces of the given roles.
pub fn nodes_with(mesh: &TetMesh, roles: &[Role]) -> Vec<bool> {
    let mut on = vec![false; mesh.points.len()];
    for (f, r) in &mesh.faces {
        if roles.contains(r) {
            for &v in f {
                on[v] = true;
            }
        }
    }
    on
}

const TOLERANCE: f64 = 1e-10;

fn iterations(n: usize) -> usize {
    (20 * n).clamp(2_000, 200_000)
}

/// The solves of one field, for the report.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FieldStats {
    pub iterations: usize,
    pub residual: f64,
    pub converged: bool,
}

impl FieldStats {
    fn add(&mut self, s: SolveStats) {
        self.iterations += s.iterations;
        self.residual = self.residual.max(s.residual);
        self.converged = s.converged && (self.converged || self.iterations == s.iterations);
    }
}

/// Heat-method distance to the wall faces with time step `(m h)^2`.
pub fn heat_distance(mesh: &TetMesh, fem: &Fem, m: f64) -> Result<(Vec<f64>, FieldStats), String> {
    let n = mesh.points.len();
    let wall = nodes_with(mesh, &[Role::Wall]);
    if !wall.iter().any(|&w| w) {
        return Err("the domain has no wall face".into());
    }
    let t = (m * fem.h).powi(2);
    let a = fem.stiffness.scaled_plus_diagonal(t, &fem.mass);
    let fixed: Vec<Option<f64>> = (0..n)
        .map(|i| {
            if !fem.used[i] {
                Some(0.0)
            } else if wall[i] {
                Some(1.0)
            } else {
                None
            }
        })
        .collect();
    let mut stats = FieldStats::default();
    let (u, s) = sparse::solve(&a, &vec![0.0; n], &fixed, None, TOLERANCE, iterations(n));
    stats.add(s);
    let grad = fem.gradient(mesh, &u);
    let mut b = vec![0.0; n];
    for (k, tet) in mesh.tets.iter().enumerate() {
        let g = grad[k];
        let len = g.length();
        if len.is_nan() || len <= 0.0 {
            continue;
        }
        let x = g * (-1.0 / len);
        for a in 0..4 {
            b[tet[a]] += fem.volume[k] * x.dot(fem.grads[k][a]);
        }
    }
    let fixed: Vec<Option<f64>> = (0..n)
        .map(|i| (!fem.used[i] || wall[i]).then_some(0.0))
        .collect();
    let (d, s) = sparse::solve(&fem.stiffness, &b, &fixed, None, TOLERANCE, iterations(n));
    stats.add(s);
    Ok((d.into_iter().map(|v| v.max(0.0)).collect(), stats))
}

/// The potential 0 on inlet faces and 1 on outlet faces.
pub fn potential(mesh: &TetMesh, fem: &Fem) -> Result<(Vec<f64>, FieldStats), String> {
    let n = mesh.points.len();
    let inlet = nodes_with(mesh, &[Role::Inlet]);
    let outlet = nodes_with(mesh, &[Role::Outlet]);
    if !inlet.iter().any(|&v| v) || !outlet.iter().any(|&v| v) {
        return Err("the potential needs an inlet face and an outlet face".into());
    }
    let fixed: Vec<Option<f64>> = (0..n)
        .map(|i| {
            if !fem.used[i] || inlet[i] {
                Some(0.0)
            } else if outlet[i] {
                Some(1.0)
            } else {
                None
            }
        })
        .collect();
    let mut stats = FieldStats::default();
    let (p, s) = sparse::solve(
        &fem.stiffness,
        &vec![0.0; n],
        &fixed,
        None,
        TOLERANCE,
        iterations(n),
    );
    stats.add(s);
    Ok((p, stats))
}

/// The Fiedler vector by inverse iteration with no flux everywhere,
/// scaled to [0, 1].
pub fn fiedler(mesh: &TetMesh, fem: &Fem, rounds: usize) -> (Vec<f64>, FieldStats) {
    let n = mesh.points.len();
    let total: f64 = fem.mass.iter().sum();
    // start: the coordinate along the longest bounding-box axis
    let (mut lo, mut hi) = (
        Vec3::new(f64::MAX, f64::MAX, f64::MAX),
        Vec3::new(f64::MIN, f64::MIN, f64::MIN),
    );
    for (i, p) in mesh.points.iter().enumerate() {
        if fem.used[i] {
            lo = Vec3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Vec3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
    }
    let e = hi - lo;
    let axis = if e.x >= e.y && e.x >= e.z {
        0
    } else if e.y >= e.z {
        1
    } else {
        2
    };
    let coord = |p: Vec3| [p.x, p.y, p.z][axis];
    let mut v: Vec<f64> = mesh.points.iter().map(|&p| coord(p)).collect();
    let fixed: Vec<Option<f64>> = (0..n).map(|i| (!fem.used[i]).then_some(0.0)).collect();
    // a tiny shift keeps the singular Neumann matrix definite
    let k_max = fem.stiffness.diagonal().iter().cloned().fold(0.0, f64::max);
    let m_max = fem.mass.iter().cloned().fold(f64::MIN_POSITIVE, f64::max);
    let shift: Vec<f64> = fem.mass.iter().map(|m| 1e-9 * k_max * m / m_max).collect();
    let a = fem.stiffness.scaled_plus_diagonal(1.0, &shift);
    let mut stats = FieldStats::default();
    let normalise = |v: &mut Vec<f64>| {
        let mean = (0..n).map(|i| fem.mass[i] * v[i]).sum::<f64>() / total;
        for i in 0..n {
            v[i] = if fem.used[i] { v[i] - mean } else { 0.0 };
        }
        let norm = (0..n)
            .map(|i| fem.mass[i] * v[i] * v[i])
            .sum::<f64>()
            .sqrt();
        if norm > 0.0 {
            for x in v.iter_mut() {
                *x /= norm;
            }
        }
    };
    normalise(&mut v);
    for _ in 0..rounds {
        let b: Vec<f64> = (0..n).map(|i| fem.mass[i] * v[i]).collect();
        let (w, s) = sparse::solve(&a, &b, &fixed, Some(&v), 1e-8, iterations(n));
        stats.add(s);
        let mut w = w;
        normalise(&mut w);
        let change = (0..n)
            .map(|i| fem.mass[i] * (w[i] - v[i]).powi(2))
            .sum::<f64>()
            .sqrt();
        let flipped = (0..n)
            .map(|i| fem.mass[i] * (w[i] + v[i]).powi(2))
            .sum::<f64>()
            .sqrt();
        v = w;
        if change.min(flipped) < 1e-6 {
            break;
        }
    }
    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
    for i in 0..n {
        if fem.used[i] {
            lo = lo.min(v[i]);
            hi = hi.max(v[i]);
        }
    }
    let span = (hi - lo).max(f64::MIN_POSITIVE);
    (
        (0..n)
            .map(|i| if fem.used[i] { (v[i] - lo) / span } else { 0.0 })
            .collect(),
        stats,
    )
}

/// The tetrahedron across each face (face `a` is opposite vertex `a`);
/// `None` on the boundary.
pub fn neighbours(tets: &[[usize; 4]]) -> Vec<[Option<usize>; 4]> {
    let mut faces: HashMap<[usize; 3], (usize, usize)> = HashMap::with_capacity(tets.len() * 2);
    let mut out = vec![[None; 4]; tets.len()];
    for (k, t) in tets.iter().enumerate() {
        for a in 0..4 {
            let key = face_key(t, a);
            match faces.remove(&key) {
                Some((other, b)) => {
                    out[k][a] = Some(other);
                    out[other][b] = Some(k);
                }
                None => {
                    faces.insert(key, (k, a));
                }
            }
        }
    }
    out
}

/// The sorted nodes of face `a` (opposite vertex `a`) of a tetrahedron.
pub fn face_key(t: &[usize; 4], a: usize) -> [usize; 3] {
    let mut f = [t[(a + 1) % 4], t[(a + 2) % 4], t[(a + 3) % 4]];
    f.sort_unstable();
    f
}

/// The tetrahedra `keep` as a mesh of their own: nodes renumbered, the
/// boundary faces they touch kept with their roles, and every face shared
/// with a dropped tetrahedron a [`Role::Cut`]. Returns the mesh and the
/// original node of each new node.
pub fn submesh(
    mesh: &TetMesh,
    keep: &[bool],
    adjacent: &[[Option<usize>; 4]],
) -> (TetMesh, Vec<usize>) {
    let mut new_id = vec![usize::MAX; mesh.points.len()];
    let mut old = Vec::new();
    let mut tets = Vec::new();
    let mut faces = Vec::new();
    let boundary: HashMap<[usize; 3], Role> = mesh
        .faces
        .iter()
        .map(|(f, r)| {
            let mut k = *f;
            k.sort_unstable();
            (k, *r)
        })
        .collect();
    let mut map = |v: usize, old: &mut Vec<usize>| {
        if new_id[v] == usize::MAX {
            new_id[v] = old.len();
            old.push(v);
        }
        new_id[v]
    };
    for (k, t) in mesh.tets.iter().enumerate() {
        if !keep[k] {
            continue;
        }
        tets.push(t.map(|v| map(v, &mut old)));
        for a in 0..4 {
            let key = face_key(t, a);
            let role = match adjacent[k][a] {
                Some(o) if keep[o] => continue,
                Some(_) => Role::Cut,
                None => boundary.get(&key).copied().unwrap_or(Role::Wall),
            };
            faces.push((key.map(|v| map(v, &mut old)), role));
        }
    }
    let points = old.iter().map(|&v| mesh.points[v]).collect();
    (
        TetMesh {
            points,
            tets,
            faces,
        },
        old,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A box [0, l] x [0, 1] x [0, 1] cut into 6 tetrahedra per cube cell.
    pub(crate) fn box_mesh(l: usize, n: usize) -> TetMesh {
        let (nx, ny, nz) = (l * n, n, n);
        let h = 1.0 / n as f64;
        let id = |i: usize, j: usize, k: usize| (i * (ny + 1) + j) * (nz + 1) + k;
        let mut points = Vec::new();
        for i in 0..=nx {
            for j in 0..=ny {
                for k in 0..=nz {
                    points.push(Vec3::new(i as f64 * h, j as f64 * h, k as f64 * h));
                }
            }
        }
        let mut tets = Vec::new();
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let c = [
                        id(i, j, k),
                        id(i + 1, j, k),
                        id(i + 1, j + 1, k),
                        id(i, j + 1, k),
                        id(i, j, k + 1),
                        id(i + 1, j, k + 1),
                        id(i + 1, j + 1, k + 1),
                        id(i, j + 1, k + 1),
                    ];
                    for t in [
                        [0, 1, 2, 6],
                        [0, 2, 3, 6],
                        [0, 3, 7, 6],
                        [0, 7, 4, 6],
                        [0, 4, 5, 6],
                        [0, 5, 1, 6],
                    ] {
                        tets.push(t.map(|a| c[a]));
                    }
                }
            }
        }
        let adjacent = neighbours(&tets);
        let mut faces = Vec::new();
        for (k, t) in tets.iter().enumerate() {
            for a in 0..4 {
                if adjacent[k][a].is_none() {
                    let f = face_key(t, a);
                    let xs: Vec<f64> = f.iter().map(|&v| points[v].x).collect();
                    let role = if xs.iter().all(|&x| x == 0.0) {
                        Role::Inlet
                    } else if xs.iter().all(|&x| (x - l as f64).abs() < 1e-12) {
                        Role::Outlet
                    } else {
                        Role::Wall
                    };
                    faces.push((f, role));
                }
            }
        }
        TetMesh {
            points,
            tets,
            faces,
        }
    }

    /// Keep: in a straight duct the potential is the axial coordinate over
    /// the length; a wrong stiffness sign or Dirichlet set breaks it, and
    /// the classification would only show it as a vaguer wrong region.
    #[test]
    fn potential_of_a_straight_duct_is_linear() {
        let m = box_mesh(4, 4);
        let fem = Fem::new(&m).unwrap();
        let (p, s) = potential(&m, &fem).unwrap();
        assert!(s.converged);
        for (i, q) in m.points.iter().enumerate() {
            assert!((p[i] - q.x / 4.0).abs() < 1e-8);
        }
    }

    /// Keep: the heat distance at the centre of a square duct of side 1 is
    /// about 0.5 and grows away from the wall; a sign error in the
    /// divergence makes it negative and clamps it to 0.
    #[test]
    fn heat_distance_peaks_on_the_axis() {
        let m = box_mesh(2, 8);
        let fem = Fem::new(&m).unwrap();
        let (d, s) = heat_distance(&m, &fem, 3.0).unwrap();
        assert!(s.converged);
        let centre = m
            .points
            .iter()
            .position(|p| {
                (p.x - 1.0).abs() < 1e-9 && (p.y - 0.5).abs() < 1e-9 && (p.z - 0.5).abs() < 1e-9
            })
            .unwrap();
        assert!(
            (d[centre] - 0.5).abs() < 0.12,
            "centre distance {}",
            d[centre]
        );
        let quarter = m
            .points
            .iter()
            .position(|p| {
                (p.x - 1.0).abs() < 1e-9 && (p.y - 0.25).abs() < 1e-9 && (p.z - 0.5).abs() < 1e-9
            })
            .unwrap();
        assert!(d[quarter] < d[centre] && d[quarter] > 0.1);
    }

    /// Keep: the Fiedler vector of a long box runs along its length; a
    /// missing deflation of the constant returns a flat field instead.
    #[test]
    fn fiedler_follows_the_long_axis() {
        let m = box_mesh(4, 2);
        let fem = Fem::new(&m).unwrap();
        let (f, _) = fiedler(&m, &fem, 40);
        let grad = fem.gradient(&m, &f);
        let mean = grad.iter().fold(Vec3::ZERO, |s, &g| s + g);
        assert!(mean.x.abs() > 10.0 * mean.y.abs().max(mean.z.abs()));
    }
}
