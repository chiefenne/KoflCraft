//! Tri-center building blocks (BS-03, MB1h): a triangle as three
//! quadrilaterals around a centre of valence 3, extruded as a prism of three
//! hexahedra, and a tetrahedron as four hexahedra around its body centre.
//!
//! **Counts.** Each corner i carries a count n_i. The piece at corner i has
//! n_j cells along its side toward corner j, so the edge (i, j) of the
//! triangle or tetrahedron gets n_i + n_j cells. Because the counts sit on
//! the corners, two triangles or tetrahedra that share an edge or a face
//! and use the same corner counts meet node to node, and the edge and face
//! points ([`TriCenter`]) are made once.
//!
//! **Shapes.** Edge points are the edge midpoints, face points the face
//! centroids, the centre the centroid of the corners. In a tetrahedron the
//! hexahedron at corner i spans i, the midpoints of its three edges, the
//! centroids of its three faces and the body centre; every corner angle of
//! a regular tetrahedron's hexahedra stays between 60° and 125° or so,
//! where one hexahedron with collapsed edges would have 0°.
//!
//! Sources: the subdivision of a simplex into hexahedra by its barycentric
//! dual (the "tri-center" map of a triangle, the hex split of a tetrahedron
//! in hex-dominant meshing), see `research/block-triangle-split.md`.

use std::collections::BTreeMap;

use mesh_core::block::{
    BLOCK_FACES, BlockSpec, BlockTopology, PatchSpec, face_corner_slot, face_key,
};
use mesh_core::geometry::Vec3;

use super::templates::{Template, counts, positive};

/// Shared edge and face points of a set of tri-center pieces, appended to
/// a vertex list once per edge and face.
pub struct TriCenter<'a> {
    pub vertices: &'a mut Vec<Vec3>,
    edges: BTreeMap<(usize, usize), usize>,
    faces: BTreeMap<[usize; 3], usize>,
}

impl<'a> TriCenter<'a> {
    pub fn new(vertices: &'a mut Vec<Vec3>) -> Self {
        Self {
            vertices,
            edges: BTreeMap::new(),
            faces: BTreeMap::new(),
        }
    }

    fn push(&mut self, p: Vec3) -> usize {
        self.vertices.push(p);
        self.vertices.len() - 1
    }

    /// The midpoint of the edge (a, b).
    pub fn edge(&mut self, a: usize, b: usize) -> usize {
        let k = (a.min(b), a.max(b));
        if let Some(&v) = self.edges.get(&k) {
            return v;
        }
        let p = (self.vertices[a] + self.vertices[b]) * 0.5;
        let v = self.push(p);
        self.edges.insert(k, v);
        v
    }

    /// The centroid of the triangle (a, b, c).
    pub fn face(&mut self, a: usize, b: usize, c: usize) -> usize {
        let mut k = [a, b, c];
        k.sort_unstable();
        if let Some(&v) = self.faces.get(&k) {
            return v;
        }
        let p = (self.vertices[a] + self.vertices[b] + self.vertices[c]) / 3.0;
        let v = self.push(p);
        self.faces.insert(k, v);
        v
    }

    /// The three hexahedra of the prism with the triangle `bottom` and the
    /// triangle `top` above it (corner i of `top` over corner i of
    /// `bottom`), `counts` per corner and `layers` cells between the two.
    /// The bottom corners may run either way round.
    pub fn prism(
        &mut self,
        bottom: [usize; 3],
        top: [usize; 3],
        counts: [usize; 3],
        layers: usize,
        domain: &str,
    ) -> Result<[BlockSpec; 3], String> {
        let (pb, pt) = (
            bottom.map(|v| self.vertices[v]),
            top.map(|v| self.vertices[v]),
        );
        let up = (pt[0] + pt[1] + pt[2]) / 3.0 - (pb[0] + pb[1] + pb[2]) / 3.0;
        let turn = (pb[1] - pb[0]).cross(pb[2] - pb[0]).dot(up);
        if !(turn.is_finite() && turn != 0.0) {
            return Err("the prism's triangle is flat or its top does not lie above it".into());
        }
        // run the corners so that bottom → top is the hexahedra's third axis
        let order = if turn > 0.0 { [0, 1, 2] } else { [0, 2, 1] };
        let (b, t, n) = (
            order.map(|i| bottom[i]),
            order.map(|i| top[i]),
            order.map(|i| counts[i]),
        );
        let ob = self.face(b[0], b[1], b[2]);
        let ot = self.face(t[0], t[1], t[2]);
        let mut out = Vec::with_capacity(3);
        for i in 0..3 {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            let (mb_ij, mb_ik) = (self.edge(b[i], b[j]), self.edge(b[i], b[k]));
            let (mt_ij, mt_ik) = (self.edge(t[i], t[j]), self.edge(t[i], t[k]));
            out.push(BlockSpec {
                vertices: [b[i], mb_ij, ob, mb_ik, t[i], mt_ij, ot, mt_ik],
                cells: [n[j], n[k], layers],
                domain: domain.into(),
            });
        }
        out.try_into()
            .map_err(|_| "three hexahedra per prism".to_string())
    }

    /// The four hexahedra of the tetrahedron `corners`, `counts` per corner.
    /// The corners may come in either orientation.
    pub fn tetrahedron(
        &mut self,
        corners: [usize; 4],
        counts: [usize; 4],
        domain: &str,
    ) -> Result<[BlockSpec; 4], String> {
        let c = corners;
        let p = c.map(|v| self.vertices[v]);
        let six_volume = (p[1] - p[0]).cross(p[2] - p[0]).dot(p[3] - p[0]);
        if !(six_volume.is_finite() && six_volume != 0.0) {
            return Err("the tetrahedron is flat".into());
        }
        let centre = self.push((p[0] + p[1] + p[2] + p[3]) / 4.0);
        let mut out = Vec::with_capacity(4);
        for i in 0..4 {
            let mut o: Vec<usize> = (0..4).filter(|&x| x != i).collect();
            // (j, k, l) such that i → j, i → k, i → l is right-handed
            let turn = |o: &[usize]| (p[o[0]] - p[i]).cross(p[o[1]] - p[i]).dot(p[o[2]] - p[i]);
            if turn(&o) < 0.0 {
                o.swap(1, 2);
            }
            let (j, k, l) = (o[0], o[1], o[2]);
            let (ci, cj, ck, cl) = (c[i], c[j], c[k], c[l]);
            let (m_ij, m_ik, m_il) = (self.edge(ci, cj), self.edge(ci, ck), self.edge(ci, cl));
            let (f_ijk, f_ijl, f_ikl) = (
                self.face(ci, cj, ck),
                self.face(ci, cj, cl),
                self.face(ci, ck, cl),
            );
            out.push(BlockSpec {
                vertices: [ci, m_ij, f_ijk, m_ik, m_il, f_ijl, centre, f_ikl],
                cells: [counts[j], counts[k], counts[l]],
                domain: domain.into(),
            });
        }
        out.try_into()
            .map_err(|_| "four hexahedra per tetrahedron".to_string())
    }
}

/// The boundary faces of `blocks`: faces used by one block only.
fn boundary_faces(blocks: &[BlockSpec]) -> Vec<[usize; 4]> {
    let mut seen: BTreeMap<[usize; 4], (usize, [usize; 4])> = BTreeMap::new();
    for b in blocks {
        for face in BLOCK_FACES {
            let f = [(0, 0), (1, 0), (1, 1), (0, 1)]
                .map(|(ca, cb)| b.vertices[face_corner_slot(face, ca, cb)]);
            seen.entry(face_key(f))
                .and_modify(|e| e.0 += 1)
                .or_insert((1, f));
        }
    }
    seen.into_values()
        .filter(|(n, _)| *n == 1)
        .map(|(_, f)| f)
        .collect()
}

fn corner_counts<const N: usize>(cells: [usize; N]) -> Result<[usize; N], String> {
    counts(&cells)?;
    Ok(cells)
}

/// A triangle extruded along +z by `span` as three hexahedra around its
/// centroid (the `triangle` template): `cells` per corner, sides `walls`,
/// ends `symmetry_low` and `symmetry_high`.
pub fn triangle(
    corners: [[f64; 2]; 3],
    cells: [usize; 3],
    span: f64,
    span_cells: usize,
) -> Result<Template, String> {
    let span = positive(span, "span")?;
    let n = corner_counts(cells)?;
    counts(&[span_cells])?;
    let mut vertices: Vec<Vec3> = Vec::new();
    for z in [0.0, span] {
        vertices.extend(corners.iter().map(|p| Vec3::new(p[0], p[1], z)));
    }
    let area = 0.5
        * ((corners[1][0] - corners[0][0]) * (corners[2][1] - corners[0][1])
            - (corners[2][0] - corners[0][0]) * (corners[1][1] - corners[0][1]))
            .abs();
    if !(area.is_finite() && area > 0.0) {
        return Err("the triangle's corners must not lie on a line".into());
    }
    let blocks = TriCenter::new(&mut vertices)
        .prism([0, 1, 2], [3, 4, 5], n, span_cells, "fluid")?
        .to_vec();
    let (mut low, mut high, mut sides) = (Vec::new(), Vec::new(), Vec::new());
    for f in boundary_faces(&blocks) {
        let z: Vec<f64> = f.iter().map(|&v| vertices[v].z).collect();
        let tol = 1e-9 * span;
        if z.iter().all(|&z| z.abs() < tol) {
            low.push(f);
        } else if z.iter().all(|&z| (z - span).abs() < tol) {
            high.push(f);
        } else {
            sides.push(f);
        }
    }
    let patch = |name: &str, faces: Vec<[usize; 4]>| PatchSpec {
        name: name.into(),
        faces,
    };
    Ok(Template {
        topology: BlockTopology {
            junctions: Vec::new(),
            vertices,
            edges: Vec::new(),
            blocks,
            patches: vec![
                patch("walls", sides),
                patch("symmetry_low", low),
                patch("symmetry_high", high),
            ],
            surfaces: Vec::new(),
        },
        analytic_volume: area * span,
        layout: None,
        unoptimised: None,
    })
}

/// A tetrahedron as four hexahedra around its centroid (the `tetrahedron`
/// template): `cells` per corner, all faces `walls`.
pub fn tetrahedron(corners: [[f64; 3]; 4], cells: [usize; 4]) -> Result<Template, String> {
    let n = corner_counts(cells)?;
    let mut vertices: Vec<Vec3> = corners
        .iter()
        .map(|p| Vec3::new(p[0], p[1], p[2]))
        .collect();
    let volume = (vertices[1] - vertices[0])
        .cross(vertices[2] - vertices[0])
        .dot(vertices[3] - vertices[0])
        .abs()
        / 6.0;
    if !(volume.is_finite() && volume > 0.0) {
        return Err("the tetrahedron's corners must not lie in a plane".into());
    }
    let blocks = TriCenter::new(&mut vertices)
        .tetrahedron([0, 1, 2, 3], n, "fluid")?
        .to_vec();
    let faces = boundary_faces(&blocks);
    Ok(Template {
        topology: BlockTopology {
            junctions: Vec::new(),
            vertices,
            edges: Vec::new(),
            blocks,
            patches: vec![PatchSpec {
                name: "walls".into(),
                faces,
            }],
            surfaces: Vec::new(),
        },
        analytic_volume: volume,
        layout: None,
        unoptimised: None,
    })
}

/// The corners of a regular tetrahedron with unit edges, one face on z = 0.
pub fn regular_tetrahedron() -> [[f64; 3]; 4] {
    let h = (2.0f64 / 3.0).sqrt();
    let r = 1.0 / 3.0f64.sqrt();
    [
        [r, 0.0, 0.0],
        [-0.5 * r, 0.5, 0.0],
        [-0.5 * r, -0.5, 0.0],
        [0.0, 0.0, h],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::block::{build, check};

    fn min_jacobian(t: &BlockTopology) -> f64 {
        build(t)
            .unwrap()
            .quality()
            .iter()
            .map(|q| q.min_jacobian)
            .fold(f64::INFINITY, f64::min)
    }

    /// Invariant (BS-03, BS-05, MB1h): a tetrahedron becomes four
    /// hexahedra whose edge (i, j) has n_i + n_j cells; the topology
    /// checks, no cell is inverted, the cells fill the tetrahedron's volume
    /// and the patch holds the 12 outer faces (three per side). A wrong
    /// corner order (a left-handed hexahedron) or a count on the wrong axis
    /// breaks it; the general block tests never build a tetrahedron.
    #[test]
    fn a_tetrahedron_becomes_four_hexahedra() {
        // corners given left-handed on purpose: the split must fix the turn
        let mut c = regular_tetrahedron();
        c.swap(1, 2);
        let t = tetrahedron(c, [3, 4, 5, 6]).unwrap();
        assert_eq!(t.topology.blocks.len(), 4);
        assert_eq!(t.topology.vertices.len(), 4 + 6 + 4 + 1);
        let report = check(&t.topology).unwrap();
        // the hexahedron at corner i holds n_j * n_k * n_l cells
        let n = [3usize, 4, 5, 6];
        let expected: usize = (0..4)
            .map(|i| (0..4).filter(|&j| j != i).map(|j| n[j]).product::<usize>())
            .sum();
        assert_eq!(report.cells, expected);
        assert!(min_jacobian(&t.topology) > 0.0);
        let grid = build(&t.topology).unwrap();
        let v: f64 = grid.volumes().iter().sum();
        assert!(
            (v - t.analytic_volume).abs() < 1e-12 * t.analytic_volume,
            "{v} {}",
            t.analytic_volume
        );
        assert_eq!(t.topology.patches[0].faces.len(), 12);
        // the edge between corners 0 and 1 runs through two hexahedra with
        // n_1 + n_0 cells
        let edge_cells: usize = t
            .topology
            .blocks
            .iter()
            .filter(|b| b.vertices[0] == 0 || b.vertices[0] == 1)
            .map(|b| {
                let other = if b.vertices[0] == 0 { 1 } else { 0 };
                let m = t
                    .topology
                    .vertices
                    .iter()
                    .position(|&p| {
                        (p - (t.topology.vertices[0] + t.topology.vertices[1]) * 0.5).length()
                            < 1e-12
                    })
                    .unwrap();
                let axis = [1usize, 3, 4]
                    .iter()
                    .position(|&s| b.vertices[s] == m)
                    .unwrap();
                assert_eq!(
                    b.cells[axis], n[other],
                    "the half next to corner {other}'s cut"
                );
                b.cells[axis]
            })
            .sum();
        assert_eq!(edge_cells, n[0] + n[1]);
    }

    /// Invariant (BS-03, MB1h): two tetrahedra on a shared face with the
    /// same corner counts meet node to node: the edge and face points are
    /// made once and the joined topology checks with the shared faces
    /// matched. Making the points per tetrahedron (duplicates) or counts
    /// per piece instead of per corner breaks it.
    #[test]
    fn two_tetrahedra_share_their_face() {
        let mut vertices = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.2, 0.3, 1.0),
            Vec3::new(0.3, 0.2, -1.0),
        ];
        let n = [3usize, 4, 5, 2, 6];
        let mut tc = TriCenter::new(&mut vertices);
        let mut blocks = tc
            .tetrahedron([0, 1, 2, 3], [n[0], n[1], n[2], n[3]], "fluid")
            .unwrap()
            .to_vec();
        blocks.extend(
            tc.tetrahedron([0, 1, 2, 4], [n[0], n[1], n[2], n[4]], "fluid")
                .unwrap(),
        );
        assert_eq!(vertices.len(), 5 + 9 + 7 + 2, "shared points made once");
        let faces = boundary_faces(&blocks);
        assert_eq!(faces.len(), 18, "three outer faces per outer side");
        let t = BlockTopology {
            junctions: Vec::new(),
            vertices,
            edges: Vec::new(),
            blocks,
            patches: vec![PatchSpec {
                name: "walls".into(),
                faces,
            }],
            surfaces: Vec::new(),
        };
        let report = check(&t).unwrap();
        assert_eq!(report.block_faces_shared, 3 * 4 + 3);
        assert!(min_jacobian(&t) > 0.0);
    }

    /// Invariant (BS-03, MB1h): the triangle template is a prism of three
    /// hexahedra whose sides have n_i + n_j cells, with positive Jacobians
    /// and the prism's volume, for either turn of the corners. A wrong turn
    /// fix or a count on the wrong side breaks it; the engine's split of an
    /// existing block (block_triangle) builds its own pieces.
    #[test]
    fn a_triangle_becomes_a_prism_of_three_hexahedra() {
        for corners in [
            [[0.0, 0.0], [2.0, 0.0], [0.7, 1.5]],
            [[0.0, 0.0], [0.7, 1.5], [2.0, 0.0]],
        ] {
            let t = triangle(corners, [4, 6, 3], 0.1, 2).unwrap();
            let report = check(&t.topology).unwrap();
            assert_eq!(report.cells, 2 * (6 * 3 + 4 * 3 + 4 * 6));
            assert!(min_jacobian(&t.topology) > 0.0);
            let v: f64 = build(&t.topology).unwrap().volumes().iter().sum();
            assert!((v - t.analytic_volume).abs() < 1e-12 * t.analytic_volume);
            assert_eq!(t.topology.patches[0].faces.len(), 6);
            assert_eq!(t.topology.patches[1].faces.len(), 3);
        }
    }
}
