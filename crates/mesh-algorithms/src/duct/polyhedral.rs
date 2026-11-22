//! Fit a triangular-to-dual transition inside each validated star region.
//! The inset surface uses a median dual shared by the incident triangles.
use super::{Cell, DuctMesh, Element};
use mesh_core::geometry::Vec3;
use std::collections::{BTreeMap, HashMap};

const INSET_FRACTION: f64 = 0.75;

pub struct Polyhedron {
    pub nodes: Vec<usize>,
    pub faces: Vec<Vec<usize>>,
    /// Strict kernel point used for geometric integration.
    pub interior: Vec3,
    /// A dual cone can be exported to Gmsh as one pyramid per base quad.
    pub cone_apex: Option<usize>,
    pub(super) expected_volume: f64,
}

impl Polyhedron {
    pub(super) fn new(
        faces: Vec<Vec<usize>>,
        interior: Vec3,
        cone_apex: Option<usize>,
        expected_volume: f64,
    ) -> Self {
        let mut nodes: Vec<_> = faces.iter().flatten().copied().collect();
        nodes.sort_unstable();
        nodes.dedup();
        Self {
            nodes,
            faces,
            interior,
            cone_apex,
            expected_volume,
        }
    }
}

#[derive(Default)]
struct DualCell {
    faces: Vec<Vec<usize>>,
    volume: f64,
}

/// Each transition cell takes its source tetrahedron's index, so all existing
/// exterior face owners and layer IDs remain valid. Dual cells are appended.
pub(super) fn convert(mesh: &mut DuctMesh) -> Result<(), String> {
    let mut corners = HashMap::new();
    let mut midpoints = HashMap::new();
    let mut dual = BTreeMap::<(usize, usize, usize), DualCell>::new();
    let points = &mut mesh.points;
    for cell in &mut mesh.cells {
        let Element::Tet([kernel, a, b, c]) = cell.element else {
            continue;
        };
        let origin = points[kernel];
        let outer = [a, b, c];
        let inner = outer.map(|v| {
            *corners.entry((kernel, v)).or_insert_with(|| {
                let index = points.len();
                points.push(origin + (points[v] - origin) * INSET_FRACTION);
                index
            })
        });
        let mid: [usize; 3] = std::array::from_fn(|i| {
            let (u, v) = (inner[i], inner[(i + 1) % 3]);
            *midpoints.entry((u.min(v), u.max(v))).or_insert_with(|| {
                let index = points.len();
                points.push(points[u] + (points[v] - points[u]) * 0.5);
                index
            })
        });
        let centroid = points[a] + ((points[b] - points[a]) + (points[c] - points[a])) / 3.0;
        let center = points.len();
        points.push(origin + (centroid - origin) * INSET_FRACTION);
        let volume = (points[b] - points[a])
            .cross(points[c] - points[a])
            .dot(points[a] - origin)
            / 6.0;
        if !volume.is_finite() || volume <= 0.0 {
            return Err("invalid source tetrahedron for dual core".into());
        }
        let inner_volume = volume * INSET_FRACTION.powi(3);
        let mut faces = vec![outer.to_vec()];
        for i in 0..3 {
            let j = (i + 1) % 3;
            let previous = (i + 2) % 3;
            let quad = [inner[i], mid[i], center, mid[previous]];
            // Opposite sides of each interface use exactly the same four IDs.
            faces.push(quad.into_iter().rev().collect());
            let part = dual.entry((cell.region, kernel, inner[i])).or_default();
            part.faces.push(quad.to_vec());
            part.faces.push(vec![kernel, center, mid[i]]);
            part.faces.push(vec![kernel, mid[previous], center]);
            part.volume += inner_volume / 3.0;

            // A shared side strip has two outer endpoints and three inset
            // points. Choose its diagonal solely from the common endpoint IDs.
            let (u, v, iu, iv, m) = (outer[i], outer[j], inner[i], inner[j], mid[i]);
            if u < v {
                faces.push(vec![u, iu, m]);
                faces.push(vec![v, u, m, iv]);
            } else {
                faces.push(vec![v, u, iu, m]);
                faces.push(vec![v, m, iv]);
            }
        }
        let interior = origin + (centroid - origin) * ((1.0 + INSET_FRACTION) * 0.5);
        cell.element = Element::Polyhedron(Box::new(Polyhedron::new(
            faces,
            interior,
            None,
            volume - inner_volume,
        )));
    }
    for ((region, kernel, corner), part) in dual {
        let interior = points[kernel] + (points[corner] - points[kernel]) * 0.5;
        mesh.cells.push(Cell {
            region,
            layer: 0,
            element: Element::Polyhedron(Box::new(Polyhedron::new(
                part.faces,
                interior,
                Some(kernel),
                part.volume,
            ))),
        });
    }
    Ok(())
}
