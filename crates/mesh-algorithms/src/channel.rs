//! Analytic rectangular channel with wall-normal prism stacks and open planar ports.

use std::collections::BTreeMap;

use mesh_core::geometry::{AxisAlignedBox, BoxSide, Vec3};
use mesh_core::mesh::{FaceAdjacency, Mesh};

/// Positive reference order: the first triangle faces the second triangle.
/// The core is a polyhedron with two quad ends and eight triangular side faces.
pub enum ChannelCell {
    Prism { nodes: [usize; 6], layer: usize },
    Core { nodes: [usize; 8] },
}

impl ChannelCell {
    pub fn layer(&self) -> usize {
        match self {
            Self::Prism { layer, .. } => *layer,
            Self::Core { .. } => 0,
        }
    }

    fn faces(&self) -> Vec<Vec<usize>> {
        match *self {
            Self::Prism {
                nodes: [a, b, c, d, e, f],
                ..
            } => vec![
                vec![a, c, b],
                vec![d, e, f],
                vec![a, b, e, d],
                vec![b, c, f, e],
                vec![c, a, d, f],
            ],
            Self::Core {
                nodes: [a, b, c, d, e, f, g, h],
            } => vec![
                vec![a, e, h, d],
                vec![b, c, g, f],
                vec![a, b, f],
                vec![a, f, e],
                vec![d, h, g],
                vec![d, g, c],
                vec![a, d, c],
                vec![a, c, b],
                vec![e, f, g],
                vec![e, g, h],
            ],
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum ChannelPatch {
    Inlet = 1,
    Outlet = 2,
    Walls = 3,
}

impl ChannelPatch {
    pub const ALL: [Self; 3] = [Self::Inlet, Self::Outlet, Self::Walls];

    pub fn from_side(side: BoxSide) -> Self {
        match side {
            BoxSide::XMin => Self::Inlet,
            BoxSide::XMax => Self::Outlet,
            _ => Self::Walls,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Inlet => "inlet",
            Self::Outlet => "outlet",
            Self::Walls => "walls",
        }
    }
}

pub struct ChannelMesh {
    bounds: AxisAlignedBox,
    mesh: Mesh,
    cells: Vec<ChannelCell>,
    insets: Vec<f64>,
}

impl ChannelMesh {
    pub fn mesh(&self) -> &Mesh {
        &self.mesh
    }
    pub fn bounds(&self) -> AxisAlignedBox {
        self.bounds
    }
    pub fn cells(&self) -> &[ChannelCell] {
        &self.cells
    }
    /// Wall distances, including zero and the layer-to-core interface.
    pub fn insets(&self) -> &[f64] {
        &self.insets
    }

    pub fn boundary(&self) -> impl Iterator<Item = (usize, usize, ChannelPatch)> + '_ {
        (0..self.mesh.face_count()).filter_map(|face| match self.mesh.face_adjacency(face)? {
            FaceAdjacency::Boundary { owner, patch } => {
                Some((face, owner, ChannelPatch::from_side(patch)))
            }
            _ => None,
        })
    }
}

/// Build a channel along +X. The transverse grid uses the same wall distances
/// on all four sides. Corner tiles use axial prisms to join the two normal stacks.
pub fn straight_channel(
    size: Vec3,
    axial_cells: usize,
    layers: usize,
    first_thickness: f64,
    growth: f64,
) -> Result<ChannelMesh, String> {
    let bounds =
        AxisAlignedBox::new(Vec3::ZERO, size).ok_or("dimensions must be finite and positive")?;
    if axial_cells == 0
        || layers == 0
        || !first_thickness.is_finite()
        || first_thickness <= 0.0
        || !growth.is_finite()
        || growth <= 0.0
    {
        return Err(
            "axial_cells, layers, first_thickness and growth must be positive and finite".into(),
        );
    }
    let n = layers
        .checked_mul(2)
        .and_then(|n| n.checked_add(2))
        .ok_or("grid size overflow")?;
    let stations = axial_cells.checked_add(1).ok_or("grid size overflow")?;
    let point_count = n
        .checked_mul(n)
        .and_then(|v| v.checked_mul(stations))
        .ok_or("grid size overflow")?;
    let cell_count = (n - 1)
        .checked_pow(2)
        .and_then(|v| v.checked_mul(2))
        .and_then(|v| (v - 1).checked_mul(axial_cells))
        .ok_or("grid size overflow")?;
    let mut insets = vec![0.0];
    let mut thickness = first_thickness;
    for _ in 0..layers {
        let previous = *insets.last().unwrap();
        let next = previous + thickness;
        if next <= previous || next >= size.y.min(size.z) * 0.5 {
            return Err("layer stack collapses or leaves no positive core".into());
        }
        insets.push(next);
        thickness *= growth;
    }
    let coordinates = |extent: f64| -> Result<Vec<f64>, String> {
        let v: Vec<_> = insets
            .iter()
            .copied()
            .chain(insets.iter().rev().map(|d| extent - d))
            .collect();
        if v.windows(2).any(|s| s[1] <= s[0]) {
            return Err("layer spacing is below coordinate precision".into());
        }
        Ok(v)
    };
    let ys = coordinates(size.y)?;
    let zs = coordinates(size.z)?;
    let mut points = Vec::new();
    points
        .try_reserve_exact(point_count)
        .map_err(|e| e.to_string())?;
    for i in 0..stations {
        let x = size.x * (i as f64 / axial_cells as f64);
        for &y in &ys {
            for &z in &zs {
                points.push(Vec3::new(x, y, z));
            }
        }
    }

    let mut cells = Vec::new();
    cells
        .try_reserve_exact(cell_count)
        .map_err(|e| e.to_string())?;
    let vertex = |i, j, k| (i * n + j) * n + k;
    for i in 0..axial_cells {
        for j in 0..n - 1 {
            for k in 0..n - 1 {
                let nodes = [
                    vertex(i, j, k),
                    vertex(i + 1, j, k),
                    vertex(i + 1, j + 1, k),
                    vertex(i, j + 1, k),
                    vertex(i, j, k + 1),
                    vertex(i + 1, j, k + 1),
                    vertex(i + 1, j + 1, k + 1),
                    vertex(i, j + 1, k + 1),
                ];
                let (base, top) = match (j == layers, k == layers) {
                    (true, true) => {
                        cells.push(ChannelCell::Core { nodes });
                        continue;
                    }
                    (false, true) => ([0, 4, 5, 1], [3, 7, 6, 2]), // normal to a Y wall
                    (true, false) => ([0, 1, 2, 3], [4, 5, 6, 7]), // normal to a Z wall
                    (false, false) => ([0, 3, 7, 4], [1, 2, 6, 5]), // axial corner transition
                };
                let layer = 1 + j.min(n - 2 - j).min(k.min(n - 2 - k));
                for [a, b, c] in [[0, 1, 2], [0, 2, 3]] {
                    cells.push(ChannelCell::Prism {
                        nodes: [base[a], base[b], base[c], top[a], top[b], top[c]]
                            .map(|v| nodes[v]),
                        layer,
                    });
                }
            }
        }
    }
    let mesh = assemble(points, &cells, size)?;
    Ok(ChannelMesh {
        bounds,
        mesh,
        cells,
        insets,
    })
}

fn assemble(points: Vec<Vec3>, cells: &[ChannelCell], size: Vec3) -> Result<Mesh, String> {
    let mut faces = BTreeMap::new();
    for (cell, shape) in cells.iter().enumerate() {
        for ring in shape.faces() {
            let mut key = ring.clone();
            key.sort_unstable();
            use std::collections::btree_map::Entry;
            match faces.entry(key) {
                Entry::Vacant(entry) => {
                    entry.insert((ring, cell, None));
                }
                Entry::Occupied(mut entry) => {
                    let (first, _, other) = entry.get_mut();
                    let start = first.iter().position(|v| *v == ring[0]).unwrap();
                    if other.replace(cell).is_some()
                        || !ring
                            .iter()
                            .enumerate()
                            .all(|(i, id)| first[(start + first.len() - i) % first.len()] == *id)
                    {
                        return Err("inconsistent incidences at a channel interface".into());
                    }
                }
            }
        }
    }
    let mut ids = Vec::new();
    let mut offsets = vec![0];
    let mut adjacency = Vec::new();
    for (ring, owner, other) in faces.into_values() {
        adjacency.push(if let Some(neighbour) = other {
            FaceAdjacency::Interior { owner, neighbour }
        } else {
            let patch = BoxSide::ALL
                .into_iter()
                .find(|side| {
                    ring.iter().all(|&id| {
                        let p = points[id];
                        let axis = side.axis();
                        [p.x, p.y, p.z][axis]
                            == if side.is_max() {
                                [size.x, size.y, size.z][axis]
                            } else {
                                0.0
                            }
                    })
                })
                .ok_or("unpaired interior face in channel")?;
            FaceAdjacency::Boundary { owner, patch }
        });
        ids.extend(ring);
        offsets.push(ids.len());
    }
    Mesh::from_faces(points, offsets, ids, adjacency, cells.len()).map_err(|e| e.to_string())
}
