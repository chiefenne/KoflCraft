//! Transverse Voronoi cells and wall-fitted prism layers for straight analytic domains.
mod section;

use mesh_core::geometry::Vec3;
use std::collections::BTreeMap;
use std::f64::consts::PI;

#[derive(Clone, Copy)]
pub enum Profile {
    BlindBlock { width: f64, height: f64 },
    Pipe { diameter: f64 },
}
impl Profile {
    pub fn blind(self) -> bool {
        matches!(self, Self::BlindBlock { .. })
    }
    pub fn area(self) -> f64 {
        match self {
            Self::BlindBlock { width, height } => width * height,
            Self::Pipe { diameter } => PI * (diameter / 2.).powi(2),
        }
    }
}

pub struct Tile {
    pub ring: Vec<usize>,
    pub layer: usize,
    pub area: f64,
    pub centroid: Vec3,
    pub neighbour_ratio: f64,
    pub non_orthogonality: f64,
    pub skewness: f64,
}
impl Tile {
    fn new(ring: Vec<usize>, layer: usize, points: &[Vec3]) -> Result<Self, String> {
        let origin = points[ring[0]];
        let mut twice_area = 0.;
        let mut moment = Vec3::ZERO;
        for pair in ring[1..].windows(2) {
            let u = points[pair[0]] - origin;
            let v = points[pair[1]] - origin;
            let det = u.cross(v).x;
            if !det.is_finite() || det <= 0. {
                return Err("degenerate or inverted section polygon".into());
            }
            twice_area += det;
            moment = moment + (u + v) * (det / 3.);
        }
        Ok(Self {
            ring,
            layer,
            area: twice_area / 2.,
            centroid: origin + moment / twice_area,
            neighbour_ratio: 1.,
            non_orthogonality: 0.,
            skewness: 0.,
        })
    }
}

pub struct Sweep {
    pub profile: Profile,
    pub points: Vec<Vec3>, // transverse stations have x = 0
    pub tiles: Vec<Tile>,
    pub outside: Vec<([usize; 2], usize)>,
    pub stations: Vec<f64>,
    pub axial_layers: Vec<usize>,
    pub axial_ratios: Vec<f64>,
    pub thicknesses: Vec<f64>,
    pub fitted_area: f64,
    pub rim_sides: usize,
}
impl Sweep {
    pub fn cells(&self) -> usize {
        self.tiles.len() * self.axial_layers.len()
    }
    pub fn layer(&self, station: usize, tile: usize) -> usize {
        let a = self.axial_layers[station];
        let b = self.tiles[tile].layer;
        match (a, b) {
            (0, _) => b,
            (_, 0) => a,
            _ => a.min(b),
        }
    }
    pub fn point_id(&self, station: usize, vertex: usize) -> usize {
        station * self.points.len() + vertex
    }
    pub fn volume(&self, station: usize, tile: usize) -> f64 {
        (self.stations[station + 1] - self.stations[station]) * self.tiles[tile].area
    }
    pub fn neighbour_ratio(&self, station: usize, tile: usize) -> f64 {
        self.axial_ratios[station].max(self.tiles[tile].neighbour_ratio)
    }
}

fn edges(ring: &[usize]) -> impl Iterator<Item = (usize, usize)> + '_ {
    ring.iter()
        .copied()
        .zip(ring.iter().copied().cycle().skip(1))
        .take(ring.len())
}

fn topology(tiles: &mut [Tile], points: &[Vec3]) -> Result<Vec<([usize; 2], usize)>, String> {
    let mut uses = BTreeMap::<[usize; 2], Vec<(usize, usize, usize)>>::new();
    for (i, t) in tiles.iter().enumerate() {
        for (a, b) in edges(&t.ring) {
            uses.entry([a.min(b), a.max(b)])
                .or_default()
                .push((i, a, b));
        }
    }
    let mut outer = Vec::new();
    for edge in uses.into_values() {
        match edge[..] {
            [(owner, a, b)] => outer.push(([a, b], owner)),
            [(i, a, b), (j, c, d)] if a == d && b == c => {
                let p = points[a];
                let q = points[b];
                let e = q - p;
                let normal = Vec3::new(0., e.z, -e.y) / e.length();
                let delta = tiles[j].centroid - tiles[i].centroid;
                let proj = delta.dot(normal);
                let to_face = (p + q) * 0.5 - tiles[i].centroid;
                let from = to_face.dot(normal);
                if !(from > 0. && proj > from) {
                    return Err("incorrect orientation at shared section edge".into());
                }
                let nonorth = delta.cross(normal).length().atan2(proj).to_degrees();
                let skew = (to_face - delta * (from / proj)).length() / delta.length();
                let ratio =
                    (tiles[i].area.max(tiles[j].area) / tiles[i].area.min(tiles[j].area)).cbrt();
                for k in [i, j] {
                    tiles[k].neighbour_ratio = tiles[k].neighbour_ratio.max(ratio);
                    tiles[k].non_orthogonality = tiles[k].non_orthogonality.max(nonorth);
                    tiles[k].skewness = tiles[k].skewness.max(skew);
                }
            }
            _ => return Err("nonmanifold or inconsistently wound section edge".into()),
        }
    }
    Ok(outer)
}

pub struct Sizing {
    pub core: f64,
    pub layers: usize,
    pub first: f64,
    pub growth: f64,
    pub ratio_limit: f64,
}

pub fn build(
    profile: Profile,
    length: f64,
    size: &Sizing,
    cell_budget: usize,
) -> Result<Sweep, String> {
    let Sizing {
        core: h,
        layers,
        first,
        growth,
        ratio_limit,
    } = *size;
    let (w, height) = match profile {
        Profile::BlindBlock { width, height } => (width, height),
        Profile::Pipe { diameter } => (diameter, diameter),
    };
    if ![w, height, length, h, first, growth, ratio_limit]
        .iter()
        .all(|x| x.is_finite() && *x > 0.)
        || layers == 0
        || growth < 1.
        || ratio_limit < 1.
    {
        return Err("positive dimensions and spacings, at least one layer, growth and size ratio >= 1 required".into());
    }
    if !(profile.area() * length).is_normal() {
        return Err("analytic volume outside representable range".into());
    }
    let mut thicknesses = Vec::new();
    let mut stack = 0.;
    let mut step = first;
    for _ in 0..layers {
        let next = stack + step;
        if next <= stack
            || next * 2. >= w.min(height)
            || (profile.blind() && next >= length)
            || thicknesses.len() >= cell_budget
        {
            return Err("layer stack collapses the core or exceeds the cell budget".into());
        }
        stack = next;
        thicknesses.push(step);
        step *= growth;
    }
    let core_length = if profile.blind() {
        length - stack
    } else {
        length
    };
    let axial = (core_length / h).ceil().max(1.) as usize;
    let stations_count = axial
        .checked_add(if profile.blind() { layers } else { 0 })
        .ok_or("axial count overflow")?;
    let section_budget = cell_budget
        .checked_div(stations_count)
        .ok_or("invalid axial count")?;
    if section_budget == 0 {
        return Err("requested sweep exceeds cell budget".into());
    }
    let section::Section {
        mut points,
        mut tiles,
        boundary: core_boundary,
    } = section::core(profile, h, stack, section_budget)?;
    let inner = topology(&mut tiles, &points)?;
    let count = inner
        .len()
        .checked_mul(layers)
        .and_then(|v| v.checked_add(tiles.len()))
        .ok_or("cell count overflow")?;
    if count > section_budget {
        return Err("requested sweep exceeds cell budget".into());
    }
    let mut previous = BTreeMap::new();
    for &([a, b], _) in &inner {
        previous.insert(a, a);
        previous.insert(b, b);
    }
    let core_perimeter: f64 = inner
        .iter()
        .map(|&([a, b], _)| (points[a] - points[b]).length())
        .sum();
    let expected_perimeter: f64 = core_boundary
        .iter()
        .zip(core_boundary.iter().cycle().skip(1))
        .map(|(&a, &b)| (a - b).length())
        .sum();
    if (core_perimeter / expected_perimeter - 1.).abs() > 1e-9 {
        return Err("unmatched interior section edge".into());
    }
    let rim_sides = core_boundary.len();
    let cos_half = (PI / rim_sides as f64).cos();
    let mut distance = 0.;
    for layer in (1..=layers).rev() {
        distance += thicknesses[layer - 1];
        let mut next = BTreeMap::new();
        for &vertex in previous.keys() {
            let q = points[vertex];
            let pos = match profile {
                Profile::Pipe { diameter } => {
                    q * (1. + distance / (diameter * 0.5 * cos_half - stack))
                }
                Profile::BlindBlock { width, height } => {
                    // Spread the tangential displacement along the entire side.
                    // Concentrating it in the corner tile creates an abrupt size jump.
                    Vec3::new(
                        0.,
                        q.y * (1. + distance / (0.5 * width - stack)),
                        q.z * (1. + distance / (0.5 * height - stack)),
                    )
                }
            };
            next.insert(vertex, points.len());
            points.push(pos);
        }
        for &([a, b], _) in &inner {
            tiles.push(Tile::new(
                vec![previous[&b], previous[&a], next[&a], next[&b]],
                layer,
                &points,
            )?);
        }
        previous = next;
    }
    let outside = topology(&mut tiles, &points)?;
    if outside.len() != inner.len() {
        return Err("unpaired wall layer edges".into());
    }
    let fitted_area = match profile {
        Profile::BlindBlock { width, height } => width * height,
        Profile::Pipe { diameter } => {
            0.5 * rim_sides as f64 * (diameter * 0.5).powi(2) * (2. * PI / rim_sides as f64).sin()
        }
    };
    let total: f64 = tiles.iter().map(|t| t.area).sum();
    if !total.is_finite() || (total / fitted_area - 1.).abs() > 1e-10 {
        return Err("section area does not match analytic profile".into());
    }
    let mut stations: Vec<_> = (0..=axial)
        .map(|i| core_length * (i as f64 / axial as f64))
        .collect();
    let mut axial_layers = vec![0; axial];
    if profile.blind() {
        for layer in (1..=layers).rev() {
            stations.push(stations.last().unwrap() + thicknesses[layer - 1]);
            axial_layers.push(layer);
        }
        *stations.last_mut().unwrap() = length;
    }
    let mut axial_ratios = vec![1.0_f64; axial_layers.len()];
    for (i, three) in stations.windows(3).enumerate() {
        let a = three[1] - three[0];
        let b = three[2] - three[1];
        if a <= 0. || b <= 0. {
            return Err("unrepresentable axial spacing".into());
        }
        let r = (a.max(b) / a.min(b)).cbrt();
        axial_ratios[i] = axial_ratios[i].max(r);
        axial_ratios[i + 1] = axial_ratios[i + 1].max(r);
    }
    let worst = tiles
        .iter()
        .map(|t| t.neighbour_ratio)
        .chain(axial_ratios.iter().copied())
        .fold(1., f64::max);
    if worst > ratio_limit + 1e-10 {
        return Err(format!(
            "measured neighbour size ratio {worst:.5} exceeds limit {ratio_limit:.5}"
        ));
    }
    Ok(Sweep {
        profile,
        points,
        tiles,
        outside,
        stations,
        axial_layers,
        axial_ratios,
        thicknesses,
        fitted_area,
        rim_sides,
    })
}

#[derive(Clone, Copy)]
pub enum Axis {
    X,
    Y,
    Z,
}

impl Sweep {
    pub fn vertices(&self, axis: Axis) -> Vec<Vec3> {
        let mut points = Vec::with_capacity(self.points.len() * self.stations.len());
        for &s in &self.stations {
            for p in &self.points {
                let (u, v) = match self.profile {
                    Profile::BlindBlock { width, height } => {
                        (p.y + 0.5 * width, p.z + 0.5 * height)
                    }
                    Profile::Pipe { .. } => (p.y, p.z),
                };
                points.push(match axis {
                    Axis::X => Vec3::new(s, u, v),
                    Axis::Y => Vec3::new(v, s, u),
                    Axis::Z => Vec3::new(u, v, s),
                });
            }
        }
        points
    }
}
