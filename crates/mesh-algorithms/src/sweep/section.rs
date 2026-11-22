use super::{Profile, Tile};
use mesh_core::geometry::Vec3;
use std::collections::HashMap;
use std::f64::consts::{PI, TAU};

fn yz(y: f64, z: f64) -> Vec3 {
    Vec3::new(0., y, z)
}
fn rectangle(y: f64, z: f64) -> Vec<Vec3> {
    vec![yz(-y, -z), yz(y, -z), yz(y, z), yz(-y, z)]
}
fn radial(r: f64, turn: f64) -> Vec3 {
    let (s, c) = (TAU * turn).sin_cos();
    yz(r * c, r * s)
}

fn sites(
    profile: Profile,
    h: f64,
    stack: f64,
    budget: usize,
) -> Result<(Vec<Vec3>, Vec<Vec3>), String> {
    let mut sites = Vec::new();
    let boundary = match profile {
        Profile::BlindBlock { width, height } => {
            let w = width - 2. * stack;
            let d = height - 2. * stack;
            let nx = (w / h).ceil().max(2.) as usize;
            let ny = (d / (h * 3.0_f64.sqrt() * 0.5)).ceil().max(2.) as usize;
            if nx.checked_mul(ny).is_none_or(|n| n > budget) {
                return Err("section sites exceed cell budget".into());
            }
            let dx = w / nx as f64;
            let dy = d / ny as f64;
            for j in 0..ny {
                for i in 0..nx {
                    let shift = if i == 0 || i == nx - 1 {
                        0.
                    } else if j % 2 == 0 {
                        -0.25
                    } else {
                        0.25
                    };
                    sites.push(yz(
                        (i as f64 + 0.5 + shift) * dx - w * 0.5,
                        (j as f64 + 0.5) * dy - d * 0.5,
                    ));
                }
            }
            rectangle(w * 0.5, d * 0.5)
        }
        Profile::Pipe { diameter } => {
            // A boundary site faces each chord, so a prism side spans a whole
            // chord. Inner rings progressively reduce their circumferential count.
            let r = diameter * 0.5;
            let n = (TAU * (r - stack) / h).ceil().max(12.) as usize;
            let n = n.checked_next_multiple_of(4).ok_or("rim count overflow")?;
            if n > budget {
                return Err("rim sites exceed cell budget".into());
            }
            let inner = r - stack / (PI / n as f64).cos();
            if inner <= 0. {
                return Err("no positive inscribed pipe core".into());
            }
            let rings = (inner / (0.8660254037844386 * h) - 0.5).round().max(1.) as usize;
            let dr = inner / (rings as f64 + 0.5);
            sites.push(Vec3::ZERO);
            for i in 1..=rings {
                let radius = dr * i as f64;
                let count = if i == rings {
                    n
                } else {
                    (TAU * radius / h).round().max(6.) as usize
                };
                if count > budget - sites.len() {
                    return Err("section sites exceed cell budget".into());
                }
                let phase = if i == rings || i % 2 == 0 { 0.5 } else { 0. };
                sites.extend((0..count).map(|j| radial(radius, (j as f64 + phase) / count as f64)));
            }
            (0..n).map(|i| radial(inner, i as f64 / n as f64)).collect()
        }
    };
    Ok((sites, boundary))
}

// Sutherland–Hodgman half-plane clipping. All coordinates are relative to
// the site while clipping neighbouring bisectors.
fn clip(ring: &mut Vec<Vec3>, origin: Vec3, normal: Vec3) {
    if ring.iter().all(|&p| (p - origin).dot(normal) <= 0.) {
        return;
    }
    let mut next = Vec::with_capacity(ring.len() + 1);
    for (&a, &b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
        let da = (a - origin).dot(normal);
        let db = (b - origin).dot(normal);
        if da <= 0. {
            next.push(a);
        }
        if (da < 0. && db > 0.) || (da > 0. && db < 0.) {
            next.push(a + (b - a) * (da / (da - db)));
        }
    }
    *ring = next;
}

struct Vertices {
    points: Vec<Vec3>,
    bins: HashMap<(i64, i64), Vec<usize>>,
    tol: f64,
}
impl Vertices {
    fn insert(&mut self, p: Vec3) -> usize {
        let x = (p.y / self.tol).floor() as i64;
        let y = (p.z / self.tol).floor() as i64;
        for i in x - 1..=x + 1 {
            for j in y - 1..=y + 1 {
                if let Some(bucket) = self.bins.get(&(i, j)) {
                    for &id in bucket {
                        if (p - self.points[id]).max_abs_component() <= self.tol {
                            return id;
                        }
                    }
                }
            }
        }
        let id = self.points.len();
        self.points.push(p);
        self.bins.entry((x, y)).or_default().push(id);
        id
    }
}

pub(super) struct Section {
    pub points: Vec<Vec3>,
    pub tiles: Vec<Tile>,
    pub boundary: Vec<Vec3>,
}

pub(super) fn core(profile: Profile, h: f64, stack: f64, budget: usize) -> Result<Section, String> {
    let (sites, boundary) = sites(profile, h, stack, budget)?;
    let bin = |p: Vec3| ((p.y / h).floor() as i64, (p.z / h).floor() as i64);
    let mut bins: HashMap<_, Vec<usize>> = HashMap::new();
    for (i, &p) in sites.iter().enumerate() {
        bins.entry(bin(p)).or_default().push(i);
    }
    let tol = boundary
        .iter()
        .map(|p| p.max_abs_component())
        .fold(h, f64::max)
        * 2e-12;
    let mut vertices = Vertices {
        points: Vec::new(),
        bins: HashMap::new(),
        tol,
    };
    let mut tiles = Vec::with_capacity(sites.len());
    for (id, &p) in sites.iter().enumerate() {
        // Any site farther than twice this square's circumradius cannot cut
        // the cell. The final check certifies the square introduced no face.
        let extent = 2. * h;
        let mut ring = rectangle(extent, extent);
        let (x, y) = bin(p);
        let mut nearby = Vec::new();
        for i in x - 6..=x + 6 {
            for j in y - 6..=y + 6 {
                if let Some(other) = bins.get(&(i, j)) {
                    for &k in other {
                        let d = sites[k] - p;
                        if k != id && d.dot(d) <= 8. * extent * extent {
                            nearby.push((d.dot(d), k));
                        }
                    }
                }
            }
        }
        nearby.sort_unstable_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        for (_, k) in nearby {
            let d = sites[k] - p;
            clip(&mut ring, d * 0.5, d);
        }
        for (&a, &b) in boundary.iter().zip(boundary.iter().cycle().skip(1)) {
            let e = b - a;
            clip(&mut ring, a - p, yz(e.z, -e.y));
        }
        if ring.iter().any(|q| q.max_abs_component() >= extent - tol) {
            return Err("insufficient site coverage at section boundary".into());
        }
        let mut ids: Vec<_> = ring.into_iter().map(|q| vertices.insert(q + p)).collect();
        ids.dedup();
        if ids.len() > 1 && ids.first() == ids.last() {
            ids.pop();
        }
        if ids.len() < 3 {
            return Err("collapsed section cell".into());
        }
        tiles.push(Tile::new(ids, 0, &vertices.points)?);
    }
    Ok(Section {
        points: vertices.points,
        tiles,
        boundary,
    })
}
