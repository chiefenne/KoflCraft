// Analytic square ducts for the region proposal (regions.propose): unions
// of unit cubes, so the surface is exact and closed, with each boundary
// face tagged wall, inlet or outlet.
// Provenance (LIC-10): written for KoflCraft by its contributors, MIT like
// the repository; no tool, no third-party geometry.

#![allow(dead_code)]

use std::collections::{HashMap, HashSet};

use mesh_core::geometry::Vec3;

/// Facet roles: 0 wall, 1 inlet, 2 outlet.
#[derive(Default)]
pub struct Duct {
    pub vertices: Vec<Vec3>,
    pub triangles: Vec<[usize; 3]>,
    pub roles: Vec<usize>,
}

/// The closed surface of a union of unit cells. `port(cell, axis, side)`
/// gives the role of the boundary face of `cell` normal to `axis` on
/// `side` (+1 or -1).
pub fn cells(cells: &[[i32; 3]], port: impl Fn([i32; 3], usize, i32) -> usize) -> Duct {
    let set: HashSet<[i32; 3]> = cells.iter().copied().collect();
    let mut index: HashMap<[i32; 3], usize> = HashMap::new();
    let mut d = Duct {
        vertices: Vec::new(),
        triangles: Vec::new(),
        roles: Vec::new(),
    };
    let mut vertex = |p: [i32; 3], d: &mut Duct| {
        *index.entry(p).or_insert_with(|| {
            d.vertices
                .push(Vec3::new(p[0] as f64, p[1] as f64, p[2] as f64));
            d.vertices.len() - 1
        })
    };
    for &c in cells {
        for axis in 0..3 {
            for side in [-1, 1] {
                let mut n = c;
                n[axis] += side;
                if set.contains(&n) {
                    continue;
                }
                let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
                let mut base = c;
                if side > 0 {
                    base[axis] += 1;
                }
                let corner = |du: i32, dv: i32| {
                    let mut p = base;
                    p[u] += du;
                    p[v] += dv;
                    p
                };
                let mut q = [corner(0, 0), corner(1, 0), corner(1, 1), corner(0, 1)];
                if side < 0 {
                    q.reverse();
                }
                let q = q.map(|p| vertex(p, &mut d));
                let role = port(c, axis, side);
                d.triangles.push([q[0], q[1], q[2]]);
                d.triangles.push([q[0], q[2], q[3]]);
                d.roles.extend([role, role]);
            }
        }
    }
    d
}

/// A T: a straight duct along x (cells 0..8) with a side branch at x = 4
/// along y (cells 1..5); inlet at x = 0, outlets at x = 8 and y = 5.
pub fn tee() -> Duct {
    let mut c: Vec<[i32; 3]> = (0..8).map(|x| [x, 0, 0]).collect();
    c.extend((1..5).map(|y| [4, y, 0]));
    cells(&c, |cell, axis, side| match (cell, axis, side) {
        ([0, 0, 0], 0, -1) => 1,
        ([7, 0, 0], 0, 1) | ([4, 4, 0], 1, 1) => 2,
        _ => 0,
    })
}

/// A Y (fork): an inlet duct along x (cells 0..4) splits at x = 4 into two
/// parallel arms (y = -1 and y = 1, cells 5..10); inlet at x = 0, outlets
/// at the ends of the arms.
pub fn fork() -> Duct {
    let mut c: Vec<[i32; 3]> = (0..5).map(|x| [x, 0, 0]).collect();
    c.extend([[4, -1, 0], [4, 1, 0]]);
    for x in 5..10 {
        c.extend([[x, -1, 0], [x, 1, 0]]);
    }
    cells(&c, |cell, axis, side| match (cell, axis, side) {
        ([0, 0, 0], 0, -1) => 1,
        ([9, _, 0], 0, 1) => 2,
        _ => 0,
    })
}

/// The STL (ASCII) and a patch map of a duct.
pub fn files(d: &Duct) -> (Vec<u8>, String) {
    let mut s = String::from("solid duct\n");
    for t in &d.triangles {
        let [a, b, c] = t.map(|i| d.vertices[i]);
        let n = (b - a).cross(c - a);
        let n = n / n.length();
        s.push_str(&format!(
            "facet normal {} {} {}\nouter loop\n",
            n.x, n.y, n.z
        ));
        for p in [a, b, c] {
            s.push_str(&format!("vertex {} {} {}\n", p.x, p.y, p.z));
        }
        s.push_str("endloop\nendfacet\n");
    }
    s.push_str("endsolid duct\n");
    let roles: Vec<String> = d.roles.iter().map(|r| r.to_string()).collect();
    let map = format!(
        "{{\"patches\": [{{\"name\": \"walls\", \"role\": \"wall\"}}, {{\"name\": \"inlet\", \"role\": \"inlet\"}}, {{\"name\": \"outlet\", \"role\": \"outlet\"}}], \"facet_patch\": [{}]}}\n",
        roles.join(", ")
    );
    (s.into_bytes(), map)
}

/// A plenum with three round pipes: the box [-1, 1]^3 with a pipe of
/// radius 0.5 and length 4 on its -x face (inlet), its +x face and its +y
/// face (one outlet patch with two caps). Each box edge has 16 segments;
/// a face with a pipe is the ring between its square rim and the pipe's
/// circle, both at the rim points' polar angles, so the pipes have 64
/// sides. The region behind the three cuts of a sweep fill has three
/// openings.
pub fn pipe_box() -> Duct {
    let mut b = Builder::default();
    b.plenum(
        Vec3::ZERO,
        [
            Some((4.0, Some(2))),
            Some((4.0, Some(1))),
            Some((4.0, Some(2))),
            None,
            None,
            None,
        ],
    );
    b.d
}

/// Two plenums like [`pipe_box`]'s joined by one round pipe of length 4:
/// an inlet pipe on the -x face of the first, the joining pipe from its +x
/// face to the -x face of the second (centre x = 6) and an outlet pipe on
/// the second's +x face. The joining pipe is cut at both ends, and the
/// region behind the cuts is two closed shells.
pub fn box_pipe_box() -> Duct {
    let mut b = Builder::default();
    b.plenum(
        Vec3::ZERO,
        [
            Some((4.0, None)),
            Some((4.0, Some(1))),
            None,
            None,
            None,
            None,
        ],
    );
    b.plenum(
        Vec3::new(6.0, 0.0, 0.0),
        [
            Some((4.0, Some(2))),
            Some((0.0, None)),
            None,
            None,
            None,
            None,
        ],
    );
    b.d
}

#[derive(Default)]
struct Builder {
    d: Duct,
    index: HashMap<[i64; 3], usize>,
}

impl Builder {
    fn vertex(&mut self, p: Vec3) -> usize {
        let key = [p.x, p.y, p.z].map(|c| (c * 1e9).round() as i64);
        let d = &mut self.d;
        *self.index.entry(key).or_insert_with(|| {
            d.vertices.push(p);
            d.vertices.len() - 1
        })
    }

    fn tri(&mut self, t: [usize; 3], role: usize) {
        self.d.triangles.push(t);
        self.d.roles.push(role);
    }

    /// The box `centre` + [-1, 1]^3; `pipes` per face (+x, -x, +y, -y, +z,
    /// -z): a pipe of radius 0.5 and the given length, ended by a cap of
    /// the given role or left open (to meet the flush ring of another
    /// plenum's face).
    fn plenum(&mut self, centre: Vec3, pipes: [Option<(f64, Option<usize>)>; 6]) {
        const M: usize = 16;
        const R: f64 = 0.5;
        const K: usize = 40;
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 1.0, 0.0);
        let z = Vec3::new(0.0, 0.0, 1.0);
        // (outward normal, u, v with u x v = normal)
        let faces = [
            (x, y, z),
            (x * -1.0, z, y),
            (y, z, x),
            (y * -1.0, x, z),
            (z, x, y),
            (z * -1.0, y, x),
        ];
        for ((n, u, v), pipe) in faces.into_iter().zip(pipes) {
            // the square rim, counter-clockwise about n from (-1, -1)
            let mut rim: Vec<(f64, f64)> = Vec::new();
            for side in 0..4 {
                for k in 0..M {
                    let t = -1.0 + 2.0 * k as f64 / M as f64;
                    rim.push(match side {
                        0 => (t, -1.0),
                        1 => (1.0, t),
                        2 => (-t, 1.0),
                        _ => (-1.0, -t),
                    });
                }
            }
            let outer: Vec<usize> = rim
                .iter()
                .map(|&(a, b)| self.vertex(centre + n + u * a + v * b))
                .collect();
            let count = outer.len();
            let Some((length, cap)) = pipe else {
                let c = self.vertex(centre + n);
                for i in 0..count {
                    self.tri([c, outer[i], outer[(i + 1) % count]], 0);
                }
                continue;
            };
            // the pipe's rings at the rim points' polar angles, out along n
            let segments = if length > 0.0 { K } else { 0 };
            let mut rings: Vec<Vec<usize>> = Vec::with_capacity(segments + 1);
            for j in 0..=segments {
                let t = if segments > 0 {
                    length * j as f64 / segments as f64
                } else {
                    0.0
                };
                let ring = rim
                    .iter()
                    .map(|&(a, b)| {
                        let s = R / (a * a + b * b).sqrt();
                        self.vertex(centre + n * (1.0 + t) + u * (a * s) + v * (b * s))
                    })
                    .collect();
                rings.push(ring);
            }
            for i in 0..count {
                let i1 = (i + 1) % count;
                let inner = rings[0].clone();
                self.tri([outer[i], outer[i1], inner[i1]], 0);
                self.tri([outer[i], inner[i1], inner[i]], 0);
                for j in 0..segments {
                    let (a, b) = (rings[j][i], rings[j][i1]);
                    let (c, e) = (rings[j + 1][i1], rings[j + 1][i]);
                    self.tri([a, b, c], 0);
                    self.tri([a, c, e], 0);
                }
                if let Some(role) = cap {
                    let end = &rings[segments];
                    let (p, q) = (end[i], end[i1]);
                    let c = self.vertex(centre + n * (1.0 + length));
                    self.tri([c, p, q], role);
                }
            }
        }
    }
}

/// Two headers joined by three parallel tubes, from unit cubes: an inlet
/// duct (x = 0) into a header column (x = 1, y = -3..3), tubes at y = -2,
/// 0 and 2 (x = 2..6), a second header (x = 7) and an outlet duct
/// (x = 8). The tubes lie on loops of the slab graph but each is one duct.
/// No tube wall is flush with a header wall (a flush square branch is a
/// known gap of the join).
pub fn header_tubes() -> Duct {
    let mut c: Vec<[i32; 3]> = vec![[0, 0, 0], [8, 0, 0]];
    for y in -3..=3 {
        c.extend([[1, y, 0], [7, y, 0]]);
    }
    for x in 2..7 {
        c.extend([[x, -2, 0], [x, 0, 0], [x, 2, 0]]);
    }
    cells(&c, |cell, axis, side| match (cell, axis, side) {
        ([0, 0, 0], 0, -1) => 1,
        ([8, 0, 0], 0, 1) => 2,
        _ => 0,
    })
}

/// A body of revolution about the x axis: `profile` lists (x, r) with x
/// increasing, `sides` the segments around; the cap at the first x is the
/// inlet, the cap at the last the outlet.
pub fn revolved(profile: &[(f64, f64)], sides: usize) -> Duct {
    let mut d = Duct {
        vertices: Vec::new(),
        triangles: Vec::new(),
        roles: Vec::new(),
    };
    for &(x, r) in profile {
        for i in 0..sides {
            let a = std::f64::consts::TAU * i as f64 / sides as f64;
            d.vertices.push(Vec3::new(x, r * a.cos(), r * a.sin()));
        }
    }
    let at = |j: usize, i: usize| j * sides + i % sides;
    for j in 0..profile.len() - 1 {
        for i in 0..sides {
            let (a, b, c, e) = (at(j, i), at(j, i + 1), at(j + 1, i + 1), at(j + 1, i));
            d.triangles.extend([[a, b, c], [a, c, e]]);
            d.roles.extend([0, 0]);
        }
    }
    let last = profile.len() - 1;
    for (j, role, x) in [(0, 1, profile[0].0), (last, 2, profile[last].0)] {
        let centre = d.vertices.len();
        d.vertices.push(Vec3::new(x, 0.0, 0.0));
        for i in 0..sides {
            let (a, b) = (at(j, i), at(j, i + 1));
            d.triangles.push(if j == 0 {
                [centre, b, a]
            } else {
                [centre, a, b]
            });
            d.roles.push(role);
        }
    }
    d
}

/// A round concentric reducer: radius 0.5 over x = 0..3, a cone to radius
/// 0.3 over x = 3..4, radius 0.3 over x = 4..7; rings every 0.1, 48 sides.
pub fn reducer() -> Duct {
    let profile: Vec<(f64, f64)> = (0..=70)
        .map(|j| {
            let x = j as f64 / 10.0;
            let r = if x <= 3.0 {
                0.5
            } else if x >= 4.0 {
                0.3
            } else {
                0.5 - 0.2 * (x - 3.0)
            };
            (x, r)
        })
        .collect();
    revolved(&profile, 48)
}
