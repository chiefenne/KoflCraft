// Analytic square ducts for the region proposal (regions.propose): unions
// of unit cubes, so the surface is exact and closed, with each boundary
// face tagged wall, inlet or outlet.
// Provenance (LIC-10): written for KoflCraft by its contributors, MIT like
// the repository; no tool, no third-party geometry.

#![allow(dead_code)]

use std::collections::{HashMap, HashSet};

use mesh_core::geometry::Vec3;

/// Facet roles: 0 wall, 1 inlet, 2 outlet.
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
