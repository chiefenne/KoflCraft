// Convex test surfaces for geometry-derived sizing (GRD-04), generated in
// code: vertices and outward-wound triangles, and an ASCII STL writer.
// Provenance (LIC-10): written for KoflCraft by its contributors from the
// analytic shapes named below, MIT like the repository; no tool, no
// third-party geometry.

#![allow(dead_code)]

use mesh_core::geometry::Vec3;

pub type Shape = (Vec<Vec3>, Vec<[usize; 3]>);

/// The box [0, max], two triangles per side.
pub fn cuboid(max: Vec3) -> Shape {
    let v = |x: f64, y: f64, z: f64| Vec3::new(x * max.x, y * max.y, z * max.z);
    let vertices = vec![
        v(0.0, 0.0, 0.0),
        v(1.0, 0.0, 0.0),
        v(1.0, 1.0, 0.0),
        v(0.0, 1.0, 0.0),
        v(0.0, 0.0, 1.0),
        v(1.0, 0.0, 1.0),
        v(1.0, 1.0, 1.0),
        v(0.0, 1.0, 1.0),
    ];
    let quads = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [2, 3, 7, 6],
        [0, 4, 7, 3],
        [1, 2, 6, 5],
    ];
    let triangles = quads
        .iter()
        .flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]])
        .collect();
    (vertices, triangles)
}

/// A sphere of radius `r` about `centre`: the icosahedron subdivided
/// `levels` times, vertices pushed onto the sphere.
pub fn icosphere(centre: Vec3, r: f64, levels: usize) -> Shape {
    let t = (1.0 + 5f64.sqrt()) / 2.0;
    let mut v: Vec<Vec3> = [
        (-1.0, t, 0.0),
        (1.0, t, 0.0),
        (-1.0, -t, 0.0),
        (1.0, -t, 0.0),
        (0.0, -1.0, t),
        (0.0, 1.0, t),
        (0.0, -1.0, -t),
        (0.0, 1.0, -t),
        (t, 0.0, -1.0),
        (t, 0.0, 1.0),
        (-t, 0.0, -1.0),
        (-t, 0.0, 1.0),
    ]
    .iter()
    .map(|&(x, y, z)| Vec3::new(x, y, z))
    .collect();
    let mut f: Vec<[usize; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    for _ in 0..levels {
        let mut mid = std::collections::BTreeMap::new();
        let mut next = Vec::with_capacity(4 * f.len());
        for tri in &f {
            let mut m = [0; 3];
            for k in 0..3 {
                let (a, b) = (tri[k], tri[(k + 1) % 3]);
                m[k] = *mid.entry((a.min(b), a.max(b))).or_insert_with(|| {
                    v.push((v[a] + v[b]) / 2.0);
                    v.len() - 1
                });
            }
            next.push([tri[0], m[0], m[2]]);
            next.push([tri[1], m[1], m[0]]);
            next.push([tri[2], m[2], m[1]]);
            next.push(m);
        }
        f = next;
    }
    let v = v
        .into_iter()
        .map(|p| centre + p * (r / p.length()))
        .collect();
    (v, f)
}

/// A latitude-longitude sphere: `rings` bands, `sectors` around, fans at
/// the poles (the layout of many CAD exports of a ball).
pub fn uv_sphere(centre: Vec3, r: f64, rings: usize, sectors: usize) -> Shape {
    let mut v = vec![
        centre + Vec3::new(0.0, 0.0, r),
        centre - Vec3::new(0.0, 0.0, r),
    ];
    for i in 1..rings {
        let theta = std::f64::consts::PI * i as f64 / rings as f64;
        for j in 0..sectors {
            let phi = std::f64::consts::TAU * j as f64 / sectors as f64;
            v.push(
                centre
                    + Vec3::new(
                        r * theta.sin() * phi.cos(),
                        r * theta.sin() * phi.sin(),
                        r * theta.cos(),
                    ),
            );
        }
    }
    let at = |i: usize, j: usize| 2 + (i - 1) * sectors + j % sectors;
    let mut f = Vec::new();
    for j in 0..sectors {
        f.push([0, at(1, j), at(1, j + 1)]);
        f.push([1, at(rings - 1, j + 1), at(rings - 1, j)]);
    }
    for i in 1..rings - 1 {
        for j in 0..sectors {
            f.push([at(i, j), at(i + 1, j), at(i + 1, j + 1)]);
            f.push([at(i, j), at(i + 1, j + 1), at(i, j + 1)]);
        }
    }
    (v, f)
}

/// An ellipsoid with semi-axes `a`: the icosphere stretched.
pub fn ellipsoid(centre: Vec3, a: Vec3, levels: usize) -> Shape {
    let (v, f) = icosphere(Vec3::ZERO, 1.0, levels);
    let v = v
        .into_iter()
        .map(|p| centre + Vec3::new(p.x * a.x, p.y * a.y, p.z * a.z))
        .collect();
    (v, f)
}

/// A closed cylinder of radius `r` and length `l` along z, its axis
/// through `(cx, cy)`, with `n` sides of long strip triangles and fan caps,
/// as CAD exporters write it.
pub fn strip_cylinder(cx: f64, cy: f64, r: f64, l: f64, n: usize) -> Shape {
    let mut v = vec![Vec3::new(cx, cy, 0.0), Vec3::new(cx, cy, l)];
    for i in 0..n {
        let a = std::f64::consts::TAU * i as f64 / n as f64;
        v.push(Vec3::new(cx + r * a.cos(), cy + r * a.sin(), 0.0));
        v.push(Vec3::new(cx + r * a.cos(), cy + r * a.sin(), l));
    }
    let (b, t) = (|i: usize| 2 + 2 * (i % n), |i: usize| 3 + 2 * (i % n));
    let mut f = Vec::new();
    for i in 0..n {
        f.push([b(i), b(i + 1), t(i + 1)]);
        f.push([b(i), t(i + 1), t(i)]);
        f.push([0, b(i + 1), b(i)]);
        f.push([1, t(i), t(i + 1)]);
    }
    (v, f)
}

/// A wedge: the triangle with its apex at the origin, opening `angle`
/// (radians) symmetrically about +x up to x = `length`, extruded along y
/// to `width`.
pub fn wedge(angle: f64, length: f64, width: f64) -> Shape {
    let h = length * (angle / 2.0).tan();
    let v = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(length, 0.0, -h),
        Vec3::new(length, 0.0, h),
        Vec3::new(0.0, width, 0.0),
        Vec3::new(length, width, -h),
        Vec3::new(length, width, h),
    ];
    let f = vec![
        [0, 1, 2],
        [3, 5, 4],
        [0, 1, 4],
        [0, 4, 3],
        [1, 2, 5],
        [1, 5, 4],
        [2, 0, 3],
        [2, 3, 5],
    ];
    outward((v, f))
}

/// Winds every triangle of a convex shape so its normal points away from
/// the vertex mean.
pub fn outward((v, mut f): Shape) -> Shape {
    let mean = v.iter().fold(Vec3::ZERO, |s, &p| s + p) / v.len() as f64;
    for t in &mut f {
        let [a, b, c] = t.map(|i| v[i]);
        if (b - a).cross(c - a).dot(a - mean) < 0.0 {
            t.swap(1, 2);
        }
    }
    (v, f)
}

/// The shape as ASCII STL, one solid.
pub fn ascii_stl((v, f): &Shape) -> Vec<u8> {
    let mut s = String::from("solid surface\n");
    for t in f {
        let [a, b, c] = t.map(|i| v[i]);
        let n = (b - a).cross(c - a);
        let n = n / n.length();
        s.push_str(&format!(
            "facet normal {:e} {:e} {:e}\nouter loop\n",
            n.x, n.y, n.z
        ));
        for p in [a, b, c] {
            s.push_str(&format!("vertex {:e} {:e} {:e}\n", p.x, p.y, p.z));
        }
        s.push_str("endloop\nendfacet\n");
    }
    s.push_str("endsolid surface\n");
    s.into_bytes()
}
