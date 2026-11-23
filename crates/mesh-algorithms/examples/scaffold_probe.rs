//! Tetrahedralize a remeshed surface MSH (triangles tagged 1..3).
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let text = std::fs::read_to_string(&a[1]).unwrap();
    let mut lines = text.lines();
    while lines.next() != Some("$Nodes") {}
    let n: usize = lines.next().unwrap().parse().unwrap();
    let mut p = Vec::new();
    for _ in 0..n { let f: Vec<f64> = lines.next().unwrap().split_whitespace().skip(1).map(|x| x.parse().unwrap()).collect(); p.push(mesh_core::geometry::Vec3::new(f[0], f[1], f[2])); }
    while lines.next() != Some("$Elements") {}
    let m: usize = lines.next().unwrap().parse().unwrap();
    let (mut t, mut g) = (Vec::new(), Vec::new());
    for _ in 0..m { let f: Vec<usize> = lines.next().unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect(); if f[1] == 2 { t.push([f[5]-1, f[6]-1, f[7]-1]); g.push(f[3]); } }
    let t0 = std::time::Instant::now();
    let (sc, st) = mesh_algorithms::scaffold::build(&p, &t, &g, Default::default()).unwrap();
    eprintln!("{st:?} in {:.1}s", t0.elapsed().as_secs_f64());
    if let Some(path) = a.get(2) {
        use std::io::Write;
        let mut out = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
        writeln!(out, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n{}", sc.points.len()).unwrap();
        for (i, p) in sc.points.iter().enumerate() { writeln!(out, "{} {:.17} {:.17} {:.17}", i + 1, p.x, p.y, p.z).unwrap(); }
        writeln!(out, "$EndNodes\n$Elements\n{}", sc.boundary.len() + sc.tetrahedra.len()).unwrap();
        let mut k = 0;
        for (t, g) in &sc.boundary { k += 1; writeln!(out, "{k} 2 2 {g} {g} {} {} {}", t[0]+1, t[1]+1, t[2]+1).unwrap(); }
        for t in &sc.tetrahedra { k += 1; writeln!(out, "{k} 4 2 4 1 {} {} {} {}", t[0]+1, t[1]+1, t[2]+1, t[3]+1).unwrap(); }
        writeln!(out, "$EndElements").unwrap();
    }
}
