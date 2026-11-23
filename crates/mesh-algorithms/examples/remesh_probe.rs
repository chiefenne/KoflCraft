use mesh_algorithms::remesh::{RemeshOptions, remesh};
use mesh_geometry::{StlOptions, import_closed_stl};
use std::collections::HashMap;
use std::io::Write;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let bytes = std::fs::read(&a[1]).unwrap();
    let tol: f64 = a[3].parse().unwrap();
    let s = import_closed_stl(&bytes, StlOptions { facet_patches: None, facet_limit: None }).unwrap();
    let tris: Vec<[usize; 3]> = s.facets().iter().map(|f| f.triangle.vertices).collect();
    let mut tag = vec![3usize; tris.len()];
    for (t, r) in [(1, 90374..90884), (2, 21168..21680)] {
        for i in r { tag[i] = t; }
    }
    let keep: Vec<bool> = tag.iter().map(|&t| t != 3).collect();
    let t0 = std::time::Instant::now();
    let r = remesh(s.vertices(), &tris, &keep, RemeshOptions { tolerance: tol, rim_size: Some(0.85), ..Default::default() }).unwrap();
    eprintln!("{:?} wall tris {} in {:.1}s", r.stats, r.source.iter().filter(|s| s.is_none()).count(), t0.elapsed().as_secs_f64());
    let mut out = std::io::BufWriter::new(std::fs::File::create(&a[2]).unwrap());
    writeln!(out, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n{}", r.vertices.len()).unwrap();
    for (i, p) in r.vertices.iter().enumerate() { writeln!(out, "{} {:.17} {:.17} {:.17}", i + 1, p.x, p.y, p.z).unwrap(); }
    // rim edges: edges between kept and new facets
    let mut owner: HashMap<(usize, usize), Vec<bool>> = HashMap::new();
    for (t, s) in r.triangles.iter().zip(&r.source) {
        for k in 0..3 { let (x, y) = (t[k], t[(k + 1) % 3]); owner.entry((x.min(y), x.max(y))).or_default().push(s.is_some()); }
    }
    let mut rim: Vec<_> = owner.iter().filter(|(_, v)| v.len() == 2 && v[0] != v[1]).map(|(k, _)| *k).collect();
    rim.sort();
    writeln!(out, "$EndNodes\n$Elements\n{}", r.triangles.len() + rim.len()).unwrap();
    for (i, (t, s)) in r.triangles.iter().zip(&r.source).enumerate() {
        let p = s.map_or(3, |f| tag[f]);
        writeln!(out, "{} 2 2 {} {} {} {} {}", i + 1, p, p, t[0] + 1, t[1] + 1, t[2] + 1).unwrap();
    }
    for (j, (x, y)) in rim.iter().enumerate() { writeln!(out, "{} 1 2 0 {} {} {}", r.triangles.len() + j + 1, j + 1, x + 1, y + 1).unwrap(); }
    writeln!(out, "$EndElements").unwrap();
    eprintln!("rim edges {}", rim.len());
    let tags: Vec<usize> = r.source.iter().map(|s| s.map_or(3, |f| tag[f])).collect();
    let t1 = std::time::Instant::now();
    let (sc, st) = mesh_algorithms::scaffold::build(&r.vertices, &r.triangles, &tags, Default::default()).unwrap();
    eprintln!("{:?} in {:.1}s; STL volume {}", st, t1.elapsed().as_secs_f64(), s.volume());
    let mut out = std::io::BufWriter::new(std::fs::File::create(&a[4]).unwrap());
    writeln!(out, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n{}", sc.points.len()).unwrap();
    for (i, p) in sc.points.iter().enumerate() { writeln!(out, "{} {:.17} {:.17} {:.17}", i + 1, p.x, p.y, p.z).unwrap(); }
    writeln!(out, "$EndNodes\n$Elements\n{}", sc.boundary.len() + sc.tetrahedra.len()).unwrap();
    let mut k = 0;
    for (t, g) in &sc.boundary { k += 1; writeln!(out, "{k} 2 2 {g} {g} {} {} {}", t[0]+1, t[1]+1, t[2]+1).unwrap(); }
    for t in &sc.tetrahedra { k += 1; writeln!(out, "{k} 4 2 4 1 {} {} {} {}", t[0]+1, t[1]+1, t[2]+1, t[3]+1).unwrap(); }
    writeln!(out, "$EndElements").unwrap();
}
