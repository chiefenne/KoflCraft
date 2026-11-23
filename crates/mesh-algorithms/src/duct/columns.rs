//! Wall-layer columns under the Voronoi core (M4).
//!
//! The core cells cut the innermost layer surface into wall polygons, one
//! per core cell and wall triangle (a *unit*). Each unit is extruded through
//! every layer along the layer nodes of its triangle: a vertex with
//! barycentric coordinates λ on the innermost triangle sits at λ on the same
//! triangle of every layer surface. Because the layer prism is affine in λ at
//! each height, the extruded units partition each coarse prism exactly, and
//! every cap of a unit is planar. Units of one core cell over adjacent wall
//! triangles that are coplanar on every layer surface are merged into one
//! column, so on a flat wall a core cell sits on one layer cell per layer,
//! whose centroid lies under its own (the poly-prism columns of the note).
//! Over a curved wall only coplanar triangles merge, so no layer cell has a
//! bent cap: FoehnCraft fans a merged face about its node average, and a
//! bent wall face would change the wall (volume) and, with a sag above the
//! first layer thickness, invert thin cells.
//!
//! A unit that covers its whole triangle and merges with nothing keeps its
//! native prism stack.
use super::voronoi_core::one_loop;
use super::*;
use std::collections::HashMap;

/// Merged triangles must be coplanar on the wall to this share of the model
/// size, so the hand-off's merged wall face is planar to roundoff.
const WALL_PLANE_TOL: f64 = 1e-10;
/// The same on the inner layer surfaces, relative to the thinnest layer.
const LAYER_PLANE_TOL: f64 = 1e-6;
/// Smallest cosine between a merged column's side face and the line to its
/// neighbour's centroid (85°); below it the column is split into its units.
const MIN_SIDE_COS: f64 = 0.087;
/// Passes that split merged columns whose contacts FoehnCraft would reject.
const SPLIT_PASSES: usize = 8;

pub(super) struct Columns {
    pub cells: Vec<Cell>,
    pub boundary: Vec<BoundaryFace>,
    /// Wall triangles that keep their native prism stack.
    pub native_triangles: usize,
    /// Poly-prism layer cells.
    pub poly_cells: usize,
    /// Poly-prism columns (one per layer is a layer cell).
    pub columns: usize,
    /// Units that share a column with another unit.
    pub merged_units: usize,
    /// Merged columns split back into units.
    pub split_columns: usize,
}

/// A side edge of a column on the innermost surface: its directed ring edge,
/// the wall triangle it lies in and the column across it (none at a port).
type Side = (usize, usize, usize, Option<usize>);

struct Unit {
    bi: usize,
    part: usize,
    rings: Vec<Vec<usize>>,
}

fn find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

/// Whether directed edges form one open chain.
fn one_chain(edges: &[(usize, usize)]) -> bool {
    let mut next = HashMap::new();
    let mut ends = std::collections::HashSet::new();
    for &(a, b) in edges {
        if next.insert(a, b).is_some() || !ends.insert(b) {
            return false;
        }
    }
    let starts: Vec<usize> = next.keys().copied().filter(|a| !ends.contains(a)).collect();
    let [start] = starts[..] else {
        return false;
    };
    let mut at = start;
    let mut count = 0;
    while let Some(&b) = next.get(&at) {
        count += 1;
        at = b;
        if count > edges.len() {
            return false;
        }
    }
    count == edges.len()
}

/// Build the layer cells of every wall triangle.
///
/// `tris` are the surface triangles, `wall` the wall facets in the order of
/// the coarse prisms (`coarse[level * wall.len() + r]`), `stack[level][v]`
/// the node of surface vertex `v` on layer surface `level`, and
/// `records[bi]` the core cells' wall polygons on the innermost surface.
#[allow(clippy::too_many_arguments)]
pub(super) fn build(
    points: &mut Vec<Vec3>,
    tris: &[[usize; 3]],
    wall: &[usize],
    coarse: &[Cell],
    old_boundary: Vec<BoundaryFace>,
    stack: &[Vec<usize>],
    rim: &HashMap<(usize, usize), usize>,
    records: &BTreeMap<usize, Vec<(usize, Vec<usize>)>>,
    scale: f64,
    thinnest: f64,
) -> Result<Columns, String> {
    let layers = stack.len() - 1;
    let rank: HashMap<usize, usize> = wall.iter().enumerate().map(|(r, &bi)| (bi, r)).collect();

    // Barycentric coordinates of the ring vertices and their images on
    // every layer surface.
    let mut bary: HashMap<(usize, usize), [f64; 3]> = HashMap::new();
    let mut image: HashMap<(usize, usize), usize> = HashMap::new();
    for &bi in wall {
        let corners = tris[bi];
        let inner = corners.map(|v| stack[layers][v]);
        let tri = inner.map(|n| points[n]);
        let normal = (tri[1] - tri[0]).cross(tri[2] - tri[0]);
        let n2 = normal.dot(normal);
        let rec = records
            .get(&bi)
            .ok_or_else(|| format!("wall triangle {bi} has no core face"))?;
        for (_, ring) in rec {
            for &g in ring {
                if bary.contains_key(&(g, bi)) {
                    continue;
                }
                let corner = inner.iter().position(|&n| n == g);
                let l = match corner {
                    Some(j) => {
                        let mut l = [0.0; 3];
                        l[j] = 1.0;
                        l
                    }
                    None => {
                        let x = points[g];
                        let la = (tri[1] - x).cross(tri[2] - x).dot(normal) / n2;
                        let lb = (tri[2] - x).cross(tri[0] - x).dot(normal) / n2;
                        let mut l = [la, lb, 1.0 - la - lb];
                        for v in &mut l {
                            if v.abs() < 1e-9 {
                                *v = 0.0;
                            }
                        }
                        let s: f64 = l.iter().sum();
                        l.map(|v| v / s)
                    }
                };
                bary.insert((g, bi), l);
                for k in 0..layers {
                    if image.contains_key(&(g, k)) {
                        continue;
                    }
                    let id = match corner {
                        Some(j) => stack[k][corners[j]],
                        None => {
                            let p = (0..3)
                                .fold(Vec3::ZERO, |s, j| s + points[stack[k][corners[j]]] * l[j]);
                            points.push(p);
                            points.len() - 1
                        }
                    };
                    image.insert((g, k), id);
                }
            }
        }
    }
    let points: &[Vec3] = points;
    let img = |g: usize, k: usize| -> usize { if k == layers { g } else { image[&(g, k)] } };

    // Units, and the unit on the left of every directed wall ring edge.
    let mut units: Vec<Unit> = Vec::new();
    let mut unit_of: HashMap<(usize, usize), usize> = HashMap::new();
    for &bi in wall {
        for (part, ring) in &records[&bi] {
            let u = *unit_of.entry((*part, bi)).or_insert_with(|| {
                units.push(Unit {
                    bi,
                    part: *part,
                    rings: Vec::new(),
                });
                units.len() - 1
            });
            units[u].rings.push(ring.clone());
        }
    }
    let mut edge_unit: HashMap<(usize, usize), usize> = HashMap::new();
    for (u, unit) in units.iter().enumerate() {
        for ring in &unit.rings {
            for k in 0..ring.len() {
                if edge_unit
                    .insert((ring[k], ring[(k + 1) % ring.len()]), u)
                    .is_some()
                {
                    return Err("a wall ring edge is used twice".into());
                }
            }
        }
    }

    // Merge units of one core cell over coplanar adjacent triangles.
    let mut plane_cache: HashMap<(usize, usize), bool> = HashMap::new();
    let mut coplanar = |t1: usize, t2: usize| -> bool {
        *plane_cache
            .entry((t1.min(t2), t1.max(t2)))
            .or_insert_with(|| {
                (0..=layers).all(|k| {
                    let p = tris[t1].map(|v| points[stack[k][v]]);
                    let q = tris[t2].map(|v| points[stack[k][v]]);
                    let n = (p[1] - p[0]).cross(p[2] - p[0]);
                    let m = (q[1] - q[0]).cross(q[2] - q[0]);
                    let len = n.length();
                    if len.is_nan() || len <= 0.0 || n.dot(m) <= 0.0 {
                        return false;
                    }
                    let unit = n / len;
                    let tol = if k == 0 {
                        WALL_PLANE_TOL * scale
                    } else {
                        LAYER_PLANE_TOL * thinnest
                    };
                    q.iter().all(|x| (*x - p[0]).dot(unit).abs() <= tol)
                })
            })
    };
    let mut parent: Vec<usize> = (0..units.len()).collect();
    for u in 0..units.len() {
        for ring in &units[u].rings {
            for k in 0..ring.len() {
                let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
                if let Some(&w) = edge_unit.get(&(b, a))
                    && w != u
                    && units[w].part == units[u].part
                    && coplanar(units[u].bi, units[w].bi)
                {
                    let (x, y) = (find(&mut parent, u), find(&mut parent, w));
                    if x != y {
                        parent[x.max(y)] = x.min(y);
                    }
                }
            }
        }
    }
    let mut by_root: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for u in 0..units.len() {
        let r = find(&mut parent, u);
        by_root.entry(r).or_default().push(u);
    }
    let mut columns: Vec<Vec<usize>> = by_root.into_values().collect();

    // FoehnCraft allows one face per cell pair: a merged column must be one
    // disk, and its contact with each neighbouring column one chain.
    let mut split_columns = 0;
    let mut col_of = vec![0; units.len()];
    let mut sides: Vec<Vec<Side>>;
    let mut pass = 0;
    loop {
        for (x, col) in columns.iter().enumerate() {
            for &u in col {
                col_of[u] = x;
            }
        }
        sides = Vec::with_capacity(columns.len());
        let mut bad = std::collections::BTreeSet::new();
        for (x, col) in columns.iter().enumerate() {
            let rings: Vec<&Vec<usize>> = col.iter().flat_map(|&u| &units[u].rings).collect();
            if col.len() > 1 && !one_loop(&rings) {
                bad.insert(x);
            }
            let mut edges = Vec::new();
            let mut by: BTreeMap<Option<usize>, Vec<(usize, usize)>> = BTreeMap::new();
            for &u in col {
                for ring in &units[u].rings {
                    for k in 0..ring.len() {
                        let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
                        let other = edge_unit.get(&(b, a)).map(|&w| col_of[w]);
                        if other == Some(x) {
                            continue;
                        }
                        edges.push((a, b, units[u].bi, other));
                        by.entry(other).or_default().push((a, b));
                    }
                }
            }
            for (other, list) in by {
                let Some(y) = other else { continue };
                if !one_chain(&list) {
                    if col.len() > 1 {
                        bad.insert(x);
                    }
                    if columns[y].len() > 1 {
                        bad.insert(y);
                    }
                    if col.len() == 1 && columns[y].len() == 1 {
                        return Err(format!(
                            "wall polygons of core cells {} and {} touch in more than one chain",
                            units[col[0]].part, units[columns[y][0]].part
                        ));
                    }
                }
            }
            sides.push(edges);
        }
        // FoehnCraft's orthogonality test (S·d > 0) between a merged column
        // and each neighbour, on every layer. A merged footprint can be
        // nonconvex (a core cell with merged or reassigned pieces), and its
        // centroid can then lie behind a side face.
        let mut centre: HashMap<(usize, usize), Vec3> = HashMap::new();
        let mut centroid = |x: usize, k: usize| -> Vec3 {
            *centre.entry((x, k)).or_insert_with(|| {
                let (mut v, mut m) = (0.0, Vec3::ZERO);
                for ring in columns[x].iter().flat_map(|&u| &units[u].rings) {
                    for i in 1..ring.len() - 1 {
                        let w = [ring[0], ring[i + 1], ring[i]];
                        let n = [
                            img(w[0], k),
                            img(w[1], k),
                            img(w[2], k),
                            img(w[0], k + 1),
                            img(w[1], k + 1),
                            img(w[2], k + 1),
                        ];
                        let vol = verify::wedge(n.map(|i| points[i])).1;
                        let c = n.iter().fold(Vec3::ZERO, |s, &i| s + points[i]) / 6.0;
                        v += vol;
                        m = m + c * vol;
                    }
                }
                m / v
            })
        };
        for x in 0..columns.len() {
            let mut by: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
            for &(a, b, _, other) in &sides[x] {
                if let Some(y) = other {
                    by.entry(y).or_default().push((a, b));
                }
            }
            for (y, list) in by {
                if columns[x].len() == 1 && columns[y].len() == 1 {
                    continue;
                }
                let skewed = (0..layers).any(|k| {
                    let area = list.iter().fold(Vec3::ZERO, |s, &(a, b)| {
                        let q =
                            [img(b, k), img(a, k), img(a, k + 1), img(b, k + 1)].map(|i| points[i]);
                        s + (q[2] - q[0]).cross(q[3] - q[1]) * 0.5
                    });
                    let d = centroid(y, k) - centroid(x, k);
                    let cos = area.dot(d) / (area.length() * d.length());
                    // NaN (a degenerate face) counts as skewed
                    cos.is_nan() || cos <= MIN_SIDE_COS
                });
                if skewed {
                    for z in [x, y] {
                        if columns[z].len() > 1 {
                            bad.insert(z);
                        }
                    }
                }
            }
        }
        if bad.is_empty() {
            break;
        }
        pass += 1;
        if pass > SPLIT_PASSES {
            return Err("merged wall columns kept failing the contact check".into());
        }
        split_columns += bad.len();
        let mut next = Vec::with_capacity(columns.len());
        for (x, col) in columns.into_iter().enumerate() {
            if bad.contains(&x) {
                next.extend(col.into_iter().map(|u| vec![u]));
            } else {
                next.push(col);
            }
        }
        columns = next;
    }

    // Cells, level by level.
    let native: Vec<bool> = columns
        .iter()
        .map(|col| {
            let unit = &units[col[0]];
            let inner = tris[unit.bi].map(|v| stack[layers][v]);
            col.len() == 1
                && unit.rings.len() == 1
                && unit.rings[0].len() == 3
                && unit.rings[0].iter().all(|g| inner.contains(g))
        })
        .collect();
    let mut cells = Vec::new();
    let mut boundary = Vec::new();
    let mut new_index = vec![usize::MAX; coarse.len()];
    let mut poly_cells = 0;
    for level in 0..layers {
        for (x, col) in columns.iter().enumerate() {
            if native[x] {
                let old = level * wall.len() + rank[&units[col[0]].bi];
                let Element::Wedge(n) = coarse[old].element else {
                    return Err("expected a layer prism".into());
                };
                new_index[old] = cells.len();
                cells.push(Cell {
                    element: Element::Wedge(n),
                    layer: level + 1,
                    region: 0,
                });
                continue;
            }
            let id = cells.len();
            let mut faces = Vec::new();
            let mut wedges = Vec::new();
            for ring in col.iter().flat_map(|&u| &units[u].rings) {
                let bottom: Vec<usize> = ring.iter().map(|&g| img(g, level)).collect();
                let top: Vec<usize> = ring.iter().rev().map(|&g| img(g, level + 1)).collect();
                if level == 0 {
                    boundary.push(BoundaryFace {
                        nodes: bottom.clone(),
                        patch: 3,
                        owner: id,
                    });
                }
                faces.push(bottom);
                faces.push(top);
                for i in 1..ring.len() - 1 {
                    let (a, b, c) = (ring[0], ring[i], ring[i + 1]);
                    wedges.push([
                        img(a, level),
                        img(c, level),
                        img(b, level),
                        img(a, level + 1),
                        img(c, level + 1),
                        img(b, level + 1),
                    ]);
                }
            }
            for &(a, b, bi, other) in &sides[x] {
                let quad = vec![
                    img(b, level),
                    img(a, level),
                    img(a, level + 1),
                    img(b, level + 1),
                ];
                if other.is_none() {
                    let (la, lb) = (bary[&(a, bi)], bary[&(b, bi)]);
                    let j = (0..3)
                        .find(|&j| la[j] == 0.0 && lb[j] == 0.0)
                        .ok_or("open wall ring edge off the triangle's rim")?;
                    let (u, v) = (tris[bi][(j + 1) % 3], tris[bi][(j + 2) % 3]);
                    let patch = *rim
                        .get(&(u.min(v), u.max(v)))
                        .ok_or("open wall ring edge that is no port rim")?;
                    boundary.push(BoundaryFace {
                        nodes: quad.clone(),
                        patch,
                        owner: id,
                    });
                }
                faces.push(quad);
            }
            let expected: f64 = wedges
                .iter()
                .map(|w| verify::wedge(w.map(|i| points[i])).1)
                .sum();
            let nodes: Vec<usize> = faces.iter().flatten().copied().collect();
            let interior =
                nodes.iter().fold(Vec3::ZERO, |s, &i| s + points[i]) / nodes.len() as f64;
            let mut poly = Polyhedron::new(faces, interior, None, expected);
            poly.prism = Some(wedges);
            cells.push(Cell {
                element: Element::Polyhedron(Box::new(poly)),
                layer: level + 1,
                region: 0,
            });
            poly_cells += 1;
        }
    }
    for f in old_boundary {
        // only native prism stacks keep their old wall and port faces
        if new_index[f.owner] != usize::MAX {
            boundary.push(BoundaryFace {
                nodes: f.nodes,
                patch: f.patch,
                owner: new_index[f.owner],
            });
        }
    }
    let merged_units = columns.iter().filter(|c| c.len() > 1).map(Vec::len).sum();
    Ok(Columns {
        cells,
        boundary,
        native_triangles: native.iter().filter(|&&n| n).count(),
        poly_cells,
        columns: native.iter().filter(|&&n| !n).count(),
        merged_units,
        split_columns,
    })
}

#[cfg(test)]
mod tests {
    use super::one_chain;

    /// Keep: two columns may touch in one open chain of edges only, since
    /// FoehnCraft allows one face per cell pair and the hand-off merges a
    /// pair's side quads into one polygon. Accepting two chains or a closed
    /// loop (a column surrounded by one neighbour) lets a merged column
    /// through whose hand-off then fails; the fixtures do not contain such a
    /// layout, so only this test sees it.
    #[test]
    fn contact_chains() {
        assert!(one_chain(&[(1, 2), (2, 3)]));
        assert!(!one_chain(&[(1, 2), (5, 6)]));
        // a closed loop around a column is no single face either
        assert!(!one_chain(&[(1, 2), (2, 3), (3, 1)]));
    }
}
