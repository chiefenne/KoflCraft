//! AABB broad phase and triangle contact checks for the input shell and layer front.
use mesh_core::geometry::Vec3;

#[derive(Clone, Copy)]
struct Bounds {
    lo: [f64; 3],
    hi: [f64; 3],
}
impl Bounds {
    fn empty() -> Self {
        Self {
            lo: [f64::INFINITY; 3],
            hi: [f64::NEG_INFINITY; 3],
        }
    }
    fn union(mut self, b: Self) -> Self {
        for i in 0..3 {
            self.lo[i] = self.lo[i].min(b.lo[i]);
            self.hi[i] = self.hi[i].max(b.hi[i]);
        }
        self
    }
    fn overlaps(self, b: Self, tol: f64) -> bool {
        (0..3).all(|i| self.lo[i] <= b.hi[i] + tol && b.lo[i] <= self.hi[i] + tol)
    }
}
struct Node {
    bounds: Bounds,
    range: std::ops::Range<usize>,
    children: Option<[usize; 2]>,
}
fn tree(ids: &mut [usize], offset: usize, boxes: &[Bounds], nodes: &mut Vec<Node>) -> usize {
    let bounds = ids.iter().fold(Bounds::empty(), |b, &i| b.union(boxes[i]));
    let index = nodes.len();
    nodes.push(Node {
        bounds,
        range: offset..offset + ids.len(),
        children: None,
    });
    if ids.len() > 8 {
        let axis = (0..3)
            .max_by(|&a, &b| {
                (bounds.hi[a] - bounds.lo[a]).total_cmp(&(bounds.hi[b] - bounds.lo[b]))
            })
            .unwrap();
        let mid = ids.len() / 2;
        ids.select_nth_unstable_by(mid, |&i, &j| {
            let b = boxes[i];
            let c = boxes[j];
            (b.lo[axis] + b.hi[axis]).total_cmp(&(c.lo[axis] + c.hi[axis]))
        });
        let (l, r) = ids.split_at_mut(mid);
        nodes[index].children = Some([
            tree(l, offset, boxes, nodes),
            tree(r, offset + mid, boxes, nodes),
        ]);
    }
    index
}

// Clip a coplanar segment against the triangle's three in-plane halfspaces.
// For a transverse segment, first reduce it to its intersection with the plane.
fn hits(mut p: Vec3, mut q: Vec3, tri: [Vec3; 3], normal: Vec3, tol: f64) -> bool {
    let a = (p - tri[0]).dot(normal);
    let b = (q - tri[0]).dot(normal);
    if a.min(b) > tol || a.max(b) < -tol {
        return false;
    }
    if a.abs().max(b.abs()) > tol {
        let t = (a / (a - b)).clamp(0., 1.);
        p = p + (q - p) * t;
        q = p;
    }
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    for i in 0..3 {
        let x = tri[i];
        let edge = tri[(i + 1) % 3] - x;
        let inward = normal.cross(edge / edge.length());
        let (dp, dq) = ((p - x).dot(inward) + tol, (q - x).dot(inward) + tol);
        if dp < 0. && dq < 0. {
            return false;
        }
        if dp < 0. {
            lo = lo.max(dp / (dp - dq));
        }
        if dq < 0. {
            hi = hi.min(dp / (dp - dq));
        }
    }
    lo <= hi
}
fn intersects(a: [usize; 3], b: [usize; 3], p: &[Vec3], tol: f64) -> bool {
    let shared = a.iter().filter(|v| b.contains(v)).count();
    let x = a.map(|i| p[i]);
    let y = b.map(|i| p[i]);
    let normal = |v: [Vec3; 3]| {
        let n = (v[1] - v[0]).cross(v[2] - v[0]);
        n / n.length()
    };
    let nx = normal(x);
    let ny = normal(y);
    if shared == 2 {
        // Distinct planes meet on the common edge. In a common plane the two
        // opposite vertices must lie on opposite sides of that edge.
        let other = b.iter().position(|v| !a.contains(v)).unwrap();
        if (y[other] - x[0]).dot(nx).abs() > tol {
            return false;
        }
        return nx.dot(ny) < 0.0;
    }
    for (nodes, other, vert, target, n) in [(a, b, x, y, ny), (b, a, y, x, nx)] {
        for i in 0..3 {
            let j = (i + 1) % 3;
            // For triangles sharing an apex any intersection extending beyond
            // it reaches the opposite edge of at least one of the triangles.
            if !other.contains(&nodes[i])
                && !other.contains(&nodes[j])
                && hits(vert[i], vert[j], target, n, tol)
            {
                return true;
            }
        }
    }
    false
}
pub(super) fn embedded(
    points: &[Vec3],
    triangles: &[[usize; 3]],
    scale: f64,
) -> Result<(), String> {
    let tol = 64.0 * f64::EPSILON * scale;
    let boxes: Vec<_> = triangles
        .iter()
        .map(|tri| {
            tri.iter().fold(Bounds::empty(), |b, &id| {
                let v = points[id];
                let coords = [v.x, v.y, v.z];
                b.union(Bounds {
                    lo: coords,
                    hi: coords,
                })
            })
        })
        .collect();
    let mut ids: Vec<_> = (0..triangles.len()).collect();
    let mut nodes = Vec::new();
    tree(&mut ids, 0, &boxes, &mut nodes);
    let mut stack = Vec::new();
    for (i, &tri) in triangles.iter().enumerate() {
        stack.push(0);
        while let Some(n) = stack.pop() {
            let node = &nodes[n];
            if !boxes[i].overlaps(node.bounds, tol) {
                continue;
            }
            if let Some(children) = node.children {
                stack.extend(children);
                continue;
            }
            for &j in &ids[node.range.clone()] {
                if j > i
                    && boxes[i].overlaps(boxes[j], tol)
                    && intersects(tri, triangles[j], points, tol)
                {
                    return Err(format!("surface collision between facets {i} and {j}"));
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rejects_transverse_coplanar_and_shared_vertex_overlaps() {
        let p = [
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0.2, 0.2, -1.],
            [0.2, 0.2, 1.],
            [0.4, 0.4, 0.],
            [0.6, 0.1, 0.],
            [0.7, 0.1, 0.],
        ]
        .map(|[x, y, z]| Vec3::new(x, y, z));
        for other in [[3, 4, 5], [5, 6, 7], [0, 5, 6], [1, 0, 5]] {
            assert!(embedded(&p, &[[0, 1, 2], other], 1.).is_err());
        }
    }
}
