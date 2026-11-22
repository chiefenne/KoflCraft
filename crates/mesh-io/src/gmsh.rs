use std::io::Write;

use mesh_core::geometry::{AxisAlignedBox, BoxSide, Vec3};
use mesh_core::mesh::{FaceAdjacency, Mesh};
use mesh_quality::{
    MeshDomain, MeshMetrics, RELATIVE_CONSERVATION_TOLERANCE, RELATIVE_GEOMETRY_TOLERANCE,
    validate_domain,
};

use crate::ExportError;
use crate::msh::{GmshVersion, Model};

/// Write a Gmsh MSH tetrahedral subdivision of the validated mesh.
/// Each volume element's elementary entity tag is its original cell index plus one.
pub fn write_gmsh(
    domain: AxisAlignedBox,
    mesh: &Mesh,
    version: GmshVersion,
    out: &mut dyn Write,
) -> Result<(), ExportError> {
    write_gmsh_domain(MeshDomain::Box(domain), mesh, version, out)
}

pub fn write_gmsh_domain(
    domain: MeshDomain<'_>,
    mesh: &Mesh,
    version: GmshVersion,
    out: &mut dyn Write,
) -> Result<(), ExportError> {
    let metrics = validate_domain(domain, mesh).map_err(ExportError::Validation)?;
    let patch_names: Vec<_> = match domain {
        MeshDomain::Box(_) => BoxSide::ALL.into_iter().map(BoxSide::name).collect(),
        MeshDomain::Surface(surface) => surface.patch_names().iter().map(String::as_str).collect(),
    };
    let mut points = mesh.vertices().to_vec();
    let scale = domain.bounds().extent().max_abs_component();
    let triangles = triangulate_faces(mesh, &metrics, &mut points, scale);
    let cell_centres = points.len();
    points.extend(metrics.cells.iter().map(|cell| cell.centroid));

    let volume_scale = scale * scale * scale;
    let mut tetrahedra = Vec::new();
    let mut boundary = Vec::new();
    for cell in 0..mesh.cell_count() {
        let centre = cell_centres + cell;
        let mut volume = 0.0;
        for &face in mesh.cell_faces(cell).unwrap() {
            let adjacency = mesh.face_adjacency(face).unwrap();
            // Both incidences use the same triangles with opposite winding.
            for &[a, mut b, mut c] in &triangles[face] {
                if adjacency.owner() != cell {
                    std::mem::swap(&mut b, &mut c);
                }
                let ids = [centre, a, b, c];
                let [p, q, r, s] = ids.map(|id| points[id]);
                let tet_volume =
                    ((q - p) / scale).dot(((r - p) / scale).cross((s - p) / scale)) / 6.0;
                if !tet_volume.is_finite() || tet_volume <= 0.0 {
                    return Err(ExportError::InvalidSubdivision {
                        cell,
                        face: Some(face),
                    });
                }
                volume += tet_volume;
                tetrahedra.push((cell, ids));
                match adjacency {
                    FaceAdjacency::Boundary { patch, .. } => {
                        boundary.push((patch as usize, [a, b, c]))
                    }
                    FaceAdjacency::Surface { patch, .. } => boundary.push((patch.0, [a, b, c])),
                    FaceAdjacency::Interior { .. } => {}
                }
            }
        }
        let Some(metric) = metrics.cells.get(cell) else {
            return Err(ExportError::InvalidSubdivision { cell, face: None });
        };
        let expected = metric.volume / volume_scale;
        if (volume / expected - 1.0).abs() > RELATIVE_CONSERVATION_TOLERANCE {
            return Err(ExportError::InvalidSubdivision { cell, face: None });
        }
    }
    let fluid = patch_names.len() + 1;
    let mut model = Model::new(points);
    for (patch, name) in patch_names.iter().enumerate() {
        model.group(2, patch + 1, name);
    }
    model.group(3, fluid, "fluid");
    for (patch, vertices) in boundary {
        model
            .element(2, patch + 1, patch + 1, &vertices)
            .map_err(ExportError::Io)?;
    }
    for (cell, vertices) in tetrahedra {
        model
            .element(4, fluid, cell + 1, &vertices)
            .map_err(ExportError::Io)?;
    }
    model.write(version, out).map_err(ExportError::Io)
}

fn triangulate_faces(
    mesh: &Mesh,
    metrics: &MeshMetrics,
    points: &mut Vec<Vec3>,
    scale: f64,
) -> Vec<Vec<[usize; 3]>> {
    let mut all = Vec::with_capacity(mesh.face_count());
    for (face, metrics) in metrics.faces.iter().enumerate() {
        let vertices = mesh.face_vertices(face).unwrap();
        let first = vertices[0];
        let origin = points[first];
        let triangles = vertices[1..].windows(2);
        if triangles.clone().all(|pair| {
            let u = (points[pair[0]] - origin) / scale;
            let v = (points[pair[1]] - origin) / scale;
            let longest_edge = u.length().max(v.length()).max((v - u).length());
            u.cross(v).dot(metrics.normal) > RELATIVE_GEOMETRY_TOLERANCE * longest_edge
        }) {
            all.push(triangles.map(|pair| [first, pair[0], pair[1]]).collect());
        } else {
            // A vertex fan can contain flat triangles when a ring has collinear vertices.
            let centroid = points.len();
            points.push(metrics.centroid);
            all.push(
                vertices
                    .iter()
                    .zip(vertices.iter().cycle().skip(1))
                    .map(|(&a, &b)| [centroid, a, b])
                    .collect(),
            );
        }
    }
    all
}
