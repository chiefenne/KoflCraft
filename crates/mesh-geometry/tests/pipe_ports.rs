use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

use mesh_core::geometry::Vec3;
use mesh_geometry::{ClosedSurface, ConvexSurface, StlOptions, SurfaceError, import_closed_stl};
use serde::Deserialize;
use sha2::{Digest, Sha256};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Fixture {
    schema_version: u32,
    source_sha256: String,
    facet_count: usize,
    vertex_count: usize,
    bounds: Bounds,
    volume: f64,
    default_patch: Patch,
    ports: [Port; 2],
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Bounds {
    min: [f64; 3],
    max: [f64; 3],
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
enum Role {
    Wall,
    Inlet,
    Outlet,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Patch {
    name: String,
    role: Role,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Port {
    name: String,
    role: Role,
    facet_range: [usize; 2],
    center: [f64; 3],
    outward_normal: [f64; 3],
    radius: f64,
    area: f64,
    rim_edges: usize,
}

#[test]
fn manifold_has_two_exterior_ports_and_feature_walls() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test-fixtures/manifold");
    let bytes = fs::read(path.join("manifold.stl")).unwrap();
    let fixture: Fixture =
        serde_json::from_slice(&fs::read(path.join("patches.json")).unwrap()).unwrap();
    assert_eq!(fixture.schema_version, 1);
    let sha256: String = Sha256::digest(&bytes)
        .into_iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    assert_eq!(sha256, fixture.source_sha256);
    assert_eq!(fixture.default_patch.role, Role::Wall);
    assert_eq!(
        fixture.ports.each_ref().map(|port| port.role),
        [Role::Inlet, Role::Outlet]
    );

    let mut names = vec![fixture.default_patch.name.clone(); fixture.facet_count];
    for port in &fixture.ports {
        let [start, end] = port.facet_range;
        assert!(start < end, "{} facet range is empty", port.name);
        for patch in &mut names[start..end] {
            assert_eq!(patch, &fixture.default_patch.name, "port ranges overlap");
            patch.clone_from(&port.name);
        }
    }
    let shell = import_closed_stl(
        &bytes,
        StlOptions {
            facet_patches: Some(&names),
            facet_limit: None,
        },
    )
    .unwrap();
    assert_eq!(shell.facets().len(), fixture.facet_count);
    assert_eq!(shell.vertices().len(), fixture.vertex_count);
    assert_eq!(shell.bounds().min(), vector(fixture.bounds.min));
    assert_eq!(shell.bounds().max(), vector(fixture.bounds.max));
    assert!((shell.volume() - fixture.volume).abs() < 1e-6);
    let expected_names: BTreeSet<_> = std::iter::once(fixture.default_patch.name.as_str())
        .chain(fixture.ports.iter().map(|port| port.name.as_str()))
        .collect();
    assert_eq!(
        shell
            .patch_names()
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        expected_names
    );

    for port in &fixture.ports {
        assert_port_facets(&shell, port);
    }
    assert_port_rims_and_wall_features(&shell, &fixture);
    assert!(matches!(
        ConvexSurface::try_from(shell),
        Err(SurfaceError::NonConvexFacet(_))
    ));
}

fn assert_port_facets(shell: &ClosedSurface, port: &Port) {
    let normal = vector(port.outward_normal);
    let center = vector(port.center);
    assert_eq!(normal.length(), 1.0);
    let [start, end] = port.facet_range;
    let scale = shell.bounds().extent().max_abs_component();
    let mut area = 0.0;
    for (index, facet) in shell.facets().iter().enumerate() {
        let is_port = shell.patch_names()[facet.triangle.patch.0] == port.name;
        let points = facet.triangle.vertices.map(|id| shell.vertices()[id]);
        let on_plane = points
            .iter()
            .all(|&point| (point - center).dot(normal) == 0.0);
        assert_eq!(is_port, (start..end).contains(&index));
        assert_eq!(
            on_plane, is_port,
            "{} has an unclassified or off-plane facet",
            port.name
        );
        if is_port {
            assert!(facet.plane.normal().dot(normal) > 1.0 - 1e-12);
            assert!(
                points
                    .iter()
                    .all(|&point| (point - center).length() <= port.radius + 2e-5)
            );
            area += facet.normalized_area * scale * scale;
        }
    }
    assert!((area - port.area).abs() < 1e-8, "{} area {area}", port.name);
}

fn assert_port_rims_and_wall_features(shell: &ClosedSurface, fixture: &Fixture) {
    let mut pending: BTreeMap<_, usize> = BTreeMap::new();
    let mut rims = [Vec::new(), Vec::new()];
    let mut sharp_wall_edges = 0;
    let mut curved_wall_edges = 0;
    let wall = &fixture.default_patch.name;
    let curved_range = 15_f64.to_radians().cos()..0.25_f64.to_radians().cos();
    let sharp_cosine = 45_f64.to_radians().cos();
    for (index, facet) in shell.facets().iter().enumerate() {
        let [a, b, c] = facet.triangle.vertices;
        let name = &shell.patch_names()[facet.triangle.patch.0];
        for (a, b) in [(a, b), (b, c), (c, a)] {
            let edge = (a.min(b), a.max(b));
            let Some(other) = pending.remove(&edge) else {
                pending.insert(edge, index);
                continue;
            };
            let other = &shell.facets()[other];
            let other_name = &shell.patch_names()[other.triangle.patch.0];
            if name == wall && other_name == wall {
                let cosine = facet.plane.normal().dot(other.plane.normal());
                sharp_wall_edges += usize::from(cosine < sharp_cosine);
                curved_wall_edges += usize::from(curved_range.contains(&cosine));
            } else if name != other_name {
                assert!(
                    name == wall || other_name == wall,
                    "the two ports share an edge"
                );
                let port_name = if name == wall { other_name } else { name };
                let port = fixture
                    .ports
                    .iter()
                    .position(|port| &port.name == port_name)
                    .unwrap();
                rims[port].push(edge);
            }
        }
    }
    assert!(pending.is_empty());
    assert!(sharp_wall_edges > 100, "expected sharp edges on the walls");
    assert!(
        curved_wall_edges > 1000,
        "expected finely faceted roundings on the walls"
    );
    for (rim, port) in rims.into_iter().zip(&fixture.ports) {
        assert_eq!(rim.len(), port.rim_edges, "{} rim edge count", port.name);
        let center = vector(port.center);
        let mut neighbours = BTreeMap::<usize, Vec<usize>>::new();
        for (a, b) in rim {
            neighbours.entry(a).or_default().push(b);
            neighbours.entry(b).or_default().push(a);
        }
        for (&vertex, incident) in &neighbours {
            assert_eq!(incident.len(), 2, "{} rim is not a loop", port.name);
            let radius = (shell.vertices()[vertex] - center).length();
            assert!(
                (radius - port.radius).abs() < 2e-5,
                "{} rim radius {radius}",
                port.name
            );
        }
        let first = *neighbours.keys().next().unwrap();
        let mut seen = BTreeSet::new();
        let mut pending = vec![first];
        while let Some(vertex) = pending.pop() {
            if seen.insert(vertex) {
                pending.extend(&neighbours[&vertex]);
            }
        }
        assert_eq!(
            seen.len(),
            neighbours.len(),
            "{} rim has multiple loops",
            port.name
        );
    }
}

fn vector([x, y, z]: [f64; 3]) -> Vec3 {
    Vec3::new(x, y, z)
}
