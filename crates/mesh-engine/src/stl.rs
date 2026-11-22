use std::time::Instant;

use mesh_algorithms::voronoi::validate_surface_sites;
use mesh_core::geometry::Vec3;
use mesh_geometry::{ConvexSurface, StlOptions, import_stl};
use serde::{Deserialize, Serialize};

use crate::artifact::{
    Point, SURFACE_ALGORITHM, generation_error, pack_artifact, point, require_schema, vector,
};
use crate::{
    BoxBounds, EngineError, ErrorCode, MeshArtifact, MeshConfig, SCHEMA_VERSION,
    generate_surface_mesh,
};

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct StlMeshConfig {
    pub schema_version: u32,
    pub sites: Vec<[f64; 3]>,
    pub fixture_id: Option<String>,
    pub seed: Option<u64>,
    /// Optional override with one name per triangle in the input file's order.
    pub facet_patches: Option<Vec<String>>,
}

#[derive(Debug, Serialize)]
pub struct StlDomainReport {
    pub schema_version: u32,
    pub bounds: BoxBounds,
    pub site_count: usize,
    pub facet_count: usize,
    pub vertex_count: usize,
    pub volume: f64,
    pub patches: Vec<StlPatchReport>,
    pub fixture_id: Option<String>,
    pub seed: Option<u64>,
    pub algorithm: &'static str,
}

#[derive(Debug, Serialize)]
pub struct StlPatchReport {
    pub name: String,
    pub facets: usize,
    pub area: f64,
}

/// Import and validate the surface and check the configured sites without building cells.
pub fn inspect_stl_domain(
    data: &[u8],
    config: &StlMeshConfig,
    facet_limit: Option<usize>,
) -> Result<StlDomainReport, EngineError> {
    let surface = import(data, config, facet_limit)?;
    validate_surface_sites(&surface, &site_vectors(&config.sites))
        .map_err(|error| EngineError::new(ErrorCode::InvalidSites, error))?;
    let scale = surface.bounds().extent().max_abs_component();
    let mut patches: Vec<_> = surface
        .patch_names()
        .iter()
        .map(|name| StlPatchReport {
            name: name.clone(),
            facets: 0,
            area: 0.0,
        })
        .collect();
    for facet in surface.facets() {
        let patch = &mut patches[facet.triangle.patch.0];
        patch.facets += 1;
        patch.area += facet.normalized_area * (scale * scale);
    }
    Ok(StlDomainReport {
        schema_version: SCHEMA_VERSION,
        bounds: bounds(&surface),
        site_count: config.sites.len(),
        facet_count: surface.facets().len(),
        vertex_count: surface.vertices().len(),
        volume: surface.volume(),
        patches,
        fixture_id: config.fixture_id.clone(),
        seed: config.seed,
        algorithm: SURFACE_ALGORITHM,
    })
}

/// Build an artifact containing both the certified volume mesh and its validated surface.
pub fn build_stl_mesh(
    data: &[u8],
    config: StlMeshConfig,
    facet_limit: Option<usize>,
) -> Result<MeshArtifact, EngineError> {
    let surface = import(data, &config, facet_limit)?;
    let started = Instant::now();
    let generated =
        generate_surface_mesh(&surface, &site_vectors(&config.sites)).map_err(generation_error)?;
    let config = MeshConfig {
        schema_version: config.schema_version,
        domain: bounds(&surface),
        sites: config.sites,
        fixture_id: config.fixture_id,
        seed: config.seed,
    };
    pack_artifact(config, &generated, Some(&surface), started)
}

fn import(
    data: &[u8],
    config: &StlMeshConfig,
    facet_limit: Option<usize>,
) -> Result<ConvexSurface, EngineError> {
    require_schema(config.schema_version)?;
    import_stl(
        data,
        StlOptions {
            facet_patches: config.facet_patches.as_deref(),
            facet_limit,
        },
    )
    .map_err(|error| EngineError::new(ErrorCode::InvalidDomain, error))
}

fn bounds(surface: &ConvexSurface) -> BoxBounds {
    let bounds = surface.bounds();
    BoxBounds {
        min: point(bounds.min()),
        max: point(bounds.max()),
    }
}

fn site_vectors(sites: &[Point]) -> Vec<Vec3> {
    sites.iter().copied().map(vector).collect()
}
