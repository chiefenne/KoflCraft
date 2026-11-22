use std::time::Instant;

use mesh_algorithms::voronoi::validate_surface_sites;
use mesh_core::geometry::Vec3;
use mesh_geometry::{ConvexSurface, StlOptions, import_stl};
use serde::{Deserialize, Serialize};

use crate::artifact::{
    Point, SURFACE_ALGORITHM, generation_error, pack_artifact, point, require_schema, vector,
};
use crate::patches::{PatchMap, PatchSurface, roles_by_name};
use crate::sizing::{SizeFieldConfig, evidence};
use crate::{
    BoxBounds, EngineError, ErrorCode, MeshArtifact, MeshConfig, SCHEMA_VERSION,
    generate_surface_mesh,
};

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct StlMeshConfig {
    pub schema_version: u32,
    /// Explicit sites; leave empty when `size_field` places them.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub sites: Vec<[f64; 3]>,
    /// The size field that places the sites (ALG-05, GRD-01 to GRD-05).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub size_field: Option<SizeFieldConfig>,
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
    /// The role from the patch map (BC-04); without a map, none.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub role: Option<mesh_io::handoff::PatchRole>,
    pub facets: usize,
    pub area: f64,
}

/// Import and validate the surface and check the configured sites without building cells.
pub fn inspect_stl_domain(
    data: &[u8],
    config: &StlMeshConfig,
    facet_limit: Option<usize>,
    patches: Option<&PatchMap>,
) -> Result<StlDomainReport, EngineError> {
    let surface = import(data, config, facet_limit, patches)?;
    let roles = patches.map(|map| roles_by_name(map, surface.patch_names()));
    let site_count = match &config.size_field {
        Some(field) => field.count(&surface)?,
        None => {
            validate_surface_sites(&surface, &site_vectors(&config.sites))
                .map_err(|error| EngineError::new(ErrorCode::InvalidSites, error))?;
            config.sites.len()
        }
    };
    let scale = surface.bounds().extent().max_abs_component();
    let mut patches: Vec<_> = surface
        .patch_names()
        .iter()
        .enumerate()
        .map(|(i, name)| StlPatchReport {
            name: name.clone(),
            role: roles.as_ref().map(|r| r[i]),
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
        site_count,
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
///
/// With a patch map (BC-01), its patch names replace the STL solid names
/// and its roles are saved with the surface for the exports (BC-04).
pub fn build_stl_mesh(
    data: &[u8],
    config: StlMeshConfig,
    facet_limit: Option<usize>,
    patches: Option<&PatchMap>,
) -> Result<MeshArtifact, EngineError> {
    let surface = import(data, &config, facet_limit, patches)?;
    let roles = patches.map(|map| roles_by_name(map, surface.patch_names()));
    let started = Instant::now();
    let (sites, sizing) = match &config.size_field {
        Some(field_config) => {
            let seed = config.seed.unwrap_or(0);
            let (field, set) = field_config.seeds(&surface, seed)?;
            (set.sites.clone(), Some((field_config, field, set, seed)))
        }
        None => (site_vectors(&config.sites), None),
    };
    let generated = generate_surface_mesh(&surface, &sites).map_err(generation_error)?;
    let evidence = sizing.map(|(field_config, field, set, seed)| {
        evidence(
            field_config,
            &field,
            &set,
            seed,
            &surface,
            generated.mesh(),
            generated.metrics(),
        )
    });
    let config = MeshConfig {
        schema_version: config.schema_version,
        domain: bounds(&surface),
        sites: sites.into_iter().map(point).collect(),
        fixture_id: config.fixture_id,
        seed: config.seed,
    };
    pack_artifact(
        config,
        &generated,
        Some((&surface, roles)),
        evidence,
        started,
    )
}

/// The number of cells a build makes: the explicit sites, or the seed count
/// of the size field (exact: the same function the build calls).
pub fn stl_site_count(
    data: &[u8],
    config: &StlMeshConfig,
    facet_limit: Option<usize>,
) -> Result<usize, EngineError> {
    match &config.size_field {
        Some(field) => field.count(&import(data, config, facet_limit, None)?),
        None => Ok(config.sites.len()),
    }
}

fn import(
    data: &[u8],
    config: &StlMeshConfig,
    facet_limit: Option<usize>,
    patches: Option<&PatchMap>,
) -> Result<ConvexSurface, EngineError> {
    require_schema(config.schema_version)?;
    let names = match patches {
        None => None,
        Some(_) if config.facet_patches.is_some() => {
            return Err(EngineError::new(
                ErrorCode::InvalidPatchMap,
                "give either facet_patches or a patch map, not both",
            ));
        }
        Some(map) => Some(PatchSurface::import(data, facet_limit)?.facet_names(map)?),
    };
    if config.size_field.is_some() && !config.sites.is_empty() {
        return Err(EngineError::new(
            ErrorCode::InvalidSites,
            "give either sites or a size_field, not both",
        ));
    }
    import_stl(
        data,
        StlOptions {
            facet_patches: names.as_deref().or(config.facet_patches.as_deref()),
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
