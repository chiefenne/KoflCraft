pub mod block;
pub mod channel;
pub mod clipping;
pub mod delaunay;
pub mod duct;
pub mod neighbours;
pub mod predicates;
pub mod remesh;
pub mod scaffold;
pub mod sizing;
pub mod sweep;
pub mod voronoi;

/// Progress lines on standard error when the environment variable
/// `KOFL_TRACE` is set (counts and lengths only, no coordinates).
pub fn trace_enabled() -> bool {
    std::env::var_os("KOFL_TRACE").is_some()
}
