use std::collections::BTreeMap;

use mesh_core::geometry::{PatchId, Vec3};

use crate::{ClosedSurface, ConvexSurface, SurfaceError, SurfaceTriangle};

#[derive(Clone, Copy, Debug, Default)]
pub struct StlOptions<'a> {
    /// Optional patch name per facet in file order.
    pub facet_patches: Option<&'a [String]>,
    pub facet_limit: Option<usize>,
}

/// Read a convex STL shell. Geometry is oriented by facet vertex winding.
pub fn import_stl(bytes: &[u8], options: StlOptions<'_>) -> Result<ConvexSurface, SurfaceError> {
    import_closed_stl(bytes, options)?.try_into()
}

/// Read and check a closed STL shell for surface analysis, including nonconvex inputs.
pub fn import_closed_stl(
    bytes: &[u8],
    options: StlOptions<'_>,
) -> Result<ClosedSurface, SurfaceError> {
    let mut surface = IndexedSurface::default();
    let binary_size = bytes
        .get(80..84)
        .map(|count| u32::from_le_bytes(count.try_into().unwrap()) as usize)
        .and_then(|count| count.checked_mul(50)?.checked_add(84));
    if binary_size == Some(bytes.len()) {
        let (records, remainder) = bytes[84..].as_chunks::<50>();
        debug_assert!(remainder.is_empty());
        for record in records {
            let triangle = std::array::from_fn(|vertex| {
                let [x, y, z] = std::array::from_fn(|axis| {
                    let offset = 12 + 12 * vertex + 4 * axis;
                    f64::from(f32::from_le_bytes(
                        record[offset..offset + 4].try_into().unwrap(),
                    ))
                });
                Vec3::new(x, y, z)
            });
            surface.add(triangle, "surface", options)?;
        }
    } else {
        let text = std::str::from_utf8(bytes)
            .map_err(|_| syntax(None, "expected ASCII STL or a matching binary facet count"))?;
        parse_ascii(text.trim_start_matches('\u{feff}'), &mut surface, options)?;
    }
    if let Some(names) = options.facet_patches
        && names.len() != surface.triangles.len()
    {
        return Err(SurfaceError::PatchCount {
            facets: surface.triangles.len(),
            names: names.len(),
        });
    }
    ClosedSurface::new(surface.vertices, surface.triangles, surface.patch_names)
}

#[derive(Default)]
struct IndexedSurface {
    vertices: Vec<Vec3>,
    triangles: Vec<SurfaceTriangle>,
    patch_names: Vec<String>,
    vertex_ids: BTreeMap<[u64; 3], usize>,
    patch_ids: BTreeMap<String, PatchId>,
}

impl IndexedSurface {
    fn add(
        &mut self,
        points: [Vec3; 3],
        name: &str,
        options: StlOptions<'_>,
    ) -> Result<(), SurfaceError> {
        let index = self.triangles.len();
        if let Some(limit) = options.facet_limit
            && index >= limit
        {
            return Err(SurfaceError::FacetLimit { limit });
        }
        let name = options
            .facet_patches
            .and_then(|names| names.get(index))
            .map_or(name, String::as_str);
        let patch = *self.patch_ids.entry(name.to_owned()).or_insert_with(|| {
            let id = PatchId(self.patch_names.len());
            self.patch_names.push(name.to_owned());
            id
        });
        let vertices = points.map(|p| {
            let floats = [p.x, p.y, p.z].map(|v| if v == 0.0 { 0.0 } else { v });
            let key = floats.map(f64::to_bits);
            *self.vertex_ids.entry(key).or_insert_with(|| {
                let index = self.vertices.len();
                self.vertices
                    .push(Vec3::new(floats[0], floats[1], floats[2]));
                index
            })
        });
        self.triangles.push(SurfaceTriangle { vertices, patch });
        Ok(())
    }
}

fn parse_ascii(
    text: &str,
    surface: &mut IndexedSurface,
    options: StlOptions<'_>,
) -> Result<(), SurfaceError> {
    let mut lines = text.lines().enumerate().filter_map(|(index, line)| {
        let line = line.trim();
        (!line.is_empty()).then_some((index + 1, line))
    });
    while let Some((line, text)) = lines.next() {
        let name = keyword(text, "solid").ok_or_else(|| syntax(Some(line), "expected solid"))?;
        let patch = if name.is_empty() { "surface" } else { name };
        loop {
            let (line, text) = next_line(&mut lines)?;
            if let Some(end_name) = keyword(text, "endsolid") {
                if !end_name.is_empty() && end_name != name {
                    return Err(syntax(Some(line), "endsolid name differs from solid name"));
                }
                break;
            }
            let normal = keyword(text, "facet")
                .and_then(|text| keyword(text, "normal"))
                .ok_or_else(|| syntax(Some(line), "expected facet normal or endsolid"))?;
            coordinates(line, normal)?;
            expect_line(&mut lines, "outer loop")?;
            let mut vertices = [Vec3::ZERO; 3];
            for vertex in &mut vertices {
                let (line, text) = next_line(&mut lines)?;
                let coords =
                    keyword(text, "vertex").ok_or_else(|| syntax(Some(line), "expected vertex"))?;
                *vertex = coordinates(line, coords)?;
            }
            expect_line(&mut lines, "endloop")?;
            expect_line(&mut lines, "endfacet")?;
            surface.add(vertices, patch, options)?;
        }
    }
    Ok(())
}

fn keyword<'a>(line: &'a str, token: &str) -> Option<&'a str> {
    let rest = line.strip_prefix(token)?;
    (rest.is_empty() || rest.starts_with(char::is_whitespace)).then(|| rest.trim())
}

fn coordinates(line: usize, text: &str) -> Result<Vec3, SurfaceError> {
    let error = syntax(Some(line), "expected three finite coordinates");
    let mut coordinates = text.split_whitespace();
    let mut next = || {
        coordinates
            .next()
            .and_then(|value| value.parse::<f64>().ok())
            .filter(|v| v.is_finite())
            .ok_or_else(|| error.clone())
    };
    let point = Vec3::new(next()?, next()?, next()?);
    if coordinates.next().is_some() {
        return Err(error);
    }
    Ok(point)
}

fn next_line<'a>(
    lines: &mut impl Iterator<Item = (usize, &'a str)>,
) -> Result<(usize, &'a str), SurfaceError> {
    lines
        .next()
        .ok_or_else(|| syntax(None, "unexpected end of ASCII surface"))
}

fn expect_line<'a>(
    lines: &mut impl Iterator<Item = (usize, &'a str)>,
    expected: &'static str,
) -> Result<(), SurfaceError> {
    let (line, text) = next_line(lines)?;
    if !text.split_whitespace().eq(expected.split_whitespace()) {
        return Err(syntax(Some(line), expected));
    }
    Ok(())
}

fn syntax(line: Option<usize>, reason: &'static str) -> SurfaceError {
    SurfaceError::InvalidStl { line, reason }
}
