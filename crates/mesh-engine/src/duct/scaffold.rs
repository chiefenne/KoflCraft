//! Read linear, physically tagged tetrahedral scaffolds in Gmsh MSH 2.2 ASCII.
use mesh_algorithms::duct::Scaffold;
use mesh_core::geometry::Vec3;
use std::collections::HashMap;

pub(super) fn read(bytes: &[u8]) -> Result<Scaffold, String> {
    let text = std::str::from_utf8(bytes).map_err(|_| "expected ASCII Gmsh MSH 2.2")?;
    let section = |name: &str| {
        let mut lines = text.lines();
        lines
            .find(|l| l.trim() == name)
            .ok_or_else(|| format!("missing {name}"))?;
        Ok::<_, String>(lines)
    };
    if section("$MeshFormat")?.next().map(str::trim) != Some("2.2 0 8") {
        return Err("expected ASCII Gmsh MSH 2.2".into());
    }
    let mut lines = section("$Nodes")?;
    let count = |value: Option<&str>| {
        value
            .and_then(|v| v.trim().parse::<usize>().ok())
            .ok_or("invalid MSH count")
    };
    let n = count(lines.next())?;
    let mut points = Vec::new();
    let mut indices = HashMap::new();
    for _ in 0..n {
        let mut fields = lines.next().ok_or("truncated nodes")?.split_whitespace();
        let tag = count(fields.next())?;
        let mut xyz = [0.0; 3];
        for x in &mut xyz {
            *x = fields
                .next()
                .and_then(|v| v.parse::<f64>().ok())
                .filter(|v| v.is_finite())
                .ok_or("invalid node coordinate")?;
        }
        if fields.next().is_some() || tag == 0 || indices.insert(tag, points.len()).is_some() {
            return Err("invalid or duplicate node record".into());
        }
        points.push(Vec3::new(xyz[0], xyz[1], xyz[2]));
    }
    if lines.next().map(str::trim) != Some("$EndNodes") {
        return Err("invalid node section".into());
    }
    let mut lines = section("$Elements")?;
    let n = count(lines.next())?;
    let mut boundary = Vec::new();
    let mut tetrahedra = Vec::new();
    for _ in 0..n {
        let mut fields = lines.next().ok_or("truncated elements")?.split_whitespace();
        let _id = count(fields.next())?;
        let kind = count(fields.next())?;
        let tags = count(fields.next())?;
        if !matches!(kind, 2 | 4) {
            return Err(
                "scaffold must contain only linear boundary triangles and volume tetrahedra".into(),
            );
        }
        let tag = count(fields.next())?;
        if tags == 0 {
            return Err("missing physical tag".into());
        }
        for _ in 1..tags {
            fields.next().ok_or("truncated element tags")?;
        }
        let mut node = || {
            indices
                .get(&count(fields.next())?)
                .copied()
                .ok_or("unknown element node")
        };
        match kind {
            2 if (1..=3).contains(&tag) => boundary.push(([node()?, node()?, node()?], tag)),
            4 if tag == 4 => tetrahedra.push([node()?, node()?, node()?, node()?]),
            _ => return Err("expected physical IDs 1 inlet, 2 outlet, 3 walls, 4 fluid".into()),
        }
        if fields.next().is_some() {
            return Err("unexpected element node".into());
        }
    }
    if lines.next().map(str::trim) != Some("$EndElements")
        || tetrahedra.is_empty()
        || boundary.is_empty()
    {
        return Err("incomplete scaffold".into());
    }
    Ok(Scaffold {
        points,
        boundary,
        tetrahedra,
    })
}
