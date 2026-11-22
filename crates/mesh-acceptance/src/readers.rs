//! Independent readers for the ASCII subsets used by the M1 fixtures.

use std::collections::{BTreeMap, BTreeSet};
use std::str::FromStr;

use mesh_core::geometry::Vec3;

use crate::{RunResult, require};

const TOLERANCE: f64 = 1e-10;
const PATCHES: [&str; 6] = ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"];

pub fn vtu(text: &str, expected: &[f64]) -> RunResult<usize> {
    let piece = between(text, "<Piece ", ">")?;
    let point_array = between(text, "<Points>", "</Points>")?;
    let points = read_points(array_values(point_array)?)?;
    require(
        xml_count(piece, "NumberOfPoints")? == points.len()
            && xml_count(piece, "NumberOfCells")? == expected.len(),
        "VTU piece counts disagree with the point array or cell count",
    )?;
    let cells = between(text, "<Cells>", "</Cells>")?;
    let array = |name| {
        let (_, tail) = cells
            .split_once(&format!("Name=\"{name}\""))
            .ok_or_else(|| format!("missing VTU array {name}"))?;
        numbers::<usize>(array_values(tail)?)
    };
    let connectivity = array("connectivity")?;
    let cell_offsets = array("offsets")?;
    let types = array("types")?;
    let faces = array("faces")?;
    let face_offsets = array("faceoffsets")?;
    require(
        types == vec![42; expected.len()],
        "VTU must contain one polyhedron per cell",
    )?;
    require(
        cell_offsets.len() == expected.len() && face_offsets.len() == expected.len(),
        "VTU offset count differs from cell count",
    )?;

    let mut previous_cell = 0;
    let mut previous_face = 0;
    let mut incidences = BTreeMap::<Vec<usize>, (Vec<usize>, usize)>::new();
    for (cell, &expected_volume) in expected.iter().enumerate() {
        let node_ids = connectivity
            .get(previous_cell..cell_offsets[cell])
            .ok_or("invalid VTU cell offset")?;
        let vertices: BTreeSet<_> = node_ids.iter().copied().collect();
        require(
            vertices.len() == node_ids.len() && vertices.len() >= 4,
            "invalid VTU cell vertex list",
        )?;
        let origin = vertices.iter().try_fold(Vec3::ZERO, |sum, &id| {
            Ok::<_, Box<dyn std::error::Error>>(sum + vertex(&points, id)?)
        })? / vertices.len() as f64;
        let stream = faces
            .get(previous_face..face_offsets[cell])
            .ok_or("invalid VTU face offset")?;
        let (&count, mut stream) = stream.split_first().ok_or("missing VTU face count")?;
        require(count >= 4, "VTU cell has too few faces")?;
        let mut edges = BTreeMap::new();
        let mut used = BTreeSet::new();
        let mut volume = 0.0;
        for _ in 0..count {
            let (&size, tail) = stream.split_first().ok_or("truncated VTU face")?;
            let (ring, tail) = tail
                .split_at_checked(size)
                .ok_or("truncated VTU vertices")?;
            stream = tail;
            require(
                ring.len() >= 3 && ring.iter().all(|id| vertices.contains(id)),
                "invalid VTU face vertices",
            )?;
            let mut key = ring.to_vec();
            key.sort_unstable();
            require(
                key.windows(2).all(|pair| pair[0] != pair[1]),
                "duplicate VTU face vertex",
            )?;
            let (first, uses) = incidences.entry(key).or_insert_with(|| (ring.to_vec(), 0));
            require(
                *uses == 0 || (*uses == 1 && opposite_rings(first, ring)),
                "invalid shared VTU face incidences",
            )?;
            *uses += 1;
            used.extend(ring.iter().copied());
            for (&a, &b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
                let counts = edges.entry((a.min(b), a.max(b))).or_insert([0, 0]);
                counts[usize::from(a > b)] += 1;
            }
            let area = face_area(&points, ring)?;
            let height = (vertex(&points, ring[0])? - origin).dot(area);
            require(height > 0.0, "VTU face is wound inward or degenerate")?;
            volume += height / 3.0;
        }
        require(
            stream.is_empty() && used == vertices,
            "incomplete VTU cell connectivity",
        )?;
        require(
            edges.values().all(|&uses| uses == [1, 1]),
            "VTU cell shell is not closed",
        )?;
        require(
            conserved(volume, expected_volume),
            &format!("VTU cell {cell} volume mismatch"),
        )?;
        previous_cell = cell_offsets[cell];
        previous_face = face_offsets[cell];
    }
    require(
        previous_cell == connectivity.len() && previous_face == faces.len(),
        "trailing VTU connectivity",
    )?;
    let mut patch_areas = [0.0; 6];
    for (ring, uses) in incidences.values() {
        if *uses == 1 {
            let (patch, area) = boundary_area(&points, ring)?;
            patch_areas[patch] += area;
        }
    }
    require(
        patch_areas.into_iter().all(|area| conserved(area, 1.0)),
        "incomplete VTU box boundary",
    )?;
    Ok(types.len())
}

pub fn gmsh(text: &str, expected: &[f64]) -> RunResult<usize> {
    let format = between(text, "$MeshFormat", "$EndMeshFormat")?.trim();
    require(
        format == "4.1 0 8" || format == "2.2 0 8",
        "expected Gmsh MSH 4.1 or 2.2 ASCII",
    )?;
    let names = rows(text, "$PhysicalNames", "$EndPhysicalNames")?;
    require(names.len() == 7, "unexpected Gmsh physical group count")?;
    for (patch, name) in PATCHES.into_iter().enumerate() {
        require(
            names
                .iter()
                .any(|line| line.trim() == format!("2 {} \"{name}\"", patch + 1)),
            "missing Gmsh boundary group",
        )?;
    }
    require(
        names.iter().any(|line| line.trim() == "3 7 \"fluid\""),
        "missing Gmsh fluid group",
    )?;
    let (points, records) = if format.starts_with("4.1") {
        msh41_records(text)?
    } else {
        msh22_records(text)?
    };
    let mut volumes = vec![0.0; expected.len()];
    let mut boundary = BTreeMap::new();
    let mut faces = BTreeMap::<[usize; 3], [u32; 2]>::new();
    let mut tetrahedra = 0;
    let mut patch_areas = [0.0; 6];
    for (index, record) in records.into_iter().enumerate() {
        require(
            record.first() == Some(&(index + 1)),
            "unexpected Gmsh element tag",
        )?;
        match record.as_slice() {
            [_, 2, 2, physical, entity, a, b, c] if physical == entity => {
                let ring = [*a, *b, *c];
                require(
                    ring.iter().all(|&id| id > 0),
                    "invalid Gmsh boundary node tag",
                )?;
                let (patch, area) = boundary_area(&points, &ring)?;
                require(*physical == patch + 1, "wrong Gmsh boundary patch")?;
                patch_areas[patch] += area;
                let (key, orientation) = triangle(ring);
                require(
                    boundary.insert(key, orientation).is_none(),
                    "duplicate Gmsh boundary triangle",
                )?;
            }
            [_, 4, 2, 7, cell, a, b, c, d] => {
                let ids = [*a, *b, *c, *d];
                require(ids.iter().all(|&id| id > 0), "invalid Gmsh node tag")?;
                let (p, q, r, s) = (
                    vertex(&points, *a)?,
                    vertex(&points, *b)?,
                    vertex(&points, *c)?,
                    vertex(&points, *d)?,
                );
                let volume = (q - p).dot((r - p).cross(s - p)) / 6.0;
                require(
                    volume > 0.0 && volume.is_finite(),
                    "nonpositive Gmsh tetrahedron",
                )?;
                *volumes
                    .get_mut(cell.checked_sub(1).ok_or("invalid Gmsh cell tag")?)
                    .ok_or("invalid Gmsh cell tag")? += volume;
                tetrahedra += 1;
                for ring in [[*b, *c, *d], [*a, *d, *c], [*a, *b, *d], [*a, *c, *b]] {
                    let (key, orientation) = triangle(ring);
                    faces.entry(key).or_default()[orientation] += 1;
                }
            }
            _ => return Err("invalid Gmsh element record".into()),
        }
    }
    for (face, uses) in faces {
        let expected_uses = match boundary.remove(&face) {
            Some(0) => [1, 0],
            Some(_) => [0, 1],
            None => [1, 1],
        };
        require(
            uses == expected_uses,
            "nonconforming Gmsh triangle incidences",
        )?;
    }
    require(
        boundary.is_empty() && patch_areas.into_iter().all(|area| conserved(area, 1.0)),
        "incomplete Gmsh box boundary",
    )?;
    for (cell, (volume, expected)) in volumes.into_iter().zip(expected).enumerate() {
        require(
            conserved(volume, *expected),
            &format!("Gmsh cell {cell} volume mismatch"),
        )?;
    }
    Ok(tetrahedra)
}

/// MSH 2.2: the points (index 0 unused, tags start at one) and the element
/// records `tag type 2 physical entity nodes…`.
fn msh22_records(text: &str) -> RunResult<(Vec<Vec3>, Vec<Vec<usize>>)> {
    let mut points = vec![Vec3::ZERO];
    for (index, line) in rows(text, "$Nodes", "$EndNodes")?.into_iter().enumerate() {
        let (tag, coordinates) = line
            .split_once(' ')
            .ok_or("missing Gmsh node coordinates")?;
        require(
            tag.parse::<usize>()? == index + 1,
            "unexpected Gmsh node tag",
        )?;
        let coordinates = read_points(coordinates)?;
        require(coordinates.len() == 1, "invalid Gmsh node record")?;
        points.push(coordinates[0]);
    }
    let records = rows(text, "$Elements", "$EndElements")?
        .into_iter()
        .map(numbers::<usize>)
        .collect::<RunResult<_>>()?;
    Ok((points, records))
}

/// MSH 4.1, read independently of the writer from the Gmsh reference manual:
/// the physical group of each entity from `$Entities`, the node blocks and
/// the element blocks, turned into the records of [`msh22_records`].
fn msh41_records(text: &str) -> RunResult<(Vec<Vec3>, Vec<Vec<usize>>)> {
    let mut lines = between(text, "$Entities", "$EndEntities")?.trim().lines();
    let head = numbers::<usize>(lines.next().ok_or("missing Gmsh entity counts")?)?;
    let [points_n, curves_n, surfaces_n, volumes_n] = head[..] else {
        return Err("invalid Gmsh entity counts".into());
    };
    require(
        points_n == 0 && curves_n == 0,
        "unexpected Gmsh point or curve entities",
    )?;
    let mut physical = BTreeMap::new();
    for (index, line) in lines.enumerate() {
        let words: Vec<&str> = line.split_whitespace().collect();
        require(words.len() >= 9, "invalid Gmsh entity record")?;
        let dim = if index < surfaces_n { 2 } else { 3 };
        require(words[7] == "1", "a Gmsh entity needs one physical group")?;
        physical.insert(
            (dim, words[0].parse::<usize>()?),
            words[8].parse::<usize>()?,
        );
    }
    require(
        physical.len() == surfaces_n + volumes_n,
        "Gmsh entity count mismatch",
    )?;

    let mut lines = between(text, "$Nodes", "$EndNodes")?.trim().lines();
    let head = numbers::<usize>(lines.next().ok_or("missing Gmsh node counts")?)?;
    let [blocks, total, ..] = head[..] else {
        return Err("invalid Gmsh node counts".into());
    };
    let mut points = vec![Vec3::ZERO; total + 1];
    let mut seen = 0;
    for _ in 0..blocks {
        let block = numbers::<usize>(lines.next().ok_or("missing Gmsh node block")?)?;
        let [_, _, 0, count] = block[..] else {
            return Err("invalid Gmsh node block".into());
        };
        let mut tags = Vec::with_capacity(count);
        for _ in 0..count {
            tags.push(
                lines
                    .next()
                    .ok_or("missing Gmsh node tag")?
                    .trim()
                    .parse::<usize>()?,
            );
        }
        for tag in tags {
            let coordinates = read_points(lines.next().ok_or("missing Gmsh node")?)?;
            require(coordinates.len() == 1, "invalid Gmsh node record")?;
            *points.get_mut(tag).ok_or("Gmsh node tag out of range")? = coordinates[0];
            seen += 1;
        }
    }
    require(seen == total, "Gmsh node count mismatch")?;

    let mut lines = between(text, "$Elements", "$EndElements")?.trim().lines();
    let head = numbers::<usize>(lines.next().ok_or("missing Gmsh element counts")?)?;
    let [blocks, total, ..] = head[..] else {
        return Err("invalid Gmsh element counts".into());
    };
    let mut records = Vec::with_capacity(total);
    for _ in 0..blocks {
        let block = numbers::<usize>(lines.next().ok_or("missing Gmsh element block")?)?;
        let [dim, entity, kind, count] = block[..] else {
            return Err("invalid Gmsh element block".into());
        };
        let group = *physical
            .get(&(dim, entity))
            .ok_or("Gmsh element block of an unknown entity")?;
        for _ in 0..count {
            let element = numbers::<usize>(lines.next().ok_or("missing Gmsh element")?)?;
            let (&tag, nodes) = element.split_first().ok_or("empty Gmsh element")?;
            let mut record = vec![tag, kind, 2, group, entity];
            record.extend_from_slice(nodes);
            records.push(record);
        }
    }
    require(records.len() == total, "Gmsh element count mismatch")?;
    records.sort_by_key(|r| r[0]);
    Ok((points, records))
}

fn opposite_rings(a: &[usize], b: &[usize]) -> bool {
    let Some(start) = a.iter().position(|id| id == &b[0]) else {
        return false;
    };
    b.iter().copied().eq(a[..=start]
        .iter()
        .rev()
        .chain(a[start + 1..].iter().rev())
        .copied())
}

fn triangle(mut ring: [usize; 3]) -> ([usize; 3], usize) {
    let [a, b, c] = ring;
    let orientation = usize::from((a > b) ^ (a > c) ^ (b > c));
    ring.sort_unstable();
    (ring, orientation)
}

fn face_area(points: &[Vec3], ring: &[usize]) -> RunResult<Vec3> {
    let a = vertex(points, ring[0])?;
    let mut area = Vec3::ZERO;
    for pair in ring[1..].windows(2) {
        area = area + (vertex(points, pair[0])? - a).cross(vertex(points, pair[1])? - a) * 0.5;
    }
    Ok(area)
}

fn boundary_area(points: &[Vec3], ring: &[usize]) -> RunResult<(usize, f64)> {
    let area = face_area(points, ring)?;
    for (axis, signed_area) in [area.x, area.y, area.z].into_iter().enumerate() {
        let high = signed_area > 0.0;
        let plane = if high { 1.0 } else { 0.0 };
        if signed_area != 0.0
            && ring.iter().all(|&id| {
                let p = points[id];
                ([p.x, p.y, p.z][axis] - plane).abs() <= TOLERANCE
            })
        {
            return Ok((2 * axis + usize::from(high), signed_area.abs()));
        }
    }
    Err("boundary face is off the box or wound inward".into())
}

fn vertex(points: &[Vec3], id: usize) -> RunResult<Vec3> {
    points
        .get(id)
        .copied()
        .ok_or_else(|| format!("unknown point {id}").into())
}

fn read_points(text: &str) -> RunResult<Vec<Vec3>> {
    let data = numbers::<f64>(text)?;
    let (points, remainder) = data.as_chunks::<3>();
    require(remainder.is_empty(), "incomplete point coordinates")?;
    require(
        data.iter()
            .all(|p| (-TOLERANCE..=1.0 + TOLERANCE).contains(p)),
        "point outside unit box",
    )?;
    Ok(points.iter().map(|&[x, y, z]| Vec3::new(x, y, z)).collect())
}

fn numbers<T: FromStr>(text: &str) -> RunResult<Vec<T>> {
    text.split_whitespace()
        .map(|s| s.parse().map_err(|_| format!("invalid number {s}").into()))
        .collect()
}

fn rows<'a>(text: &'a str, start: &str, end: &str) -> RunResult<Vec<&'a str>> {
    let mut lines = between(text, start, end)?.trim().lines();
    let count: usize = lines
        .next()
        .ok_or("missing MSH record count")?
        .trim()
        .parse()?;
    let rows: Vec<_> = lines.collect();
    require(rows.len() == count, "MSH record count mismatch")?;
    Ok(rows)
}

fn between<'a>(text: &'a str, start: &str, end: &str) -> RunResult<&'a str> {
    let (_, tail) = text
        .split_once(start)
        .ok_or_else(|| format!("missing {start}"))?;
    tail.split_once(end)
        .map(|(body, _)| body)
        .ok_or_else(|| format!("missing {end}").into())
}

fn array_values(text: &str) -> RunResult<&str> {
    between(text, ">", "</DataArray>")
}

fn xml_count(tag: &str, name: &str) -> RunResult<usize> {
    Ok(between(tag, &format!("{name}=\""), "\"")?.parse()?)
}

fn conserved(actual: f64, expected: f64) -> bool {
    expected > 0.0 && (actual / expected - 1.0).abs() <= TOLERANCE
}
