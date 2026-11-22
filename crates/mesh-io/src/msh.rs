//! One element model for every Gmsh writer, serialized as MSH 4.1 (the
//! default, the version FoehnCraft reads, its D33) or as MSH 2.2 (kept as an
//! option, OPEN-11).
//!
//! Written from the public description of the format in the Gmsh reference
//! manual, section "MSH file format". The writers of this crate collect
//! nodes, physical groups and elements into a [`Model`]; only the
//! serialization differs between the versions.
//!
//! In MSH 4.1 every element belongs to an elementary entity. The model keeps
//! the elementary tag each writer gives an element (the original cell index
//! plus one for the box, the region for the duct) as the entity, and each
//! entity carries the one physical group of its elements. All nodes are
//! classified on the first volume entity; KoflCraft writes no geometry, so
//! the entities have no bounding entities.

use std::collections::BTreeMap;
use std::io::{self, Write};

use mesh_core::geometry::Vec3;

/// The MSH version a writer serializes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum GmshVersion {
    /// MSH 4.1 ASCII, the default (FC-01).
    #[default]
    V4_1,
    /// MSH 2.2 ASCII, kept as an option (OPEN-11).
    V2_2,
}

impl GmshVersion {
    /// The version as the `$MeshFormat` section names it.
    pub fn label(self) -> &'static str {
        match self {
            Self::V4_1 => "4.1",
            Self::V2_2 => "2.2",
        }
    }
}

/// Gmsh element types written by KoflCraft and their dimension and node count.
fn element_shape(kind: u8) -> Option<(u8, usize)> {
    match kind {
        2 => Some((2, 3)), // triangle
        3 => Some((2, 4)), // quadrangle
        4 => Some((3, 4)), // tetrahedron
        5 => Some((3, 8)), // hexahedron
        6 => Some((3, 6)), // prism
        7 => Some((3, 5)), // pyramid
        _ => None,
    }
}

struct Element {
    kind: u8,
    physical: usize,
    entity: usize,
    /// Range in `Model::nodes`.
    start: usize,
    /// Parent cell of a piece of a subdivided polyhedron (`ParentCellId`).
    parent: Option<usize>,
}

/// Nodes, physical groups and elements of one mesh.
pub(crate) struct Model {
    points: Vec<Vec3>,
    groups: Vec<(u8, usize, String)>,
    elements: Vec<Element>,
    nodes: Vec<usize>,
}

fn invalid(message: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}

impl Model {
    pub(crate) fn new(points: Vec<Vec3>) -> Self {
        Self {
            points,
            groups: Vec::new(),
            elements: Vec::new(),
            nodes: Vec::new(),
        }
    }

    /// Add a physical group of dimension 2 (boundary) or 3 (volume).
    pub(crate) fn group(&mut self, dim: u8, tag: usize, name: &str) {
        self.groups.push((dim, tag, name.to_string()));
    }

    /// Add an element with zero-based node ids.
    pub(crate) fn element(
        &mut self,
        kind: u8,
        physical: usize,
        entity: usize,
        nodes: &[usize],
    ) -> io::Result<()> {
        self.element_of(kind, physical, entity, nodes, None)
    }

    /// Add a piece of the polyhedron `parent` with zero-based node ids.
    pub(crate) fn element_of(
        &mut self,
        kind: u8,
        physical: usize,
        entity: usize,
        nodes: &[usize],
        parent: Option<usize>,
    ) -> io::Result<()> {
        let (_, count) =
            element_shape(kind).ok_or_else(|| invalid(format!("element type {kind}")))?;
        if nodes.len() != count {
            return Err(invalid(format!(
                "element type {kind} needs {count} nodes, got {}",
                nodes.len()
            )));
        }
        if let Some(&node) = nodes.iter().find(|&&n| n >= self.points.len()) {
            return Err(invalid(format!("element node {node} is not a point")));
        }
        if physical == 0 || entity == 0 {
            return Err(invalid("Gmsh tags start at one".to_string()));
        }
        self.elements.push(Element {
            kind,
            physical,
            entity,
            start: self.nodes.len(),
            parent,
        });
        self.nodes.extend_from_slice(nodes);
        Ok(())
    }

    fn nodes_of(&self, element: &Element) -> &[usize] {
        let count = element_shape(element.kind).map_or(0, |(_, n)| n);
        &self.nodes[element.start..element.start + count]
    }

    fn has_parents(&self) -> bool {
        self.elements.iter().any(|e| e.parent.is_some())
    }

    pub(crate) fn write(&self, version: GmshVersion, out: &mut dyn Write) -> io::Result<()> {
        match version {
            GmshVersion::V2_2 => self.write_v2(out),
            GmshVersion::V4_1 => self.write_v4(out),
        }
    }

    fn write_names(&self, out: &mut dyn Write) -> io::Result<()> {
        writeln!(out, "$PhysicalNames\n{}", self.groups.len())?;
        for (dim, tag, name) in &self.groups {
            writeln!(out, "{dim} {tag} \"{name}\"")?;
        }
        writeln!(out, "$EndPhysicalNames")
    }

    fn write_point(out: &mut dyn Write, p: Vec3) -> io::Result<()> {
        writeln!(out, "{:e} {:e} {:e}", p.x, p.y, p.z)
    }

    fn write_parents(
        &self,
        out: &mut dyn Write,
        tags: impl Iterator<Item = (usize, usize)>,
    ) -> io::Result<()> {
        let rows: Vec<_> = tags.collect();
        writeln!(
            out,
            "$ElementData\n1\n\"ParentCellId\"\n1\n0\n3\n0\n1\n{}",
            rows.len()
        )?;
        for (tag, parent) in rows {
            writeln!(out, "{tag} {parent}")?;
        }
        writeln!(out, "$EndElementData")
    }

    fn write_v2(&self, out: &mut dyn Write) -> io::Result<()> {
        writeln!(out, "$MeshFormat\n2.2 0 8\n$EndMeshFormat")?;
        self.write_names(out)?;
        writeln!(out, "$Nodes\n{}", self.points.len())?;
        for (index, &point) in self.points.iter().enumerate() {
            write!(out, "{} ", index + 1)?;
            Self::write_point(out, point)?;
        }
        writeln!(out, "$EndNodes\n$Elements\n{}", self.elements.len())?;
        for (index, e) in self.elements.iter().enumerate() {
            write!(
                out,
                "{} {} 2 {} {}",
                index + 1,
                e.kind,
                e.physical,
                e.entity
            )?;
            for node in self.nodes_of(e) {
                write!(out, " {}", node + 1)?;
            }
            writeln!(out)?;
        }
        writeln!(out, "$EndElements")?;
        if self.has_parents() {
            let tags = self
                .elements
                .iter()
                .enumerate()
                .filter_map(|(i, e)| e.parent.map(|p| (i + 1, p)));
            self.write_parents(out, tags)?;
        }
        Ok(())
    }

    fn write_v4(&self, out: &mut dyn Write) -> io::Result<()> {
        // Entities: (dim, tag) -> physical group and bounding box.
        let mut entities: BTreeMap<(u8, usize), (usize, Vec3, Vec3)> = BTreeMap::new();
        for e in &self.elements {
            let (dim, _) = element_shape(e.kind).unwrap_or((0, 0));
            let mut low = Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
            let mut high = Vec3::new(-f64::INFINITY, -f64::INFINITY, -f64::INFINITY);
            for &n in self.nodes_of(e) {
                let p = self.points[n];
                low = Vec3::new(low.x.min(p.x), low.y.min(p.y), low.z.min(p.z));
                high = Vec3::new(high.x.max(p.x), high.y.max(p.y), high.z.max(p.z));
            }
            let entry = entities
                .entry((dim, e.entity))
                .or_insert((e.physical, low, high));
            if entry.0 != e.physical {
                return Err(invalid(format!(
                    "entity {} of dimension {dim} has elements of the physical groups {} and {}",
                    e.entity, entry.0, e.physical
                )));
            }
            entry.1 = Vec3::new(
                entry.1.x.min(low.x),
                entry.1.y.min(low.y),
                entry.1.z.min(low.z),
            );
            entry.2 = Vec3::new(
                entry.2.x.max(high.x),
                entry.2.y.max(high.y),
                entry.2.z.max(high.z),
            );
        }
        let Some(&(_, node_entity)) = entities.keys().find(|(dim, _)| *dim == 3) else {
            return Err(invalid("the mesh has no volume elements".to_string()));
        };
        let count = |d: u8| entities.keys().filter(|(dim, _)| *dim == d).count();

        writeln!(out, "$MeshFormat\n4.1 0 8\n$EndMeshFormat")?;
        self.write_names(out)?;
        writeln!(out, "$Entities\n0 0 {} {}", count(2), count(3))?;
        for (&(_, tag), &(physical, low, high)) in &entities {
            writeln!(
                out,
                "{tag} {:e} {:e} {:e} {:e} {:e} {:e} 1 {physical} 0",
                low.x, low.y, low.z, high.x, high.y, high.z
            )?;
        }
        writeln!(out, "$EndEntities")?;

        let n = self.points.len();
        writeln!(out, "$Nodes\n1 {n} 1 {n}\n3 {node_entity} 0 {n}")?;
        for tag in 1..=n {
            writeln!(out, "{tag}")?;
        }
        for &point in &self.points {
            Self::write_point(out, point)?;
        }
        writeln!(out, "$EndNodes")?;

        // Elements in blocks of one entity and type. Tags follow the block
        // order, so they rise through the file.
        let mut blocks: BTreeMap<(u8, usize, u8), Vec<usize>> = BTreeMap::new();
        for (index, e) in self.elements.iter().enumerate() {
            let (dim, _) = element_shape(e.kind).unwrap_or((0, 0));
            blocks
                .entry((dim, e.entity, e.kind))
                .or_default()
                .push(index);
        }
        let total = self.elements.len();
        writeln!(out, "$Elements\n{} {total} 1 {total}", blocks.len())?;
        let mut tag = 0;
        let mut parents = Vec::new();
        for (&(dim, entity, kind), members) in &blocks {
            writeln!(out, "{dim} {entity} {kind} {}", members.len())?;
            for &index in members {
                tag += 1;
                let e = &self.elements[index];
                write!(out, "{tag}")?;
                for node in self.nodes_of(e) {
                    write!(out, " {}", node + 1)?;
                }
                writeln!(out)?;
                if let Some(parent) = e.parent {
                    parents.push((tag, parent));
                }
            }
        }
        writeln!(out, "$EndElements")?;
        if !parents.is_empty() {
            self.write_parents(out, parents.into_iter())?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Read either version back into (node coordinates, physical names,
    /// sorted element records `(kind, physical, entity, nodes, parent)`).
    type Record = (u8, usize, usize, Vec<usize>, Option<usize>);

    fn section<'a>(text: &'a str, name: &str) -> Option<&'a str> {
        let start = text.find(&format!("${name}\n"))? + name.len() + 2;
        let end = text.find(&format!("$End{name}"))?;
        Some(&text[start..end])
    }

    fn numbers(line: &str) -> Vec<f64> {
        line.split_whitespace()
            .map(|w| w.parse().unwrap())
            .collect()
    }

    fn read(text: &str) -> (Vec<[f64; 3]>, String, Vec<Record>) {
        let format = section(text, "MeshFormat").unwrap().trim().to_string();
        let names = section(text, "PhysicalNames").unwrap().to_string();
        let mut points = Vec::new();
        let mut elements = Vec::new();
        let mut parent_of = BTreeMap::new();
        if let Some(data) = section(text, "ElementData") {
            for line in data.lines().skip(9) {
                let v = numbers(line);
                parent_of.insert(v[0] as usize, v[1] as usize);
            }
        }
        if format.starts_with("2.2") {
            for line in section(text, "Nodes").unwrap().lines().skip(1) {
                let v = numbers(line);
                assert_eq!(v[0] as usize, points.len() + 1);
                points.push([v[1], v[2], v[3]]);
            }
            for line in section(text, "Elements").unwrap().lines().skip(1) {
                let v: Vec<usize> = numbers(line).into_iter().map(|x| x as usize).collect();
                assert_eq!(v[2], 2);
                let nodes = v[5..].iter().map(|n| n - 1).collect();
                elements.push((v[1] as u8, v[3], v[4], nodes, parent_of.get(&v[0]).copied()));
            }
        } else {
            assert_eq!(format, "4.1 0 8");
            let mut physical = BTreeMap::new();
            let entities = section(text, "Entities").unwrap();
            let mut lines = entities.lines();
            let head = numbers(lines.next().unwrap());
            let surfaces = head[2] as usize;
            for (i, line) in lines.enumerate() {
                let v = numbers(line);
                let dim = if i < surfaces { 2 } else { 3 };
                assert_eq!(v[7], 1.0);
                physical.insert((dim, v[0] as usize), v[8] as usize);
            }
            let nodes = section(text, "Nodes").unwrap();
            let mut lines = nodes.lines();
            let head = numbers(lines.next().unwrap());
            let n = head[1] as usize;
            let block = numbers(lines.next().unwrap());
            assert_eq!(block[3] as usize, n);
            let rows: Vec<_> = lines.collect();
            for (i, tag) in rows[..n].iter().enumerate() {
                assert_eq!(tag.trim().parse::<usize>().unwrap(), i + 1);
            }
            for line in &rows[n..] {
                let v = numbers(line);
                points.push([v[0], v[1], v[2]]);
            }
            let body = section(text, "Elements").unwrap();
            let mut lines = body.lines();
            let head = numbers(lines.next().unwrap());
            let mut expected_tag = 1;
            for _ in 0..head[0] as usize {
                let b = numbers(lines.next().unwrap());
                let (dim, entity, kind) = (b[0] as u8, b[1] as usize, b[2] as u8);
                for _ in 0..b[3] as usize {
                    let v: Vec<usize> = numbers(lines.next().unwrap())
                        .into_iter()
                        .map(|x| x as usize)
                        .collect();
                    assert_eq!(v[0], expected_tag, "element tags rise");
                    expected_tag += 1;
                    let nodes = v[1..].iter().map(|n| n - 1).collect();
                    let phys = physical[&(dim, entity)];
                    elements.push((kind, phys, entity, nodes, parent_of.get(&v[0]).copied()));
                }
            }
            assert_eq!(expected_tag - 1, head[1] as usize);
        }
        elements.sort();
        (points, names, elements)
    }

    fn sample() -> Model {
        let points = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 1.0, 1.0),
        ];
        let mut m = Model::new(points);
        m.group(2, 1, "inlet");
        m.group(2, 2, "walls");
        m.group(3, 3, "fluid");
        m.element(2, 1, 1, &[0, 2, 1]).unwrap();
        m.element(2, 2, 2, &[0, 1, 3]).unwrap();
        m.element_of(4, 3, 1, &[0, 1, 2, 3], Some(0)).unwrap();
        m.element(2, 1, 1, &[1, 2, 4]).unwrap();
        m.element_of(4, 3, 2, &[1, 2, 3, 4], Some(1)).unwrap();
        m
    }

    /// Invariant: both versions hold the same nodes, groups, elements with
    /// their physical group and entity, and parent cells. A change to one
    /// serializer alone (lost element, wrong tag, reordered parent data)
    /// breaks it; the per-format tests elsewhere read one version each, so
    /// only this test ties the two together.
    #[test]
    fn both_versions_carry_the_same_model() {
        let m = sample();
        let mut v2 = Vec::new();
        let mut v4 = Vec::new();
        m.write(GmshVersion::V2_2, &mut v2).unwrap();
        m.write(GmshVersion::V4_1, &mut v4).unwrap();
        let a = read(std::str::from_utf8(&v2).unwrap());
        let b = read(std::str::from_utf8(&v4).unwrap());
        assert_eq!(a, b);
        assert_eq!(a.2.len(), 5);
    }

    /// Invariant: an entity names one physical group, as MSH 4.1 requires.
    /// Mixing them would silently put elements into the wrong group.
    #[test]
    fn an_entity_with_two_physical_groups_is_refused() {
        let mut m = sample();
        m.element(2, 2, 1, &[0, 1, 4]).unwrap();
        let mut out = Vec::new();
        assert!(m.write(GmshVersion::V4_1, &mut out).is_err());
    }
}
