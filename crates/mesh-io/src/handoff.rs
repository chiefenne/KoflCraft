//! The face-based mesh KoflCraft hands to FoehnCraft (requirements section
//! 16): points, faces with owner and neighbour, named boundary patches with
//! a role, and named cell zones. Every mesh type of KoflCraft (box, STL,
//! channel, sweep, duct) converts into it, and the FoehnCraft writers (the
//! single-file VTU of FC-02 and the OpenFOAM `polyMesh` of FC-03) and the
//! round-trip check of FC-06 work on it alone.
//!
//! Faces are stored in OpenFOAM order: the interior faces first, ordered by
//! owner and then neighbour with the owner the lower cell, then the
//! boundary faces patch by patch. A face is wound outward from its owner.

use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::io::{self, Write};

use mesh_core::geometry::Vec3;
use mesh_core::mesh::{FaceAdjacency, Mesh};
use serde::{Deserialize, Serialize};

/// The role of a boundary patch, mapped to a FoehnCraft boundary condition
/// in `docs/foehncraft.md` (FC-04, GEO-06).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PatchRole {
    Wall,
    Inlet,
    Outlet,
    Symmetry,
    Interface,
    Unspecified,
}

impl PatchRole {
    /// The role a patch name suggests: `inlet`, `outlet`, `wall`,
    /// `symmetry` anywhere in the lower-case name, else unspecified.
    pub fn from_name(name: &str) -> Self {
        let n = name.to_ascii_lowercase();
        if n.contains("inlet") {
            Self::Inlet
        } else if n.contains("outlet") {
            Self::Outlet
        } else if n.contains("wall") {
            Self::Wall
        } else if n.contains("symmetry") {
            Self::Symmetry
        } else if n.contains("interface") {
            Self::Interface
        } else {
            Self::Unspecified
        }
    }

    /// The FoehnCraft boundary-condition type a case would start from.
    pub fn foehncraft_type(self) -> Option<&'static str> {
        match self {
            Self::Wall => Some("wall"),
            Self::Inlet => Some("velocity-inlet"),
            Self::Outlet => Some("static-pressure-outlet"),
            Self::Symmetry => Some("symmetry"),
            Self::Interface | Self::Unspecified => None,
        }
    }

    /// The OpenFOAM patch type of the `boundary` file.
    pub fn openfoam_type(self) -> &'static str {
        match self {
            Self::Wall => "wall",
            Self::Symmetry => "symmetry",
            _ => "patch",
        }
    }
}

/// A named boundary patch.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Patch {
    /// The name in the files: a TOML bare key, as FoehnCraft names
    /// boundaries in `[boundary.<name>]` (see [`handoff_name`]).
    pub name: String,
    /// The name the patch had in KoflCraft's input, when it differs.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_name: Option<String>,
    pub role: PatchRole,
}

/// The name FoehnCraft gets for a KoflCraft patch or zone name: ASCII
/// letters, digits, `_` and `-` are kept, every other character becomes
/// `_`; an empty name becomes `patch`.
pub fn handoff_name(name: &str) -> String {
    let n: String = name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' || c == '-' {
                c
            } else {
                '_'
            }
        })
        .collect();
    if n.is_empty() { "patch".to_string() } else { n }
}

/// Patches from KoflCraft names: hand-off names made unique by a suffix,
/// roles from the names.
pub fn patches_from_names<S: AsRef<str>>(names: &[S]) -> Vec<Patch> {
    let mut used: HashMap<String, usize> = HashMap::new();
    names
        .iter()
        .map(|source| {
            let source = source.as_ref();
            let mut name = handoff_name(source);
            let n = used.entry(name.clone()).or_insert(0);
            *n += 1;
            if *n > 1 {
                name = format!("{name}_{n}");
            }
            Patch {
                source_name: (name != source).then(|| source.to_string()),
                role: PatchRole::from_name(source),
                name,
            }
        })
        .collect()
}

/// Why a set of cells is not a conforming face-based mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct HandoffError(pub String);

impl fmt::Display for HandoffError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for HandoffError {}

fn fail<T>(message: impl Into<String>) -> Result<T, HandoffError> {
    Err(HandoffError(message.into()))
}

/// One cell as its outward face rings, and its zone.
pub struct CellFaces {
    pub faces: Vec<Vec<usize>>,
    pub zone: usize,
}

impl CellFaces {
    /// Reverse every ring if the cell's signed volume is negative, for
    /// elements given in a node order whose orientation is not known.
    pub fn outward(mut self, points: &[Vec3]) -> Self {
        let origin = self
            .faces
            .first()
            .and_then(|f| f.first())
            .map_or(Vec3::ZERO, |&v| points[v]);
        let v6: f64 = self
            .faces
            .iter()
            .map(|r| ring_volume6(points, r, origin))
            .sum();
        if v6 < 0.0 {
            for r in &mut self.faces {
                r.reverse();
            }
        }
        self
    }
}

/// The face-based hand-off mesh.
#[derive(Clone, Debug, PartialEq)]
pub struct HandoffMesh {
    pub points: Vec<Vec3>,
    face_offsets: Vec<usize>,
    face_nodes: Vec<usize>,
    /// Owner of every face.
    pub owner: Vec<usize>,
    /// Neighbour of every interior face (the first `neighbour.len()` faces).
    pub neighbour: Vec<usize>,
    /// First face of each patch, and one past the last face.
    patch_starts: Vec<usize>,
    pub patches: Vec<Patch>,
    pub cell_zone: Vec<usize>,
    pub zones: Vec<String>,
    /// Faces of KoflCraft's mesh that were merged into their neighbours
    /// because two cells shared more than one face (see [`merge_shared`]).
    pub merged_faces: usize,
}

/// The faces of every cell as one list (see [`HandoffMesh::cell_face_list`]):
/// a few bytes per face instead of a vector per cell, which matters at
/// millions of cells.
pub struct CellFaceList {
    start: Vec<usize>,
    /// `2 f` for a face the cell owns (wound outward), `2 f + 1` for one
    /// it is the neighbour of (wound inward).
    faces: Vec<usize>,
}

impl CellFaceList {
    /// The faces of cell `c`: (face, wound outward from `c`).
    pub fn cell(&self, c: usize) -> impl Iterator<Item = (usize, bool)> + '_ {
        self.faces[self.start[c]..self.start[c + 1]]
            .iter()
            .map(|&k| (k / 2, k % 2 == 0))
    }

    /// The number of faces of cell `c`.
    pub fn len(&self, c: usize) -> usize {
        self.start[c + 1] - self.start[c]
    }
}

/// Counts and conserved sums FoehnCraft must reproduce (FC-06).
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct HandoffSummary {
    pub cells: usize,
    pub faces: usize,
    pub interior_faces: usize,
    pub boundary_faces: usize,
    pub points: usize,
    pub total_volume: f64,
    pub min_volume: f64,
    /// Largest |Σ n·dS| of a cell relative to its surface area.
    pub max_relative_closure: f64,
    /// Interior faces merged so that two cells share at most one face.
    pub merged_faces: usize,
    /// Interior faces with S·d ≤ 0 (see [`HandoffMesh::inverted_faces`]).
    pub inverted_faces: usize,
    pub patches: Vec<PatchSummary>,
    pub zones: Vec<ZoneSummary>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct PatchSummary {
    pub name: String,
    pub role: PatchRole,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub foehncraft_type: Option<&'static str>,
    pub faces: usize,
    pub area: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ZoneSummary {
    pub name: String,
    pub cells: usize,
    pub volume: f64,
}

/// Vector area of a closed polygon (exact for any ring, planar or not).
pub(crate) fn ring_area(points: &[Vec3], ring: &[usize]) -> Vec3 {
    let mut a = Vec3::ZERO;
    for (i, &v) in ring.iter().enumerate() {
        let w = ring[(i + 1) % ring.len()];
        a = a + points[v].cross(points[w]);
    }
    a * 0.5
}

/// Six times the signed volume a ring contributes, by the fan of triangles
/// from the ring's vertex average, relative to `origin`.
fn ring_volume6(points: &[Vec3], ring: &[usize], origin: Vec3) -> f64 {
    let mut c = Vec3::ZERO;
    for &v in ring {
        c = c + points[v];
    }
    let c = c * (1.0 / ring.len() as f64) - origin;
    let mut v6 = 0.0;
    for (i, &a) in ring.iter().enumerate() {
        let b = ring[(i + 1) % ring.len()];
        v6 += c.dot((points[a] - origin).cross(points[b] - origin));
    }
    v6
}

/// FoehnCraft (like most finite-volume solvers) allows one face between two
/// cells. Where KoflCraft's mesh has several (the pieces of a non-planar
/// interface of the polyhedral dual), their union becomes one polygon: the
/// edges shared by two pieces cancel and the rest must form one loop.
/// Anything else (a ring with a hole, two separate patches) is an error.
/// Input sorted by cell pair; returns the faces and how many were merged
/// away.
#[allow(clippy::type_complexity)]
fn merge_shared(
    interior: Vec<(usize, usize, Vec<usize>)>,
) -> Result<(Vec<(usize, usize, Vec<usize>)>, usize), HandoffError> {
    let mut out: Vec<(usize, usize, Vec<usize>)> = Vec::with_capacity(interior.len());
    let mut merged = 0;
    let mut i = 0;
    while i < interior.len() {
        let mut j = i + 1;
        while j < interior.len() && (interior[j].0, interior[j].1) == (interior[i].0, interior[i].1)
        {
            j += 1;
        }
        if j == i + 1 {
            out.push(interior[i].clone());
        } else {
            let (lo, hi) = (interior[i].0, interior[i].1);
            let mut edges: HashMap<(usize, usize), usize> = HashMap::new();
            for (_, _, ring) in &interior[i..j] {
                for (k, &a) in ring.iter().enumerate() {
                    let b = ring[(k + 1) % ring.len()];
                    if let Some(n) = edges.get_mut(&(b, a)) {
                        *n -= 1;
                        if *n == 0 {
                            edges.remove(&(b, a));
                        }
                    } else {
                        *edges.entry((a, b)).or_insert(0) += 1;
                    }
                }
            }
            let mut next: HashMap<usize, usize> = HashMap::new();
            for (&(a, b), &n) in &edges {
                if n != 1 || next.insert(a, b).is_some() {
                    return fail(format!(
                        "cells {lo} and {hi} share {} faces that do not merge into one polygon",
                        j - i
                    ));
                }
            }
            let Some(&start) = next.keys().min() else {
                return fail(format!("cells {lo} and {hi} share faces that cancel out"));
            };
            let mut ring = vec![start];
            let mut at = next[&start];
            while at != start {
                if ring.len() > next.len() {
                    break;
                }
                ring.push(at);
                at = next[&at];
            }
            if ring.len() != next.len() {
                return fail(format!(
                    "cells {lo} and {hi} share {} faces whose union is not one polygon",
                    j - i
                ));
            }
            merged += j - i - 1;
            out.push((lo, hi, ring));
        }
        i = j;
    }
    Ok((out, merged))
}

fn key(ring: &[usize]) -> Vec<usize> {
    let mut k = ring.to_vec();
    k.sort_unstable();
    k
}

impl HandoffMesh {
    /// Build from cells given as outward face rings. Two cells that share a
    /// face must list the same vertices (in opposite winding); a face of one
    /// cell must be one of `boundary` (ring and patch index). Anything else,
    /// a face shared by three cells, a boundary face no cell has, is an
    /// error: the result is always conforming and fully classified.
    pub fn from_cells(
        points: Vec<Vec3>,
        cells: Vec<CellFaces>,
        zones: Vec<String>,
        boundary: &[(Vec<usize>, usize)],
        patches: Vec<Patch>,
    ) -> Result<Self, HandoffError> {
        let mut boundary_patch: HashMap<Vec<usize>, usize> = HashMap::new();
        for (ring, patch) in boundary {
            if *patch >= patches.len() {
                return fail(format!("boundary face of the unknown patch {patch}"));
            }
            if boundary_patch.insert(key(ring), *patch).is_some() {
                return fail("a boundary face is listed twice");
            }
        }
        let mut open: HashMap<Vec<usize>, (usize, usize)> = HashMap::new();
        let mut interior: Vec<(usize, usize, Vec<usize>)> = Vec::new();
        for (c, cell) in cells.iter().enumerate() {
            if cell.zone >= zones.len() {
                return fail(format!("cell {c} has the unknown zone {}", cell.zone));
            }
            if cell.faces.len() < 4 {
                return fail(format!("cell {c} has {} faces", cell.faces.len()));
            }
            for (f, ring) in cell.faces.iter().enumerate() {
                if ring.len() < 3 || ring.iter().any(|&v| v >= points.len()) {
                    return fail(format!("cell {c} face {f} is not a polygon of the points"));
                }
                match open.remove(&key(ring)) {
                    None => {
                        open.insert(key(ring), (c, f));
                    }
                    Some((o, of)) => {
                        let first = &cells[o].faces[of];
                        let (lo, hi, ring) = if o < c {
                            (o, c, first.clone())
                        } else {
                            (c, o, ring.clone())
                        };
                        if lo == hi {
                            return fail(format!("cell {c} has a face twice"));
                        }
                        interior.push((lo, hi, ring));
                    }
                }
            }
        }
        // A face seen a third time reopens; it is caught here as an
        // unmatched face that is no boundary face.
        let mut bfaces: Vec<(usize, usize, Vec<usize>)> = Vec::with_capacity(open.len());
        for (k, (c, f)) in open {
            let Some(patch) = boundary_patch.remove(&k) else {
                return fail(format!(
                    "cell {c} face {f} has no neighbour and is no boundary face"
                ));
            };
            bfaces.push((patch, c, cells[c].faces[f].clone()));
        }
        if !boundary_patch.is_empty() {
            return fail(format!(
                "{} boundary faces belong to no cell",
                boundary_patch.len()
            ));
        }
        interior.sort_by_key(|a| (a.0, a.1));
        let (interior, merged_faces) = merge_shared(interior)?;
        bfaces.sort_by(|a, b| (a.0, a.1, &a.2).cmp(&(b.0, b.1, &b.2)));

        let mut mesh = Self {
            points,
            face_offsets: vec![0],
            face_nodes: Vec::new(),
            owner: Vec::new(),
            neighbour: Vec::new(),
            patch_starts: Vec::new(),
            patches,
            cell_zone: cells.iter().map(|c| c.zone).collect(),
            zones,
            merged_faces,
        };
        for (lo, hi, ring) in interior {
            mesh.push_face(&ring, lo);
            mesh.neighbour.push(hi);
        }
        let mut next = 0;
        for (patch, c, ring) in bfaces {
            while next <= patch {
                mesh.patch_starts.push(mesh.owner.len());
                next += 1;
            }
            mesh.push_face(&ring, c);
        }
        while next <= mesh.patches.len() {
            mesh.patch_starts.push(mesh.owner.len());
            next += 1;
        }
        mesh.check_cells()?;
        Ok(mesh)
    }

    /// Build from KoflCraft's canonical mesh (box and STL meshes, the
    /// channel): every face from its owner, and boundary faces by patch
    /// index (`BoxSide` order for boxes, the surface patch for STL).
    pub fn from_mesh(mesh: &Mesh, patches: Vec<Patch>) -> Result<Self, HandoffError> {
        Self::from_mesh_with(mesh, patches, &|adjacency| match adjacency {
            FaceAdjacency::Boundary { patch, .. } => Some(patch as usize),
            FaceAdjacency::Surface { patch, .. } => Some(patch.0),
            FaceAdjacency::Interior { .. } => None,
        })
    }

    /// As [`Self::from_mesh`], with the patch index of a boundary face
    /// given by `patch_of`.
    pub fn from_mesh_with(
        mesh: &Mesh,
        patches: Vec<Patch>,
        patch_of: &dyn Fn(FaceAdjacency) -> Option<usize>,
    ) -> Result<Self, HandoffError> {
        let mut cells: Vec<CellFaces> = (0..mesh.cell_count())
            .map(|_| CellFaces {
                faces: Vec::new(),
                zone: 0,
            })
            .collect();
        let mut boundary = Vec::new();
        for face in 0..mesh.face_count() {
            let (Some(ring), Some(adjacency)) =
                (mesh.face_vertices(face), mesh.face_adjacency(face))
            else {
                return fail(format!("face {face} is missing"));
            };
            cells[adjacency.owner()].faces.push(ring.to_vec());
            if let FaceAdjacency::Interior { neighbour, .. } = adjacency {
                cells[neighbour]
                    .faces
                    .push(ring.iter().rev().copied().collect());
            } else {
                let patch = patch_of(adjacency)
                    .ok_or_else(|| HandoffError(format!("face {face} has no patch")))?;
                boundary.push((ring.to_vec(), patch));
            }
        }
        Self::from_cells(
            mesh.vertices().to_vec(),
            cells,
            vec!["fluid".to_string()],
            &boundary,
            patches,
        )
    }

    fn push_face(&mut self, ring: &[usize], owner: usize) {
        self.face_nodes.extend_from_slice(ring);
        self.face_offsets.push(self.face_nodes.len());
        self.owner.push(owner);
    }

    pub fn cell_count(&self) -> usize {
        self.cell_zone.len()
    }

    pub fn face_count(&self) -> usize {
        self.owner.len()
    }

    /// The vertices of a face, wound outward from its owner.
    pub fn face(&self, face: usize) -> &[usize] {
        &self.face_nodes[self.face_offsets[face]..self.face_offsets[face + 1]]
    }

    /// The faces of patch `p`.
    pub fn patch_faces(&self, p: usize) -> std::ops::Range<usize> {
        self.patch_starts[p]..self.patch_starts[p + 1]
    }

    /// The cells as lists of (face, outward from this cell).
    pub fn cell_faces(&self) -> Vec<Vec<(usize, bool)>> {
        let list = self.cell_face_list();
        (0..self.cell_count())
            .map(|c| list.cell(c).collect())
            .collect()
    }

    /// The faces of every cell in one flat list, in the order of
    /// [`Self::cell_faces`] (owned faces, then neighbour faces, each by
    /// face number), without a vector per cell.
    pub fn cell_face_list(&self) -> CellFaceList {
        let mut start = vec![0usize; self.cell_count() + 1];
        for &c in self.owner.iter().chain(&self.neighbour) {
            start[c + 1] += 1;
        }
        for c in 0..self.cell_count() {
            start[c + 1] += start[c];
        }
        let mut next = start.clone();
        let mut faces = vec![0usize; start[self.cell_count()]];
        for (f, &o) in self.owner.iter().enumerate() {
            faces[next[o]] = 2 * f;
            next[o] += 1;
        }
        for (f, &n) in self.neighbour.iter().enumerate() {
            faces[next[n]] = 2 * f + 1;
            next[n] += 1;
        }
        CellFaceList { start, faces }
    }

    /// Per cell: volume, |Σ n dS|, Σ |dS| and centroid, all from the
    /// tetrahedra (cell origin, face node average, edge), the decomposition
    /// FoehnCraft uses too (its `geometry` module).
    fn cell_measures(&self) -> Vec<(f64, f64, f64, Vec3)> {
        let cells = self.cell_face_list();
        (0..self.cell_count())
            .map(|c| {
                let first = cells.cell(c).next().map_or(0, |(f, _)| f);
                let origin = self.points[self.face(first)[0]];
                let mut v6 = 0.0;
                let mut moment = Vec3::ZERO;
                let mut sum = Vec3::ZERO;
                let mut total = 0.0;
                for (f, out) in cells.cell(c) {
                    let s = if out { 1.0 } else { -1.0 };
                    let ring = self.face(f);
                    let mut c = Vec3::ZERO;
                    for &v in ring {
                        c = c + self.points[v];
                    }
                    let c = c * (1.0 / ring.len() as f64) - origin;
                    for (k, &a) in ring.iter().enumerate() {
                        let b = ring[(k + 1) % ring.len()];
                        let (pa, pb) = (self.points[a] - origin, self.points[b] - origin);
                        let t6 = s * c.dot(pa.cross(pb));
                        v6 += t6;
                        moment = moment + (c + pa + pb) * (t6 / 4.0);
                    }
                    let a = ring_area(&self.points, ring);
                    sum = sum + a * s;
                    total += a.length();
                }
                let centroid = origin + moment * (1.0 / v6);
                (v6 / 6.0, sum.length(), total, centroid)
            })
            .collect()
    }

    /// Interior faces whose area vector does not point from the owner's
    /// centroid towards the neighbour's (S·d ≤ 0). FoehnCraft counts them
    /// as inverted and calls the mesh invalid.
    pub fn inverted_faces(&self) -> usize {
        let m = self.cell_measures();
        (0..self.neighbour.len())
            .filter(|&f| {
                let d = m[self.neighbour[f]].3 - m[self.owner[f]].3;
                ring_area(&self.points, self.face(f)).dot(d) <= 0.0
            })
            .count()
    }

    fn check_cells(&self) -> Result<(), HandoffError> {
        for (c, (v, closure, total, _)) in self.cell_measures().into_iter().enumerate() {
            if !(v > 0.0 && v.is_finite()) {
                return fail(format!("cell {c} has the volume {v:e}"));
            }
            if closure.is_nan() || closure > 1e-9 * total {
                return fail(format!("cell {c} is not closed (|Σ n dS| = {closure:e})"));
            }
        }
        Ok(())
    }

    /// Counts, volumes and areas.
    pub fn summary(&self) -> HandoffSummary {
        let measures = self.cell_measures();
        let mut zones: Vec<ZoneSummary> = self
            .zones
            .iter()
            .map(|n| ZoneSummary {
                name: n.clone(),
                cells: 0,
                volume: 0.0,
            })
            .collect();
        for (c, &(v, _, _, _)) in measures.iter().enumerate() {
            zones[self.cell_zone[c]].cells += 1;
            zones[self.cell_zone[c]].volume += v;
        }
        let patches = self
            .patches
            .iter()
            .enumerate()
            .map(|(p, patch)| PatchSummary {
                name: patch.name.clone(),
                role: patch.role,
                foehncraft_type: patch.role.foehncraft_type(),
                faces: self.patch_faces(p).len(),
                area: self
                    .patch_faces(p)
                    .map(|f| ring_area(&self.points, self.face(f)).length())
                    .sum(),
            })
            .collect();
        HandoffSummary {
            cells: self.cell_count(),
            faces: self.face_count(),
            interior_faces: self.neighbour.len(),
            boundary_faces: self.face_count() - self.neighbour.len(),
            points: self.points.len(),
            total_volume: measures.iter().map(|m| m.0).sum(),
            min_volume: measures.iter().map(|m| m.0).fold(f64::INFINITY, f64::min),
            max_relative_closure: measures.iter().map(|m| m.1 / m.2).fold(0.0, f64::max),
            merged_faces: self.merged_faces,
            inverted_faces: self.inverted_faces(),
            patches,
            zones,
        }
    }

    /// Patches that are empty; FoehnCraft drops them, so they are reported.
    pub fn empty_patches(&self) -> Vec<&str> {
        (0..self.patches.len())
            .filter(|&p| self.patch_faces(p).is_empty())
            .map(|p| self.patches[p].name.as_str())
            .collect()
    }

    /// Number of faces per cell, for statistics.
    pub fn faces_per_cell(&self) -> BTreeMap<usize, usize> {
        let cells = self.cell_face_list();
        let mut counts = BTreeMap::new();
        for c in 0..self.cell_count() {
            *counts.entry(cells.len(c)).or_insert(0) += 1;
        }
        counts
    }
}

// ---------------------------------------------------------------------
// writers

/// Text gathered in a buffer and handed on in large writes: the writers
/// below print hundreds of millions of numbers for a large mesh, and one
/// formatted write per number costs more than the numbers.
struct Text<'a> {
    out: &'a mut dyn Write,
    buf: Vec<u8>,
}

impl<'a> Text<'a> {
    const SPILL: usize = 1 << 16;

    fn new(out: &'a mut dyn Write) -> Self {
        Self {
            out,
            buf: Vec::with_capacity(Self::SPILL + 4096),
        }
    }

    fn str(&mut self, s: &str) {
        self.buf.extend_from_slice(s.as_bytes());
    }

    /// A non-negative integer in decimal, as `{}` prints it.
    fn uint(&mut self, mut v: usize) {
        let mut digits = [0u8; 20];
        let mut i = digits.len();
        loop {
            i -= 1;
            digits[i] = b'0' + (v % 10) as u8;
            v /= 10;
            if v == 0 {
                break;
            }
        }
        self.buf.extend_from_slice(&digits[i..]);
    }

    fn byte(&mut self, b: u8) {
        self.buf.push(b);
    }

    /// Hand the buffer on once it is large; call between lines.
    fn spill(&mut self) -> io::Result<()> {
        if self.buf.len() >= Self::SPILL {
            self.out.write_all(&self.buf)?;
            self.buf.clear();
        }
        Ok(())
    }

    fn finish(self) -> io::Result<()> {
        self.out.write_all(&self.buf)
    }
}

fn string_array(t: &mut Text<'_>, name: &str, values: &[String]) -> io::Result<()> {
    t.str("      <DataArray type=\"String\" Name=\"");
    t.str(name);
    t.str("\" NumberOfTuples=\"");
    t.uint(values.len());
    t.str("\" format=\"ascii\">");
    for v in values {
        for b in v.bytes() {
            t.uint(usize::from(b));
            t.byte(b' ');
        }
        t.str("0 ");
        t.spill()?;
    }
    t.str("</DataArray>\n");
    Ok(())
}

/// The sorted distinct nodes of cell `c` (its VTK connectivity) into
/// `nodes`.
fn cell_nodes(mesh: &HandoffMesh, cells: &CellFaceList, c: usize, nodes: &mut Vec<usize>) {
    nodes.clear();
    for (f, _) in cells.cell(c) {
        nodes.extend_from_slice(mesh.face(f));
    }
    nodes.sort_unstable();
    nodes.dedup();
}

/// Write the single-file VTU FoehnCraft reads (FC-02, its D192): every cell a
/// polyhedron (type 42), every boundary face a polygon (type 7), the cell
/// array `region` (zone index for cells, patch index for faces) and the
/// field-data string arrays `zone_names` and `boundary_names` (`id=name`).
/// `metadata` adds string arrays such as the mesh digest (FC-05).
///
/// The file is produced in one pass into `out` (no copy of it in memory);
/// besides the mesh it needs the cell-face list and the cell offsets.
pub fn write_foehncraft_vtu(
    mesh: &HandoffMesh,
    metadata: &[(&str, String)],
    out: &mut dyn Write,
) -> io::Result<()> {
    let cells = mesh.cell_face_list();
    let n_cells = mesh.cell_count();
    let boundary = mesh.neighbour.len()..mesh.face_count();
    let total = n_cells + boundary.len();
    let mut t = Text::new(out);
    t.str("<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n  <UnstructuredGrid>\n    <FieldData>\n");
    let ids = |names: Vec<&str>| -> Vec<String> {
        names
            .iter()
            .enumerate()
            .map(|(i, n)| format!("{i}={n}"))
            .collect()
    };
    string_array(
        &mut t,
        "zone_names",
        &ids(mesh.zones.iter().map(String::as_str).collect()),
    )?;
    string_array(
        &mut t,
        "boundary_names",
        &ids(mesh.patches.iter().map(|p| p.name.as_str()).collect()),
    )?;
    string_array(
        &mut t,
        "boundary_roles",
        &ids(mesh
            .patches
            .iter()
            .map(|p| match p.role {
                PatchRole::Wall => "wall",
                PatchRole::Inlet => "inlet",
                PatchRole::Outlet => "outlet",
                PatchRole::Symmetry => "symmetry",
                PatchRole::Interface => "interface",
                PatchRole::Unspecified => "unspecified",
            })
            .collect()),
    )?;
    for (name, value) in metadata {
        string_array(&mut t, name, std::slice::from_ref(value))?;
    }
    t.str("    </FieldData>\n    <Piece NumberOfPoints=\"");
    t.uint(mesh.points.len());
    t.str("\" NumberOfCells=\"");
    t.uint(total);
    t.str("\">\n      <Points>\n        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for p in &mesh.points {
        writeln!(t.buf, "{:e} {:e} {:e}", p.x, p.y, p.z)?;
        t.spill()?;
    }
    t.str("        </DataArray>\n      </Points>\n      <Cells>\n");

    t.str("        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    let mut offsets = Vec::with_capacity(total);
    let mut nodes = Vec::new();
    let mut n = 0usize;
    for c in 0..n_cells {
        cell_nodes(mesh, &cells, c, &mut nodes);
        for &v in &nodes {
            t.uint(v);
            t.byte(b' ');
        }
        t.byte(b'\n');
        t.spill()?;
        n += nodes.len();
        offsets.push(n);
    }
    for f in boundary.clone() {
        for &v in mesh.face(f) {
            t.uint(v);
            t.byte(b' ');
        }
        t.byte(b'\n');
        t.spill()?;
        n += mesh.face(f).len();
        offsets.push(n);
    }
    t.str("        </DataArray>\n        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    for o in offsets {
        t.uint(o);
        t.byte(b'\n');
        t.spill()?;
    }
    t.str("        </DataArray>\n        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for _ in 0..n_cells {
        t.str("42\n");
        t.spill()?;
    }
    for _ in boundary.clone() {
        t.str("7\n");
        t.spill()?;
    }
    t.str("        </DataArray>\n        <DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">\n");
    let mut face_offsets = Vec::with_capacity(n_cells);
    let mut n = 0usize;
    for c in 0..n_cells {
        t.uint(cells.len(c));
        n += 1;
        for (f, outward) in cells.cell(c) {
            let ring = mesh.face(f);
            t.byte(b' ');
            t.uint(ring.len());
            let mut put = |v: usize| {
                t.byte(b' ');
                t.uint(v);
            };
            if outward {
                ring.iter().for_each(|&v| put(v));
            } else {
                ring.iter().rev().for_each(|&v| put(v));
            }
            n += 1 + ring.len();
        }
        t.byte(b'\n');
        t.spill()?;
        face_offsets.push(n);
    }
    t.str("        </DataArray>\n        <DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">\n");
    for o in face_offsets {
        t.uint(o);
        t.byte(b'\n');
        t.spill()?;
    }
    for _ in boundary.clone() {
        t.str("-1\n");
        t.spill()?;
    }
    t.str("        </DataArray>\n      </Cells>\n      <CellData Scalars=\"region\">\n        <DataArray type=\"Int32\" Name=\"region\" format=\"ascii\">\n");
    for &z in &mesh.cell_zone {
        t.uint(z);
        t.byte(b'\n');
        t.spill()?;
    }
    for p in 0..mesh.patches.len() {
        for _ in mesh.patch_faces(p) {
            t.uint(p);
            t.byte(b'\n');
            t.spill()?;
        }
    }
    t.str(
        "        </DataArray>\n      </CellData>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n",
    );
    t.finish()
}

fn foam_header(out: &mut dyn Write, class: &str, object: &str, note: &str) -> io::Result<()> {
    writeln!(
        out,
        "FoamFile\n{{\n    version     2.0;\n    format      ascii;\n    class       {class};\n    note        \"{note}\";\n    location    \"constant/polyMesh\";\n    object      {object};\n}}\n"
    )
}

/// The files of an OpenFOAM `constant/polyMesh` directory (FC-03), in the
/// order they are written.
pub const POLYMESH_FILES: [&str; 6] = [
    "points",
    "faces",
    "owner",
    "neighbour",
    "boundary",
    "cellZones",
];

/// Write one file of [`POLYMESH_FILES`] (FC-03): `points`, `faces`, `owner`,
/// `neighbour`, `boundary` or `cellZones`. Patch types follow the roles
/// (`wall`, `symmetry`, otherwise `patch`); empty patches are left out, as
/// OpenFOAM tools do.
pub fn write_polymesh_file(mesh: &HandoffMesh, file: &str, out: &mut dyn Write) -> io::Result<()> {
    let note = format!(
        "nPoints:{} nCells:{} nFaces:{} nInternalFaces:{}",
        mesh.points.len(),
        mesh.cell_count(),
        mesh.face_count(),
        mesh.neighbour.len()
    );
    let mut t = Text::new(out);
    let list = |t: &mut Text<'_>, values: &[usize]| -> io::Result<()> {
        t.uint(values.len());
        t.str("\n(\n");
        for &x in values {
            t.uint(x);
            t.byte(b'\n');
            t.spill()?;
        }
        t.str(")\n");
        Ok(())
    };
    match file {
        "points" => {
            foam_header(&mut t.buf, "vectorField", "points", &note)?;
            t.uint(mesh.points.len());
            t.str("\n(\n");
            for p in &mesh.points {
                writeln!(t.buf, "({:e} {:e} {:e})", p.x, p.y, p.z)?;
                t.spill()?;
            }
            t.str(")\n");
        }
        "faces" => {
            foam_header(&mut t.buf, "faceList", "faces", &note)?;
            t.uint(mesh.face_count());
            t.str("\n(\n");
            for f in 0..mesh.face_count() {
                let r = mesh.face(f);
                t.uint(r.len());
                t.byte(b'(');
                for (i, &v) in r.iter().enumerate() {
                    if i > 0 {
                        t.byte(b' ');
                    }
                    t.uint(v);
                }
                t.str(")\n");
                t.spill()?;
            }
            t.str(")\n");
        }
        "owner" => {
            foam_header(&mut t.buf, "labelList", "owner", &note)?;
            list(&mut t, &mesh.owner)?;
        }
        "neighbour" => {
            foam_header(&mut t.buf, "labelList", "neighbour", &note)?;
            list(&mut t, &mesh.neighbour)?;
        }
        "boundary" => {
            foam_header(&mut t.buf, "polyBoundaryMesh", "boundary", &note)?;
            let used: Vec<usize> = (0..mesh.patches.len())
                .filter(|&p| !mesh.patch_faces(p).is_empty())
                .collect();
            writeln!(t.buf, "{}\n(", used.len())?;
            for p in used {
                let patch = &mesh.patches[p];
                let range = mesh.patch_faces(p);
                let group = if patch.role == PatchRole::Wall {
                    "\n        inGroups    List<word> 1(wall);"
                } else {
                    ""
                };
                writeln!(
                    t.buf,
                    "    {}\n    {{\n        type        {};{group}\n        nFaces      {};\n        startFace   {};\n    }}",
                    patch.name,
                    patch.role.openfoam_type(),
                    range.len(),
                    range.start
                )?;
            }
            t.str(")\n");
        }
        "cellZones" => {
            foam_header(&mut t.buf, "regIOobject", "cellZones", &note)?;
            let mut members: Vec<Vec<usize>> = vec![Vec::new(); mesh.zones.len()];
            for (c, &z) in mesh.cell_zone.iter().enumerate() {
                members[z].push(c);
            }
            writeln!(t.buf, "{}\n(", mesh.zones.len())?;
            for (name, members) in mesh.zones.iter().zip(&members) {
                writeln!(
                    t.buf,
                    "{name}\n{{\n    type cellZone;\ncellLabels      List<label> {}\n(",
                    members.len()
                )?;
                for &c in members {
                    t.uint(c);
                    t.byte(b'\n');
                    t.spill()?;
                }
                t.str(")\n;\n}\n");
            }
            t.str(")\n");
        }
        other => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("{other} is not a polyMesh file"),
            ));
        }
    }
    t.finish()
}

/// The files of [`POLYMESH_FILES`] in memory, for tests and small meshes.
pub fn write_polymesh(mesh: &HandoffMesh) -> io::Result<Vec<(&'static str, Vec<u8>)>> {
    POLYMESH_FILES
        .iter()
        .map(|&name| {
            let mut bytes = Vec::new();
            write_polymesh_file(mesh, name, &mut bytes)?;
            Ok((name, bytes))
        })
        .collect()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// A unit cube split into two hexahedra along x, with boundary patches
    /// `xmin`, `xmax` and `sides`.
    pub(crate) fn two_hexes() -> HandoffMesh {
        let mut points = Vec::new();
        for z in [0.0, 1.0] {
            for y in [0.0, 1.0] {
                for x in [0.0, 0.5, 1.0] {
                    points.push(Vec3::new(x, y, z));
                }
            }
        }
        let id = |i: usize, j: usize, k: usize| i + 3 * j + 6 * k;
        let hex = |i: usize| -> Vec<Vec<usize>> {
            vec![
                vec![id(i, 0, 0), id(i, 1, 0), id(i + 1, 1, 0), id(i + 1, 0, 0)], // z-
                vec![id(i, 0, 1), id(i + 1, 0, 1), id(i + 1, 1, 1), id(i, 1, 1)], // z+
                vec![id(i, 0, 0), id(i + 1, 0, 0), id(i + 1, 0, 1), id(i, 0, 1)], // y-
                vec![id(i, 1, 0), id(i, 1, 1), id(i + 1, 1, 1), id(i + 1, 1, 0)], // y+
                vec![id(i, 0, 0), id(i, 0, 1), id(i, 1, 1), id(i, 1, 0)],         // x-
                vec![
                    id(i + 1, 0, 0),
                    id(i + 1, 1, 0),
                    id(i + 1, 1, 1),
                    id(i + 1, 0, 1),
                ], // x+
            ]
        };
        let cells: Vec<CellFaces> = (0..2)
            .map(|i| CellFaces {
                faces: hex(i),
                zone: 0,
            })
            .collect();
        let mut boundary = Vec::new();
        for (c, cell) in cells.iter().enumerate() {
            for (f, ring) in cell.faces.iter().enumerate() {
                let patch = match (c, f) {
                    (0, 4) => 0,
                    (1, 5) => 1,
                    (_, 0..=3) => 2,
                    _ => continue,
                };
                boundary.push((ring.clone(), patch));
            }
        }
        HandoffMesh::from_cells(
            points,
            cells,
            vec!["fluid".into()],
            &boundary,
            patches_from_names(&["xmin", "xmax", "side walls"]),
        )
        .unwrap()
    }

    /// Invariant: a conforming cell set becomes one interior face per
    /// shared pair, boundary faces grouped by patch, exact volume and
    /// areas. Breaks if faces are matched by winding instead of vertex set
    /// or patches lose faces; the format writers rely on this order.
    #[test]
    fn two_hexes_share_one_face_and_keep_their_patches() {
        let m = two_hexes();
        let s = m.summary();
        assert_eq!((s.cells, s.faces, s.interior_faces), (2, 11, 1));
        assert!((s.total_volume - 1.0).abs() < 1e-15);
        let areas: Vec<_> = s.patches.iter().map(|p| (p.faces, p.area)).collect();
        assert_eq!(areas, vec![(1, 1.0), (1, 1.0), (8, 4.0)]);
        assert_eq!(m.patches[2].name, "side_walls");
        assert_eq!(m.patches[2].role, PatchRole::Wall);
        assert_eq!(m.patches[2].source_name.as_deref(), Some("side walls"));
        assert!(s.max_relative_closure < 1e-15);
    }

    /// Invariant: an unclassified face is an error, never a silent extra
    /// patch (AGENTS.md: no partial result as success).
    #[test]
    fn a_face_without_neighbour_or_patch_is_refused() {
        let m = two_hexes();
        let cells: Vec<CellFaces> = m
            .cell_faces()
            .iter()
            .map(|faces| CellFaces {
                faces: faces
                    .iter()
                    .map(|&(f, out)| {
                        let mut r = m.face(f).to_vec();
                        if !out {
                            r.reverse();
                        }
                        r
                    })
                    .collect(),
                zone: 0,
            })
            .collect();
        let boundary: Vec<_> = (1..m.face_count())
            .filter(|&f| f != 1)
            .map(|f| (m.face(f).to_vec(), 0))
            .collect();
        let error = HandoffMesh::from_cells(
            m.points.clone(),
            cells,
            vec!["fluid".into()],
            &boundary,
            patches_from_names(&["all"]),
        )
        .unwrap_err();
        assert!(error.0.contains("no boundary face"), "{error}");
    }

    #[test]
    fn handoff_names_are_bare_keys_and_unique() {
        let p = patches_from_names(&["west inlet", "west/inlet", "walls"]);
        let names: Vec<_> = p.iter().map(|p| p.name.as_str()).collect();
        assert_eq!(names, ["west_inlet", "west_inlet_2", "walls"]);
        assert_eq!(p[0].role, PatchRole::Inlet);
    }
}

#[cfg(test)]
mod writer_tests {
    use super::tests::two_hexes;
    use super::*;

    /// Invariant: two cells that share two faces get one merged face whose
    /// area vector is the sum of the pieces (FoehnCraft refuses two faces
    /// between one pair). Breaks if the merge drops or duplicates an edge;
    /// only this test builds such a pair.
    #[test]
    fn two_faces_between_one_pair_merge_into_one() {
        let m = two_hexes();
        // Split the shared face x = 0.5 into two triangles.
        let shared = m.face(0).to_vec();
        let mut cells: Vec<CellFaces> = m
            .cell_faces()
            .iter()
            .map(|faces| CellFaces {
                faces: faces
                    .iter()
                    .flat_map(|&(f, out)| {
                        let r = m.face(f).to_vec();
                        let mut pieces = if f == 0 {
                            vec![vec![r[0], r[1], r[2]], vec![r[0], r[2], r[3]]]
                        } else {
                            vec![r]
                        };
                        if !out {
                            pieces.iter_mut().for_each(|p| p.reverse());
                        }
                        pieces
                    })
                    .collect(),
                zone: 0,
            })
            .collect();
        cells.truncate(2);
        let boundary: Vec<_> = (1..m.face_count())
            .map(|f| {
                let p = (0..3).find(|&p| m.patch_faces(p).contains(&f)).unwrap();
                (m.face(f).to_vec(), p)
            })
            .collect();
        let merged = HandoffMesh::from_cells(
            m.points.clone(),
            cells,
            vec!["fluid".into()],
            &boundary,
            m.patches.clone(),
        )
        .unwrap();
        assert_eq!(merged.merged_faces, 1);
        assert_eq!(merged.neighbour.len(), 1);
        let mut a = merged.face(0).to_vec();
        let mut b = shared;
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
        assert_eq!(merged.summary().inverted_faces, 0);
    }

    /// The VTU writer as it was before it streamed: one formatted write
    /// per number, the cells' faces in vectors.
    fn reference_string_array(
        out: &mut dyn Write,
        name: &str,
        values: &[String],
    ) -> io::Result<()> {
        write!(
            out,
            "      <DataArray type=\"String\" Name=\"{name}\" NumberOfTuples=\"{}\" format=\"ascii\">",
            values.len()
        )?;
        for v in values {
            for b in v.bytes() {
                write!(out, "{b} ")?;
            }
            write!(out, "0 ")?;
        }
        writeln!(out, "</DataArray>")
    }

    fn reference_vtu(
        mesh: &HandoffMesh,
        metadata: &[(&str, String)],
        out: &mut dyn Write,
    ) -> io::Result<()> {
        // the cells' faces as the writer built them before the cell-face list
        let mut cells = vec![Vec::new(); mesh.cell_count()];
        for (f, &o) in mesh.owner.iter().enumerate() {
            cells[o].push((f, true));
        }
        for (f, &n) in mesh.neighbour.iter().enumerate() {
            cells[n].push((f, false));
        }
        let boundary = mesh.neighbour.len()..mesh.face_count();
        let total = cells.len() + boundary.len();
        writeln!(
            out,
            "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n  <UnstructuredGrid>\n    <FieldData>"
        )?;
        let ids = |names: Vec<&str>| -> Vec<String> {
            names
                .iter()
                .enumerate()
                .map(|(i, n)| format!("{i}={n}"))
                .collect()
        };
        reference_string_array(
            out,
            "zone_names",
            &ids(mesh.zones.iter().map(String::as_str).collect()),
        )?;
        reference_string_array(
            out,
            "boundary_names",
            &ids(mesh.patches.iter().map(|p| p.name.as_str()).collect()),
        )?;
        reference_string_array(
            out,
            "boundary_roles",
            &ids(mesh
                .patches
                .iter()
                .map(|p| match p.role {
                    PatchRole::Wall => "wall",
                    PatchRole::Inlet => "inlet",
                    PatchRole::Outlet => "outlet",
                    PatchRole::Symmetry => "symmetry",
                    PatchRole::Interface => "interface",
                    PatchRole::Unspecified => "unspecified",
                })
                .collect()),
        )?;
        for (name, value) in metadata {
            reference_string_array(out, name, std::slice::from_ref(value))?;
        }
        writeln!(
            out,
            "    </FieldData>\n    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{total}\">\n      <Points>\n        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">",
            mesh.points.len()
        )?;
        for p in &mesh.points {
            writeln!(out, "{:e} {:e} {:e}", p.x, p.y, p.z)?;
        }
        writeln!(out, "        </DataArray>\n      </Points>\n      <Cells>")?;

        let ring = |f: usize, outward: bool| -> Vec<usize> {
            let mut r = mesh.face(f).to_vec();
            if !outward {
                r.reverse();
            }
            r
        };
        writeln!(
            out,
            "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">"
        )?;
        let mut offsets = Vec::with_capacity(total);
        let mut n = 0usize;
        for faces in &cells {
            let mut nodes: Vec<usize> = faces
                .iter()
                .flat_map(|&(f, _)| mesh.face(f).iter().copied())
                .collect();
            nodes.sort_unstable();
            nodes.dedup();
            for v in &nodes {
                write!(out, "{v} ")?;
            }
            writeln!(out)?;
            n += nodes.len();
            offsets.push(n);
        }
        for f in boundary.clone() {
            for v in mesh.face(f) {
                write!(out, "{v} ")?;
            }
            writeln!(out)?;
            n += mesh.face(f).len();
            offsets.push(n);
        }
        writeln!(
            out,
            "        </DataArray>\n        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"
        )?;
        for o in offsets {
            writeln!(out, "{o}")?;
        }
        writeln!(
            out,
            "        </DataArray>\n        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
        )?;
        for _ in 0..cells.len() {
            writeln!(out, "42")?;
        }
        for _ in boundary.clone() {
            writeln!(out, "7")?;
        }
        writeln!(
            out,
            "        </DataArray>\n        <DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">"
        )?;
        let mut face_offsets = Vec::with_capacity(total);
        let mut n = 0i64;
        for faces in &cells {
            write!(out, "{}", faces.len())?;
            n += 1;
            for &(f, outward) in faces {
                let r = ring(f, outward);
                write!(out, " {}", r.len())?;
                for v in &r {
                    write!(out, " {v}")?;
                }
                n += 1 + r.len() as i64;
            }
            writeln!(out)?;
            face_offsets.push(n);
        }
        writeln!(
            out,
            "        </DataArray>\n        <DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">"
        )?;
        for o in face_offsets {
            writeln!(out, "{o}")?;
        }
        for _ in boundary.clone() {
            writeln!(out, "-1")?;
        }
        writeln!(
            out,
            "        </DataArray>\n      </Cells>\n      <CellData Scalars=\"region\">\n        <DataArray type=\"Int32\" Name=\"region\" format=\"ascii\">"
        )?;
        for &z in &mesh.cell_zone {
            writeln!(out, "{z}")?;
        }
        for p in 0..mesh.patches.len() {
            for _ in mesh.patch_faces(p) {
                writeln!(out, "{p}")?;
            }
        }
        writeln!(
            out,
            "        </DataArray>\n      </CellData>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>"
        )
    }

    /// Invariant: the streamed writer prints the same bytes as the plain
    /// one it replaced, and its integer printing equals `{}`. Breaks if the
    /// buffered digits, the cell-face order or a section header drift;
    /// FoehnCraft and MurmelCraft read these files and their digests are
    /// recorded, and the layout test below checks only a few markers.
    #[test]
    fn streamed_vtu_equals_the_plain_writer() {
        let m = two_hexes();
        let meta = [
            ("koflcraft_mesh_digest", "abc".to_string()),
            ("koflcraft_algorithm", "x y".to_string()),
        ];
        for metadata in [&meta[..], &[]] {
            let (mut new, mut old) = (Vec::new(), Vec::new());
            write_foehncraft_vtu(&m, metadata, &mut new).unwrap();
            reference_vtu(&m, metadata, &mut old).unwrap();
            assert_eq!(
                String::from_utf8(new).unwrap(),
                String::from_utf8(old).unwrap()
            );
        }
        let mut out = Vec::new();
        let mut t = Text::new(&mut out);
        let values = [0, 7, 9, 10, 99, 100, 12_345, u32::MAX as usize, usize::MAX];
        for v in values {
            t.uint(v);
            t.byte(b' ');
        }
        t.finish().unwrap();
        let expected: String = values.iter().map(|v| format!("{v} ")).collect();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    /// Invariant: the VTU holds every cell as a polyhedron and every
    /// boundary face as a polygon with its patch in `region` and the names
    /// in the field data, the layout FoehnCraft's reader expects (its
    /// D192). FC-06 reads the files with FoehnCraft itself; this test keeps
    /// the layout checkable without it.
    #[test]
    fn foehncraft_vtu_holds_cells_faces_and_names() {
        let m = two_hexes();
        let mut out = Vec::new();
        write_foehncraft_vtu(&m, &[("koflcraft_mesh_digest", "abc".into())], &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("NumberOfCells=\"12\""));
        let types = text.split("Name=\"types\"").nth(1).unwrap();
        let types = types.split("</DataArray>").next().unwrap();
        assert_eq!(types.matches("42").count(), 2);
        assert_eq!(types.lines().filter(|l| l.trim() == "7").count(), 10);
        let codes = |s: &str| {
            s.bytes()
                .map(|b| b.to_string())
                .collect::<Vec<_>>()
                .join(" ")
        };
        assert!(text.contains(&codes("2=side_walls")));
        assert!(text.contains(&codes("0=fluid")));
        assert!(text.contains(&codes("abc")));
    }

    /// Invariant: OpenFOAM `boundary` start faces and counts tile the
    /// boundary faces after the interior ones, in patch order.
    #[test]
    fn polymesh_boundary_tiles_the_faces() {
        let m = two_hexes();
        let files = write_polymesh(&m).unwrap();
        let boundary = String::from_utf8(files[4].1.clone()).unwrap();
        for (name, n, start) in [("xmin", 1, 1), ("xmax", 1, 2), ("side_walls", 8, 3)] {
            let block = boundary.split(&format!("    {name}\n")).nth(1).unwrap();
            assert!(block.contains(&format!("nFaces      {n};")));
            assert!(block.contains(&format!("startFace   {start};")));
        }
        assert!(boundary.contains("type        wall;"));
    }
}
