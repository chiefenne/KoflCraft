use super::*;
use std::collections::{BTreeMap, BTreeSet, HashMap};

pub struct Measurements {
    pub volume: f64,
    /// Minimum, fifth percentile and median over core polyhedra.
    pub core_sphericity: Option<[f64; 3]>,
    pub minimum_volume: f64,
    pub minimum_wedge_jacobian: Option<f64>,
    /// Boundary area of every patch id, indexed by id - 1.
    pub areas: Vec<f64>,
    pub prism_core_faces: usize,
    pub port_layer_areas: [Vec<f64>; 2],
    pub minimum_tet_mean_ratio: f64,
    pub tetrahedra_below_mean_ratio_0_1: usize,
    pub normal_thickness_ranges: Vec<[f64; 2]>,
    pub minimum_kernel_clearance: f64,
    /// The solver rules, counted beside the strict checks.
    pub solver_rules: SolverRules,
}

/// Cells counted against the validity rules of production finite-volume
/// codes (OpenFOAM's mesh description and checkMesh, the STAR-CCM+ and
/// Fluent cell checks; research/robust-polyhedral-meshing): a cell is
/// usable when its shell is closed, its volume is positive, every face
/// pyramid on the cell centre is positive and every face points away from
/// the cell (owner-to-neighbour vector, or centre-to-face vector at the
/// boundary, at less than 90° to the face normal). Faces need not be flat:
/// a face's area vector, centre and pyramids come from the triangle fan on
/// its vertex mean, as in those codes. Report only; the strict checks
/// below still decide.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SolverRules {
    pub cells: usize,
    /// Cells failing at least one rule.
    pub failing_cells: usize,
    /// Failing cells of the wall layers (layer > 0); the rest are core cells.
    pub failing_layer_cells: usize,
    /// Shell not closed: an edge not used by exactly two faces of the cell,
    /// or the face area vectors not summing to zero (1e-6 of their sum).
    pub open_cells: usize,
    pub nonpositive_volume: usize,
    /// Cells with a face pyramid of nonpositive volume.
    pub negative_pyramids: usize,
    /// Cells with a face at 90° or more (non-orthogonality).
    pub faces_turned_away: usize,
    /// Largest non-orthogonality over all faces, degrees.
    pub maximum_non_orthogonality: f64,
    /// Largest face warping, degrees: the angle between a triangle of a
    /// face's fan and the face's area vector (a graded quality measure, as
    /// the warpage angle of STAR-CCM+; no limit).
    pub maximum_warping: f64,
}

/// Area vector and centre of a face as finite-volume codes take them:
/// triangles fanned from the vertex mean, centres weighted by area.
fn fv_face(p: &[Vec3], ring: &[usize]) -> (Vec3, Vec3) {
    let mean = ring.iter().fold(Vec3::ZERO, |s, &v| s + p[v]) / ring.len() as f64;
    if ring.len() == 3 {
        let [a, b, c] = [p[ring[0]], p[ring[1]], p[ring[2]]];
        return ((b - a).cross(c - a) * 0.5, mean);
    }
    let mut area = Vec3::ZERO;
    let mut moment = Vec3::ZERO;
    let mut weight = 0.0;
    for (i, &v) in ring.iter().enumerate() {
        let (a, b) = (p[v], p[ring[(i + 1) % ring.len()]]);
        let n = (a - mean).cross(b - mean) * 0.5;
        area = area + n;
        let w = n.length();
        moment = moment + (a + b + mean) * (w / 3.0);
        weight += w;
    }
    let centre = if weight > 0.0 { moment / weight } else { mean };
    (area, centre)
}

/// Counts the cells of a mesh against [`SolverRules`].
#[cfg(test)]
fn solver_rules(p: &[Vec3], cells: &[Cell]) -> SolverRules {
    match solver_rules_and_volumes(p, cells) {
        Ok((rules, _, _)) => rules,
        Err(e) => panic!("{e}"),
    }
}

/// Largest angle between a fan triangle of a face and its area vector.
fn warping(p: &[Vec3], ring: &[usize], area: Vec3) -> f64 {
    if ring.len() == 3 {
        return 0.0;
    }
    let mean = ring.iter().fold(Vec3::ZERO, |s, &v| s + p[v]) / ring.len() as f64;
    let mut worst = 0.0_f64;
    for (i, &v) in ring.iter().enumerate() {
        let n = (p[v] - mean).cross(p[ring[(i + 1) % ring.len()]] - mean);
        let cos = n.dot(area) / (n.length() * area.length());
        if cos.is_finite() {
            worst = worst.max(cos.clamp(-1.0, 1.0).acos().to_degrees());
        }
    }
    worst
}

/// A cell's shell as the solver sees it.
pub(super) struct Shell {
    pub open: bool,
    pub volume: f64,
    pub centre: Vec3,
    pub pyramids_ok: bool,
}

impl Shell {
    /// Closed, positive volume and every face pyramid on the centre positive
    /// (the solver rules that need no neighbour).
    pub fn valid(&self) -> bool {
        !self.open && self.volume > 0.0 && self.pyramids_ok
    }
}

/// Face rings stored one after another, so the faces of a cell need no
/// allocation per face.
#[derive(Default)]
pub(super) struct Rings {
    flat: Vec<usize>,
    ends: Vec<usize>,
}

impl Rings {
    fn clear(&mut self) {
        self.flat.clear();
        self.ends.clear();
    }
    fn push(&mut self, ring: &[usize]) {
        self.flat.extend_from_slice(ring);
        self.ends.push(self.flat.len());
    }
    fn len(&self) -> usize {
        self.ends.len()
    }
    fn get(&self, i: usize) -> &[usize] {
        let start = if i == 0 { 0 } else { self.ends[i - 1] };
        &self.flat[start..self.ends[i]]
    }
    fn iter(&self) -> impl Iterator<Item = &[usize]> + '_ {
        (0..self.len()).map(|i| self.get(i))
    }
    /// The faces of an element, in the order of [`Element::faces`].
    fn load(&mut self, element: &Element) {
        self.clear();
        match *element {
            Element::Polyhedron(ref poly) => poly.faces.iter().for_each(|r| self.push(r)),
            Element::Tet([a, b, c, d]) => {
                for ring in [[a, c, b], [a, b, d], [b, c, d], [c, a, d]] {
                    self.push(&ring);
                }
            }
            Element::Wedge([a, b, c, d, e, f]) => {
                self.push(&[a, c, b]);
                self.push(&[d, e, f]);
                self.push(&[a, b, e, d]);
                self.push(&[b, c, f, e]);
                self.push(&[c, a, d, f]);
            }
        }
    }
}

/// The shell of a cell given by its outward face rings.
pub(super) fn cell_shell(p: &[Vec3], rings: &[Vec<usize>]) -> Shell {
    let mut flat = Rings::default();
    rings.iter().for_each(|r| flat.push(r));
    let geometry: Vec<(Vec3, Vec3)> = flat.iter().map(|r| fv_face(p, r)).collect();
    shell(p, &flat, &geometry, &mut Vec::new())
}

fn shell(
    p: &[Vec3],
    rings: &Rings,
    geometry: &[(Vec3, Vec3)],
    edges: &mut Vec<(usize, usize)>,
) -> Shell {
    // closed shell: every edge twice, area vectors summing to zero
    edges.clear();
    for r in rings.iter() {
        for (i, &a) in r.iter().enumerate() {
            let b = r[(i + 1) % r.len()];
            edges.push((a.min(b), a.max(b)));
        }
    }
    edges.sort_unstable();
    let twice = edges.chunk_by(|x, y| x == y).all(|uses| uses.len() == 2);
    let sum = geometry.iter().fold(Vec3::ZERO, |s, g| s + g.0);
    let total: f64 = geometry.iter().map(|g| g.0.length()).sum();
    let open = !twice
        || sum.length().partial_cmp(&(1e-6 * total)) == Some(std::cmp::Ordering::Greater)
        || !total.is_finite();
    let (volume, centre) = fan_volume_centre(p, rings);
    let pyramids_ok = geometry.iter().all(|(area, c)| area.dot(*c - centre) > 0.0);
    Shell {
        open,
        volume,
        centre,
        pyramids_ok,
    }
}

/// The union of rings that share their inner edges, as one ring (the
/// hand-off merges the faces between two cells this way); `None` when the
/// union is not one loop.
pub(super) fn merge_rings<R: AsRef<[usize]>>(rings: &[R]) -> Option<Vec<usize>> {
    if let [ring] = rings {
        return Some(ring.as_ref().to_vec());
    }
    // directed edges left after the shared ones cancel, with their counts
    let mut edges: Vec<((usize, usize), u32)> = Vec::new();
    for ring in rings {
        let ring = ring.as_ref();
        for (k, &a) in ring.iter().enumerate() {
            let b = ring[(k + 1) % ring.len()];
            if let Some(i) = edges.iter().position(|&(e, _)| e == (b, a)) {
                edges[i].1 -= 1;
                if edges[i].1 == 0 {
                    edges.swap_remove(i);
                }
            } else if let Some(i) = edges.iter().position(|&(e, _)| e == (a, b)) {
                edges[i].1 += 1;
            } else {
                edges.push(((a, b), 1));
            }
        }
    }
    if edges.iter().any(|&(_, n)| n != 1) {
        return None;
    }
    let mut next: Vec<(usize, usize)> = edges.into_iter().map(|(e, _)| e).collect();
    next.sort_unstable();
    if next.windows(2).any(|w| w[0].0 == w[1].0) {
        return None;
    }
    let follow = |at: usize| {
        next.binary_search_by_key(&at, |e| e.0)
            .ok()
            .map(|i| next[i].1)
    };
    let &(start, _) = next.first()?;
    let mut ring = vec![start];
    let mut at = follow(start)?;
    while at != start && ring.len() <= next.len() {
        ring.push(at);
        at = follow(at)?;
    }
    (ring.len() == next.len()).then_some(ring)
}

/// Volume and centroid of a closed shell as the hand-off and FoehnCraft
/// take them: tetrahedra on the first vertex, each face fanned from its
/// vertex mean.
fn fan_volume_centre(p: &[Vec3], rings: &Rings) -> (f64, Vec3) {
    let Some(origin) = rings.iter().next().and_then(|r| r.first()).map(|&v| p[v]) else {
        return (0.0, Vec3::ZERO);
    };
    let mut v6 = 0.0;
    let mut moment = Vec3::ZERO;
    for ring in rings.iter() {
        let c = ring.iter().fold(Vec3::ZERO, |s, &v| s + p[v]) / ring.len() as f64 - origin;
        for (k, &a) in ring.iter().enumerate() {
            let (pa, pb) = (p[a] - origin, p[ring[(k + 1) % ring.len()]] - origin);
            let t6 = c.dot(pa.cross(pb));
            v6 += t6;
            moment = moment + (c + pa + pb) * (t6 / 4.0);
        }
    }
    (v6 / 6.0, origin + moment / v6)
}

/// A face's vertex set, sorted, with its hash: every cell that uses the face
/// makes the same key. Triangles and quadrilaterals need no allocation.
#[derive(Clone, PartialEq, Eq)]
pub struct FaceKey {
    hash: u64,
    vertices: KeyVertices,
}

#[derive(Clone, PartialEq, Eq)]
enum KeyVertices {
    Small([u32; 4]),
    Large(Box<[usize]>),
}

impl std::hash::Hash for FaceKey {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        state.write_u64(self.hash);
    }
}

fn mix(h: u64, v: u64) -> u64 {
    let x = (h ^ v).wrapping_mul(0x9e37_79b9_7f4a_7c15);
    x ^ (x >> 29)
}

/// The splitmix64 finalizer: every bit of the hash depends on every vertex.
fn avalanche(mut x: u64) -> u64 {
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

impl FaceKey {
    pub fn new(ring: &[usize]) -> Self {
        let small = ring.len() <= 4 && ring.iter().all(|&v| v < u32::MAX as usize);
        if small {
            let mut vertices = [u32::MAX; 4];
            for (slot, &v) in vertices.iter_mut().zip(ring) {
                *slot = v as u32;
            }
            vertices[..ring.len()].sort_unstable();
            let h = vertices
                .iter()
                .fold(ring.len() as u64, |h, &v| mix(h, u64::from(v)));
            Self {
                hash: avalanche(h),
                vertices: KeyVertices::Small(vertices),
            }
        } else {
            let mut vertices = ring.to_vec();
            vertices.sort_unstable();
            let h = vertices
                .iter()
                .fold(ring.len() as u64, |h, &v| mix(h, v as u64));
            Self {
                hash: avalanche(h),
                vertices: KeyVertices::Large(vertices.into_boxed_slice()),
            }
        }
    }
}

/// Hands the precomputed [`FaceKey`] hash to the map.
#[derive(Default)]
pub struct KeyHasher(u64);

impl std::hash::Hasher for KeyHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = mix(self.0, u64::from(b));
        }
    }
    fn write_u64(&mut self, v: u64) {
        self.0 = v;
    }
}

/// A map by face vertex set.
pub type KeyMap<V> = HashMap<FaceKey, V, std::hash::BuildHasherDefault<KeyHasher>>;
/// A set of face vertex sets.
pub type KeySet = std::collections::HashSet<FaceKey, std::hash::BuildHasherDefault<KeyHasher>>;

/// The cells that use one face, ascending (a cell that lists the face twice
/// appears twice): the order in which a pass over the cells meets them.
#[derive(Clone)]
enum Users {
    Few([usize; 2], u8),
    Many(Vec<usize>),
}

impl Users {
    fn as_slice(&self) -> &[usize] {
        match self {
            Users::Few(cells, n) => &cells[..usize::from(*n)],
            Users::Many(cells) => cells,
        }
    }
    fn insert(&mut self, cell: usize) {
        if let Users::Few(cells, n @ 1) = self {
            cells[1] = cell;
            *n = 2;
            if cells[0] > cells[1] {
                cells.swap(0, 1);
            }
            return;
        }
        let mut cells = self.as_slice().to_vec();
        let at = cells.partition_point(|&c| c <= cell);
        cells.insert(at, cell);
        *self = Users::Many(cells);
    }
    /// Removes one use by `cell`; true when no cell uses the face any more.
    fn remove(&mut self, cell: usize) -> bool {
        let mut cells = self.as_slice().to_vec();
        if let Some(i) = cells.iter().position(|&c| c == cell) {
            cells.remove(i);
        }
        *self = match cells[..] {
            [] => return true,
            [a] => Users::Few([a, 0], 1),
            [a, b] => Users::Few([a, b], 2),
            _ => Users::Many(cells),
        };
        false
    }
}

/// Cores used by the checks of the whole mesh.
fn threads() -> usize {
    std::thread::available_parallelism()
        .map_or(1, |n| n.get())
        .clamp(1, 64)
}

/// Runs `f` on the chunks of `0..count` on all cores (a shared queue, so
/// costly cells anywhere in the range balance out); results in chunk order.
fn chunked<T: Send>(
    count: usize,
    chunk: usize,
    f: impl Fn(std::ops::Range<usize>) -> T + Sync,
) -> Result<Vec<T>, String> {
    let chunk = chunk.max(1);
    let chunks = count.div_ceil(chunk);
    let next = std::sync::atomic::AtomicUsize::new(0);
    let workers = threads().min(chunks).max(1);
    let (f, next) = (&f, &next);
    let results = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                scope.spawn(move || {
                    let mut done = Vec::new();
                    loop {
                        let k = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        if k >= chunks {
                            break done;
                        }
                        let start = k * chunk;
                        done.push((k, f(start..(start + chunk).min(count))));
                    }
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join()).collect::<Vec<_>>()
    });
    let mut slots: Vec<Option<T>> = (0..chunks).map(|_| None).collect();
    for result in results {
        for (k, value) in result.map_err(|_| "a check worker panicked".to_string())? {
            slots[k] = Some(value);
        }
    }
    slots
        .into_iter()
        .map(|s| s.ok_or_else(|| "a check chunk is missing".to_string()))
        .collect()
}

/// Face `at` of an element (a polyhedron's own ring, no copy of the rest).
fn ring_of<'a>(element: &'a Element, at: usize, buffer: &'a mut Rings) -> &'a [usize] {
    match element {
        Element::Polyhedron(poly) => &poly.faces[at],
        _ => {
            buffer.load(element);
            buffer.get(at)
        }
    }
}

/// The number of faces of an element.
fn face_count(element: &Element) -> usize {
    match element {
        Element::Polyhedron(poly) => poly.faces.len(),
        Element::Tet(_) => 4,
        Element::Wedge(_) => 5,
    }
}

/// The faces of a mesh by vertex set, with the cells that use each one. The
/// map is split into shards by hash, built in parallel.
pub(super) struct FaceTable {
    shards: Vec<KeyMap<Users>>,
}

impl FaceTable {
    pub(super) fn new(cells: &[Cell]) -> Result<Self, String> {
        let count = threads();
        let rings: usize = cells.iter().map(|c| face_count(&c.element)).sum();
        let shards = chunked(count, 1, |shard| {
            let mut map = KeyMap::with_capacity_and_hasher(
                rings / 2 / count + 16,
                std::hash::BuildHasherDefault::default(),
            );
            let mut buffer = Rings::default();
            for (id, cell) in cells.iter().enumerate() {
                buffer.load(&cell.element);
                for ring in buffer.iter() {
                    let key = FaceKey::new(ring);
                    if Self::shard_of(key.hash, count) == shard.start {
                        map.entry(key)
                            .and_modify(|u: &mut Users| u.insert(id))
                            .or_insert(Users::Few([id, 0], 1));
                    }
                }
            }
            map
        })?;
        Ok(Self { shards })
    }

    fn shard_of(hash: u64, count: usize) -> usize {
        (hash >> 32) as usize % count
    }

    fn users(&self, ring: &[usize]) -> &[usize] {
        let key = FaceKey::new(ring);
        self.shards[Self::shard_of(key.hash, self.shards.len())]
            .get(&key)
            .map_or(&[], Users::as_slice)
    }

    /// The first cell other than `id` that uses the face, as the hand-off
    /// pairs it.
    fn other(&self, ring: &[usize], id: usize) -> Option<usize> {
        self.users(ring).iter().copied().find(|&c| c != id)
    }

    /// The cells that share a face with the cell, ascending, and itself not.
    fn neighbours(&self, id: usize, element: &Element, buffer: &mut Rings) -> Vec<usize> {
        buffer.load(element);
        let mut out: Vec<usize> = buffer
            .iter()
            .flat_map(|r| self.users(r).iter().copied())
            .filter(|&c| c != id)
            .collect();
        out.sort_unstable();
        out.dedup();
        out
    }

    /// Adds (`add`) or removes the uses of every face of a cell; returns the
    /// other cells that use those faces.
    fn update(&mut self, id: usize, element: &Element, add: bool) -> Vec<usize> {
        let mut buffer = Rings::default();
        buffer.load(element);
        let count = self.shards.len();
        let mut others = Vec::new();
        for ring in buffer.iter() {
            let key = FaceKey::new(ring);
            let shard = &mut self.shards[Self::shard_of(key.hash, count)];
            if add {
                let users = shard
                    .entry(key)
                    .and_modify(|u| u.insert(id))
                    .or_insert(Users::Few([id, 0], 1));
                others.extend(users.as_slice().iter().filter(|&&c| c != id));
            } else if let Some(users) = shard.get_mut(&key) {
                others.extend(users.as_slice().iter().filter(|&&c| c != id));
                if users.remove(id) {
                    shard.remove(&key);
                }
            }
        }
        others
    }
}

/// One cell as the solver sees it after the hand-off: the faces it shares
/// with one neighbour merged into one polygon. Buffers are reused from cell
/// to cell.
#[derive(Default)]
struct Exported {
    rings: Rings,
    tagged: Vec<(usize, usize)>,
    merged: Rings,
    /// Per merged face: the cell across, or none on the boundary.
    across: Vec<Option<usize>>,
    /// Per merged face: area vector and centre.
    geometry: Vec<(Vec3, Vec3)>,
    edges: Vec<(usize, usize)>,
}

/// The rules of one cell that need no neighbour, and its largest warping.
struct Verdict {
    open: bool,
    volume: f64,
    centre: Vec3,
    pyramids_ok: bool,
    warping: f64,
}

impl Verdict {
    fn nonpositive(&self) -> bool {
        self.volume.is_nan() || self.volume <= 0.0
    }
    fn bad(&self) -> bool {
        self.open || self.nonpositive() || !self.pyramids_ok
    }
}

impl Exported {
    fn evaluate(&mut self, p: &[Vec3], table: &FaceTable, id: usize, element: &Element) -> Verdict {
        self.rings.load(element);
        self.tagged.clear();
        self.merged.clear();
        self.across.clear();
        self.geometry.clear();
        for (i, ring) in self.rings.iter().enumerate() {
            match table.other(ring, id) {
                Some(other) => self.tagged.push((other, i)),
                None => {
                    self.merged.push(ring);
                    self.across.push(None);
                }
            }
        }
        // by neighbour, ascending; each neighbour's faces in cell order
        self.tagged.sort_unstable();
        let mut unmergeable = false;
        for group in self.tagged.chunk_by(|a, b| a.0 == b.0) {
            let other = group[0].0;
            if let [(_, i)] = group {
                self.merged.push(self.rings.get(*i));
                self.across.push(Some(other));
                continue;
            }
            let pieces: Vec<&[usize]> = group.iter().map(|&(_, i)| self.rings.get(i)).collect();
            match merge_rings(&pieces) {
                Some(ring) => {
                    self.merged.push(&ring);
                    self.across.push(Some(other));
                }
                None => {
                    unmergeable = true;
                    for ring in pieces {
                        self.merged.push(ring);
                        self.across.push(Some(other));
                    }
                }
            }
        }
        let mut worst = 0.0_f64;
        for ring in self.merged.iter() {
            let g = fv_face(p, ring);
            worst = worst.max(warping(p, ring, g.0));
            self.geometry.push(g);
        }
        let Shell {
            open,
            volume,
            centre,
            pyramids_ok,
        } = shell(p, &self.merged, &self.geometry, &mut self.edges);
        Verdict {
            open: open || unmergeable,
            volume,
            centre,
            pyramids_ok,
            warping: worst,
        }
    }
}

/// The angle in degrees between the owner-to-neighbour (or centre-to-face)
/// vector and the face normal, and whether the face points back (90° or
/// more, or no angle at all).
fn non_orthogonality(d: Vec3, area: Vec3) -> (f64, bool) {
    let cos = d.dot(area) / (d.length() * area.length());
    let degrees = if cos.is_finite() {
        cos.clamp(-1.0, 1.0).acos().to_degrees()
    } else {
        180.0
    };
    (degrees, cos.is_nan() || cos <= 0.0)
}

/// Cells per chunk of the parallel checks.
const CHECK_CHUNK: usize = 2048;
/// Cells per chunk of a re-check after splits (a few thousand cells).
const RECHECK_CHUNK: usize = 256;

/// The solver rules and each cell's volume as the solver takes it (the sum
/// of its face pyramids), with the failing cells and the cells they share a
/// face with.
#[allow(clippy::type_complexity)]
fn solver_rules_and_volumes(
    p: &[Vec3],
    cells: &[Cell],
) -> Result<(SolverRules, Vec<f64>, Vec<(usize, Vec<usize>)>), String> {
    let table = FaceTable::new(cells)?;
    solver_rules_on(p, cells, &table)
}

#[allow(clippy::type_complexity)]
fn solver_rules_on(
    p: &[Vec3],
    cells: &[Cell],
    table: &FaceTable,
) -> Result<(SolverRules, Vec<f64>, Vec<(usize, Vec<usize>)>), String> {
    // Each cell as exported: the faces between two cells merged into one
    // polygon (the hand-off and FoehnCraft allow one face per pair).
    struct Part {
        bad: Vec<bool>,
        volumes: Vec<f64>,
        centres: Vec<Vec3>,
        counts: [usize; 3],
        warping: f64,
        // (owner, neighbour, owner's area vector) once per pair
        shared: Vec<(usize, usize, Vec3)>,
        boundary: Vec<(usize, Vec3, Vec3)>,
    }
    let parts = chunked(cells.len(), CHECK_CHUNK, |range| {
        let mut part = Part {
            bad: Vec::with_capacity(range.len()),
            volumes: Vec::with_capacity(range.len()),
            centres: Vec::with_capacity(range.len()),
            counts: [0; 3],
            warping: 0.0,
            shared: Vec::new(),
            boundary: Vec::new(),
        };
        let mut cell = Exported::default();
        for id in range {
            let verdict = cell.evaluate(p, table, id, &cells[id].element);
            part.warping = part.warping.max(verdict.warping);
            part.counts[0] += usize::from(verdict.open);
            part.counts[1] += usize::from(verdict.nonpositive());
            part.counts[2] += usize::from(!verdict.pyramids_ok);
            part.bad.push(verdict.bad());
            part.volumes.push(verdict.volume);
            part.centres.push(verdict.centre);
            for (n, &(area, c)) in cell.across.iter().zip(&cell.geometry) {
                match *n {
                    Some(other) if other > id => part.shared.push((id, other, area)),
                    Some(_) => {}
                    None => part.boundary.push((id, area, c)),
                }
            }
        }
        part
    })?;
    let mut out = SolverRules {
        cells: cells.len(),
        ..SolverRules::default()
    };
    let mut bad = Vec::with_capacity(cells.len());
    let mut volumes = Vec::with_capacity(cells.len());
    let mut centres = Vec::with_capacity(cells.len());
    for part in &parts {
        out.maximum_warping = out.maximum_warping.max(part.warping);
        out.open_cells += part.counts[0];
        out.nonpositive_volume += part.counts[1];
        out.negative_pyramids += part.counts[2];
        bad.extend_from_slice(&part.bad);
        volumes.extend_from_slice(&part.volumes);
        centres.extend_from_slice(&part.centres);
    }
    let centres = &centres;
    let turned_parts = chunked(parts.len(), 1, |range| {
        let mut worst = 0.0_f64;
        let mut turned = Vec::new();
        for part in &parts[range] {
            for &(owner, neighbour, area) in &part.shared {
                let (degrees, back) = non_orthogonality(centres[neighbour] - centres[owner], area);
                worst = worst.max(degrees);
                if back {
                    turned.extend([owner, neighbour]);
                }
            }
            for &(owner, area, c) in &part.boundary {
                let (degrees, back) = non_orthogonality(c - centres[owner], area);
                worst = worst.max(degrees);
                if back {
                    turned.push(owner);
                }
            }
        }
        (worst, turned)
    })?;
    let mut turned = vec![false; cells.len()];
    for (worst, ids) in turned_parts {
        out.maximum_non_orthogonality = out.maximum_non_orthogonality.max(worst);
        for id in ids {
            turned[id] = true;
        }
    }
    out.faces_turned_away = turned.iter().filter(|&&t| t).count();
    let mut failing = Vec::new();
    let mut buffer = Rings::default();
    for (id, ((b, t), cell)) in bad.iter().zip(&turned).zip(cells).enumerate() {
        if *b || *t {
            out.failing_cells += 1;
            out.failing_layer_cells += usize::from(cell.layer > 0);
            failing.push((id, table.neighbours(id, &cell.element, &mut buffer)));
        }
    }
    Ok((out, volumes, failing))
}

/// A merged face: the cell across (none on the boundary), area vector and
/// centre.
type MergedFace = (Option<usize>, Vec3, Vec3);

/// The solver rules on the mesh as exported, kept up to date while cells
/// are split: after a change only the cells whose verdict can change are
/// checked again (those within two face neighbours of a changed cell), with
/// the same verdicts as a check of the whole mesh.
pub(super) struct ExportedCheck {
    table: FaceTable,
    /// The failing cells, each with the cells it shares a face with.
    failing: BTreeMap<usize, Vec<usize>>,
    /// Cells whose faces changed, and the cells that share those faces.
    touched: BTreeSet<usize>,
}

impl ExportedCheck {
    pub(super) fn new(p: &[Vec3], cells: &[Cell]) -> Result<Self, String> {
        let table = FaceTable::new(cells)?;
        let (rules, _, failing) = solver_rules_on(p, cells, &table)?;
        if crate::trace_enabled() && rules.failing_cells > 0 {
            eprintln!(
                "kofl trace: first check: {} of {} cells fail ({} in layers; open {}, volume {}, pyramids {}, turned away {})",
                rules.failing_cells,
                rules.cells,
                rules.failing_layer_cells,
                rules.open_cells,
                rules.nonpositive_volume,
                rules.negative_pyramids,
                rules.faces_turned_away
            );
        }
        Ok(Self {
            table,
            failing: failing.into_iter().collect(),
            touched: BTreeSet::new(),
        })
    }

    pub(super) fn failing(&self) -> &BTreeMap<usize, Vec<usize>> {
        &self.failing
    }

    /// Call with the old element before cell `id` changes.
    pub(super) fn remove(&mut self, id: usize, element: &Element) {
        let others = self.table.update(id, element, false);
        self.touched.insert(id);
        self.touched.extend(others);
    }

    /// Call with the new element after cell `id` changed or was added.
    pub(super) fn insert(&mut self, id: usize, element: &Element) {
        let others = self.table.update(id, element, true);
        self.touched.insert(id);
        self.touched.extend(others);
    }

    /// Checks again every cell whose verdict the changes since the last call
    /// can have changed, on all cores.
    pub(super) fn recheck(&mut self, p: &[Vec3], cells: &[Cell]) -> Result<(), String> {
        let touched: Vec<usize> = std::mem::take(&mut self.touched).into_iter().collect();
        let table = &self.table;
        // the face neighbours of each cell of `ids`, in that order
        let around = |ids: &[usize]| -> Result<Vec<Vec<usize>>, String> {
            let parts = chunked(ids.len(), RECHECK_CHUNK, |range| {
                let mut buffer = Rings::default();
                ids[range]
                    .iter()
                    .map(|&s| table.neighbours(s, &cells[s].element, &mut buffer))
                    .collect::<Vec<_>>()
            })?;
            Ok(parts.into_iter().flatten().collect())
        };
        let mut dirty: BTreeSet<usize> = touched.iter().copied().collect();
        dirty.extend(around(&touched)?.into_iter().flatten());
        let dirty: Vec<usize> = dirty.into_iter().collect();
        let next = around(&dirty)?;
        let mut needed: BTreeSet<usize> = dirty.iter().copied().collect();
        needed.extend(next.iter().flatten().copied());
        let needed: Vec<usize> = needed.into_iter().collect();
        let parts = chunked(needed.len(), RECHECK_CHUNK, |range| {
            let mut exported = Exported::default();
            needed[range]
                .iter()
                .map(|&id| {
                    let v = exported.evaluate(p, table, id, &cells[id].element);
                    let faces: Vec<MergedFace> = exported
                        .across
                        .iter()
                        .zip(&exported.geometry)
                        .map(|(&n, &(area, c))| (n, area, c))
                        .collect();
                    (v, faces)
                })
                .collect::<Vec<_>>()
        })?;
        let verdicts: Vec<(Verdict, Vec<MergedFace>)> = parts.into_iter().flatten().collect();
        let verdict = |id: usize| needed.binary_search(&id).ok().map(|i| &verdicts[i]);
        for (&x, neighbours) in dirty.iter().zip(next) {
            let Some((own, faces)) = verdict(x) else {
                return Err(format!("cell {x} missing from the re-check"));
            };
            let mut fails = own.bad();
            for &(n, area, c) in faces {
                fails |= match n {
                    None => non_orthogonality(c - own.centre, area).1,
                    Some(y) if y > x => verdict(y).is_none_or(|(theirs, _)| {
                        non_orthogonality(theirs.centre - own.centre, area).1
                    }),
                    Some(_) => false,
                };
            }
            // faces whose owner is a neighbour with a lower id
            for &y in neighbours.iter().filter(|&&y| y < x) {
                fails |= verdict(y).is_none_or(|(theirs, their_faces)| {
                    their_faces.iter().any(|&(n, area, _)| {
                        n == Some(x) && non_orthogonality(own.centre - theirs.centre, area).1
                    })
                });
            }
            if fails {
                self.failing.insert(x, neighbours);
            } else {
                self.failing.remove(&x);
            }
        }
        Ok(())
    }
}

fn det(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    a.cross(b).dot(c)
}

/// At fixed height the linear-wedge determinant is affine over the reference
/// triangle. On each of its three vertical edges it is quadratic in height.
/// Integrate that polynomial exactly to obtain the physical wedge volume.
pub(super) fn wedge(p: [Vec3; 6]) -> (f64, f64, [f64; 2]) {
    let [a, b, c, d, e, f] = p;
    let (u, v) = (b - a, c - a);
    let (du, dv) = ((e - d) - u, (f - d) - v);
    let normal = u.cross(v) / u.cross(v).length();
    let mut minimum = f64::INFINITY;
    let mut volume = 0.0;
    let mut thickness = [f64::INFINITY, f64::NEG_INFINITY];
    for w in [d - a, e - b, f - c] {
        let (q, l, h) = (det(du, dv, w), det(du, v, w) + det(u, dv, w), det(u, v, w));
        minimum = minimum.min(h).min(q + l + h);
        if q > 0.0 {
            let t = -l / (2.0 * q);
            if (0.0..1.0).contains(&t) {
                minimum = minimum.min((q * t + l) * t + h);
            }
        }
        volume += (h + l * 0.5 + q / 3.0) / 6.0;
        thickness[0] = thickness[0].min(w.dot(normal));
        thickness[1] = thickness[1].max(w.dot(normal));
    }
    (minimum, volume, thickness)
}

fn same_ring(a: &[usize], b: &[usize], reverse: bool) -> bool {
    let n = a.len();
    b.len() == n
        && a.iter().position(|v| *v == b[0]).is_some_and(|start| {
            (0..n).all(|i| a[(start + if reverse { n - i } else { i }) % n] == b[i])
        })
}
/// A face counts as planar when its corners are within this share of its
/// size of the fitted plane. Such a face is fanned from its vertex mean, as
/// FoehnCraft does with every face, so the checked volume is the solver's.
/// The core welds points up to 1e-10 of the model size apart, which moves a
/// corner of a small face by more than 1e-10 of that face's size.
const PLANAR_SHARE: f64 = 1e-7;

/// Two-point Gauss rule on [0, 1].
const GAUSS: [f64; 2] = [0.211_324_865_405_187_1, 0.788_675_134_594_812_9];

/// A face integrated from `origin`: (six times the cone volume, area), and
/// whether every part of it is seen strictly from the front.
///
/// Planar polygons are fanned from their vertex mean, as FoehnCraft does. A
/// quadrilateral that is not planar is the bilinear patch of a layer prism
/// side (as in a linear wedge); its cone volume is integrated exactly by the
/// 2×2 Gauss rule (the integrand has degree two per parameter).
fn face_integral(
    points: &[Vec3],
    face: &[usize],
    origin: Vec3,
) -> Result<(f64, f64, bool), &'static str> {
    let (newell, centre) = super::convex::ring_plane(points, face);
    let norm = newell.length();
    if norm.is_nan() || norm <= 0.0 {
        return Err("degenerate polyhedron face");
    }
    let unit = super::convex::fitted_normal(points, face, newell, centre);
    let diameter = face
        .iter()
        .map(|&v| (points[v] - centre).length())
        .fold(0.0, f64::max);
    let coordinate_scale = centre.max_abs_component().max(diameter);
    // a triangle is planar; on a sliver its Newell normal is too inaccurate
    // for the offset test
    let planar = face.len() == 3
        || face.iter().all(|&v| {
            (points[v] - centre).dot(unit).abs()
                <= PLANAR_SHARE * diameter + 256.0 * f64::EPSILON * coordinate_scale
        });
    if !planar {
        if face.len() != 4 {
            return Err("nonplanar polyhedron face");
        }
        let [p0, p1, p2, p3] = [0, 1, 2, 3].map(|i| points[face[i]]);
        let mut v6 = 0.0;
        let mut area = 0.0;
        let mut front = true;
        let at = |u: f64, v: f64| {
            let x = p0 * ((1. - u) * (1. - v))
                + p1 * (u * (1. - v))
                + p2 * (u * v)
                + p3 * ((1. - u) * v);
            let xu = (p1 - p0) * (1. - v) + (p2 - p3) * v;
            let xv = (p3 - p0) * (1. - u) + (p2 - p1) * u;
            (x, xu.cross(xv))
        };
        for &u in &GAUSS {
            for &v in &GAUSS {
                let (x, n) = at(u, v);
                // 2 = 6 / 3 for the cone, 1/4 for the Gauss weights
                v6 += 0.5 * (x - origin).dot(n);
                area += 0.25 * n.length();
                front &= (x - origin).dot(n) > 0.0;
            }
        }
        for (u, v) in [(0., 0.), (1., 0.), (1., 1.), (0., 1.)] {
            let (x, n) = at(u, v);
            front &= (x - origin).dot(n) > 0.0;
        }
        return Ok((v6, area, front));
    }
    let mut v6 = 0.0;
    let mut area = 0.0;
    let mut front = true;
    for k in 0..face.len() {
        let (a, b) = (points[face[k]], points[face[(k + 1) % face.len()]]);
        let cross = (a - centre).cross(b - centre);
        let signed = cross.dot(centre - origin);
        front &= signed > 0.0 && cross.dot(newell) > 0.0;
        v6 += signed;
        area += 0.5 * cross.length();
    }
    Ok((v6, area, front))
}

/// Integrate the oriented shell from its certified interior point. The expected
/// volume is computed independently from the source tetrahedron partitions.
/// Every edge used twice in opposite directions, Euler characteristic 2.
fn shell_closed(poly: &Polyhedron) -> Result<(), &'static str> {
    let mut edges = HashMap::<(usize, usize), (usize, i32)>::new();
    for face in &poly.faces {
        for (&a, &b) in face.iter().zip(face.iter().cycle().skip(1)) {
            let e = edges.entry((a.min(b), a.max(b))).or_default();
            e.0 += 1;
            e.1 += if a < b { 1 } else { -1 };
        }
    }
    if edges.values().any(|&use_count| use_count != (2, 0))
        || poly.nodes.len() + poly.faces.len() != edges.len() + 2
    {
        return Err("open or nonmanifold polyhedron shell");
    }
    Ok(())
}

fn poly_volume(poly: &Polyhedron, points: &[Vec3]) -> Result<(f64, f64), &'static str> {
    let mut volume = 0.0;
    let mut area = 0.0;
    for face in &poly.faces {
        let (v6, a, front) = face_integral(points, face, poly.interior)?;
        if !v6.is_finite() || !front {
            return Err("inverted face or invalid polyhedron kernel");
        }
        volume += v6 / 6.0;
        area += a;
    }
    shell_closed(poly)?;
    if !volume.is_finite() || (volume / poly.expected_volume - 1.0).abs() > 1e-8 {
        return Err("polyhedron shell volume differs from its source partition");
    }
    Ok((
        volume,
        (36. * std::f64::consts::PI).cbrt() * volume.cbrt().powi(2) / area,
    ))
}

pub(super) fn check(
    surface: &ClosedSurface,
    patches: &[usize],
    p: &[Vec3],
    cells: &[Cell],
    boundary: &[BoundaryFace],
    layers: usize,
) -> Result<Measurements, String> {
    check_with_interface_count(
        surface,
        patches,
        p,
        cells,
        boundary,
        layers,
        patches.iter().filter(|&&p| super::is_wall(p)).count(),
    )
}

/// Area of a boundary face from its vector area. For a planar face this is
/// the area also when the face is not convex, where a sum of unsigned fan
/// triangles counts a notch twice.
fn boundary_area(p: &[Vec3], nodes: &[usize]) -> f64 {
    let a = p[nodes[0]];
    nodes[1..]
        .windows(2)
        .fold(Vec3::ZERO, |s, w| s + (p[w[0]] - a).cross(p[w[1]] - a))
        .length()
        * 0.5
}

/// Pairs every face: each boundary face with its owner and the same winding,
/// each interior face with exactly one other cell and the opposite winding.
/// Returns the number of faces between a layer cell and a core cell. The
/// faces are split into shards by hash, one thread each.
fn pair_faces(cells: &[Cell], boundary: &[BoundaryFace]) -> Result<usize, String> {
    let shards = threads();
    // sized for every face up front: growing the maps cost more than the
    // lookups on a million-cell mesh
    let faces: usize = cells.iter().map(|c| face_count(&c.element)).sum();
    let results = chunked(shards, 1, |range| -> Result<usize, String> {
        let shard = range.start;
        let mine = |k: &FaceKey| FaceTable::shard_of(k.hash, shards) == shard;
        let mut exterior: KeyMap<&BoundaryFace> = KeyMap::default();
        for face in boundary {
            let k = FaceKey::new(&face.nodes);
            if mine(&k) && exterior.insert(k, face).is_some() {
                return Err("duplicate outer boundary face".into());
            }
        }
        // (cell, face index, layer) of a face met once so far
        let mut pending: KeyMap<(usize, usize, usize)> =
            KeyMap::with_capacity_and_hasher(faces / 8 / shards + 16, Default::default());
        let mut complete =
            KeySet::with_capacity_and_hasher(faces / 2 / shards + 16, Default::default());
        let (mut rings, mut earlier) = (Rings::default(), Rings::default());
        let mut prism_core_faces = 0;
        for (id, cell) in cells.iter().enumerate() {
            rings.load(&cell.element);
            for (index, ring) in rings.iter().enumerate() {
                let k = FaceKey::new(ring);
                if !mine(&k) {
                    continue;
                }
                if complete.contains(&k) {
                    return Err("too many incidences on one face".into());
                }
                if let Some(face) = exterior.remove(&k) {
                    complete.insert(k);
                    if id != face.owner || !same_ring(&face.nodes, ring, false) {
                        return Err("incorrect boundary owner or winding".into());
                    }
                } else if let Some((previous, at, other_layer)) = pending.remove(&k) {
                    if (other_layer == 0) != (cell.layer == 0) {
                        prism_core_faces += 1;
                    }
                    complete.insert(k);
                    if !same_ring(
                        ring_of(&cells[previous].element, at, &mut earlier),
                        ring,
                        true,
                    ) {
                        return Err("interior face has inconsistent winding".into());
                    }
                } else {
                    pending.insert(k, (id, index, cell.layer));
                }
            }
        }
        if !pending.is_empty() || !exterior.is_empty() {
            return Err("unpaired interior or missing exterior face".into());
        }
        Ok(prism_core_faces)
    })?;
    results.into_iter().sum()
}

pub(super) fn check_with_interface_count(
    surface: &ClosedSurface,
    patches: &[usize],
    p: &[Vec3],
    cells: &[Cell],
    boundary: &[BoundaryFace],
    layers: usize,
    expected_interfaces: usize,
) -> Result<Measurements, String> {
    let started = std::time::Instant::now();
    let (solver_rules, solver_volumes, _) = solver_rules_and_volumes(p, cells)?;
    if crate::trace_enabled() {
        eprintln!(
            "kofl trace: solver rules ({:.2} s): {} of {} cells fail ({} in layers; open {}, volume {}, pyramids {}, turned away {}), max non-orthogonality {:.2}°",
            started.elapsed().as_secs_f64(),
            solver_rules.failing_cells,
            solver_rules.cells,
            solver_rules.failing_layer_cells,
            solver_rules.open_cells,
            solver_rules.nonpositive_volume,
            solver_rules.negative_pyramids,
            solver_rules.faces_turned_away,
            solver_rules.maximum_non_orthogonality
        );
    }
    if solver_rules.failing_cells > 0 {
        // Name the likely cause in plain words when layer cells fail.
        let hint = if solver_rules.failing_layer_cells > 0 {
            "; wall layer collision or fold: thinner or fewer layers help"
        } else {
            ""
        };
        return Err(format!(
            "{} of {} cells fail the solver rules ({} in layers; open {}, volume {}, pyramids {}, turned away {}){hint}",
            solver_rules.failing_cells,
            solver_rules.cells,
            solver_rules.failing_layer_cells,
            solver_rules.open_cells,
            solver_rules.nonpositive_volume,
            solver_rules.negative_pyramids,
            solver_rules.faces_turned_away
        ));
    }
    let ruled = std::time::Instant::now();
    let ids = patches.iter().copied().max().unwrap_or(0).max(super::WALL);
    let mut areas = vec![0.; ids];
    let mut port_layer_areas = [vec![0.; layers + 1], vec![0.; layers + 1]];
    for face in boundary {
        let area = boundary_area(p, &face.nodes);
        *areas
            .get_mut(face.patch.wrapping_sub(1))
            .ok_or("boundary face with an unknown patch")? += area;
        if face.patch <= 2 {
            port_layer_areas[face.patch - 1][cells[face.owner].layer] += area;
        }
    }
    let mut expected = vec![0.; ids];
    let scale = surface.bounds().extent().max_abs_component();
    for (facet, &patch) in surface.facets().iter().zip(patches) {
        expected[patch - 1] += facet.normalized_area * scale * scale;
    }
    for (i, (&a, &e)) in areas.iter().zip(&expected).enumerate() {
        if a == 0. && e == 0. {
            continue; // an id no facet uses (no inlet cap, say)
        }
        if !a.is_finite() || !e.is_finite() || e <= 0. || (a / e - 1.).abs() >= 1e-9 {
            return Err(format!(
                "patch {} area {a} differs from STL area {e}",
                i + 1
            ));
        }
    }
    // The solver rules decided above; the element measures are quality
    // reports only, taken on all cores.
    struct Report {
        minimum_j: f64,
        minimum_kernel_clearance: f64,
        normal_thickness_ranges: Vec<[f64; 2]>,
        minimum_tet_mean_ratio: f64,
        tetrahedra_below_mean_ratio_0_1: usize,
        sphericities: Vec<f64>,
    }
    let reports = chunked(cells.len(), CHECK_CHUNK, |range| {
        let mut r = Report {
            minimum_j: f64::INFINITY,
            minimum_kernel_clearance: f64::INFINITY,
            normal_thickness_ranges: vec![[f64::INFINITY, 0.0_f64]; layers],
            minimum_tet_mean_ratio: 1.0_f64,
            tetrahedra_below_mean_ratio_0_1: 0,
            sphericities: Vec::new(),
        };
        for cell in &cells[range] {
            match cell.element {
                Element::Polyhedron(ref poly) if poly.prism.is_some() => {
                    for w in poly.prism.iter().flatten() {
                        let (jac, _, range) = wedge(w.map(|v| p[v]));
                        r.minimum_j = r.minimum_j.min(jac);
                        if let Some(t) = cell
                            .layer
                            .checked_sub(1)
                            .and_then(|l| r.normal_thickness_ranges.get_mut(l))
                        {
                            t[0] = t[0].min(range[0]);
                            t[1] = t[1].max(range[1]);
                        }
                    }
                }
                Element::Polyhedron(ref poly) => {
                    if cell.layer == 0
                        && let Ok((_, sphericity)) = poly_volume(poly, p)
                    {
                        r.sphericities.push(sphericity);
                    }
                }
                Element::Tet(n) => {
                    let [o, a, b, c] = n.map(|v| p[v]);
                    let determinant = det(b - a, c - a, a - o);
                    let clearance = determinant / (b - a).cross(c - a).length();
                    r.minimum_kernel_clearance = r.minimum_kernel_clearance.min(clearance);
                    let squared_edges = [b - a, c - b, a - c, a - o, b - o, c - o]
                        .iter()
                        .map(|e| e.dot(*e))
                        .sum::<f64>();
                    let mean_ratio = 12.0 * (determinant / 2.0).cbrt().powi(2) / squared_edges;
                    r.minimum_tet_mean_ratio = r.minimum_tet_mean_ratio.min(mean_ratio);
                    r.tetrahedra_below_mean_ratio_0_1 += usize::from(mean_ratio < 0.1);
                }
                Element::Wedge(n) => {
                    let (jac, _, range) = wedge(n.map(|v| p[v]));
                    r.minimum_j = r.minimum_j.min(jac);
                    let t = &mut r.normal_thickness_ranges[cell.layer - 1];
                    t[0] = t[0].min(range[0]);
                    t[1] = t[1].max(range[1]);
                }
            }
        }
        r
    })?;
    let mut minimum_j = f64::INFINITY;
    let mut minimum_kernel_clearance = f64::INFINITY;
    let mut normal_thickness_ranges = vec![[f64::INFINITY, 0.0_f64]; layers];
    let mut minimum_tet_mean_ratio = 1.0_f64;
    let mut tetrahedra_below_mean_ratio_0_1 = 0;
    let mut sphericities = Vec::new();
    for r in reports {
        minimum_j = minimum_j.min(r.minimum_j);
        minimum_kernel_clearance = minimum_kernel_clearance.min(r.minimum_kernel_clearance);
        for (all, part) in normal_thickness_ranges
            .iter_mut()
            .zip(&r.normal_thickness_ranges)
        {
            all[0] = all[0].min(part[0]);
            all[1] = all[1].max(part[1]);
        }
        minimum_tet_mean_ratio = minimum_tet_mean_ratio.min(r.minimum_tet_mean_ratio);
        tetrahedra_below_mean_ratio_0_1 += r.tetrahedra_below_mean_ratio_0_1;
        sphericities.extend(r.sphericities);
    }
    let measured = std::time::Instant::now();
    let mut minimum_volume = f64::INFINITY;
    let mut volume = 0.;
    for (id, &value) in solver_volumes.iter().enumerate() {
        if !value.is_finite() || value <= 0. {
            return Err(format!(
                "nonpositive or unrepresentable volume at cell {id}"
            ));
        }
        volume += value;
        minimum_volume = minimum_volume.min(value);
    }
    let prism_core_faces = pair_faces(cells, boundary)?;
    if crate::trace_enabled() {
        eprintln!(
            "kofl trace: measures {:.2} s, face pairing {:.2} s",
            (measured - ruled).as_secs_f64(),
            measured.elapsed().as_secs_f64()
        );
    }
    if !volume.is_finite() || (volume / surface.volume() - 1.).abs() >= 1e-9 {
        return Err(format!(
            "mesh volume {volume} differs from STL volume {}",
            surface.volume()
        ));
    }
    // every layer reaches each cap the surface has
    let has_cap = |id: usize| patches.contains(&id);
    if port_layer_areas
        .iter()
        .enumerate()
        .any(|(i, a)| has_cap(i + 1) && a.iter().any(|v| *v <= 0.))
    {
        return Err("a layer does not reach both ports".into());
    }
    if layers > 0 && prism_core_faces != expected_interfaces {
        return Err("incomplete prism-to-core interface".into());
    }
    sphericities.sort_unstable_by(f64::total_cmp);
    let core_sphericity = (!sphericities.is_empty()).then(|| {
        let n = sphericities.len() - 1;
        [sphericities[0], sphericities[n / 20], sphericities[n / 2]]
    });
    Ok(Measurements {
        core_sphericity,
        prism_core_faces,
        minimum_tet_mean_ratio,
        tetrahedra_below_mean_ratio_0_1,
        volume,
        minimum_volume,
        minimum_wedge_jacobian: minimum_j.is_finite().then_some(minimum_j),
        areas,
        port_layer_areas,
        normal_thickness_ranges,
        minimum_kernel_clearance,
        solver_rules,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cell(element: Element) -> Cell {
        Cell {
            element,
            layer: 0,
            region: 0,
        }
    }

    #[test]
    fn warped_faces_pass_the_solver_rules() {
        // a wedge whose top corner is lifted: its three side quads are not
        // flat, which the solver rules accept
        let p = [
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [1., 0., 1.3],
            [0., 1., 1.],
        ]
        .map(|[x, y, z]| Vec3::new(x, y, z));
        let r = solver_rules(&p, &[cell(Element::Wedge([0, 1, 2, 3, 4, 5]))]);
        assert_eq!(r.failing_cells, 0, "{r:?}");
        assert!(r.maximum_non_orthogonality < 90.0);
    }

    #[test]
    fn two_tetrahedra_pass_and_an_inverted_one_fails() {
        let p = [
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [0.6, 0.6, 0.6],
        ]
        .map(|[x, y, z]| Vec3::new(x, y, z));
        let good = [
            cell(Element::Tet([0, 1, 2, 3])),
            cell(Element::Tet([1, 2, 3, 4])),
        ];
        // the second tetrahedron must be oriented positively
        let r = solver_rules(&p, &good);
        assert_eq!((r.cells, r.failing_cells), (2, 0), "{r:?}");
        let inverted = [cell(Element::Tet([0, 2, 1, 3]))];
        let r = solver_rules(&p, &inverted);
        assert_eq!(r.failing_cells, 1);
        assert_eq!(r.nonpositive_volume, 1);
    }

    #[test]
    fn wedge_with_positive_end_jacobians_can_invert_halfway() {
        let points = [
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [-2., 0., 1.],
            [0., -1., 1.],
        ]
        .map(|[x, y, z]| Vec3::new(x, y, z));
        // J(t) = (1 - 3t)(1 - 2t), with a negative minimum at t = 5/12.
        let (minimum, volume, _) = wedge(points);
        assert!((minimum + 1.0 / 24.0).abs() < 1e-14);
        assert!(volume > 0.0);
    }

    /// Keep rule: a long thin planar face far from the origin passes the
    /// planarity check. The ring is a core face of the Manifold as built
    /// with conforming wall seeds (2.6 long, 0.00007 wide, 263 from the
    /// origin): two wall cells that barely touch over a wall vertex. Its
    /// vertices are on one plane to 1e-13, but a Newell normal from
    /// absolute coordinates is tilted by 1e-10 over the length and the
    /// offset test rejected the whole mesh. Taking the normal from a plain
    /// Newell sum again, or dropping the least-squares plane, breaks this;
    /// the fixtures of the fast run have no such face, and the Manifold
    /// build is a full-run test.
    #[test]
    fn thin_face_far_from_the_origin_is_planar() {
        let points = [
            [6.414802894970417, 262.9597716879253, 18.751259654212134],
            [6.859459004309515, 262.6533452053385, 18.034802806686436],
            [6.983285982556163, 262.5680121477317, 17.835285281547556],
            [6.997576348350637, 262.5581642080918, 17.812259778890425],
            [7.023771642070697, 262.54011220892545, 17.770052335510094],
            [7.023766730368811, 262.54005614619194, 17.77009088156666],
            [6.997572759971385, 262.5581082208852, 17.812295683907582],
            [6.983282513930192, 262.56795661671845, 17.8353207160207],
            [6.859456505831417, 262.6532936738063, 18.034834272658962],
            [6.414835399147537, 262.9597127964614, 18.751226084952386],
            [5.4468412280837315, 263.62682445664535, 20.310898102481982],
        ]
        .map(|[x, y, z]| Vec3::new(x, y, z));
        let ring: Vec<usize> = (0..points.len()).collect();
        let (newell, centre) = crate::duct::convex::ring_plane(&points, &ring);
        let unit = crate::duct::convex::fitted_normal(&points, &ring, newell, centre);
        let flatness = points
            .iter()
            .map(|&p| (p - centre).dot(unit).abs())
            .fold(0.0, f64::max);
        assert!(
            flatness < 1e-12,
            "offset from the fitted plane {flatness:e}"
        );
        // seen from a point off the plane the face integrates as planar
        let origin = centre - unit;
        let result = face_integral(&points, &ring, origin);
        let (v6, area, _) = result.expect("the face is planar");
        assert!(area > 5e-5 && area < 2e-4, "area {area:e}");
        assert!((v6.abs() / 6.0 - area / 3.0).abs() < 1e-9 * area);
    }

    /// Keep rule: the fitted plane does not hide a face that is bent. A
    /// pentagon with one corner lifted by a millionth of its size is not
    /// planar, with or without the fit; a fit that absorbed the lift (a
    /// wrong sign or axis in the tilt) would accept warped polyhedron faces
    /// everywhere, and no mesh test builds one.
    #[test]
    fn bent_face_is_still_rejected() {
        let mut points = [[0., 0.], [1., 0.], [1.4, 0.8], [0.5, 1.5], [-0.4, 0.8]]
            .map(|[x, y]| Vec3::new(100.0 + x, 200.0 + y, 50.0));
        points[2].z += 1e-6;
        let ring: Vec<usize> = (0..points.len()).collect();
        let origin = Vec3::new(100.5, 200.7, 49.0);
        assert_eq!(
            face_integral(&points, &ring, origin).err(),
            Some("nonplanar polyhedron face")
        );
    }

    /// Invariant: the area of a notched planar boundary face is its true
    /// area. Summing unsigned fan triangles from the first corner (the code
    /// before) counts the notch twice, and the patch area check then fails
    /// a valid mesh; the fixtures only have convex boundary faces.
    #[test]
    fn notched_boundary_face_has_its_true_area() {
        let v = |x: f64, y: f64| Vec3::new(x, y, 0.0);
        let points = [v(0., 0.), v(4., 0.), v(4., 4.), v(2., 1.), v(0., 4.)];
        // square of 16 without the triangle (4,4) (2,1) (0,4) of area 6
        let area = boundary_area(&points, &[0, 1, 2, 3, 4]);
        assert!((area - 10.0).abs() < 1e-12);
    }

    /// Eight tetrahedra around the centre of an octahedron, all positive.
    fn octahedron() -> (Vec<Vec3>, Vec<Cell>) {
        let p = vec![
            Vec3::new(0., 0., 0.),
            Vec3::new(1., 0., 0.),
            Vec3::new(-1., 0., 0.),
            Vec3::new(0., 1., 0.),
            Vec3::new(0., -1., 0.),
            Vec3::new(0., 0., 1.),
            Vec3::new(0., 0., -1.),
        ];
        let mut cells = Vec::new();
        for x in [1, 2] {
            for y in [3, 4] {
                for z in [5, 6] {
                    let mut t = [0, x, y, z];
                    let [o, a, b, c] = t.map(|i| p[i]);
                    if (a - o).cross(b - o).dot(c - o) < 0.0 {
                        t.swap(2, 3);
                    }
                    cells.push(cell(Element::Tet(t)));
                }
            }
        }
        (p, cells)
    }

    /// Invariant: the faces the checks read are the faces of the element,
    /// in the same order (the order decides the fan origin of the volume).
    /// Reordering a wedge's faces in `Rings::load` breaks this; the mesh
    /// tests would only see slightly different volumes.
    #[test]
    fn flat_rings_are_the_element_faces() {
        let poly = Polyhedron::new(vec![vec![0, 1, 2], vec![2, 1, 3, 4]], Vec3::ZERO, None, 1.0);
        let mut rings = Rings::default();
        for element in [
            Element::Tet([3, 1, 4, 1]),
            Element::Wedge([5, 9, 2, 6, 5, 3]),
            Element::Polyhedron(Box::new(poly)),
        ] {
            rings.load(&element);
            let flat: Vec<Vec<usize>> = rings.iter().map(<[usize]>::to_vec).collect();
            assert_eq!(flat, element.faces());
        }
    }

    /// Invariant: two faces that share an edge merge into one ring starting
    /// at the smallest vertex, and faces that do not form one loop do not
    /// merge. A merge that kept a shared edge, or started elsewhere, would
    /// change the exported faces; the duct tests see that only as changed
    /// cell counts.
    #[test]
    fn rings_merge_into_one_loop() {
        let a = vec![4, 5, 6, 7];
        let b = vec![7, 6, 8, 9];
        assert_eq!(merge_rings(&[&a, &b]), Some(vec![4, 5, 6, 8, 9, 7]));
        let apart = vec![10, 11, 12];
        assert_eq!(merge_rings(&[&a, &apart]), None);
    }

    /// Invariant: after cells change, the local re-check reports exactly the
    /// failing cells (and their neighbours) of a check of the whole mesh.
    /// Re-checking only the changed cells, not the cells within two face
    /// neighbours, breaks this: a neighbour whose merged face or centre
    /// moved keeps its old verdict. The duct builds would show it only as a
    /// slightly different cell count, or a final check that fails.
    #[test]
    fn local_recheck_matches_the_whole_mesh() {
        let (mut p, mut cells) = octahedron();
        let mut check = ExportedCheck::new(&p, &cells).expect("check");
        let whole = |p: &[Vec3], cells: &[Cell]| -> Vec<(usize, Vec<usize>)> {
            solver_rules_and_volumes(p, cells).expect("check").2
        };
        let local = |check: &ExportedCheck| -> Vec<(usize, Vec<usize>)> {
            check
                .failing()
                .iter()
                .map(|(&k, v)| (k, v.clone()))
                .collect()
        };
        assert_eq!(local(&check), whole(&p, &cells));
        assert!(local(&check).is_empty());
        // turn one tetrahedron inside out
        let Element::Tet([o, a, b, c]) = cells[3].element else {
            panic!("tetrahedron expected");
        };
        check.remove(3, &cells[3].element);
        cells[3].element = Element::Tet([o, b, a, c]);
        check.insert(3, &cells[3].element);
        check.recheck(&p, &cells).unwrap();
        assert!(!local(&check).is_empty());
        assert_eq!(local(&check), whole(&p, &cells));
        // split another at an edge midpoint: a new cell, faces that no
        // longer match their neighbours
        let Element::Tet([o, a, b, c]) = cells[5].element else {
            panic!("tetrahedron expected");
        };
        p.push((p[o] + p[a]) * 0.5);
        let m = p.len() - 1;
        check.remove(5, &cells[5].element);
        cells[5].element = Element::Tet([o, m, b, c]);
        check.insert(5, &cells[5].element);
        cells.push(cell(Element::Tet([m, a, b, c])));
        check.insert(cells.len() - 1, &cells[cells.len() - 1].element);
        check.recheck(&p, &cells).unwrap();
        assert_eq!(local(&check), whole(&p, &cells));
        // and back: the inverted cell restored
        check.remove(3, &cells[3].element);
        let Element::Tet(restored) = octahedron().1[3].element else {
            panic!("tetrahedron expected");
        };
        cells[3].element = Element::Tet(restored);
        check.insert(3, &cells[3].element);
        check.recheck(&p, &cells).unwrap();
        assert_eq!(local(&check), whole(&p, &cells));
    }
}
