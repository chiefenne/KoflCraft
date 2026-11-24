//! The binary hand-off to FoehnCraft and MurmelCraft: a blocked VTKHDF
//! file (layout 1, see `docs/foehncraft.md`), which a viewer can read block
//! by block and ParaView opens as an ordinary partitioned data set.
//!
//! - `/VTKHDF` is a VTKHDF 2.4 `PartitionedDataSetCollection`: the block
//!   `internal` holds the cells, one partition per block of
//!   [`BLOCK_CELLS`] cells in Hilbert order, and each boundary patch is a
//!   `PolyData` block of its own.
//! - Each partition repeats the points it uses and one layer of ghost
//!   cells (the cells of other blocks that share a point with it), so a
//!   block can be read alone.
//! - Every cell is a polyhedron (type 42) whose faces are listed per cell
//!   and wound outward from it.
//!
//! Written from the VTK documentation ("VTK File Formats", section VTKHDF,
//! version 2.4), the layout spec agreed with MurmelCraft, and Skilling,
//! "Programming the Hilbert curve", AIP Conf. Proc. 707 (2004), for the
//! order. The HDF5 container is [`crate::hdf5`].

use std::io::{self, Write};

use mesh_core::geometry::Vec3;

use crate::handoff::{CellFaceList, HandoffMesh, PatchRole};
use crate::hdf5::{DataWriter, Element, H5File, Obj};

/// Owned cells per block (the last block holds the rest).
pub const BLOCK_CELLS: usize = 65_536;
/// Bits per axis of the Hilbert index (63 bits in all).
const HILBERT_BITS: u32 = 21;
/// The version of the blocked layout, the root attribute `murmelcraft_blocks`.
pub const LAYOUT_VERSION: i64 = 1;
const NONE: u32 = u32::MAX;

/// The Hilbert index of a point of the 2²¹ grid: Skilling's transpose
/// algorithm (undo the excess work, Gray-encode), then the transposed
/// bits read from the most significant, x before y before z.
pub fn hilbert_index(mut x: [u32; 3]) -> u64 {
    let m = 1u32 << (HILBERT_BITS - 1);
    let mut q = m;
    while q > 1 {
        let p = q - 1;
        for i in 0..3 {
            if x[i] & q != 0 {
                x[0] ^= p;
            } else {
                let t = (x[0] ^ x[i]) & p;
                x[0] ^= t;
                x[i] ^= t;
            }
        }
        q >>= 1;
    }
    for i in 1..3 {
        x[i] ^= x[i - 1];
    }
    let mut t = 0;
    let mut q = m;
    while q > 1 {
        if x[2] & q != 0 {
            t ^= q - 1;
        }
        q >>= 1;
    }
    for v in &mut x {
        *v ^= t;
    }
    let mut h = 0u64;
    for b in (0..HILBERT_BITS).rev() {
        for v in x {
            h = (h << 1) | u64::from((v >> b) & 1);
        }
    }
    h
}

/// The distinct points of cell `c` in the order its faces first use them,
/// each face wound outward from the cell: the `Connectivity` of a
/// polyhedron in the layout.
fn cell_points(mesh: &HandoffMesh, cells: &CellFaceList, c: usize, points: &mut Vec<usize>) {
    points.clear();
    for (f, outward) in cells.cell(c) {
        let ring = mesh.face(f);
        let mut take = |p: usize| {
            if !points.contains(&p) {
                points.push(p);
            }
        };
        if outward {
            ring.iter().for_each(|&p| take(p));
        } else {
            ring.iter().rev().for_each(|&p| take(p));
        }
    }
}

/// The bounding box of points.
fn bounds(points: impl Iterator<Item = Vec3>) -> [f64; 6] {
    let mut b = [
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    ];
    for p in points {
        for (k, v) in [p.x, p.y, p.z].into_iter().enumerate() {
            b[k] = b[k].min(v);
            b[k + 3] = b[k + 3].max(v);
        }
    }
    b
}

/// How the cells are cut into blocks, and what each partition holds.
struct Plan {
    block_cells: usize,
    /// KoflCraft's cell number of each file cell id.
    order: Vec<usize>,
    /// The file cell id of each KoflCraft cell.
    file_id: Vec<usize>,
    /// The distinct points of each KoflCraft cell, counted.
    node_count: Vec<u32>,
    /// Ghost file cell ids of each block (`ghost_start[b]..ghost_start[b+1]`).
    ghost_start: Vec<usize>,
    ghosts: Vec<usize>,
    /// Points of each partition, as mesh point ids, in local id order.
    point_start: Vec<usize>,
    points: Vec<usize>,
    /// The lowest block whose owned cells use each point.
    point_owner: Vec<u32>,
    /// Per partition: connectivity ids, faces, face connectivity ids.
    connectivity: Vec<usize>,
    faces: Vec<usize>,
    face_connectivity: Vec<usize>,
    /// x, y, z minimum then maximum over the points of the owned cells.
    block_bounds: Vec<[f64; 6]>,
}

impl Plan {
    fn blocks(&self) -> usize {
        self.ghost_start.len() - 1
    }

    fn owned(&self, b: usize) -> std::ops::Range<usize> {
        b * self.block_cells..((b + 1) * self.block_cells).min(self.order.len())
    }

    /// The file cell ids of partition `b`: owned cells, then ghosts.
    fn cells(&self, b: usize) -> impl Iterator<Item = usize> + '_ {
        self.owned(b).chain(
            self.ghosts[self.ghost_start[b]..self.ghost_start[b + 1]]
                .iter()
                .copied(),
        )
    }

    fn partition_cells(&self, b: usize) -> usize {
        self.owned(b).len() + self.ghost_start[b + 1] - self.ghost_start[b]
    }

    fn partition_points(&self, b: usize) -> &[usize] {
        &self.points[self.point_start[b]..self.point_start[b + 1]]
    }

    fn new(mesh: &HandoffMesh, cells: &CellFaceList, block_cells: usize) -> Self {
        let n = mesh.cell_count();
        let block_cells = block_cells.max(1);
        let nodes = &mut Vec::new();

        // Hilbert order of the centroids (the mean of each cell's points)
        let bb = bounds(mesh.points.iter().copied());
        let scale = |v: f64, k: usize| -> u32 {
            let (lo, hi) = (bb[k], bb[k + 3]);
            if hi > lo {
                let max = f64::from((1u32 << HILBERT_BITS) - 1);
                ((v - lo) / (hi - lo) * max).clamp(0.0, max) as u32
            } else {
                0
            }
        };
        let mut node_count = Vec::with_capacity(n);
        let mut keys: Vec<(u64, usize)> = Vec::with_capacity(n);
        for c in 0..n {
            cell_points(mesh, cells, c, nodes);
            node_count.push(nodes.len() as u32);
            let mut s = Vec3::ZERO;
            for &v in nodes.iter() {
                s = s + mesh.points[v];
            }
            let m = s * (1.0 / nodes.len().max(1) as f64);
            keys.push((
                hilbert_index([scale(m.x, 0), scale(m.y, 1), scale(m.z, 2)]),
                c,
            ));
        }
        keys.sort_unstable();
        let order: Vec<usize> = keys.into_iter().map(|(_, c)| c).collect();
        let mut file_id = vec![0usize; n];
        for (f, &c) in order.iter().enumerate() {
            file_id[c] = f;
        }
        let blocks = n.div_ceil(block_cells);

        // The blocks whose owned cells use each point: the lowest in
        // `point_owner`, the others (ascending) as (point, block) pairs.
        let np = mesh.points.len();
        let mut point_owner = vec![NONE; np];
        let mut last = vec![NONE; np];
        let mut extra: Vec<(usize, u32)> = Vec::new();
        for (f, &c) in order.iter().enumerate() {
            let b = (f / block_cells) as u32;
            cell_points(mesh, cells, c, nodes);
            for &p in nodes.iter() {
                if point_owner[p] == NONE {
                    point_owner[p] = b;
                } else if last[p] != b {
                    extra.push((p, b));
                }
                last[p] = b;
            }
        }
        drop(last);
        extra.sort_unstable();
        let mut extra_start = vec![0usize; np + 1];
        for &(p, _) in &extra {
            extra_start[p + 1] += 1;
        }
        for p in 0..np {
            extra_start[p + 1] += extra_start[p];
        }

        // Ghosts: a cell is a ghost of every other block that uses one of
        // its points.
        let mut pairs: Vec<(u32, usize)> = Vec::new();
        for (f, &c) in order.iter().enumerate() {
            let b = (f / block_cells) as u32;
            cell_points(mesh, cells, c, nodes);
            for &p in nodes.iter() {
                let others = std::iter::once(point_owner[p]).chain(
                    extra[extra_start[p]..extra_start[p + 1]]
                        .iter()
                        .map(|e| e.1),
                );
                for o in others {
                    if o != b {
                        pairs.push((o, f));
                    }
                }
            }
        }
        drop(extra);
        drop(extra_start);
        pairs.sort_unstable();
        pairs.dedup();
        let mut ghost_start = vec![0usize; blocks + 1];
        for &(b, _) in &pairs {
            ghost_start[b as usize + 1] += 1;
        }
        for b in 0..blocks {
            ghost_start[b + 1] += ghost_start[b];
        }
        let ghosts: Vec<usize> = pairs.into_iter().map(|(_, f)| f).collect();

        let mut plan = Self {
            block_cells,
            order,
            file_id,
            node_count,
            ghost_start,
            ghosts,
            point_start: vec![0],
            points: Vec::new(),
            point_owner,
            connectivity: Vec::with_capacity(blocks),
            faces: Vec::with_capacity(blocks),
            face_connectivity: Vec::with_capacity(blocks),
            block_bounds: Vec::with_capacity(blocks),
        };

        // Points of each partition in order of first use, and its counts
        let mut seen = vec![false; np];
        for b in 0..blocks {
            let first = plan.points.len();
            let (mut k, mut nf, mut fk) = (0, 0, 0);
            let owned = plan.owned(b).len();
            // the end of the owned cells' points in `plan.points`
            let mut owned_points = first;
            let cells_b: Vec<usize> = plan.cells(b).collect();
            for (i, &f) in cells_b.iter().enumerate() {
                let c = plan.order[f];
                cell_points(mesh, cells, c, nodes);
                for &p in nodes.iter() {
                    if !seen[p] {
                        seen[p] = true;
                        plan.points.push(p);
                    }
                }
                if i + 1 == owned {
                    owned_points = plan.points.len();
                }
                k += nodes.len();
                for (face, _) in cells.cell(c) {
                    nf += 1;
                    fk += mesh.face(face).len();
                }
            }
            for &p in &plan.points[first..] {
                seen[p] = false;
            }
            plan.block_bounds.push(bounds(
                plan.points[first..owned_points]
                    .iter()
                    .map(|&p| mesh.points[p]),
            ));
            plan.point_start.push(plan.points.len());
            plan.connectivity.push(k);
            plan.faces.push(nf);
            plan.face_connectivity.push(fk);
        }
        plan
    }
}

/// The role attribute of a patch block (none when unspecified).
fn role_name(role: PatchRole) -> Option<&'static str> {
    match role {
        PatchRole::Wall => Some("wall"),
        PatchRole::Inlet => Some("inlet"),
        PatchRole::Outlet => Some("outlet"),
        PatchRole::Symmetry => Some("symmetry"),
        PatchRole::Interface => Some("interface"),
        PatchRole::Unspecified => None,
    }
}

/// The datasets of the `internal` block, in the order of their data.
struct Internal {
    number_of_points: Obj,
    number_of_cells: Obj,
    number_of_connectivity_ids: Obj,
    number_of_faces: Obj,
    number_of_face_connectivity_ids: Obj,
    number_of_polyhedron_to_face_ids: Obj,
    points: Obj,
    types: Obj,
    offsets: Obj,
    connectivity: Obj,
    face_offsets: Obj,
    face_connectivity: Obj,
    polyhedron_offsets: Obj,
    polyhedron_to_faces: Obj,
    zone: Obj,
    global_cell_ids: Obj,
    original_cell_ids: Obj,
    cell_ghosts: Obj,
    global_point_ids: Obj,
    point_ghosts: Obj,
    block_owned_cells: Obj,
    block_bounds: Obj,
    zone_names: Obj,
}

/// The datasets of a patch block.
struct PatchBlock {
    patch: usize,
    points: usize,
    connectivity: usize,
    number_of_points: Obj,
    xyz: Obj,
    /// NumberOfCells, NumberOfConnectivityIds, Offsets, Connectivity of
    /// Vertices, Lines, Polygons and Strips.
    topology: [[Obj; 4]; 4],
    global_point_ids: Obj,
    owner_cell: Obj,
}

const TOPOLOGY: [&str; 4] = ["Vertices", "Lines", "Polygons", "Strips"];
const POLYGONS: usize = 2;

fn len(n: usize) -> u64 {
    n as u64
}

/// Write the blocked VTKHDF file of a hand-off mesh. `metadata` becomes
/// string attributes of `/VTKHDF` (`koflcraft_version`, digests).
pub fn write_foehncraft_vtkhdf(
    mesh: &HandoffMesh,
    metadata: &[(&str, String)],
    out: &mut dyn Write,
) -> io::Result<()> {
    write_blocked_vtkhdf(mesh, metadata, BLOCK_CELLS, out)
}

/// [`write_foehncraft_vtkhdf`] with another number of owned cells per
/// block (the layout allows 16,384 to 1,048,576; tests use fewer).
pub fn write_blocked_vtkhdf(
    mesh: &HandoffMesh,
    metadata: &[(&str, String)],
    block_cells: usize,
    out: &mut dyn Write,
) -> io::Result<()> {
    let patches: Vec<usize> = (0..mesh.patches.len())
        .filter(|&p| !mesh.patch_faces(p).is_empty())
        .collect();
    for &p in &patches {
        let name = &mesh.patches[p].name;
        if name == "internal" || name == "Assembly" {
            return Err(io::Error::other(format!(
                "the patch name {name} is taken by the VTKHDF layout"
            )));
        }
    }
    let cells = mesh.cell_face_list();
    let plan = Plan::new(mesh, &cells, block_cells);
    let blocks = plan.blocks();
    let total_cells = len(plan.order.len() + plan.ghosts.len());
    let total_points = len(plan.points.len());
    let total_faces = len(plan.faces.iter().sum());
    let nb = len(blocks);

    let mut file = H5File::new();
    let root = file.root();
    let vtkhdf = file.group(root, "VTKHDF");
    file.attr_i64s(vtkhdf, "Version", &[2, 4]);
    file.attr_str(vtkhdf, "Type", "PartitionedDataSetCollection");
    file.attr_i64(vtkhdf, "murmelcraft_blocks", LAYOUT_VERSION);
    file.attr_str(vtkhdf, "block_order", "hilbert");
    for (name, value) in metadata {
        file.attr_str(vtkhdf, name, value);
    }
    let assembly = file.group(vtkhdf, "Assembly");
    for name in
        std::iter::once("internal").chain(patches.iter().map(|&p| mesh.patches[p].name.as_str()))
    {
        let node = file.group(assembly, name);
        file.soft_link(node, name, &format!("/VTKHDF/{name}"));
    }

    let g = file.group(vtkhdf, "internal");
    file.attr_i64(g, "Index", 0);
    file.attr_str(g, "Type", "UnstructuredGrid");
    file.attr_i64s(g, "Version", &[2, 4]);
    file.attr_i64(g, "block_cells", plan.block_cells as i64);
    let i64s = |file: &mut H5File, parent: Obj, name: &str, n: u64| {
        file.dataset(parent, name, Element::I64, &[n])
    };
    let number_of_points = i64s(&mut file, g, "NumberOfPoints", nb);
    let number_of_cells = i64s(&mut file, g, "NumberOfCells", nb);
    let number_of_connectivity_ids = i64s(&mut file, g, "NumberOfConnectivityIds", nb);
    let number_of_faces = i64s(&mut file, g, "NumberOfFaces", nb);
    let number_of_face_connectivity_ids = i64s(&mut file, g, "NumberOfFaceConnectivityIds", nb);
    let number_of_polyhedron_to_face_ids = i64s(&mut file, g, "NumberOfPolyhedronToFaceIds", nb);
    let points = file.dataset(g, "Points", Element::F64, &[total_points, 3]);
    let types = file.dataset(g, "Types", Element::U8, &[total_cells]);
    let offsets = i64s(&mut file, g, "Offsets", total_cells + nb);
    let connectivity = i64s(
        &mut file,
        g,
        "Connectivity",
        len(plan.connectivity.iter().sum()),
    );
    let face_offsets = i64s(&mut file, g, "FaceOffsets", total_faces + nb);
    let face_connectivity = i64s(
        &mut file,
        g,
        "FaceConnectivity",
        len(plan.face_connectivity.iter().sum()),
    );
    let polyhedron_offsets = i64s(&mut file, g, "PolyhedronOffsets", total_cells + nb);
    let polyhedron_to_faces = i64s(&mut file, g, "PolyhedronToFaces", total_faces);
    let cell_data = file.group(g, "CellData");
    let zone = file.dataset(cell_data, "zone", Element::U32, &[total_cells]);
    let global_cell_ids = i64s(&mut file, cell_data, "GlobalCellIds", total_cells);
    let original_cell_ids = i64s(&mut file, cell_data, "OriginalCellIds", total_cells);
    let cell_ghosts = file.dataset(cell_data, "vtkGhostType", Element::U8, &[total_cells]);
    let point_data = file.group(g, "PointData");
    let global_point_ids = i64s(&mut file, point_data, "GlobalPointIds", total_points);
    let point_ghosts = file.dataset(point_data, "vtkGhostType", Element::U8, &[total_points]);
    let field_data = file.group(g, "FieldData");
    let block_owned_cells = i64s(&mut file, field_data, "BlockOwnedCells", nb);
    let block_bounds = file.dataset(field_data, "BlockBounds", Element::F64, &[nb, 6]);
    let zone_entries: Vec<String> = mesh
        .zones
        .iter()
        .enumerate()
        .map(|(i, z)| format!("{i}={z}"))
        .collect();
    let zone_width = zone_entries
        .iter()
        .map(String::len)
        .max()
        .unwrap_or(0)
        .max(1);
    let zone_names = file.dataset(
        field_data,
        "zone_names",
        Element::Str(zone_width),
        &[len(zone_entries.len())],
    );
    let internal = Internal {
        number_of_points,
        number_of_cells,
        number_of_connectivity_ids,
        number_of_faces,
        number_of_face_connectivity_ids,
        number_of_polyhedron_to_face_ids,
        points,
        types,
        offsets,
        connectivity,
        face_offsets,
        face_connectivity,
        polyhedron_offsets,
        polyhedron_to_faces,
        zone,
        global_cell_ids,
        original_cell_ids,
        cell_ghosts,
        global_point_ids,
        point_ghosts,
        block_owned_cells,
        block_bounds,
        zone_names,
    };

    // Patch blocks: their points in order of first use
    let mut local = vec![NONE; mesh.points.len()];
    let mut patch_blocks = Vec::with_capacity(patches.len());
    for (index, &p) in patches.iter().enumerate() {
        let patch = &mesh.patches[p];
        let g = file.group(vtkhdf, &patch.name);
        file.attr_i64(g, "Index", index as i64 + 1);
        file.attr_str(g, "Type", "PolyData");
        file.attr_i64s(g, "Version", &[2, 4]);
        if let Some(role) = role_name(patch.role) {
            file.attr_str(g, "role", role);
        }
        let mut used = Vec::new();
        let mut k = 0;
        for f in mesh.patch_faces(p) {
            for &v in mesh.face(f) {
                if local[v] == NONE {
                    local[v] = 0;
                    used.push(v);
                }
            }
            k += mesh.face(f).len();
        }
        for &v in &used {
            local[v] = NONE;
        }
        let faces = len(mesh.patch_faces(p).len());
        let number_of_points = i64s(&mut file, g, "NumberOfPoints", 1);
        let xyz = file.dataset(g, "Points", Element::F64, &[len(used.len()), 3]);
        let mut topology = [[root; 4]; 4];
        for (t, name) in TOPOLOGY.iter().enumerate() {
            let tg = file.group(g, name);
            let (cells, ids) = if t == POLYGONS {
                (faces, len(k))
            } else {
                (0, 0)
            };
            topology[t] = [
                i64s(&mut file, tg, "NumberOfCells", 1),
                i64s(&mut file, tg, "NumberOfConnectivityIds", 1),
                i64s(&mut file, tg, "Offsets", cells + 1),
                i64s(&mut file, tg, "Connectivity", ids),
            ];
        }
        let pd = file.group(g, "PointData");
        let global_point_ids = i64s(&mut file, pd, "GlobalPointIds", len(used.len()));
        let cd = file.group(g, "CellData");
        let owner_cell = i64s(&mut file, cd, "OwnerCell", faces);
        file.group(g, "FieldData");
        patch_blocks.push(PatchBlock {
            patch: p,
            points: used.len(),
            connectivity: k,
            number_of_points,
            xyz,
            topology,
            global_point_ids,
            owner_cell,
        });
    }

    let mut w = file.write(out)?;
    write_internal(
        mesh,
        &cells,
        &plan,
        &internal,
        &zone_entries,
        zone_width,
        &mut local,
        &mut w,
    )?;
    for block in &patch_blocks {
        write_patch(mesh, &plan, block, &mut local, &mut w)?;
    }
    w.finish()
}

#[allow(clippy::too_many_arguments)]
fn write_internal(
    mesh: &HandoffMesh,
    cells: &CellFaceList,
    plan: &Plan,
    d: &Internal,
    zone_entries: &[String],
    zone_width: usize,
    local: &mut [u32],
    w: &mut DataWriter<'_>,
) -> io::Result<()> {
    let blocks = plan.blocks();
    let per_block = |w: &mut DataWriter<'_>, ds: Obj, value: &dyn Fn(usize) -> usize| {
        w.begin(ds)?;
        (0..blocks).try_for_each(|b| w.id(value(b)))
    };
    per_block(w, d.number_of_points, &|b| plan.partition_points(b).len())?;
    per_block(w, d.number_of_cells, &|b| plan.partition_cells(b))?;
    per_block(w, d.number_of_connectivity_ids, &|b| plan.connectivity[b])?;
    per_block(w, d.number_of_faces, &|b| plan.faces[b])?;
    per_block(w, d.number_of_face_connectivity_ids, &|b| {
        plan.face_connectivity[b]
    })?;
    per_block(w, d.number_of_polyhedron_to_face_ids, &|b| plan.faces[b])?;

    w.begin(d.points)?;
    for &p in &plan.points {
        let x = mesh.points[p];
        w.f64(x.x)?;
        w.f64(x.y)?;
        w.f64(x.z)?;
    }
    w.begin(d.types)?;
    for _ in 0..plan.order.len() + plan.ghosts.len() {
        w.u8(42)?;
    }
    w.begin(d.offsets)?;
    for b in 0..blocks {
        let mut n = 0;
        w.id(0)?;
        for f in plan.cells(b) {
            n += plan.node_count[plan.order[f]] as usize;
            w.id(n)?;
        }
    }

    // The point ids of a partition become local ids while it is written.
    let localise = |local: &mut [u32], b: usize| {
        for (i, &p) in plan.partition_points(b).iter().enumerate() {
            local[p] = i as u32;
        }
    };
    let release = |local: &mut [u32], b: usize| {
        for &p in plan.partition_points(b) {
            local[p] = NONE;
        }
    };
    let nodes = &mut Vec::new();
    w.begin(d.connectivity)?;
    for b in 0..blocks {
        localise(local, b);
        for f in plan.cells(b) {
            cell_points(mesh, cells, plan.order[f], nodes);
            for &p in nodes.iter() {
                w.i64(i64::from(local[p]))?;
            }
        }
        release(local, b);
    }
    w.begin(d.face_offsets)?;
    for b in 0..blocks {
        let mut n = 0;
        w.id(0)?;
        for f in plan.cells(b) {
            for (face, _) in cells.cell(plan.order[f]) {
                n += mesh.face(face).len();
                w.id(n)?;
            }
        }
    }
    w.begin(d.face_connectivity)?;
    for b in 0..blocks {
        localise(local, b);
        for f in plan.cells(b) {
            for (face, outward) in cells.cell(plan.order[f]) {
                let ring = mesh.face(face);
                if outward {
                    for &p in ring {
                        w.i64(i64::from(local[p]))?;
                    }
                } else {
                    for &p in ring.iter().rev() {
                        w.i64(i64::from(local[p]))?;
                    }
                }
            }
        }
        release(local, b);
    }
    w.begin(d.polyhedron_offsets)?;
    for b in 0..blocks {
        let mut n = 0;
        w.id(0)?;
        for f in plan.cells(b) {
            n += cells.len(plan.order[f]);
            w.id(n)?;
        }
    }
    w.begin(d.polyhedron_to_faces)?;
    for &faces in &plan.faces {
        for i in 0..faces {
            w.id(i)?;
        }
    }

    let each_cell =
        |w: &mut DataWriter<'_>,
         ds: Obj,
         put: &dyn Fn(&mut DataWriter<'_>, usize, usize, bool) -> io::Result<()>| {
            w.begin(ds)?;
            for b in 0..blocks {
                let owned = plan.owned(b);
                for f in plan.cells(b) {
                    put(w, f, plan.order[f], !owned.contains(&f))?;
                }
            }
            Ok::<(), io::Error>(())
        };
    each_cell(w, d.zone, &|w, _, c, _| {
        let z = u32::try_from(mesh.cell_zone[c]).map_err(|_| io::Error::other("too many zones"))?;
        w.u32(z)
    })?;
    each_cell(w, d.global_cell_ids, &|w, f, _, _| w.id(f))?;
    each_cell(w, d.original_cell_ids, &|w, _, c, _| w.id(c))?;
    each_cell(w, d.cell_ghosts, &|w, _, _, ghost| w.u8(u8::from(ghost)))?;

    w.begin(d.global_point_ids)?;
    for &p in &plan.points {
        w.id(p)?;
    }
    w.begin(d.point_ghosts)?;
    for b in 0..blocks {
        for &p in plan.partition_points(b) {
            w.u8(u8::from(plan.point_owner[p] as usize != b))?;
        }
    }

    per_block(w, d.block_owned_cells, &|b| plan.owned(b).len())?;
    w.begin(d.block_bounds)?;
    for bb in &plan.block_bounds {
        for &v in bb {
            w.f64(v)?;
        }
    }
    w.begin(d.zone_names)?;
    for z in zone_entries {
        w.str(z, zone_width)?;
    }
    Ok(())
}

fn write_patch(
    mesh: &HandoffMesh,
    plan: &Plan,
    d: &PatchBlock,
    local: &mut [u32],
    w: &mut DataWriter<'_>,
) -> io::Result<()> {
    let faces = mesh.patch_faces(d.patch);
    let mut used = Vec::with_capacity(d.points);
    for f in faces.clone() {
        for &v in mesh.face(f) {
            if local[v] == NONE {
                local[v] = used.len() as u32;
                used.push(v);
            }
        }
    }
    w.begin(d.number_of_points)?;
    w.id(used.len())?;
    w.begin(d.xyz)?;
    for &v in &used {
        let x = mesh.points[v];
        w.f64(x.x)?;
        w.f64(x.y)?;
        w.f64(x.z)?;
    }
    for (t, [cells, ids, offsets, connectivity]) in d.topology.iter().enumerate() {
        let polygons = t == POLYGONS;
        w.begin(*cells)?;
        w.id(if polygons { faces.len() } else { 0 })?;
        w.begin(*ids)?;
        w.id(if polygons { d.connectivity } else { 0 })?;
        w.begin(*offsets)?;
        w.id(0)?;
        if polygons {
            let mut n = 0;
            for f in faces.clone() {
                n += mesh.face(f).len();
                w.id(n)?;
            }
        }
        w.begin(*connectivity)?;
        if polygons {
            for f in faces.clone() {
                for &v in mesh.face(f) {
                    w.i64(i64::from(local[v]))?;
                }
            }
        }
    }
    w.begin(d.global_point_ids)?;
    for &v in &used {
        w.id(v)?;
    }
    // a boundary face is wound outward from its owner, out of the domain
    w.begin(d.owner_cell)?;
    for f in faces {
        w.id(plan.file_id[mesh.owner[f]])?;
    }
    for &v in &used {
        local[v] = NONE;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::collections::{BTreeMap, BTreeSet};

    use super::*;
    use crate::handoff::{CellFaces, patches_from_names};
    use crate::hdf5::read::{Item, read};

    /// An n × n × n grid of hexahedra (stretched per axis), two zones split
    /// at x = n/2, patches `inlet` (x−), `outlet` (x+) and `walls`.
    fn hex_grid(n: usize) -> HandoffMesh {
        let np = n + 1;
        let id = |i: usize, j: usize, k: usize| i + np * j + np * np * k;
        let mut points = Vec::new();
        for k in 0..np {
            for j in 0..np {
                for i in 0..np {
                    points.push(Vec3::new(i as f64, j as f64 * 0.5, k as f64 * 2.0));
                }
            }
        }
        let mut cells = Vec::new();
        let mut boundary = Vec::new();
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    let faces = vec![
                        vec![
                            id(i, j, k),
                            id(i, j + 1, k),
                            id(i + 1, j + 1, k),
                            id(i + 1, j, k),
                        ],
                        vec![
                            id(i, j, k + 1),
                            id(i + 1, j, k + 1),
                            id(i + 1, j + 1, k + 1),
                            id(i, j + 1, k + 1),
                        ],
                        vec![
                            id(i, j, k),
                            id(i + 1, j, k),
                            id(i + 1, j, k + 1),
                            id(i, j, k + 1),
                        ],
                        vec![
                            id(i, j + 1, k),
                            id(i, j + 1, k + 1),
                            id(i + 1, j + 1, k + 1),
                            id(i + 1, j + 1, k),
                        ],
                        vec![
                            id(i, j, k),
                            id(i, j, k + 1),
                            id(i, j + 1, k + 1),
                            id(i, j + 1, k),
                        ],
                        vec![
                            id(i + 1, j, k),
                            id(i + 1, j + 1, k),
                            id(i + 1, j + 1, k + 1),
                            id(i + 1, j, k + 1),
                        ],
                    ];
                    let on = [k == 0, k == n - 1, j == 0, j == n - 1, i == 0, i == n - 1];
                    for (s, ring) in faces.iter().enumerate() {
                        if on[s] {
                            let patch = [2, 2, 2, 2, 0, 1][s];
                            boundary.push((ring.clone(), patch));
                        }
                    }
                    cells.push(CellFaces {
                        faces,
                        zone: usize::from(2 * i >= n),
                    });
                }
            }
        }
        HandoffMesh::from_cells(
            points,
            cells,
            vec!["fluid".into(), "rotor".into()],
            &boundary,
            patches_from_names(&["inlet", "outlet", "walls"]),
        )
        .unwrap()
    }

    /// Invariant: the index is a Hilbert curve: the 64 points of the
    /// 4 × 4 × 4 corner of the grid get the indices 0 to 63, and
    /// consecutive indices are neighbouring points.
    /// Breaks if: a step of Skilling's transpose or the bit interleaving
    /// changes. No other test: the files stay valid with any order, only
    /// regions read more blocks.
    #[test]
    fn hilbert_index_walks_the_grid_one_step_at_a_time() {
        let mut by_index = BTreeMap::new();
        for x in 0..4u32 {
            for y in 0..4u32 {
                for z in 0..4u32 {
                    by_index.insert(hilbert_index([x, y, z]), [x, y, z]);
                }
            }
        }
        assert_eq!(
            by_index.keys().copied().collect::<Vec<_>>(),
            (0..64).collect::<Vec<_>>()
        );
        let walk: Vec<[u32; 3]> = by_index.into_values().collect();
        assert_eq!(walk[0], [0, 0, 0]);
        for pair in walk.windows(2) {
            let step: u32 = (0..3).map(|k| pair[0][k].abs_diff(pair[1][k])).sum();
            assert_eq!(step, 1, "{:?} to {:?}", pair[0], pair[1]);
        }
    }

    /// A partition cell read back: file cell id, KoflCraft cell, ghost,
    /// zone, nodes, faces.
    type Cell = (usize, usize, u8, u32, Vec<usize>, Vec<Vec<usize>>);

    /// One partition read back, in mesh point ids.
    struct Part {
        owned: usize,
        points: Vec<usize>,
        point_ghost: Vec<u8>,
        cells: Vec<Cell>,
    }

    fn parts(items: &BTreeMap<String, Item>) -> Vec<Part> {
        let get = |p: &str| &items[&format!("/VTKHDF/internal/{p}")];
        let ints =
            |p: &str| -> Vec<usize> { get(p).i64s().into_iter().map(|v| v as usize).collect() };
        let (np, nc, nk) = (
            ints("NumberOfPoints"),
            ints("NumberOfCells"),
            ints("NumberOfConnectivityIds"),
        );
        let (nf, nfk) = (ints("NumberOfFaces"), ints("NumberOfFaceConnectivityIds"));
        assert_eq!(nf, ints("NumberOfPolyhedronToFaceIds"));
        let owned = ints("FieldData/BlockOwnedCells");
        let (offsets, conn) = (ints("Offsets"), ints("Connectivity"));
        let (face_offsets, face_conn) = (ints("FaceOffsets"), ints("FaceConnectivity"));
        let (poly_offsets, poly_faces) = (ints("PolyhedronOffsets"), ints("PolyhedronToFaces"));
        let (gcell, ocell, gpoint) = (
            ints("CellData/GlobalCellIds"),
            ints("CellData/OriginalCellIds"),
            ints("PointData/GlobalPointIds"),
        );
        let (cghost, pghost) = (
            &get("CellData/vtkGhostType").data,
            &get("PointData/vtkGhostType").data,
        );
        let zone = get("CellData/zone").u32s();
        assert!(get("Types").data.iter().all(|&t| t == 42));
        let (mut p0, mut c0, mut k0, mut f0, mut fk0) = (0, 0, 0, 0, 0);
        let mut out = Vec::new();
        for b in 0..np.len() {
            let points: Vec<usize> = gpoint[p0..p0 + np[b]].to_vec();
            let mut cells = Vec::new();
            for j in 0..nc[b] {
                let c = c0 + j;
                let o = &offsets[c0 + b..];
                let nodes = conn[k0 + o[j]..k0 + o[j + 1]]
                    .iter()
                    .map(|&l| points[l])
                    .collect();
                let po = &poly_offsets[c0 + b..];
                let faces = poly_faces[f0 + po[j]..f0 + po[j + 1]]
                    .iter()
                    .map(|&lf| {
                        let fo = &face_offsets[f0 + b..];
                        face_conn[fk0 + fo[lf]..fk0 + fo[lf + 1]]
                            .iter()
                            .map(|&l| points[l])
                            .collect()
                    })
                    .collect();
                cells.push((gcell[c], ocell[c], cghost[c], zone[c], nodes, faces));
            }
            out.push(Part {
                owned: owned[b],
                point_ghost: pghost[p0..p0 + np[b]].to_vec(),
                points,
                cells,
            });
            p0 += np[b];
            c0 += nc[b];
            k0 += nk[b];
            f0 += nf[b];
            fk0 += nfk[b];
        }
        assert_eq!(
            (p0, c0, k0, f0, fk0),
            (
                gpoint.len(),
                gcell.len(),
                conn.len(),
                poly_faces.len(),
                face_conn.len()
            )
        );
        out
    }

    /// Invariant: each cell is owned by exactly one block, in file order;
    /// every partition holds its owned cells with their faces wound
    /// outward, then exactly its point neighbours from other blocks as
    /// ghosts, and marks the points a lower block owns; patches keep
    /// their faces and owners; the HDF5 checksums hold.
    /// Breaks if: the ghost layer, the local point ids, the per-partition
    /// offsets, a winding or the owner mapping change. No other test: a
    /// reader would load such a file without complaint and show wrong
    /// cells at block borders.
    #[test]
    fn every_cell_is_owned_once_and_each_block_holds_its_point_neighbours() {
        let mesh = hex_grid(5);
        let block = 16;
        let mut out = Vec::new();
        write_blocked_vtkhdf(
            &mesh,
            &[("koflcraft_version", "9.9".into())],
            block,
            &mut out,
        )
        .unwrap();
        let items = read(&out);
        let root = &items["/VTKHDF"];
        assert_eq!(
            root.attrs["Version"],
            [2i64.to_le_bytes(), 4i64.to_le_bytes()].concat()
        );
        assert_eq!(root.attrs["Type"], b"PartitionedDataSetCollection");
        assert_eq!(root.attrs["koflcraft_version"], b"9.9");
        assert_eq!(root.attrs["murmelcraft_blocks"], 1i64.to_le_bytes());
        assert_eq!(
            root.links,
            ["Assembly", "internal", "inlet", "outlet", "walls"]
        );
        assert_eq!(
            items["/VTKHDF/Assembly/walls/walls"].soft.as_deref(),
            Some("/VTKHDF/walls")
        );
        assert_eq!(
            items["/VTKHDF/internal"].attrs["block_cells"],
            16i64.to_le_bytes()
        );

        let n = mesh.cell_count();
        let list = mesh.cell_face_list();
        let parts = parts(&items);
        assert_eq!(parts.len(), n.div_ceil(block));
        // the owned cells, block by block, are the file cell ids in order
        let mut orig_of = Vec::new();
        let mut block_of = Vec::new();
        for (b, part) in parts.iter().enumerate() {
            for (j, cell) in part.cells.iter().enumerate() {
                assert_eq!(cell.2, u8::from(j >= part.owned), "owned cells come first");
                if j < part.owned {
                    assert_eq!(cell.0, orig_of.len());
                    orig_of.push(cell.1);
                    block_of.push(b);
                }
            }
        }
        assert_eq!(orig_of.iter().copied().collect::<BTreeSet<_>>().len(), n);
        let file_of: BTreeMap<usize, usize> =
            orig_of.iter().enumerate().map(|(f, &c)| (c, f)).collect();
        let nodes_of = |c: usize| -> BTreeSet<usize> {
            list.cell(c)
                .flat_map(|(f, _)| mesh.face(f).to_vec())
                .collect()
        };
        let bb = items["/VTKHDF/internal/FieldData/BlockBounds"].f64s();
        for (b, part) in parts.iter().enumerate() {
            for (f, c, _, zone, nodes, faces) in &part.cells {
                assert_eq!(file_of[c], *f, "a ghost carries its original's id");
                assert_eq!(*zone as usize, mesh.cell_zone[*c]);
                assert_eq!(nodes.iter().copied().collect::<BTreeSet<_>>(), nodes_of(*c));
                let expected: Vec<Vec<usize>> = list
                    .cell(*c)
                    .map(|(face, outward)| {
                        let mut ring = mesh.face(face).to_vec();
                        if !outward {
                            ring.reverse();
                        }
                        ring
                    })
                    .collect();
                assert_eq!(faces, &expected, "faces wound outward from the cell");
                let mut first_use: Vec<usize> = Vec::new();
                for p in expected.iter().flatten() {
                    if !first_use.contains(p) {
                        first_use.push(*p);
                    }
                }
                assert_eq!(nodes, &first_use, "points in the order the faces use them");
            }
            // the ghosts are exactly the point neighbours of other blocks
            let used: BTreeSet<usize> = part.cells[..part.owned]
                .iter()
                .flat_map(|c| nodes_of(c.1))
                .collect();
            let expected: Vec<usize> = (0..n)
                .filter(|&f| {
                    block_of[f] != b && nodes_of(orig_of[f]).iter().any(|p| used.contains(p))
                })
                .collect();
            let ghosts: Vec<usize> = part.cells[part.owned..].iter().map(|c| c.0).collect();
            assert_eq!(ghosts, expected, "block {b}");
            for (p, &ghost) in part.points.iter().zip(&part.point_ghost) {
                let owner = (0..n)
                    .filter(|&f| nodes_of(orig_of[f]).contains(p))
                    .map(|f| block_of[f])
                    .min();
                assert_eq!(ghost, u8::from(owner != Some(b)));
                assert_eq!(part.points.iter().filter(|&&q| q == *p).count(), 1);
            }
            let tight = bounds(used.iter().map(|&p| mesh.points[p]));
            assert_eq!(&bb[6 * b..6 * b + 6], &tight);
        }

        for (k, patch) in mesh.patches.iter().enumerate() {
            let at = |p: &str| &items[&format!("/VTKHDF/{}/{p}", patch.name)];
            assert_eq!(
                items[&format!("/VTKHDF/{}", patch.name)].attrs["Index"],
                (k as i64 + 1).to_le_bytes()
            );
            let gp = at("PointData/GlobalPointIds").i64s();
            let offsets = at("Polygons/Offsets").i64s();
            let conn = at("Polygons/Connectivity").i64s();
            let owner = at("CellData/OwnerCell").i64s();
            let faces: Vec<usize> = mesh.patch_faces(k).collect();
            assert_eq!(owner.len(), faces.len());
            for (i, &f) in faces.iter().enumerate() {
                let ring: Vec<usize> = conn[offsets[i] as usize..offsets[i + 1] as usize]
                    .iter()
                    .map(|&l| gp[l as usize] as usize)
                    .collect();
                assert_eq!(ring, mesh.face(f), "wound out of the domain");
                assert_eq!(orig_of[owner[i] as usize], mesh.owner[f]);
            }
            assert!(at("Lines/Connectivity").data.is_empty());
            assert_eq!(at("Lines/Offsets").i64s(), [0]);
        }
    }
}
