//! Export of block-structured grids as unstructured hexahedra (BS-14): a
//! VTU with the block and the i, j, k index of every cell, MSH 4.1 with
//! hexahedra and boundary quadrangles, and the FoehnCraft hand-off mesh.

use std::io::{self, Write};

use mesh_algorithms::block::BlockGrid;
use mesh_algorithms::block::grid::HEX_FACES;

use crate::handoff::{CellFaces, HandoffError, HandoffMesh, Patch};
use crate::msh::{GmshVersion, Model};
use crate::{end_array, end_vtu, start_array, vtu_points, write_array, write_ids};

/// The domains of the grid in order of first appearance, and each block's
/// domain index.
pub fn domains(grid: &BlockGrid) -> (Vec<String>, Vec<usize>) {
    let mut names: Vec<String> = Vec::new();
    let of = grid
        .blocks
        .iter()
        .map(|b| match names.iter().position(|n| *n == b.domain) {
            Some(i) => i,
            None => {
                names.push(b.domain.clone());
                names.len() - 1
            }
        })
        .collect();
    (names, of)
}

/// The volume as VTK hexahedra (type 12) with cell data `BlockId`, `I`,
/// `J`, `K`.
pub fn write_volume(grid: &BlockGrid, out: &mut dyn Write) -> io::Result<()> {
    vtu_points(&grid.points, grid.hexes.len(), out)?;
    start_array(out, "Int64", "connectivity")?;
    for h in &grid.hexes {
        write!(out, "          ")?;
        write_ids(out, h, false)?;
        writeln!(out)?;
    }
    end_array(out)?;
    let offsets: Vec<usize> = (1..=grid.hexes.len()).map(|c| 8 * c).collect();
    write_array(out, "offsets", &offsets)?;
    start_array(out, "UInt8", "types")?;
    for _ in &grid.hexes {
        writeln!(out, "          12")?;
    }
    end_array(out)?;
    writeln!(out, "      </Cells>\n      <CellData Scalars=\"BlockId\">")?;
    let n = grid.hexes.len();
    let (mut block, mut i, mut j, mut k) = (
        Vec::with_capacity(n),
        Vec::with_capacity(n),
        Vec::with_capacity(n),
        Vec::with_capacity(n),
    );
    for r in &grid.blocks {
        let [ni, nj, nk] = r.cells;
        for kk in 0..nk {
            for jj in 0..nj {
                for ii in 0..ni {
                    block.push(r.block);
                    i.push(ii);
                    j.push(jj);
                    k.push(kk);
                }
            }
        }
    }
    write_array(out, "BlockId", &block)?;
    write_array(out, "I", &i)?;
    write_array(out, "J", &j)?;
    write_array(out, "K", &k)?;
    writeln!(out, "      </CellData>")?;
    end_vtu(out)
}

/// MSH (4.1 by default): one physical surface per patch (tags 1…P), one
/// physical volume per domain (tags P+1…), one volume entity per block.
pub fn write_gmsh(
    grid: &BlockGrid,
    patches: &[Patch],
    version: GmshVersion,
    out: &mut dyn Write,
) -> io::Result<()> {
    let (names, domain_of) = domains(grid);
    let mut model = Model::new(grid.points.clone());
    for (p, patch) in patches.iter().enumerate() {
        model.group(2, p + 1, &patch.name);
    }
    for (d, name) in names.iter().enumerate() {
        model.group(3, patches.len() + 1 + d, name);
    }
    for &(q, p, _) in &grid.boundary {
        model.element(3, p + 1, p + 1, &q)?;
    }
    for (b, r) in grid.blocks.iter().enumerate() {
        let count = r.cells.iter().product::<usize>();
        for h in &grid.hexes[r.first_cell..r.first_cell + count] {
            model.element(5, patches.len() + 1 + domain_of[b], b + 1, h)?;
        }
    }
    model.write(version, out)
}

/// The FoehnCraft hand-off mesh: patches as given, one zone per domain.
pub fn handoff(grid: &BlockGrid, patches: Vec<Patch>) -> Result<HandoffMesh, HandoffError> {
    let (names, domain_of) = domains(grid);
    let mut cells = Vec::with_capacity(grid.hexes.len());
    for (b, r) in grid.blocks.iter().enumerate() {
        let count = r.cells.iter().product::<usize>();
        for h in &grid.hexes[r.first_cell..r.first_cell + count] {
            cells.push(CellFaces {
                faces: HEX_FACES.iter().map(|f| f.map(|s| h[s]).to_vec()).collect(),
                zone: domain_of[b],
            });
        }
    }
    let boundary: Vec<(Vec<usize>, usize)> = grid
        .boundary
        .iter()
        .map(|&(q, p, _)| (q.to_vec(), p))
        .collect();
    HandoffMesh::from_cells(grid.points.clone(), cells, names, &boundary, patches)
}
