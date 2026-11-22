use std::error::Error;
use std::io::{self, BufWriter, Write};

use koflcraft::generate_mesh;
use koflcraft::geometry::{AxisAlignedBox, Vec3};
use koflcraft::mesh_export::write_vtu;

fn main() -> Result<(), Box<dyn Error>> {
    let domain = AxisAlignedBox::new(Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let sites = [
        Vec3::new(0.15, 0.22, 0.33),
        Vec3::new(0.72, 0.28, 0.11),
        Vec3::new(0.60, 0.72, 0.54),
        Vec3::new(0.30, 0.50, 0.82),
    ];
    let mesh = generate_mesh(domain, &sites)?;
    let mut output = BufWriter::new(io::stdout().lock());
    write_vtu(domain, mesh.mesh(), &mut output)?;
    output.flush()?;
    Ok(())
}
