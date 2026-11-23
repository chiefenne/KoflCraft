"""Build the synthetic cooling-jacket fixture (fluid volume) as a binary STL.

A clean-room stand-in for an industrial motor cooling jacket: a thin annular
fluid shell whose inner side carries six ribs (rings with a gap, so the
coolant passes between the channels), two blocks on the outer side and two
bent pipes, each ending in a flat circular cap (the ports). Lengths in mm.

    python3 -m pip install manifold3d numpy
    python3 make_jacket.py cooling-jacket.stl
"""
import struct
import sys

import numpy as np
from manifold3d import CrossSection, Manifold, set_circular_segments

set_circular_segments(96)

R_IN, R_OUT, HEIGHT = 60.0, 66.0, 60.0
RIB_DEPTH, RIB_WIDTH, RIBS, GAP_DEG = 4.0, 2.0, 6, 24.0
PIPE_R = 4.0


def ring(r_in, r_out, z0, z1):
    return Manifold.cylinder(z1 - z0, r_out).translate([0, 0, z0]) - Manifold.cylinder(
        z1 - z0 + 2, r_in
    ).translate([0, 0, z0 - 1])


def jacket():
    fluid = ring(R_IN, R_OUT, 0.0, HEIGHT)
    pitch = HEIGHT / (RIBS + 1)
    gap = (
        Manifold.cube([R_OUT + 10, 2 * (R_OUT + 10) * np.tan(np.radians(GAP_DEG / 2)), HEIGHT + 4])
        .translate([0, -(R_OUT + 10) * np.tan(np.radians(GAP_DEG / 2)), -2])
    )
    for k in range(1, RIBS + 1):
        z = k * pitch
        rib = ring(R_IN - 1, R_IN + RIB_DEPTH, z - RIB_WIDTH / 2, z + RIB_WIDTH / 2) - gap
        fluid = fluid - rib
    return fluid


def pipe(z, out_dir):
    """A block on the outer wall, a radial pipe and a bent end leg."""
    block = Manifold.cube([10.0, 12.0, 12.0]).translate([R_OUT - 2.0, -6.0, z - 6.0])
    start = np.array([R_OUT + 6.0, 0.0, z])
    knee = start + np.array([30.0, 0.0, 0.0])
    end = knee + 25.0 * np.asarray(out_dir) / np.linalg.norm(out_dir)
    leg1 = segment(start, knee)
    leg2 = segment(knee, end)
    joint = Manifold.sphere(PIPE_R, 96).translate(list(knee))
    return block + leg1 + leg2 + joint


def segment(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = b - a
    length = np.linalg.norm(d)
    cyl = Manifold.cylinder(length, PIPE_R)
    z = d / length
    # rotate +z onto d
    v = np.cross([0, 0, 1], z)
    s, c = np.linalg.norm(v), z[2]
    if s < 1e-12:
        m = np.eye(3) if c > 0 else np.diag([1, -1, -1])
    else:
        k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        m = np.eye(3) + k + k @ k * ((1 - c) / s**2)
    t = np.hstack([m, a.reshape(3, 1)])
    return cyl.transform(t.tolist())


def main(path):
    body = jacket() + pipe(15.0, [1.0, 0.0, -1.0]) + pipe(45.0, [1.0, 0.0, 1.0])
    mesh = body.to_mesh()
    v = np.asarray(mesh.vert_properties)[:, :3]
    f = np.asarray(mesh.tri_verts)
    with open(path, "wb") as out:
        out.write(b"KoflCraft synthetic cooling jacket".ljust(80, b" "))
        out.write(struct.pack("<I", len(f)))
        for t in f:
            p = v[t]
            n = np.cross(p[1] - p[0], p[2] - p[0])
            n = n / np.linalg.norm(n)
            out.write(struct.pack("<12fH", *n, *p.ravel(), 0))
    print(f"{len(f)} facets, genus {body.genus()}, volume {body.volume():.3f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "cooling-jacket.stl")
