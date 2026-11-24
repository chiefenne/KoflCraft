"""Build the synthetic cooling-jacket fixture (fluid volume) as a binary STL.

A clean-room stand-in for an industrial motor cooling jacket: a thin annular
fluid shell whose inner side carries six ribs (rings with a gap, so the
coolant passes between the channels), two blocks on the outer side and two
bent pipes, each ending in a flat circular cap (the ports). Lengths in mm.

    python3 -m pip install manifold3d numpy
    python3 make_jacket.py cooling-jacket.stl
    python3 make_jacket.py cooling-jacket-thin.stl thin

The `thin` variant adds the two thin features of the industrial part: rib 3
reaches to 0.5 mm below the outer wall (a thin fluid gap over its top) and
rib 5 is 0.5 mm wide (a thin solid wall between two channels).
"""
import struct
import sys

import numpy as np
from manifold3d import Manifold, Mesh, set_circular_segments

set_circular_segments(96)

R_IN, R_OUT, HEIGHT = 60.0, 66.0, 60.0
RIB_DEPTH, RIB_WIDTH, RIBS, GAP_DEG = 4.0, 2.0, 6, 24.0
PIPE_R, BEND_R = 4.0, 12.0


def ring(r_in, r_out, z0, z1):
    return Manifold.cylinder(z1 - z0, r_out).translate([0, 0, z0]) - Manifold.cylinder(
        z1 - z0 + 2, r_in
    ).translate([0, 0, z0 - 1])


THIN_GAP, THIN_WALL = 0.5, 0.5


def jacket(thin=False):
    fluid = ring(R_IN, R_OUT, 0.0, HEIGHT)
    pitch = HEIGHT / (RIBS + 1)
    gap = (
        Manifold.cube([R_OUT + 10, 2 * (R_OUT + 10) * np.tan(np.radians(GAP_DEG / 2)), HEIGHT + 4])
        .translate([0, -(R_OUT + 10) * np.tan(np.radians(GAP_DEG / 2)), -2])
    )
    for k in range(1, RIBS + 1):
        z = k * pitch
        depth = R_OUT - R_IN - THIN_GAP if thin and k == 3 else RIB_DEPTH
        width = THIN_WALL if thin and k == 5 else RIB_WIDTH
        rib = ring(R_IN - 1, R_IN + depth, z - width / 2, z + width / 2) - gap
        fluid = fluid - rib
    return fluid


def pipe(z, down):
    """A block on the outer wall and a swept tube: a radial leg, a 45° bend
    of centreline radius BEND_R tangent to both legs, and an end leg turning
    down or up, closed by flat caps (the far one is the port)."""
    block = Manifold.cube([10.0, 12.0, 12.0]).translate([R_OUT - 2.0, -6.0, z - 6.0])
    k = -1.0 if down else 1.0
    start = np.array([R_OUT + 6.0, 0.0, z])
    knee = start + np.array([30.0, 0.0, 0.0])
    centre = knee + np.array([0.0, 0.0, k * BEND_R])
    frames = []  # (point, tangent) along the centreline
    for t in np.linspace(0.0, 1.0, 16):
        frames.append((start + t * (knee - start), np.array([1.0, 0.0, 0.0])))
    for a in np.radians(np.linspace(0.0, 45.0, 25))[1:]:
        frames.append((centre + BEND_R * np.array([np.sin(a), 0.0, -k * np.cos(a)]),
                       np.array([np.cos(a), 0.0, k * np.sin(a)])))
    end, d = frames[-1]
    for t in np.linspace(0.0, 25.0, 13)[1:]:
        frames.append((end + t * d, d))
    return block + tube(frames)


def tube(frames, n=96):
    """A closed tube of radius PIPE_R through the frames, flat end caps."""
    side = np.array([0.0, 1.0, 0.0])  # the bend stays in the xz plane
    verts = []
    for p, t in frames:
        u = np.cross(t, side)
        u /= np.linalg.norm(u)
        for j in range(n):
            a = 2 * np.pi * j / n
            verts.append(p + PIPE_R * (np.cos(a) * u + np.sin(a) * side))
    rings = len(frames)
    tris = []
    for i in range(rings - 1):
        for j in range(n):
            a, b = i * n + j, i * n + (j + 1) % n
            c, d = a + n, b + n
            tris += [[a, d, b], [a, c, d]]
    c0, c1 = len(verts), len(verts) + 1
    verts += [frames[0][0], frames[-1][0]]
    last = (rings - 1) * n
    for j in range(n):
        tris.append([c0, j, (j + 1) % n])
        tris.append([c1, last + (j + 1) % n, last + j])
    mesh = Mesh(vert_properties=np.array(verts, np.float32), tri_verts=np.array(tris, np.uint32))
    body = Manifold(mesh)
    if body.volume() <= 0:
        raise RuntimeError(f"tube: {body.status()} volume {body.volume()}")
    return body


def main(path, thin=False):
    body = jacket(thin) + pipe(15.0, True) + pipe(45.0, False)
    if body.genus() != 1:
        raise RuntimeError(f"genus {body.genus()}")
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
    main(sys.argv[1] if len(sys.argv) > 1 else "cooling-jacket.stl", "thin" in sys.argv[2:])
