"""Build the box-furnace fixture (fluid volume) as a binary STL and its
koflcraft.patch-map.

The geometry of FoehnCraft's radiation gallery furnace (example 259): a
chamber 4 m long (x), 2 m deep (y) and 2.5 m high (z) with horizontal
tubes of 0.16 m diameter through it under the roof (z = 2.1 m, from
x = 0.9 m in steps of 0.4 m), a burner disc of 0.5 m diameter in the front
wall (x = 0) and a rectangular opening in the back wall (x = 4 m) under the
roof. The burner and the opening are flush with their walls. Patches:
`walls` (floor, roof, side walls, rest of the back wall), `front` (rest of
the front wall), `tubes`, `burner` (role inlet) and `outlet` (role outlet).

    python3 -m pip install gmsh numpy
    python3 make_furnace.py furnace.stl 1      # the fixture: one tube
    python3 make_furnace.py furnace-7.stl 7    # the gallery furnace

The patch map goes beside the STL (`furnace.patches.json`).
"""
import hashlib
import json
import struct
import sys

import gmsh
import numpy as np

LENGTH, DEPTH, HEIGHT = 4.0, 2.0, 2.5
TUBE_R, TUBE_Z = 0.08, 2.1
BURNER = (DEPTH / 2, 0.9, 0.25)  # y, z of the centre, radius
OUTLET = (0.6, 1.4, 1.75, 2.35)  # y0, y1, z0, z1 in the back wall
RIM = 0.02  # STL spacing on the burner and opening rims


def main():
    out = sys.argv[1]
    tubes_n = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    tube_x = [0.9 + 0.4 * i for i in range(tubes_n)]
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    occ = gmsh.model.occ
    box = occ.addBox(0, 0, 0, LENGTH, DEPTH, HEIGHT)
    tubes = [occ.addCylinder(x, -0.1, TUBE_Z, 0, DEPTH + 0.2, 0, TUBE_R) for x in tube_x]
    chamber = occ.cut([(3, box)], [(3, t) for t in tubes])[0] if tubes else [(3, box)]
    yb, zb, rb = BURNER
    disc = occ.addDisk(0, yb, zb, rb, rb, zAxis=[1, 0, 0], xAxis=[0, 1, 0])
    y0, y1, z0, z1 = OUTLET
    pts = [occ.addPoint(LENGTH, y, z) for y, z in [(y0, z0), (y1, z0), (y1, z1), (y0, z1)]]
    lines = [occ.addLine(pts[i], pts[(i + 1) % 4]) for i in range(4)]
    opening = occ.addPlaneSurface([occ.addCurveLoop(lines)])
    occ.fragment(chamber, [(2, disc), (2, opening)])
    occ.synchronize()
    vol = gmsh.model.getEntities(3)
    eps = 1e-6
    names = ["walls", "front", "tubes", "burner", "outlet"]
    groups = {n: [] for n in names}
    for _, tag in gmsh.model.getBoundary(vol, oriented=False):
        x0, ya, za, x1, yb_, zb_ = gmsh.model.getBoundingBox(2, tag)
        if x1 < eps:
            r = max(ya - (yb - rb), (yb + rb) - yb_, za - (zb - rb), (zb + rb) - zb_)
            burner = abs(r) < 1e-3 and yb_ - ya < 2 * rb + 1e-3
            groups["burner" if burner else "front"].append(tag)
        elif x0 > LENGTH - eps:
            inside = ya > y0 - eps and yb_ < y1 + eps and za > z0 - eps and zb_ < z1 + eps
            groups["outlet" if inside else "walls"].append(tag)
        elif zb_ - za < 2 * TUBE_R + 1e-3 and za > TUBE_Z - TUBE_R - 1e-3 and x1 - x0 < 2 * TUBE_R + 1e-3:
            groups["tubes"].append(tag)
        else:
            groups["walls"].append(tag)
    for name in ("burner", "outlet"):
        ends = gmsh.model.getBoundary([(2, t) for t in groups[name]], recursive=True)
        gmsh.model.mesh.setSize(ends, RIM)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.25)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
    gmsh.model.mesh.generate(2)
    ntag, xyz, _ = gmsh.model.mesh.getNodes()
    nodes = dict(zip(ntag, xyz.reshape(-1, 3)))
    facets, patch = [], []
    for gi, name in enumerate(names):
        for tag in groups[name]:
            _, _, en = gmsh.model.mesh.getElements(2, tag)
            for a, b, c in np.array(en[0]).reshape(-1, 3):
                facets.append([nodes[a], nodes[b], nodes[c]])
                patch.append(gi)
    gmsh.finalize()
    f = np.array(facets)
    n = np.cross(f[:, 1] - f[:, 0], f[:, 2] - f[:, 0])
    c = f.mean(axis=1)
    outward = np.zeros_like(c)
    for i, (cc, pi) in enumerate(zip(c, patch)):
        if names[pi] == "tubes":
            xt = min(tube_x, key=lambda x: abs(x - cc[0]))
            outward[i] = [xt - cc[0], 0, TUBE_Z - cc[2]]  # into the tube
        else:
            d = np.concatenate([cc, np.array([LENGTH, DEPTH, HEIGHT]) - cc])
            k = int(np.argmin(d))
            outward[i][k % 3] = -1 if k < 3 else 1
    flip = (n * outward).sum(1) < 0
    f[flip] = f[flip][:, [0, 2, 1]]
    n = np.cross(f[:, 1] - f[:, 0], f[:, 2] - f[:, 0])
    data = bytearray(b"box furnace".ljust(80, b" ")) + struct.pack("<I", len(f))
    for tri, nn in zip(f, n):
        data += struct.pack("<12fH", *(nn / np.linalg.norm(nn)), *tri.flatten(), 0)
    with open(out, "wb") as fh:
        fh.write(data)
    roles = {"burner": "inlet", "outlet": "outlet"}
    used = [k for i, k in enumerate(names) if i in patch]
    patch_map = {
        "format": "koflcraft.patch-map",
        "format_version": "1.0",
        "source_sha256": hashlib.sha256(data).hexdigest(),
        "facet_count": len(f),
        "patches": [{"name": k, "role": roles.get(k, "wall")} for k in used],
        "facet_patch": [used.index(names[p]) for p in patch],
    }
    with open(out.rsplit(".", 1)[0] + ".patches.json", "w") as fh:
        json.dump(patch_map, fh)
    print(f"{out}: {len(f)} facets", {k: patch.count(i) for i, k in enumerate(names)})


if __name__ == "__main__":
    main()
