"""Write the shell-side fluid of the heat exchanger as a closed STL.

    python3 make_shell_side.py SHELL_TUBE_HX.step shell_side.stl [facet size in mm]

KoflCraft has no STEP reader yet (STP-01 to STP-13), so this script
prepares the fluid with Gmsh (pip install gmsh): the fluid is the empty
space of the model, so it closes the housing's inner surface with a disc
over each nozzle opening, fills it, cuts away the tubes, the tube sheets
and the baffles, keeps the piece that touches the two shell-side nozzles
(air) and meshes its surface. The STL's triangles are made consistently
outward by edge propagation, so KoflCraft reads it as one closed shell.
"""

import collections
import math
import sys

import gmsh

step, out = sys.argv[1], sys.argv[2]
size = float(sys.argv[3]) if len(sys.argv) > 3 else 1.5
occ = gmsh.model.occ
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
occ.importShapes(step)
occ.synchronize()
solids = [t for _, t in gmsh.model.getEntities(3)]
# the housing is the solid with the widest bounding box along the axis
housing = max(solids, key=lambda t: gmsh.model.getBoundingBox(3, t)[4] - gmsh.model.getBoundingBox(3, t)[1])


def edges(face):
    return [t for _, t in gmsh.model.getBoundary([(2, face)], oriented=False)]


faces = [t for _, t in gmsh.model.getBoundary([(3, housing)], oriented=False)]
# the annular end faces of the four nozzles: planar, at |x| = 100 or z = 100
ends = []
for f in faces:
    bb = gmsh.model.getBoundingBox(2, f)
    if (abs(bb[0] - bb[3]) < 1e-6 and abs(abs(bb[0]) - 100) < 1e-6) or (
        abs(bb[2] - bb[5]) < 1e-6 and abs(bb[2] - 100) < 1e-6
    ):
        ends.append(f)
rest = [f for f in faces if f not in ends]
owners = collections.defaultdict(list)
for f in rest:
    for e in edges(f):
        owners[e].append(f)
adj = collections.defaultdict(set)
for fs in owners.values():
    for a in fs:
        adj[a].update(b for b in fs if b != a)
groups, seen = [], set()
for f in rest:
    if f in seen:
        continue
    stack, group = [f], []
    seen.add(f)
    while stack:
        x = stack.pop()
        group.append(x)
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    groups.append(group)


def is_bore(f):  # the cylinder of radius 45 inside the shell
    bb = gmsh.model.getBoundingBox(2, f)
    return abs(bb[0] + 45) < 1e-3 and abs(bb[3] - 45) < 1e-3 and bb[4] - bb[1] > 100


inner = next(g for g in groups if any(is_bore(f) for f in g))
free = [e for e, fs in owners.items() if len([f for f in fs if f in inner]) == 1]
caps = [occ.addPlaneSurface([occ.addCurveLoop([e])]) for e in free]
cavity = occ.addVolume([occ.addSurfaceLoop(inner + caps)])
occ.synchronize()
occ.cut([(3, cavity)], [(3, s) for s in solids if s != housing], removeObject=True, removeTool=False)
occ.synchronize()
# the shell side: the piece reaching both nozzles along x, between the tube sheets
air = None
for _, t in gmsh.model.getEntities(3):
    bb = gmsh.model.getBoundingBox(3, t)
    if bb[0] < -99.9 and bb[3] > 99.9 and abs(bb[1] - 6) < 1e-3:
        air = t
volume = occ.getMass(3, air)
boundary = gmsh.model.getBoundary([(3, air)], oriented=True)
gmsh.model.removeEntities([(3, t) for _, t in gmsh.model.getEntities(3) if t != air], recursive=True)
keep = {abs(t) for _, t in boundary}
gmsh.model.removeEntities([(2, t) for _, t in gmsh.model.getEntities(2) if t not in keep])
gmsh.option.setNumber("Mesh.MeshSizeMax", size)
gmsh.option.setNumber("Mesh.MeshSizeMin", size / 2)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
gmsh.model.mesh.generate(2)
ids, coords, _ = gmsh.model.mesh.getNodes()
at = {int(i): tuple(coords[3 * k : 3 * k + 3]) for k, i in enumerate(ids)}
tris = []
for _, t in boundary:
    types, _, nodes = gmsh.model.mesh.getElements(2, abs(t))
    for ty, nd in zip(types, nodes):
        if ty == 2:
            tris += [[int(x) for x in nd[i : i + 3]] for i in range(0, len(nd), 3)]
gmsh.finalize()

# consistent winding: neighbours must run their shared edge in opposite directions
edge = collections.defaultdict(list)
for i, t in enumerate(tris):
    for k in range(3):
        a, b = t[k], t[(k + 1) % 3]
        edge[(min(a, b), max(a, b))].append(i)
done = [False] * len(tris)
for s in range(len(tris)):
    if done[s]:
        continue
    done[s] = True
    stack = [s]
    while stack:
        i = stack.pop()
        t = tris[i]
        for k in range(3):
            a, b = t[k], t[(k + 1) % 3]
            for j in edge[(min(a, b), max(a, b))]:
                if j == i or done[j]:
                    continue
                u = tris[j]
                if any((u[m], u[(m + 1) % 3]) == (a, b) for m in range(3)):
                    tris[j] = [u[0], u[2], u[1]]
                done[j] = True
                stack.append(j)


def cross(p, q):
    return (p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0])


signed = sum(sum(a * b for a, b in zip(at[t[0]], cross(at[t[1]], at[t[2]]))) for t in tris) / 6
if signed < 0:  # inward: turn every triangle round
    tris = [[t[0], t[2], t[1]] for t in tris]
lines = ["solid shell_side"]
for t in tris:
    p, q, r = (at[i] for i in t)
    n = cross([q[k] - p[k] for k in range(3)], [r[k] - p[k] for k in range(3)])
    length = math.sqrt(sum(x * x for x in n)) or 1.0
    lines.append("facet normal %.6e %.6e %.6e\nouter loop" % tuple(x / length for x in n))
    lines += ["vertex %.9e %.9e %.9e" % x for x in (p, q, r)]
    lines.append("endloop\nendfacet")
lines.append("endsolid shell_side")
open(out, "w").write("\n".join(lines) + "\n")
print(f"{out}: {len(tris)} facets, volume {abs(signed):.0f} mm3 (Gmsh: {volume:.0f})")
