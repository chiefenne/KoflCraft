"""Rebuild this fixture's tagged MSH 2.2 scaffold using Gmsh and Python's stdlib.

Usage: python3 prepare_scaffold.py /path/to/new-scaffold.msh
"""

import hashlib
import json
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

CASE = Path(__file__).resolve().parent


def ranges_to_labels(count, groups, default):
    labels = [default] * count
    for label, ranges in groups:
        for start, end in ranges:
            labels[start:end] = [label] * (end - start)
    return labels


def chart(vertices, region, patch):
    if patch != 3:
        return (patch,)
    # Split each pipe at its equator so a smooth half tube has a single chart.
    if region != 2:
        height = 40 if region in (0, 1) else 0
        return (patch, region, sum(p[2] - height for p in vertices) >= 0)
    # Retain distinct planar chamber faces across all sharp edges.
    a, b, c = vertices
    u = [b[i] - a[i] for i in range(3)]
    v = [c[i] - a[i] for i in range(3)]
    normal = [u[(i+1)%3]*v[(i+2)%3] - u[(i+2)%3]*v[(i+1)%3] for i in range(3)]
    axis = max(range(3), key=lambda i: abs(normal[i]))
    return (patch, region, axis, normal[axis] > 0)


def write_charts(path):
    data = (CASE / "SolidCraft_vs_FoehnCraft.stl").read_bytes()
    plan = json.loads((CASE / "partition.json").read_bytes())
    patch_map = json.loads((CASE / "patches.json").read_bytes())
    if hashlib.sha256(data).hexdigest() != plan["source_sha256"]:
        raise ValueError("unexpected source STL")
    count, = struct.unpack_from("<I", data, 80)
    regions = ranges_to_labels(count, ((i, r["facets"]) for i, r in enumerate(plan["regions"])), 0)
    patches = ranges_to_labels(count, ((i+1, [p["facet_range"]]) for i, p in enumerate(patch_map["ports"])), 3)
    parent = list(range(count))

    def root(i):
        while i != parent[i]:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    points, pending, surfaces = {}, {}, {}
    triangles, groups = [], []
    rim = set()
    for i, record in enumerate(struct.iter_unpack("<12fH", memoryview(data)[84:])):
        vertices = [record[j:j+3] for j in (3, 6, 9)]
        nodes = tuple(points.setdefault(p, len(points)+1) for p in vertices)
        triangles.append(nodes)
        groups.append(chart(vertices, regions[i], patches[i]))
        for j, a in enumerate(nodes):
            edge = tuple(sorted((a, nodes[(j+1)%3])))
            other = pending.pop(edge, None)
            if other is None:
                pending[edge] = i
            elif groups[i] == groups[other]:
                parent[root(other)] = root(i)
            if patches[i] != 3:
                rim.symmetric_difference_update([edge])
    if pending or len(rim) != 480:
        raise ValueError("expected closed shell and both original 240-edge rims")
    with path.open("w") as out:
        print("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes", len(points), sep="\n", file=out)
        for pos, index in points.items():
            print(index, *pos, file=out)
        print("$EndNodes\n$Elements", count + len(rim), sep="\n", file=out)
        for i, nodes in enumerate(triangles):
            tag = surfaces.setdefault(root(i), len(surfaces)+1)
            print(i+1, 2, 2, patches[i], tag, *nodes, file=out)
        for index, edge in enumerate(sorted(rim), 1):
            print(count+index, 1, 2, 0, index, *edge, file=out)
        print("$EndElements", file=out)


def main(output):
    output = output.resolve()
    if output.exists():
        raise FileExistsError(output)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        write_charts(work / "charts.msh")
        geometry = work / "mesh.geo"
        geometry.write_bytes((CASE / "prepare-scaffold.geo").read_bytes())
        subprocess.run(["gmsh", geometry, "-3", "-nt", "1", "-o", output], check=True)
    print("SHA-256", hashlib.sha256(output.read_bytes()).hexdigest())


if __name__ == "__main__":
    main(Path(sys.argv[1]))
