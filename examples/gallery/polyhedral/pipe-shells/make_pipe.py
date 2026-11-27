"""Write pipe.stl: a round pipe of radius 5 and length 40 along z, its wall a
96-gon, its two ends flat fan-triangulated discs (ASCII STL, 384 facets)."""
import math
import sys

R, L, N = 5.0, 40.0, 96
v = [(0.0, 0.0, 0.0), (0.0, 0.0, L)]
for i in range(N):
    a = 2 * math.pi * i / N
    v += [(R * math.cos(a), R * math.sin(a), 0.0), (R * math.cos(a), R * math.sin(a), L)]
b = lambda i: 2 + 2 * (i % N)
t = lambda i: 3 + 2 * (i % N)
faces = []
for i in range(N):
    faces += [(b(i), b(i + 1), t(i + 1)), (b(i), t(i + 1), t(i)), (0, b(i + 1), b(i)), (1, t(i), t(i + 1))]


def sub(p, q):
    return [p[k] - q[k] for k in range(3)]


out = ["solid pipe"]
for f in faces:
    p, q, r = (v[k] for k in f)
    u, w = sub(q, p), sub(r, p)
    n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
    s = math.sqrt(sum(c * c for c in n))
    out.append("facet normal %.9e %.9e %.9e\nouter loop" % tuple(c / s for c in n))
    out += ["vertex %.9e %.9e %.9e" % x for x in (p, q, r)]
    out.append("endloop\nendfacet")
out.append("endsolid pipe")
open(sys.argv[1] if len(sys.argv) > 1 else "pipe.stl", "w").write("\n".join(out) + "\n")
