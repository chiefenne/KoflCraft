# Pitfalls: lessons and experiments not to repeat

This page keeps what earlier work found out the hard way (PROC-02), so that
nobody spends time finding it again. The rules that follow from it are in
[AGENTS.md](../AGENTS.md). Add a lesson in the pull request that learns it:
what was tried, what it measured, and what to do instead.

## Clipping against nonconvex domains

Clipping Voronoi cells against a nonconvex domain or a triangle soup does
not give valid CFD cells. The M1 clipper assumes a convex box or a convex
validated STL shell. For nonconvex shapes, use the partitioned path
([nonconvex-duct.md](nonconvex-duct.md)) or a boundary-conforming method
(section 13.3 of the requirements).

## Local cells are not a global mesh

A set of locally valid convex polyhedra is not a conforming mesh. Shared
faces and vertices must be reconciled, and conservation checked, after
construction (ALG-12). GPU papers that skip this report volume errors up to
0.2 % (Liu et al. 2022): not acceptable here.

## Gmsh scaffolds are not reproducible across preparations

Gmsh optimisation can change node positions and counts between two
preparations of the same `.geo` file. Pin the checked-in scaffold and its
SHA-256 in the partition file; regenerate only on purpose, and update the
digest in the same commit.

## A fine cell size does not mean a fine wall

The polyhedral duct takes its wall from the volume scaffold, not from the
STL. The Gmsh scaffold of Manifold remeshed the 1-unit STL wall into
triangles of about 6 units (curvature sizing off), and every finer core or
layer cell sits on those flat triangles: 0.54 units off the STL, about seven
first-layer thicknesses. Gmsh's curvature sizing had no effect on the
discrete charts, and a chord size field made it worse (0.72: the
reparametrisation of long pipe charts stretches triangles to 12 units). Using
the STL triangles directly inverts tetrahedra at the wall offset (28-unit
slivers). Measure the wall against the source both ways
(`report.scaffold.wall_deviation`) and bound it with
`max_wall_deviation`; build the scaffold in Rust (`volume_scaffold.build`),
which remeshes the wall to a chord tolerance.

## Delaunay wall slivers survive refinement

A Delaunay tetrahedralization of a curved, finely remeshed wall has flat
tetrahedra whose four corners are all on the wall (two wall triangles and
the opposite diagonal). Their circumspheres bulge outward, so no interior
point removes them, and inserting points just inside the wall only adds
points. Remove them by edge removal around their diagonal (best
triangulation of the ring by dynamic programming); on Manifold that takes
the worst mean ratio from 0.002 to 0.05.

## Fitted volume differs from the STL

The polyhedral duct mode remeshes the wall along curved pipes and fillets.
On the Manifold fixture the fitted fluid volume is 0.183 % below the STL
volume. Report both values and the difference as a limitation; never call
the fitted mesh "volume conserving" with respect to the STL.

## Gmsh format versions

FoehnCraft reads only MSH 4.1, so KoflCraft writes 4.1 by default and 2.2
only on request (`gmsh-2.2`). A file that Gmsh itself opens is not proof
that FoehnCraft will. Test the hand-off with FoehnCraft's reader (FC-06).
In 4.1 every element belongs to an entity and each entity names one
physical group; a writer that gives one entity elements of two groups
gets an error, not a silently regrouped file.

## Counts a file declares

A reader that allocates from a count declared in the header can abort on a
broken file. FoehnCraft found this in four of its mesh readers. Compare the
declared count with the remaining file size before allocating.

## Whole files in memory

The writers once built each file in a `Vec<u8>` and the bundle wrote them
all at the end, so a 2.5 million cell sweep held its 3.4 GB of files next
to its meshes (peak 6.5 GB; about 12 GB of `foehncraft.vtu` alone at 10
million polyhedra). Hashing for the evidence wanted the bytes, which is why
they were kept. Stream instead: `MeshDirectory::stream_new` and
`Bundle::stream` hash while writing and return the digest. Inside a writer,
gather a line in a buffer and write it at once: a `write!` per number into
a `&mut dyn Write` costs more than the number. Keep per-cell data flat
(`HandoffMesh::cell_face_list`), not a `Vec` per cell.

## Tetrahedral subdivision is not the polyhedral mesh

The Gmsh export subdivides polyhedra into tetrahedra (and the duct export
into pyramids and tetrahedra). A solver reading it sees a different cell
count and different face connectivity. To hand true polyhedra to a solver,
use VTU or OpenFOAM `polyMesh` (OPEN-12).

## Hash iteration order makes a mesh differ between runs

`std::collections::HashMap` and `HashSet` seed their hasher randomly per process. Code that iterates one (a Delaunay cavity, a set of faces to split) and lets the order decide tetrahedron indices or insertion order gives a different mesh on every run of the same input. The Manifold scaffold's smallest core sphericity moved between 0.03 and 0.17 that way. The scaffold modules use the fixed-hasher aliases `delaunay::HashMap` and `delaunay::HashSet`; use them, or a `BTreeMap`, wherever iteration order reaches the output.
