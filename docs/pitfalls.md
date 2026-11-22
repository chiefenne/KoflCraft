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

## Fitted volume differs from the STL

The polyhedral duct mode remeshes the wall along curved pipes and fillets.
On the SolidCraft fixture the fitted fluid volume is 0.183 % below the STL
volume. Report both values and the difference as a limitation; never call
the fitted mesh "volume conserving" with respect to the STL.

## Gmsh format versions

KoflCraft writes MSH 2.2; FoehnCraft reads only MSH 4.1. A file that Gmsh
itself opens is not proof that FoehnCraft will. Test the hand-off with
FoehnCraft's reader (FC-06).

## Counts a file declares

A reader that allocates from a count declared in the header can abort on a
broken file. FoehnCraft found this in four of its mesh readers. Compare the
declared count with the remaining file size before allocating.

## Tetrahedral subdivision is not the polyhedral mesh

The Gmsh export subdivides polyhedra into tetrahedra (and the duct export
into pyramids and tetrahedra). A solver reading it sees a different cell
count and different face connectivity. To hand true polyhedra to a solver,
use VTU or OpenFOAM `polyMesh` (OPEN-12).
