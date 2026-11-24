# Symbolic vertices for the Voronoi core

Design note for the first hardening step of the Voronoi core ([m4-voronoi-core](m4-voronoi-core.md)) on thin CAD-like parts. Question: why does the core stop on the cooling jackets with "piece of site S in tetrahedron T has no continuation in O", and what makes the cuts of neighbouring scaffold tetrahedra agree? The code is [`duct/symbolic.rs`](../crates/mesh-algorithms/src/duct/symbolic.rs) and the final diagram in [`duct/voronoi_core.rs`](../crates/mesh-algorithms/src/duct/voronoi_core.rs).

## The failure

Measured on the standard cooling jacket at `size_factor` 2 (a confidential part, so only counts and lengths are given):

- 2 of 855,210 pieces had no continuation, and 1,300 pieces were dropped as collapsed.
- Both missing pieces are needles of mirrored sites at reentrant edges: 0.5 to 1.0 long, 1e-8 wide, volume 4e-18 to 5e-17. The bisector of a mirrored pair is the plane that halves the wall angle, and it is only accurate to about 1e-8 because STL normals are single precision.
- The piece exists in both tetrahedra when clipped with tolerance zero. With the tolerances (clip 1e-12 × scale, weld 1e-10 × scale = 2e-8) one tetrahedron dropped it after the weld collapsed all its faces, while its neighbour kept a face with three distinct vertices.
- The thin jacket stopped on main with the same error (2 missing continuations among 2,039,597 pieces, shared face 2.8e-8 across, neighbour piece volume 6e-25).

So the reading was right: each tetrahedron decided on its own rounded coordinates whether a thin piece exists, and two neighbours decided differently. The clip tolerance and the weld both took part.

## Candidates

| | Method | Source |
|---|---|---|
| Classical | Keep coordinates as the identity of a vertex and tune the tolerances, or patch the weld so that a face seen from one side survives on the other | |
| Graphics | Name every vertex by the planes that define it, decide from the name, and compute coordinates once per name | Lévy & Liu (2010), *Lp centroidal Voronoi tessellation and its applications* (symbolic vertices of the restricted diagram); Yan, Wang, Lévy & Liu (2013); Lévy (2016), *Robustness and efficiency of geometric programs: the Predicate Construction Kit (PCK)* |

No source code was read (LIC-12 to LIC-15). The construction below is built from the papers and from first principles.

## Construction

1. **Names.** A vertex of a piece of site s is named by its three defining planes other than "s itself": a scaffold node (`V`), a scaffold edge and one other site (`E`), a scaffold face and two other sites (`F`), or three other sites (`X`). The name together with s gives four sites or nodes, sorted, so two tetrahedra and two sites build the same name for the same point.
2. **Position from the name.** The coordinates are a function of the name only: the point of the scaffold edge at equal distance from two sites, the point of the face at equal distance from three, the circumcentre of four. Two tetrahedra therefore get bit-identical coordinates on their common face.
3. **Sides from the name.** Whether a vertex is inside, on or outside the bisector towards another site o is decided from the name. A vertex whose name contains o is on it. A scaffold node within the clip tolerance of the bisector is on it. A vertex on a scaffold edge or face whose nodes are all on the bisector is on it.
4. **Weld.** All named points of the final diagram are welded once, globally, by single linkage (union-find over a grid hash) with the old tolerance. Scaffold nodes anchor their class. The result does not depend on the order of the points. A piece that loses all faces is dropped, and it is then gone on every side.
5. **Leftovers of the weld.** A part with fewer than four exposed faces (a film) joins the neighbour with the largest contact. A wall polygon below 1e-6 of the largest polygon its cell has on the same wall triangle is merged into that polygon. A wall polygon whose first fan triangle is degenerate is rotated to its best corner, because the wall columns fan from the first corner. A pinched ring is cut into simple loops and not dropped.

Lloyd still uses the floating-point pieces. Only the final diagram of each round is symbolic.

### Two findings on the way

- **A pure position is not always usable.** The circumcentre of four nearly cocircular sites lands far from the piece. The position from the name is therefore used when it lies within 1e4 × the clip tolerance of the crossing of the cut edge, and the edge crossing otherwise. In that case two tetrahedra can differ by rounding, and the weld joins them.
- **The wall seeds are ties by design.** The conforming wall seeds ([m4-wall-columns](m4-wall-columns.md)) put bisectors exactly through wall vertices and wall edges. The floating-point sign of such a test is noise of 1e-15. Rule 3 exists for this: the tie is detected at the scaffold nodes and inherited by every vertex named on their edges and faces.

## Results

Release build, Linux x86_64.

| Case | Before | With symbolic vertices |
|---|---|---|
| Elbow | 1,004 polyhedra | 1,004 polyhedra |
| Manifold | valid, 95,666 polyhedra, 71,687 core cells, sphericity minimum 0.531 | valid, 95,938 polyhedra, 71,752 core cells, sphericity minimum 0.558, no welded vertex, no film, no inverted face |
| Standard jacket, `size_factor` 2 | core stops: "no continuation" | core completes in all four rounds; the run stops in the final checks of the wall columns |
| Thin jacket, `size_factor` 2 | core stops: "no continuation" | not finished when this note was written |

The standard jacket does not yet produce a mesh. What stops it now is downstream of the core: 23 layer sub-wedges are inverted where a wall polygon has a speck of about 1e-5 near a wall vertex and is not convex, and the wall patch area is 0.22 too large in the check (3e-6 relative) because the check sums unsigned fan areas of such a polygon. 600 wall polygons keep a straight vertex.

### Variants tried and rejected

| Variant | Outcome |
|---|---|
| Remove flat faces and zip the edges across them | Open wall ring edges. Three points 2 to 3e-8 apart have no consistent order along an edge |
| Wider tie band (1e-9 × scale) | "bisector face has no partner": real crossings are then called ties on one side |
| Larger weld (1e-8 × scale) | No effect on the failing cells, and it would join distinct points that the planarity check (1e-10 × diameter) tells apart |
| Merge all wall polygons of a cell on one triangle | Area excess and straight vertices on Manifold; only needles are merged now |

### Follow-up: wall polygons that are not convex

The wall columns and the patch area check fanned every wall polygon from its first corner, which is only valid for a convex polygon. The columns now triangulate each polygon by ear clipping (`ears` in [`duct/columns.rs`](../crates/mesh-algorithms/src/duct/columns.rs)): straight corners (below 1e-9 of the longest edge squared) are taken out first and put back by splitting the triangle on their edge, and among the free ears the roundest is cut first. The area check uses the vector area of a boundary face. With that all layer cells of the standard jacket pass (23 sub-wedges were inverted), and the elbow and Manifold are unchanged in their counts.

The jacket then stops in the last check, on core cells. A census with the check made non-fatal found 61 core faces beyond the planarity limit of 1e-10 × diameter (32 beyond 1e-9, 4 beyond 1e-8, largest 2.3e-8), only 8 of them with an edge below 1e-6, 2 faces seen from behind by a rounding amount and 1 cell with an open shell. A weld of 1e-12 × scale instead of 1e-10 makes it worse (needle wall polygons, open wall rings). The cause of the off-plane points is not yet established.

The maintainer then chose a planarity limit of 1e-7 of the face size (`PLANAR_SHARE` in [`duct/verify.rs`](../crates/mesh-algorithms/src/duct/verify.rs)); a face within it is fanned from its vertex mean, as FoehnCraft fans every face. That covers all 61 faces. What remains on the jacket, measured with the checks made non-fatal: 5 faces with one fan triangle seen from behind by 1e-11 to 1e-9 (each at an edge of 3e-8 to 2.6e-5 whose end points are a few 1e-9 off the plane), 3 pairs of cells with a missing triangle or quadrilateral between them (an open shell), and at least one cell whose volume differs from its pieces.

## Step 2: exact sides

After the follow-up above the standard jacket still stopped in the last check. The measured cause: for about 18 core faces the two cells held different rings for the same face, because the side of a vertex was still decided by tolerance from each cell's own view. Taking the rings of both cells is worse (edges used three times).

**Method** (Lévy's side predicates of the PCK, 2016, written for clipping per tetrahedron; from the paper, no code read). A named vertex q lies on a scaffold simplex (a node, an edge, a face or the tetrahedron) with corners p_j and on the bisectors of its site s with one to three sites k_i. With the offset f_k(p) = |p − k|² − |p − s|², which is affine in p, the barycentric weights λ of q solve Σ λ = 1 and Σ λ f_{k_i}(p_j) = 0. The side of q of the bisector of s and o is the sign of Σ λ f_o(p_j), by Cramer's rule the quotient of two determinants (up to 4 × 4) of offsets at the corners. `offset_side` in [`predicates.rs`](../crates/mesh-algorithms/src/predicates.rs) evaluates both signs with a floating-point filter first and with Shewchuk's expansion arithmetic when the filter cannot decide. The answer is in, out or exactly on. It uses only the name, so every tetrahedron gives the same answer, and it is the same from every site of the name, because the linear systems differ by row operations only.

**Wall seed ties.** Seeds on the ball of one wall vertex are equidistant from it by construction, but only to rounding in stored coordinates. Declaring single offsets zero is not enough: two tied seeds would still see different distances to a third site. Each squared distance from a scaffold vertex is therefore taken to a stand-in site: the lowest seed on the ball of that vertex for every seed on it, the site itself otherwise. All offsets then come from one table of squared distances, which is all the argument above needs. The seeds are recognised in every rebuild round by their coordinates, which never change.

**Security radius.** The walk over the neighbouring sites stops when the next site is farther than twice the reach of the piece. The reach came from rounded vertex coordinates, and a vertex where three planes nearly share a line can be far from its computed place. Two pieces of a mirrored pair were then cut by different third sites ("bisector face has no partner"). Each vertex now carries an upper bound on its distance from the site: no more than the farther end of the edge it was cut on, and no more than its computed distance plus the rounding of the cut and the doubt of the two ends.

**Results** (release build; jacket figures with the final checks made non-fatal to count every failure):

| Standard jacket, `size_factor` 2 | Tolerance sides | Exact sides |
|---|---|---|
| Cell pairs with an open shell | 3 | 0 |
| Core faces with a fan triangle seen from behind | 5 | 2 |
| Cells with a volume mismatch | several, behind the open shells | 1 (5e-6) |
| Inverted layer sub-wedges | 0 | 7 (three columns, at a wall polygon edge of 6e-6) |

The elbow has 1,004 cells as before. Manifold is valid with the same counts as before (95,938 polyhedra, 71,752 core cells), so on a part without thin features the exact tests reproduce the rule-based result. The jacket does not finish yet.

**Cost.** The long Manifold test (elbow and Manifold builds with their checks) takes 233 s against 158 s before; the jacket build 5.5 minutes against 2.5. Part of it is the exact arithmetic, part the more careful stop rule, which clips against more sites. The two have not been separated.

## Where the work stopped (2026-11-24)

The maintainer closed the hardening for now. State on the standard jacket at size factor 2 (about 335,000 cells, part about 200 mm long, cells 0.2 to 1 mm), counted with every final check made non-fatal:

| Finding | Count | Size |
|---|---|---|
| Wall columns with inverted layer sub-wedges | 3 (7 sub-wedges) | at wall polygon edges of 6.3e-6 and 2.9e-8 |
| Core faces with a fan triangle seen from behind | 2 | signed area −8e-12 and −8e-10, at edges of 3.9e-8 and 2.6e-5 |
| Volume mismatch | 1 | −5.4e-6 relative |

The mesh is not exported. Elbow and Manifold are valid.

**Cause, traced for the wall columns.** The names of the vertices and the side tests are exact, the coordinates are not. `Frame::position` evaluates the name in floating point and, where that lands more than 1e4 tolerances from the crossing of the clipped edge, takes the crossing, which itself comes from rounded points. At a vertex that its bisectors fix badly the coordinates are off by up to about 1e-5. Two true vertices of a cell can be closer than that (four or five sites almost on one sphere), and the polygon then folds: the two wall polygons that share the 6.3e-6 edge each hold one vertex on the wrong side of it, and ear clipping has no positive triangle left for that corner.

**Experiment, not merged.** Coordinates from exact weights: the barycentric weights of a named vertex are λ_j = C_j / Σ C with C_j the cofactors of the row of ones in the determinant of the side test; with expansions for the C_j, Σ C_j p_j and Σ C_j are exact and one division rounds. A floating point filter keeps the old formula where the weights are good to 1e-11. Result on the jacket: no inverted sub-wedge is left (3 columns before), but 10 core faces have a fan triangle seen from behind (2 before; edges 2e-8 to 6e-6, signed area −8e-12 to −4e-9) and one cell is off in volume by −1.2e-8. Why the core faces got worse is not known; candidates are the global weld (1e-10 of the model size, about 2e-8 here, the length of most of these edges) and the fan from the vertex mean on faces with such an edge. Manifold and elbow were not run with it. The patch is kept outside the repository.

**Open options.**

1. Exact coordinates as above, after finding the second cause in the core faces.
2. Merge vertices of a cell closer than a limit. Removes the tiny edges themselves, but it is a tolerance decision that changes which cells touch, so it has to be taken once per vertex, not per cell.
3. Open, not checked here: a jacket set up afresh (`detect_ports` and `setup_duct`, tolerance 0.026 instead of the 0.0442 of the runs above) was reported not to leave the wall remesh within 15 minutes.

Also still open: the thin jacket runs out of memory in the core stage under 3.5 GB, and the jacket build takes 5.5 minutes with exact sides (2.5 before).

## What the literature offers beyond this step

Read on 2026-11-24 from the papers (the arXiv texts were cut at about 40,000 characters, so later sections were read as summaries; none of the texts read names a code licence).

- **Exact sides.** Lévy's PCK gives side tests for points defined as bisectors meeting an edge, a triangle or a tetrahedron, with an arithmetic filter, Shewchuk expansions and a globally consistent symbolic perturbation; exact arithmetic costs 10 to 50 times more and is normally needed in under 1% of calls ([preprint](https://core.ac.uk/download/49461268.pdf)). Attene's indirect predicates keep constructed points implicit and evaluate predicates through a filter cascade without division ([arXiv 2105.09772](https://arxiv.org/abs/2105.09772)). Either fits the names of this step, and [`predicates.rs`](../crates/mesh-algorithms/src/predicates.rs) already has the expansions. At the wall seeds nearly every test is a tie, so the filter would fail there; rule 3 or a perturbation has to stay.
- **No clipping at the wall.** VoroCrust places seeds in inside/outside pairs so that the wall consists of bisector faces ([arXiv 1902.08767](https://arxiv.org/abs/1902.08767)). That would remove the designed ties, at the price of a new core boundary. The paper notes that narrow regions need many cells. The patent check LIC-08 named in [m4-voronoi-core](m4-voronoi-core.md) applies.
- **Boundary recovery.** Diazzi, Panozzo, Vaxman & Attene (2023) recover constrained segments with implicit Steiner points, with a convergence guarantee for segments, and report no failure on 4,408 valid Thingi10k models ([arXiv 2309.09805](https://arxiv.org/abs/2309.09805)). It is a candidate for the recovery guarantee the scaffold still lacks ([thin-gaps](thin-gaps.md)).
- **Rounding stays open.** Lévy's mesh CSG paper ([arXiv 2405.12949](https://arxiv.org/abs/2405.12949)) and Diazzi & Attene ([arXiv 2109.14434](https://arxiv.org/abs/2109.14434)) both leave the conversion of exact points to floating point unsolved. A weld of some kind stays necessary.

## Limits

- Since step 2 the side decisions are exact. Vertex coordinates are still rounded, and the weld, the film merge and the planarity limit deal with what rounding leaves.
- A position falls back to the edge crossing for nearly cocircular sites, so the weld still carries part of the consistency.
- The jackets do not finish (see the follow-up above): a few fan triangles seen from behind and three missing faces between cells remain.
- Manifold's counts changed slightly, as they do with any change of the site layout.

## Regression

`neighbouring_tetrahedra_share_names_and_coordinates`, `two_sites_name_their_common_face_alike` and `a_scaffold_vertex_on_the_bisector_belongs_to_both_pieces` in `symbolic.rs`; `the_weld_is_the_same_for_every_order` and `a_pinched_ring_is_cut_not_dropped` in `voronoi_core.rs`; `side_of_a_voronoi_vertex_is_exact`, `both_cells_of_a_face_agree` and `declared_tie_is_on_the_bisector` in `predicates.rs`; `ears_of_notched_and_straight_rings_are_positive` in `columns.rs`; `notched_boundary_face_has_its_true_area` in `verify.rs`. The jacket is confidential, so no test in the repository reproduces the original failure on a whole part; the unit tests hold the properties it violated.
