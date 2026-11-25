# Block layout optimisation and wall rings for the multi-element airfoil (MB1g)

**Status:** built, 2026-11-25. Requested by the maintainer: bring the
multi-element template's largest non-orthogonality down from 86°, study
Sauer & Morsbach (2023) and newer methods, redesign the gaps and the slat
nose. The fix was chosen on a decision card: "Ring + layout optimizer
(paper)", general for every block topology, not only this template.
Code: [`block::layout`](../crates/mesh-algorithms/src/block/layout.rs) and
[`block::multi_element`](../crates/mesh-algorithms/src/block/multi_element.rs).

## 1. Diagnosis of the first template

Reference case (case A below): NACA 4410 slat 25° (gap 0.05, overlap
-0.04), NACA 4412 main, NACA 4412 flap 30° (gap 0.03, overlap -0.02),
64 surface cells, wall spacing 5e-4.

- **Spacing runs out along H-lines.** The H-topology carried every wall
  vertex's spacing straight to the box: the column line from a nose corner
  to the box bottom carried the corner's tangential spacing, the row lines
  to the inlet the wall spacing. Far from any wall the cells were 0.87 long
  and 4e-4 wide.
- **Taper plus aspect ratio.** On such a sliver a 7 % taper moves the cell
  centroid 5e-3 along the cell while the neighbour's centroid is 3e-4
  away, so the finite-volume non-orthogonality (face normal against the
  line between the two cell centroids) read 86° where the grid lines cross
  at right angles.
- **Gap fans.** Three row lines left the previous trailing edge within a
  few millimetres and had to reach three vertices spread over the next
  element's nose; TFI sheared these wedge blocks (75°).
- **The smoother cannot help.** The elliptic smoother freezes the TFI
  residual (MB1d), so away from the walls the TFI grid is a fixed point:
  TFI and elliptic gave the same 85.96°, and elliptic took 6 minutes
  (40,000 sweeps) for 20,752 cells. The fix has to be in the layout.

## 2. Candidates (ALG-16)

| | Method | Source | Outcome |
|---|---|---|---|
| CFD, chosen idea | Block corners and edge shapes found by optimising cell angles and growth on a coarse control grid | Sauer & Morsbach (2023) | Idea taken, built on our TFI grid (section 4) |
| Graphics | Quad layouts from integer-grid maps; cross fields | Bommes, Zimmer & Kobbelt (2009); Bommes et al. (2013); Lyon, Campen & Kobbelt (2021); Knöppel, Crane, Pinkall & Schröder (2013); Ray, Vallet, Li & Lévy (2008) | Deferred: they answer "which topology", not "where the corners of a given topology go"; candidates for automatic layouts later |
| CFD, field based | Cross-field block layouts with wall layers | Deng, Wang & Qin (2023); Roche et al. (2023) | Deferred with the graphics row; Roche's wall-layer-first structure is the ring below |
| Classical | Hand-placed corners, then elliptic smoothing | Thompson, Thames & Mastin (1974); the template as it was | Baseline of the benchmark |
| Multi-element practice | A thick wall wrap per element, joined to H-blocks or overset | Jones, Biedron & Whitlock (1995); Chi et al. (2002) | Confirms the ring; overset is not an option for point-matched blocks |

**What we take from Sauer & Morsbach.** An O-ring around each blade is
where wall orthogonality and spacing live; singular vertices stay off the
wall; the layout is the minimiser of an angle and growth criterion, not a
hand placement. **What we leave.** Their degree-2 B-spline intermediate
layer with singularity blending replaces TFI wholesale; the optimiser below
evaluates our own TFI grid instead, so it keeps the TFI + elliptic pipeline
and works on every topology the checker accepts. Their Levenberg–Marquardt
solve needs a smooth objective; the worst angle is not smooth in the
parameters (section 4). Only the paper was used; no source code was read.

## 3. Wall rings and the cove

Each element sits in an O-ring one block thick: six blocks whose inner
vertices are the six wall vertices (lower and upper trailing-edge corners,
middle of the blunt trailing edge, lower and upper corner, leading edge)
and whose outer loop is the wall offset along its normals by
`RING` = 5 % of the element's chord, at most `RING_GAP_SHARE` = 30 % of a
gap and at least `2 · RING_MIN_CELLS` wall spacings. The cells across the
ring grow by at most 1.2. So the near-wall cells are orthogonal by
construction, and the H-blocks meet a smooth loop instead of the wall's
corners and spacings (the "thick wrap-around" of the multi-element
practice).

**Cove.** Where the next element's nose sits below the previous element's
chord (a flap under the main element), the lower trailing-edge vertex
slides forward along the lower side, by the next element's lower nose
length and at most `COVE_REACH` = 40 % of the lower side. The gap blocks
then face the nose across the slot instead of wrapping round the trailing
edge. The cove wall takes the facing nose's spacings, scaled by the cove
over the nose length; with uniform spacing the ring cells there tilted by
74°.

## 4. Layout optimisation

`block::layout::optimize` works on any topology extruded along z (two
vertex layers):

- **Rigid blocks.** A block with a face on a rigid patch (a wall, here the
  rings) keeps its vertices and edges.
- **Parameters.** A vertex on the outer boundary slides along it where the
  boundary is straight; any other free vertex moves in the plane. An edge
  inside the domain bends by two Hermite modes, 6.75·s(1−s)² and
  6.75·s²(1−s) times the chord normal, on top of its template shape; that
  sets the angles at which it leaves its ends. The end spacings of Tanh
  and Geometric laws scale in log space within a factor of 20.
- **Cost per face** of the TFI grid: (θ/45°)² + (θ/40°)¹² + (ln g / ln
  2.5)², θ the non-orthogonality, g the growth across the face. The
  twelfth power lets the worst faces dominate without making the cost a
  maximum. A folded corner costs 1e4, and a corner sharper than sine 0.1
  costs 100 · ((0.1 − sin) / 0.1)²; without that continuous term the search
  could not see a way out of a fold, since the fold count is discrete.
- **Search.** Pattern search (Hooke & Jeeves 1961), one parameter at a
  time, steps halving from 2 % of the vertex's shortest edge down by a
  factor 100. A trial re-grids only the blocks the parameter touches and
  the faces between them, so a sweep costs about one grid of the whole
  topology per parameter. Every accepted move keeps the corners
  right-handed, checked with the edge end tangents as the topology check
  does. 40 sweeps at most.
- **Output.** The result is written back into the topology: vertices on
  both layers, bent edges as splines, scaled laws. The report gives the
  before and after values the written topology grids to; the template
  passes it to `block.template`'s result under `layout`.

The pattern search was chosen over a gradient method because θ and the
growth are piecewise smooth in the parameters only (the worst face changes
from trial to trial), and a trial is cheap: 0.5 s to 2.5 s for 40 to 96
parameters and 2,800 to 24,000 cells.

## 5. Benchmark

TFI grid of the template, 2D finite-volume non-orthogonality on the
one-cell extrusion (face normal against the centroid line), the same
measure as FoehnCraft's checkMesh-style report. Before: the template on
main (commit 5220e23). After: rings, cove and layout optimisation. Time is
the wall time of the template build (with the optimisation) and the TFI
grid, release build, one core.

Main NACA 4412, slat NACA 4410 (0.16 chord), flap NACA 4412 (0.3 chord);
deflection in degrees, gap and overlap in main chords; 64 surface cells
and wall spacing 5e-4 unless noted.

| Case | Elements | Before: max / mean | Before: cells > 70° | After: max / mean | After: cells > 70° | Cells after | Time |
|---|---|---|---|---|---|---|---|
| A | slat 25° (0.05, -0.04), flap 30° (0.03, -0.02) | 85.96° / 17.40° | 240 | 54.72° / 14.46° | 0 | 23,782 | 2.2 s |
| B | flap 20° (0.025, -0.015) | 76.11° / 14.98° | 89 | 51.22° / 14.37° | 0 | 19,232 | 1.7 s |
| C | slat 15° (0.03, -0.025) | 86.35° / 14.02° | 200 | 59.92° / 14.55° | 0 | 18,528 | 1.8 s |
| D | slat 20° (0.04, -0.03) | refused (folds in the slat gap) | | 65.11° / 14.70° | 0 | 18,688 | 1.9 s |
| E | main only | 76.11° / 12.15° | 58 | 47.38° / 9.57° | 0 | 15,520 | 1.4 s |
| F | NACA 0012 main, slat 25° (0.05, -0.04); 24 cells, spacing 2e-3 | 89.45° / 16.32° | 66 | 64.63° / 12.73° | 0 | 2,808 | 0.5 s |
| G | as A; 32 cells, spacing 1e-3 | 89.14° / 18.54° | 117 | 61.95° / 15.35° | 0 | 5,652 | 0.8 s |
| H | slat 30° (0.02, -0.015), flap 40° (0.015, -0.01) | refused (folds in the slat gap) | | 68.22° / 15.99° | 0 | 22,912 | 2.3 s |
| I | slat 10° (0.08, -0.06), flap 10° (0.05, -0.04) | 86.21° / 13.92° | 168 | 54.95° / 13.27° | 0 | 24,306 | 2.5 s |

The rings add cells (A: 20,752 to 23,782). Of the improvement, the
optimiser does most: with the rings and the cove but before the
optimisation, case A reads 89.75° with 65 cells over 70° and 10 folded
corners, case H 89.95° with 122 cells over 70° and 46 folded corners; the
optimiser unfolds them all. Every case runs 40 sweeps of 2,300 to 5,100
trials over 40 to 96 parameters.

**Elliptic, case A.** The template's elliptic smoothing now converges in
1,310 sweeps (34 s) instead of running into the 40,000 sweep limit (6
minutes): max 54.85°, mean 14.23° (TFI 54.72°, 14.46°); wall angle max
44.9°, mean 4.6° (TFI 44.1°, 16.8°).

**Where the rest sits.** The remaining maxima are at the box top and in
the gap blocks, away from the walls. The largest growth across a face
(17 to 28) sits in the ring at the corner cells of the blunt trailing
edge, where the trailing-edge face cells meet the ring's wall-normal
spacing; the rings are rigid, so the optimiser does not touch it.

## 6. Tried and dropped

- **Offset-consistent outer laws** (the ring's outer loop carrying the
  wall law, scaled by the offset): 1,267 cells over 70° in the main-only
  case, because a Tanh law with a large first cell (the corner arc) cannot
  follow it.
- **Unclamped nose corners** (corner vertices free on the ring): failed the
  spacing and right-handedness checks.
- **Trailing-edge spacing clustering** on the H-lines: the mean got worse.
- **Uniform cove wall spacing:** the 74° ring tilt of section 3.

## 7. Open

- The triangle (tri-center) split, three quads around a valence-3 vertex,
  2D and 3D, with an automatic trigger: the maintainer's next PR, for every
  block topology.
- The gap blocks still carry the worst angles of the multi-element cases;
  the split above is meant for them.
- Overlapping elements (a slot topology) are still refused; the fold check
  stays, and its message names the gap to change (slat 10°, gap 0.02,
  overlap -0.005 still folds).
- The optimiser treats extruded topologies only; a 3D version needs face
  shape parameters.

## Sources

- Sauer, M. & Morsbach, C. (2023). An optimization based
  multi-block-structured grid generation method. *International Journal
  for Numerical Methods in Engineering*. doi:10.1002/nme.7308.
- Hooke, R. & Jeeves, T. A. (1961). "Direct search" solution of numerical
  and statistical problems. *Journal of the ACM* 8(2), 212–229.
- Thompson, J. F., Thames, F. C. & Mastin, C. W. (1974). Automatic
  numerical generation of body-fitted curvilinear coordinate system.
  *Journal of Computational Physics* 15, 299–319.
- Bommes, D., Zimmer, H. & Kobbelt, L. (2009). Mixed-integer
  quadrangulation. *ACM Transactions on Graphics* 28(3).
- Bommes, D., Campen, M., Ebke, H.-C., Alliez, P. & Kobbelt, L. (2013).
  Integer-grid maps for reliable quad meshing. *ACM Transactions on
  Graphics* 32(4).
- Lyon, M., Campen, M. & Kobbelt, L. (2021). Quad layouts via constrained
  T-mesh quantization. *Computer Graphics Forum* 40(2).
- Knöppel, F., Crane, K., Pinkall, U. & Schröder, P. (2013). Globally
  optimal direction fields. *ACM Transactions on Graphics* 32(4).
- Ray, N., Vallet, B., Li, W. C. & Lévy, B. (2008). N-symmetry direction
  field design. *ACM Transactions on Graphics* 27(2).
- Deng, X., Wang, Y. & Qin, N. (2023). Cross-field based block layouts
  for airfoils with wall layers. *AIAA Journal* (White Rose eprint
  204784).
- Roche, J. et al. (2023). Advancing-front block layouts from a blended
  wall-normal field. *International Meshing Roundtable* 2023, paper 11.
- Jones, W. T., Biedron, R. T. & Whitlock, M. E. (1995). Multi-element
  airfoil grids, NASA report, NTRS 19970013280.
- Chi, X. et al. (2002). Multi-element airfoil grids, NASA report, NTRS
  20030007845.

The last four are named by subject, not by their exact titles; the
identifiers find them.
