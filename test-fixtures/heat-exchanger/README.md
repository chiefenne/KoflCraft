# Shell-and-tube heat exchanger (STEP)

[`SHELL_TUBE_HX.step`](SHELL_TUBE_HX.step) is the maintainer's test case for
the STEP reader (section 4b, STP-01 to STP-13) and the multi-domain mesher
(MD-18). It is stored byte for byte.

**Version 1 (2026-11-23), geometry only.** The maintainer did not set the
materials and names in this version and will replace the file with a
corrected one, which becomes the shipped test example. Until then, the
names, colours and materials below describe this file only. They are not
acceptance targets, but the reader must still handle them, including the
repeated body names. When the update arrives, replace the file, update its
SHA-256 and the inventory, and keep the geometry checks if the geometry is
unchanged.

**Two shipped variants (maintainer decision 2026-11-23).** The corrected
example comes as two files: an **open** one, whose nozzles end in annular
wall faces as in this file and which tests the automatic port caps (MD-19),
and a **closed** one, in which each nozzle opening is closed by a named
planar cap face inside the inner radius (`tube_inlet`, `tube_outlet`,
`shell_inlet`, `shell_outlet`). The closed one lets the STEP reader and the
fluid domains be tested before MD-19 is built, and both must give the same
mesh domains and port patches.

## Provenance (LIC-10)

| Item | Value |
|---|---|
| Supplied by | the maintainer (chiefenne), 2026-11-23, as the test file for the STEP reader and the multi-domain mesher |
| Created | 2026-11-22 in Autodesk Fusion (`originating_system` Autodesk Translation Framework v15.17.0.0, written with ST-Developer v20.2) |
| Schema | AP214 (`AUTOMOTIVE_DESIGN { 1 0 10303 214 3 1 1 }`), ISO 10303-21 clear text |
| Size, SHA-256 | 246,554 bytes, `4e073c630cc8d03b5cde7c37d30a6b83f5653df8c81d32396bee54493d28cb9a` |
| Licence | supplied for use as a KoflCraft test fixture; the geometry is the maintainer's own model |

## What the file contains

Inspected on 2026-11-23 with a throw-away entity count, not with a KoflCraft
reader. The numbers below are what the reader must reproduce.

- **Units:** length in millimetres (`SI_UNIT(.MILLI.,.METRE.)`), angles in radians,
  mass in kilograms. Distance accuracy 0.01 mm (`UNCERTAINTY_MEASURE_WITH_UNIT`).
- **Assembly:** one root product `SHELL_TUBE_HX` with three parts, each used once
  with an identity placement: `TUBES:1`, `HOUSING:1`, `BAFFLES:1`.
- **Bodies:** 13 manifold solid B-reps (13 closed shells, 117 faces):
  `TUBES` has 7 bodies, `HOUSING` has 1 and `BAFFLES` has 5. The body names are
  Fusion's German defaults `Körper1` to `Körper7`, written with the escape
  `K\X\F6rper…`. They are **not unique**: `Körper1`, `Körper2`, `Körper4`,
  `Körper5` and `Körper7` occur in more than one part.
- **Colours** (one per part, on every body): `TUBES` grey
  `Kunststoff - glänzend (Grau)`, `HOUSING` red `Farbe - Emaille glänzend (Rot)`,
  `BAFFLES` yellow `Farbe - Emaille glänzend (Gelb)`.
- **Material:** each part carries an AP214 material property with the name
  `Kunststoff` (plastic) and density 1290 kg/m³. All three parts have the **same**
  material in the file.
- **Surfaces:** 41 planes, 64 cylinders, 4 cones, 4 tori, 4 B-spline surfaces.
  **Curves:** 95 lines, 136 circles, 16 B-spline curves.
- **No fluid bodies.** The two fluids are the empty space of the model.

## Geometry (y is the axis of the exchanger, mm)

| Part | Bodies | Shape |
|---|---|---|
| `TUBES` | 7 | Straight tubes, inner radius 7.5, outer 10.5 (wall 3), from y = 0 to 300; one on the axis and six on a circle of radius 30 |
| `BAFFLES` | 5 | Two tube sheets (y = 0–6 and 294–300) and three segmental baffles (y = 71–77, 147–153, 223–229), all 6 thick, outer radius 45, seven tube holes of radius 10.5. The baffle cuts are at z = −20, +20, −20, so the windows alternate |
| `HOUSING` | 1 | Shell, inner radius 45, outer 51 (wall 6), with dished heads ending at y = −87 and 387. Two tube-side nozzles along z (inner radius 20, outer 26, at y = −30 and 330, ending at z = 100) and two shell-side nozzles along x (inner radius 15, outer 21, at y = 28, z = 23 ending at x = −100, and at y = 266, z = −23 ending at x = 100) |

## Domains and interfaces the mesher must find

- **Fluid `tube_side`:** the two heads plus the inside of the seven tubes, picked by
  a material point in a head (MD-06). Inlet and outlet: the two z-nozzles.
- **Fluid `shell_side`:** the space between the tube sheets, the shell and the tubes,
  around the baffles, picked by a material point between the tube sheets. Inlet and
  outlet: the two x-nozzles.
- The two fluids never touch; they meet only through the tube walls and the tube
  sheets (MD-08).
- **Solids:** 13 bodies in three parts. Domains may be one per body or one per part
  (MD-06); the material tag starts from the file (`Kunststoff`) and is set per part
  through the domain tools (MD-07, MD-15). The corrected file is expected to name
  different materials per part, which then become the expected tags.
- **Ports:** the nozzle ends are annular faces of the housing wall; the fluid
  openings inside them (radius 20 and 15) are not faces of the file and are closed
  by planar caps from their rim circles (MD-19).
- **Contacts to imprint (MD-10):** tube outer faces in the tube-sheet and baffle
  holes (coincident cylinders of radius 10.5, partly cut by the baffle windows),
  tube-sheet and baffle rims against the shell bore (radius 45), tube ends flush
  with the outer tube-sheet faces.
- **Thin parts and gaps (MD-11):** tube walls 3, sheets, baffles and shell 6; the
  gaps between neighbouring tubes are 9, and each shell-side nozzle bore starts
  7 from the nearest tube-sheet face.

## Analytic checks

Each tube wall has volume π(10.5² − 7.5²)·300 = 16,200π ≈ 50,893.80 mm³, and the
tube-side fluid inside each tube π·7.5²·300 = 16,875π ≈ 53,014.38 mm³. The
inner wall of each tube (area 2π·7.5·300) is a fluid–solid interface of
`tube_side`. These relations take the place of the synthetic tube-in-tube case
planned before this file arrived.
