# Shell-and-tube heat exchanger (STEP)

The maintainer's test case for the STEP reader (section 4b, STP-01 to STP-13)
and the multi-domain mesher (MD-18), shipped in two variants (maintainer
decision 2026-11-23). Both files are stored byte for byte.

| File | Variant | What it tests |
|---|---|---|
| [`SHELL_TUBE_HX.step`](SHELL_TUBE_HX.step) | **Open:** the nozzles end in annular wall faces and the fluid openings are not closed | Automatic port caps (MD-19); the reader on solids only |
| [`SHELL_TUBE_HX_closed_fluid_domains.step`](SHELL_TUBE_HX_closed_fluid_domains.step) | **Closed:** four named planar cap faces close the nozzle openings | The reader including named surface models (STP-01, STP-04), the fluid domains and named ports before MD-19 is built |

Both variants must give the same domains and port patches. Version 2
(2026-11-23) replaces the version 1 file of the same day, which had Fusion's
default body names (`Körper1` to `Körper7`, repeated across parts); its
geometry is unchanged (same entity counts).

## Provenance (LIC-10)

| Item | Open | Closed |
|---|---|---|
| Supplied by | the maintainer (chiefenne), 2026-11-23 | same |
| Written | 2026-11-23 09:53 +01:00 | 2026-11-23 09:59 +01:00 |
| Created with | Autodesk Fusion (`originating_system` Autodesk Translation Framework v15.17.0.0), written with ST-Developer v20.2 | same |
| Schema | AP242 (`AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 2 1 4 }`), ISO 10303-21 clear text | same |
| Size | 246,727 bytes | 250,462 bytes |
| SHA-256 | `538b0c4e2a4895d0eaee62305f1701fba692abcf31349b65541152e6203fe522` | `327ecb1b32edf324489f53fb351a16870915609c60754597e1a2a5a975450c2f` |
| Licence | the maintainer's own model, supplied as a KoflCraft test fixture | same |

## What the files contain

Inspected on 2026-11-23 with a throw-away entity count, not with a KoflCraft
reader. These are the values the reader must reproduce.

- **Units:** length in millimetres (`SI_UNIT(.MILLI.,.METRE.)`), angles in radians,
  mass in kilograms. Distance accuracy 0.01 mm.
- **Assembly:** one root product `SHELL_TUBE_HX` with three parts, each used once
  with an identity placement: `TUBES:1`, `HOUSING:1`, `BAFFLES:1`.
- **Solids:** 13 manifold solid B-reps, 117 faces, all with unique names:

| Part | Bodies |
|---|---|
| `TUBES` | `tube_center`, `tube_mid_left`, `tube_mid_right`, `tube_top_left`, `tube_top_right`, `tube_bottom_left`, `tube_bottom_right` |
| `HOUSING` | `housing` |
| `BAFFLES` | `tube_sheet_inlet` (y = 0–6), `tube_sheet_outlet` (y = 294–300), `baffle_1` (y = 71–77), `baffle_2` (y = 147–153), `baffle_3` (y = 223–229) |

- **Cap faces (closed file only):** four `SHELL_BASED_SURFACE_MODEL` entities at the
  root of the assembly, outside the three parts. Each holds one `OPEN_SHELL` with one
  planar disc face:

| Name | Nozzle | Plane, centre | Radius | Normal in the file |
|---|---|---|---|---|
| `water_inlet` | tube side, at y = −30 | z = 100, (0, −30, 100) | 20 | −z (into the fluid) |
| `water_outlet` | tube side, at y = 330 | z = 100, (0, 330, 100) | 20 | +z (out of the fluid) |
| `air_inlet` | shell side, at y = 28 | x = −100, (−100, 28, 23) | 15 | −x (out of the fluid) |
| `air_outlet` | shell side, at y = 266 | x = 100, (100, 266, −23) | 15 | −x (into the fluid) |

  The cap normals are not consistent, so a reader orients a cap by its domain,
  not by the file. `water_inlet` also carries a red style on its face.

- **Materials:** in **both** files every part carries the AP214/AP242 material
  property `Kunststoff` (plastic) with density 1290 kg/m³. The materials the
  maintainer assigned in Fusion did not reach the material properties of the
  closed file. Only the appearances changed:

| Part | Appearance, open file | Appearance, closed file |
|---|---|---|
| `TUBES` | `Kunststoff - glänzend (Grau)` | `Farbe - Emaille glänzend (Rot)` |
| `HOUSING` | `Farbe - Emaille glänzend (Rot)` | `Stahl - satiniert` |
| `BAFFLES` | `Farbe - Emaille glänzend (Gelb)` | `Farbe - Emaille glänzend (Gelb)` |
| cap faces | — | `Kunststoff - glänzend (Grau)` |

  So the reader must report `Kunststoff` as the file's material for all parts
  (STP-13). The intended material tags of the solids are still to be confirmed
  by the maintainer and are recorded here when they are. The fluids have no
  material in the file; the cap names suggest `water` on the tube side and `air`
  on the shell side.

- **Surfaces:** planes, cylinders, cones, tori and B-spline surfaces (41, 64, 4, 4, 4
  in the open file; the closed file adds the 4 cap planes). **Curves:** lines,
  circles and B-spline curves (95, 136, 16; the closed file adds 4 cap circles).
- **No fluid bodies.** The two fluids are the empty space of the model.

## Geometry (y is the axis of the exchanger, mm)

| Part | Shape |
|---|---|
| `TUBES` | Straight tubes, inner radius 7.5, outer 10.5 (wall 3), from y = 0 to 300; `tube_center` on the axis, six on a circle of radius 30 |
| `BAFFLES` | Two tube sheets and three segmental baffles, all 6 thick, outer radius 45, seven tube holes of radius 10.5. The cuts of `baffle_1`, `baffle_2`, `baffle_3` are at z = −20, +20, −20, so the windows alternate |
| `HOUSING` | Shell, inner radius 45, outer 51 (wall 6), with dished heads ending at y = −87 and 387. Two tube-side nozzles along z (inner radius 20, outer 26, at y = −30 and 330, ending at z = 100) and two shell-side nozzles along x (inner radius 15, outer 21, at y = 28, z = 23 ending at x = −100, and at y = 266, z = −23 ending at x = 100) |

## Expected domains, ports and interfaces

| Domain | Kind | Picked by | Ports (patch names) |
|---|---|---|---|
| `tube_side` | fluid | material point in a head, or the region bounded by the water caps | `water_inlet`, `water_outlet` |
| `shell_side` | fluid | material point between the tube sheets, or the region bounded by the air caps | `air_inlet`, `air_outlet` |
| one per solid body, named as in the file (or one per part, MD-06) | solid | body | — |

- In the closed file the port names come from the cap faces. In the open file the
  caps are proposed from the rim circles of the nozzle openings (MD-19) and given
  the same names, so both variants produce the same patches.
- The two fluids never touch; they meet only through the tube walls and the tube
  sheets (MD-08).
- **Contacts to imprint (MD-10):** tube outer faces in the tube-sheet and baffle
  holes (coincident cylinders of radius 10.5, partly cut by the baffle windows),
  tube-sheet and baffle rims against the shell bore (radius 45), tube ends flush
  with the outer tube-sheet faces.
- **Thin parts and gaps (MD-11):** tube walls 3, sheets, baffles and shell 6; the
  gaps between neighbouring tubes are 9, and each shell-side nozzle bore starts 7
  from the nearest tube-sheet face.

## Analytic checks

Each tube wall has volume π(10.5² − 7.5²)·300 = 16,200π ≈ 50,893.80 mm³, and the
tube-side fluid inside each tube π·7.5²·300 = 16,875π ≈ 53,014.38 mm³. The
inner wall of each tube (area 2π·7.5·300) is a fluid–solid interface of
`tube_side`. The cap areas are 400π (water) and 225π (air) mm².
