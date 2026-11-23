# Heat exchanger (placeholder)

This folder will hold the maintainer's heat-exchanger example, announced on
2026-11-23 and not yet supplied. It is the M5b acceptance fixture for
multi-domain meshing ([MD-18](../../docs/requirements.md#multi-domain-meshes)):
several fluid passes and several solid materials, meshed for conjugate heat
transfer in FoehnCraft.

When the file arrives:

1. Store it byte for byte and record its provenance beside it (author, CAD
   tool and version, licence, SHA-256), as LIC-10 requires.
2. List its bodies here with the intended domain name, kind (`fluid`,
   `solid`, `porous`) and material tag (MD-06, MD-07), and the expected
   interfaces between them (MD-09).
3. If it is a STEP file, it is also the first STEP fixture (STP-01 to
   STP-05) and decides OPEN-14; record the file's schema and length unit.
4. Name the acceptance checks: one zone per fluid pass and solid material,
   conformal interfaces named on both sides that pass MD-14, and a
   `foehn mesh` read-back with one region per domain (FC-06, FC-12).

Until then, a synthetic tube-in-tube case (inner fluid, tube wall, annulus
fluid, outer shell of a second material) is the first multi-domain fixture,
because its volumes and interface areas are analytic.
