"""Run the application gallery entries and render their figures.

    python3 tools/gallery/app_gallery.py <kofl binary> <work dir> [entry ...]

The same machinery as poly_gallery.py (bash scripts printed on
docs/applications.md, MurmelCraft's headless export.image, Pillow); the
figures go to docs/figures/applications/. Entries: heat-exchanger (needs
gmsh: pip install gmsh), furnace-7 (needs gmsh and numpy) and airfoil
(the block mesher's three-element airfoil, from examples/gallery/block).
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import block_gallery as bg  # noqa: E402
import poly_gallery as pg  # noqa: E402

FIGURES = bg.ROOT / "docs" / "figures" / "applications"
bg.FIGURES = FIGURES
pg.FIGURES = FIGURES

pg.SCRIPTS["heat-exchanger"] = """
mkdir -p out/heat-exchanger && cd out/heat-exchanger
python3 ../../examples/gallery/applications/heat-exchanger/make_shell_side.py \\
  ../../test-fixtures/heat-exchanger/SHELL_TUBE_HX.step shell_side.stl
cp ../../examples/gallery/applications/heat-exchanger/*.json .
kofl call patch.detect_ports --input detect.json     # air_inlet and air_outlet are the two caps found
kofl call patch.rename --input rename-inlet.json
kofl call patch.rename --input rename-outlet.json    # -> patches.json
kofl call mesh.setup_duct --input setup.json         # -> partition.json
kofl generate-duct shell_side.stl patches.json partition.json polyhedral mesh --handoff vtu
"""
pg.OUTPUTS["heat-exchanger"] = ["heat-exchanger"]

PATCH, TYPE, std, FIT, zoom, clip = pg.PATCH, pg.TYPE, pg.std, pg.FIT, pg.zoom, pg.clip
pg.SHOTS["heat-exchanger"] = ("app_heat_exchanger", [
    ("heat-exchanger/mesh/foehncraft.vtu", [PATCH, std("front_right_top"), FIT]),
    ("heat-exchanger/mesh/foehncraft.vtu", [TYPE, std("top"), clip("z", 0.0, "flat"), FIT]),
    ("heat-exchanger/mesh/foehncraft.vtu", [TYPE, std("top"), clip("z", 0.0, "flat"), FIT, zoom(3.5, (-0.45, 0.2))]),
])
pg.SHOTS["furnace-7"] = ("app_furnace_7", pg.SHOTS["furnace-7"][1])

BLOCKS = bg.BLOCKS
NONORTH = {"tool": "display.color_by", "input": {"by": "quality", "metric": "non_orthogonality"}}


def airfoil(kofl, work, render_only):
    d = work / "multi-element"
    if not render_only:
        bg.run(kofl, work, "multi-element")
    shots = [
        [NONORTH, std("top"), clip("z", 0.0, "flat"), FIT, zoom(2.0)],
        [NONORTH, std("top"), clip("z", 0.0, "flat"), FIT, zoom(7.0, (-0.02, 0.0))],
    ]
    from PIL import Image

    pngs = []
    for k, steps in enumerate(shots):
        png = work / "shots" / f"app_airfoil-{k}.png"
        png.parent.mkdir(exist_ok=True)
        bg.shoot(d, "three-mesh", steps, png)
        pngs.append(Image.open(png).convert("RGB"))
    w = 800
    ims = [im.resize((w, round(im.height * w / im.width)), Image.LANCZOS) for im in pngs]
    sheet = Image.new("RGB", (w * len(ims), max(im.height for im in ims)), "white")
    for k, im in enumerate(ims):
        sheet.paste(im, (k * w, 0))
    FIGURES.mkdir(parents=True, exist_ok=True)
    sheet.save(FIGURES / "app_airfoil.jpg", quality=88)


def main():
    kofl = Path(sys.argv[1]).resolve()
    work = Path(sys.argv[2]).resolve()
    work.mkdir(parents=True, exist_ok=True)
    entries = sys.argv[3:] or ["heat-exchanger", "furnace-7", "airfoil"]
    render_only = os.environ.get("RENDER_ONLY") == "1"
    for entry in entries:
        if entry == "airfoil":
            airfoil(kofl, work, render_only)
        else:
            if not render_only:
                pg.run(kofl, work, entry)
            pg.render(entry, work)
        print(entry, "done", flush=True)


if __name__ == "__main__":
    main()
