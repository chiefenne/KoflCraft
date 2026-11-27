"""Render the six hero pictures of the README.

    MURMEL=/path/to/murmel python3 tools/readme-hero/make_heroes.py <work dir>

<work dir> holds the meshes of the gallery runs (see tools/gallery/README.md):
  <work>/gallery/out/  poly_gallery.py manifold and app_gallery.py heat-exchanger furnace-7
  <work>/gblock/       block_gallery.py elbow naca multi-element
Each picture is MurmelCraft's headless export.image (1600 x 1600, transparent,
no legend or axes), cropped to the mesh, centred on white and written as a
600 x 600 JPEG to assets/readme/hero_<name>.jpg.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "gallery"))
import block_gallery as bg  # noqa: E402
import poly_gallery as pg  # noqa: E402
from PIL import Image  # noqa: E402

std, FIT, zoom, clip, BLOCKS, explode = bg.std, bg.FIT, bg.zoom, pg.clip, bg.BLOCKS, bg.explode
TYPE = pg.TYPE
NONORTH = {"tool": "display.color_by", "input": {"by": "quality", "metric": "non_orthogonality"}}
IMAGE = {"width": 1600, "height": 1600, "scale": 2, "legend": False, "axes": False, "transparent": True}


def shots(work):
    g, b = work / "gallery" / "out", work / "gblock"
    return {
        "furnace_tube": (g / "furnace-7/foehncraft.vtu",
                         [TYPE, std("front"), clip("y"), FIT, zoom(6.0, (0.0, 0.62))]),
        "naca_cgrid": (b / "naca/naca-mesh/volume.vtu", [BLOCKS, std("top"), FIT, zoom(8.0, (-0.05, 0.0))]),
        "manifold_bend": (g / "manifold/foehncraft.vtu",
                          [NONORTH, std("top"), clip("z", 0.0, "flat"), FIT, zoom(3.0, (-0.35, 0.0))]),
        "elbow_blocks": (b / "elbow/elbow-mesh/volume.vtu", [BLOCKS, explode(0.25), std("front_top"), FIT, zoom(1.1)]),
        "heat_exchanger": (g / "heat-exchanger/mesh/foehncraft.vtu",
                           [TYPE, std("front"), clip("y", 0.0, "flat"), FIT, zoom(0.6)]),
        "three_element": (b / "multi-element/three-mesh/volume.vtu",
                          [BLOCKS, std("top"), FIT, zoom(7.0, (-0.02, 0.0))]),
    }


def main():
    work = Path(sys.argv[1]).resolve()
    tmp = work / "hero"
    tmp.mkdir(parents=True, exist_ok=True)
    for name, (mesh, steps) in shots(work).items():
        png = tmp / f"{name}.png"
        bg.shoot_file(mesh, steps, png, dict(IMAGE))
        im = Image.open(png).convert("RGBA")
        im = im.crop(im.getchannel("A").point(lambda a: 255 if a > 8 else 0).getbbox())
        side, margin = 1200, 0.06
        f = side * (1 - 2 * margin) / max(im.width, im.height)
        im = im.resize((round(im.width * f), round(im.height * f)), Image.LANCZOS)
        sheet = Image.new("RGBA", (side, side), "white")
        sheet.alpha_composite(im, ((side - im.width) // 2, (side - im.height) // 2))
        out = ROOT / "assets" / "readme" / f"hero_{name}.jpg"
        sheet.convert("RGB").resize((600, 600), Image.LANCZOS).save(out, quality=88)
        print(name, "done", flush=True)


if __name__ == "__main__":
    main()
