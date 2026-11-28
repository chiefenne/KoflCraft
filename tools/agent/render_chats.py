"""Render the pictures of docs/chats.md from the meshes of the recorded chats.

    MURMEL=/path/to/murmel python3 tools/agent/render_chats.py <chats dir> [chat ...]

<chats dir> holds one folder per chat with the outputs of its calls
(pipe/mesh, naca/naca-tfi, manifold/mesh). The pictures are written by
MurmelCraft's headless export.image, as for the gallery
(tools/gallery/block_gallery.py), into docs/figures/chats/.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "gallery"))
import block_gallery as bg  # noqa: E402
import poly_gallery as pg  # noqa: E402

FIGURES = bg.ROOT / "docs" / "figures" / "chats"

# chat -> (figure, [(mesh file below the chat's folder, steps[, export.image options])])
SHOTS = {
    "pipe": ("chat_pipe_baffle", [
        ("mesh/foehncraft.vtu", [pg.PATCH, pg.hide("patch:walls"), bg.std("front_right_top"), bg.FIT]),
        ("mesh/foehncraft.vtu", [pg.TYPE, bg.std("front"), pg.clip("y", mode="flat"), bg.FIT,
                                 bg.zoom(2.2, (0.0, 0.2))]),
    ]),
    "naca": ("chat_naca_cgrid", [
        ("NACA_MESH/volume.vtu", [bg.BLOCKS, bg.std("top"), bg.FIT], bg.NO_LEGEND),
        ("NACA_MESH/volume.vtu", [bg.BLOCKS, bg.std("top"), bg.FIT, bg.zoom(30.0, (-0.012, 0.0))], bg.NO_LEGEND),
    ]),
    "manifold": ("chat_manifold_sweep", [
        ("mesh/foehncraft.vtu", [pg.TYPE, bg.std("top"), pg.clip("z", 40.0, mode="flat"), bg.FIT]),
        ("mesh/foehncraft.vtu", [pg.TYPE, bg.std("top"), pg.clip("z", 40.0, mode="flat"), bg.FIT,
                                 bg.zoom(3.5, (0.04, 0.62))]),
    ]),
}


def main():
    from PIL import Image

    chats = Path(sys.argv[1]).resolve()
    naca = (chats / "naca" / "mesh-name").read_text().strip() if (chats / "naca" / "mesh-name").exists() else "naca-tfi"
    FIGURES.mkdir(parents=True, exist_ok=True)
    for chat in sys.argv[2:] or list(SHOTS):
        figure, shots = SHOTS[chat]
        ims = []
        for k, (mesh, steps, *extra) in enumerate(shots):
            png = chats / chat / f"shot-{k}.png"
            bg.shoot_file(chats / chat / mesh.replace("NACA_MESH", naca), steps, png, extra[0] if extra else None)
            ims.append(Image.open(png).convert("RGB"))
        w = 800
        ims = [im.resize((w, round(im.height * w / im.width)), Image.LANCZOS) for im in ims]
        sheet = Image.new("RGB", (w * len(ims), max(im.height for im in ims)), "white")
        for k, im in enumerate(ims):
            sheet.paste(im, (k * w, 0))
        sheet.save(FIGURES / f"{figure}.jpg", quality=88)
        print(chat, "done")


if __name__ == "__main__":
    main()
