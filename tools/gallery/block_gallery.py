"""Run the block-mesher gallery entries and render their figures.

    python3 tools/gallery/block_gallery.py <kofl binary> <work dir> [entry ...]

Each entry copies its inputs from examples/gallery/block/<entry>/ into
<work dir>/<entry>/, runs the kofl calls listed below in that order and
renders the figures into docs/figures/gallery/ with MurmelCraft's
headless export.image (MURMEL: the murmel binary, default murmel on the
PATH; a Vulkan driver, Mesa lavapipe on a headless Linux machine) and
Pillow (see tools/gallery/README.md).
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path

import os

ROOT = Path(__file__).resolve().parents[2]
INPUTS = ROOT / "examples" / "gallery" / "block"
FIGURES = ROOT / "docs" / "figures" / "gallery"

# entry -> list of (tool, input file)
CALLS = {
    "channel": [("block.template", "template.json"), ("block.grid", "grid.json")],
    "pipe": [("block.template", "template.json"), ("block.grid", "grid.json")],
    "elbow": [("block.template", "template.json"), ("block.grid", "grid.json")],
    "naca": [
        ("block.template", "template.json"),
        ("block.grid", "grid-tfi.json"),
        ("block.grid", "grid.json"),
    ],
    "pipe-split": [("block.split_block", "split.json"), ("block.grid", "grid.json")],
    "multi-element": [("block.template", "template.json"), ("block.grid", "grid.json")],
    "triangle-split": [
        ("block.grid", "grid-before.json"),
        ("block.optimize_layout", "optimize.json"),
        ("block.grid", "grid-after.json"),
    ],
    "tri-center": [
        ("block.template", "template-triangle.json"),
        ("block.grid", "grid-triangle.json"),
        ("block.template", "template-tetrahedron.json"),
        ("block.grid", "grid-tetrahedron.json"),
    ],
}
# entries whose inputs come from another entry's outputs
NEEDS = {"pipe-split": ("pipe", ["pipe.json"])}

MURMEL = os.environ.get("MURMEL", "murmel")
SIZE = (1600, 1000)


def std(view):
    return {"tool": "view.standard", "input": {"view": view}}


FIT = {"tool": "view.fit", "input": {}}
BLOCKS = {"tool": "display.color_by", "input": {"by": "category", "field": "BlockId"}}


def explode(f):
    return {"tool": "display.explode", "input": {"factor": f}}


def clip(axis, mode="crinkle"):
    return {"tool": "clip.set", "input": {"axis": axis, "offset": 0.0, "mode": mode, "enabled": True}}


def zoom(factor, at=(0.0, 0.0)):
    """Zoom in by factor about a point of the view (normalised device coordinates)."""
    return {"tool": "view.zoom", "input": {"factor": 1.0 / factor, "cursor_ndc": list(at)}}


NO_LEGEND = {"legend": False}


# entry -> list of (mesh folder in the entry's work dir, steps); the shots
# of an entry are put side by side into docs/figures/gallery/<figure>.png
SHOTS = {
    "channel": ("block_channel", [
        ("channel-mesh", [BLOCKS, std("front_right"), FIT]),
        ("channel-mesh", [BLOCKS, std("left"), FIT]),
    ]),
    "pipe": ("block_pipe_ogrid", [
        ("pipe-mesh", [BLOCKS, std("left"), FIT]),
        ("pipe-mesh", [BLOCKS, explode(0.35), std("front_right"), FIT]),
    ]),
    "elbow": ("block_elbow", [
        ("elbow-mesh", [BLOCKS, std("top"), clip("z", "flat"), FIT], NO_LEGEND),
        ("elbow-mesh", [BLOCKS, explode(0.25), std("front_top"), FIT], NO_LEGEND),
    ]),
    "naca": ("block_naca_cgrid", [
        ("naca-mesh", [BLOCKS, std("top"), FIT]),
        ("naca-mesh", [BLOCKS, std("top"), FIT, zoom(8.0, (-0.05, 0.0))]),
    ]),
    "pipe-split": ("block_pipe_split", [("pipe-split-mesh", [BLOCKS, explode(0.35), std("front_right"), FIT])]),
    "multi-element": ("block_multi_element", [
        ("three-mesh", [BLOCKS, std("top"), FIT, zoom(7.0, (-0.02, 0.0))], NO_LEGEND),
        ("three-mesh", [BLOCKS, explode(0.6), std("top"), FIT, zoom(3.0)], NO_LEGEND),
    ]),
    "triangle-split": ("block_triangle_split", [
        ("wedge-before", [BLOCKS, std("top"), FIT]),
        ("wedge-after", [BLOCKS, std("top"), FIT]),
        ("wedge-after", [BLOCKS, explode(0.25), std("top"), FIT]),
    ]),
    "tri-center": ("block_tri_center", [
        ("triangle-mesh", [BLOCKS, explode(0.25), std("top"), FIT]),
        ("tetrahedron-mesh", [BLOCKS, explode(0.5), std("front_right"), FIT]),
    ]),
}


def run(kofl, work, entry):
    d = work / entry
    if d.exists():
        shutil.rmtree(d)
    shutil.copytree(INPUTS / entry, d)
    if entry in NEEDS:
        src, files = NEEDS[entry]
        for f in files:
            shutil.copy(work / src / f, d / f)
    results = []
    for tool, inp in CALLS[entry]:
        out = subprocess.run(
            [kofl, "call", tool, "--input", inp], cwd=d, capture_output=True, text=True
        )
        if out.returncode != 0:
            sys.exit(f"{entry}: {tool} {inp} failed\n{out.stdout}\n{out.stderr}")
        results.append(json.loads(out.stdout)["result"])
    (d / "results.json").write_text(json.dumps(results, indent=1))
    return d, results


def shoot(d, folder, steps, png, extra=None):
    """One picture of <d>/<folder>/volume.vtu by MurmelCraft's headless export.image."""
    shoot_file(d / folder / "volume.vtu", steps, png, extra)


def shoot_file(mesh, steps, png, extra=None):
    """One picture of the mesh file <mesh> by MurmelCraft's headless export.image."""
    image = {"path": str(png), "width": SIZE[0], "height": SIZE[1], "scale": 2, "theme": "light"}
    image.update(extra or {})
    script = ([{"tool": "mesh.open", "input": {"path": str(mesh)}}]
              + steps + [{"tool": "export.image", "input": image}])
    path = png.with_suffix(".json")
    path.write_text(json.dumps(script, indent=1))
    run = subprocess.run([MURMEL, "run", "--script", str(path)], capture_output=True, text=True)
    calls = json.loads(run.stdout)["calls"] if run.stdout else []
    failed = [(s["tool"], c["diagnostics"]) for s, c in zip(script, calls) if c["exit_code"] or c["diagnostics"]]
    if run.returncode or failed or len(calls) != len(script):
        sys.exit(f"{png.name}: {failed or run.stderr}")


def render(entry, d):
    from PIL import Image

    figure, shots = SHOTS[entry]
    pngs = []
    for k, (folder, steps, *extra) in enumerate(shots):
        png = d / f"shot-{k}.png"
        shoot(d, folder, steps, png, extra[0] if extra else None)
        pngs.append(Image.open(png).convert("RGB"))
    # side by side, each scaled to 800 pixels wide
    w = 800
    ims = [im.resize((w, round(im.height * w / im.width)), Image.LANCZOS) for im in pngs]
    sheet = Image.new("RGB", (w * len(ims), max(im.height for im in ims)), "white")
    for k, im in enumerate(ims):
        sheet.paste(im, (k * w, 0))
    FIGURES.mkdir(parents=True, exist_ok=True)
    sheet.save(FIGURES / f"{figure}.jpg", quality=88)


def main():
    kofl = str(Path(sys.argv[1]).resolve())
    work = Path(sys.argv[2]).resolve()
    work.mkdir(parents=True, exist_ok=True)
    entries = sys.argv[3:] or list(CALLS)
    for entry in entries:
        if entry in NEEDS and not (work / NEEDS[entry][0]).exists():
            run(kofl, work, NEEDS[entry][0])
        d, results = run(kofl, work, entry)
        render(entry, d)
        print(entry, "done")


if __name__ == "__main__":
    main()
