"""Run the polyhedral gallery entries and render their figures.

    python3 tools/gallery/poly_gallery.py <kofl binary> <work dir> [entry ...]

Each entry is the shell script printed on docs/gallery/polyhedral.md. It
runs with bash in <work dir>, which holds links to the repository's
test-fixtures/ and examples/ and the outputs in out/, with `kofl` on the
PATH. The figures are rendered into docs/figures/gallery/ with
MurmelCraft's headless export.image (MURMEL: the murmel binary, default
murmel on the PATH) and Pillow, as for block_gallery.py. The 7-tube
furnace needs the gmsh and numpy Python modules.
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import block_gallery as bg  # noqa: E402

ROOT = bg.ROOT
FIGURES = bg.FIGURES

# entry -> the commands of its gallery entry, run from <work dir>
SCRIPTS = {
    "box": """
kofl generate test-fixtures/box-4.json out/box.json
kofl export out/box.json vtu out/box.vtu
""",
    "channel": """
kofl generate-channel test-fixtures/straight-channel.json out/channel --handoff vtu
""",
    "sweep": """
kofl generate-sweep test-fixtures/analytic/pipe.json out/sweep-pipe --handoff vtu
kofl generate-sweep test-fixtures/analytic/block.json out/sweep-block --handoff vtu
""",
    "convex-stl": """
kofl generate-stl test-fixtures/unit-box.stl test-fixtures/stl-size-field.json out/stl.json
kofl export out/stl.json vtu out/stl.vtu
""",
    "elbow": """
E=test-fixtures/elbow
kofl generate-duct $E/elbow.stl $E/patches.json $E/partition.json bulk out/elbow-bulk --handoff vtu
kofl generate-duct $E/elbow.stl $E/patches.json $E/partition.json layered out/elbow-layered --handoff vtu
kofl generate-duct $E/elbow.stl $E/patches.json $E/partition.json polyhedral out/elbow-polyhedral --handoff vtu
""",
    "manifold": """
M=test-fixtures/manifold
kofl generate-duct $M/manifold.stl $M/patches.json $M/partition.json polyhedral out/manifold --handoff vtu
""",
    "bare-stl": """
mkdir -p out/bare && cp test-fixtures/manifold/manifold.stl out/bare/ && cd out/bare
echo '{"surface": "manifold.stl", "output": "ports.json"}' > detect.json
kofl call patch.detect_ports --input detect.json
echo '{"surface": "manifold.stl", "patches": "ports.json", "output": "partition.json", "size_factor": 2}' > setup.json
kofl call mesh.setup_duct --input setup.json
kofl generate-duct manifold.stl ports.json partition.json polyhedral mesh --handoff vtu
""",
    "furnace": """
F=test-fixtures/furnace
kofl generate-duct $F/furnace.stl $F/furnace.patches.json $F/partition.json polyhedral out/furnace --handoff vtu
""",
    "furnace-variants": """
F=test-fixtures/furnace V=examples/gallery/polyhedral/furnace-variants
kofl generate-duct $F/furnace.stl $F/furnace.patches.json $V/partition-yield.json polyhedral out/furnace-yield --handoff vtu
kofl generate-duct $F/furnace.stl $F/furnace.patches.json $V/partition-hex-core.json polyhedral out/furnace-hex-core --handoff vtu
""",
    "furnace-7": """
mkdir -p out/furnace-7-input && cd out/furnace-7-input
sed 's/"Mesh.MeshSizeFromCurvature", 24/"Mesh.MeshSizeFromCurvature", 48/' \\
  ../../test-fixtures/furnace/make_furnace.py > make_furnace7.py
python3 make_furnace7.py furnace-7.stl 7
python3 ../../examples/gallery/polyhedral/furnace-7/set_sha256.py furnace-7.stl \\
  ../../examples/gallery/polyhedral/furnace-7/partition-7.json partition-7.json
kofl generate-duct furnace-7.stl furnace-7.patches.json partition-7.json polyhedral ../furnace-7 --handoff vtu
""",
    "pipe-shells": """
mkdir -p out/pipe-shells && cp examples/gallery/polyhedral/pipe-shells/* out/pipe-shells/ && cd out/pipe-shells
python3 make_pipe.py pipe.stl
kofl call patch.propose --input propose.json      # opening_1, opening_2, walls
kofl call patch.rename --input name-in.json       # opening_1 -> inlet
kofl call patch.rename --input name-out.json      # opening_2 -> outlet
kofl call patch.set_role --input role-in.json
kofl call patch.set_role --input role-out.json    # -> patches.json
kofl call bc.skeleton --input skeleton.json       # -> case.toml for FoehnCraft
kofl call mesh.setup_duct --input setup.json      # -> partition.json
kofl call shell.add --input baffle.json           # -> partition-1.json
kofl call shell.add --input plate.json            # -> partition-2.json
kofl call grading.add --input probe.json          # -> partition-3.json
kofl generate-duct pipe.stl patches.json partition-3.json polyhedral mesh --handoff vtu
""",
}


def color(by, **k):
    return {"tool": "display.color_by", "input": {"by": by, **k}}


def clip(axis, offset=0.0, mode="crinkle"):
    return {"tool": "clip.set", "input": {"axis": axis, "offset": offset, "mode": mode, "enabled": True}}


def hide(*items):
    return {"tool": "display.visibility", "input": {"action": "hide", "items": list(items)}}


std, FIT, zoom, NO_LEGEND = bg.std, bg.FIT, bg.zoom, bg.NO_LEGEND
PATCH, TYPE = color("patch"), color("cell_type")
FACES = color("quality", metric="face_count")

# entry -> (figure, [(mesh file in out/, steps[, export.image options])])
SHOTS = {
    "box": ("poly_box", [
        ("box.vtu", [FACES, std("front_right_top"), FIT]),
    ]),
    "channel": ("poly_channel", [
        ("channel/foehncraft.vtu", [PATCH, std("front_right_top"), FIT]),
        ("channel/foehncraft.vtu", [TYPE, std("front_right_top"), clip("y"), FIT]),
    ]),
    "sweep": ("poly_sweep", [
        ("sweep-pipe/foehncraft.vtu", [TYPE, std("left"), clip("x", mode="flat"), FIT, zoom(1.8, (0.0, 0.45))]),
        ("sweep-block/foehncraft.vtu", [TYPE, std("front_right_top"), clip("y"), FIT]),
    ]),
    "convex-stl": ("poly_convex_stl", [
        ("stl.vtu", [FACES, std("front_right_top"), clip("z"), FIT]),
        ("stl.vtu", [FACES, std("top"), clip("z", 0.25 - 0.5, "flat"), FIT]),
    ]),
    "elbow": ("poly_elbow", [
        ("elbow-bulk/foehncraft.vtu", [TYPE, std("top"), clip("z", mode="flat"), FIT]),
        ("elbow-layered/foehncraft.vtu", [TYPE, std("top"), clip("z", mode="flat"), FIT]),
        ("elbow-polyhedral/foehncraft.vtu", [TYPE, std("top"), clip("z", mode="flat"), FIT]),
    ]),
    "manifold": ("poly_manifold", [
        ("manifold/foehncraft.vtu", [PATCH, std("front_right_top"), FIT]),
        ("manifold/foehncraft.vtu", [TYPE, std("front_right_top"), clip("z"), FIT]),
    ]),
    "bare-stl": ("poly_bare_stl", [
        ("bare/mesh/foehncraft.vtu", [PATCH, std("front_right_top"), FIT]),
        ("bare/mesh/foehncraft.vtu", [TYPE, std("front_right_top"), clip("z"), FIT]),
    ]),
    "furnace": ("poly_furnace", [
        ("furnace/foehncraft.vtu", [PATCH, std("front_right_top"), clip("y"), FIT]),
        ("furnace/foehncraft.vtu", [TYPE, std("front"), clip("y"), FIT, zoom(6.0, (-0.5, 0.62))]),
    ]),
    "furnace-variants": ("poly_furnace_variants", [
        ("furnace/foehncraft.vtu", [TYPE, std("top"), clip("z", 0.85, "flat"), FIT, zoom(3.0, (-0.55, -0.75))]),
        ("furnace-yield/foehncraft.vtu", [TYPE, std("top"), clip("z", 0.85, "flat"), FIT, zoom(3.0, (-0.55, -0.75))]),
        ("furnace-hex-core/foehncraft.vtu", [TYPE, std("front"), clip("y"), FIT]),
    ]),
    "furnace-7": ("poly_furnace_7", [
        ("furnace-7/foehncraft.vtu", [PATCH, std("front_right_top"), clip("y"), FIT]),
        ("furnace-7/foehncraft.vtu", [TYPE, std("front"), clip("y"), FIT, zoom(3.0, (0.0, 0.62))]),
    ]),
    "pipe-shells": ("poly_pipe_shells", [
        ("pipe-shells/mesh/foehncraft.vtu", [PATCH, hide("patch:walls"), std("front_right_top"), FIT]),
        ("pipe-shells/mesh/foehncraft.vtu", [TYPE, std("front"), clip("y", mode="flat"), FIT, zoom(3.0, (0.0, -0.4))]),
    ]),
}
# entries whose pictures use another entry's outputs
NEEDS = {"furnace-variants": "furnace"}
# what an entry writes into out/; removed before it runs again
OUTPUTS = {
    "box": ["box.json", "box.vtu"],
    "channel": ["channel"],
    "sweep": ["sweep-pipe", "sweep-block"],
    "convex-stl": ["stl.json", "stl.vtu"],
    "elbow": ["elbow-bulk", "elbow-layered", "elbow-polyhedral"],
    "manifold": ["manifold"],
    "bare-stl": ["bare"],
    "furnace": ["furnace"],
    "furnace-variants": ["furnace-yield", "furnace-hex-core"],
    "furnace-7": ["furnace-7-input", "furnace-7"],
    "pipe-shells": ["pipe-shells"],
}


def run(kofl, work, entry):
    for link in ("test-fixtures", "examples"):
        if not (work / link).exists():
            (work / link).symlink_to(ROOT / link)
    bin_dir = work / "bin"
    bin_dir.mkdir(exist_ok=True)
    if not (bin_dir / "kofl").exists():
        (bin_dir / "kofl").symlink_to(kofl)
    (work / "out").mkdir(exist_ok=True)
    for name in OUTPUTS[entry]:
        path = work / "out" / name
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()
    env = dict(os.environ, PATH=f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    out = subprocess.run(["bash", "-e", "-c", SCRIPTS[entry]], cwd=work, env=env, capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"{entry} failed\n{out.stdout[-3000:]}\n{out.stderr[-3000:]}")


def render(entry, work):
    from PIL import Image

    figure, shots = SHOTS[entry]
    pngs = []
    for k, (mesh, steps, *extra) in enumerate(shots):
        png = work / "shots" / f"{figure}-{k}.png"
        png.parent.mkdir(exist_ok=True)
        bg.shoot_file(work / "out" / mesh, steps, png, extra[0] if extra else None)
        pngs.append(Image.open(png).convert("RGB"))
    w = 800
    ims = [im.resize((w, round(im.height * w / im.width)), Image.LANCZOS) for im in pngs]
    sheet = Image.new("RGB", (w * len(ims), max(im.height for im in ims)), "white")
    for k, im in enumerate(ims):
        sheet.paste(im, (k * w, 0))
    FIGURES.mkdir(parents=True, exist_ok=True)
    sheet.save(FIGURES / f"{figure}.jpg", quality=88)


def main():
    kofl = Path(sys.argv[1]).resolve()
    work = Path(sys.argv[2]).resolve()
    work.mkdir(parents=True, exist_ok=True)
    entries = sys.argv[3:] or list(SCRIPTS)
    render_only = os.environ.get("RENDER_ONLY") == "1"
    for entry in entries:
        if not render_only:
            if entry in NEEDS and not (work / "out" / NEEDS[entry]).exists():
                run(kofl, work, NEEDS[entry])
            run(kofl, work, entry)
        render(entry, work)
        print(entry, "done", flush=True)


if __name__ == "__main__":
    main()
