"""Render the README slideshow of KoflCraft meshes with MurmelCraft.

    python3 -I tools/readme-tour/tour.py MESH_DIR SHOTS_DIR

MESH_DIR holds the mesh made by make-tour.sh (manifold/foehncraft.vtu). For each slide this script
writes a `murmel run --script` file (open the mesh, set the colouring, the
clip plane and the camera, then `export.image` at 2560 x 1600, scale 2) and
runs it headless; the PNGs land in SHOTS_DIR as NN-name.png. Set MURMEL to
the murmel binary (default: murmel on the PATH).
"""
import json
import os
import pathlib
import subprocess
import sys

MURMEL = os.environ.get("MURMEL", "murmel")


def cell_types():
    return {"tool": "display.color_by", "input": {"by": "cell_type"}}


def clip(axis, mode):
    return {"tool": "clip.set", "input": {"axis": axis, "offset": 0.0, "mode": mode, "enabled": True}}


def camera(target, distance, yaw, pitch):
    return {"tool": "view.set_camera",
            "input": {"target": target, "distance": distance, "yaw_deg": yaw, "pitch_deg": pitch,
                      "projection": "perspective"}}


# (file name, mesh folder, steps, export.image extras)
SLIDES = [
    ("01-manifold-patches", "manifold",
     [{"tool": "display.color_by", "input": {"by": "patch"}}, camera([150.0, 200.0, 20.0], 900.0, -45.0, 35.0)], {}),
    ("02-manifold-cell-types", "manifold",
     [cell_types(), camera([150.0, 200.0, 20.0], 900.0, -45.0, 35.0), clip("z", "crinkle")], {}),
    ("03-manifold-quality", "manifold",
     [{"tool": "display.color_by", "input": {"by": "quality", "metric": "non_orthogonality"}},
      camera([150.0, 200.0, 20.0], 900.0, -45.0, 35.0), clip("z", "flat")], {"theme": "light"}),
]


def main():
    meshes, shots = pathlib.Path(sys.argv[1]).resolve(), pathlib.Path(sys.argv[2]).resolve()
    only = set(sys.argv[3:])
    shots.mkdir(parents=True, exist_ok=True)
    for name, mesh, steps, extra in SLIDES:
        if only and name not in only:
            continue
        image = {"path": str(shots / f"{name}.png"), "width": 2560, "height": 1600, "scale": 2}
        image.update(extra)
        script = ([{"tool": "mesh.open", "input": {"path": str(meshes / mesh / "foehncraft.vtu")}}]
                  + steps + [{"tool": "export.image", "input": image}])
        path = shots / f"{name}.json"
        path.write_text(json.dumps(script, indent=1))
        run = subprocess.run([MURMEL, "run", "--script", str(path)], capture_output=True, text=True, check=False)
        calls = json.loads(run.stdout)["calls"] if run.stdout else []
        failed = [(s["tool"], c["diagnostics"]) for s, c in zip(script, calls) if c["exit_code"] or c["diagnostics"]]
        if run.returncode or failed or len(calls) != len(script):
            sys.exit(f"{name}: {failed or run.stderr}")
        print(f"{name}.png")


if __name__ == "__main__":
    main()
