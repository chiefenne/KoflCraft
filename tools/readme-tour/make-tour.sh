#!/usr/bin/env bash
# Make the README slideshow (assets/readme/tour.webp and the stills in
# assets/readme/tour/) from the KoflCraft Manifold mesh, pictured by MurmelCraft's
# headless export.image. Needs a release build of mesh-cli here, a release
# build of MurmelCraft's murmel (MURMEL, default murmel on the PATH), a
# Vulkan driver (Mesa lavapipe is enough) and Python 3 with Pillow.
#
#   tools/readme-tour/make-tour.sh [WORK_DIR]
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work=$(realpath -m "${1:-$repo/target/readme-tour}")
kofl=$repo/target/release/mesh-cli
mkdir -p "$work/meshes"
cd "$repo"

mesh() { # folder, then the mesh-cli arguments before the output folder
  local out=$work/meshes/$1; shift
  [ -f "$out/foehncraft.vtu" ] || "$kofl" "$@" "$out" > "$out.json"
}
m=test-fixtures/manifold
mesh manifold generate-duct $m/manifold.stl $m/patches.json $m/partition.json polyhedral

python3 -I "$here/tour.py" "$work/meshes" "$work/shots"
python3 -I "$here/make_webp.py" "$work/shots" "$repo/assets/readme"
