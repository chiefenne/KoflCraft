#!/usr/bin/env bash
# Make the README slideshow (assets/readme/tour.webp and the stills in
# assets/readme/tour/) from KoflCraft meshes, pictured by MurmelCraft's
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
m=test-fixtures/manifold e=test-fixtures/elbow
mesh manifold generate-duct $m/manifold.stl $m/patches.json $m/partition.json polyhedral
mesh elbow generate-duct $e/elbow.stl $e/patches.json $e/partition.json polyhedral
mesh sweep generate-sweep test-fixtures/analytic/block.json
mesh pipe generate-sweep test-fixtures/analytic/pipe.json
mesh channel generate-channel test-fixtures/straight-channel.json

python3 -I "$here/tour.py" "$work/meshes" "$work/shots"
python3 -I "$here/make_webp.py" "$work/shots" "$repo/assets/readme"
