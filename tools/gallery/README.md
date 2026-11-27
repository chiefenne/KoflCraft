# Gallery pictures

The pictures of [docs/gallery](../../docs/gallery.md) are KoflCraft meshes
pictured by [MurmelCraft](https://github.com/chiefenne/MurmelCraft)'s
headless `export.image`, made again by one script per gallery page.

1. `cargo build --release -p mesh-cli` here, and in a MurmelCraft checkout
   `cargo build --release -p murmel-cli` (with `display.explode`, MurmelCraft
   PR #124 or later).
2. A Vulkan driver for the offscreen render (`apt install
   mesa-vulkan-drivers` gives lavapipe) and Python 3 with Pillow.
3. `MURMEL=/path/to/murmel python3 tools/gallery/block_gallery.py target/release/kofl target/gallery [entry ...]`

[block_gallery.py](block_gallery.py) copies each entry's inputs from
`examples/gallery/block/<entry>/` into the work folder, runs the `kofl
call` sequence the gallery page lists (`CALLS`), then for each shot
(`SHOTS`) writes a `murmel run --script` file (open `volume.vtu`, colour
by `BlockId`, explode, camera, `export.image` 1600 x 1000 at scale 2,
light theme) and puts the shots of an entry side by side into
`docs/figures/gallery/<figure>.jpg` (800 pixels wide each, JPEG quality
88). Each run's tool results are in `<work>/<entry>/results.json`; the
numbers on the gallery page come from them.

[poly_gallery.py](poly_gallery.py) makes the page
[Polyhedral and layered meshes](../../docs/gallery/polyhedral.md) the same
way: `MURMEL=/path/to/murmel python3 tools/gallery/poly_gallery.py target/release/kofl target/gallery [entry ...]`.
Each entry (`SCRIPTS`) is the command block printed on the page, run by
bash in the work folder, which links `test-fixtures/` and `examples/`
and collects the outputs in `out/`; the shots (`SHOTS`) open the meshes'
`foehncraft.vtu`, colour by patch or cell type and cut them with clip
planes. The seven-tube furnace needs `pip install gmsh numpy` (and, on a
bare Linux, `libglu1-mesa`). `RENDER_ONLY=1` renders from the outputs of
an earlier run without meshing again.
