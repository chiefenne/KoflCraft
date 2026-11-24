# The README slideshow

`assets/readme/tour.webp` (the slideshow near the top of the README) and
the stills of [docs/gallery.md](../../docs/gallery.md) are pictures of
KoflCraft meshes written by MurmelCraft's headless `export.image`, on
Linux without a screen. The recipe matches MurmelCraft's own README
slideshow, so the two read the same.

To make them again:

1. `cargo build --locked --release -p mesh-cli` here, and in a
   [MurmelCraft](https://github.com/chiefenne/MurmelCraft) checkout
   `cargo build --release -p murmel-cli`.
2. A Vulkan driver for the offscreen render (`apt install
   mesa-vulkan-drivers` gives lavapipe) and Python 3 with Pillow (WebP
   support).
3. `MURMEL=/path/to/murmel tools/readme-tour/make-tour.sh` makes the five
   meshes in `target/readme-tour/meshes` (the Manifold takes about four
   minutes on four cores, the rest seconds), renders the slides and writes
   `assets/readme/tour.webp` and `assets/readme/tour/*.jpg`.

How it works:

- [make-tour.sh](make-tour.sh) meshes the Manifold and elbow
  (`generate-duct … polyhedral`), the block and pipe (`generate-sweep`)
  and the straight channel (`generate-channel`), and skips a mesh that is
  already there.
- [tour.py](tour.py) holds the slides: for each one a `murmel run
  --script` file opens the mesh, sets the colouring, the clip plane and the
  camera with the same tools as the window's controls, and ends with
  `export.image` at 2560 x 1600, scale 2 (lines, legend and text grown to
  match). Change a slide there and run the script again.
- [make_webp.py](make_webp.py) scales the pictures to 1280 x 800 and
  writes an animated WebP with Pillow: each slide holds 4 s, then five
  blended frames of 60 ms crossfade to the next (smoothstep), quality 70,
  method 6, looping; about 1.7 MB for eight slides. The stills are
  1440 x 900 JPEG, quality 88.

What to know:

- A clip plane's "half nearer the eye" is judged from the camera when the
  plane is set: set the camera first, or give the plane an origin and a
  normal.
- GIF would cut the colour maps to 256 colours; GitHub plays animated WebP
  in an `<img>`, and a viewer that cannot animate it shows the first
  slide.
