# Colour maps

KoflCraft, FoehnCraft and MurmelCraft use one set of colour maps
(maintainer's decision, 2026-11-24), so a mesh, its quality and the flow
on it read the same in every picture of the three projects.

- **Rainbow Desaturated** is the map of every field picture: a quality
  measure drawn as a field, a mesh size, a wall distance, a flow result.
- **Inferno** is the second map, for a second field next to the first.
- A diverging map (cool to warm) stays for signed values centred on zero.

KoflCraft draws no pictures of its own: its meshes are shown by
MurmelCraft (`murmel run --script` with `export.image`, as in the README
pictures), which has these maps built in and Rainbow Desaturated as its
default. Any picture or plot added to this repository takes them too.

## The shared table

Rainbow Desaturated is the preset of that name in ParaView
(`Remoting/Views/ColorMaps.json` in github.com/Kitware/ParaView): eight
colours at the positions below, mixed linearly in RGB. Copyright (c) 2000
Kitware Inc., BSD 3-Clause licence. Every value of the source is a whole
number of 255ths, so the colours are exact.

| Position | Colour |
|---|---|
| 0.0 | `#4747DB` |
| 0.143 | `#00005C` |
| 0.285 | `#00FFFF` |
| 0.429 | `#008000` |
| 0.571 | `#FFFF00` |
| 0.714 | `#FF6100` |
| 0.857 | `#6B0000` |
| 1.0 | `#E04D4D` |

The same table is in MurmelCraft (`assets/colormaps.toml`, test
`rainbow_desaturated_is_paraviews`) and FoehnCraft
(`tools/plot/colormaps.py`, `python/tests/test_colormaps.py`, D261);
`tests/colour_maps.rs` here checks that this page still holds it.

Inferno is matplotlib's map (Nathaniel J. Smith, Stefan van der Walt and
Eric Firing, CC0 1.0); ParaView calls it "Inferno (matplotlib)".
