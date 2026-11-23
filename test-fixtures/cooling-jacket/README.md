# Synthetic cooling jacket

A clean-room stand-in for a motor cooling jacket, the fluid volume only (mm): a thin annular shell (radii 60 and 66, height 60) whose inner side carries six ribs 4 deep and 2 wide, each a ring with a 24° gap so the coolant passes between the channels; two blocks on the outer wall; two pipes of radius 4 leaving radially and bending by 45°, each ending in a flat circular cap (96 sides). The shell is closed, genus 1 (the hole through the ring) and nonconvex. 8,120 facets.

`make_jacket.py` builds it with [manifold3d](https://github.com/elalish/manifold) (Apache-2.0): `python3 -m pip install manifold3d numpy && python3 make_jacket.py cooling-jacket.stl`. The ports are not in contiguous facet ranges; `patch.detect_ports` finds them ([workflow](../../docs/nonconvex-duct.md#from-a-bare-stl)).
