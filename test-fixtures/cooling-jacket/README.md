# Synthetic cooling jacket

A clean-room stand-in for a motor cooling jacket, the fluid volume only (mm): a thin annular shell (radii 60 and 66, height 60) whose inner side carries six ribs 4 deep and 2 wide, each a ring with a 24° gap so the coolant passes between the channels; two blocks on the outer wall; two pipes of radius 4 leaving radially and bending by 45° (torus bends of centreline radius 12, swept as one tube each), each ending in a flat circular cap (96 sides). The shell is closed, genus 1 (the hole through the ring) and nonconvex. 24,732 facets.

`make_jacket.py` builds it with [manifold3d](https://github.com/elalish/manifold) (Apache-2.0): `python3 -m pip install manifold3d numpy && python3 make_jacket.py cooling-jacket.stl`. The ports are not in contiguous facet ranges; `patch.detect_ports` finds them ([workflow](../../docs/nonconvex-duct.md#from-a-bare-stl)).

`cooling-jacket-thin.stl` (`python3 make_jacket.py cooling-jacket-thin.stl thin`) carries the two thin features of the industrial part: rib 3 reaches to 0.5 mm below the outer wall (a thin fluid gap) and rib 5 is 0.5 mm wide (a thin solid wall). It reproduces the boundary recovery blowup that the thickness size source and the remesh rules of [research/thin-gaps.md](../../research/thin-gaps.md) remove.
