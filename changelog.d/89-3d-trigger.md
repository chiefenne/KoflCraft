- **3D trigger of the prism split.** `block.optimize_layout` on a topology
  that is not extruded along z (or has T-junctions) now finds prism-shaped
  blocks by geometry, tries the split with the T-junction stop and with the
  sheet, measures the 3D grid (non-orthogonality, growth, scaled Jacobian)
  and keeps a split only when its blocks get at least 2° better with a worst
  scaled Jacobian of 0.4 or more, while the whole grid gets no worse and
  grows by at most 10 %. The report (`mode: 3d`) names every split kept and
  every candidate left, with the reason.
