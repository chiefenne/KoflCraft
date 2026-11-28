- **H-grid for rectangular ducts** (`regions.fill`). A section with four
  sharp corners is now meshed as the wall layers around a tensor core
  that is clustered to all four walls, in place of the O-grid whose
  columns fan towards a smaller core. On the elbow the worst
  non-orthogonality drops from 59.1° to 49.5° (at the mitre) and from 45°
  to 26.6° in the straight legs, and the smallest scaled Jacobian rises
  from 0.37 to 0.49, for about twice the cells. Round pipes keep the
  O-grid. The cut caps of a square duct are now Delaunay, so the square T
  joins to the box mesh with two wall layers.
