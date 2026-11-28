- **More than two openings behind the cuts** (`regions.fill` with join).
  The polyhedral mesher now gives each connected cap of a port id its own
  layer similarity, so several cuts can share one id, and the join tells
  them apart by their planes. A box with three round pipes joins into one
  closed, valid mesh of 105,638 cells, and FoehnCraft reads it as valid.
  A cut that lies next to a ring of STL nodes no longer leaves a strip of
  needles that the scaffold could not recover. Square sections (the T)
  still fail at the cut.
