- **Parallel tubes and reducers are swept** (`regions.propose`,
  `regions.fill`). A run whose level sets meet the wall in one closed loop
  is a sweep, also between two headers; a second loop (a body inside, as
  in the furnace) makes it an obstacle. Boxes are found by section area
  instead of volume, so the wide pipe of a reducer is no longer a box. The
  fill sweeps tubes cut at both ends and meshes every closed shell behind
  the cuts on its own: the heat exchanger's tube side now joins its seven
  tubes to both bonnets. A sweep's cuts stay
  in one unbroken run of one-loop level sets, so a tube never ends in a
  cut across a box.
