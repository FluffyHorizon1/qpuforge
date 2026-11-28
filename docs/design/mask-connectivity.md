# Design note: mask connectivity (SQUID landing, flux launcher, ground bumps and straps, control-line coupling)

Status: **redesign in progress; not ready to merge.** The geometry and the checks below are implemented and pass on
AF-17 with re-synthesised cells in a development run. The committed cell files, masks, previews, reviews and specs are
still those of the previous geometry, and BF-49 and GF-97 have not been taken through the gates. The numbers in
"Results" are from that development run, not from committed review outputs; they will be replaced by the regenerated
reviews.

## Problem

Release 0.4 added a metal-connectivity check (`ConnectivityReport`, LVS-lite on the merged etch layer) and reported it
without gating. On the three committed masks it found:

| (0.4 masks, from the committed reviews) | AF-17 | BF-49 | GF-97 |
|---|---|---|---|
| floating ground fragments | 143 | 405 | 821 |
| shorted lines reaching the ground plane | 16 of 58 | 33 of 178 | 49 of 362 |
| qubit flux lines reaching the ground plane | 0 of 17 | 0 of 49 | 0 of 97 |
| junctions landing on metal at both ends | 24 of 41 | 80 of 129 | 168 of 265 |
| DRC space violations (ground necks kept by the connectivity guard) | 4 | 0 | 0 |

Three causes:

1. **Qubit SQUIDs landed in etch.** The SQUID stub sits in the (−x, −y) inner corner of the cross, inside the arm
   cut-outs. The junction region ran from the stub across 30 µm of etch and ended in the −x arm's cut-out: no ground
   within reach.
2. **Enclosed ground.** With T = 120 µm T-bars and 30 µm gaps, the T-bar cut-outs of neighbouring arms overlap at their
   tips (T/2 + gap = 90 µm > L − t_w − gap ≈ 85 µm), so the three free inner corners of every qubit are sealed pockets.
   The coupler and bus run in parallel bands from qubit to qubit; the strip of ground between a coupler body and its
   bus, and the interior of every plaquette, are sealed by those bands. Resonators shorted into a sealed plaquette, and
   flux lines shorted into a sealed pocket, had no DC return, which is why the flux lines had been drawn open.
3. **Control lines too close.** Extracting the XY and flux lines with the qubit (new: `line_coupling`) showed that the
   XY line, which ran diagonally into the (+x, −y) pocket, coupled 3.0 fF to the qubit, and the flux line 2.1 fF if
   open. For a line matched at its far end, 1/T1 = ω² C_c² Z0 / C_Σ: at 5.5 GHz that is a T1 limit of 0.2 µs (XY) and
   0.4 µs (flux), far below the design T1 (75–157 µs).


### What the first fix got wrong (review of the previous revision)

The previous revision landed the SQUID on a tongue of ground inside the (−x, −y) pocket, shorted the flux line into the
same pocket and tied the pockets to the plane with a spanning tree of single ground straps. Checked with the gates of
this note, on the masks of that revision:

| (previous revision) | AF-17 | BF-49 | GF-97 |
|---|---|---|---|
| qubit SQUID grounds on plane metal without a strap | 0 of 17 | 0 of 49 | 0 of 97 |
| deepest chain of straps to plane metal | 12 | 22 | 34 |
| flux-line mutual to the SQUID loop, thin-film model | ≈ 0 | ≈ 0 | ≈ 0 |
| C(XY, coupler) against C(XY, qubit) ≈ 0.03 fF | 1.61 fF | 1.42 fF | 1.42 fF |
| XY drive on a neighbouring qubit, of its own | 109 % | 123 % | 123 % |

So the qubit's mode current and the flux line's return shared one airbridge, the flux line could not tune the qubit,
and the XY line drove the couplers. The connectivity check called this clean because it treated any strapped net as
ground and any non-plane net as an island.

## Design

All dimensions are in the cell files (`layout_cells.json`) or are named constants in `layout.hpp`.

- **SQUID on the end of a T-bar.** The qubit's SQUID sits on a finger at the −x end of the −y arm's T-bar and lands on a
  4 µm collar of ground (`kCollarW`) on the edge of the (−x, −y) plaquette. The collar is plaquette metal: no strap is in
  the SQUID's ground path. The other inner corners of the cross keep a relief (`corner_relief_gap_um`) so that their
  pockets join the ground outside through necks instead of being sealed.
- **Flux launcher.** The flux line comes from a bump pad about 200 µm away in the same plaquette as a staircase CPW. Ten
  squares before the SQUID the gap on the qubit's side stops and the centre conductor joins the ground there; only one
  slot, on the side away from the qubit, runs on to a 5 µm collar slot behind the collar (`flux_slot_up_um`,
  `flux_slot_down_um`, `flux_slot_feed_um`). The line current returns along the collar past the loop: this is the
  tuning mutual M. The unequal slot arms set how much of the line's current also flows in the qubit's own return path
  (the common-mode mutual M_cm, below).
- **Coupler SQUID and launcher.** Each coupler SQUID sits on the edge of its pad nearest its flux pad, with the same
  kind of launcher (`flux_slot_short_side_um`, `flux_slot_other_side_um`).
- **Ground bumps (GDS layer 6).** The couplers and buses cut the ground into plaquettes that are sealed in 2D. A sealed
  plaquette takes its reference from flip-chip ground bumps to the top chip: two beside every signal bump (flux,
  readout), and a lattice of ground bumps at 50 µm pitch on free ground metal, kept 40 µm from qubits, couplers and
  buses and 20 µm from lines. On metal joined to the chip border every third lattice point is used.
- **Ground straps (GDS layers 4, 5).** Modelled airbridges remain only where 2D metal cannot do the job: stitching
  straps in pairs across each coupler body and bus (nine per coupling), and whatever else is needed for small ground
  fragments. A strap's feet have 2 µm of ground metal around them, a strap keeps 10 µm from qubit metal, and every pair
  of fragments is tied by two straps at least 25 µm apart. No SQUID ground and no flux-line short has a strap in its
  path.
- **Control lines: two modes** (`control_lines` in the cell file). `combined` (default): one line per qubit carries
  the flux bias and the microwave drive; the drive couples through M_cm. `separate`: an XY line per qubit as well,
  ending open in front of a corner pocket. Both are drawn and checked by the same code. Which one the chips use is an
  open decision; the wiring sections of the specs still describe separate XY lines.

## Checks and gates

Connectivity (`analyze_connectivity`, gated by the `gds` command, the tests and the reviews):

- nets: no two islands share metal, every shorted line reaches the ground net, no floating ground, every junction region
  has an owner island and lands owner-to-ground with nothing else on its boundary, exactly one junction per qubit and
  coupler;
- ground quality: every SQUID ground and every qubit flux-line short is on plane metal (the border component or a
  component under a ground bump) with no strap in the path, and the flux line's return is continuous to its bump pad;
  no fragment hangs on a single strap; no chain of straps deeper than 2;
- ground bumps: the report separates SQUID grounds on metal joined to the chip border in 2D from those on plaquette
  metal grounded by bumps alone; a bump-grounded plaquette has at least 4 ground bumps; every SQUID ground on one has
  a bump within 400 µm; every signal bump has a ground bump within 75 µm; bumps are 30 µm squares at least 50 µm apart,
  a ground bump has 3 µm of ground metal around it and keeps 10 µm from every island and strap.

Flux coupling (`flux_coupling`). A thin superconducting film carrying sheet currents is, by Babinet duality, the
electrostatic problem of conductors shaped like the etched apertures. The existing method-of-moments solver is run on
the apertures of a window of the mask; the flux through the SQUID loop per ampere in the line is M, and
I(Φ0/2) = Φ0 / 2M. The same solution gives M_cm: the line's flux linked by the qubit's own island-to-ground current
path. The method reproduces the inductance of a CPW to −1.8 % against conformal mapping.

Line budget (`line_budget`, `line_budget_crosstalk`), evaluated for a measure qubit, a data qubit and a coupler of each
chip, at every qubit's idle frequency:

| gate | limit |
|---|---|
| I(Φ0/2), qubits and couplers | ≤ 2 mA |
| T1 limit of each channel (flux line capacitive, through the loop at asymmetry ≥ 0.1, common mode; XY line; coupler) | ≥ 10 × design T1 |
| drive power for a 25 MHz Rabi rate, at the chip | ≤ −55 dBm |
| drive on a neighbouring qubit, as a fraction of its own drive | ≤ 10 % |

A galvanic path shared by the qubit and a line through a strap is excluded by the connectivity gate, not budgeted.

## Assumptions

These are assumptions of the model, not things it verifies.

- **Top chip.** The top chip carries a continuous ground plane on which all ground bumps land. Above every qubit,
  coupler, bus and resonator, and for 40 µm around them (the lattice keep-out), the top chip is assumed free of metal.
  The capacitance extraction has no top-chip metal; with metal 10 µm above a qubit island its E_C would shift by
  several per cent.
- **Bump inductance.** 2 pH per bump (`kBumpInductance_pH`): the partial self-inductance of a post 30 µm across and
  10 µm tall is about 1.1 pH; the rest allows for current spreading into the bump. The plaquette metal is taken as one
  node on its N bumps in parallel. A line's return current then moves that node by jω(L_b/N)·I against the top ground;
  a qubit grounded on it has three of its four quadrants on metal that does not move, so 0.75·L_b/N is added to its
  M_cm (worst sign), and a neighbour with one quadrant on it sees 0.25·L_b/N. The quadrant fractions are taken as
  equal, the straps between plaquettes are ignored, and the sheet inductance of the plaquette (which makes distant
  bumps less effective than the nearest ones) is ignored.
- **Thin film.** The duality model has no film thickness and no kinetic inductance. For aluminium or niobium films of
  100–200 nm the Pearl length 2λ²/t is of the order of 0.1 µm against features of 4–10 µm, so the current
  distribution changes at the level of 1–3 %: of that order in M, and 0.02–0.05 pH in M_cm, which is a difference of
  terms of about 1.5 pH. For high-kinetic-inductance films (NbTiN, TiN, granular aluminium) the model does not apply.
- **Combined line.** Driving through the flux line needs the bias and the 4–6 GHz drive combined off the chip
  (a diplexer or bias tee at the mixing chamber). The flux line then cannot carry the usual low-pass filter below the
  qubit frequency: the drive arm must be attenuated like an XY line, and the T1 limits here assume a matched 50 Ω line
  with no excess noise at the qubit frequency. Coupler flux lines keep a pure flux role and can be filtered as usual.
- **Separate XY line.** The model is a single metal plane. In a flip-chip assembly the top chip's ground would cut the
  long-range coplanar coupling between an XY line and distant couplers; the crosstalk and coupler-T1 numbers of the
  separate mode are what this model can verify, not a statement about the assembled device.
- **Junctions.** Josephson asymmetry is taken as at least 0.1 where the design gives less.

## Results

AF-17, development run with re-synthesised cells (fine mesh; not yet the committed review):

RESULTS_TABLE

## Not modelled / open

- Full-wave electromagnetics of any part of the chip; the top chip and the bumps beyond the lumped estimate above.
- The capacitance of the ground straps to the couplers and buses they span (about 0.2 fF per strap over a 6 µm strip
  at 3 µm height; not in the extraction cells).
- Any foundry DRC deck (airbridge and bump rules included); the in-house checks are the 4 µm / 4 µm merged-etch
  width/space check, the placement checks, and the connectivity and bump rules above.
- The choice between combined and separate control lines, and with it the line counts in the specs.
- The masks are design-study outputs. They are not fabrication-ready.
