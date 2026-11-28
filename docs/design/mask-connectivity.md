# Design note: mask connectivity (SQUID landing, flux-line return, ground straps, control-line coupling)

Status: **proposed, not ready to merge.** The cell files, masks, previews and reviews were regenerated with it, but the
adversarial review found three high-severity electrical problems (below, "Open review findings").

**Open review findings (PR #28 adversarial review; these block merging):**
1. *Shared ground bridge (high).* The SQUID ground lead and the flux-line short sit in the same sealed (−x, −y) pocket, which reaches the plane through one modelled airbridge (~100 µm), chained up to 12 / 22 / 34 bridges deep on AF-17 / BF-49 / GF-97. Qubit plasma current and the 50 Ω flux-line return share that bridge's inductance: a galvanic T1 channel estimated at T1 ≈ 0.1–1 µs, which `line_coupling` cannot see (it treats the pocket as ideal ground). Fix: land the SQUID lead and the flux short on continuous plane metal; gate that each junction's ground side is on the plane component itself, not reached through a strap.
2. *Flux line cannot tune the qubit (high).* The short at (−61, −61) µm gives a mutual inductance of ~0.004–0.011 pH to the SQUID loop, i.e. ~90–290 mA for Φ0/2. The flux lines are therefore **not** fixed. Fix: run the line current along the loop and return it on continuous plane; compute M and gate I(Φ0/2) ≲ 1–2 mA.
3. *XY drive is coupler-dominated (high).* The XY bump pad sits ~30 µm from two coupler pads (the sliver fill removes the 2 µm ground between their cut-outs); `line_coupling` draws couplers without pads. A window extraction of the assembled mask gives C(XY, coupler) ≈ 1.4–1.6 fF against C(XY, qubit) ≈ 0.036 fF: neighbour crosstalk 60–160 %. **The 0.03 fF / 23.8× control-line figures do not hold.** Fix: move the bump or keep a ≥ 4 µm ground strip the fill leaves alone; extract from a mask window; add coupler-mediated drive, crosstalk and a coupler T1 budget.
4. *Connectivity as ideal ground (medium).* The connectivity check treats a spanning tree of single bridges as ideal ground (one 36 µm bridge carries the returns of up to 77 circuits). The junction-landing rule accepts any non-plane net at the island end (a mutant with SQUIDs on the XY pad passes `layout_presets`). The control-line model and the XY move have no tests. `line_budget` uses the 5.5 GHz flux coupling at every f01, overstating the flux T1 limit by (f01/5.5)² (BF-49 9.7× → 8.2×).


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

## Changes

- **SQUID ground lead.** `xmon` cuts a 10 µm wide tongue of ground (`kTongueW`) out of the qubit's cut-outs, from the
  (−x, −y) pocket up beside the stub, 10 µm (`kSquidGap`) from the stub's −x face and 10 µm below the −x arm. The SQUID
  region bridges stub → tongue, i.e. island → ground. The tongue is part of the extraction cells (it is drawn by `xmon`),
  so the cells were re-synthesised (below).
- **Flux line shorted to ground.** The flux line is shorted again into the (−x, −y) pocket next to the SQUID; the pocket
  is tied to the plane by ground straps.
- **Ground straps (modelled airbridges).** `detail::place_ground_straps` ties every component that may be grounded
  (floating ground, and components holding only shorted lines with their pocket) to the plane: 10 µm wide bridges with
  10 µm × 10 µm feet on ground metal of two different components, spans up to 150 µm, chosen by Kruskal over candidate
  straps (cost: span + 40 µm per line crossed + 80 µm per coupler body or bus strip crossed). A strap never crosses a
  qubit, and a foot never touches island metal. Straps are written to GDS layers 4 (span) and 5 (feet).
  Crossing a coupler or a bus is unavoidable for the plaquette interiors and the coupler–bus strips. A bridge adds
  about ε0 · 10 µm · w / h of capacitance to what it spans (≈ 0.2 fF over a 6 µm bus at 3 µm height, against
  C_bus ≈ 65 fF); this is **not** in the extraction cells.
- **Small fragments removed.** Ground fragments that no strap can reach and smaller than 2500 µm² (slivers of a pocket
  cut by a line, ≈ 150 µm²) are etched away.
- **Unguarded sliver fill.** The mask clean-up no longer keeps narrow ground necks to protect connectivity
  (`guarded_fill` is still in the library): ground a fill cuts off is strapped back or removed. AF-17's four necks are
  gone.
- **XY line moved out.** The XY line now ends open outside the qubit's cut-out, at (+247, −279) µm from the qubit centre
  (`kXYEnd`, `kXYLevel`), in the free corner between the two coupler pads, with its bump pad beside it.
- **Control-line budget.** `line_coupling` extracts the qubit with its flux and XY lines (star-cell geometry plus both
  lines up to their bump pads). The flux line is shorted at the SQUID: a point at path distance z from the short of a
  line matched at its far end sees Re Z = Z0 sin²(βz), so the qubit's coupling to the flux line is its capacitance to each
  panel of the line weighted by sin(βz) (panel charges with the qubit at 1 V; the line is drawn open without its last
  two squares so that it is a separate conductor). `line_budget` evaluates the T1 limit of each line at every qubit's
  idle frequency and asks for at least 10× the qubit's design T1 (`kLineT1Margin`). The mutual inductance of the flux
  line to the SQUID loop (the flux coupling itself, and its own T1 channel) is not modelled.
- **Junction check.** A junction now counts as landed only if it joins an island to the ground plane (both are grounded
  transmons), with straps taken into account.

## Results

| (this revision, from the reviews) | AF-17 | BF-49 | GF-97 |
|---|---|---|---|
| mask DRC, width + space violations | 0 + 0 | 0 + 0 | 0 + 0 |
| floating ground fragments | 0 | 0 | 0 |
| shorted lines reaching the plane | 58 of 58 | 178 of 178 | 362 of 362 |
| SQUIDs joining an island to the plane | 41 of 41 | 129 of 129 | 265 of 265 |
| ground straps (over lines / over coupler or bus strips) | 122 (47 / 41) | 373 (129 / 146) | 774 (265 / 315) |
| XY / flux coupling, measure qubit (fF) | 0.029 / 0.030 | 0.030 / 0.032 | 0.030 / 0.032 |
| control-line T1 limit / design T1, worst qubit (XY, flux) | 23.8×, 22.4× | 11.4×, **9.7×** | 11.4×, **9.7×** |

Re-synthesised cells (fine mesh, against the design targets): AF-17 E_C within 0.12 %, couplings within 5.0 %, readout
capacitance −0.5 %; Beta-1 (also GF-97) E_C within 0.08 %, couplings within 5.4 %, readout capacitance +0.3 %. Both
re-null the parked exchange.

Gating: the mask DRC and metal connectivity pass on all three chips and are gated. The control-line budget is missed
by the flux line on the Beta-1 cells (9.7× against 10×), so it is reported, not gated, until all three chips meet it.
Closing that gap needs the flux line's approach to the qubit reworked (for example a shorter run beside the T-bar tips).

## Not modelled / open

- The capacitance of the ground straps to what they span (above).
- Inductive coupling of the flux line to the SQUID loop and the resulting flux-noise and T1 channel.
- Any foundry DRC deck (airbridge rules included); the in-house checks are the 4 µm / 4 µm merged-etch width/space
  check, the placement checks and the connectivity check.
