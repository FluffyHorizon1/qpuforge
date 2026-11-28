# Alpha-1 v1.1 (AF-17) — design specification

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

| | |
|---|---|
| **Status** | Simulation gate **passed**: every criterion holds for every pair, including CZ leakage, and the layout passes the in-house mask checks: placement and spacing, the 4 µm width/space DRC (0 violations), metal connectivity with modelled ground straps (clean) and the control-line T1 budget (§5, §6). **Not yet checked:** a foundry DRC deck (airbridge rules included), the straps' capacitance to the coupler and bus strips they cross, and the flux lines' inductive coupling to the SQUIDs (§7). |
| **Supersedes** | Alpha-1 v1.0 (`designs/alpha1/ALPHA-1.md`, qpuforge v0.2) |
| **Purpose** | First error-correction-ready chip: 17-qubit distance-3 rotated surface code (memory experiment) |
| **Tooling** | qpuforge v0.4: `alpha1_v11_review` (399 s), `qpuforge gds`, `qpuforge cz`, `qpuforge hwgate` |

All performance numbers are **simulation results** under the assumptions in §8.

## 1. What changed from v1.0

| Area | v1.0 | v1.1 |
|---|---|---|
| Tunable coupler | Abstract, g_max = 40 MHz. Idle ZZ: 5 kHz typed (24 Nov v1.0 review); since PR #2 computed by `alpha1()` from a coupler transmon at the exchange null, −4 to −14 kHz | **Physical coupler transmon** per pair: f_max 8.5 GHz (park at the sweet spot), α_c −130 MHz, g_qc 160 MHz |
| Direct coupling | Not modelled in the 24 Nov v1.0 review (4 MHz); since PR #2 `alpha1()` has 1 MHz (0.035 fF), cancelled at the coupler's exchange null | **Per-pair g12 = 9.4–9.8 MHz**, chosen so the parked coupler cancels the qubit–qubit exchange: no avoided crossings when a measure qubit sweeps past its other neighbours |
| CZ gate | 23.3 ns, 2e-3 control error (assumed) | **Pulse-level**: sin² flux ramps (7 ns) + flat top, 3-knob calibration per pair. 35.6–49.2 ns, simulated coherent error 7.9e-5, leakage 4.8e-5 |
| Leakage | Folded into control error | **Simulated**: CZ leakage of the moving qubit, leaked-partner phase scrambling, transport (5 %), seepage at 2/T1 |
| Reset | Active reset 200 ns, 0.3 % | Same, made **leakage-removing** (multi-level reset: \|2⟩ → \|0⟩) |
| Layout | Floor plan (SVG) | **MoM-synthesised cells, assembled chip, DRC, GDSII**, post-layout circuit check |

The design file `alpha1_v11.json` keeps the 4 MHz bare direct coupling (0.14 fF per link) it was generated with:
`alpha1_v11()` is built with `alpha1_with_coupling` at 4 MHz, whereas `alpha1()` (v1.0) has used 1 MHz since PR #2.
The per-pair g12 in the table above is set by the pulse-level calibration itself (`pulse.hpp`).

## 2. Chip summary

| | |
|---|---|
| Qubits | 17 flux-tunable transmons (SQUID), all idling at their sweet spot: 9 data (5.11–5.29 GHz), 8 measure (5.81–5.99 GHz); α ≈ −220 MHz (also post-layout) |
| Couplers | 24 grounded tunable-coupler transmons (one per surface-code CZ pair) |
| Readout | 17 λ/4 resonators at 6.90–7.20 GHz (60 MHz steps), κ = 4 MHz, 2χ = κ, one Purcell filter per feedline, 3 feedlines; 500 ns readout, 0.5 % separation error |
| Reset | Active, 200 ns, 0.3 % residual excitation, leakage-removing |
| Coherence (model) | T1 65–76 µs (dielectric Q = 3e6, Purcell), Tφ 120–170 µs, T_eff 45 mK |
| Chip | 5.8 × 5.8 mm, qubit pitch 900 µm |
| Lines | 17 XY + 17 qubit flux + 24 coupler flux + 3 readout in + 3 readout out = **64**. The 58 control lines and the 17 resonator–feedline couplings end on flip-chip bump pads (75); the feedlines run on the interposer |

## 3. CZ gate (pulse level)

**Protocol.** The measure qubit (the higher frequency) is the *mover*. For each CZ:

1. Its flux pulse brings |11⟩ into resonance with |20⟩ (f_mover + α = f_data). The interaction point is 5.33–5.50 GHz, an excursion of 308–665 MHz.
2. At the same time the coupler swings from 8.5 GHz down to ≈ 6.03–6.14 GHz, which turns on the |11⟩–|20⟩ coupling (gap 36 MHz for the 5.99/5.29 GHz test pair with the coupler at 6.10 GHz, `results/test_output.txt`).
3. One full |11⟩ → |20⟩ → |11⟩ cycle gives the conditional phase π.

Three knobs are calibrated per pair:
- the mover's interaction flux;
- the flat-top length;
- the coupler bias, which also closes the off-resonant |01⟩–|10⟩ swap.

The calibration minimises 1 − F_avg. If a pair then fails the 2e-4 leakage limit, it is calibrated again with leakage
weighted ten times; the repeat is kept only if it meets the limit, its 1 − F_avg is no worse and it stays within 5e-4
(CZ model version 4, `../../docs/design/cz-leakage.md`). One of the 16 frequency combinations needs this: Q8–Q11 goes
from 2.3e-4 to 1.6e-5 leakage (1 − F_avg 4.1e-4 → 3.6e-4, gate 3.5 ns longer, total CZ error +0.5 %).

The simulation covers:
- exact transmon charge-basis spectra for each element, truncated to 4/4/4 levels (qubit/coupler/qubit; qpuforge v0.3 release, was 4/3/4, see CHANGELOG);
- a smooth sign gauge for the charge matrix elements;
- piecewise-exact propagation with dt = 0.1 ns (converged to 6e-7 in F);
- the single-qubit rotating frame, so idle ZZ counts as a physical error.

| Quantity (24 pairs) | Value |
|---|---|
| Duration | 35.6–49.2 ns (mean 43.7 ns) |
| Coherent error (1 − F_avg, closed system) | 7.9e-5 mean, 3.6e-4 max (Q8–Q11) |
| Leakage per CZ (moving qubit) | 4.8e-5 mean, 1.2e-4 max (Q8–Q10) |
| 1/f flux noise (local dephasing at the excursion + quasi-static sensitivity) | 3.2e-5 |
| Residual idle ZZ at park | 3.4–4.5 kHz |
| Spectator crossings (Landau–Zener, nulled exchange) | < 1e-15 |
| Residual control allowance (assumed: distortion, drift, TLS, crosstalk) | 1.0e-3 |
| **Total CZ error incl. T1/Tφ during the gate** | **1.86e-3** |

The per-pair table is `alpha1_v11_review.json → cz_calibration.per_pair`. The calibration cache is `results/cz_cache.json`.

## 4. Error-correction performance (distance-3 memory, 5–25 rounds, Pauli-frame Monte Carlo + exact MWPM)

| | Logical error per round | Round |
|---|---|---|
| v1.0 model | 0.454 % ± 0.006 % | 847 ns |
| **v1.1 model (pulse CZ + leakage + leakage-removing reset), memory-Z** | **0.377 % ± 0.004 %** | 928 ns |
| v1.1 model, memory-X (\|+⟩ logical state, phase flips) | 0.360 % ± 0.005 % | 928 ns |
| Logical lifetime | **246 µs** = 3.3 × best physical T1 (75.6 µs) | |
| Distance-5 projection (same physics, 49 qubits) | 0.123 % ± 0.007 % → **Λ = 2.95 ± 0.18** | |

The v1.0, lifetime and distance-5 rows are memory-Z. The two bases differ by 0.017 %, which is 2.7 standard errors.
The ± values are 1σ Monte Carlo errors of the fits. The distance-5 projection compares distance 3 and distance 5 over
5–15 rounds, a shorter fit than the headline rows. Since qpuforge 0.4 the Monte Carlo is thread-count independent
(`docs/design/deterministic-mc.md`): a rerun gives the same numbers bit for bit at any `OMP_NUM_THREADS`.

The v1.0 row is the abstract v1.0 model run inside this review with the current `alpha1()`. Since PR #2 that is 1 MHz
direct coupling and a coupler transmon parked at the exchange null (idle ZZ −4 to −14 kHz, computed), not the 4 MHz
coupling and typed 5 kHz idle ZZ of the 24 Nov v1.0 review (`designs/alpha1/ALPHA-1.md`).

v1.1 is better than v1.0 despite the slower CZ. The physically modelled control error (≈ 1.06e-3) is about half the
2e-3 that v1.0 assumed.

**Leakage handling** (eps over 5–25 rounds; P_L after 50 rounds):

| | eps / round | P_L (50 rounds) | leaked data / measure |
|---|---|---|---|
| Leakage not modelled | 0.366 % | 0.150 | — |
| No leakage removal | 0.395 % | 0.163 | 1.0e-3 / 1.4e-3 |
| **Leakage-removing reset (chosen)** | **0.366 %** | **0.153** | 4.9e-5 / 9.7e-5 |
| + DQLR each round (30 ns, 5e-4 error) | 0.402 % | 0.168 | 5.5e-6 / 8.8e-5 |
| Stress, 10× CZ leakage: LRR / LRR + DQLR | — | 0.183 / 0.191 | 3.9e-4 / 4.5e-5 |

**Decision:** the leakage-removing reset is enough for Alpha-1. DQLR removes data-qubit leakage but costs more than it
saves at distance 3, even at 10× leakage. It is carried to Beta as an option to re-test at distance 5.

**Error budget** (eps change when one source is removed): T1 −80 %, CZ/1q control −40 %, Tφ −32 %, readout −18 %,
reset −13 %, leakage −5 %, thermal −3 %, ZZ +1 %. Changes of a few per cent are within the Monte Carlo noise.

**Sensitivity** (single change):

| Change | eps change |
|---|---|
| Dielectric Q ×2 | −36 % |
| T_eff 35 mK | −27 % |
| Readout 300 ns | −22 % |
| CZ residual 5e-4 | −19 % |
| Readout separation 0.25 % | −10 % |
| Flux noise ×3 | +9 % |
| CZ residual 2e-3 | +37 % |

These ranks set the Beta-1 upgrade list.

## 5. Physical layout

**Cells.** An Xmon-style cross with a T-bar at each arm end faces each coupler's T-bar across a 5 µm gap. The readout
claw wraps a pad in the cross's inner corner. The coupler is a grounded transmon: two T-bars, a 10 µm body and a square
mid pad; its SQUID sits on the pad edge with an inductive flux line. A thin floating strip (the bus) runs parallel to each coupler and sets the direct coupling.
The target direct coupling is weak, so the bus has no T-bar: it is a 6 µm strip whose ends stop 17.5 µm short of the
qubit T-bars. (The 0.4 synthesis shrank a bus T-bar to 75 nm instead; the cell generator no longer draws metal narrower
than the strip.) It adds 1.2 fF of loading.

The dimensions were synthesised by Gauss-Newton against the circuit targets, using method-of-moments extraction
(quasi-static, ε_eff = (11.45 + 1)/2 on silicon) with a coarse-to-fine mesh correction:

| Cell | Dimension | Value |
|---|---|---|
| Measure / data qubit | arm length (centre to T-bar edge) | 123.4 / 125.7 µm |
| | arm width / gap / T-bar | 16 / 30 / 120 × 12 µm |
| | readout pad / claw gap | 38.2 / 5 µm |
| | SQUID: stub / span to the ground lead / lead width | 16 × 16 / 10 / 10 µm |
| Coupler | T-bar (facing the qubit) | 42.7 µm long, 5 µm gap |
| | body / trench / mid pad | 10 / 20 / 246.1 × 246.1 µm |
| Bus | strip / trench / end gap to the qubit T-bar | 6 / 15 / 14.1 µm |

| Extracted (fine mesh) | Value | Target | Error |
|---|---|---|---|
| E_C measure / data | 202.4 / 200.0 MHz | 202.2 / 199.7 MHz | +0.09 % / +0.12 % |
| E_C coupler | 125.6 MHz | 125.6 MHz | −0.03 % |
| G measure–coupler / coupler–data (GHz per n n) | 0.0605 / 0.0585 | 0.0576 / 0.0612 | +5.0 % / −4.3 % |
| G measure–data (direct) | 0.00566 | 0.00556 | +2.0 % (aimed 1.5 % high, see below) |
| C qubit–readout | 9.10 fF | 9.15 fF | −0.5 % |
| C_Σ qubit / C coupler / C_qc / C_qq | 96.0–97.1 / 155 / 5.61 / 0.099 fF | | |

**Post-layout circuit check.** The junctions are retargeted to the design frequencies and the representative pair is
rebuilt from the extracted circuit:
- α = −220 MHz; coupler α −130 MHz; g_qc = 168 MHz; readout g = 155 MHz (design 156).
- g12 = 9.72 MHz against the 9.54 MHz needed for an exchange-free park at the sweet spot. This would leave a 0.37 MHz parked gap (3.3e-5 transfer per spectator crossing).
- Because g12 is *above* the null value, parking each coupler 102 MHz below its sweet spot (8.398 GHz, df/dΦ = 3.0 GHz/Φ0) nulls it exactly. ZZ is then −4.4 kHz. The synthesis aims the direct coupling 1.5 % above target for this reason: a g12 just below the null value could not be re-nulled.
- Calibration plan: per-pair park flux = exchange null.

**GDSII** (`AF17_v11.gds`, cell `AF17_V11`, 1 nm database unit):

| Layer | Content |
|---|---|
| 1 | Etch: ground-plane cut-outs minus metal islands (11,662 polygons) |
| 2 | Josephson junction / SQUID regions (41) |
| 3 | Flip-chip bump pads (75): XY, qubit flux, coupler flux, resonator–feedline coupling |
| 4 / 5 | Ground straps (modelled airbridges, 122): span / feet |
| 10 | Chip outline |
| 20 | Labels |

- Placement checks (10 µm metal spacing, 4 µm at designed coupling gaps, metal inside the outline): **clean**.
- Mask DRC: minimum width and space of 4 µm on the merged etch layer, projection metric: **0 width and 0 space violations** (34 + 202 in 0.4, 0 + 4 in the previous revision). A clean-up pass merges cut-outs closer than 4 µm (ground slivers where the diagonal lines pass the corners of the qubit cut-out); it removes ground only, never island metal. Ground it cuts off the plane is strapped back (below). Corner-to-corner (Euclidean) spacing is not checked; the foundry's DRC deck must still be run.
- Metal connectivity (LVS-lite, gated): **clean** (`alpha1_v11_review.json → layout.connectivity`; design note `docs/design/mask-connectivity.md`).
  - Each qubit SQUID bridges the stub to a 10 µm ground lead cut out of the (−x, −y) inner corner: island to plane (41 of 41 junctions, couplers included).
  - Each qubit flux line is shorted into that corner next to the SQUID (17 of 17 reach the plane; 58 of 58 shorted lines in all).
  - The inner corners of each qubit, the strip between each coupler and its bus, and each plaquette are enclosed by cut-outs. 122 ground straps (modelled airbridges, 10 µm wide, ≤ 150 µm span, GDS layers 4/5) tie them to the plane: 47 cross a line, 41 a coupler body or bus strip, none a qubit. No floating ground remains.
  - Not modelled: a strap's capacitance to the coupler or bus it crosses (about 0.2 fF over a 6 µm bus at 3 µm bridge height, against C_bus ≈ 67 fF).
- Control lines (gated): the XY line ends open outside the qubit cut-out, 247 µm right of and 279 µm below the qubit centre. Extracted with the qubit (`line_coupling`): XY 0.029 / 0.030 fF (measure / data qubit); flux line 0.030 fF effective, i.e. the capacitance to each piece of the shorted line weighted by sin(βz) from the short (2.0 fF if it were open). T1 limits at the idle frequencies: at least 1.56 ms (XY) and 1.47 ms (flux), 23.8× and 22.4× the design T1 of the worst qubit (Q16); budget 10×. In the previous revision the XY line ran into the qubit's corner and coupled 3.0 fF (T1 limit 0.2 µs). The flux line's mutual inductance to the SQUID loop is not modelled.
- All 17 resonators are drawn at their λ/4 length (4.17–4.35 mm).
- Previews: `AF17_v11_layout.png`/`.svg`, `AF17_v11_detail.png`/`.svg`.
- Independent reader: `tools/gds_preview.py`.

## 6. Simulation gate (ROADMAP.md)

| Criterion | Threshold | Result |
|---|---|---|
| d=3 logical error per round | ≤ 0.6 % | 0.377 % (memory-Z), 0.360 % (memory-X) ✔ |
| Logical lifetime vs best T1 | > 2× | 3.3× ✔ |
| CZ coherent error / leakage, every pair | ≤ 5e-4 / ≤ 2e-4 | largest pair 3.6e-4 / 1.2e-4 (means 7.9e-5 / 4.8e-5) ✔ |
| Idle ZZ at park | ≤ 10 kHz | 3.4–4.5 kHz, mean 3.9 kHz (4.3 kHz re-nulled post-layout) ✔ (computed: `alpha1_v11_review.json → exit_criteria`) |
| Mask | placement/spacing checks and 4 µm width/space DRC clean | clean; 0 + 0 ✔ |
| Metal connectivity | conductors separate; shorted lines and all ground on the plane; SQUIDs island-to-plane | clean (122 ground straps) ✔ |
| Control lines | each qubit's XY and flux line T1 limits ≥ 10× its design T1 | 23.8× (XY), 22.4× (flux) ✔ |
| Layout | E_C within 3 %, couplings within 6 %, exchange re-nullable | 0.12 %; 5.0 %; yes ✔ |
| Fabrication Monte Carlo, 2 % junction spread | ≥ 90 % collision-free | 98 % (v1.0 review; same frequency plan) ✔ |

## 7. Before tape-out (Alpha fabrication checklist)

1. Full-wave/FEM cross-check of the qubit, coupler and readout cells (Palace / Sonnet / Q3D). Then a second synthesis iteration on the residual errors (up to 5.0 % on the qubit–coupler couplings).
2. Interposer design: feedlines, Purcell filters, edge launchers, wiring of the 75 bump pads.
3. Airbridge process for the 122 modelled ground straps (and their capacitance to the coupler and bus strips they cross, in the extraction); the foundry DRC deck; the flux lines' mutual inductance to the SQUIDs; dicing lanes, alignment marks, test structures (junction test arrays, resonator Q).
4. Junction process: target Rn per qubit from the design JSON, with ≤ 2 % spread (or trimming).
5. Calibration plan: idle frequencies; coupler exchange-null park per pair; the 3-knob CZ; leakage-removing reset tune-up.
6. Characterisation: follow `CHARACTERISATION-PLAN.md`. It sets the order of measurements and the data file, and `qpuforge hwgate` turns that file into the hardware-gate verdict of ROADMAP.md, including model versus measurement. The most important single number is the CZ residual control error.

## 8. Assumptions and limitations

- The CZ residual control allowance (1e-3) is assumed, not simulated, and it dominates the CZ error.
- Pulse simulation is closed-system. Decoherence during the gate is applied by the emulators; 1/f flux noise is applied analytically.
- Leakage is a frame-level model. The decoder is leakage-unaware.
- Capacitances come from a quasi-static extraction on a thick substrate. Not included: box modes, kinetic inductance, crosstalk, participation-ratio TLS loss.
- The mask passes the in-house checks only (§5): no foundry DRC deck has been run; the ground straps are modelled airbridges whose capacitance is not in the extraction; the flux lines' inductive coupling to the SQUIDs is not modelled. Control wiring and feedlines on the interposer are not drawn.
- Coherence comes from the device model (Q = 3e6 dielectric loss, Purcell, 2 µΦ0 flux noise, photon shot noise), not from measurement.

## 9. Files

| File | Content |
|---|---|
| `alpha1_v11.json` | Chip design (qubits, couplings, coupler element, leakage spec) |
| `alpha1_v11_review.json/.txt` | Full review output |
| `layout_cells.json` | Synthesised cell geometry, achieved vs target circuit, synthesis log |
| `AF17_v11.gds` | GDSII layout |
| `AF17_v11_layout.png`/`.svg`, `AF17_v11_detail.png`/`.svg` | Previews rendered from the GDSII by an independent reader |
| `CHARACTERISATION-PLAN.md` | What to measure on the fabricated chip, the data-file format, and how the hardware gate is evaluated |
| `example_characterisation_synthetic.json`, `example_gate_report_synthetic.json` | A virtual chip's data file and its gate report (synthetic; for rehearsing the tool chain) |
| `../../results/cz_cache.json` | Pulse-level CZ calibrations (coherent part) |
