# Alpha-1 v1.1 (AF-17) — design specification

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

| | |
|---|---|
| **Status** | Simulation gate **passed**: every criterion holds for every pair, including CZ leakage (§6). Layout (GDSII): placement and spacing checks clean; **the mask is not fabrication-ready**: four DRC space violations remain and the metal connectivity check is not clean (floating ground, flux-line return, junction landing; §5). |
| **Supersedes** | Alpha-1 v1.0 (`designs/alpha1/ALPHA-1.md`, qpuforge v0.2) |
| **Purpose** | First error-correction-ready chip: 17-qubit distance-3 rotated surface code (memory experiment) |
| **Tooling** | qpuforge v0.4: `alpha1_v11_review` (353 s), `qpuforge gds`, `qpuforge cz`, `qpuforge hwgate` |

All performance numbers are **simulation results** under the assumptions in §8.

## 1. What changed from v1.0

| Area | v1.0 | v1.1 |
|---|---|---|
| Tunable coupler | Abstract (g_max = 40 MHz, 5 kHz residual ZZ) | **Physical coupler transmon** per pair: f_max 8.5 GHz (park at the sweet spot), α_c −130 MHz, g_qc 160 MHz |
| Direct coupling | Not modelled | **Per-pair g12 = 9.4–9.8 MHz**, chosen so the parked coupler cancels the qubit–qubit exchange: no avoided crossings when a measure qubit sweeps past its other neighbours |
| CZ gate | 23.3 ns, 2e-3 control error (assumed) | **Pulse-level**: sin² flux ramps (7 ns) + flat top, 3-knob calibration per pair. 39.7–49.2 ns, simulated coherent error 8.1e-5, leakage 4.4e-5 |
| Leakage | Folded into control error | **Simulated**: CZ leakage of the moving qubit, leaked-partner phase scrambling, transport (5 %), seepage at 2/T1 |
| Reset | Active reset 200 ns, 0.3 % | Same, made **leakage-removing** (multi-level reset: \|2⟩ → \|0⟩) |
| Layout | Floor plan (SVG) | **MoM-synthesised cells, assembled chip, DRC, GDSII**, post-layout circuit check |

## 2. Chip summary

| | |
|---|---|
| Qubits | 17 flux-tunable transmons (SQUID), all idling at their sweet spot: 9 data (5.11–5.29 GHz), 8 measure (5.81–5.99 GHz); α ≈ −220 MHz (also post-layout) |
| Couplers | 24 grounded tunable-coupler transmons (one per surface-code CZ pair) |
| Readout | 17 λ/4 resonators at 6.90–7.20 GHz (60 MHz steps), κ = 4 MHz, 2χ = κ, one Purcell filter per feedline, 3 feedlines; 500 ns readout, 0.5 % separation error |
| Reset | Active, 200 ns, 0.3 % residual excitation, leakage-removing |
| Coherence (model) | T1 65–76 µs (dielectric Q = 3e6, Purcell), Tφ 120–170 µs, T_eff 45 mK |
| Chip | 5.8 × 5.8 mm, qubit pitch 900 µm |
| Lines | 17 XY + 17 qubit flux + 24 coupler flux + 3 readout in + 3 readout out = **64**, all ending on flip-chip bump pads |

## 3. CZ gate (pulse level)

**Protocol.** The measure qubit (the higher frequency) is the *mover*. For each CZ:

1. Its flux pulse brings |11⟩ into resonance with |20⟩ (f_mover + α = f_data). The interaction point is 5.33–5.50 GHz, an excursion of 308–665 MHz.
2. At the same time the coupler swings from 8.5 GHz down to ≈ 6.04–6.14 GHz, which turns on an effective coupling with a 36 MHz |11⟩–|20⟩ gap.
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
| Duration | 39.7–49.2 ns (mean 43.9 ns) |
| Coherent error (1 − F_avg, closed system) | 8.1e-5 mean |
| Leakage per CZ (moving qubit) | 4.4e-5 mean, 8.4e-5 max |
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
| **v1.1 model (pulse CZ + leakage + leakage-removing reset), memory-Z** | **0.374 % ± 0.004 %** | 928 ns |
| v1.1 model, memory-X (\|+⟩ logical state, phase flips) | 0.368 % ± 0.005 % | 928 ns |
| Logical lifetime | **248 µs** = 3.3 × best physical T1 (75.6 µs) | |
| Distance-5 projection (same physics, 49 qubits) | 0.117 % ± 0.007 % → **Λ = 3.07 ± 0.19** | |

The v1.0, lifetime and distance-5 rows are memory-Z. The two bases differ by 0.006 %, which is 1.0 standard error.
The ± values are 1σ Monte Carlo errors of the fits. The distance-5 projection compares distance 3 and distance 5 over
5–15 rounds, a shorter fit than the headline rows.

v1.1 is better than v1.0 despite the slower CZ. The physically modelled control error (≈ 1.07e-3) is about half the
2e-3 that v1.0 assumed.

**Leakage handling** (eps over 5–25 rounds; P_L after 50 rounds):

| | eps / round | P_L (50 rounds) | leaked data / measure |
|---|---|---|---|
| Leakage not modelled | 0.363 % | 0.151 | — |
| No leakage removal | 0.401 % | 0.163 | 8.7e-4 / 1.3e-3 |
| **Leakage-removing reset (chosen)** | **0.374 %** | **0.154** | 4.6e-5 / 8.9e-5 |
| + DQLR each round (30 ns, 5e-4 error) | 0.406 % | 0.166 | 5.3e-6 / 7.9e-5 |
| Stress, 10× CZ leakage: LRR / LRR + DQLR | — | 0.180 / 0.189 | 3.8e-4 / 4.4e-5 |

**Decision:** the leakage-removing reset is enough for Alpha-1. DQLR removes data-qubit leakage but costs more than it
saves at distance 3, even at 10× leakage. It is carried to Beta as an option to re-test at distance 5.

**Error budget** (eps change when one source is removed): T1 −79 %, CZ/1q control −43 %, Tφ −30 %, readout −17 %,
reset −11 %, leakage −3 %, thermal −2 %, ZZ +1 %. Changes of a few per cent are within the Monte Carlo noise.

**Sensitivity** (single change):

| Change | eps change |
|---|---|
| Dielectric Q ×2 | −37 % |
| T_eff 35 mK | −27 % |
| Readout 300 ns | −25 % |
| CZ residual 5e-4 | −20 % |
| Readout separation 0.25 % | −11 % |
| Flux noise ×3 | +8 % |
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
| Measure / data qubit | arm length (centre to T-bar edge) | 127.4 / 129.8 µm |
| | arm width / gap / T-bar | 16 / 30 / 120 × 12 µm |
| | readout pad / claw gap | 38.4 / 5 µm |
| Coupler | T-bar (facing the qubit) | 42.0 µm long, 5 µm gap |
| | body / trench / mid pad | 10 / 20 / 248.4 × 248.4 µm |
| Bus | strip / trench / end gap to the qubit T-bar | 6 / 15 / 17.5 µm |

| Extracted (fine mesh) | Value | Target | Error |
|---|---|---|---|
| E_C measure / data | 202.3 / 199.8 MHz | 202.2 / 199.7 MHz | +0.05 % |
| E_C coupler | 125.6 MHz | 125.6 MHz | −0.01 % |
| G measure–coupler / coupler–data (GHz per n n) | 0.0603 / 0.0585 | 0.0576 / 0.0612 | +4.6 % / −4.4 % |
| G measure–data (direct) | 0.00563 | 0.00556 | +1.4 % (aimed 1.5 % high, see below) |
| C qubit–readout | 9.10 fF | 9.15 fF | −0.5 % |
| C_Σ qubit / C coupler / C_qc / C_qq | 96.0–97.2 / 155 / 5.61 / 0.104 fF | | |

**Post-layout circuit check.** The junctions are retargeted to the design frequencies and the representative pair is
rebuilt from the extracted circuit:
- α = −220 MHz; coupler α −130 MHz; g_qc = 167 MHz; readout g = 155 MHz (design 156).
- g12 = 9.66 MHz against the 9.54 MHz needed for an exchange-free park at the sweet spot. This would leave a 0.34 MHz parked gap (2.9e-5 transfer per spectator crossing).
- Because g12 is *above* the null value, parking each coupler 95 MHz below its sweet spot (8.405 GHz, df/dΦ = 2.9 GHz/Φ0) nulls it exactly. ZZ is then −4.3 kHz. The synthesis aims the direct coupling 1.5 % above target for this reason: a g12 just below the null value could not be re-nulled.
- Calibration plan: per-pair park flux = exchange null.

**GDSII** (`AF17_v11.gds`, cell `AF17_V11`, 1 nm database unit):

| Layer | Content |
|---|---|
| 1 | Etch: ground-plane cut-outs minus metal islands (14,689 polygons) |
| 2 | Josephson junction / SQUID regions (41) |
| 3 | Flip-chip bump pads (75): XY, qubit flux, coupler flux, resonator–feedline coupling |
| 10 | Chip outline |
| 20 | Labels |

- Placement checks (10 µm metal spacing, 4 µm at designed coupling gaps, metal inside the outline): **clean**.
- Mask DRC: minimum width and space of 4 µm on the merged etch layer, projection metric: **0 width and 4 space violations** (34 + 202 in 0.4). The 0.4 violations came from three places:
  - the 75 nm bus T-bar (now never narrower than the 6 µm strip);
  - a 1 µm etch slit beside the "shorted" end of each qubit flux line (the line now ends open with a full trench; see below);
  - 0.3–3.3 µm ground slivers where the diagonal readout, XY and flux lines pass the corners of the qubit cut-out. A clean-up pass merges cut-outs closer than 4 µm. It removes ground only, never island metal, and a connectivity guard undoes any merge that would cut ground off the plane. The four remaining violations are such necks (0.42 µm at Q7 and Q13, 1.59 µm at Q10 and Q16): each is the only link between a ground region and the plane.
  - Corner-to-corner (Euclidean) spacing is not checked; the foundry's DRC deck must still be run.
- Metal connectivity (LVS-lite, reported, not gated): **not clean.** The ground between the arms of each qubit is enclosed by the arm and T-bar cut-outs (arm lengths below about 132 µm seal the inner corners) and floats; so does ground around the buses and in pockets near the resonators and coupler flux lines. Shorted lines that end in such a pocket do not reach the plane; the qubit flux lines end open (no DC return); the qubit SQUID junctions end in etch on their ground side. The extraction treats every enclosed ground fragment as grounded. The counts are in `alpha1_v11_review.json → layout.connectivity`.
- All 17 resonators are drawn at their λ/4 length (4.17–4.35 mm).
- Previews: `AF17_v11_layout.png`, `AF17_v11_detail.png`.
- Independent reader: `tools/gds_preview.py`.

## 6. Simulation gate (ROADMAP.md)

| Criterion | Threshold | Result |
|---|---|---|
| d=3 logical error per round | ≤ 0.6 % | 0.374 % (memory-Z), 0.368 % (memory-X) ✔ |
| Logical lifetime vs best T1 | > 2× | 3.3× ✔ |
| CZ coherent error / leakage, every pair | ≤ 5e-4 / ≤ 2e-4 | largest pair 3.6e-4 / 8.4e-5 (means 8.1e-5 / 4.4e-5) ✔ |
| Idle ZZ at park | ≤ 10 kHz | 3.4–4.5 kHz, mean 3.9 kHz (4.3 kHz re-nulled post-layout) ✔ |
| Layout | placement/spacing checks clean; E_C within 3 %, couplings within 6 %, exchange re-nullable (mask DRC and connectivity reported, not gated) | clean; 0.05 %; 4.6 %; yes ✔ (mask DRC 0 + 4, connectivity not clean) |
| Fabrication Monte Carlo, 2 % junction spread | ≥ 90 % collision-free | 98 % (v1.0 review; same frequency plan) ✔ |

## 7. Before tape-out (Alpha fabrication checklist)

1. Full-wave/FEM cross-check of the qubit, coupler and readout cells (Palace / Sonnet / Q3D). Then a second synthesis iteration on the residual errors (up to 4.6 % on the qubit–coupler couplings).
2. Interposer design: feedlines, Purcell filters, edge launchers, wiring of the 64 bump pads.
3. Airbridges across the coupler and flux-line trenches; dicing lanes, alignment marks, test structures (junction test arrays, resonator Q).
4. Junction process: target Rn per qubit from the design JSON, with ≤ 2 % spread (or trimming).
5. Calibration plan: idle frequencies; coupler exchange-null park per pair; the 3-knob CZ; leakage-removing reset tune-up.
6. Characterisation: follow `CHARACTERISATION-PLAN.md`. It sets the order of measurements and the data file, and `qpuforge hwgate` turns that file into the hardware-gate verdict of ROADMAP.md, including model versus measurement. The most important single number is the CZ residual control error.

## 8. Assumptions and limitations

- The CZ residual control allowance (1e-3) is assumed, not simulated, and it dominates the CZ error.
- Pulse simulation is closed-system. Decoherence during the gate is applied by the emulators; 1/f flux noise is applied analytically.
- Leakage is a frame-level model. The decoder is leakage-unaware.
- Capacitances come from a quasi-static extraction on a thick substrate. Not included: box modes, kinetic inductance, crosstalk, participation-ratio TLS loss.
- The mask is not fabrication-ready (§5): four DRC space violations, floating ground, no flux-line return, qubit junctions not landing; the coupling of the XY and flux lines to the qubits is not modelled. Control wiring and feedlines on the interposer are not drawn.
- Coherence comes from the device model (Q = 3e6 dielectric loss, Purcell, 2 µΦ0 flux noise, photon shot noise), not from measurement.

## 9. Files

| File | Content |
|---|---|
| `alpha1_v11.json` | Chip design (qubits, couplings, coupler element, leakage spec) |
| `alpha1_v11_review.json/.txt` | Full review output |
| `layout_cells.json` | Synthesised cell geometry, achieved vs target circuit, synthesis log |
| `AF17_v11.gds` | GDSII layout |
| `AF17_v11_layout.{png,svg}`, `AF17_v11_detail.{png,svg}` | Previews rendered from the GDSII by an independent reader |
| `CHARACTERISATION-PLAN.md` | What to measure on the fabricated chip, the data-file format, and how the hardware gate is evaluated |
| `example_characterisation_synthetic.json`, `example_gate_report_synthetic.json` | A virtual chip's data file and its gate report (synthetic; for rehearsing the tool chain) |
| `../../results/cz_cache.json` | Pulse-level CZ calibrations (coherent part) |
