# Alpha-1 v1.1 (AF-17) — design specification

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

| | |
|---|---|
| **Status** | Simulation gate **passed**. Layout first cut (GDSII), DRC clean. Ready for an external EM cross-check and fabrication planning. |
| **Supersedes** | Alpha-1 v1.0 (`designs/alpha1/ALPHA-1.md`, qpuforge v0.2) |
| **Purpose** | First error-correction-ready chip: 17-qubit distance-3 rotated surface code (memory experiment) |
| **Tooling** | qpuforge v0.3: `alpha1_v11_review` (1,305 s), `qpuforge gds`, `qpuforge cz` |

All performance numbers are **simulation results** under the assumptions in §8.

## 1. What changed from v1.0

| Area | v1.0 | v1.1 |
|---|---|---|
| Tunable coupler | Abstract (g_max = 40 MHz, 5 kHz residual ZZ) | **Physical coupler transmon** per pair: f_max 8.5 GHz (park at the sweet spot), α_c −130 MHz, g_qc 160 MHz |
| Direct coupling | Not modelled | **Per-pair g12 = 9.4–9.8 MHz**, chosen so the parked coupler cancels the qubit–qubit exchange: no avoided crossings when a measure qubit sweeps past its other neighbours |
| CZ gate | 23.3 ns, 2e-3 control error (assumed) | **Pulse-level**: sin² flux ramps (7 ns) + flat top, 3-knob calibration per pair. 35.6–49.2 ns, simulated coherent error 7.5e-5, leakage 5.2e-5 |
| Leakage | Folded into control error | **Simulated**: CZ leakage of the moving qubit, leaked-partner phase scrambling, transport (5 %), seepage at 2/T1 |
| Reset | Active reset 200 ns, 0.3 % | Same, made **leakage-removing** (multi-level reset: \|2⟩ → \|0⟩) |
| Layout | Floor plan (SVG) | **MoM-synthesised cells, assembled chip, DRC, GDSII**, post-layout circuit check |

## 2. Chip summary

| | |
|---|---|
| Qubits | 17 flux-tunable transmons (SQUID), all idling at their sweet spot: 9 data (5.11–5.29 GHz), 8 measure (5.81–5.99 GHz); α ≈ −220 MHz (−215 MHz post-layout) |
| Couplers | 24 grounded tunable-coupler transmons (one per surface-code CZ pair) |
| Readout | 17 λ/4 resonators at 6.90–7.20 GHz (60 MHz steps), κ = 4 MHz, 2χ = κ, one Purcell filter per feedline, 3 feedlines; 500 ns readout, 0.5 % separation error |
| Reset | Active, 200 ns, 0.3 % residual excitation, leakage-removing |
| Coherence (model) | T1 65–76 µs (dielectric Q = 3e6, Purcell), Tφ 120–170 µs, T_eff 45 mK |
| Chip | 5.8 × 5.8 mm, qubit pitch 900 µm |
| Lines | 17 XY + 17 qubit flux + 24 coupler flux + 3 readout in + 3 readout out = **64**, all ending on flip-chip bump pads |

## 3. CZ gate (pulse level)

**Protocol.** The measure qubit (the higher frequency) is the *mover*. For each CZ:

1. Its flux pulse brings |11⟩ into resonance with |20⟩ (f_mover + α = f_data). The interaction point is 5.33–5.50 GHz, an excursion of 310–665 MHz.
2. At the same time the coupler swings from 8.5 GHz down to ≈ 6.03–6.13 GHz, which turns on an effective coupling with a 36 MHz |11⟩–|20⟩ gap.
3. One full |11⟩ → |20⟩ → |11⟩ cycle gives the conditional phase π.

Three knobs are calibrated per pair:
- the mover's interaction flux;
- the flat-top length;
- the coupler bias, which also closes the off-resonant |01⟩–|10⟩ swap.

The simulation covers:
- exact transmon charge-basis spectra for each element, truncated to 4/3/4 levels;
- a smooth sign gauge for the charge matrix elements;
- piecewise-exact propagation with dt = 0.1 ns (converged to 6e-7 in F);
- the single-qubit rotating frame, so idle ZZ counts as a physical error.

| Quantity (24 pairs) | Value |
|---|---|
| Duration | 35.6–49.2 ns (mean 43.9 ns) |
| Coherent error (1 − F_avg, closed system) | 7.5e-5 mean |
| Leakage per CZ (moving qubit) | 5.2e-5 mean, 1.3e-4 max |
| 1/f flux noise (local dephasing at the excursion + quasi-static sensitivity) | 3.2e-5 |
| Residual idle ZZ at park | 2.4–3.3 kHz |
| Spectator crossings (Landau–Zener, nulled exchange) | < 1e-15 |
| Residual control allowance (assumed: distortion, drift, TLS, crosstalk) | 1.0e-3 |
| **Total CZ error incl. T1/Tφ during the gate** | **1.86e-3** |

The per-pair table is `alpha1_v11_review.json → cz_calibration.per_pair`. The calibration cache is `results/cz_cache.json`.

## 4. Error-correction performance (distance-3 memory-Z, 5–25 rounds, Pauli-frame Monte Carlo + exact MWPM)

| | Logical error per round | Round |
|---|---|---|
| v1.0 model | 0.445 % ± 0.005 % | 847 ns |
| **v1.1 model (pulse CZ + leakage + leakage-removing reset)** | **0.370 % ± 0.003 %** | 928 ns |
| Logical lifetime | **251 µs** = 3.3 × best physical T1 (75.6 µs) | |
| Distance-5 projection (same physics, 49 qubits) | 0.128 % → **Λ = 2.89** | |

v1.1 is better than v1.0 despite the slower CZ. The physically modelled control error (≈ 1.06e-3) is about half the
2e-3 that v1.0 assumed.

**Leakage handling** (eps over 5–25 rounds; P_L after 50 rounds):

| | eps / round | P_L (50 rounds) | leaked data / measure |
|---|---|---|---|
| Leakage not modelled | 0.365 % | 0.152 | — |
| No leakage removal | 0.409 % | 0.168 | 1.1e-3 / 1.5e-3 |
| **Leakage-removing reset (chosen)** | **0.379 %** | **0.155** | 6.9e-5 / 1.1e-4 |
| + DQLR each round (30 ns, 5e-4 error) | 0.399 % | 0.167 | 7.3e-6 / 9.7e-5 |
| Stress, 10× CZ leakage: LRR / LRR + DQLR | — | 0.190 / 0.197 | 5.0e-4 / 5.9e-5 |

**Decision:** the leakage-removing reset is enough for Alpha-1. DQLR removes data-qubit leakage but costs more than it
saves at distance 3, even at 10× leakage. It is carried to Beta as an option to re-test at distance 5.

**Error budget** (eps change when one source is removed): T1 −78 %, CZ/1q control −40 %, Tφ −31 %, readout −16 %,
reset −9 %, leakage −4 %, thermal −2 %, ZZ −1 %.

**Sensitivity** (single change):

| Change | eps change |
|---|---|
| Dielectric Q ×2 | −34 % |
| T_eff 35 mK | −25 % |
| Readout 300 ns | −22 % |
| CZ residual 5e-4 | −16 % |
| Readout separation 0.25 % | −7 % |
| Flux noise ×3 | +10 % |
| CZ residual 2e-3 | +43 % |

These ranks set the Beta-1 upgrade list.

## 5. Physical layout

**Cells.** An Xmon-style cross with a T-bar at each arm end faces each coupler's T-bar across a 5 µm gap. The readout
claw wraps a pad in the cross's inner corner. The coupler is a grounded transmon: two T-bars, a 10 µm body and a square
mid pad; its SQUID sits on the pad edge with an inductive flux line. A thin floating strip runs parallel to each coupler.
The synthesis drove its T-bar to ~0, so it only adds 1.5 fF of loading and could be dropped in v1.2.

The dimensions were synthesised by Gauss-Newton against the circuit targets, using method-of-moments extraction
(quasi-static, ε_eff = (11.45 + 1)/2 on silicon) with a coarse-to-fine mesh correction:

| Cell | Dimension | Value |
|---|---|---|
| Measure / data qubit | arm length (centre to T-bar edge) | 129.3 / 131.9 µm |
| | arm width / gap / T-bar | 16 / 30 / 120 × 12 µm |
| | readout pad / claw gap | 40.9 / 5 µm |
| Coupler | T-bar (facing the qubit) | 43.5 µm long, 5 µm gap |
| | body / trench / mid pad | 10 / 20 / 252.6 × 252.6 µm |

| Extracted (fine mesh) | Value | Target | Error |
|---|---|---|---|
| E_C measure / data | 198.2 / 195.8 MHz | 202.2 / 199.7 MHz | −2.0 % |
| E_C coupler | 124.4 MHz | 125.6 MHz | −0.9 % |
| G measure–coupler / coupler–data (GHz per n n) | 0.0598 / 0.0579 | 0.0576 / 0.0612 | +3.8 % / −5.3 % |
| G measure–data (direct) | 0.00568 | 0.00555 | +2.2 % |
| C qubit–readout | 9.36 fF | 9.15 fF | +2.3 % |
| C_Σ qubit / C coupler / C_qc / C_qq | 98.0–99.2 / 157 / 5.72 / 0.105 fF | | |

**Post-layout circuit check.** The junctions are retargeted to the design frequencies and the representative pair is
rebuilt from the extracted circuit:
- α = −215 MHz; coupler α −129 MHz; g_qc = 168 MHz; readout g = 158 MHz (design 156).
- g12 = 9.94 MHz against the 9.54 MHz needed for an exchange-free park at the sweet spot. This would leave a 0.77 MHz parked gap (1.5e-4 transfer per spectator crossing).
- Because g12 is *above* the null value, parking each coupler 205 MHz below its sweet spot (8.295 GHz, df/dΦ = 4.2 GHz/Φ0) nulls it exactly. ZZ is then −3.7 kHz.
- Calibration plan: per-pair park flux = exchange null.

**GDSII** (`AF17_v11.gds`, cell `AF17_V11`, 1 nm database unit):

| Layer | Content |
|---|---|
| 1 | Etch: ground-plane cut-outs minus metal islands (14,380 polygons) |
| 2 | Josephson junction / SQUID regions (41) |
| 3 | Flip-chip bump pads (75): XY, qubit flux, coupler flux, resonator–feedline coupling |
| 10 | Chip outline |
| 20 | Labels |

- DRC: 10 µm metal spacing, 4 µm at designed coupling gaps — **clean**.
- All 17 resonators are drawn at their λ/4 length (4.17–4.35 mm).
- Previews: `AF17_v11_layout.png`, `AF17_v11_detail.png`.
- Independent reader: `tools/gds_preview.py`.

## 6. Simulation gate (ROADMAP.md)

| Criterion | Threshold | Result |
|---|---|---|
| d=3 logical error per round | ≤ 0.6 % | 0.370 % ✔ |
| Logical lifetime vs best T1 | > 2× | 3.3× ✔ |
| CZ coherent error / leakage | ≤ 5e-4 / ≤ 2e-4 | 7.5e-5 / 5.2e-5 (max 1.3e-4) ✔ |
| Idle ZZ at park | ≤ 10 kHz | 2.8 kHz (3.7 kHz re-nulled post-layout) ✔ |
| Layout | DRC clean; E_C within 3 %, couplings within 6 %, exchange re-nullable | clean; 2.0 %; 5.3 %; yes ✔ |
| Fabrication Monte Carlo, 2 % junction spread | ≥ 90 % collision-free | 98 % (v1.0 review; same frequency plan) ✔ |

## 7. Before tape-out (Alpha fabrication checklist)

1. Full-wave/FEM cross-check of the qubit, coupler and readout cells (Palace / Sonnet / Q3D). Then a second synthesis iteration on the residual 2–5 % errors.
2. Interposer design: feedlines, Purcell filters, edge launchers, wiring of the 64 bump pads.
3. Airbridges across the coupler and flux-line trenches; dicing lanes, alignment marks, test structures (junction test arrays, resonator Q).
4. Junction process: target Rn per qubit from the design JSON, with ≤ 2 % spread (or trimming).
5. Calibration plan: idle frequencies; coupler exchange-null park per pair; the 3-knob CZ; leakage-removing reset tune-up.
6. Characterisation plan to measure the parameters behind the hardware gate in ROADMAP.md. The most important single number is the CZ residual control error.

## 8. Assumptions and limitations

- The CZ residual control allowance (1e-3) is assumed, not simulated, and it dominates the CZ error.
- Pulse simulation is closed-system. Decoherence during the gate is applied by the emulators; 1/f flux noise is applied analytically.
- Leakage is a frame-level model. The decoder is leakage-unaware.
- Capacitances come from a quasi-static extraction on a thick substrate. Not included: box modes, kinetic inductance, crosstalk, participation-ratio TLS loss.
- The layout is a first cut. Control wiring and feedlines on the interposer are not drawn.
- Coherence comes from the device model (Q = 3e6 dielectric loss, Purcell, 2 µΦ0 flux noise, photon shot noise), not from measurement.

## 9. Files

| File | Content |
|---|---|
| `alpha1_v11.json` | Chip design (qubits, couplings, coupler element, leakage spec) |
| `alpha1_v11_review.json/.txt` | Full review output |
| `layout_cells.json` | Synthesised cell geometry, achieved vs target circuit, synthesis log |
| `AF17_v11.gds` | GDSII layout |
| `AF17_v11_layout.png`, `AF17_v11_detail.png` | Previews rendered from the GDSII by an independent reader |
| `../../results/cz_cache.json` | Pulse-level CZ calibrations (coherent part) |
