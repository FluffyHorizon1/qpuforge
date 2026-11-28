# Gamma-1 (GF-97) — design specification

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

| | |
|---|---|
| **Status** | Memory review complete (9 of 9 memory-gate criteria pass). Layout (GDSII) is DRC- and connectivity-clean in the model (2D LVS-lite connectivity with modelled airbridges; airbridges, EM and inductance not verified); control-line T1 budget **not met** (flux line); foundry DRC not yet run (§7). Logical operations not yet reviewed (§9). |
| **Purpose** | A distance-7 logical memory on 97 qubits: the third point on the suppression curve, and the first chip where decoding and correlated errors need their own answers |
| **Builds on** | Beta-1 rev C (BF-49): same physics, same frequency plan, same synthesised cells |
| **Tooling** | qpuforge v0.4: `gamma1_review` (2205 s; 100,000 shots per point, twice that at distance 7) |

All performance numbers are **simulation results**. The assumptions are in §10.

## 1. Design

Gamma-1 changes nothing in the physics of Beta-1. It is the same cell, repeated to distance 7.

| | Beta-1 (BF-49) | Gamma-1 (GF-97) |
|---|---|---|
| Code | distance 5 | **distance 7** |
| Qubits / couplers | 49 / 80 | **97 (49 data + 48 measure) / 168** |
| Readout | 6 feedlines, up to 9 resonators each | **11 feedlines, up to 9 resonators each**, 6.90–7.38 GHz |
| Coaxial lines | 190 | **384** |
| Chip | 9.4 × 9.4 mm | **13.0 × 13.0 mm** |
| Decoder | exact matching | exact matching, any number of defects (new in v0.4) |

**Device model:** mean T1 147 µs (best 157 µs), T2 202 µs, readout error 0.55 % at 300 ns.
CZ: 168 pairs, all calibrated at pulse level. Mean 43.3 ns, mean error 9.08e-04 including decoherence,
leakage 5.0e-05 per gate (largest pair 1.2e-04). Frequency collisions: 0.

## 2. Error suppression with code distance

Memory experiments of 14 rounds in both bases, Pauli-frame Monte Carlo, exact matching. The logical error per
round is taken from the failure probability of that one run length, the same way at every distance.

| | d = 3 | d = 5 | d = 7 |
|---|---|---|---|
| Memory-Z (\|0⟩, bit flips) | 0.0843 % | 0.0137 % | **0.00318 %** |
| Memory-X (\|+⟩, phase flips) | 0.0832 % | 0.0145 % | **0.00264 %** |
| Failures / shots (Z; X) | 1167 / 100,000; 1152 / 100,000 | 191 / 100,000; 202 / 100,000 | 89 / 200,000; 74 / 200,000 |

- **Λ(5→7) = 4.30 ± 0.55 (Z) and 5.47 ± 0.74 (X).** Λ(3→5) = 6.17 ± 0.48 (Z) and 5.75 ± 0.44 (X). The ± values are 1σ Monte Carlo errors.
- Distance-7 round 746 ns → logical lifetime **25.6 ms**, 163× the best physical T1.
- The statistical error on the distance-7 values is 12 % (one standard error, relative). It dominates the error on Λ(5→7).
- Λ is a matched-model best case: the decoder is given the same noise model that generated the errors (`../beta1/BETA-1.md` §8).

## 3. Decoding

Minimum-weight perfect matching, exact at any number of defects: a dynamic programme up to 8 defects and Edmonds'
blossom algorithm beyond.

- Distance 7, 14 rounds: **11.5 defects per shot** (largest 37); 1.53 % of shots have more than 22.
- Software decoding takes about **2.54 µs per round** on one CPU core, against a round of 0.75 µs. The
  figure moves with the machine and its load. On average it is the same order as the chip's own pace; the worst case
  is far slower, and a batch decoder is not a streaming one. Real-time decoding stays an FPGA task (§9).
- **What exactness buys.** Before v0.4 the decoder fell back to greedy matching above 22 defects. On a 40-round
  distance-7 run (4,000 shots, same syndromes, 31.9 defects per shot, 90 % above 22):

| Decoder | Logical failure probability | Failures | Decode time per round |
|---|---|---|---|
| Exact matching (v0.4) | **0.0008** | 3 | 7.5 µs |
| Greedy above 22 defects (v0.3) | 0.0805 | 322 | 839 µs |

## 4. High-energy impacts

A cosmic-ray or gamma impact fills the substrate with phonons and the qubits with quasiparticles. Error rates rise
across the whole chip at once and relax over milliseconds, which breaks the assumption of independent errors that the
code relies on. The model adds an excess decay or dephasing rate to every qubit, lets it relax exponentially, and
leaves the decoder with its quiet-chip priors.

| Scenario | Impact rate | Recovery | Peak excess | P(logical flip per impact) | Unprotected for | Visible in the syndrome for | Floor per round | Share of the d = 7 error |
|---|---|---|---|---|---|---|---|---|
| No mitigation (T1 burst; memory-Z) | 0.1 /s | 25 ms | T1 → 1 µs | 0.50 | 147 ms | 103 ms | 3.7e-08 | 0.13 % |
| Gap-engineered junctions (dephasing burst; memory-X) | 0.1 /s | 1 ms | Tφ → 0.1 µs | 0.50 | 8.4 ms | 6.1 ms | 3.7e-08 | 0.13 % |

- **An impact costs the logical qubit.** The burst outlasts thousands of rounds, so the logical state is randomised
  whatever the distance. The floor is therefore set by the impact rate alone: rate × round time ÷ 2.
- **It does not limit Gamma-1.** At distance 7 the floor is 0.13 % of the logical error. At this Λ, impacts
  would start to dominate from about distance 15: before the distances a useful machine needs.
- **It is detectable.** The detection-event count stays above 1.5× its quiet value for 103 ms of the
  147 ms unprotected window. A control system that watches that rate can flag the window and discard or
  restart what ran in it. That costs about 1.5 % of running time and turns a silent logical error into a known one.
- **Gap engineering shortens the exposure,** from 147 ms to 8.4 ms per impact in this model. It does
  not remove the logical flip: a millisecond of strong dephasing is still more than a thousand rounds.

Sources for the scenarios are in §10.

## 5. Long experiment (d = 7, 40 rounds)

| Reset | P_L after 40 rounds (failures / shots) | Detection rate, round 2 → round 40 | Drift (first 5 vs last 5 rounds) | Leaked data / measure |
|---|---|---|---|---|
| **Leakage-removing (Gamma-1)** | 0.00110 (55 / 50,000) | 0.0330 → 0.0319 | -1.1 % | 6.8e-05 / 1.0e-04 |
| Plain (no leakage removal) | 0.02246 (1123 / 50,000) | 0.0330 → 0.0410 | +23.0 % | 1.3e-03 / 1.8e-03 |

Without leakage removal the detection rate climbs through the run and the logical failure probability is 20×
higher. The rising detection rate is the signature to look for on hardware.

## 6. Fabrication yield (97 qubits, frequency-collision rules)

| Junction resistance spread | Collision-free chips | Frequency scatter (rms) |
|---|---|---|
| 1 % | 91/100 (91 %) | 29 MHz |
| **0.5 % (trimmed)** | **100/100 (100 %)** | 14 MHz |

## 7. Wiring and layout

- **Lines:** 97 XY + 97 qubit flux + 168 coupler flux + 11 readout in + 11 readout out = **384 coaxial lines**.
  - 11 feedlines with 8.8 resonators each on average, 6.90–7.38 GHz in 60 MHz steps.
  - Beta-1 needs 190 lines and Alpha-1 64. This count is past what coaxial wiring in one refrigerator carries comfortably: Gamma-1 needs high-density flex lines or multiplexed flux control. That is a programme risk, not a chip risk.
- **Layout:** `GF97.gds` (cell `GF97_GAMMA1`), built from the Beta-1 rev C cells.
  - 13.0 × 13.0 mm, 67,134 etch polygons, 265 junction regions, 459 bump pads.
  - DRC: placement and island-spacing checks clean. Minimum width/space on the merged etch layer (4 µm / 4 µm, projection metric): clean.
  - Metal connectivity: **clean**: 0 floating ground fragments; 362 of 362 shorted lines reach the ground plane (qubit flux lines 97 of 97); 265 of 265 junctions join an island to the plane; 0 conductors share metal with another island; 774 ground straps (modelled airbridges; 265 over lines, 315 over coupler or bus strips), 0 small ground fragments removed.
  - Control lines (reported, not gated): XY line 0.030 / 0.030 fF (measure / data qubit), flux line 0.032 / 0.032 fF effective (2.37 / 2.32 fF if it were open); T1 limits at least 1.56 ms (XY) and 1.32 ms (flux), 11.4× and 9.7× the design T1 of the worst qubit (budget 10×: **not met**).
  - Extracted against target: E_C within 0.1 %, couplings within 5.4 %, qubit–readout capacitance +0.4 %.
  - Parked exchange re-nulled by parking the couplers at 8.415 GHz (residual ZZ -4.31 kHz).

## 8. Gamma memory gate

| Criterion | Result |
| --- | --- |
| Λ(5→7) ≥ 2 in both bases | ✔ pass |
| distance-7 logical error per round ≤ 0.01 % in both bases | ✔ pass |
| distance-7 logical lifetime > 50× best physical T1 | ✔ pass |
| logical-error floor from impacts below 10 % of the distance-7 error (both scenarios) | ✔ pass |
| detection rate flat over 40 rounds with the leakage-removing reset (drift < 5 %) | ✔ pass |
| trimmed fabrication yield ≥ 80 % collision-free | ✔ pass |
| layout passes placement, island-spacing and mask width/space checks (§7) | ✔ pass |
| metal connectivity clean: conductors separate, shorted lines and all ground on the plane, SQUIDs island-to-plane (§7) | ✔ pass |
| extracted E_C within 3 %, couplings and readout capacitance within 6 %; parked exchange re-nullable | ✔ pass |

## 9. What Phase Gamma still needs

This review covers the distance-7 memory. The rest of the phase is open:

- **Logical operations.** Two distance-5 patches with lattice surgery (a logical CNOT). qpuforge needs multi-patch
  circuits with merge and split schedules, and detectors that follow them.
- **Leakage-aware decoding.** Three-state readout heralds a leaked measure qubit. Feeding that to the decoder should
  recover part of the leakage cost in the budget.
- **Real-time decoding.** A streaming decoder under 1 µs per round on an FPGA. §3 gives the software reference.
- **Impact veto.** A detector on the detection-event rate and a restart protocol (§4).

## 10. Assumptions and limitations

- Every Beta-1 assumption applies (`../beta1/BETA-1.md` §8): the upgrades over Alpha-1 are targets until Alpha hardware confirms the baseline.
- **Mask.** The layout passes the placement and island-spacing checks and the mask-level width/space check (4 µm on the merged etch layer, projection metric) with no violations; the metal connectivity check is clean with its ground straps; the control-line T1 budget (reported, not gated) is **not met**: XY line 0.030 / 0.030 fF (measure / data qubit), flux line 0.032 / 0.032 fF effective (2.37 / 2.32 fF if it were open); T1 limits at least 1.56 ms (XY) and 1.32 ms (flux), 11.4× and 9.7× the design T1 of the worst qubit (budget 10×: **not met**) (§7). Not verified: the airbridges (their capacitance to the coupler and bus strips they cross is not in the extraction cells), EM (full-wave) extraction, inductance (the flux lines' coupling to the SQUIDs), a foundry DRC deck.
- **CZ calibrations** use pulse model version 4: four levels each for the two qubits and the coupler. The 168 pairs fall into 16 frequency combinations, the same 16 as on Alpha-1 v1.1 and Beta-1. The largest per-pair leakage, 1.2e-04 per gate, is within the 2e-4 limit of the Alpha simulation gate (ROADMAP.md).
- **Impact model.** Each scenario is sampled at a ladder of severities, each half a decade below the last, until the excess logical error is lost in the statistical noise; below the last point the excess is taken as proportional to the severity. Runs in which the logical qubit is randomised within the 14 rounds count as fully failed. The unmitigated scenario follows McEwen et al., *Nature Physics* 18, 107 (2022): chip-wide loss of T1 by one to two orders of magnitude, recovering with a time constant of about 25 ms. The rate of 0.1 per second is a round figure for a centimetre-scale chip; it scales with substrate area and shielding and must be measured on the actual chip. The gap-engineered scenario follows McEwen et al., *Phys. Rev. Lett.* 133, 240601 (2024) and Kurilovich et al., arXiv:2506.18228 (2025): no T1 degradation, but frequency shifts of up to a few MHz lasting about 1 ms, modelled here as a dephasing burst. The model treats the burst as chip-wide and uniform. Gap engineering does not address the scrambling of two-level-system defects by an impact (same PRL); that is not modelled either.
- **Estimator.** One run length per distance gives the logical error per round without a fit; it folds the time-boundary effects of the first and last rounds into the value. The same estimator at every distance keeps Λ consistent.
- **Decode time** is a software figure per decoding thread on the build machine, measured while the other OpenMP threads ran the same experiment. It covers matching only, not the syndrome extraction from the simulator.
- The frequency-collision yield uses the static (sweet-spot) collision rules.

## 11. Files

| File | Content |
|---|---|
| `gamma1.json` | Chip design (97 qubits, 168 couplers) |
| `gamma1_review.json/.txt` | Full review output |
| `GF97.gds`, `GF97_layout.png` | GDSII layout and preview |
| `../beta1/layout_cells.json` | Cell geometry (shared with Beta-1 rev C) |
| `../../results/cz_cache.json` | Pulse-level CZ calibrations (shared) |
