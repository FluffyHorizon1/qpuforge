# Gamma-1 (GF-97) — design specification

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

| | |
|---|---|
| **Status** | Memory review complete ({{GATE_SUMMARY}}). Layout (GDSII) {{MASK_STATUS}} (§7). Logical operations are simulated outside this review (lattice surgery at d = 3 and 5, `../../results/surgery_study.txt`) but not reviewed on this chip (§9). |
| **Purpose** | A distance-7 logical memory on 97 qubits: the third point on the suppression curve, and the first chip where decoding and correlated errors need their own answers |
| **Builds on** | Beta-1 rev C (BF-49): same physics, same frequency plan, same synthesised cells |
| **Tooling** | qpuforge v0.4: `gamma1_review` ({{RUNTIME}} s; {{SHOTS}} shots per point, twice that at distance 7) |

All performance numbers are **simulation results**. The assumptions are in §10.

## 1. Design

Gamma-1 changes nothing in the physics of Beta-1. It is the same cell, repeated to distance 7.

| | Beta-1 (BF-49) | Gamma-1 (GF-97) |
|---|---|---|
| Code | distance 5 | **distance 7** |
| Qubits / couplers | 49 / 80 | **97 (49 data + 48 measure) / 168** |
| Readout | 6 feedlines, up to 9 resonators each | **{{FEED}} feedlines, up to 9 resonators each**, {{RO_LO}}–{{RO_HI}} GHz |
| Coaxial lines | {{BETA_LINES}} | **{{W_TOT}}** |
| Chip | 9.4 × 9.4 mm | **{{CHIP_W}} × {{CHIP_H}} mm** |
| Decoder | exact matching | exact matching, any number of defects (new in v0.4) |

**Device model:** mean T1 {{T1_MEAN}} µs (best {{T1_BEST}} µs), T2 {{T2_MEAN}} µs, readout error {{RO_ERR}} % at 300 ns.
CZ: {{N_CZ}} pairs, all calibrated at pulse level. Mean {{TCZ}} ns, mean error {{CZ_ERR}} including decoherence,
leakage {{CZ_LEAK}} per gate (largest pair {{CZ_LEAK_MAX}}). Frequency collisions: {{COLLISIONS}}.

## 2. Error suppression with code distance

Memory experiments of {{ROUNDS}} rounds in both bases, Pauli-frame Monte Carlo, exact matching. The logical error per
round is taken from the failure probability of that one run length, the same way at every distance.

| | d = 3 | d = 5 | d = 7 |
|---|---|---|---|
| Memory-Z (\|0⟩, bit flips) | {{EZ3}} | {{EZ5}} | **{{EZ7}}** |
| Memory-X (\|+⟩, phase flips) | {{EX3}} | {{EX5}} | **{{EX7}}** |
| Failures / shots (Z; X) | {{F3}} | {{F5}} | {{F7}} |

- **Λ(5→7) = {{L57Z}} (Z) and {{L57X}} (X).** Λ(3→5) = {{L35Z}} (Z) and {{L35X}} (X). The ± values are 1σ Monte Carlo errors.
- Distance-7 round {{ROUND7}} ns → logical lifetime **{{LIFE7}} ms**, {{LIFE_X}}× the best physical T1.
- The statistical error on the distance-7 values is {{E7ERR}} (one standard error, relative). It dominates the error on Λ(5→7).
- Λ is a matched-model best case: the decoder is given the same noise model that generated the errors (`../beta1/BETA-1.md` §8).

## 3. Decoding

Minimum-weight perfect matching, exact at any number of defects: a dynamic programme up to 8 defects and Edmonds'
blossom algorithm beyond.

- Distance 7, {{ROUNDS}} rounds: **{{DEF_MEAN}} defects per shot** (largest {{DEF_MAX}}); {{DEF_OVER}} of shots have more than 22.
- Software decoding takes about **{{DEC_US}} µs per round** on one CPU core, against a round of {{ROUND7_US}} µs. The
  figure moves with the machine and its load. On average it is the same order as the chip's own pace; the worst case
  is far slower, and a batch decoder is not a streaming one. Real-time decoding stays an FPGA task (§9).
- **What exactness buys.** Before v0.4 the decoder fell back to greedy matching above 22 defects. On a {{LONG_R}}-round
  distance-7 run ({{LONG_N}} shots, same syndromes, {{LONG_DEF}} defects per shot, {{LONG_OVER}} above 22):

| Decoder | Logical failure probability | Failures | Decode time per round |
|---|---|---|---|
| Exact matching (v0.4) | **{{LONG_PL_EX}}** | {{LONG_F_EX}} | {{LONG_US_EX}} µs |
| Greedy above 22 defects (v0.3) | {{LONG_PL_LG}} | {{LONG_F_LG}} | {{LONG_US_LG}} µs |

## 4. High-energy impacts

A cosmic-ray or gamma impact fills the substrate with phonons and the qubits with quasiparticles. Error rates rise
across the whole chip at once and relax over milliseconds, which breaks the assumption of independent errors that the
code relies on. The model adds an excess decay or dephasing rate to every qubit, lets it relax exponentially, and
leaves the decoder with its quiet-chip priors.

| Scenario | Impact rate | Recovery | Peak excess | P(logical flip per impact) | Unprotected for | Visible in the syndrome for | Floor per round | Share of the d = 7 error |
|---|---|---|---|---|---|---|---|---|
| No mitigation (T1 burst; memory-Z) | {{B1_RATE}} /s | {{B1_TAU}} ms | T1 → {{B1_T1}} µs | {{B1_PF}} | {{B1_DEAD}} ms | {{B1_VIS}} ms | {{B1_FLOOR}} | {{B1_SHARE}} |
| Gap-engineered junctions (dephasing burst; memory-X) | {{B2_RATE}} /s | {{B2_TAU}} ms | Tφ → {{B2_TP}} µs | {{B2_PF}} | {{B2_DEAD}} ms | {{B2_VIS}} ms | {{B2_FLOOR}} | {{B2_SHARE}} |

- **An impact costs the logical qubit.** The burst outlasts thousands of rounds, so the logical state is randomised
  whatever the distance. The floor is therefore set by the impact rate alone: rate × round time ÷ 2.
- **It does not limit Gamma-1.** At distance 7 the floor is {{B1_SHARE}} of the logical error. At this Λ, impacts
  would start to dominate from about distance {{B1_DSTAR}}: before the distances a useful machine needs.
- **It is detectable.** The detection-event count stays above 1.5× its quiet value for {{B1_VIS}} ms of the
  {{B1_DEAD}} ms unprotected window. A control system that watches that rate can flag the window and discard or
  restart what ran in it. That costs about {{B1_LOST}} of running time and turns a silent logical error into a known one.
- **Gap engineering shortens the exposure,** from {{B1_DEAD}} ms to {{B2_DEAD}} ms per impact in this model. It does
  not remove the logical flip: a millisecond of strong dephasing is still more than a thousand rounds.

Sources for the scenarios are in §10.

## 5. Long experiment (d = 7, 40 rounds)

| Reset | P_L after 40 rounds (failures / shots) | Detection rate, round 2 → round 40 | Drift (first 5 vs last 5 rounds) | Leaked data / measure |
|---|---|---|---|---|
| **Leakage-removing (Gamma-1)** | {{LR_PL}} ({{LR_F}}) | {{LR_D0}} → {{LR_D1}} | {{LR_DRIFT}} | {{LR_LD}} / {{LR_LM}} |
| Plain (no leakage removal) | {{PR_PL}} ({{PR_F}}) | {{PR_D0}} → {{PR_D1}} | {{PR_DRIFT}} | {{PR_LD}} / {{PR_LM}} |

Without leakage removal the detection rate climbs through the run and the logical failure probability is {{PR_RATIO}}×
higher. The rising detection rate is the signature to look for on hardware.

## 6. Fabrication yield (97 qubits, frequency-collision rules)

| Junction resistance spread | Collision-free chips | Frequency scatter (rms) |
|---|---|---|
| 1 % | {{Y10}} | {{S10}} MHz |
| **0.5 % (trimmed)** | **{{Y05}}** | {{S05}} MHz |

## 7. Wiring and layout

- **Lines:** {{W_XY}} XY + {{W_ZQ}} qubit flux + {{W_ZC}} coupler flux + {{W_RI}} readout in + {{W_RO}} readout out = **{{W_TOT}} coaxial lines**.
  - {{FEED}} feedlines with {{QPF}} resonators each on average, {{RO_LO}}–{{RO_HI}} GHz in 60 MHz steps.
  - Beta-1 needs {{BETA_LINES}} lines and Alpha-1 {{ALPHA_LINES}}. This count is past what coaxial wiring in one refrigerator carries comfortably: Gamma-1 needs high-density flex lines or multiplexed flux control. That is a programme risk, not a chip risk.
- **Layout:** `GF97.gds` (cell `GF97_GAMMA1`), built from the Beta-1 rev C cells.
  - {{CHIP_W}} × {{CHIP_H}} mm, {{ETCH}} etch polygons, {{JJ}} junction regions, {{BUMPS}} bump pads.
  - DRC: {{DRC}}.
  - Metal connectivity (reported, not gated): {{CONN}}.
  - Extracted against target: E_C within {{EC_ERR}}, couplings within {{G_ERR}}, qubit–readout capacitance {{RO_C_ERR}}.
  - Parked exchange re-nulled by parking the couplers at {{PARK}} GHz (residual ZZ {{PARK_ZZ}} kHz).

## 8. Gamma memory gate

{{GATES}}

## 9. What Phase Gamma still needs

This review covers the distance-7 memory. Each other part of the phase now has software outside this review, which
the review does not use; what is still open is named in each item:

- **Logical operations.** Lattice surgery is simulated (`surgery.hpp`, `qpuforge surgery`,
  `../../docs/design/lattice-surgery.md`): logical ZZ and XX measurement between two patches and the three-patch CNOT,
  at d = 3 and 5 (`../../results/surgery_study.txt`). With Beta-1 physics (also Gamma-1's) at d = 5 the ZZ outcome is
  wrong in 0.14 % of measurements and a CNOT on |00⟩ gets some output wrong in 0.71 % of runs. Open: a chip layout for
  the patches and their routing region, real-time feed-forward of the corrections, a joint decoder for all observables.
- **Leakage-aware decoding.** Implemented (`HeraldSpec`, `qpuforge qec --herald`,
  `../../docs/design/heralded-leakage.md`): three-state readout heralds a leaked measure qubit and the matching is
  re-weighted. At d = 7 with the leakage-removing reset and exact matching, heralds lower the logical error by 10 % at
  the chip's own CZ leakage and by 46 % at ten times that (`../../results/decoder_study.txt` §3).
- **Real-time decoding.** In software, union-find (`qpuforge qec --decoder uf`) gives 0.98–1.10× the logical error of
  exact matching at d = 3 to 7. Its sliding-window version (`--decoder uf-stream`) comes within 1.16× of batch
  union-find with a 10-round window committing 5, and takes 1.5–1.8 µs per committed round at d = 7, against a round of
  0.73 µs in that 50-round run (`../../results/decoder_study.txt` §2, §5). §3 gives the exact-matching reference. Open:
  a decoder under 1 µs per round on an FPGA, and a faster, pipelined window scheme.
- **Impact veto.** Implemented in software (`qpuforge qec --bursts --veto`, `../../docs/design/impact-veto.md`): at
  d = 7 a 3-round × radius-2 window catches 90–95 % of impacts and removes 98–99.5 % of the failures in impacted
  shots. Shot-by-shot discard leaves the floor of §4 at 3.7e-08 per round; holding the chip for 7 recovery times (T1
  burst) or 10 (dephasing burst) after a veto brings it below what the study resolves, at 1.8 % and 0.1 % dead time
  (`../../results/decoder_study.txt` §4). Open: the restart protocol on hardware.

## 10. Assumptions and limitations

- Every Beta-1 assumption applies (`../beta1/BETA-1.md` §8): the upgrades over Alpha-1 are targets until Alpha hardware confirms the baseline.
- {{MASK_NOTE}}
- **CZ calibrations** use pulse model version {{CZ_VERSION}}: four levels each for the two qubits and the coupler. The {{N_CZ}} pairs fall into 16 frequency combinations, the same 16 as on Alpha-1 v1.1 and Beta-1. {{CZ_LEAK_NOTE}}
- **Impact model.** Each scenario is sampled at a ladder of severities, each half a decade below the last, until the excess logical error is lost in the statistical noise; below the last point the excess is taken as proportional to the severity. Runs in which the logical qubit is randomised within the {{ROUNDS}} rounds count as fully failed. The unmitigated scenario follows McEwen et al., *Nature Physics* 18, 107 (2022): chip-wide loss of T1 by one to two orders of magnitude, recovering with a time constant of about 25 ms. The rate of 0.1 per second is a round figure for a centimetre-scale chip; it scales with substrate area and shielding and must be measured on the actual chip. The gap-engineered scenario follows McEwen et al., *Phys. Rev. Lett.* 133, 240601 (2024) and Kurilovich et al., arXiv:2506.18228 (2025): no T1 degradation, but frequency shifts of up to a few MHz lasting about 1 ms, modelled here as a dephasing burst. The model treats the burst as chip-wide and uniform. Gap engineering does not address the scrambling of two-level-system defects by an impact (same PRL); that is not modelled either.
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
