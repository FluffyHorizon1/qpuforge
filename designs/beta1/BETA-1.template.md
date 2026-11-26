# Beta-1 (BF-49) — design specification

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

| | |
|---|---|
| **Status** | Simulation review complete ({{GATE_SUMMARY}}). Layout first cut (GDSII). |
| **Purpose** | Show that the logical error falls with code size: a distance-5 rotated surface code on 49 qubits |
| **Builds on** | Alpha-1 v1.1 (AF-17): same architecture, same frequency plan and synthesised cells |
| **Tooling** | qpuforge v0.3: `beta1_review` ({{RUNTIME}} s, {{SHOTS}} shots per point) |

All performance numbers are **simulation results**. The assumptions are in §7.

## 1. Design: Alpha-1 v1.1 + the upgrades the reviews ranked highest

| Area | Alpha-1 v1.1 | Beta-1 | Alpha-1 v1.1 sensitivity | How |
|---|---|---|---|---|
| Size | 17 qubits, d = 3 | **49 qubits (25 data + 24 measure), 80 couplers, d = 5** | — | same cells, 900 µm pitch |
| Dielectric loss | Q = 3e6 | **Q = 6e6** (T1 ≈ ×2) | −34 % | tantalum base layer, improved substrate cleaning / surface treatment |
| Effective temperature | 45 mK | **35 mK** | −25 % | IR-tight sample box, Eccosorb filtering, better line thermalisation |
| Readout | 500 ns, κ = 4 MHz | **300 ns, κ = 2χ = 8 MHz** | −22 % | stronger readout coupling behind **individual Purcell filters** (suppression 1000) |
| CZ residual control error | 1e-3 | **5e-4** | −16 % | flux-pulse predistortion (step-response calibration), closed-loop recalibration |
| Leakage | leakage-removing reset | **leakage-removing reset** (DQLR evaluated, not adopted) | — | re-tested at d = 5 in long experiments (§3) |
| Junction spread | 2 % | **0.5 %** (trimmed) | yield | post-fabrication laser annealing of junctions |

**Review history.** The first review of this design (rev A) caught three problems, all fixed in this revision (rev B):

1. **Faster readout cancelled the materials upgrade.** With κ = 8 MHz behind the per-feedline Purcell filters, Purcell decay became the T1 limit: mean T1 was 62 µs, below Alpha-1. Individual Purcell filters (suppression 1000) restore a dielectric-limited T1 of {{T1_MEAN}} µs.
2. **DQLR did not pay.** Data-qubit leakage removal raised the 50-round logical failure from 0.147 to 0.166 at d = 5. It is now an option, not the baseline.
3. **Two DRC violations in the 49-qubit layout.** Two resonator coupling pads ran into neighbouring coupler buses. The pads now leave the serpentine upwards, and the placement/spacing checks are clean. The merged-etch width/space check added in v0.3 still reports violations (§6).

Rev A still met the distance-scaling gate despite its T1 problem (d=3 0.224 %, d=5 0.067 %, Λ = 3.34).

**Device model (rev B):** mean T1 {{T1_MEAN}} µs (best {{T1_BEST}} µs), T2 {{T2_MEAN}} µs, readout error {{RO_ERR}} % at 300 ns.
CZ: {{N_CZ}} pairs, all calibrated at pulse level. Mean {{TCZ}} ns, mean error {{CZ_ERR}} including decoherence,
leakage {{CZ_LEAK}} per gate.

## 2. Error suppression with code distance

Memory-Z experiments, Pauli-frame Monte Carlo, exact MWPM decoding. Logical error per round, fitted over the round
counts shown:

| | d = 3 (5–20 rounds) | d = 5 (5–20 rounds) | d = 7 (9 rounds, {{N7}} shots) |
|---|---|---|---|
| **Beta-1** | {{EPS3}} | **{{EPS5}}** | {{EPS7}} |
| Alpha-1 v1.1 physics (reference) | {{A_EPS3}} | {{A_EPS5}} | — |

- **Λ(3→5) = {{L35}}.** Λ(5→7) = {{L57}}, comparing d = 5 ({{EPS5_9}} per round, {{N5_9}} shots) and d = 7 at the same 9 rounds. The ± values are 1σ Monte Carlo errors. At d = 7 the decoder is exact up to 22 defects and greedy beyond, which can only add logical errors at d = 7.
- **Λ is a matched-model best case.** The decoder's error model comes from the same noise model that generates the errors. Its priors are therefore exact: no drift, no calibration mismatch and no unmodelled correlated errors. Hardware Λ will be lower.
- Distance-5 round {{ROUND5}} ns → logical lifetime **{{LIFE5}} ms**, i.e. {{LIFE_X}}× the best physical T1.
- Against Alpha-1 v1.1 at distance 3, the Beta-1 distance-5 logical qubit is **{{BEAT}}× better**.

## 3. Long experiments and leakage (d = 5, 50 rounds)

| | P_L after 50 rounds | leaked data / measure | round |
|---|---|---|---|
| **Leakage-removing reset (Beta-1 baseline)** | {{LR_PL}} | {{LR_LD}} / {{LR_LM}} | {{LR_RN}} ns |
| + DQLR every round (option) | {{LD_PL}} | {{LD_LD}} / {{LD_LM}} | {{LD_RN}} ns |

{{DQLR_DECISION}}

## 4. Error budget and remaining sensitivity (Beta-1 physics at d = 3)

Leave-one-out (eps change when the source is removed): {{BUDGET}}.

Next-step sensitivity: {{SENS}}.

## 5. Fabrication yield (d = 5, frequency-collision rules)

| Junction resistance spread | Collision-free chips | Frequency scatter (rms) |
|---|---|---|
| 2 % (untrimmed) | {{Y20}} | {{S20}} MHz |
| 1 % | {{Y10}} | {{S10}} MHz |
| **0.5 % (trimmed)** | **{{Y05}}** | {{S05}} MHz |

## 6. Wiring and layout

- **Lines:** {{W_XY}} XY + {{W_ZQ}} qubit flux + {{W_ZC}} coupler flux + {{W_RI}} readout in + {{W_RO}} readout out = **{{W_TOT}} coaxial lines**, {{FEED}} feedlines (~{{QPF}} qubits per feedline).
  - Readout frequencies are unique on each feedline: 6.90–7.86 GHz in 60 MHz steps, with spacing/κ = 7.5.
  - The wiring needs high-density flex cabling and a large-format dilution refrigerator. Each readout output chain is a TWPA plus a HEMT.
  - Phase Gamma will need more feedlines.
- **Layout:** `BF49.gds` (cell `BF49_BETA1`), built from the synthesised Alpha-1 v1.1 cells. The circuit targets are the same because the frequency plan is unchanged.
  - {{CHIP_W}} × {{CHIP_H}} mm, {{ETCH}} etch polygons, {{BUMPS}} bump pads.
  - DRC: {{DRC}}.
- **Interposer** (feedlines, Purcell filters, fan-out to edge launchers): designed in Beta; not drawn here.

## 7. Beta simulation gate (ROADMAP.md)

{{GATES}}

## 8. Assumptions and limitations

- The upgrades are **targets**: the Beta-1 numbers assume they are achieved (Q = 6e6, 35 mK, 300 ns readout at 0.5 % separation error, CZ residual 5e-4). Alpha hardware characterisation must confirm the baseline first. The sensitivity table shows what each target is worth.
- All Alpha-1 v1.1 limitations apply, as listed in `../alpha1-v11/ALPHA-1-v1.1.md` §8: CZ residual allowance, frame-level leakage, quasi-static capacitances, first-cut layout.
- **`residual_ctrl_err` is an assumed allowance, not a simulation result.** It covers pulse distortion, calibration drift, TLS defects and crosstalk, none of which are simulated: 5e-4 per CZ for Beta-1 and 1e-3 for Alpha-1 v1.1. It is added on top of each pair's simulated coherent, flux-noise and decoherence error. The 5e-4 target assumes predistortion and closed-loop recalibration work as planned (§1).
- **Leakage is modelled, not decoded.** The frame simulator tracks leakage, but the decoder gets no leakage herald: there is no 3-state readout, no erasure flags and no reweighting. A leaked qubit's outcomes reach the decoder as ordinary detection events. A leakage-aware decoder is future work (ROADMAP Phase Gamma).
- **Λ is a matched-model best case** (§2): the decoder is built from the true error model.
- The frequency-collision yield uses the v0.2 collision rules with static (sweet-spot) frequencies. Tunable qubits and per-pair CZ calibration give extra margin that the rules do not credit.
- d = 7 numbers come from the frame model with a decoder that falls back to greedy matching beyond 22 defects.

## 9. Files

| File | Content |
|---|---|
| `beta1.json` | Chip design (49 qubits, 80 couplers) |
| `beta1_review.json/.txt` | Full review output |
| `BF49.gds`, `BF49_layout.png` | GDSII layout and preview |
| `../../results/cz_cache.json` | Pulse-level CZ calibrations, shared with Alpha-1 v1.1. Model version 2, generated by the review runs |
| `../../results/PROVENANCE.md` | Source commit, compiler, seeds and date of the committed review outputs |
