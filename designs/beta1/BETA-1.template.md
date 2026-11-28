# Beta-1 (BF-49) — design specification

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

| | |
|---|---|
| **Status** | Rev C. Simulation review complete ({{GATE_SUMMARY}}). Layout (GDSII) {{MASK_STATUS}} (§6). |
| **Purpose** | Show that the logical error falls with code size: a distance-5 rotated surface code on 49 qubits |
| **Builds on** | Alpha-1 v1.1 (AF-17): same architecture and frequency plan; cells re-synthesised for the stronger readout coupling |
| **Tooling** | qpuforge v0.4: `beta1_review` ({{RUNTIME}} s, {{SHOTS}} shots per point) |

All performance numbers are **simulation results**. The assumptions are in §8.

## 1. Design: Alpha-1 v1.1 + the upgrades the reviews ranked highest

| Area | Alpha-1 v1.1 | Beta-1 | Alpha-1 v1.1 sensitivity | How |
|---|---|---|---|---|
| Size | 17 qubits, d = 3 | **49 qubits (25 data + 24 measure), 80 couplers, d = 5** | — | same cell topology, 900 µm pitch |
| Dielectric loss | Q = 3e6 | **Q = 6e6** (T1 ≈ ×2) | {{A_SENS_Q}} | tantalum base layer, improved substrate cleaning / surface treatment |
| Effective temperature | 45 mK | **35 mK** | {{A_SENS_T}} | IR-tight sample box, Eccosorb filtering, better line thermalisation |
| Readout | 500 ns, κ = 4 MHz | **300 ns, κ = 2χ = 8 MHz** | {{A_SENS_RO}} | stronger readout coupling behind **individual Purcell filters** (suppression 1000) |
| Readout plan | 3 feedlines, up to 6 resonators each | **{{FEED}} feedlines, up to 9 resonators each, {{RO_LO}}–{{RO_HI}} GHz** | — | one TWPA + HEMT chain per feedline |
| CZ residual control error | 1e-3 | **5e-4** | {{A_SENS_CZ}} | flux-pulse predistortion (step-response calibration), closed-loop recalibration |
| Leakage | leakage-removing reset | **leakage-removing reset** (DQLR evaluated, not adopted) | — | re-tested at d = 5 in long experiments (§3) |
| Junction spread | 2 % | **0.5 %** (trimmed) | yield | post-fabrication laser annealing of junctions |

**Review history.**

Rev A → rev B (first review):

1. **Faster readout cancelled the materials upgrade.** With κ = 8 MHz behind the per-feedline Purcell filters, Purcell decay became the T1 limit: mean T1 was 62 µs, below Alpha-1. Individual Purcell filters (suppression 1000) restore a dielectric-limited T1.
2. **DQLR did not pay.** Data-qubit leakage removal raised the 50-round logical failure at d = 5. It is an option, not the baseline.
3. **Two placement violations in the 49-qubit layout.** Two resonator coupling pads ran into neighbouring coupler buses. The pads now leave the serpentine upwards, and the placement/spacing checks are clean. The merged-etch width/space check added in v0.3 still reported violations at rev B; the current rev C layout passes it (§6).

Rev B → rev C (second review, with qpuforge v0.4):

4. **Too many resonators per feedline.** Three feedlines for 49 qubits put 17 resonators on each line across almost 1 GHz. Readout chains are run at 6–8 per line today. Rev C uses {{FEED}} feedlines with at most 9 resonators each, in {{RO_LO}}–{{RO_HI}} GHz. Cost: {{EXTRA_LINES}} more coaxial lines.
5. **The reused Alpha-1 cells could not deliver the Beta-1 readout.** Beta-1 doubles the dispersive shift, which needs a qubit–readout coupling of {{RO_G_T}} MHz. The Alpha-1 cell gives about 160 MHz (9.4 fF where 13.3 fF is needed): half the intended shift, and a slower or noisier readout than the model assumed. Rev B passed because its layout gate checked design rules only. Rev C re-synthesises the cells (readout pad {{RO_PAD_A}} → {{RO_PAD_B}} µm in the current cells, arms shortened to hold E_C) and the gate now checks the extracted circuit.
6. **Decoding above 22 defects is now exact.** Rev B's decoder fell back to greedy matching there, and its long runs suffered for it. On identical syndromes of a {{DEC_R}}-round distance-5 run ({{DEC_N}} shots, {{DEC_DEF}} defects per shot, {{DEC_OVER}} of shots above 22), the fallback fails in {{DEC_PL_LG}} of shots ({{DEC_F_LG}} failures) and exact matching in {{DEC_PL_EX}} ({{DEC_F_EX}}). Rev B's long-run numbers were pessimistic; the leakage decision is re-tested with the exact decoder in §3.

**Device model (rev C):** mean T1 {{T1_MEAN}} µs (best {{T1_BEST}} µs), T2 {{T2_MEAN}} µs, readout error {{RO_ERR}} % at 300 ns.
CZ: {{N_CZ}} pairs, all calibrated at pulse level. Mean {{TCZ}} ns, mean error {{CZ_ERR}} including decoherence,
leakage {{CZ_LEAK}} per gate (largest pair {{CZ_LEAK_MAX}}).

## 2. Error suppression with code distance

Memory experiments in both bases, Pauli-frame Monte Carlo, exact matching decoder. Logical error per round, fitted
over the round counts shown:

| | d = 3 (5–20 rounds) | d = 5 (5–20 rounds) | d = 7 (9 rounds, {{N7}} shots) |
|---|---|---|---|
| **Beta-1, memory-Z** | {{EPS3}} | **{{EPS5}}** | {{EPS7}} |
| **Beta-1, memory-X** | {{EPS3X}} | **{{EPS5X}}** | — |
| Alpha-1 v1.1 physics, memory-Z (reference) | {{A_EPS3}} | {{A_EPS5}} | — |

- **Λ(3→5) = {{L35}} (Z) and {{L35X}} (X).** The ± values are 1σ Monte Carlo errors.
- Λ(5→7) = {{L57}}, comparing d = 5 ({{EPS5_9}} per round, {{N5_9}} shots) and d = 7 ({{F7}} failures) at the same 9 rounds. A 9-round run is short for distance 7. The distance-7 chip is reviewed in full, at 14 rounds and in both bases, as Gamma-1 (`../gamma1/GAMMA-1.md`).
- **Λ is a matched-model best case.** The decoder's error model comes from the same noise model that generates the errors. Its priors are therefore exact: no drift, no calibration mismatch and no unmodelled correlated errors. Hardware Λ will be lower.
- Distance-5 round {{ROUND5}} ns → logical lifetime **{{LIFE5}} ms**, i.e. {{LIFE_X}}× the best physical T1.
- Against Alpha-1 v1.1 at distance 3, the Beta-1 distance-5 logical qubit is **{{BEAT}}× better**.

## 3. Long experiments and leakage (d = 5)

| | P_L after 50 rounds | leaked data / measure | round |
|---|---|---|---|
| **Leakage-removing reset (Beta-1 baseline)** | {{LR_PL}} | {{LR_LD}} / {{LR_LM}} | {{LR_RN}} ns |
| + DQLR every round (option) | {{LD_PL}} | {{LD_LD}} / {{LD_LM}} | {{LD_RN}} ns |

{{DQLR_DECISION}}

**Stability over 100 rounds** (the simulation counterpart of the hardware gate "no growth of the detection-event rate"):

| Reset | Detection-rate drift, first 10 vs last 10 rounds | P_L after 100 rounds | Defects per shot (largest) |
|---|---|---|---|
| **Leakage-removing (Beta-1)** | {{ST_L_DRIFT}} | {{ST_L_PL}} | {{ST_L_DEF}} |
| Plain (no leakage removal) | {{ST_P_DRIFT}} | {{ST_P_PL}} | {{ST_P_DEF}} |

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

- **Lines:** {{W_XY}} XY + {{W_ZQ}} qubit flux + {{W_ZC}} coupler flux + {{W_RI}} readout in + {{W_RO}} readout out = **{{W_TOT}} coaxial lines**.
  - {{FEED}} feedlines with {{QPF}} resonators each on average, {{RO_LO}}–{{RO_HI}} GHz in 60 MHz steps (spacing/κ = 7.5).
  - The wiring needs high-density flex cabling and a large-format dilution refrigerator. Each readout output chain is a TWPA plus a HEMT.
- **Cells:** `layout_cells.json`, synthesised against the Beta-1 circuit targets by method-of-moments extraction.
  - Extracted against target: E_C within {{EC_ERR}}, couplings within {{G_ERR}}, qubit–readout capacitance {{RO_C_ERR}} ({{RO_C}} fF; readout g {{RO_G}} MHz for a design value of {{RO_G_T}} MHz).
  - Parked exchange re-nulled by parking the couplers at {{PARK}} GHz (residual ZZ {{PARK_ZZ}} kHz).
- **Layout:** `BF49.gds` (cell `BF49_BETA1`).
  - {{CHIP_W}} × {{CHIP_H}} mm, {{ETCH}} etch polygons, {{JJ}} junction regions, {{BUMPS}} bump pads.
  - DRC: {{DRC}}.
  - Metal connectivity: {{CONN}}.
  - Control lines (reported, not gated): {{LINES}}.
- **Interposer** (feedlines, Purcell filters, fan-out to edge launchers): designed in Beta; not drawn here.

## 7. Beta simulation gate (ROADMAP.md)

{{GATES}}

## 8. Assumptions and limitations

- The upgrades are **targets**: the Beta-1 numbers assume they are achieved (Q = 6e6, 35 mK, 300 ns readout at 0.5 % separation error, CZ residual 5e-4). Alpha hardware characterisation must confirm the baseline first. The sensitivity table shows what each target is worth.
- All Alpha-1 v1.1 limitations apply, as listed in `../alpha1-v11/ALPHA-1-v1.1.md` §8: CZ residual allowance, frame-level leakage, quasi-static capacitances, no foundry DRC deck yet.
- **`residual_ctrl_err` is an assumed allowance, not a simulation result.** It covers pulse distortion, calibration drift, TLS defects and crosstalk, none of which are simulated: 5e-4 per CZ for Beta-1 and 1e-3 for Alpha-1 v1.1. It is added on top of each pair's simulated coherent, flux-noise and decoherence error. The 5e-4 target assumes predistortion and closed-loop recalibration work as planned (§1).
- **Leakage is modelled, not decoded, in this review.** The frame simulator tracks leakage, but the review's decoder gets no leakage herald, so a leaked qubit's outcomes reach it as ordinary detection events. Leakage-aware decoding now exists for Phase Gamma (3-state readout heralds the leaked measure qubits and the matching is re-weighted: `HeraldSpec`, `qpuforge qec --herald`, `../../docs/design/heralded-leakage.md`); this review does not use it.
- **Λ is a matched-model best case** (§2): the decoder is built from the true error model.
- {{MASK_NOTE}}
- **CZ leakage.** {{CZ_LEAK_NOTE}} Beta-1 uses the same 16 qubit-frequency combinations as Alpha-1 v1.1, so this is the pair class reported there (`../alpha1-v11/ALPHA-1-v1.1.md` §6).
- The frequency-collision yield uses the v0.2 collision rules with static (sweet-spot) frequencies. Tunable qubits and per-pair CZ calibration give extra margin that the rules do not credit.
- The couplers park {{PARK_BELOW}} MHz below their sweet spot to null the exchange left by the layout. Their flux noise there is not in the error-correction model.
- The cell synthesis meets the layout tolerances with little margin on the qubit–coupler couplings: it fits their mean, and the measure and data arms differ in length, so one lies above target and the other below. Junction sizes are set from the extracted E_C before tape-out.

## 9. Files

| File | Content |
|---|---|
| `beta1.json` | Chip design (49 qubits, 80 couplers) |
| `beta1_review.json/.txt` | Full review output |
| `layout_cells.json` | Cell geometry, synthesis log, extracted circuit |
| `BF49.gds`, `BF49_layout.png` | GDSII layout and preview |
| `../../results/cz_cache.json` | Pulse-level CZ calibrations, shared with Alpha-1 v1.1 and Gamma-1. CZ model version 4, generated by `qpuforge cz` (`../../results/PROVENANCE.md`); the review reuses it (0 simulated, 16 cached) |
| `../../results/PROVENANCE.md` | Source commit, compiler, seeds and date of the committed review outputs |
