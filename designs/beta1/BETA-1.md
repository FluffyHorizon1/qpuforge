# Beta-1 (BF-49) — design specification

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

| | |
|---|---|
| **Status** | Rev C. Simulation review complete (9 of 9 simulation-gate criteria pass). Layout (GDSII) passes the mask width/space DRC (§6). |
| **Purpose** | Show that the logical error falls with code size: a distance-5 rotated surface code on 49 qubits |
| **Builds on** | Alpha-1 v1.1 (AF-17): same architecture and frequency plan; cells re-synthesised for the stronger readout coupling |
| **Tooling** | qpuforge v0.4: `beta1_review` (995 s, 100,000 shots per point) |

All performance numbers are **simulation results**. The assumptions are in §8.

## 1. Design: Alpha-1 v1.1 + the upgrades the reviews ranked highest

| Area | Alpha-1 v1.1 | Beta-1 | Alpha-1 v1.1 sensitivity | How |
|---|---|---|---|---|
| Size | 17 qubits, d = 3 | **49 qubits (25 data + 24 measure), 80 couplers, d = 5** | — | same cell topology, 900 µm pitch |
| Dielectric loss | Q = 3e6 | **Q = 6e6** (T1 ≈ ×2) | −37 % | tantalum base layer, improved substrate cleaning / surface treatment |
| Effective temperature | 45 mK | **35 mK** | −27 % | IR-tight sample box, Eccosorb filtering, better line thermalisation |
| Readout | 500 ns, κ = 4 MHz | **300 ns, κ = 2χ = 8 MHz** | −25 % | stronger readout coupling behind **individual Purcell filters** (suppression 1000) |
| Readout plan | 3 feedlines, up to 6 resonators each | **6 feedlines, up to 9 resonators each, 6.90–7.38 GHz** | — | one TWPA + HEMT chain per feedline |
| CZ residual control error | 1e-3 | **5e-4** | −20 % | flux-pulse predistortion (step-response calibration), closed-loop recalibration |
| Leakage | leakage-removing reset | **leakage-removing reset** (DQLR evaluated, not adopted) | — | re-tested at d = 5 in long experiments (§3) |
| Junction spread | 2 % | **0.5 %** (trimmed) | yield | post-fabrication laser annealing of junctions |

**Review history.**

Rev A → rev B (first review):

1. **Faster readout cancelled the materials upgrade.** With κ = 8 MHz behind the per-feedline Purcell filters, Purcell decay became the T1 limit: mean T1 was 62 µs, below Alpha-1. Individual Purcell filters (suppression 1000) restore a dielectric-limited T1.
2. **DQLR did not pay.** Data-qubit leakage removal raised the 50-round logical failure at d = 5. It is an option, not the baseline.
3. **Two placement violations in the 49-qubit layout.** Two resonator coupling pads ran into neighbouring coupler buses. The pads now leave the serpentine upwards, and the placement/spacing checks are clean. The merged-etch width/space check added in v0.3 still reports violations (§6).

Rev B → rev C (second review, with qpuforge v0.4):

4. **Too many resonators per feedline.** Three feedlines for 49 qubits put 17 resonators on each line across almost 1 GHz. Readout chains are run at 6–8 per line today. Rev C uses 6 feedlines with at most 9 resonators each, in 6.90–7.38 GHz. Cost: 6 more coaxial lines.
5. **The reused Alpha-1 cells could not deliver the Beta-1 readout.** Beta-1 doubles the dispersive shift, which needs a qubit–readout coupling of 229 MHz. The Alpha-1 cell gives about 160 MHz (9.4 fF where 13.3 fF is needed): half the intended shift, and a slower or noisier readout than the model assumed. Rev B passed because its layout gate checked design rules only. Rev C re-synthesises the cells (readout pad 40.9 → 67.4 µm, arms shortened to hold E_C) and the gate now checks the extracted circuit.
6. **Decoding above 22 defects is now exact.** Rev B's decoder fell back to greedy matching there, and its long runs suffered for it. On identical syndromes of a 50-round distance-5 run (4,000 shots, 19.3 defects per shot, 29 % of shots above 22), the fallback fails in 5.3 % of shots (214 failures) and exact matching in 0.8 % (32). Rev B's long-run numbers were pessimistic; the leakage decision is re-tested with the exact decoder in §3.

**Device model (rev C):** mean T1 147 µs (best 157 µs), T2 201 µs, readout error 0.55 % at 300 ns.
CZ: 80 pairs, all calibrated at pulse level. Mean 43.6 ns, mean error 9.12e-04 including decoherence,
leakage 4.4e-05 per gate (largest pair 8.4e-05).

## 2. Error suppression with code distance

Memory experiments in both bases, Pauli-frame Monte Carlo, exact matching decoder. Logical error per round, fitted
over the round counts shown:

| | d = 3 (5–20 rounds) | d = 5 (5–20 rounds) | d = 7 (9 rounds, 100,000 shots) |
|---|---|---|---|
| **Beta-1, memory-Z** | 0.0753 % | **0.0131 %** | 0.0022 % |
| **Beta-1, memory-X** | 0.0744 % | **0.0140 %** | — |
| Alpha-1 v1.1 physics, memory-Z (reference) | 0.3585 % | 0.1238 % | — |

- **Λ(3→5) = 5.76 ± 0.57 (Z) and 5.32 ± 0.49 (X).** The ± values are 1σ Monte Carlo errors.
- Λ(5→7) = 7.81 ± 1.96, comparing d = 5 (0.0174 % per round, 50,000 shots) and d = 7 (20 failures) at the same 9 rounds. A 9-round run is short for distance 7. The distance-7 chip is reviewed in full, at 14 rounds and in both bases, as Gamma-1 (`../gamma1/GAMMA-1.md`).
- **Λ is a matched-model best case.** The decoder's error model comes from the same noise model that generates the errors. Its priors are therefore exact: no drift, no calibration mismatch and no unmodelled correlated errors. Hardware Λ will be lower.
- Distance-5 round 734 ns → logical lifetime **5.61 ms**, i.e. 35.7× the best physical T1.
- Against Alpha-1 v1.1 at distance 3, the Beta-1 distance-5 logical qubit is **27.4× better**.

## 3. Long experiments and leakage (d = 5)

| | P_L after 50 rounds | leaked data / measure | round |
|---|---|---|---|
| **Leakage-removing reset (Beta-1 baseline)** | 0.0073 | 5.8e-05 / 8.8e-05 | 728 ns |
| + DQLR every round (option) | 0.0086 | 6.1e-06 / 7.9e-05 | 758 ns |

**Decision confirmed:** the leakage-removing reset alone does at least as well over 50 rounds (the two differ by 2.4 standard errors). DQLR clears data-qubit leakage, but its time and gate cost are not repaid at this leakage rate, so it stays an option.

**Stability over 100 rounds** (the simulation counterpart of the hardware gate "no growth of the detection-event rate"):

| Reset | Detection-rate drift, first 10 vs last 10 rounds | P_L after 100 rounds | Defects per shot (largest) |
|---|---|---|---|
| **Leakage-removing (Beta-1)** | +0.6 % | 0.0146 | 38.1 (104) |
| Plain (no leakage removal) | +37.7 % | 0.0736 | 47.7 (268) |

## 4. Error budget and remaining sensitivity (Beta-1 physics at d = 3)

Leave-one-out (eps change when the source is removed): T1 -76 %, control -57 %, Tphi -21 %, leakage -14 %, readout -10 %, zz -7 %, reset -5 %, thermal +1 %.

Next-step sensitivity: dielectric Q 6e6 → 1.2e7 -38 %, readout 300 → 200 ns -20 %, CZ residual 5e-4 → 2.5e-4 -17 %, readout separation 0.5 % → 0.25 % +0 %.

## 5. Fabrication yield (d = 5, frequency-collision rules)

| Junction resistance spread | Collision-free chips | Frequency scatter (rms) |
|---|---|---|
| 2 % (untrimmed) | 192/200 (96 %) | 57 MHz |
| 1 % | 185/200 (92 %) | 29 MHz |
| **0.5 % (trimmed)** | **199/200 (100 %)** | 14 MHz |

## 6. Wiring and layout

- **Lines:** 49 XY + 49 qubit flux + 80 coupler flux + 6 readout in + 6 readout out = **190 coaxial lines**.
  - 6 feedlines with 8.2 resonators each on average, 6.90–7.38 GHz in 60 MHz steps (spacing/κ = 7.5).
  - The wiring needs high-density flex cabling and a large-format dilution refrigerator. Each readout output chain is a TWPA plus a HEMT.
- **Cells:** `layout_cells.json`, synthesised against the Beta-1 circuit targets by method-of-moments extraction.
  - Extracted against target: E_C within 0.1 %, couplings within 5.5 %, qubit–readout capacitance -0.0 % (13.29 fF; readout g 229 MHz for a design value of 229 MHz).
  - Parked exchange re-nulled by parking the couplers at 8.417 GHz (residual ZZ -4.36 kHz).
- **Layout:** `BF49.gds` (cell `BF49_BETA1`).
  - 9.4 × 9.4 mm, 41,573 etch polygons, 129 junction regions, 227 bump pads.
  - DRC: placement and island-spacing checks clean. Minimum width/space on the merged etch layer (4 µm / 4 µm, projection metric): clean.
- **Interposer** (feedlines, Purcell filters, fan-out to edge launchers): designed in Beta; not drawn here.

## 7. Beta simulation gate (ROADMAP.md)

| Criterion | Result |
| --- | --- |
| Λ(3→5) ≥ 2 in both bases | ✔ pass |
| logical error per round at d = 5 ≤ 0.1 % in both bases | ✔ pass |
| d = 5 at least 3× below Alpha-1 v1.1 at d = 3 | ✔ pass |
| d = 5 logical lifetime > 5× best physical T1 | ✔ pass |
| trimmed fabrication yield ≥ 80 % collision-free | ✔ pass |
| chosen leakage scheme at least as good as the alternative over 50 rounds | ✔ pass |
| detection rate flat over 100 rounds with the leakage-removing reset (drift < 5 %) | ✔ pass |
| layout passes placement, island-spacing and mask width/space DRC (4 µm / 4 µm on the merged etch layer, see §6) | ✔ pass |
| extracted E_C within 3 %, couplings and readout capacitance within 6 %; parked exchange re-nullable | ✔ pass |

## 8. Assumptions and limitations

- The upgrades are **targets**: the Beta-1 numbers assume they are achieved (Q = 6e6, 35 mK, 300 ns readout at 0.5 % separation error, CZ residual 5e-4). Alpha hardware characterisation must confirm the baseline first. The sensitivity table shows what each target is worth.
- All Alpha-1 v1.1 limitations apply, as listed in `../alpha1-v11/ALPHA-1-v1.1.md` §8: CZ residual allowance, frame-level leakage, quasi-static capacitances, no foundry DRC deck yet.
- **`residual_ctrl_err` is an assumed allowance, not a simulation result.** It covers pulse distortion, calibration drift, TLS defects and crosstalk, none of which are simulated: 5e-4 per CZ for Beta-1 and 1e-3 for Alpha-1 v1.1. It is added on top of each pair's simulated coherent, flux-noise and decoherence error. The 5e-4 target assumes predistortion and closed-loop recalibration work as planned (§1).
- **Leakage is modelled, not decoded.** The frame simulator tracks leakage, but the decoder gets no leakage herald: there is no 3-state readout, no erasure flags and no reweighting. A leaked qubit's outcomes reach the decoder as ordinary detection events. A leakage-aware decoder is future work (ROADMAP Phase Gamma).
- **Λ is a matched-model best case** (§2): the decoder is built from the true error model.
- **Mask DRC.** The layout passes the placement and island-spacing checks and the mask-level check (minimum width and space of 4 µm on the merged etch layer, projection metric) with no violations, after 392 ground slivers narrower than the minimum were merged into the neighbouring cut-outs (`layout.hpp`, `fill_ground_slivers`; the extraction cells get the same clean-up). Corner-to-corner (Euclidean) spacing is not checked: the foundry's DRC deck must still be run before tape-out.
- **CZ leakage.** The largest per-pair leakage, 8.4e-05 per gate, is within the 2e-4 limit of the Alpha simulation gate (ROADMAP.md). Beta-1 uses the same 16 qubit-frequency combinations as Alpha-1 v1.1, so this is the pair class reported there (`../alpha1-v11/ALPHA-1-v1.1.md` §6).
- The frequency-collision yield uses the v0.2 collision rules with static (sweet-spot) frequencies. Tunable qubits and per-pair CZ calibration give extra margin that the rules do not credit.
- The couplers park 83 MHz below their sweet spot to null the exchange left by the layout. Their flux noise there is not in the error-correction model.
- The cell synthesis meets the layout tolerances with little margin on the qubit–coupler couplings: it fits their mean, and the measure and data arms differ in length, so one lies above target and the other below. Junction sizes are set from the extracted E_C before tape-out.

## 9. Files

| File | Content |
|---|---|
| `beta1.json` | Chip design (49 qubits, 80 couplers) |
| `beta1_review.json/.txt` | Full review output |
| `layout_cells.json` | Cell geometry, synthesis log, extracted circuit |
| `BF49.gds`, `BF49_layout.png` | GDSII layout and preview |
| `../../results/cz_cache.json` | Pulse-level CZ calibrations, shared with Alpha-1 v1.1. Model version 2, generated by the review runs |
| `../../results/PROVENANCE.md` | Source commit, compiler, seeds and date of the committed review outputs |
