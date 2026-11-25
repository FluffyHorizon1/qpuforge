# Beta-1 (BF-49) — design specification

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

| | |
|---|---|
| **Status** | Simulation review complete (7 of 7 simulation-gate criteria pass). Layout first cut (GDSII). |
| **Purpose** | Show that the logical error falls with code size: a distance-5 rotated surface code on 49 qubits |
| **Builds on** | Alpha-1 v1.1 (AF-17): same architecture, same frequency plan and synthesised cells |
| **Tooling** | qpuforge v0.3: `beta1_review` (2263 s, 60,000 shots per point) |

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

1. **Faster readout cancelled the materials upgrade.** With κ = 8 MHz behind the per-feedline Purcell filters, Purcell decay became the T1 limit: mean T1 was 62 µs, below Alpha-1. Individual Purcell filters (suppression 1000) restore a dielectric-limited T1 of 147 µs.
2. **DQLR did not pay.** Data-qubit leakage removal raised the 50-round logical failure from 0.147 to 0.166 at d = 5. It is now an option, not the baseline.
3. **Two DRC violations in the 49-qubit layout.** Two resonator coupling pads ran into neighbouring coupler buses. The pads now leave the serpentine upwards, and the layout is DRC clean.

Rev A still met the distance-scaling gate despite its T1 problem (d=3 0.224 %, d=5 0.067 %, Λ = 3.34).

**Device model (rev B):** mean T1 147 µs (best 158 µs), T2 218 µs, readout error 0.55 % at 300 ns.
CZ: 80 pairs, all calibrated at pulse level. Mean 43.5 ns, mean error 8.92e-04 including decoherence,
leakage 5.5e-05 per gate.

## 2. Error suppression with code distance

Memory-Z experiments, Pauli-frame Monte Carlo, exact MWPM decoding. Logical error per round, fitted over the round
counts shown:

| | d = 3 (5–20 rounds) | d = 5 (5–20 rounds) | d = 7 (9 rounds, 100,000 shots) |
|---|---|---|---|
| **Beta-1** | 0.0815 % | **0.0172 %** | 0.0057 % |
| Alpha-1 v1.1 physics (reference) | 0.3741 % | 0.1546 % | — |

- **Λ(3→5) = 4.75.** Λ(5→7) = 3.10, comparing d = 5 and d = 7 at the same 9 rounds. At d = 7 the decoder is exact up to 22 defects and greedy beyond, so this value is conservative.
- Distance-5 round 734 ns → logical lifetime **4.28 ms**, i.e. 27.1× the best physical T1.
- Against Alpha-1 v1.1 at distance 3, the Beta-1 distance-5 logical qubit is **21.8× better**.

## 3. Long experiments and leakage (d = 5, 50 rounds)

| | P_L after 50 rounds | leaked data / measure | round |
|---|---|---|---|
| **Leakage-removing reset (Beta-1 baseline)** | 0.0467 | 7.9e-05 / 1.1e-04 | 728 ns |
| + DQLR every round (option) | 0.0548 | 8.1e-06 / 1.0e-04 | 758 ns |

**Decision confirmed:** the leakage-removing reset alone gives the lower 50-round failure. DQLR clears data-qubit leakage, but its time and gate cost outweigh the benefit at this leakage rate, so it stays an option for hardware that turns out to leak more.

## 4. Error budget and remaining sensitivity (Beta-1 physics at d = 3)

Leave-one-out (eps change when the source is removed): T1 -73 %, control -54 %, Tphi -22 %, leakage -15 %, readout -13 %, zz -10 %, reset -7 %, thermal +5 %.

Next-step sensitivity: dielectric Q 6e6 → 1.2e7 -46 %, CZ residual 5e-4 → 2.5e-4 -24 %, readout 300 → 200 ns -9 %, readout separation 0.5 % → 0.25 % -5 %.

## 5. Fabrication yield (d = 5, frequency-collision rules)

| Junction resistance spread | Collision-free chips | Frequency scatter (rms) |
|---|---|---|
| 2 % (untrimmed) | 192/200 (96 %) | 57 MHz |
| 1 % | 185/200 (92 %) | 29 MHz |
| **0.5 % (trimmed)** | **199/200 (100 %)** | 14 MHz |

## 6. Wiring and layout

- **Lines:** 49 XY + 49 qubit flux + 80 coupler flux + 3 readout in + 3 readout out = **184 coaxial lines**, 3 feedlines (~16 qubits per feedline).
  - Readout frequencies are unique on each feedline: 6.90–7.86 GHz in 60 MHz steps, with spacing/κ = 7.5.
  - The wiring needs high-density flex cabling and a large-format dilution refrigerator. Each readout output chain is a TWPA plus a HEMT.
  - Phase Gamma will need more feedlines.
- **Layout:** `BF49.gds` (cell `BF49_BETA1`), built from the synthesised Alpha-1 v1.1 cells. The circuit targets are the same because the frequency plan is unchanged.
  - 9.4 × 9.4 mm, 42,196 etch polygons, 227 bump pads.
  - DRC: clean.
- **Interposer** (feedlines, Purcell filters, fan-out to edge launchers): designed in Beta; not drawn here.

## 7. Beta simulation gate (ROADMAP.md)

| Criterion | Result |
| --- | --- |
| Λ(3→5) ≥ 2 | ✔ pass |
| logical error per round at d = 5 ≤ 0.1 % | ✔ pass |
| d = 5 at least 3× below Alpha-1 v1.1 at d = 3 | ✔ pass |
| d = 5 logical lifetime > 5× best physical T1 | ✔ pass |
| trimmed fabrication yield ≥ 80 % collision-free | ✔ pass |
| chosen leakage scheme at least as good as the alternative over 50 rounds | ✔ pass |
| layout DRC clean | ✔ pass |

## 8. Assumptions and limitations

- The upgrades are **targets**: the Beta-1 numbers assume they are achieved (Q = 6e6, 35 mK, 300 ns readout at 0.5 % separation error, CZ residual 5e-4). Alpha hardware characterisation must confirm the baseline first. The sensitivity table shows what each target is worth.
- All Alpha-1 v1.1 limitations apply, as listed in `../alpha1-v11/ALPHA-1-v1.1.md` §8: CZ residual allowance, frame-level leakage, quasi-static capacitances, first-cut layout.
- The frequency-collision yield uses the v0.2 collision rules with static (sweet-spot) frequencies. Tunable qubits and per-pair CZ calibration give extra margin that the rules do not credit.
- d = 7 numbers come from the frame model with a decoder that falls back to greedy matching beyond 22 defects.

## 9. Files

| File | Content |
|---|---|
| `beta1.json` | Chip design (49 qubits, 80 couplers) |
| `beta1_review.json/.txt` | Full review output |
| `BF49.gds`, `BF49_layout.png` | GDSII layout and preview |
| `../../results/cz_cache.json` | Pulse-level CZ calibrations (shared with Alpha-1 v1.1) |
