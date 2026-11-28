# Secured Cyber Solutions — superconducting QPU programme roadmap

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

This roadmap takes the programme from emulation to a fault-tolerant processor. It is organised as phases with
measurable gates. Every gate has two parts:

- a **simulation gate**: the qpuforge reviews in this repository;
- a **hardware gate**: chip characterisation after fabrication.

A phase is complete only when both gates pass. All numbers marked "sim" come from the qpuforge models, with the
assumptions listed in each design spec. They are design targets, not measured performance.

| Phase | Chip | Purpose | Simulation gate | Hardware gate |
|---|---|---|---|---|
| 0 — Tooling | atomforge, qpuforge v0.1–v0.4 | Materials and processor emulation, validated physics | done (500 validation checks on `main`, 504 with the mask-connectivity work) | n/a |
| **Alpha** | **AF-17 (Alpha-1 v1.1)** | First error-correction-ready chip: distance-3 surface code, 17 qubits | passed (designs/alpha1-v11), including mask DRC and metal connectivity; foundry DRC not yet run | open: fabricate and characterise |
| **Beta** | **BF-49 (Beta-1 rev C)** | Show logical error falling with code size: distance 5, 49 qubits | criteria met in simulation, including mask DRC and metal connectivity; control-line budget missed by the flux line (9.7× vs 10×); foundry DRC not yet run (designs/beta1) | open |
| **Gamma** | **GF-97 (Gamma-1)** + 2×d5 | Distance 7; logical operations (lattice surgery); real-time decoding | memory criteria met in simulation, including mask DRC and metal connectivity; control-line budget missed by the flux line (9.7× vs 10×); foundry DRC not yet run (designs/gamma1). Lattice-surgery ZZ, XX and CNOT simulated at d = 3, 5 (results/surgery_study.txt) | planned |
| Delta | modular | Several logical qubits linked between modules; logical-level benchmarks | planned | planned |

---

## Phase 0 — Tooling (complete)

- **atomforge**: atomistic alloy and structure search.
  - Embedded-atom potentials (EAM), structure relaxation, molecular dynamics (MD).
  - Convex-hull and Monte Carlo (MC) phase searches.
  - DFT probe via PySCF.
  - Use in the programme: materials screening for junction and base-layer metals.
- **qpuforge**: the processor emulator.
  - Device physics; compiler; state-vector and density-matrix emulators.
  - Surface-code QEC with an exact MWPM decoder.
  - Pulse-level tunable-coupler CZ (v0.3).
  - Leakage-modelled QEC (v0.3): leakage is simulated; the decoder is not leakage-aware.
  - Method-of-moments (MoM) capacitance extraction, layout synthesis and GDSII output (v0.3).
  - Exact matching decoder for any number of defects; memory experiments in both bases (v0.4).
  - Correlated error bursts from high-energy impacts (v0.4).
  - Hardware-gate evaluation from measured data, with a virtual chip to rehearse it (v0.4).

## Phase Alpha — AF-17

**Objective.** Fabricate a 17-qubit distance-3 surface-code chip whose measured behaviour can be tied to the model.
After this phase the model is calibrated against our own hardware, not just against the literature.

### Alpha design iterations

- **v1.0** (qpuforge v0.2): architecture and frequency plan. Tunable couplers (abstract model), fast readout and active reset. Simulation: 0.45% logical error per round.
- **v1.1** (qpuforge v0.3): the open modelling gaps of v1.0 are closed.
  - Physical coupler transmons with pulse-level CZ calibration.
  - Leakage model, with a leakage-removing (multi-level) reset chosen.
  - MoM-synthesised layout and GDSII.
  - Post-layout circuit check.
  - Results are in `designs/alpha1-v11/ALPHA-1-v1.1.md`.

### Alpha simulation gate (all must hold)

| Criterion | Threshold |
|---|---|
| Distance-3 logical error per round (memory-Z, 25 rounds) | ≤ 0.6% |
| Logical lifetime vs best physical T1 | > 2× |
| CZ coherent error / leakage per gate (pulse simulation), every pair | ≤ 5e-4 / ≤ 2e-4 |
| Idle ZZ at the coupler park point | ≤ 10 kHz |
| Layout | Placement/spacing checks and mask width/space DRC clean; extracted E_C within 3%, couplings within 6% of target; residual parked exchange re-nullable |
| Metal connectivity | Conductors separate; shorted lines and all ground on the plane (ground straps allowed); SQUIDs island-to-plane |
| Control lines | Each qubit's XY and flux line T1 limits ≥ 10× its design T1 (reported; to be gated once all three chips meet it) |
| Fabrication Monte Carlo at 2% junction spread | ≥ 90% collision-free |

### Alpha hardware gate (after fabrication)

| Criterion | Threshold | Why |
|---|---|---|
| Median T1 / T2 | ≥ 50 µs / ≥ 40 µs | QEC budget assumes ~65–75 µs |
| Single-qubit gate error (simultaneous RB) | ≤ 1e-3 | |
| CZ error (simultaneous, interleaved RB / XEB) | ≤ 5e-3 | the model uses ~2e-3 incl. decoherence |
| CZ leakage per gate | ≤ 1e-3 | leakage study shows tolerance to ~10× the simulated rate with LRR |
| Readout assignment error at 500 ns | ≤ 2% | |
| Idle ZZ with couplers parked | ≤ 20 kHz | |
| Distance-3 memory: detection-event rate stable over 25 rounds | drift < 10% | stability shows leakage is under control |
| Distance-3 logical error per round | ≤ 1.5% | |
| Model vs measurement | ε within 2× of the re-run model fed with measured parameters | validates the model for Beta |

**Deliverables:** GDSII (`AF17_v11.gds`), design JSON, CZ calibration table, a characterisation plan, and an updated
model fitted to measured T1/T2, gate errors and leakage.

- The characterisation plan is `designs/alpha1-v11/CHARACTERISATION-PLAN.md`. It defines the data file
  (`qpuforge.characterisation.v0`) and the order of measurements.
- `qpuforge hwgate` evaluates every row of the table above from that file, including the last one: it puts the measured
  parameters into the model, re-runs the memory experiments and compares. The chain has been rehearsed on a virtual chip;
  it has not seen hardware data.

**Main risks.**
- TLS defects lowering T1 → mitigated by materials work: tantalum or improved surface treatment.
- Flux-pulse distortion raising the CZ error → mitigated by predistortion and closed-loop calibration.
- Frequency collisions after fabrication → mitigated by junction trimming and the flux-tunable qubits.
- Wiring and packaging (~70 lines) → mitigated by flip-chip with an interposer.

## Phase Beta — BF-49

**Objective.** Demonstrate error suppression with code size: a distance-5 logical qubit beating distance 3 by Λ ≥ 2.

**Design.** The Alpha-1 v1.1 architecture at distance 5, plus the upgrades the Alpha reviews ranked highest:
- dielectric Q ×2;
- effective temperature 35 mK;
- 300 ns readout;
- CZ residual control error 5e-4;
- individual Purcell filters, to keep T1 limited by materials at the faster readout;
- leakage-removing reset (data-qubit leakage removal evaluated, kept as an option);
- junction trimming to 0.5% spread.

Rev C (qpuforge v0.4) corrects two things the second review found:
- the readout plan: at most 9 resonators per feedline (6 feedlines), where rev B had up to 17 on each of 3 (17, 17 and 15);
- the cells: re-synthesised for the Beta-1 readout coupling, which the reused Alpha-1 cells missed by 30%.

Details are in `designs/beta1/BETA-1.md`.

### Beta simulation gate (all must hold)

| Criterion | Threshold |
|---|---|
| Λ(3→5), both memory bases | ≥ 2 |
| Distance-5 logical error per round, both bases | ≤ 0.1% |
| Distance-5 vs Alpha-1 v1.1 distance-3 | ≥ 3× lower |
| Distance-5 logical lifetime | > 5× best physical T1 |
| Trimmed fabrication yield (0.5% spread) | ≥ 80% collision-free |
| Leakage removal over 50 rounds | chosen scheme at least as good as the alternative |
| Detection rate over 100 rounds with the leakage-removing reset | drift < 5% |
| Layout | Placement/spacing checks and mask width/space DRC clean; extracted E_C within 3%, couplings and readout capacitance within 6% of target; residual parked exchange re-nullable. Foundry DRC deck clean before tape-out |
| Metal connectivity | as Alpha |
| Control lines | as Alpha |

### Beta hardware gate

| Criterion | Threshold |
|---|---|
| Measured Λ(3→5) | ≥ 1.5 (stretch: ≥ 2) |
| Distance-5 logical error per round | ≤ 0.3% |
| Median T1 | ≥ 100 µs |
| CZ error (simultaneous) | ≤ 3e-3 |
| Readout error at 300 ns | ≤ 1% |
| Long runs (≥ 100 rounds) | no growth of the detection-event rate (drift < 5%) |
| Model vs measurement | ε within 2× of the re-run model fed with measured parameters |

`qpuforge hwgate --phase beta` evaluates this table. The measured Λ(3→5) needs distance-3 memory runs on sub-patches of
the chip in the same data file.

**Main risks.**
- Wiring scale (~190 lines): needs high-density flex cabling and a large dilution refrigerator.
- Decoder latency for real-time operation.
- Run-to-run calibration drift across 80 couplers.

## Phase Gamma — distance 7 and logical operations

**Objective.** A distance-7 memory with Λ(5→7) ≥ 2, then logical operations between two distance-5 patches.

**Gamma-1 (GF-97).** Beta-1 physics and cells at distance 7: 97 qubits, 168 couplers, 11 feedlines. Details are in
`designs/gamma1/GAMMA-1.md`.

### Gamma memory gate (simulation; all must hold)

| Criterion | Threshold |
|---|---|
| Λ(5→7), both memory bases | ≥ 2 |
| Distance-7 logical error per round, both bases | ≤ 0.01% |
| Distance-7 logical lifetime | > 50× best physical T1 |
| Logical-error floor from high-energy impacts | < 10% of the distance-7 error |
| Detection rate over 40 rounds with the leakage-removing reset | drift < 5% |
| Trimmed fabrication yield (0.5% spread) | ≥ 80% collision-free |
| Layout, metal connectivity, control lines | as Beta |

### Still open in Phase Gamma

- Two distance-5 patches with lattice surgery, giving a logical CNOT. **Simulated** (`surgery.hpp`,
  `docs/design/lattice-surgery.md`): logical ZZ and XX measurement between two patches and the three-patch CNOT, at
  d = 3 and 5. With Beta-1 physics at d = 5 the ZZ outcome is wrong in 0.14 % of measurements and a CNOT on |00⟩ gets
  some output wrong in 0.71 % of runs
  (`results/surgery_study.txt`). Still open: a chip layout for two or three d = 5 patches and their routing region,
  real-time feed-forward of the corrections, and a joint decoder for all observables.
- Real-time decoding at < 1 µs per round: FPGA first, then an ASIC. The exact software decoder is the reference.
  Done in software: a weighted union-find decoder (`DecodeMethod::UnionFind`), within 1.1× of the exact decoder's
  logical error at d = 3 and 5 (`results/decoder_study.txt` for Gamma-1), and its streaming (sliding-window)
  version (`StreamingUnionFind`): at d = 7 it is within 1.16x of batch union-find once the window extends about d
  rounds past the commit region. In software it takes 1.5-1.8 us per committed round against a 0.73 us round. Open:
  the hardware (and a faster, pipelined window scheme).
- Wiring: about 380 coaxial lines at distance 7. High-density flex lines or multiplexed flux control are needed before
  this chip can be operated.

### Done in Phase Gamma (software)

- Leakage-aware decoding from heralded leakage (3-state readout of the measure qubits): `HeraldSpec`, with the
  leakage-removing reset. See `docs/design/heralded-leakage.md` and `results/decoder_study.txt`.
- An impact veto: per-shot impact injection, a sliding space-time window veto calibrated to a false-veto rate, and the
  floor with shot-by-shot discard or a latched hold. See `docs/design/impact-veto.md`.
- Linear-time fault enumeration: setup cost linear in the number of rounds (Gamma-1, d = 7, 50 rounds: 1.4 s against
  44 s for the 0.4 routine on four threads).

**Main risks.**
- Wiring scale (above).
- High-energy impacts set a floor that distance does not remove. It is far below the distance-7 error, and becomes
  the limit from about distance 15 at the simulated Λ. Gap-engineered junctions shorten each event; they should be in
  the process by Delta, together with an impact veto.
- Calibration time for 168 couplers.

## Phase Delta — modular scale-out (planned)

- Several logical qubits.
- Inter-module links: microwave links or transduction.
- Logical-level benchmarks.
- Service integration with the SCS product line.

---

## Programme dependencies

- **Fabrication.** Superconducting foundry or cleanroom access.
  - Alpha: Al/AlOx/Al junctions on Nb or Ta base layer.
  - Beta: adds trimming.
- **Cryogenics.** Dilution refrigerator.
  - Alpha: ≥ 70 lines.
  - Beta: ≥ 190 lines.
  - Gamma: ≥ 384 lines, or multiplexed control.
  - Both: TWPA plus HEMT readout chains.
- **Control electronics.** One flux line per qubit and per coupler. Microwave drive and readout multiplexing (up to 9 qubits per feedline).
- **Software.** qpuforge model updates after every characterisation run (the model-vs-measurement gate).
