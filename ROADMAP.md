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
| 0 — Tooling | atomforge, qpuforge v0.1–v0.3 | Materials and processor emulation, validated physics | done (96 validation checks) | n/a |
| **Alpha** | **AF-17 (Alpha-1 v1.1)** | First error-correction-ready chip: distance-3 surface code, 17 qubits | passed (designs/alpha1-v11) | open: fabricate and characterise |
| **Beta** | **BF-49 (Beta-1)** | Show logical error falling with code size: distance 5, 49 qubits | see designs/beta1 | open |
| Gamma | GF-97 + 2×d5 | Distance 7; logical operations (lattice surgery); real-time decoding | planned | planned |
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

## Phase Alpha — AF-17

**Objective.** Fabricate a 17-qubit distance-3 surface-code chip whose measured behaviour can be tied to the model.
After this phase the model is calibrated against our own hardware, not just against the literature.

### Alpha design iterations

- **v1.0** (qpuforge v0.2): architecture and frequency plan. Tunable couplers (abstract model), fast readout and active reset. Simulation: 0.45% logical error per round.
- **v1.1** (qpuforge v0.3, this release): the open modelling gaps of v1.0 are closed.
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
| CZ coherent error / leakage per gate (pulse simulation) | ≤ 5e-4 / ≤ 2e-4 |
| Idle ZZ at the coupler park point | ≤ 10 kHz |
| Layout | Placement/spacing checks clean (mask width/space DRC clean before fabrication); extracted E_C within 3%, couplings within 6% of target; residual parked exchange re-nullable |
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

**Main risks.**
- TLS defects lowering T1 → mitigated by materials work: tantalum or improved surface treatment.
- Flux-pulse distortion raising the CZ error → mitigated by predistortion and closed-loop calibration.
- Frequency collisions after fabrication → mitigated by junction trimming and the flux-tunable qubits.
- Wiring and packaging (~70 lines) → mitigated by flip-chip with an interposer.

## Phase Beta — BF-49

**Objective.** Demonstrate error suppression with code size: a distance-5 logical qubit beating distance 3 by Λ ≥ 2.

**Design.** The Alpha-1 v1.1 architecture and cells at distance 5, plus the upgrades the Alpha reviews ranked highest:
- dielectric Q ×2;
- effective temperature 35 mK;
- 300 ns readout;
- CZ residual control error 5e-4;
- individual Purcell filters, to keep T1 limited by materials at the faster readout;
- leakage-removing reset (data-qubit leakage removal evaluated, kept as an option);
- junction trimming to 0.5% spread.

Details are in `designs/beta1/BETA-1.md`.

### Beta simulation gate (all must hold)

| Criterion | Threshold |
|---|---|
| Λ(3→5) | ≥ 2 |
| Distance-5 logical error per round | ≤ 0.1% |
| Distance-5 vs Alpha-1 v1.1 distance-3 | ≥ 3× lower |
| Distance-5 logical lifetime | > 5× best physical T1 |
| Trimmed fabrication yield (0.5% spread) | ≥ 80% collision-free |
| Leakage removal over 50 rounds | chosen scheme at least as good as the alternative |
| Layout | Placement/spacing checks clean in simulation; mask width/space DRC clean before tape-out |

### Beta hardware gate

| Criterion | Threshold |
|---|---|
| Measured Λ(3→5) | ≥ 1.5 (stretch: ≥ 2) |
| Distance-5 logical error per round | ≤ 0.3% |
| Median T1 | ≥ 100 µs |
| CZ error (simultaneous) | ≤ 3e-3 |
| Readout error at 300 ns | ≤ 1% |
| Long runs (≥ 100 rounds) | no growth of the detection-event rate |

**Main risks.**
- Wiring scale (~180 lines): needs high-density flex cabling and a large dilution refrigerator.
- Decoder latency for real-time operation.
- Run-to-run calibration drift across 80 couplers.

## Phase Gamma — distance 7 and logical operations (planned)

- Distance 7 (97 qubits), targeting Λ(5→7) ≥ 2.
- Two distance-5 patches with lattice surgery, giving a logical CNOT.
- Real-time decoding at < 1 µs per round: FPGA first, then an ASIC.
- Leakage-aware decoding using heralded leakage from 3-state readout.
- qpuforge work this phase needs:
  - multi-patch circuits;
  - correlated-error models, including cosmic-ray bursts;
  - a decoder for more than 22 defects (sparse blossom / union-find).

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
  - Beta: ≥ 180 lines.
  - Both: TWPA plus HEMT readout chains.
- **Control electronics.** One flux line per qubit and per coupler. Microwave drive and readout multiplexing (~6–8 qubits per feedline).
- **Software.** qpuforge model updates after every characterisation run (the model-vs-measurement gate).
