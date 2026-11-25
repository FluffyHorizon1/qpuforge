# qpuforge — superconducting QPU chip design, whole-processor emulation and error correction

Header-only C++17. You describe a chip (junctions, capacitors, couplers, readout resonators, materials). qpuforge works out
the device physics of every qubit and coupler. It then runs real quantum circuits on a noisy model of the **whole processor**:
OpenQASM 2.0, built-in benchmarks, or full **surface-code error-correction cycles** with mid-circuit measurement and reset.
A change to the chip design shows up as a change in algorithm performance and logical error rate.

```
chip.json ──► device physics ─────────► noise model ──► compiled circuit ──► noisy emulation ───────────► metrics
 (Rn, C, Cc,    transmon spectra,         T1, Tphi,      layout, routing,    statevector trajectories      fidelity, XEB,
  couplers,     C-matrix couplings,       thermal, ZZ,   RZ/SX/X/CZ/         (exact noise) or              GHZ, error budgets,
  resonators,   ZZ, mediated exchange,    exchange,      measure/reset,      Pauli frames (fast, QEC)      logical error per
  materials)    CZ timing, collisions     control, RO,   per-qubit timelines                               round, fab yield
                                          reset
```

**Designs** (each folder has a spec, design JSON, review data and, from v1.1, a GDSII layout):

* [`designs/alpha1/ALPHA-1.md`](designs/alpha1/ALPHA-1.md): Alpha-1 v1.0, the first error-correction-ready design
  (17 qubits, distance-3 surface code; v0.2 models).
* [`designs/alpha1-v11/ALPHA-1-v1.1.md`](designs/alpha1-v11/ALPHA-1-v1.1.md): **Alpha-1 v1.1 (AF-17)**. Physical tunable
  couplers with pulse-level CZ calibration, leakage-aware QEC, a synthesised layout and GDSII.
* [`designs/beta1/BETA-1.md`](designs/beta1/BETA-1.md): **Beta-1 (BF-49)**, a 49-qubit distance-5 processor.
* [`ROADMAP.md`](ROADMAP.md): programme phases (Alpha → Beta → Gamma → Delta) with simulation and hardware exit criteria.

### What's new in v0.3

| Layer | Module | What it does |
|---|---|---|
| Pulse-level CZ | `pulse.hpp` | 3-body qubit–coupler–qubit charge Hamiltonian (exact transmon spectra, smooth sign gauge), piecewise-exact propagation, 3-knob calibration (interaction flux, flat-top length, coupler bias), exchange-free coupler park at the sweet spot, flux-noise and spectator-crossing errors, per-pair calibration of a whole chip with a cache |
| Leakage-aware QEC | `qec.hpp` | Per-qubit leakage flag beside the Pauli frame: CZ leakage injection, leaked-partner phase scrambling, transport, seepage (2/T1), \|2⟩ reads as 1, leakage-removing reset, data-qubit leakage removal (DQLR) |
| Physical layout | `layout.hpp`, `gds.hpp` | Method-of-moments capacitance extraction (validated to 0.04 % / 1.3 %), parametric Xmon/coupler/readout cells, Gauss-Newton synthesis of cell dimensions against circuit targets, chip assembly with readout resonators and control lines, DRC, GDSII writer, post-layout circuit check |
| Fast eigensolver | `core.hpp` | Householder + implicit QL (`sym_eigen_tri`), 5× faster than Jacobi at n = 64 |

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
./build/test_qpuforge          # 96 validation checks (~6 min on 2 cores)
./build/alpha1_v11_review      # Alpha-1 v1.1 review  -> designs/alpha1-v11/alpha1_v11_review.json
./build/beta1_review           # Beta-1 review        -> designs/beta1/beta1_review.json
```

## CLI

```bash
qpuforge new chip.json --preset alpha1                      # Alpha-1 (add --distance 5 for the 49-qubit extension)
qpuforge new chip.json --preset alpha1v11                   # Alpha-1 v1.1 (physical couplers, leakage model)
qpuforge new chip.json --preset beta1                       # Beta-1 (distance 5, 49 qubits)
qpuforge cz chip.json --cache results/cz_cache.json --out cal.json      # pulse-level CZ calibration of every pair
qpuforge qec chip.json --leakage-removal lrr|lrr+dqlr|none --cache results/cz_cache.json
qpuforge gds chip.json chip.gds --cells cells.json [--synthesize]       # layout synthesis + DRC + GDSII
qpuforge new chip.json --surface 3 --arch tunable           # any surface-code chip
qpuforge new chip.json --grid 3x3 --arch tunable --g-MHz 4  # grid / line / ring chips
qpuforge analyze chip.json                                  # device report + design-rule check
qpuforge layout chip.json chip.svg                          # floor-plan
qpuforge run chip.json examples/ghz5.qasm --traj 400        # emulate a circuit (mid-circuit measure/reset supported)
qpuforge run chip.json xeb:10                               # also ghz:N qft:N bv:N
qpuforge budget chip.json xeb:10                            # which noise source dominates a circuit
qpuforge qec chip.json --statevector 300                    # surface-code memory: logical error per round
qpuforge fab chip.json --samples 300                        # fabrication-yield Monte Carlo
qpuforge run chip.json qft:6 --off zz,control               # ablation: switch noise sources off
```

`new` options: `--grid RxC | --line N | --ring N | --surface D`, `--f-hi/--f-lo` (GHz), `--alpha`, `--g-MHz`, `--arch fixed|tunable`,
`--Q`, `--plan 2|4|8`, `--stagger-MHz`, `--feedlines`, `--kappa-MHz`, `--chi-MHz`, `--purcell`, `--readout-ns`, `--reset-ns`,
`--reset-err`, `--active-reset`, `--err2q`, `--residual-zz-kHz`. Designs are plain JSON, so you can edit any qubit or coupler by hand.

## Physics and emulation model

| Layer | What is computed | How |
|---|---|---|
| Junctions | E_J from normal-state resistance | Ambegaokar–Baratoff; SQUID flux tuning with asymmetry |
| Capacitance | E_C of every qubit, every coupling, parasitics | Full Maxwell capacitance matrix, inverted |
| Transmon | f01, anharmonicity, charge dispersion, ⟨i\|n\|j⟩ | Exact diagonalisation in the charge basis |
| Couplers | Static ZZ, dressed frequencies | Exact diagonalisation of each coupled pair (5 levels each) |
| Crowding | Exchange between near-degenerate qubits | 2nd-order Schrieffer–Wolff through every other qubit; parked tunable couplers block it |
| Coherence | T1 = dielectric (Q) + Purcell; Tφ = flux 1/f + charge + thermal photons | Koch 2007, Clerk–Utami; thermal population at T_eff |
| CZ gate | Duration, interaction point, feasibility | t = shape / (2√2 g) + ramps; tunable couplers via g_max and residual ZZ |
| Design rules | Degenerate pairs, \|11⟩–\|20⟩ resonance, ZZ > 1 MHz, infeasible CZ, parasitic resonance | Flags + collision count, used for fab yield |
| Inverse design | Targets (f01, α, g, χ) → C_shunt, R_n, C_c, resonator coupling | Iterated against the full capacitance-matrix model |
| Compiler | Layout, SWAP routing, 1q merging, {RZ, SX, X, CZ, measure, reset} | ≤ 2 physical pulses per 1q gate; per-qubit ASAP timelines cut into exact time slices |
| Processor | Every slice: native ops → control Pauli errors → ZZ conditional phases → resonant exchange → amplitude damping (T1 + thermal) → dephasing | Quantum trajectories, or exact density matrix (≤ 10 qubits) |
| Dynamic circuits | Mid-circuit measurement (collapse + assignment error), reset (+ reset error) | Classical record per shot |
| Error correction | Rotated surface code of any odd d; memory experiment; detector error model; decoder | Pauli-frame Monte Carlo with state-aware decoherence and coherent-ZZ windows; DEM from single-fault enumeration; exact minimum-weight matching |

## Validation (`tests/test_qpuforge.cpp`, 96 checks, all passing)

* 1q native decomposition: 3,016 unitaries, worst error 8e-16.
* Transmon spectra vs an independent numpy phase-grid solver (≤ 2e-4 GHz). Charge dispersion and ⟨0|n|1⟩ vs Koch et al.
  asymptotics. E_J(10 kΩ, 180 µeV) = 14.04 GHz.
* Two-transmon ZZ vs 2nd-order perturbation theory (within 5 %), and ZZ ∝ g². Inverse design hits f01 within 1 MHz,
  α within 2 MHz and g within 0.05 MHz.
* Compiler: 48 random circuits (CX, CP, CCX, SWAP, arbitrary U) on 4 chip shapes reproduce the logical state to 1e-14 after
  routing and time-slice scheduling.
* Noise: T1, Ramsey T2 and thermal equilibrium exact vs analytic. Exchange matches the Rabi formula to 1e-9. Trajectories match
  the exact density matrix with every noise source on.
* Benchmarks: GHZ metric = max over phase of ⟨GHZ_φ|ρ|GHZ_φ⟩. Sampled XEB matches the exact density-matrix XEB.
* Dynamic circuits: measure/reset/Bell semantics exact; reset error rate statistically correct.
* Surface code (d = 3, 5, 7): stabilisers commute, logicals anticommute, CZ schedule conflict-free. Noiseless memory
  experiments give zero detection events.
* Circuit distance: every single fault is corrected at d = 3 (6,434 faults). At d = 5 every single fault (12,975) and 20,000
  random fault pairs are corrected.
* Pauli-frame model vs exact statevector trajectories: detection rates and logical error agree within statistics.
* Eigensolver: Householder + QL vs Jacobi on random symmetric matrices (eigenvalues to 1e-13, orthogonality 4e-15).
* Pulse-level CZ: charge matrix elements stay continuous over a 0–0.32 Φ0 flux sweep (no sign flips). Idle evolution
  reproduces the static ZZ phase to 1e-14 rad. The |11⟩–|20⟩ gap with the coupler on matches 2nd-order perturbation theory
  (36.2 vs 34.9 MHz). The designed g12 nulls the parked exchange (< 1e-5 MHz). Calibrated CZ: 40 ns, 1−F = 1.0e-4,
  leakage 6.8e-5, |φ−π| < 2e-3. Halving dt changes F by 6e-7.
* Leakage frame model: bit-identical to the leakage-free model when all rates are zero. Gate-leakage plus seepage
  populations match the two-state Markov solution (within 3.5 %). A leakage-removing reset clears measure-qubit leakage
  (5× lower) and DQLR clears data-qubit leakage (3× lower).
* Per-pair calibration reaches the noise model and the frame program: overridden pairs carry their own control error
  and leakage.
* Layout: unit square plate capacitance vs the accepted value 0.3667874·4πε0·a (−0.04 %). CPW capacitance per length vs
  conformal mapping (−1.3 % on the design mesh, −0.6 % refined). Overlapping conductors are rejected. The AF-17 chip
  passes DRC, all 17 resonators sit at their quarter-wave lengths, and the GDSII stream parses back record by record.
  `tools/gds_preview.py` is an independent Python GDSII reader and renderer.

## Example: 3×3 chip design study (`examples/design_study.cpp`)

| Design (3×3 chip) | t_CZ | mean \|ZZ\| | GHZ-9 | random circuits (10 cycles) | QFT-6 |
|---|---|---|---|---|---|
| fixed coupling g = 2 MHz | 275 ns | 13 kHz | 0.73 | 0.71 | 0.32 |
| fixed coupling g = 6 MHz | 98 ns | 116 kHz | 0.27 | 0.58 | 0.02 |
| fixed coupling g = 12 MHz | 54 ns | 462 kHz | 0.06 | 0.10 | 0.01 |
| **tunable coupler** (g_max 40 MHz, 5 kHz residual) | 23 ns | 5 kHz | **0.93** | **0.86** | **0.77** |

Error correction (`examples/alpha1_review.cpp`, distance-3 surface code, logical error per round): Gen-0 fixed coupling
1.23 % → tunable couplers 0.54 % → plus fast readout 0.54 % → plus active reset **0.45 % (Alpha-1)**. The same physics at
distance 5 gives 0.18 % (Λ ≈ 2.6). Full data is in `results/` and `designs/alpha1/`.

The tunable coupler wins clearly. With fixed coupling, faster gates and always-on ZZ crosstalk trade against each other.

**v0.3 results** (simulation; see the design specs for assumptions):

| Design | Code distance | Logical error per round | Notes |
|---|---|---|---|
| Alpha-1 v1.1 (AF-17) | 3 | **0.370 %** | pulse-level CZ (36–49 ns), leakage model, leakage-removing reset; logical lifetime 251 µs = 3.3× best T1 |
| Alpha-1 v1.1 physics | 5 | 0.128 % | Λ(3→5) = 2.89 |
| Beta-1 (BF-49) | 3 / 5 / 7 | 0.0815 % / **0.0172 %** / 0.0057 % | Λ(3→5) = 4.75, Λ(5→7) = 3.1; d = 5 logical lifetime 4.3 ms = 27× best T1 |

## Limitations (read before trusting numbers)

* **CZ residual control error is an assumption.** The pulse simulation gives the coherent error (≈ 1e-4), leakage, and
  flux-noise terms. Pulse distortion, calibration drift, TLS defects and crosstalk are not simulated. They enter as one
  allowance: 1e-3 for Alpha-1 v1.1 and 5e-4 for Beta-1. This allowance dominates the CZ error and is the single most
  important number to measure on hardware.
* Pulse simulation is closed-system: T1 and Tφ during the gate are applied by the processor and QEC emulators, and the
  1/f flux noise at the excursion point is added analytically (quasi-static). Spectator crossings use a Landau-Zener
  estimate, not a many-qubit pulse simulation.
* Leakage is modelled at the frame level (flags, phase scrambling, transport, seepage, removal), not with 3-level
  dynamics in the QEC circuits. The decoder is leakage-unaware: no heralding from 3-state readout.
* **Capacitances come from a quasi-static MoM on an infinitely thick substrate.** Box modes, airbridges, the
  flip-chip interposer, kinetic inductance, flux and readout crosstalk, and participation ratios (TLS loss) are not
  extracted. Before tape-out, cross-check the cells with a full-wave/FEM tool (Palace, Sonnet, Q3D).
* The extraction assembles each qubit's self-capacitance from a star cell and the pair couplings from a pair cell,
  i.e. nearest-neighbour locality. Boundary qubits with fewer couplers have up to ~5 % different E_C. Their junctions are
  retargeted, but α shifts slightly.
* The GDSII layout is a first cut. It has no wiring to the chip edge: all control lines end on flip-chip bump pads,
  and feedlines live on the interposer, which is not drawn. There are no airbridges, test structures, alignment marks or
  dicing lanes.
* **No dynamical decoupling.** Static-ZZ penalties for fixed-coupling designs are therefore upper bounds.
* Not modelled: correlated errors (quasiparticle bursts, cosmic rays), measurement-induced dephasing, readout and flux
  crosstalk.
* Pauli-frame QEC model: accurate for incoherent noise. Strong coherent ZZ (tens of kHz and up) is under-counted, so use
  `--statevector` to cross-check such designs. The MWPM decoder is exact up to 22 defects and greedy beyond. Distance-7
  results with many defects are therefore pessimistic.
* The statevector backend is practical to about 26–28 qubits; the Pauli-frame backend handles d = 5–7 easily.

## Roadmap

See [`ROADMAP.md`](ROADMAP.md) for programme phases and exit criteria. Next on the tooling side:

1. Full-wave/FEM cross-check of the cells; participation ratios → TLS-limited T1 per element.
2. Multi-qubit pulse simulation of spectator crossings; flux-pulse distortion and predistortion model.
3. Leakage heralding (3-state readout) and a leakage-aware decoder; a decoder for large defect counts (union-find /
   sparse blossom) for d ≥ 7; correlated-error models.
4. Multi-patch circuits (lattice surgery) for Phase Gamma.
5. Materials loop: atomforge → DFT for junction barriers, oxides and TLS defects. Loss tangents → Q → T1.
6. Dear ImGui/ImPlot front-end; GPU kernels for 30+ qubit statevector emulation.

## Licence

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

This is proprietary software of Secured Cyber Solutions Ltd, licensed under the
[Secured Cyber Solutions Proprietary Software Licence v1.0](LICENSE). It is **not open source**:
no use, copying, modification or distribution is permitted without written permission from SCS.
Third-party components are listed in [NOTICE](NOTICE).

Licensing enquiries: https://securedcybersolutions.co.uk
