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
  couplers with pulse-level CZ calibration, leakage-modelled QEC, a synthesised layout and GDSII.
* [`designs/alpha1-v11/CHARACTERISATION-PLAN.md`](designs/alpha1-v11/CHARACTERISATION-PLAN.md): what to measure on the
  fabricated chip, the data-file format, and how the hardware gate is evaluated from it.
* [`designs/beta1/BETA-1.md`](designs/beta1/BETA-1.md): **Beta-1 (BF-49), rev C**, a 49-qubit distance-5 processor.
* [`designs/gamma1/GAMMA-1.md`](designs/gamma1/GAMMA-1.md): **Gamma-1 (GF-97)**, a 97-qubit distance-7 memory.
* [`ROADMAP.md`](ROADMAP.md): programme phases (Alpha → Beta → Gamma → Delta) with simulation and hardware exit criteria.

### Unreleased: flux lines

| Layer | Module | What it does |
|---|---|---|
| Flux-line model | `fluxline.hpp` | What a flux line does to the CZ pulse, and what the CZ needs from it. A line is a cascade of first-order sections (step response 1 + a exp(−t/τ)), a low-pass and a delay, evaluated in closed form; the AWG is a sample-and-hold with a DAC; the predistortion filter inverts a fit of the line. `simulate_cz_drive` runs the CZ of `pulse.hpp` on the flux that arrives, for a gate on its own or inside the pulse train of the memory cycle. `cz_fidelity_in_frame` scores a gate whose virtual Z gates were calibrated elsewhere, and `cz_offset_response` gives its error per unit of flux offset. The calibration does not use the model, so every review is unchanged |

Study: [`results/flux_distortion_study.txt`](results/flux_distortion_study.txt) (`examples/flux_distortion_study.cpp`); design
note: [`docs/design/flux-line-predistortion.md`](docs/design/flux-line-predistortion.md). A residual tail of 0.1 % with a time
constant of 30–100 ns costs each CZ of the Alpha-1 v1.1 cycle 1.3e-3 when the gates are used as calibrated in isolation, and
1.2e-4 when their virtual Z gates are calibrated in the cycle (the residual control-error allowance is 1e-3).

### Unreleased: Phase Gamma decoding

| Layer | Module | What it does |
|---|---|---|
| Union-find decoder | `unionfind.hpp`, `qec.hpp` | Weighted union-find (Delfosse-Nickerson growth, peeling), `DecodeMethod::UnionFind`, `qpuforge qec --decoder uf`: the real-time candidate; exact matching stays the reference. On Gamma-1 (14 rounds, same syndromes) its logical error per round is 0.98-1.10x the exact decoder's at d = 3, 5 and 7 |
| Streaming union-find | `unionfind.hpp`, `qec.hpp` | Sliding-window decoding over rounds (`StreamingUnionFind`, `qpuforge qec --decoder uf-stream --window W --commit C`): each window is decoded with the later rounds as a boundary, and its first C rounds are committed. At d = 7 (50 rounds) it gives 1.16x batch union-find's failures for 10/5 windows and 1.05x for 20/10. The window must extend about d rounds past the commit region |
| Heralded leakage | `qec.hpp` | 3-state readout of the measure qubits (`HeraldSpec`: herald fidelity, false-herald rate); the decoder lowers the weights of the edges of the heralded qubit's round, with exact matching on the per-shot re-weighted graph (`decode_scaled`). `qpuforge qec --herald`. With the leakage-removing reset on Gamma-1 d = 7: -10% logical error at the chip's own leakage, -46% at 10x |
| Impact veto | `qec.hpp` | Per-shot impact injection at a random time and place (`ImpactSpec`, deterministic per (seed, shot)), a sliding space-time window veto calibrated to a false-veto rate (`VetoSpec`, `calibrate_veto`), and the floor over the burst profile with the veto (`burst_floor_veto`). `qpuforge qec --bursts --veto` |
| Linear-time setup | `qec.hpp` | Single-fault enumeration is linear in the number of rounds (per-qubit tails from each measurement cut): Gamma-1 d = 7 at 50 rounds in 1.4 s against 44 s for the 0.4 routine, identical output |

Study: [`results/decoder_study.txt`](results/decoder_study.txt) (`examples/decoder_study.cpp`). Design notes:
[`docs/design/`](docs/design/).

### Unreleased: lattice surgery

| Layer | Module | What it does |
|---|---|---|
| Lattice surgery | `surgery.hpp` | Several surface-code patches on one lattice with merge and split schedules: logical ZZ and XX measurement between two patches and the three-patch logical CNOT, with detectors and observables across the merge and split boundaries, Pauli-frame corrections, a stabiliser-tableau check of the noiseless circuit, one exact matching decoder per observable and the memory experiments' noise model. With Beta-1 physics the ZZ outcome is wrong in 0.74 % of measurements at d = 3 and 0.14 % at d = 5; a CNOT on \|00⟩ gets some output wrong in 2.7 % and 0.71 % of runs (`results/surgery_study.txt`, [design note](docs/design/lattice-surgery.md)) |
| Joint and correlated decoding | `surgery.hpp` | One decoder pass for all observables of an operation: exact matching or union-find on one graph whose edges carry an observable mask (`SurgeryDecoder::Joint`, `JointUF`; `qpuforge surgery --decoder joint`). On 1,600,000 shots it gives the same prediction as the separate decoders in every shot, in 1.9–3.4× less decode time. A two-pass correlated decoder (`Correlated`) uses the detectors of both Pauli types: at d = 5 with Beta-1 physics it lowers the CNOT's "some output wrong" from 0.71 % to 0.42 %, at 12 times the decode time of joint matching (184 µs per shot against 15 µs), so it is an offline reference, not a real-time candidate |
| Routing strips | `surgery.hpp` | ZZ, XX and the CNOT between patches that are not adjacent, through a straight routing strip of chosen length that is merged with both patches and read out at the split (`gap`; `qpuforge surgery --gap 3`). The strip adds no rounds and 2d qubits per row. The error of the merge outcome grows about linearly with the strip: at d = 5 with Beta-1 physics the ZZ outcome is wrong in 0.14 % of measurements for adjacent patches and 0.45 % through a strip of 7, and a CNOT on \|00⟩ gets some output wrong in 0.71 % and 1.82 % of runs |
| Chip account | `surgery_study` | Two d = 5 patches need 109 qubits, 188 couplers, 13 feedlines and 432 lines when adjacent and 169 qubits with a strip of 7; the three patches of a CNOT need 169 to 289 qubits. The 8-frequency plan and the 16 CZ calibrations of Beta-1 cover every lattice without a frequency collision. No layout or mask |

### What's new in v0.4

| Layer | Module | What it does |
|---|---|---|
| Exact decoding at any size | `matching.hpp`, `qec.hpp` | Edmonds' blossom algorithm (primal-dual, O(n³)) for minimum-weight perfect matching of any number of defects. The dynamic programme is kept for up to 8 defects, where it is faster. The greedy fallback above 22 defects is gone. On a 30-round distance-5 run (20 defects per shot on average) decoding takes about 39 µs per shot, against 35 ms for the v0.3 decoder (wall-clock, from `results/test_output.txt`) |
| Both memory bases | `qec.hpp` | Memory experiments in the X basis (\|+⟩ logical state, X-type detectors, logical X) beside the Z basis |
| Correlated bursts | `qec.hpp` | High-energy impacts as a chip-wide excess decay or dephasing rate that relaxes exponentially; the decoder keeps its quiet-chip priors. Gives the probability that one impact flips the logical qubit, the unprotected time, how long the burst is visible in the syndrome, and the logical-error floor |
| Hardware gate | `characterise.hpp` | A data-file format for measured chip parameters, a model built from them, and automatic evaluation of the roadmap's hardware gates, including model versus measurement. A virtual chip generates test data until hardware exists |
| Decoder and syndrome statistics | `qec.hpp` | Detection fraction per round, defects per shot, decode time per shot; `qpuforge qec --decoder legacy` replays the old decoder on the same syndromes |
| Error bars | `qec.hpp` | `fit_eps_per_round_err`: the 1σ error of a fitted logical error per round, propagated through the fit. v0.3 took it from the last round count alone, which was 20–40 % too small |
| Faster setup | `qec.hpp` | Single-fault enumeration follows only the qubits a fault touches, skips the time slices that do not act on them, and runs in parallel. On the Gamma-1 chip (distance 7, 14 rounds) it takes 3.9 s on two cores against 33 s for the v0.3 routine, with identical output |
| Designs | `presets.hpp` | Readout plan that scales with the chip (at most 9 resonators per feedline); Beta-1 rev C; Gamma-1 (GF-97) |
| Reproducible Monte Carlo | `qec.hpp`, `emulator.hpp` | Every frame-simulation shot is seeded from (seed, shot) and partial sums are added in a fixed order, so every Monte Carlo result is bit-identical for any `OMP_NUM_THREADS`. Committed outputs, except the historical Alpha-1 v1.0 review (`results/alpha1_review.*`), no longer depend on the machine's core count |
| Hardened serialisation | `json.hpp`, `gds.hpp` | RFC 8259 JSON parsing with a nesting limit; GDSII writes fail loudly on I/O errors; JSON and GDSII round-trip tests |

### What's new in v0.3

| Layer | Module | What it does |
|---|---|---|
| Pulse-level CZ | `pulse.hpp` | 3-body qubit–coupler–qubit charge Hamiltonian (exact transmon spectra, smooth sign gauge), piecewise-exact propagation, 3-knob calibration (interaction flux, flat-top length, coupler bias), exchange-free coupler park at the sweet spot, flux-noise and spectator-crossing errors, per-pair calibration of a whole chip with a cache |
| Leakage-modelled QEC | `qec.hpp` | Per-qubit leakage flag beside the Pauli frame (the decoder does not use leakage information): CZ leakage injection, leaked-partner phase scrambling, transport, seepage (2/T1), \|2⟩ reads as 1, leakage-removing reset, data-qubit leakage removal (DQLR) |
| Physical layout | `layout.hpp`, `gds.hpp` | Method-of-moments capacitance extraction (validated to 0.04 % / 1.3 %), parametric Xmon/coupler/readout cells, Gauss-Newton synthesis of cell dimensions against circuit targets, chip assembly with readout resonators and control lines, placement/spacing checks and merged-layer width/space DRC, GDSII writer, post-layout circuit check |
| Fast eigensolver | `core.hpp` | Householder + implicit QL (`sym_eigen_tri`), 5× faster than Jacobi at n = 64 |

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
./build/test_qpuforge          # 541 validation checks (~15 min on 4 cores); or: ctest --test-dir build -LE smoke
./build/test_qpuforge matching bursts   # selected groups by name (--list prints the names)
./build/alpha1_v11_review      # Alpha-1 v1.1 review  -> designs/alpha1-v11/alpha1_v11_review.json
./build/beta1_review           # Beta-1 review        -> designs/beta1/beta1_review.json
./build/gamma1_review          # Gamma-1 review       -> designs/gamma1/gamma1_review.json
./build/surgery_study          # lattice surgery      -> results/surgery_study.{txt,json}
./build/flux_distortion_study  # flux-line requirements -> results/flux_distortion_study.{txt,json}
python3 tools/fill_spec.py beta    # BETA-1.md from its template and review data (also: gamma)
```

`ctest` runs one entry per test, labelled `unit`, `physics`, `pulse`, `qec`, `layout` and `smoke`.

The `json` and `gds` unit groups (64 checks) cover the JSON parser/writer (including a lossless design-file round trip for the `gen0`, `alpha1`, `alpha1v11` and `beta1` presets) and the GDSII stream writer (record framing, UNITS encoding, coordinate/record-length guards, reproducible `SOURCE_DATE_EPOCH` dates).

- The five `smoke` entries run the review programs and the surgery study end to end with 200 shots, and the flux-line study on a reduced set. They write into the build tree and
  leave the committed outputs untouched.
- Each review program also takes `--out DIR`, `--cache FILE` and `--smoke`; see `examples/review_io.hpp`.
- Committed review outputs, except the historical Alpha-1 v1.0 review, record their source commit, compiler, seeds and
  date (`results/PROVENANCE.md`).

## CLI

```bash
qpuforge new chip.json --preset alpha1                      # Alpha-1 (add --distance 5 for the 49-qubit extension)
qpuforge new chip.json --preset alpha1v11                   # Alpha-1 v1.1 (physical couplers, leakage model)
qpuforge new chip.json --preset beta1                       # Beta-1 (distance 5, 49 qubits)
qpuforge new chip.json --preset gamma1                      # Gamma-1 (distance 7, 97 qubits)
qpuforge qec chip.json --basis both --rounds 14 --cache results/cz_cache.json   # memory in the Z and the X basis
qpuforge surgery --op cnot --d 3,5 --input 00 --cache results/cz_cache.json     # lattice surgery (also --op zz, --op xx)
qpuforge surgery --op zz --d 3 --gap 3 --decoder all --cache results/cz_cache.json   # through a routing strip, every decoder
qpuforge qec chip.json --bursts --cache results/cz_cache.json                   # + logical-error floor from impacts
qpuforge synth-chip chip.json data.json --seed 3 --cache results/cz_cache.json  # virtual chip -> characterisation file
qpuforge hwgate chip.json data.json --phase alpha --out report.json             # hardware gate from measured data
qpuforge cz chip.json --cache results/cz_cache.json --out cal.json      # pulse-level CZ calibration of every pair
qpuforge qec chip.json --leakage-removal lrr+dqlr --cache results/cz_cache.json  # leakage removal (also lrr, dqlr, none)
qpuforge gds chip.json chip.gds --cells cells.json --synthesize         # layout synthesis + DRC + GDSII (without --synthesize, an existing cells.json is reused)
qpuforge new chip.json --surface 3 --arch tunable           # any surface-code chip
qpuforge new chip.json --grid 3x3 --arch tunable --g-MHz 4  # grid / line / ring chips
qpuforge analyze chip.json                                  # device report + design-rule check
qpuforge analyze chip.json --materials sidecar.json         # apply an atomforge forge.materials.v0 sidecar first (add --loss-tangent to use Q = 1/tanδ)
qpuforge qec chip.json --materials sidecar.json             # --materials works on run, budget, fab, qec, cz, synth-chip, hwgate too
qpuforge materials apply chip.json sidecar.json --out chip2.json   # save the applied values, with a materials_source block
qpuforge analyze designs/alpha1-v11/alpha1_v11.json --materials designs/alpha1-v11/materials-alpha1-v11-handwritten.json   # hand-written sidecar (assumed values)
qpuforge layout chip.json chip.svg                          # floor-plan
qpuforge run chip.json examples/ghz5.qasm --traj 400        # emulate a circuit (mid-circuit measure/reset supported)
qpuforge run chip.json xeb:10                               # also ghz:N qft:N bv:N
qpuforge budget chip.json xeb:10                            # which noise source dominates a circuit
qpuforge qec chip.json --statevector 300                    # surface-code memory: logical error per round
qpuforge fab chip.json --samples 300                        # fabrication-yield Monte Carlo
qpuforge run chip.json qft:6 --off zz,control               # ablation: switch noise sources off
```

`new` options: `--grid RxC | --line N | --ring N | --surface D`, `--f-hi/--f-lo` (GHz), `--alpha`, `--g-MHz`, `--arch fixed|tunable`,
`--Q`, `--plan 2|4|8`, `--stagger-MHz`, `--name`, `--feedlines`, `--kappa-MHz`, `--chi-MHz`, `--purcell`, `--readout-ns`,
`--reset-ns`, `--reset-err`, `--active-reset`, `--err2q`, `--residual-zz-kHz`; a `--preset` design takes only `--distance`.
Designs are plain JSON, so you can edit any qubit or coupler by hand.

Unknown options are errors: each command accepts only the options its usage lists (`qpuforge` without arguments), and an
unknown, misspelt or repeated option, a missing value or a malformed number stops it before it reads any file. Exit
status: 0 done, 1 usage error, 2 error while running (also `gds`: layout not DRC- and connectivity-clean), 3 `hwgate`: gate not passed,
`surgery`: a noiseless check failed.
`qpuforge --help` (or `-h`, `help`) prints the usage and `qpuforge <command> --help` prints one command's usage, on stdout, exit 0.

## Physics and emulation model

| Layer | What is computed | How |
|---|---|---|
| Junctions | E_J from normal-state resistance | Ambegaokar–Baratoff; SQUID flux tuning with asymmetry |
| Capacitance | E_C of every qubit, every coupling, parasitics | Full Maxwell capacitance matrix, inverted |
| Transmon | f01, anharmonicity, charge dispersion, ⟨i\|n\|j⟩ | Exact diagonalisation in the charge basis |
| Couplers | Static ZZ, dressed frequencies; Alpha-1 idle ZZ | Pair diagonalisation (5 levels). With `model_coupler`, a third transmon is parked at the exchange null (3 levels) |
| Crowding | Exchange between near-degenerate qubits | 2nd-order Schrieffer–Wolff through every other qubit; parked tunable couplers block it |
| Coherence | T1 = dielectric (Q) + Purcell; Tφ = flux 1/f + charge + thermal photons | Koch 2007, Clerk–Utami; thermal population at T_eff |
| CZ gate | Duration, interaction point, feasibility, flux-pulse error | t = shape / (2√2 g) + ramps; off-sweet-spot 1/f noise uses the short-time t² infidelity |
| Design rules | Degenerate pairs, \|11⟩–\|20⟩ resonance, ZZ > 1 MHz, infeasible CZ, parasitic resonance | Flags + collision count, used for fab yield |
| Inverse design | Targets (f01, α, g, χ) → C_shunt, R_n, C_c, resonator coupling | Iterated against the full capacitance-matrix model |
| Compiler | Layout, SWAP routing, 1q merging, {RZ, SX, X, CZ, measure, reset} | ≤ 2 physical pulses per 1q gate; per-qubit ASAP timelines cut into exact time slices |
| Processor | Every slice: native ops → control Pauli errors → ZZ conditional phases → resonant exchange → amplitude damping (T1 + thermal) → dephasing | Quantum trajectories, or exact density matrix (≤ 10 qubits) |
| Dynamic circuits | Mid-circuit measurement (collapse + assignment error), reset (+ reset error) | Classical record per shot |
| Error correction | Rotated surface code of any odd d; memory experiment in either basis; detector error model; decoder | Pauli-frame Monte Carlo with state-aware decoherence and coherent-ZZ windows; DEM from single-fault enumeration; exact minimum-weight matching at any number of defects |
| Correlated errors | Logical-error floor from high-energy impacts | Chip-wide excess decay / dephasing with exponential recovery, sampled at a ladder of severities and integrated over the burst |
| Hardware gate | Roadmap criteria from a measured-data file; model versus measurement | Measured parameters replace predicted ones; gate errors are split into decoherence, leakage and control parts; the memory experiment is re-simulated |

## Validation (`tests/test_qpuforge.cpp`, 541 checks, all passing)

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
* JSON: parser/writer round trips, typed accessors, error paths, and a lossless `design_to_json` / `design_from_json` cycle for the `gen0`, `alpha1`, `alpha1v11` and `beta1` presets.
* GDSII: record framing, UNITS encoding, closed-ring coordinates, record-length and int32 guards, `SOURCE_DATE_EPOCH` stamps.
* Pauli-frame model vs exact statevector trajectories (1,000): detection rates and logical error agree within statistics
  (P_L 0.0240 ± 0.0048 vs 0.0237).
* Eigensolver: Householder + QL vs Jacobi on random symmetric matrices (eigenvalues to 1e-13, orthogonality 4e-15).
* Pulse-level CZ: charge matrix elements stay continuous over a 0–0.32 Φ0 flux sweep (no sign flips). Idle evolution
  reproduces the static ZZ phase to better than 1e-12 rad. The |11⟩–|20⟩ gap with the coupler on matches 2nd-order perturbation theory
  (36.2 vs 34.9 MHz). The designed g12 nulls the parked exchange (< 1e-5 MHz). Calibrated CZ: 40 ns, 1−F = 1.0e-4,
  leakage 6.6e-5, |φ−π| < 2e-3. Halving dt changes F by 6e-7. With 4/4/4 levels, one more coupler level changes the CZ
  phase by 6e-5 rad and one more qubit level by 2e-4 rad. The earlier default of 3 coupler levels was off by 4e-3 rad.
* Flux lines: the partial fractions of a cascade with a low-pass match a step-by-step integration of the sections (8e-7), and
  a section followed by its inverse is the identity. A pulse keeps its area through a line, and two pulses give the sum of
  each. Predistortion at 1 GS/s is exact for one section. The CZ on a drive reproduces `simulate_cz` (3e-16). A slow 1 %
  tail equals a pulse scaled by 1.01, and the amplitude calibration removes it. After an earlier pulse, the flux offset
  matches a quadrature of the tail and the measure qubit's phase shift matches the bare transmon's frequency shift (−0.057
  against −0.059 rad). The response to constant offsets matches the calibration's curvature recipe and a direct simulation.
* Leakage frame model: bit-identical to the leakage-free model when all rates are zero. Gate-leakage plus seepage
  populations match the two-state Markov solution (within 4 %). A leakage-removing reset clears measure-qubit leakage
  (5× lower) and DQLR clears data-qubit leakage (3× lower).
* Per-pair calibration reaches the noise model and the frame program: overridden pairs carry their own control error
  and leakage.
* Layout: unit square plate capacitance vs the accepted value 0.3667874·4πε0·a (−0.04 %). CPW capacitance per length vs
  conformal mapping (−1.3 % on the design mesh, −0.6 % refined). Overlapping conductors are rejected. The AF-17 chip
  passes the placement and island-spacing checks, all 17 resonators sit at their quarter-wave lengths, and the GDSII stream parses back record by record.
  A synthetic layer exercises the merged-etch width/space DRC; the ground-sliver clean-up and its connectivity guard, the
  metal-connectivity map, the ground-strap placer, the bus knob and the shorted staircase end have their own checks
  (`layout_geometry`). AF-17, BF-49 and GF-97, assembled from the committed cells, pass the placement checks and the
  mask DRC (0 + 0), are connectivity-clean with their ground straps, meet the extraction targets and re-null the parked
  exchange (`layout_presets`). No foundry DRC deck has been run: see CHANGELOG, Known issues.
  `tools/gds_preview.py` is an independent Python GDSII reader and renderer.
* Matching: the blossom matcher agrees with brute force on 600 random graphs (complete and sparse, even and odd). In
  the decoder it agrees with the exact dynamic programme on 2,075 defect sets (weights to 3e-14, same logical
  prediction every time) and is never heavier than greedy on large sets. On a 30-round distance-5 run it fails
  3× less often than the old greedy fallback on the same syndromes.
* Fault enumeration: the sparse and the linear-time versions both reproduce the reference fault by fault on 188,932
  faults (detectors, logical flip, probability, order; d = 3 and 5, both bases, 1 to 8 rounds, DQLR, equal and
  unequal CZ durations). The linear-time version's work per fault does not grow with the number of rounds (d = 5:
  0.85 time slices per fault at both 10 and 40 rounds, against 33.7 and 141.9 for the sparse version).
* Union-find decoder: every correction clears its syndrome (12,000 syndromes), every single fault the priors allow is
  decoded correctly (55,917 mechanisms), and at d = 3 and 5 its logical error is 1.01-1.08x the exact decoder's on
  the same syndromes (limit 1.5x).
* Heralded leakage: matching on a re-weighted graph equals brute force on 300 random cases (to 9e-16); at d = 5 with
  stressed measure-qubit leakage, heralds cut failures from 928 to 694 (exact) and from 984 to 763 (union-find) on
  the same syndromes.
* Impact veto: impacts without excess noise leave every shot unchanged; a calibrated threshold holds its false-veto
  rate on fresh data; the veto catches 69% of chip-wide impacts at d = 5 and removes most of their excess error.
* Memory-X: noiseless runs are deterministic, the detector error model is graph-like, and every single fault is
  corrected at d = 3 and d = 5. On the exact state-vector emulator a logical-Z chain flips the X memory and fires no
  detector, and leaves the Z memory alone; the same holds with X and Z exchanged.
* Fit error bar: for a four-point fit it equals the scatter of 4,000 synthetic fits (5.51e-5 against 5.52e-5) and the
  linear propagation of the point errors. The last point alone would give 3.7e-5.
* Bursts: the excess rates are exact, the logical error rises with severity, a 25 ms burst randomises the logical
  qubit, and the expected number of flips scales linearly with the recovery time.
* Hardware gate: a scatter-free "measurement" of a design returns that design's noise model (to 1e-9 or better in
  every parameter). A virtual chip passes model versus measurement (ratio 1.06); the same chip with a hidden dephasing
  source fails it (ratio 6.1). Threshold failures and missing data are reported.
* Designs: feedline plans of Beta-1 and Gamma-1, no frequency collisions, and the CZ interaction point solved at the
  tuned anharmonicity (1–2.5 MHz above the idle-anharmonicity estimate).
* Lattice surgery: on a stabiliser tableau, no detector fires and every logical outcome is right for ZZ and XX
  measurement on |00⟩ |01⟩ |10⟩ |11⟩ |++⟩ |+−⟩ |−+⟩ |−−⟩ (d = 3, 5) and for the CNOT on Z-, X-basis and Bell-pair inputs.
  A reference-qubit process check passes for Z- and X-basis patches. A single patch through the surgery engine
  reproduces the memory experiment shot for shot, with and without leakage, and its chip equals Beta-1 at d = 3. The
  d = 3 detector error models have no undetectable logical faults, and every single fault is corrected except a few
  coherent-ZZ crosstalk terms that straddle two rounds (total probability at most 2.1e-7).
  The joint decoder's graph agrees with the per-observable graphs pair by pair, and on the same shots it reproduces
  the separate decoders. Through routing strips (lengths 3 and 5) the tableau checks pass for ZZ, XX and the CNOT,
  including the reference-qubit process check; removing the strip's closing detectors or a strip qubit from a logical
  is caught.
* Materials sidecar (`forge.materials.v0`, `materials.hpp`, [`docs/DATA_CONTRACT.md`](docs/DATA_CONTRACT.md)): the
  placeholder file exactly as atomforge writes it leaves Alpha-1 v1.1 unchanged. Sidecar values of Q_dielectric and the
  junction gap reach `analyze()`: T1_diel matches 2πf/Q to 1e-12 and E_J scales with the gap. The loss tangent becomes
  Q = 1/tanδ only on request. Wrong schemas, units, values, unknown quantities and unknown or missing provenance are
  rejected. `materials apply --out` writes a chip equal to the in-memory application, every `Rn_kOhm` kept, with a
  `materials_source` block (SHA-256 of the sidecar, checked against the FIPS 180-4 vectors).
* The new groups also pass under AddressSanitizer and UndefinedBehaviorSanitizer.

## Example: 3×3 chip design study (`examples/design_study.cpp`)

| Design (3×3 chip) | t_CZ | mean \|ZZ\| | GHZ-9 | random circuits (10 cycles) | QFT-6 |
|---|---|---|---|---|---|
| fixed coupling g = 2 MHz | 275 ns | 13 kHz | 0.74 | 0.71 | 0.32 |
| fixed coupling g = 6 MHz | 98 ns | 116 kHz | 0.28 | 0.58 | 0.02 |
| fixed coupling g = 12 MHz | 54 ns | 462 kHz | 0.05 | 0.10 | 0.01 |
| **tunable coupler** (g_max 40 MHz, 5 kHz residual) | 23 ns | 5 kHz | **0.93** | **0.86** | **0.76** |

Error correction (`examples/alpha1_review.cpp`, distance-3 surface code, logical error per round): Gen-0 fixed coupling
1.23 % → tunable couplers 0.54 % → plus fast readout 0.54 % → plus active reset **0.45 % (Alpha-1)**. The same physics at
distance 5 gives 0.18 % (Λ ≈ 2.6). This is the historical Alpha-1 v1.0 review (pre-PR #2 model, not regenerated; see
`designs/alpha1/README.md`); full data is in `results/` and `designs/alpha1/`.

The tunable coupler wins clearly. With fixed coupling, faster gates and always-on ZZ crosstalk trade against each other.

**Results of the committed reviews** (simulation with CZ model version 4, `results/PROVENANCE.md`; see the design specs for
assumptions). Logical error per round, memory-Z / memory-X:

| Design | Code distance | Logical error per round | Notes |
|---|---|---|---|
| Alpha-1 v1.1 (AF-17) | 3 | **0.377 %** / 0.360 % | pulse-level CZ (36–49 ns), leakage model, leakage-removing reset; logical lifetime 246 µs = 3.3× best T1 |
| Alpha-1 v1.1 physics | 5 | 0.123 % (Z) | Λ(3→5) = 2.9 ± 0.2, from a fit over 5–15 rounds |
| Beta-1 rev C (BF-49) | 3 | 0.0772 % / 0.0718 % | |
| Beta-1 rev C (BF-49) | 5 | **0.0147 %** / 0.0128 % | Λ(3→5) = 5.3 ± 0.5 (Z), 5.6 ± 0.6 (X); logical lifetime 5.0 ms = 32× best T1 |
| Gamma-1 (GF-97) | 7 | **0.0032 %** / 0.0026 % | Λ(5→7) = 4.3 ± 0.6 (Z), 5.5 ± 0.7 (X); logical lifetime 26 ms = 163× best T1; impact floor 3.7e-08 per round |

Every Λ is a matched-model best case: the decoder is given the noise model that generated the errors. The ± values
are 1σ Monte Carlo errors. The Beta-1 and Gamma-1 rows use different estimators (a fit over 5–20 rounds, and one
14-round run), so compare within a spec, not across them.

## Limitations (read before trusting numbers)

* **CZ residual control error is an assumption.** The pulse simulation gives the coherent error (≈ 1e-4), leakage, and
  flux-noise terms. Pulse distortion, calibration drift, TLS defects and crosstalk are not simulated. They enter as one
  allowance: 1e-3 for Alpha-1 v1.1 and 5e-4 for Beta-1. This allowance dominates the CZ error and is the single most
  important number to measure on hardware. `fluxline.hpp` models the flux lines and says what the allowance asks of
  them (`results/flux_distortion_study.txt`), but the calibration does not use it: the allowance still stands for all of it.
* Pulse simulation is closed-system: T1 and Tφ during the gate are applied by the processor and QEC emulators, and the
  1/f flux noise at the excursion point is added analytically (quasi-static). Spectator crossings use a Landau-Zener
  estimate, not a many-qubit pulse simulation.
* Leakage is modelled at the frame level (flags, phase scrambling, transport, seepage, removal), not with 3-level
  dynamics in the QEC circuits. Heralded leakage (3-state readout) covers the measure qubits only, and the decoder
  uses fixed weight factors for heralded edges, not a posterior from the leakage model.
* **No layout in this repository is fabrication-ready.** All three chips are DRC- and connectivity-clean in the model;
  airbridges, EM and inductance are not verified, and no foundry DRC deck has been run. A review of the connectivity
  work found open electrical problems (shared ground bridge, flux lines that cannot tune the SQUIDs, coupler-dominated
  XY drive: `docs/design/mask-connectivity.md`). The in-house checks are:
  placement and island spacing, the 4 µm / 4 µm width/space check on the merged etch layer, and metal connectivity with
  modelled ground straps (airbridges). The control-line T1 budget (10× the design T1) is met on AF-17 and missed by the
  flux line on BF-49 and GF-97 (9.7×). The straps' capacitance to the coupler and bus
  strips they cross and the flux lines' inductive coupling to the SQUIDs are not modelled
  (`docs/design/mask-connectivity.md`).
* **The hardware-gate chain has only seen synthetic data.** The virtual chip is generated by the same simulator that
  evaluates it, so the rehearsal tests the plumbing and the gate logic, not the model's agreement with a real chip.
* **Impacts are a uniform burst with one recovery time.** `burst_floor` treats them as chip-wide; per-shot injection
  (`ImpactSpec`) adds a Gaussian spatial profile. The spread of phonons and the scrambling of two-level-system defects
  by an impact are not modelled. The impact rate is an assumed round figure. The veto discards whole shots; on
  Gamma-1, shot-by-shot discard does not lower the impact floor (the weak tail of a 25 ms burst still randomises the
  logical qubit). A latched hold of 5 recovery times leaves it unchanged too; from 7 recovery times on it falls below
  what the study resolves at 20,000 shots per severity (`results/decoder_study.txt`, `docs/design/impact-veto.md`).
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
* Not modelled: measurement-induced dephasing, readout and flux crosstalk.
* Pauli-frame QEC model: accurate for incoherent noise. Strong coherent ZZ (tens of kHz and up) is under-counted, so use
  `--statevector` to cross-check such designs.
* Setting up a memory experiment is linear in the number of rounds for the fault enumeration, but the matching
  decoder's all-pairs shortest paths cost the cube of the detector count (Gamma-1, d = 7: whole setup 1.4 s at 25
  rounds, 3.9 s at 50, on four threads).
* Results from before v0.4 with more than 22 defects in a shot used a greedy decoder fallback and were pessimistic. The
  design reviews in this repository have been re-run with the exact decoder, except the historical Alpha-1 v1.0 review.
* The statevector backend is practical to about 26–28 qubits; the Pauli-frame backend handles d = 5–7 easily.

## Roadmap

See [`ROADMAP.md`](ROADMAP.md) for programme phases and exit criteria. Next on the tooling side:

1. Full-wave/FEM cross-check of the cells; participation ratios → TLS-limited T1 per element.
2. Multi-qubit pulse simulation of spectator crossings; a net-zero (bipolar) CZ and a flux-line term in the calibration
   (the line model and its study are in `fluxline.hpp`).
3. An operation-count model of the streaming union-find decoder for the FPGA; for lattice surgery, routes with bends
   and a router for more than three patches, a correlated decoder fast enough for real time, and a layout for the
   patches and their routing region.
4. Data-qubit leakage heralding, and herald weights derived from the leakage model.
5. Materials loop: atomforge → DFT for junction barriers, oxides and TLS defects. Loss tangents → Q → T1.
6. Dear ImGui/ImPlot front-end; GPU kernels for 30+ qubit statevector emulation.

## Licence

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

This is proprietary software of Secured Cyber Solutions Ltd, licensed under the
[Secured Cyber Solutions Proprietary Software Licence v1.0](LICENSE). It is **not open source**:
no use, copying, modification or distribution is permitted without written permission from SCS.
Third-party components are listed in [NOTICE](NOTICE).

Licensing enquiries: https://securedcybersolutions.co.uk
