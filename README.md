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

**Alpha-1**, the first error-correction-ready design (17 qubits, distance-3 surface code), is in `designs/alpha1/`. Start with
[`ALPHA-1.md`](designs/alpha1/ALPHA-1.md). The junction physics and the tuned-anharmonicity correction are written up in
[`docs/WHITEPAPER.md`](docs/WHITEPAPER.md).

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
./build/test_qpuforge          # 71 validation checks (~10 min on 2 cores)
```

## CLI

```bash
qpuforge new chip.json --preset alpha1                      # Alpha-1 (add --distance 5 for the 49-qubit extension)
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
| Couplers | Static ZZ, dressed frequencies; Alpha-1 idle ZZ | Pair diagonalisation (5 levels). With `model_coupler`, a third transmon parked at the exchange null (3 levels) |
| Crowding | Exchange between near-degenerate qubits | 2nd-order Schrieffer–Wolff through every other qubit; parked tunable couplers block it |
| Coherence | T1 = dielectric (Q) + Purcell; Tφ = flux 1/f + charge + thermal photons | Koch 2007, Clerk–Utami; thermal population at T_eff |
| CZ gate | Duration, interaction point, feasibility, flux-pulse error | t = shape / (2√2 g) + ramps; \|11⟩–\|20⟩ uses α at the tuned EJ. Off-sweet-spot 1/f noise uses the short-time t² infidelity |
| Design rules | Degenerate pairs, \|11⟩–\|20⟩ resonance, ZZ > 1 MHz, infeasible CZ, parasitic resonance | Flags + collision count, used for fab yield |
| Inverse design | Targets (f01, α, g, χ) → C_shunt, R_n, C_c, resonator coupling | Iterated against the full capacitance-matrix model |
| Compiler | Layout, SWAP routing, 1q merging, {RZ, SX, X, CZ, measure, reset} | ≤ 2 physical pulses per 1q gate; per-qubit ASAP timelines cut into exact time slices |
| Processor | Every slice: native ops → control Pauli errors → ZZ conditional phases → resonant exchange → amplitude damping (T1 + thermal) → dephasing | Quantum trajectories, or exact density matrix (≤ 10 qubits) |
| Dynamic circuits | Mid-circuit measurement (collapse + assignment error), reset (+ reset error) | Classical record per shot |
| Error correction | Rotated surface code of any odd d; memory experiment; detector error model; decoder | Pauli-frame Monte Carlo with state-aware decoherence and coherent-ZZ windows; DEM from single-fault enumeration; exact minimum-weight matching |

## Validation (`tests/test_qpuforge.cpp`, 154 checks, all passing)

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

## Limitations (read before trusting numbers)

* **Two-level qubits** in the processor and QEC emulation: leakage to |2⟩ during CZ is folded into the control error, not
  simulated. CZ is modelled at gate level (duration plus errors), not as a pulse-level flux trajectory.
* **No dynamical decoupling.** Static-ZZ penalties for fixed-coupling designs are therefore upper bounds.
* Tunable couplers on Alpha-1 are a third transmon: direct g = 1 MHz, g_qc = 40 MHz, idle ZZ −4 to −14 kHz at the exchange
  null (3.79–3.97 GHz). The AF-9 study still treats residual ZZ as an input. g_max is the on-state scale, not a simulated pulse.
* Capacitances come from the lumped model. For a real layout they must come from EM simulation (Q3D, Palace, Sonnet).
* Not modelled: correlated errors (TLS drift, quasiparticle bursts, cosmic rays), measurement-induced dephasing, readout
  and flux crosstalk. Flux noise during the CZ is in the reported two-qubit error (about 2–5 × 10⁻⁵); the logical-error
  review of 24 November was not re-run for it.
* Pauli-frame QEC model: accurate for incoherent noise. Strong coherent ZZ (tens of kHz and up) is under-counted, so use
  `--statevector` to cross-check such designs.
* The statevector backend is practical to about 26–28 qubits; the Pauli-frame backend handles d = 5–7 easily.

## Roadmap toward tape-out

1. **Pulse level:** a multi-level transmon and coupler Lindblad simulation of CZ and DRAG pulses, giving leakage and coherent
   errors. Feed the resulting process matrices into the processor and QEC layers.
2. **EM loop:** import capacitance matrices and participation ratios from FEM of a GDS layout (Palace / KLayout).
3. **Leakage-aware QEC:** leakage tracking and leakage-reduction units in the memory experiment, plus correlated-error models.
4. **Cross-resonance architecture** (fixed-frequency) with its collision rules.
5. **Materials loop:** atomforge → DFT for junction barriers, oxides and TLS defects. Loss tangents → Q → T1 → logical error.
6. Dear ImGui/ImPlot front-end; GPU kernels for 30+ qubit statevector emulation.
