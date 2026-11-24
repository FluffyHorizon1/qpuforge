# qpuforge — superconducting QPU chip design & whole-processor emulation

Header-only C++17. You describe a chip (junctions, capacitors, couplers, readout resonators, materials).
qpuforge works out the device physics of every qubit and coupler. It then runs real quantum circuits
(OpenQASM 2.0 or built-in benchmarks) on a noisy model of the **whole processor**, so a change to the
chip design shows up directly as a change in algorithm performance.

```
chip.json ──► device physics ──────────────► noise model ──► compiled circuit ──► noisy emulation ──► metrics
 (Rn, C, Cc,    exact transmon spectra,         T1, Tphi,      layout, SWAP         statevector +       state fidelity,
  resonators,   C-matrix couplings, ZZ,         thermal, ZZ,   routing, RZ/SX/X/    quantum             GHZ fidelity,
  materials)    CZ time, coherence, rules       exchange,      CZ, time-scheduled   trajectories        XEB, counts,
                                                control, RO    moments                                  error budget
```

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
# or: g++ -std=c++17 -O3 -march=native -fopenmp -Iinclude src/qpuforge.cpp -o qpuforge
./build/test_qpuforge          # 54 validation checks
```

## CLI

```bash
qpuforge new chip.json --grid 3x3 --arch tunable --g-MHz 4      # inverse-design a chip from targets
qpuforge analyze chip.json                                       # device report + design-rule check
qpuforge run chip.json examples/ghz5.qasm --traj 400             # emulate a circuit
qpuforge run chip.json xeb:10                                    # also ghz:N qft:N bv:N
qpuforge budget chip.json xeb:10                                 # which noise source dominates?
qpuforge fab chip.json --samples 300                             # fabrication yield Monte Carlo
qpuforge run chip.json qft:6 --off zz,control                    # ablation: switch sources off
```

`new` options: `--grid RxC | --line N | --ring N`, `--f-hi/--f-lo` (GHz), `--alpha` (GHz), `--g-MHz`,
`--arch fixed|tunable`, `--Q` (dielectric quality factor), `--plan 2|4|8` (frequency allocation),
`--stagger-MHz`. Generated files are plain JSON: edit any qubit or coupler by hand to test a variant.

## Physics model

| Layer | What is computed | How |
|---|---|---|
| Junctions | E_J from normal-state resistance | Ambegaokar–Baratoff, SQUID flux tuning with asymmetry |
| Capacitance | E_C of every qubit, every coupling | Full Maxwell capacitance matrix, inverted (gives parasitic couplings too) |
| Transmon | f01, anharmonicity, charge dispersion, ⟨i\|n\|j⟩ | Exact diagonalisation in the charge basis (implicit QL) |
| Couplers | ZZ (static crosstalk), dressed frequencies | Exact diagonalisation of each coupled pair (5 levels each) |
| Crowding | Effective exchange between near-degenerate qubits | 2nd-order Schrieffer–Wolff through every other qubit (mediated coupling) |
| Coherence | T1 = dielectric (Q) + Purcell; Tφ = flux 1/f + charge + thermal photons | Standard formulas (Koch 2007, Clerk–Utami); thermal population at T_eff |
| CZ gate | Duration from coupling (\|11⟩↔\|20⟩), interaction point, feasibility | t = shape / (2√2 g) + ramps; tunable couplers parametrised by g_max and residual ZZ |
| Design rules | Degenerate pairs, \|11⟩–\|20⟩ resonance, ZZ > 1 MHz, infeasible CZ, parasitic resonance | Flags plus collision count; used for fab yield |
| Processor | Every moment: ideal native ops → control Pauli errors → ZZ conditional phases → resonant exchange → generalised amplitude damping (T1 + thermal) → dephasing | Quantum trajectories (exact unravelling), or exact density matrix for ≤ 10 qubits |
| Readout | Assignment errors (T1 decay during readout + separation error) | Applied to sampled bits |

The compiler does interaction-aware placement, shortest-path SWAP routing, merging of single-qubit gates, and
decomposition to {virtual RZ, SX, X, CZ} (at most 2 physical pulses per 1q gate). It then builds an ASAP
schedule in time "moments". Idle qubits decohere for the whole moment, the way they do on hardware.

## Validation (`tests/test_qpuforge.cpp`, all passing)

* 1q native decomposition: 3,016 unitaries, worst error 8e-16.
* Transmon spectra vs an independent numpy phase-grid solver (≤ 2e-4 GHz). Charge dispersion and
  ⟨0|n|1⟩ vs the Koch et al. asymptotics. E_J(10 kΩ, 180 µeV) = 14.04 GHz, matching a hand calculation.
* Two-transmon ZZ vs 2nd-order perturbation theory (within 5%). ZZ ∝ g² (ratio 3.997).
* Inverse design hits f01 within 1 MHz, α within 2 MHz and g within 0.05 MHz.
* Compiler: 48 random circuits (CX, CP, CCX, SWAP, arbitrary U) on 4 chip shapes reproduce the logical
  state to 1e-14 after routing and scheduling.
* Noise: T1 decay, Ramsey T2 and thermal equilibrium exact vs analytic. Resonant and detuned exchange
  match the Rabi formula to 1e-9 across 40 moments of random length.
* Trajectories vs exact density matrix on GHZ-5 and QFT-4 with every noise source on (within statistical error).
* GHZ metric = max over relative phase of ⟨GHZ_φ|ρ|GHZ_φ⟩. The sampled XEB estimator equals the exact density-matrix XEB.

## Example design study (`examples/design_study.cpp`, 3×3 chip)

| Design | t_CZ | mean \|ZZ\| | GHZ-9 | random circuits (10 cycles) | QFT-6 |
|---|---|---|---|---|---|
| fixed coupling g = 2 MHz | 275 ns | 13 kHz | 0.73 | 0.68 | 0.28 |
| fixed coupling g = 6 MHz | 98 ns | 116 kHz | 0.27 | 0.43 | 0.02 |
| fixed coupling g = 12 MHz | 54 ns | 462 kHz | 0.06 | 0.04 | 0.01 |
| **tunable coupler** (g_max 40 MHz, 5 kHz residual) | 23 ns | 5 kHz | **0.93** | **0.86** | **0.77** |

Other findings: the 2-frequency checkerboard has 14 resonance collisions and the 8-frequency plan has 0.
On short circuits, though, the 8-frequency plan's compressed detunings cost some ZZ. Dielectric Q above
about 3e6 (T1 ≈ 90 µs) gives diminishing returns, because control error then dominates. Fabrication with 2%
junction spread yields 97% collision-free chips; 5% spread yields only 63%. Full data is in `results/`.

## Limitations (read before trusting numbers)

* The processor emulator uses **two-level qubits**. Leakage to |2⟩ during CZ is not simulated; it is folded into the
  control error. The CZ is modelled at gate level (duration plus errors), not as a pulse-level flux trajectory.
* **No dynamical decoupling**, so static-ZZ penalties for fixed-coupling designs are upper bounds.
* Tunable couplers are parametric (g_max, residual ZZ), not a full coupler-transmon circuit.
* Capacitances come from the lumped model. For a real layout they must come from EM simulation (Q3D, Palace, Sonnet).
* Terminal measurement only (no mid-circuit measurement or reset). Unused spectator qubits are assumed in |0⟩.
* Flux-noise dephasing is an order-of-magnitude 1/f estimate. Dephasing during flux excursions is ignored.
* The statevector backend is practical to about 26–28 qubits on a large-memory machine. Validated here up to 9.

## Roadmap toward tape-out

1. **Pulse level:** multi-level transmon + coupler Lindblad simulation of CZ/DRAG pulses to get leakage and
   coherent errors. Feed per-gate process matrices to the processor layer instead of Pauli channels.
2. **Cross-resonance architecture** (fixed-frequency, IBM-style) with its collision rules.
3. **QEC:** mid-circuit measurement and reset, then a distance-3 surface-code memory experiment (17 qubits).
   The logical error rate vs design is the figure of merit that matters. Add a stabilizer simulator for larger distances.
4. **EM loop:** import capacitance matrices and participation ratios from FEM of a GDS layout (Palace / KLayout),
   so the flow runs design → layout → physics → processor metrics.
5. **Materials loop:** atomforge → DFT for Al/AlOx junction barriers, Ta/Nb surface oxides and TLS defects.
   Loss tangents → Q → T1 → processor metrics.
6. Dear ImGui/ImPlot front-end, and GPU kernels for 30+ qubit emulation.
