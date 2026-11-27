# Changelog

All notable changes to qpuforge. Versions follow `MAJOR.MINOR.PATCH`; before 1.0, a minor release may change the API, and
every such change is listed under **API changes**.

## [0.4.0] (unreleased)

Finishes the software side of Phases Alpha and Beta and starts Phase Gamma: an exact decoder at any number of defects,
memory experiments in both bases, correlated error bursts, hardware-gate evaluation from measured data, Beta-1 rev C and
the Gamma-1 (GF-97) distance-7 design. See the README "What's new in v0.4" table.

### Added
- `matching.hpp`: `BlossomMatcher`, Edmonds' blossom algorithm for minimum-weight perfect matching (dense, primal-dual,
  O(n³), integer weights).
- `qec.hpp`:
  - Exact decoding at any number of defects. `DecodeMethod` selects `Auto` (dynamic programme up to 8 defects, blossom
    above), `DP`, `Blossom`, `Greedy` or `Legacy` (the 0.3 behaviour: exact up to 22 defects, greedy above).
  - Memory experiments in the X basis: `memory_experiment(S, rounds, basis)`; `prepare_memory(..., basis)`.
  - Correlated bursts from high-energy impacts: `BurstSpec`, `burst_noise`, `burst_floor`.
  - Syndrome and decoder statistics in `MemoryResult`: `detection_by_round`, `mean_defects`, `max_defects`,
    `frac_over_22`, `decode_us_per_shot`.
  - `FrameRunOptions`: run a memory experiment with another setup's decoder, or with a chosen decode method.
  - `enumerate_faults` follows only the qubits a fault touches, skips time slices without an operation on them, and
    runs in parallel: 8.5× faster on the Gamma-1 chip at 14 rounds on two cores. The 0.3 routine is kept as
    `enumerate_faults_reference`; the two agree fault by fault.
- `characterise.hpp`: the `qpuforge.characterisation.v0` data file, `model_from_measurement`,
  `evaluate_hardware_gate` (Alpha and Beta hardware gates from ROADMAP.md) and `synthesize_characterisation` (a virtual
  chip for rehearsal).
- `presets.hpp`: `plan_readout` (at most 9 resonators per feedline) and `gamma1()`. `device.hpp`: `f01_for_11_20`.
- CLI: `synth-chip`, `hwgate`, `qec --basis Z|X|both`, `qec --bursts`, `qec --decoder exact|legacy|greedy`,
  `new --preset gamma1`, `gds --synthesize --from`.
- Designs: Beta-1 rev C, Gamma-1 (GF-97), and the Alpha-1 v1.1 characterisation plan.
- `examples/gamma1_review.cpp`; `examples/review_util.hpp` (shared by the Beta and Gamma reviews); `tools/fill_spec.py`.
- Tests: six new groups (`matching`, `fault_enumeration`, `memory_bases`, `bursts`, `characterisation`, `designs_v04`),
  235 checks in all (186 in 0.3.0), and a third review smoke test.

### API changes
- **`MatchingDecoder::decode(defects, DecodeScratch&, DecodeMethod = Auto)`** replaces
  `decode(defects, std::vector<double>& cost, std::vector<int>& choice)`. `DecodeScratch` holds the work arrays of both
  matchers, and reports the matching weight and the method used.
- `prepare_memory(D, R, nm, rounds, basis = 'Z', build_decoder = true)` and
  `run_memory_frame(M, shots, seed, FrameRunOptions = {})` gain optional arguments. 0.3 calls compile unchanged.
- `memory_z_experiment(S, rounds)` is now a wrapper for `memory_experiment(S, rounds, 'Z')`. `MemoryExperiment` gains
  `detector_round` and `basis`.
- **Behaviour change:** shots with more than 22 defects are now decoded exactly. Logical error rates of long or large
  experiments drop; nothing changes where shots stay below 22 defects (distance 3 up to 25 rounds). Pass
  `DecodeMethod::Legacy` (`qpuforge qec --decoder legacy`) to reproduce 0.3 results on the same syndromes.
- **Behaviour change:** `beta1()` is rev C. Its readout plan has 6 feedlines (was 3), with resonators from 6.90 GHz in
  60 MHz steps. `gds` no longer overwrites a cell file it was asked to load.
- **Behaviour change:** `ChipReport` pairs: `cz_interaction_GHz` is solved at the anharmonicity the moving qubit has at
  the interaction point (`f01_for_11_20`), 1–2.5 MHz above the 0.3 estimate.
- `tools/fill_beta_spec.py` is replaced by `tools/fill_spec.py beta` (and `gamma`).

### Changed
- All review outputs were regenerated with CZ model version 2 and the exact decoder (`results/PROVENANCE.md`).
- The Beta and Gamma reviews report both memory bases, the largest per-pair CZ leakage, and the mask-level DRC counts.
- The Beta gate adds a 100-round detection-rate stability criterion and an extraction criterion for the layout.
- Review smoke runs shorten the long experiments (setup time is quadratic in the number of rounds).

### Fixed
- **Error bars of fitted logical error rates were too small.** The reviews took the error of a fit over several round
  counts from the last round count alone. The slope of a free-intercept fit is less certain than that: by about 20 % for
  five round counts and 40 % for four. `fit_eps_per_round_err` propagates the binomial errors of all points through
  the fit; it is checked against the scatter of 4,000 synthetic fits. Every ± on a fitted value, and on a Λ derived
  from one, is wider than in 0.3. Single-run estimates were already correct.

### Known issues
- CZ leakage: the worst frequency combination leaks 2.3e-04 per gate, above the 2e-4 limit of the Alpha
  simulation gate (mean 5.7e-05). Alpha-1, Beta-1 and Gamma-1 share the same 16 frequency combinations, so the
  pair class appears on all three chips (Q8–Q11 on AF-17).
- No layout passes the merged-etch width/space DRC (4 µm / 4 µm, projection metric): AF-17 has 34 width + 202 space violations,
  BF-49 rev C 98 width + 504 space violations, and GF-97 194 width + 1,016 space violations. The
  masks are not fabrication-ready.
- Setting up a memory experiment takes time quadratic in the number of rounds (single-fault enumeration).
- The hardware-gate chain has only been run on synthetic data from the virtual chip.

## [0.3.0] - YYYY-MM-DD


First tagged release. Adds pulse-level CZ calibration, a leakage-modelled QEC frame simulator, layout synthesis with
GDSII output, and the Alpha-1 v1.1 (AF-17) and Beta-1 (BF-49) designs. See the README "What's new in v0.3" table.

### Added
- `pulse.hpp`: a qubit–coupler–qubit charge Hamiltonian, piecewise-exact CZ propagation, and 3-knob calibration. Also
  flux-noise and spectator-crossing error terms, plus per-pair calibration of a whole chip with an on-disk cache.
- `qec.hpp`: a per-qubit leakage flag next to the Pauli frame. It models CZ leakage injection, transport, seepage,
  leakage-removing reset and data-qubit leakage removal (DQLR).
- `layout.hpp`, `gds.hpp`: method-of-moments capacitance extraction, parametric cells, Gauss–Newton cell synthesis,
  chip assembly, DRC and a GDSII writer.
- `core.hpp`: Householder + implicit-QL tridiagonal eigensolver (`sym_eigen_tri`).
- Circuits: mid-circuit `reset` (OpenQASM `reset`, `Circuit::reset`, `Op::RESET`) and dynamic-circuit detection.
- CLI: the `cz` and `gds` sub-commands, plus `qec --leakage-removal` and `--version`.

### API changes
- **`Op::RESET` is a new enumerator, appended after `Op::X`.** All existing `Op` values keep their 0.2 numbering. A
  pre-release build of 0.3 put `RESET` between `MEASURE` and `BARRIER`, which renumbered `BARRIER`, `RZ`, `SX` and `X`.
  If you serialised `Op` as an integer with that pre-release, re-generate the data. Code that `switch`es over `Op`
  must handle `RESET`.
- **`compiler.hpp`: `lower(const Circuit&, meas&, bool& dynamic)`** is new. It reports whether the circuit has
  mid-circuit measurement feed-forward or resets. The 0.2 signature `lower(const Circuit&, meas&)` is still available
  as an overload that discards the flag.
- `NoiseModel::only()` accepts `"reset"` and `"leakage"`.
- `ChipDesign` gains `coupler` (`CouplerElement`) and `leakage` (`LeakageSpec`). Old chip JSON files still load, with
  physical couplers off.
- **Behaviour change:** `LeakageSpec::model` defaults to `true`. As a result, the Pauli-frame QEC model applies
  `leak_1q` = 1e-6 per single-qubit gate, plus seepage, to old chips too. To get 0.2 behaviour, set
  `"leakage": {"model": false}`.

### Changed (release clean-up)
- **CZ model:** the default coupler truncation is now `Lc = 4` (was 3); qubit/coupler/qubit levels are 4/4/4.
  - Lc = 3 was not converged: it left a 4e-3 rad conditional-phase error. At Lc = 4, one more level changes the phase
    by 6e-5 rad.
  - Per-pair calibration is about 1.6× slower.
  - All CZ calibrations and review outputs were regenerated. Headline changes are within Monte Carlo error. One
    exception: Alpha-1 v1.1 pair Q8–Q11 now leaks 2.3e-4 per CZ (was 1.3e-4), above the 2e-4 gate limit; the mean is
    5.7e-5.
- **CZ cache:** the cache is versioned (`CZModel::kVersion = 2`).
  - Keys include the time step, optimiser budget, level truncations and charge cutoffs.
  - The file is wrapped as `{schema, model_version, entries}`; unversioned or other-version files are ignored with a
    warning and re-simulated.
  - `results/cz_cache.json` was regenerated by the review runs.
- **Wording:** "leakage-aware QEC" is now "leakage-modelled QEC". Leakage is simulated, but the decoder does not use it.
- **DRC:** new minimum width/space check on the merged etch layer (4 µm, projection metric).
  - The AF-17 and BF-49 first-cut layouts do not pass it. "DRC clean" claims were reworded.
  - The Beta exit criterion now covers placement and island spacing only.
- **QEC:** `summarize()` warns when fault mechanisms flipping more than two detectors are missing from the matching
  graph (new `hyperedge_p`), and when faults flip the logical without any detector.
- **Specs:** BETA-1 caveats added: no leakage herald in the decoder; Λ is a matched-model best case;
  `residual_ctrl_err` is an assumed allowance.
- **Review programs:** they take `--out`, `--cache` and `--smoke`, and write a provenance block (see
  `results/PROVENANCE.md`). Review outputs are generated by code only; there are no hand-written additions.
- **Tests:** 186 checks (was 96), split into labelled CTest entries, plus two review smoke tests.

### Fixed
- `gds`: record lengths over 65,530 data bytes (polygons with more than 8,190 vertices) now throw instead of silently
  wrapping. Coordinates are range-checked. BGNLIB/BGNSTR carry the real date (or `SOURCE_DATE_EPOCH`).
- clang `-Wc++20-extensions` (structured binding captured in a lambda) in `calibrate_cz_pair`.
- `-Wunknown-pragmas` flood in builds without OpenMP.
- The test suite no longer writes a fixed `/tmp/qpuforge_test.gds`.

### Known issues
- Alpha-1 v1.1 pair Q8–Q11 leaks 2.3e-4 per CZ, above the 2e-4 gate limit.
- The AF-17 and BF-49 first-cut layouts fail the merged-etch width/space DRC (4 µm / 4 µm, projection metric):
  AF-17 has 34 width + 202 space violations, and BF-49 has 98 width + 602 space violations. The masks are not
  fabrication-ready.
