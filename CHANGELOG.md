# Changelog

All notable changes to qpuforge. Versions follow `MAJOR.MINOR.PATCH`; before 1.0, a minor release may change the API, and
every such change is listed under **API changes**.

## [0.4.0] (unreleased)

Finishes the software side of Phases Alpha and Beta and starts Phase Gamma: an exact decoder at any number of defects,
memory experiments in both bases, correlated error bursts, hardware-gate evaluation from measured data, Beta-1 rev C and
the Gamma-1 (GF-97) distance-7 design. It also adopts thread-count-independent Monte Carlo and folds in the 0.3
follow-up fixes: JSON and GDSII hardening, and the scope of two PR #2 changes (the Alpha-1 v1.1 direct coupling and the
CZ flux-noise term). See the README "What's new in v0.4" table.

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
- Thread-count-independent Monte Carlo (`docs/design/deterministic-mc.md`). Every frame-simulation shot draws from its
  own generator, seeded by `frame_shot_seed(seed, shot)`, and partial sums are added in a fixed block order; the
  trajectory emulator adds its sampled XEB in trajectory order. Results are bit-identical for any `OMP_NUM_THREADS` and
  OpenMP schedule.
- `presets.hpp`: `alpha1_with_coupling(d, g_GHz)`. `device.hpp`: `ChipDesign::cz_flux_noise` (design JSON
  `gates.cz_flux_noise`).
- Tests: eleven new groups (`matching`, `fault_enumeration`, `memory_bases`, `bursts`, `characterisation`,
  `designs_v04`, and the serialisation groups `json`, `gds`, `json_rfc`, `json_depth`, `gds_write`), 321 checks in all
  (186 in 0.3.0), and a third review smoke test. Outside the validation checks: the `determinism` (1, 2, 3, 4 and 8
  threads in one process) and `determinism_env` (separate processes under `OMP_NUM_THREADS` = 1, 2, 4, 8) guards, the
  `fp_contract` and `fp_contract_flags` build-flag guards, and the `cli_version`, `cli_version_short` and `cli_usage`
  CLI guards.

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
- **Behaviour change:** the frame Monte Carlo seeds every shot separately. The same seed therefore gives a different,
  statistically equivalent sample than 0.3, and the same sample at any thread count. `openmp_threads` in the review
  provenance is informational.
- **Behaviour change:** `alpha1_v11()`, and through it `beta1()` and `gamma1()`, keep the 4 MHz direct qubit–qubit
  coupling (`alpha1_with_coupling`). Since PR #2 they had inherited 1 MHz from `alpha1()`, so the presets no longer
  reproduced the committed design files; they do again. `alpha1()` (Alpha-1 v1.0) stays at 1 MHz.
- **Behaviour change:** the CZ flux-noise term is added to `PairReport::err_2q` only when `ChipDesign::cz_flux_noise`
  is on: on for `alpha1()`, off for `alpha1_v11()`, whose pulse calibration carries its own flux terms. Design JSON
  without `gates.cz_flux_noise` follows `tunable_coupler.model_coupler`. Fixed-coupling and generic tunable designs get
  their pre-PR #2 `e_2q` back. The term is still computed and printed for every design. Neither change moves a logical
  error rate or circuit fidelity: the emulators take the control error and their own decoherence model, not `err_2q`.
- `Json::parse` rejects input that is not RFC 8259 JSON (see Fixed).

### Changed
- All review outputs were regenerated from the integrated code with CZ model version 2, the exact decoder and the
  thread-count-independent Monte Carlo (`results/PROVENANCE.md`).
- The Beta and Gamma reviews report both memory bases, the largest per-pair CZ leakage, and the mask-level DRC counts.
- The Beta gate adds a 100-round detection-rate stability criterion and an extraction criterion for the layout.
- Review smoke runs shorten the long experiments (setup time is quadratic in the number of rounds).

### Fixed
- **Error bars of fitted logical error rates were too small.** The reviews took the error of a fit over several round
  counts from the last round count alone. The slope of a free-intercept fit is less certain than that: by about 20 % for
  five round counts and 40 % for four. `fit_eps_per_round_err` propagates the binomial errors of all points through
  the fit; it is checked against the scatter of 4,000 synthetic fits. Every ± on a fitted value, and on a Λ derived
  from one, is wider than in 0.3. Single-run estimates were already correct.
- JSON: `Json::parse` accepts only RFC 8259 numbers and string escapes (`\uXXXX` is decoded; unknown escapes and raw
  control characters are rejected), and `dump` escapes control characters. Nesting deeper than 256 levels throws
  ("json: nesting too deep") instead of overflowing the stack.
- GDSII: `gds::write` throws if the output file cannot be opened or written.
- Build: the `fp_contract` guard is built with the same flags as the real targets (`QPUFORGE_COMPILE_OPTIONS`), so a
  build that drops `-ffp-contract=off` fails it.
- CLI: the usage text lists `--version` / `-V`.

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
- `device.hpp` (PR #2): an Alpha-1 three-transmon coupler model. With `ChipDesign::tc_model` (design JSON
  `tunable_coupler.model_coupler`, plus `g_qc_MHz` and `coupler_EC_GHz`), `analyze()` parks a coupler transmon at the
  exchange null (`coupler_exchange_null()`) and takes the idle ZZ from a three-level diagonalisation of qubit, coupler
  and qubit, instead of the typed `tc_residual_zz_kHz`. New `PairReport` fields: `coupler_f_GHz`, `coupler_J_MHz`.
- `device.hpp` (PR #2): a CZ flux-noise term. `analyze()` finds the flux of the higher qubit at the CZ interaction
  point and adds the short-time 1/f infidelity there (`flux_dephasing_rate()`, `flux_pulse_infidelity()`) to `e_2q`.
  New `PairReport` fields: `cz_flux`, `cz_df_dflux`, `cz_flux_err`. For Alpha-1 it is 2.0–4.6e-5 per CZ.
- CLI (PR #2): `qpuforge analyze` prints two new pair columns, `flux_e` (CZ flux-noise error) and `c_GHz` (coupler
  park frequency).
- `docs/WHITEPAPER.md` (PR #2): "Transmon anharmonicity and the Alpha-1 processor", with the AF17 layout figures
  `docs/af17-*.png`.
- Tests (PR #2): new `coupler_and_flux` CTest group (label `physics`).

### API changes
- **`Op::RESET` is a new enumerator, appended after `Op::X`.** All existing `Op` values keep their 0.2 numbering. A
  pre-release build of 0.3 put `RESET` between `MEASURE` and `BARRIER`, which renumbered `BARRIER`, `RZ`, `SX` and `X`.
  If you serialised `Op` as an integer with that pre-release, re-generate the data. Code that `switch`es over `Op`
  must handle `RESET`.
- **`compiler.hpp`: `lower(const Circuit&, meas&, bool& dynamic)`** is new. It reports whether the circuit has
  mid-circuit measurement feed-forward or resets. The 0.2 signature `lower(const Circuit&, meas&)` is still available
  as an overload that discards the flag.
- `NoiseModel::only()` accepts `"reset"` and `"leakage"`.
- `ChipDesign` gains `tc_model` and `tc_coupler_EC` (PR #2). Design JSON without `model_coupler` loads with the
  coupler model off.
- `ChipDesign` gains `coupler` (`CouplerElement`) and `leakage` (`LeakageSpec`). Old chip JSON files still load, with
  physical couplers off.
- **Behaviour change:** `LeakageSpec::model` defaults to `true`. As a result, the Pauli-frame QEC model applies
  `leak_1q` = 1e-6 per single-qubit gate, plus seepage, to old chips too. To get 0.2 behaviour, set
  `"leakage": {"model": false}`.
- **Behaviour change (PR #2): Alpha-1 direct coupling is 1 MHz (was 4 MHz).** `alpha1()` sets `g = 0.001`, so the
  link capacitance is 0.035 fF (was 0.14 fF), and it turns the coupler model on. Alpha-1 idle ZZ is now computed,
  −4.3 to −13.8 kHz, not the typed 5 kHz. `designs/alpha1/alpha1.json` and `examples/alpha1.json` were updated to
  match. `alpha1_v11()` and `beta1()` are built from `alpha1()`, so they inherit the 1 MHz coupling too (it was 4 MHz
  when their v0.3 review outputs were generated); see the open issue on Alpha-1 v1.1 / Beta-1 coupling.
- **Behaviour change (PR #2): `e_2q` from `analyze()` includes the CZ flux-noise term for every design whose higher
  qubit is tunable**, not just Alpha-1. Pulse-calibrated pairs (v1.1, Beta-1 reviews) overwrite `e_2q` with the
  calibrated value, so their review CZ errors are not affected, but abstract and fixed-coupling models are: mean `e_2q`
  of `gen0_surface17()` rises from 4.09e-3 to 5.98e-3, and of a 3×3 fixed-coupling grid at g = 3 MHz from 3.87e-3 to
  5.15e-3. Committed outputs generated before PR #2 (`results/design_study.*`, the Gen-0 row of the Alpha-1 ladder) do
  not include it.

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
- **Tests:** 186 checks (was 96; includes the PR #2 `coupler_and_flux` group), split into labelled CTest entries, plus two review smoke tests.

### Fixed
- `gds`: record lengths over 65,530 data bytes (polygons with more than 8,190 vertices) now throw instead of silently
  wrapping. Coordinates are range-checked. BGNLIB/BGNSTR carry the real date (or `SOURCE_DATE_EPOCH`).
- clang `-Wc++20-extensions` (structured binding captured in a lambda) in `calibrate_cz_pair`.
- `-Wunknown-pragmas` flood in builds without OpenMP.
- The test suite no longer writes a fixed `/tmp/qpuforge_test.gds`.

### Known issues
- PR #2 changes that are not yet reflected in committed outputs: the Alpha-1 v1.0 review (`designs/alpha1/alpha1_review.*`,
  `results/alpha1_review.*`) predates the coupler model, the 1 MHz coupling and the flux term; the v1.1 / Beta-1 outputs
  were generated at 4 MHz coupling and without the flux term in `analyze()`.
- `docs/WHITEPAPER.md` describes a tuned-α |11⟩–|20⟩ solve that `analyze()` does not do; `analyze()` uses the idle-α
  estimate `f_partner − α_idle`, about 2 MHz low.
- Alpha-1 v1.1 pair Q8–Q11 leaks 2.3e-4 per CZ, above the 2e-4 gate limit.
- The AF-17 and BF-49 first-cut layouts fail the merged-etch width/space DRC (4 µm / 4 µm, projection metric):
  AF-17 has 34 width + 202 space violations, and BF-49 has 98 width + 602 space violations. The masks are not
  fabrication-ready.
