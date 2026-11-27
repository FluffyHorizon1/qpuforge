# Changelog

All notable changes to qpuforge. Versions follow `MAJOR.MINOR.PATCH`; before 1.0, a minor release may change the API, and
every such change is listed under **API changes**.

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
