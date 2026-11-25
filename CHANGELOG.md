# Changelog

All notable changes to qpuforge. Versions follow `MAJOR.MINOR.PATCH`; before 1.0, a minor release may change the API, and
every such change is listed under **API changes**.

## [0.3.0] (unreleased)

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
- CLI: the `cz` and `gds` sub-commands, plus `qec --leakage-removal`.

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
