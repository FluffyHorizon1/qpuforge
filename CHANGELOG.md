# Changelog

All notable changes to qpuforge. Versions follow `MAJOR.MINOR.PATCH`; before 1.0, a minor release may change the API, and
every such change is listed under **API changes**.

## [Unreleased]

Clears the merged-etch width/space DRC of all three masks, makes their metal connectivity clean and gates it, and
brings every CZ pair class within the 2e-4 leakage limit of the Alpha gate. Each qubit SQUID now lands on a ground
lead, each qubit flux line is shorted to the plane, modelled airbridges ("ground straps") tie the enclosed ground to the
plane, and the XY line is moved out of the qubit's cut-out after an extraction of the control-line couplings showed a
sub-microsecond T1 limit (docs/design/mask-connectivity.md). The masks are DRC- and connectivity-clean in the model (2D
LVS-lite connectivity with modelled airbridges, treated as ideal ground); airbridges, EM and inductance are not
verified, and no foundry DRC deck has been run. **A review found three high-severity electrical problems with this
approach (shared ground bridge, flux lines that cannot tune the qubits, coupler-dominated XY drive): see Known issues.**

The decoding items of Phase Gamma: a union-find decoder, leakage-aware decoding from heralded leakage, per-shot impact
injection with an impact veto, and linear-time fault enumeration. Study: `examples/decoder_study.cpp` ->
`results/decoder_study.{txt,json}`.

A streaming union-find decoder for real-time decoding: it decodes a sliding window of rounds, commits the corrections
in the window's first rounds and moves on (`qpuforge qec --decoder uf-stream`). `decoder_study` section 5 compares it
with the batch decoder on the same syndromes.

Phase Gamma logical operations: lattice surgery between surface-code patches (logical Z_A Z_B and X_A X_B measurement,
and a three-patch CNOT). Each operation is checked noiselessly with a stabiliser tableau and simulated with the
Pauli-frame Monte Carlo. Study: `examples/surgery_study.cpp` -> `results/surgery_study.{txt,json}` (CZ model version 4).

Lattice surgery, next steps: one decoder pass that returns all observables of an operation (exact matching or
union-find on one graph whose edges carry an observable mask); a two-pass correlated decoder that uses the detectors
of both Pauli types; operations between patches that are not adjacent, through a routing strip; and a chip account for
two and three d = 5 patches with their routing region. `surgery_study` sections 4 to 6.
The next decoding items of Phase Gamma: herald weights derived from the leakage model instead of fixed factors,
heralding of data-qubit leakage through DQLR (measured, and found to cost more than it recovers), and an
operation-count model of the streaming union-find decoder with two stated hardware cost models. The hardware part is a
model: nothing was synthesised or run on an FPGA. Study: `examples/decoder_study2.cpp` ->
`results/decoder_study2.{txt,json}`.

A strict command line: each command accepts only its own options, and an unknown, misspelt or repeated option, a
missing value or a malformed number is a usage error (exit 1), reported before any file is read.

The qpuforge side of the atomforge → qpuforge materials contract (`docs/DATA_CONTRACT.md`): a reader for the
`forge.materials.v0` sidecar that atomforge writes, and a CLI flag that applies one before the device report.

A model of the flux lines and of their digital correction (`fluxline.hpp`), and a study of what the CZ gates of
Alpha-1 v1.1 and Beta-1 need from them (`examples/flux_distortion_study.cpp` -> `results/flux_distortion_study.{txt,json}`,
`docs/design/flux-line-predistortion.md`). The calibration does not use the model: no cache entry or review number changes.

### Added
- `fluxline.hpp`:
  - `FluxLine` (first-order sections, `FluxSection`, in cascade with the step response 1 + a exp(-t / tau), a
    single-pole low-pass, a delay), `FluxLineFilter` (its response to a held drive in closed form), `FluxAwg` (sample-and-hold, DAC resolution),
    `FluxChannel`, `Predistorter` (digital first-order inverses of a fit of the line), `flux_line_output` (the flux a
    channel delivers for a train of pulses).
  - `CZDrive`, `simulate_cz_drive` (the propagator and fidelity of `simulate_cz` for any mover and coupler flux),
    `cz_ideal_drive` (reproduces `simulate_cz`), `CZFrame`, `cz_frame`, `cz_fidelity_in_frame` (the fidelity with virtual Z
    gates calibrated in another context), `cz_offset_response` (added error per unit of constant flux offset, with free
    and with fixed virtual Z gates).
  - `FluxSetup`, `cz_gate_drive`, `cz_train_drive` (a gate on two channels, on its own or inside pulse trains),
    `cz_matched_command` and `flux_amplitude_scale` (the amplitude calibration on a distorted line), `cz_system_of` (the
    pair system of a stored calibration), `cz_starts` (CZ start times of a compiled program).
- `pulse.hpp`: `pair_system` (the system `design_pair_system` builds for a given direct coupling; no change of arithmetic).
- `examples/flux_distortion_study.cpp` and its smoke test (`smoke_flux_distortion_study`).
- Tests: group `flux_predistortion` (label `pulse`, 32 checks, about 3 s). 541 validation checks in all; the committed
  `results/test_output.txt` predates the group and holds 500.
- q17, FEM reference cross-check scaffold (decided by Finnian Bond 2026-11-28: both Elmer and Palace; deviation
  reported first, tolerance enforced later). `fem.hpp` (not in `qpuforge.hpp`; include it with `layout.hpp`):
  `FemReference`, the `qpuforge.fem_crosscheck.v0` reader / writer (`reference_from_json`, `reference_to_json`,
  `load_reference`), `import_elmer` (`CapacitanceMatrix.dat`, lumped → Maxwell), `import_palace` (`terminal-C.csv`),
  `FemTolerance` (`enforced = false`), `compare_fem`, `circuit_of`, `comparison_summary`, `comparison_json`.
  `docs/DATA_CONTRACT.md` §6. Test group `fem` (label layout) with synthetic fixtures in `tests/data/fem/`. No solver
  is a build dependency; no committed number changes.
- q18, opt-in flux-line distortion term (docs/design/flux-line-distortion.md; decided by Finnian Bond 2026-11-28:
  additive to `residual_ctrl_err`, off by default, no regeneration). `device.hpp`: `FluxTail`, `FluxLineSpec`,
  `ChipDesign::flux_line` (design JSON `tunable_coupler.flux_line`, written only when enabled). `pulse.hpp`:
  `flux_line_eps_eff`, `flux_distortion_err`, `CZCalibration::distortion_modeled`, `flux_eps_eff`, `err_distortion`
  (added to `err_ctrl`, hence `err_2q_ctrl`, when on; `cz_calibration_json` writes the two values only when on).
  `qpuforge cz` prints the term when on. Off, every output is unchanged; the CZ cache key and `CZModel::kVersion` are
  unchanged. Test group `flux_line` (label pulse; its checks are not in `results/test_output.txt`, which is not
  regenerated).
- `layout.hpp`:
  - `detail::place_ground_straps` (`detail::Strap`, `detail::IslandIndex`, `detail::island_index`,
    `detail::foot_comp`, `detail::UnionFindLite`; `kStrapW`, `kStrapMaxSpan`, `kStrapLinePenalty`,
    `kStrapCouplerPenalty`): ground straps (modelled airbridges, GDS layers 4 = span, 5 = feet) that tie every floating
    ground fragment and every shorted-line pocket to the ground plane, chosen by Kruskal over candidate straps. A strap
    never spans a qubit; it may span lines, coupler bodies and bus strips. `ChipLayout::straps`.
  - Unreachable ground fragments smaller than `kRemoveFragmentUm2` (2500 µm²) are etched away
    (`ConnectivityReport::fragments_removed`).
  - `ConnectivityReport::straps`, `straps_landed`, `strap_lines_crossed`, `strap_couplers_crossed`; `clean()` now also
    requires every strap to land.
  - `qubit_control_lines` (`kXYEnd`, `kXYLevel`), `line_coupling` (`LineCoupling`), `line_budget` (`LineBudget`,
    `kLineT1Margin`), `line_budget_summary`, `line_budget_json`: the qubit's capacitance to its XY line and to its
    shorted flux line (each panel weighted by sin(βz) from the short), and the T1 limit of each line at every qubit's
    idle frequency against its design T1. The reviews record it under `layout.control_lines`.
  - `extract_capacitance(..., excite, panel_charge)`: optional per-panel charges (`PanelCharge`).
  - `kSquidGap`, `kTongueW`: the qubit SQUID's span and ground-lead width.
- Gate criterion `connectivity` in the three reviews; `drc` now includes the mask width/space check; the control-line
  budget is printed and recorded (`layout.control_lines`), not gated.
- Tests: `layout_geometry` checks the strap placer (a pocket gets one strap, never over a qubit; a pocket ringed by a
  qubit gets none); `layout_presets` requires each committed chip to be connectivity-clean, with one junction per qubit
  and coupler and one flux line per qubit, and DRC 0 + 0. Four checks more than before (one and three).
- `docs/design/mask-connectivity.md`.
- `layout.hpp` (earlier in this cycle):
  - `detail::fill_ground_slivers`: mask clean-up that merges two cut-outs separated by a strip of ground narrower than
    the minimum (projection metric). It removes ground only; a gap that touches island metal is left to the DRC.
  - `detail::guarded_fill`: the clean-up with a connectivity guard: a fill that would cut ground off the ground plane
    is undone, and the DRC reports the neck it leaves. (`assemble_chip` no longer uses it: see Changed.)
  - `detail::metal_map` (`MetalMap`): exact connected components of the metal (outline minus etch) by scanline and
    union-find.
  - `ConnectivityReport` (`ChipLayout::conn`), `connectivity_summary`, `connectivity_json` (LVS-lite): conductors
    sharing metal, conductors on the plane, shorted lines (resonators, qubit and coupler flux lines) reaching the plane,
    floating ground fragments, junctions joining an island to the plane, ground straps. The reviews record it under
    `layout.connectivity` and gate it.
  - `detail::etch_of` (cut-outs minus metal, factored out of `assemble_chip`); `cell_ground_holes` (the same clean-up
    on the ground cut-outs of the extraction cells).
  - `set_bus_knob` / `bus_knob`, `kBusRetract`, `CouplerGeom::g_b` (`bus_end_gap_um` in the cell file),
    `CouplerGeom::bus_gap()`, `CouplerGeom::bus_tbar()`, `LayoutTargets::G_md_renull_margin`.
- `pulse.hpp`: `CZModel::leak_limit`, `CZModel::coh_limit`, `CZModel::leak_weight`;
  `optimize_cz(..., leak_weight, leak_limit, coh_limit)`.
- `alpha1_v11_review`: the Alpha simulation gate computed from per-pair maxima (`exit_criteria`), the extraction errors,
  the re-null check and the connectivity report.
- Tests: `layout_geometry` (bus knob, shorted staircase end, ground-sliver fill, metal map, guarded fill) and
  `layout_presets` (the three committed chips: placement clean, DRC 0 + 0, connectivity clean, extraction on target,
  exchange re-nullable; label `layout`, a few minutes). `pulse_cz` checks the two pair classes with the most leakage with and without the fallback.
- `docs/design/cz-leakage.md`.
- `unionfind.hpp`: `UnionFindDecoder`, Delfosse-Nickerson union-find with weighted growth and peeling
  (docs/design/union-find-decoder.md), with its work arrays `UFScratch`. `DecodeMethod::UnionFind`,
  `MatchingDecoder::union_find()`, `DecodeScratch::uf`; `qpuforge qec --decoder uf`.
- `qec.hpp`:
  - Heralded leakage (docs/design/heralded-leakage.md): `HeraldSpec` in `FrameRunOptions` (3-state readout of the
    measure qubits with a herald fidelity and a false-herald rate), `MatchingDecoder::decode_scaled` (exact matching on
    a per-shot re-weighted graph), `herald_touch_edges`, `herald_footprints`, `MemorySetup::herald_edges`.
    `qpuforge qec --herald [--herald-fidelity F] [--herald-false P] [--herald-factor X] [--herald-rule footprint|detector]`.
    The herald options need `--herald`; `--herald-fidelity`, `--herald-false` and `--herald-factor` are in [0, 1].
  - Impact injection and veto (docs/design/impact-veto.md): `ImpactSpec` and `VetoSpec` in `FrameRunOptions`,
    `BurstSpec::radius`, `impact_shot_seed`, `calibrate_veto`, `burst_floor_veto`. `qpuforge qec --bursts --veto
    [--veto-window R] [--veto-radius L] [--veto-false P] [--impact-p P] [--burst-radius L]`. `--veto` needs `--bursts`
    and these options need `--veto`; `--veto-false` and `--impact-p` are in [0, 1], `--veto-window` >= 1, and
    `--veto-radius` and `--burst-radius` >= 0.
  - `MemoryResult`: herald counts, impact and discard counts, `p_logical_kept`, `eps_kept`, `veto_hist`.
  - `MemorySetup::sim_xy`, `det_xy` (lattice positions); `FrameProgram::FMoment::known` (known basis state per slice).
  - `enumerate_faults_direct`: the 0.4 sparse enumeration, kept for timing comparisons.
- Tests: groups `union_find`, `heralds`, `impact_veto` (label qec); `fault_enumeration` now also checks the linear-time
  version and its scaling (a deterministic work count: the optional `steps` output of both enumerations); the
  thread-count determinism digest covers union-find, heralds, and impacts with the veto.
- `unionfind.hpp`: `StreamingUnionFind` and `StreamScratch`, sliding-window decoding over rounds with a commit region
  (docs/design/streaming-union-find.md); `FrameRunOptions::stream`; `DecodeScratch::stream`;
  `MemoryResult::stream_window_us_mean` and `stream_window_us_max`. `qpuforge qec --decoder uf-stream [--window W]
  [--commit C]` (6 and 3 rounds by default; both need `--decoder uf-stream`, and `--commit` must be smaller than
  `--window`). Test group `streaming_uf`; `decoder_study` section 5 (d = 7, 50 rounds).
- `surgery.hpp`:
  - Multi-patch geometry: `SurgeryLattice` (patches and merged regions on one qubit lattice), `region_faces`.
  - Schedules: `SurgerySchedule` / `SurgeryPhase`; `build_surgery` writes the circuit and derives the detectors across
    merges and splits, the reference bits and the classical-bit bookkeeping.
  - Protocols: `surgery_parity` (logical Z_A Z_B or X_A X_B measurement), `surgery_cnot` (three-patch CNOT with
    Pauli-frame corrections), `surgery_memory` (one patch), `surgery_experiment` (by name).
  - `StabilizerTableau`, `tableau_records` and `check_surgery_noiseless`: noiseless verification of detectors and
    observables on random measurement branches.
  - `surgery_chip`: a chip for the lattice with Alpha-1 v1.1 or Beta-1 physics.
  - `prepare_surgery`, `surgery_fault_signatures`, `run_surgery_frame`: one detector error model and exact matching
    decoder per observable, and the Pauli-frame Monte Carlo (thread-count independent).
- CLI: `qpuforge surgery [--op zz|xx|cnot] [--d 3,5] [--input ..] [--readout Z|X] [--chip beta1|alpha1_v11] ...`.
  `--op`, `--readout` and `--chip` are choices, `--d` takes odd distances from 3 to 15, and `--readout` applies to
  `--op cnot` only. Exit status 3 when a noiseless check fails.
- `examples/surgery_study.cpp` and `results/surgery_study.{txt,json}`; `docs/design/lattice-surgery.md`.
- Tests: `surgery` and `surgery_cnot` groups (82 checks); the determinism guard also covers
  `run_surgery_frame` (the ZZ measurement and the three-patch CNOT); a smoke entry for the study; `cli_*surgery*` guards
  for the command line (outside the validation checks).
- `surgery.hpp`, joint decoding: `SurgeryFault` and `surgery_faults` (each fault mechanism with its decoded detectors,
  a mask of the observables it flips and, on request, its detectors of the other Pauli type); `JointEdge`, `JointDEM`
  and `build_joint_dem` (one edge per detector pair: the probability of all its mechanisms, the mask of the most
  likely class, and counts of the edges and of the probability that one mask cannot represent); `SurgerySetup::jdem`
  and `SurgerySetup::joint`, a `MatchingDecoder` whose `decode()` returns the mask of all observables from one
  matching, or from one union-find pass with `DecodeMethod::UnionFind`.
- `surgery.hpp`, correlated decoding: `prepare_surgery(..., correlated = true)` builds a union-find decoder on the
  detectors of the other Pauli type and, for each of its edges, the joint-graph edges that share mechanisms with it.
  `SurgeryDecoder::Correlated` decodes the other type first, lowers the weights of the correlated edges to their
  conditional probability and runs exact matching on the joint graph.
- `surgery.hpp`: `SurgeryDecoder` (`Separate`, `Joint`, `SeparateUF`, `JointUF`, `Correlated`),
  `surgery_decoder_name`, and `run_surgery_decoders`, which decodes the same shots with several decoders and returns
  for each the error per observable, the shots with any observable wrong, the shots that differ from, are fixed or are
  broken relative to the first decoder, and the decode time per shot. `SurgeryResult` gains `decoder`, `decode_ns`,
  `differ`, `any_fixed` and `any_broken`.
- `surgery.hpp`, routing: `surgery_parity`, `surgery_cnot`, `surgery_experiment` and `surgery_cnot_reference` take a
  `gap` (odd, default 1). With `gap` > 1 the patches are joined through a routing strip, d wide and `gap` long, which
  is merged with both patches for the measurement and read out at the split. Strips are straight.
- CLI: `qpuforge surgery --gap G` (odd, 1 to 15) and `--decoder separate|joint|uf|joint-uf|correlated|all`; `all`
  decodes the same shots with every decoder and prints decode times. Both are in the option table, the usage and the
  `cli_*` guards (`cli_surgery_even_gap`, `cli_surgery_bad_decoder`, `cli_ok_surgery_routing`, `cli_exit_status`).
- `examples/surgery_study.cpp`: section 4 (the shots of section 2 decoded with every decoder), section 5 (noiseless
  check and logical error against the strip length, Beta-1 physics, strips of 1, 3, 5 and 7 at d = 3 and 5) and
  section 6 (chip account for d = 5 patches and their routing region: qubits, couplers, feedlines, lines, design
  frequencies, collisions, fabrication Monte Carlo). `results/surgery_study.{txt,json}` regenerated; sections 1 to 3
  are unchanged.
- Tests: `surgery_joint` (57 checks) and `surgery_routing` (51 checks) groups. The determinism guard also covers
  `run_surgery_decoders` with every decoder and a routed XX measurement.
- Tests: `cli_*` guards (label `unit`) for unknown options, missing values, malformed numbers and choices, repeated and
  extra arguments, one accepted invocation per command with every option given (two for `new`, a topology and
  `--preset`; three for `qec`: the burst options, then the decoding options and the streaming decoder at non-default
  values, matched against what the run echoes), the materials options, failed writes and the exit statuses;
  `cli_exit_status` also checks every option dependency and value range of `qec` and `analyze` on a chip file that does
  not exist (exit 1 before any file is read).
- `materials.hpp`: `load_materials` / `materials_from_json` read and check a `forge.materials.v0` sidecar, and
  `apply_materials` writes it into a `ChipDesign`. The sidecar's `Q_dielectric` (unit `dimensionless`) sets
  `ChipDesign::Q_dielectric` and its `junction_gap_ueV` (unit `ueV`) sets `ChipDesign::gap_ueV`; a value must be
  finite and > 0. With `MaterialsApplyOptions::use_loss_tangent`, a `loss_tangent` sets
  `Q_dielectric = 1 / loss_tangent` (the uniform-loss model of `analyze()`) when the sidecar's `Q_dielectric` is null
  or absent. A null value, such as each of atomforge's placeholders, changes nothing. `Rn_spread_pct` and the
  environment fields are never taken from a sidecar. The returned `MaterialsReport` has one line per quantity, applied
  or skipped, with the reason. A wrong schema, an unknown quantity, a wrong unit, a bad value or a mismatched `maps_to`
  throws `std::runtime_error` before the design is changed.
- CLI: `analyze --materials <sidecar.json> [--loss-tangent]` applies a sidecar, prints one line per quantity, then
  prints the device report. The chip JSON file is not changed. `--loss-tangent` without `--materials`, and an empty
  `--materials` value, are usage errors (exit 1), reported before any file is read.
- CLI: `--materials <sidecar.json> [--loss-tangent]` is also accepted by `run`, `budget`, `fab`, `qec`, `cz`, `synth-chip`
  and `hwgate`, through one shared loader, so a sidecar value reaches the emulators and the logical-error runs and is no
  longer silently dropped by a sub-command that did not know the flag. `layout` and `gds` declare it only to refuse it
  (a usage error, exit 1; they draw the chip from the chip JSON). When a sidecar changes `junction_gap_ueV`, the report line says that each qubit keeps its
  `Rn_kOhm`, so `E_J` and `f01` scale with the gap.
- Tests: the `materials` group (label `unit`, 25 checks), with sidecar fixtures in `tests/data/`, among them
  atomforge's placeholder sidecar byte for byte. Outside the validation checks: the `cli_materials`,
  `cli_materials_placeholder`, `cli_materials_loss_tangent`, `cli_materials_loss_tangent_opt_in`,
  `cli_materials_refused` and `cli_usage_materials` CLI guards, plus `cli_materials_fab`, `cli_materials_run`,
  `cli_materials_budget`, `cli_materials_run_refused`, `cli_materials_layout_refused`, `cli_materials_gds_refused`,
  `cli_materials_qec`, `cli_materials_cz_missing_chip`, `cli_materials_hwgate_missing_chip` and `cli_loss_tangent_alone`
  for the other sub-commands. `test_qpuforge` is built with `QPUFORGE_TEST_DATA`, the path of `tests/data`.
- `qec.hpp`, herald weights derived from the leakage model (docs/design/derived-herald-weights.md):
  - `HeraldSpec::derived`; `HeraldModel` (`MemorySetup::herald_model`), `build_herald_model`, `prepare_herald_model`.
    Per heralded readout: the probability that the qubit was leaked, and the probability of a fault on each
    decoder-graph edge given that it was. Weights are log-likelihood ratios; an edge more likely than not is applied.
  - `frame_detail::leak_chain` (`LeakChain`): one leakage excitation followed through the circuit with the sampler's
    transitions (forward pass, backward smoothing per readout, forward continuation).
  - `fault_signatures(F, E, list, steps)`: the linear-time enumeration for any list of mechanisms.
  - `qpuforge qec --herald --herald-weights fixed|derived` (`--herald-factor` and `--herald-rule` apply to `fixed`).
- Data-qubit leakage heralding (docs/design/data-leakage-herald.md): `LeakageSpec::dqlr_transfer` (`dqlr_transfer` in
  the chip JSON, written only when set), `NoiseModel::dqlr_transfer`, `FrameProgram::dqlr_transfer`,
  `FrameProgram::FMoment::dqlr_to`. The DQLR pulse leaves the leakage on the measure qubit, and that qubit's next
  3-state readout heralds it. `qpuforge qec --leakage-removal lrr+dqlr --herald --dqlr-herald`.
- `unionfind.hpp`, for the hardware model (docs/design/streaming-union-find-fpga.md):
  - `UFOpCounts` and `UnionFindDecoder::decode_counted`: loads, growth steps, edge visits and updates, finds, pointer
    hops, unions, peeling steps, corrections. The plain `decode()` is the same template with an empty counter.
  - `StreamOptions` for `StreamingUnionFind::decode`: `scaled` (heralds through the streaming decoder), `predecode`
    (an isolated-pair pre-decoder), `count` (per-window counts in `StreamScratch::counts` and `cleared`).
    `StreamingUnionFind::global()`, `window_vertices`, `window_edges`, `committed_rounds`.
  - `UFCostModel`, with `serial_bram()` and `parallel_growth()`: cycles per window from the counts.
- `qec.hpp`: `run_stream_profile` and `StreamProfile`: histograms per window of every counter, of the cycles per
  committed round and of the queueing latency, for a list of cost models.
- `examples/decoder_study2.cpp` and `results/decoder_study2.{txt,json}`; `docs/design/derived-herald-weights.md`,
  `data-leakage-herald.md`, `streaming-union-find-fpga.md`.
- Tests: groups `herald_weights` (7 checks), `dqlr_herald` (6) and `stream_ops` (5), label `qec`: 18 checks on top of
  the 500 below (`results/test_output.txt` predates them). The thread-count determinism digest also covers the herald
  model, derived weights through the batch and the streaming decoder, DQLR with heralding, and the operation profile.
  CLI guards for the two new options, and `cli_herald_stream` (below).
- 500 validation checks in all (`results/test_output.txt`): 346 on `main` before the decoding work (the 321 of 0.4.0
  plus 25 in `materials`), plus 21 from the decoding work (#20), 48 from fabrication readiness (#19), 3 from the
  streaming decoder (#27) and 82 from lattice surgery (#21). The 0.4 groups now give 328: fabrication readiness added
  four `pulse_cz` checks and one `layout` check, the decoding work two `fault_enumeration` checks. The strict CLI (#24)
  adds CTest guards, not validation checks.
- `docs/DATA_CONTRACT.md`: an "Implementation status" section. The proposal text is unchanged.
- `materials.hpp`: `provenance` is an enumerated field on `forge.materials.v0` (decision D4, Finnian Bond,
  2026-11-28): `placeholder`, `assumed`, `literature`, `measured`, `dft-loss-tangent`, `participation-weighted`. A
  number must carry one, and `literature` / `measured` must name their source in `reference` (new member
  `MaterialsQuantity::reference`, printed on the report line). Values must lie in a physical range (`Q_dielectric`
  1e3–1e10, `loss_tangent` 1e-10–1e-3, `junction_gap_ueV` 10–3000 µeV). `produced_by.tool` is required, every
  `placeholders` entry must name a quantity in the file, and `contains_placeholders` must agree with the quantities.
  atomforge's writer already meets every rule. The `materials` group has 29 checks (4 new), 504 validation checks in
  all; `results/test_output.txt` still records the 500-check run of `main` named in `results/PROVENANCE.md`.
- CLI: `qpuforge materials apply <chip.json> <sidecar.json> --out <new-chip.json> [--loss-tangent]` (decision D2) writes
  a copy of the chip with the sidecar's accepted values and a `materials_source` block: the sidecar's file name and
  SHA-256, its `produced_by` and `system`, the loss-tangent opt-in, the gap policy and one entry per quantity (applied,
  before, after, provenance, reference). Each qubit keeps its `Rn_kOhm` (D1), so a new gap moves `f01`, as it does in
  memory. A chip that already has a `materials_source` keeps it as `previous`. `design_from_json` ignores the block, so
  the copy loads anywhere; `--out` naming an input file, a missing `--out` or another action than `apply` is refused,
  and nothing is written when a check fails. `materials.hpp` gains `materials_apply_file`, `materials_source_json`,
  `MaterialsApplyResult` and `materials_detail::sha256_hex` / `read_file`. The `materials` group has 34 checks (5 new),
  509 in all; CLI guards `cli_materials_apply*` and new `cli_exit_status` / `cli_help_command` lines.
- The `layout` / `gds` refusal of `--materials` now points at `qpuforge materials apply`.
- `designs/alpha1-v11/materials-alpha1-v11-handwritten.json` (decision D5, first half): the first sidecar with values,
  written by hand. `Q_dielectric` 3e6 and `junction_gap_ueV` 180 are the values Alpha-1 v1.1 already assumes, each with
  `provenance: "assumed"` and a `reason`; no number in it is presented as sourced, computed or measured, and
  `loss_tangent` is null. Applying it changes no Alpha-1 v1.1 number; it exercises the reader, the provenance checks and
  `materials apply --out` on a real design now, before the DFT bridge writes sidecars. Two `materials` checks (36 in the
  group, 511 in all) and the `cli_materials_alpha1_handwritten` guard.

### Fixed
- CLI: each command accepts only the options it reads. An unknown or misspelt option (`--shot 1000`, `-grid`, or another
  command's option such as `hwgate --hidden-dephasing-us`) used to be ignored without a word; it is now an error that
  names the option and the closest match, prints the command's usage and exits 1 before any file is read. So are an
  option given twice or without its value, a missing or extra positional argument, a malformed number (`--shots lots`;
  `--shots 1e5`, which `std::stoi` read as 1), a count out of range (`fab --samples 0` crashed, `run --traj 0` reported
  a fidelity of 1), a negative seed (`std::stoul` wrapped it) and an unknown choice (`--arch tunabel` gave a
  fixed-coupling chip; `--off contrl` and `--leakage-removal lrrr` did nothing). `new --preset` rejects the generator
  options it ignored, `new` takes one topology, `qec` rejects the `--burst-*` options without `--bursts` and a negative
  `--burst-rate` or `--burst-tau-ms`, and `budget` no longer accepts `--shots`, which it never used. For a misspelling
  that is a unique prefix of one option, the error names that option (`--burst-rad`: `--burst-radius`, not
  `--burst-rate`).
- CLI: `hwgate --out`, `cz --out`, the `synth-chip` data file and the `gds --cells` write-back fail with "cannot write"
  (exit 2) when the file cannot be opened or written, as `gds::write` does. `save_design` (`qpuforge new`) and
  `write_layout_svg` (`qpuforge layout`) throw likewise. All of them used to report "wrote" for a missing directory.
- Review and study programs (`alpha1_v11_review`, `beta1_review`, `gamma1_review`, `decoder_study`, `decoder_sweeps`,
  `surgery_study`): `provenance.source_commit` and `provenance.source_tree_dirty` now describe the tree the program was
  built from. They were compile definitions (`QPUFORGE_GIT_SHA`, `QPUFORGE_GIT_DIRTY`; no longer read) computed when
  CMake configured the build directory, so a directory configured at one commit and rebuilt after later commits kept
  writing the old commit into its outputs. The stamp is now taken on every build: the target `qpuforge_git_stamp_update`
  runs `cmake/git_stamp.cmake`, which rewrites the one source file of a small static library (`qpuforge_git_stamp`) only
  when the commit or the dirty flag changes. A new commit costs one small compile and six relinks, and a build directory
  does not have to be configured again before a recorded run. The Makefile rules `<program>/fast` skip a program's
  dependencies and therefore do not take the stamp. The dirty flag means what it did (tracked files modified, untracked
  files not counted) but is now read at every build: a build made after a recorded run has rewritten its tracked output
  files records `true`. Outside a git checkout or without git the stamp is `"unknown"` / `false`, as it is for a program
  compiled by hand (`examples/git_stamp.hpp`), and the build still works; a source tree unpacked inside another
  repository no longer takes that repository's commit. Two new CTest guards (label `unit`, reported as skipped where git
  is not installed): `git_stamp` runs the script on a scratch repository, and `git_stamp_wiring` configures a copy of
  `CMakeLists.txt` once, with a few lines in place of each program source, and builds it again after each commit. Every
  program that includes `git_stamp.hpp` must then print the new commit, and nothing but the stamp may be compiled.
  `fp_contract_flags` now requires the `QPUFORGE_COMPILE_OPTIONS` line for every `add_executable` and `add_library` in
  `CMakeLists.txt` (the target helper, `test_determinism` and the new library).

- `qpuforge qec --decoder uf-stream --herald` decoded with the batch union-find decoder in two places. The
  "decoded without heralds" comparison copied the options before the streaming decoder was attached. And
  `run_memory_frame` sent every shot that had a herald to `decode_scaled`, whatever `FrameRunOptions::stream` said, so
  the run itself was a mixture of the two decoders. Heralds now go through the streaming decoder
  (`StreamOptions::scaled`), and the comparison copies the run's own options. `cli_herald_stream` pins it: with heralds
  that never fire, the run and its comparison must agree, and must differ from `--decoder uf`.

### Changed
- **Qubit SQUID.** `xmon` cuts a 10 µm ground lead out of the qubit's (−x, −y) inner corner, 10 µm from the SQUID stub;
  the SQUID region bridges stub → lead. Before, it ended in etch on its ground side (0 of 17 qubit SQUIDs landed). A
  junction now counts as landed only if it joins an island to the ground plane.
- **Qubit flux line** shorted again into the (−x, −y) corner next to the SQUID (the corner is strapped to the plane).
- **XY line** ends open outside the qubit's cut-out, 247 µm right of and 279 µm below the qubit centre (it ran into the
  qubit's (+x, −y) corner and coupled 3.0 fF: a T1 limit of 0.2 µs at 5.5 GHz). Now 0.029–0.030 fF on AF-17.
- **Sliver fill** in `assemble_chip` is no longer guarded (`detail::guarded_fill` stays in the library): ground a
  fill cuts off is strapped back or removed. AF-17's four ground necks are gone.
- **Cells re-synthesised** (warm start) with the SQUID ground lead, and every mask, preview and review regenerated.
  AF-17: arms 123.4 / 125.7 µm (127.4 / 129.8 before); E_C within 0.12 %, couplings within 5.0 %, readout C −0.5 %.
  Beta-1 cells (also GF-97): E_C within 0.08 %, couplings within 5.4 %, readout C +0.3 %. Both re-null the parked
  exchange (AF-17 couplers 102 MHz below their sweet spot).
- **Masks**: DRC 0 + 0 on all three; metal connectivity clean on all three: 0 floating ground fragments (143 / 405 /
  821 before), every shorted line and every SQUID on the plane, with 122 / 373 / 774 ground straps (AF-17 / BF-49 /
  GF-97).
- **Gates.** The mask width/space DRC (`drc`) and metal connectivity (`connectivity`) are gated in the Alpha, Beta and
  Gamma reviews (ROADMAP.md). The control-line T1 budget (10× the design T1) is reported, not gated: AF-17 23.8× (XY) /
  22.4× (flux); the Beta-1 cells (BF-49, GF-97) 11.4× / **9.7×**. Alpha-1 v1.1: 7 of 7 criteria; Beta-1: 10 of 10;
  Gamma-1: 9 of 9.
- `qpuforge gds` prints the connectivity summary and exits 2 unless both the DRC and the connectivity are clean.
- `tools/gds_preview.py` draws the ground straps.
- **Layout cells.** The bus is never drawn narrower than its strip: below that, the direct-coupling knob of the
  synthesis pulls the bus ends back from the qubits instead of shrinking the T-bar to a sliver (0.075 µm on the 0.4
  AF-17 cells, 0.25 µm on BF-49). The synthesis aims the direct coupling 1.5 % above target so that the parked exchange
  can always be re-nulled by parking the coupler below its sweet spot.
- `stair(..., short_end)` trims the trench beside a shorted end, so the short meets the ground without the 1 µm etch slit
  0.4 left beside each qubit flux-line "short" (a DRC width violation at every qubit).
- **Mask DRC** (width + space violations), 0.4 → now: AF-17 34 + 202 → 0 + 0, BF-49 98 + 504 → 0 + 0,
  GF-97 194 + 1,016 → 0 + 0.
- **CZ calibration (CZ model version 4).** A pair whose plain `1 - F_avg` calibration fails the 2e-4 leakage limit is
  re-calibrated with the objective `1 - F_avg + 10 · leakage`; the repeat is kept only if it meets the limit, its
  1 − F_avg is no worse and it stays within 5e-4. 15 of the 16 pair classes keep their version-2 calibration bit for bit.
  5.81/5.29 GHz (Q8–Q11 on AF-17) goes from 2.33e-4 to 1.60e-5 leakage and from 4.06e-4 to 3.55e-4 coherent error; its
  gate is 3.5 ns longer and its total CZ error 0.5 % higher. Largest leakage over the 16 classes: 1.20e-4 (5.87/5.29).
  `results/cz_cache.json` regenerated; caches of other model versions are rejected.
- Regenerated on the merged tree with CZ model version 4 and the new cells (`results/PROVENANCE.md`): the three reviews
  at 100000 shots, `BETA-1.md` and `GAMMA-1.md`, the Alpha-1 v1.1 characterisation rehearsal, the three GDS files and
  their previews, and the test log. Alpha-1 v1.1: memory-Z 0.377 % ± 0.004 % per round, memory-X 0.360 %, Λ(3→5)
  projection 2.95 ± 0.18; the Alpha gate, now computed per pair, passes. Beta-1 and Gamma-1 pass every criterion (counts under Gates above).
- `enumerate_faults` is linear in the number of rounds (docs/design/linear-fault-enumeration.md): faults are propagated
  to the next measurement cut and finished with precomputed per-qubit tails. Output unchanged, fault by fault.
- The matching decoder's DP, blossom and greedy paths read distances through an oracle, so that `decode_scaled` can
  reuse them; results of `decode` are unchanged.
- **Behaviour change:** a malformed number, an unknown circuit (`ghz:x`) and an unknown `--basis`, `--decoder`,
  `--preset` or `--phase` value are usage errors (exit 1); they used to exit 2, numbers with "error: stoi". The usage
  text and README document the exit statuses: 0 done, 1 usage error, 2 error while running (also `gds`: layout not
  DRC- and connectivity-clean), 3 `hwgate`: gate not passed, `surgery`: a noiseless check failed.
- CLI: the usage text lists every option of every command (`budget --seed`, `qec --seed --cache --leakage-removal`,
  `synth-chip --cache`, `hwgate --seed --cache` and most `new` options were missing), and a usage error prints the
  usage of its command only. Options may come before the positional arguments.
- `build_surgery`: a face that ends at a split and whose whole support is read out there in its basis closes with a
  detector (the faces inside a routing strip), and a face that ends is forgotten, so that a face returning in a later
  phase is not compared with a measurement from before its data were read out. Neither case occurs in the schedules of
  PR #21 (adjacent patches): their circuits, detectors and results are unchanged.

### Breaking CLI changes

The strict option table (#24) and the follow-up for `--help` change what scripts can rely on. The list is limited to what
the entries under Fixed do not already state; unknown, misspelt, repeated and misplaced options, `budget --shots` and the
`qec` option dependencies are described there.

- `new --plan` accepts only 2, 4 or 8 (`--plan expects 2, 4 or 8, got '16'`, exit 1), as the usage line always said and
  `device.hpp` handles. Other values used to be accepted.
- `--help`, `-h` and `help` print the usage to stdout and exit 0, and `<command> --help` or `-h` (also `help <command>`)
  prints that command's usage only, exit 0. Before the strict table (#24), `--help` was an unknown command: the usage went to
  stdout with exit status 1. Since #24 it was an unknown option: a message and the usage on stderr, exit 1. Scripts that
  treated exit 1 as "help shown" must now expect 0. `--help` as the value of an option that takes one (`run --off --help`)
  is still that option's bad value, a usage error.
- The `cli_help`, `cli_help_short`, `cli_help_word` and `cli_help_command` guards (label `unit`) run `tests/check_cli_help.cmake`,
  which asserts exit status 0, an empty stderr and the usage text on stdout; a pass regex alone cannot see the exit status.

### API changes
- `surgery.hpp`: `prepare_surgery` has a sixth parameter `correlated` (default `false`), `run_surgery_frame` a fourth
  parameter `decoder` (default `SurgeryDecoder::Separate`), and `surgery_parity`, `surgery_cnot`, `surgery_experiment`
  and `surgery_cnot_reference` a last parameter `gap` (default 1). Calls written for the earlier signatures compile
  and give the same results. `SurgerySetup` and `SurgeryResult` have new members; `SurgerySetup` holds a
  `std::unique_ptr` more and stays move-only. New names in `qpuforge`: `SurgeryFault`, `JointEdge`, `JointDEM`,
  `SurgeryDecoder`, `build_joint_dem`, `surgery_faults`, `surgery_fault_signatures_of`, `surgery_decoder_name`,
  `run_surgery_decoders`; `surgery_fault_signatures` has an overload taking a fault list.
- `qec.hpp`: `MatchingDecoder` has a second constructor, `MatchingDecoder(n_det, edges)`, from any list of edges with
  members `a`, `b`, `obs` and `p`. `obs` may be a mask of up to 8 observables; `decode()` and `decode_scaled()` then
  return the XOR of the masks along the matching instead of 0 or 1. The constructor from a `DetectorErrorModel` is
  unchanged.
- `UnionFindDecoder::decode` is now a wrapper around a private template; its signature and results are unchanged.
  `StreamingUnionFind::decode(defects, ws)` is unchanged; `decode(defects, ws, StreamOptions)` is new. `StreamScratch`
  gained members. A `StreamingUnionFind` now also keeps the whole detector graph (`global()`).
- `FrameRunOptions::stream` together with heralds used in the decoder: the shots with a herald are decoded by the
  streaming decoder (they were decoded by the batch one). The streaming decoder must then be built from the setup's
  own detector error model; `run_memory_frame` throws otherwise.
- `unionfind.hpp` adds `UFOpCounts`, `UFCostModel`, `StreamOptions` and the `uf_detail` namespace to `qpuforge`;
  `qec.hpp` adds `HeraldModel`, `StreamProfile`, `fault_signatures`, `build_herald_model`, `prepare_herald_model`,
  `run_stream_profile`, and `frame_detail::LeakChain` / `leak_chain`.
- `BurstSpec` is now defined before `FrameRunOptions` (same members, plus `radius`, default 0).
- `frame_detail::FaultDesc` (internal) gained `owner`.
- `enumerate_faults(F, E, long* steps = nullptr)` and `optimize_cz(S, fluxc_on, t_rise, dt, max_eval, leak_weight = 0,
  leak_limit = 1, coh_limit = 1)` gain optional arguments. 0.4 calls compile unchanged; `optimize_cz` called without
  them has no leakage fallback (`calibrate_cz_pair` passes the `CZModel` values).
- `CZModel::kVersion` is 4 (was 2): a CZ cache of another model version, such as one written by 0.4, is ignored with a
  warning and its calibrations are re-simulated.
- New header `include/qpuforge/materials.hpp`, included by `qpuforge.hpp`. It adds `MaterialsQuantity`,
  `MaterialsProducer`, `MaterialsSidecar`, `MaterialsApplyOptions`, `MaterialsReportLine`, `MaterialsReport`,
  `materials_from_json`, `load_materials`, `apply_materials` and the `materials_detail` namespace to `qpuforge`. Code
  that includes `qpuforge.hpp`, defines one of these names itself and uses `using namespace qpuforge` can now hit an
  ambiguity.
- New headers `include/qpuforge/unionfind.hpp`, included by `qec.hpp`, and `include/qpuforge/surgery.hpp`, included by
  `qpuforge.hpp`. They add `UFScratch`, `UnionFindDecoder`, `StreamScratch`, `StreamingUnionFind` and the
  lattice-surgery names (`Surgery*`, `StabilizerTableau`, `make_surgery_lattice`, `region_faces`, `build_surgery`,
  `surgery_*`, `tableau_records`, `check_surgery_noiseless`, `prepare_surgery`, `run_surgery_frame` and the
  `surgery_detail` namespace) to `qpuforge`, with the same possible ambiguity; for the union-find names, also in code
  that includes only `qec.hpp`.
- `save_design` and `write_layout_svg` (`device.hpp`) throw `std::runtime_error` when the file cannot be opened or
  written; they used to return as if it had been written.
- The materials reader changes nothing else: chip JSON files and the output of `qpuforge analyze` without
  `--materials` are as before.
- Stricter sidecar validation (D4): a hand-written sidecar whose numbers have no `provenance`, an unknown one, or a
  `literature` / `measured` value without `reference` is now rejected (`std::runtime_error`), as is a value outside the
  ranges above, a sidecar without `produced_by.tool`, or inconsistent `placeholders` / `contains_placeholders`.

### Known issues
- **Open review findings on the mask connectivity now on `main`** (details in docs/design/mask-connectivity.md). The
  approach is known to be electrically wrong on the points below; a redesign is in progress in PR #28:
  - the SQUID ground lead and the flux-line short share one modelled airbridge to the plane (galvanic T1 channel,
    estimated 0.1–1 µs; the pocket is treated as ideal ground);
  - the flux-line short couples only ~0.004–0.011 pH to the SQUID loop (~90–290 mA for Φ0/2): the flux lines do not
    yet tune the qubits;
  - the XY bump pad is ~30 µm from two coupler pads (C ≈ 1.4–1.6 fF against 0.036 fF to the qubit): the 0.03 fF /
    23.8× control-line figures do not hold;
  - the connectivity gate counts single-bridge chains as ideal ground; the junction rule and the control-line model
    lack mutation-resistant tests; `line_budget` overstates the flux T1 limit by (f01/5.5 GHz)².
- **Control-line budget missed on BF-49 and GF-97.** The flux line's T1 limit is 9.7× the design T1 of the worst
  qubit (budget 10×), so the budget is reported, not gated. The XY line meets it (11.4×).
- **No foundry DRC deck has been run.** The masks pass the in-house checks (placement, 4 µm width/space on the merged
  etch layer, metal connectivity, control-line budget). Corner-to-corner spacing and airbridge rules are not checked.
- **Ground straps are modelled, not extracted.** A strap's capacitance to the coupler body or bus strip it crosses
  (about 0.2 fF over a 6 µm bus at 3 µm bridge height) is not in the extraction cells.
- **The flux lines' mutual inductance to the SQUID loops** (the flux coupling itself and its own T1 channel) is not
  modelled; the control-line budget covers their capacitive coupling only.

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
  `designs_v04`, and the serialisation groups `json`, `gds`, `json_rfc`, `json_depth`, `gds_write`), checks of the
  opt-in flux term and the 4 MHz direct coupling in `coupler_and_flux`, and the tuned-anharmonicity |11⟩–|20⟩ test in
  `transmon`: 321 checks in all (186 in 0.3.0), and a third review smoke test. Outside the validation checks: the
  `determinism` (1, 2, 3, 4 and 8 threads in one process) and `determinism_env` (separate processes under
  `OMP_NUM_THREADS` = 1, 2, 4, 8) guards, the `fp_contract` and `fp_contract_flags` build-flag guards, and the
  `cli_version`, `cli_version_short` and `cli_usage` CLI guards.

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
- Every review output except the historical Alpha-1 v1.0 review (`results/alpha1_review.*`: pre-PR #2 model, not
  regenerated), and `results/design_study.*` and `results/test_output.txt`, were regenerated from `802ceb7` with CZ
  model version 2, the exact decoder and the thread-count-independent Monte Carlo (`results/PROVENANCE.md`).
- Build: every target is compiled with `-ffp-contract=off` (`QPUFORGE_COMPILE_OPTIONS`), so `-march` / FMA choices no
  longer change results; `fp_contract` and `fp_contract_flags` check it.
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

## [0.3.0] - 2026-11-27


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
  `docs/af17-*.png` (an earlier drawing that does not match `AF17_v11.gds`).
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
