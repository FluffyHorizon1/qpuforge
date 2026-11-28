# Feasibility note: FEM cross-check (q17) and flux-pulse distortion (q18)

Status: **decided, implementation planned.** Decided by Finnian Bond, 2026-11-28 (relayed by Forge PM); see
[Decisions](#decisions-finnian-bond-2026-11-28) and [Implementation plan](#implementation-plan). Scope note for issues q17 (backlog19) and q18 (backlog20), i.e. the two Roadmap tooling items
in the README. Written against `main` at 7380b58 (0.4.0). Function names are used instead of line numbers.

**Recommendation.** Do both as small PRs that leave every committed number unchanged until Finnian chooses to regenerate.
q17: fix a solver-agnostic exchange format, keep the solver out of CI, and validate the FEM setup against the MoM under
*matching* assumptions before comparing against the realistic stack. q18: start with a sensitivity study that needs no new
physics code. It shows that sub-per-mille flux error is already the size of the 1e-3 residual allowance.

## q17: FEM cross-check of the cells

**What it checks and where it hooks in.** `extract_capacitance` is a MoM with zero-thickness metal on an infinitely thick
substrate, no airbridges, no interposer. It feeds `star_cell` / `pair_cell` → `triple_capacitance` / `triple_circuit` →
E_C, g_qc, g12, coupler E_C, and from there `synthesize_cells`, `evaluate_cells` and `post_layout` (α, g, residual
|11>-|20> gap, ZZ at park). The check compares the cells' Maxwell matrices and the derived E_C and g. The second half of the
issue, participation ratios, hooks into `analyze()`, which uses one chip-wide `G_diel = 2πf / Q_dielectric`. The
`forge.materials.v0` reader (`materials.hpp`) carries only `Q_dielectric` or `loss_tangent`.

**Options.**
- *External solver via files.* Elmer (`Calculate Capacitance Matrix` writes `CapacitanceMatrix.dat` in the lumped
  convention, so flip the off-diagonals) or Palace (`Electrostatic` writes the Maxwell matrix to `terminal-C.csv`). It gives an
  independent number with finite substrate and metal thickness. Palace's eigenmode
  `Boundaries/Postprocessing/Dielectric` (MA, MS, SA types) also writes interface participation and Q to `surface-Q.csv`,
  which is what the TLS item needs. FEniCS means owning the weak form and meshing for no gain.
- *Analytic / approximate, in tree.* More closed forms (the conformal-map CPW exists), a `MeshOptions` convergence sweep, a
  finite-substrate image-series Green's function. This is cheap and guards regressions, but it uses the same engine, so it
  is a sensitivity check, not an independent one.

Use both: the analytic side as the guard. **Decided: support both Elmer and Palace from the start** for the
capacitance cross-check (one exchange format, an importer per solver's native output); Palace also carries the later
participation work.

**Keeping CI free of heavy dependencies.**
1. Exchange format = JSON like the materials sidecar: `qpuforge.fem_crosscheck.v0` with conductor names, Maxwell C in fF,
   solver and version, mesh and domain parameters, cell hash. A C++ reader and `compare_cells()` report deltas in C, E_C, g, α.
2. The reference is produced offline and committed. An ordinary CTest compares the MoM with it at a stated tolerance.
   No `find_package`, no solver at build time.
   **Decided:** the FEM-vs-model tolerance is *reported* first (per cell and per quantity, in the test log and the
   review JSON) and not enforced; a pass/fail tolerance is fixed later, once stage-one agreement has been measured.
3. The driver (`tools/fem_crosscheck.py`: conductor rectangles → gmsh → solver → JSON) runs only if the solver is present
   and otherwise returns CTest's skip code. A container job may run on `workflow_dispatch` or nightly, never as a required check.

**Effort (PR-sized).**
1. Format, reader, `compare_cells()`, synthetic fixture, `docs/DATA_CONTRACT.md` section. C++ only. Small.
2. Geometry export and driver, validated first on known-answer cases (CPW against conformal mapping, a simple coplanar
   capacitor) with matching assumptions, so a setup error cannot pass as physics. Medium.
3. Run the three cells in two stages (matching assumptions, then the realistic stack). Commit the reference and a
   `PROVENANCE.md` row, and replace the README "before tape-out" sentence with measured deltas. Medium; someone must run the solver.
4. Participation ratios: Palace EPR → per-qubit p_MA, p_MS, p_SA, an optional `participation` block, `1/Q_i = Σ p·tanδ`,
   falling back to `Q_dielectric` when absent. Needs per-interface tanδ in the sidecar, which is an atomforge contract
   change. Separate and later.

**Risks.** 1 µm edge features on ~1 mm cells give large meshes. Domain size can dominate the difference, so a
domain-convergence run is needed. The MoM matrix includes the ground conductor and FEM uses a boundary, so the reduction must be
written down. C_Σ errors map ~1:1 onto E_C (1 % ≈ 2 MHz of α at E_C ≈ 0.22 GHz). Participation values depend on the assumed
interface thickness and tanδ, so they are only as good as the materials inputs.

**First step (decided):** PR 1 with importers for both Elmer (`CapacitanceMatrix.dat`) and Palace (`terminal-C.csv`),
deviation reported, not enforced; then PR 2 with the known-answer validation as acceptance. No claim about the real stack
until stage-one agreement is shown.

## q18: flux-pulse distortion and predistortion

**What it models and where it hooks in.** The flux line is a linear system between the AWG and the SQUID: finite bandwidth,
ringing, exponential tails from ns to µs. The CZ is `CZPulse` (sin² ramps, flat top) run by `simulate_cz`, whose `step` moves
qubit-1 and coupler flux with one shape `s(t)`. Distortion acts on the flux waveform *before* the flux-to-frequency map
(`FluxTransmon::f01`), so it belongs inside `simulate_cz`, not as a frequency correction.

**Error budget.** For calibrated chips, `calibrate_chip_cz` sets `err_ctrl = (err_coherent − leakage) + err_flux_local +
err_flux_quasistatic + residual_ctrl_err` and writes `err_2q_ctrl`, `leak_2q` and `err_2q`. The QEC and circuit emulators read
`err_2q_ctrl` and `leak_2q`, **not** `analyze()`'s `e_2q` (`err_2q_control` plus the optional `cz_flux_err`). A distortion term
must therefore be added in `calibrate_cz_pair` / `calibrate_chip_cz`, as a separate `err_distortion` plus its leakage share,
to reach any logical-error rate or fidelity. A proxy in `analyze()` alone would change the reported `e_2q` and nothing
downstream (as for q03's flux term).

**Sensitivity from the committed data** (`results/cz_cache.json`: Alpha-1 v1.1 coupler, 16 unique pairs, t_cz 36–49 ns).
The stored `curv_flux` is already the curvature of 1−F in flux (~1e5 per Φ0²). Take a constant fractional offset ε of both flux
steps over the whole gate. This is an upper bound: hardware calibration re-tunes amplitudes, and tails decay. Then
infidelity ≈ ½·curvature·(εΔΦ)²:

| ε | median over 16 pairs | worst pair |
|---|---|---|
| 1e-4 | 1.4e-5 – 4.8e-5 | 9.8e-5 |
| 3e-4 | 1.3e-4 – 4.3e-4 | 8.8e-4 |
| 1e-3 | 1.4e-3 – 4.8e-3 | 9.8e-3 |

(The range comes from how the stored curvature, which is the mover's plus the coupler's, splits between the two lines. The worst-pair column uses the upper end.)
A residual offset of a few 1e-4 of the step therefore uses the whole 1e-3 allowance. A model turns that allowance into a
requirement on predistortion accuracy.

**Options.**
- *Line model:* `FluxLine` per line (mover, coupler) = bandwidth limit plus Σ aₖ·exp(−t/τₖ) tails, as an exact-exponential IIR.
- *Predistortion:* inverse of the fitted tails (IIR/FIR), limited by AWG sample rate, bit depth and FIR length. The residual is
  set by fit error in aₖ, τₖ, which is the quantity to sweep; a perfect inverse would be an unrealistic zero.
- *Effects on the CZ:* ramp distortion changes adiabaticity (|11>→|20> leakage, spectator LZ transfer). Flat-top tails detune
  the interaction point during the hold (conditional phase). Coupler-line distortion changes the coupling strength.
  Single-qubit Z phases are absorbed by the existing free-Z fidelity metric.
- *Fast path* (no `simulate_cz` change): effective offset = tail averaged over the hold, then the stored curvature (slow tails
  only). *Full path:* a waveform-driven `simulate_cz` with the flat top sampled (~0.5 ns). Measured here: 1.6 ms per propagation
  step (dim 64), 141 steps per call, so sampling a 30–60 ns flat top costs roughly 1.4–1.9× per call, before the 260-evaluation
  optimiser. An ad hoc 5.0/4.6 GHz pair took ~90 s to calibrate, so measure real pairs first.

**Effort (PR-sized).**
1. The sensitivity table above as a test, from existing calibrations. No model code. Small.
2. `FluxLine` (IIR, inverse filter), closed-form tests, not wired in. Small.
3. Waveform-driven CZ path; an identity filter must reproduce `simulate_cz`. Medium.
4. `err_distortion` in `CZCalibration`, a design-JSON `flux_line` block, default off, its own budget line. **Decided:**
   additive on top of `residual_ctrl_err` (not carved out of it), so Alpha-1 v1.1 and Beta-1 results are not regenerated. Extend the cz cache key
   only when enabled so existing entries and outputs stay valid. Medium.
5. Predistortion study: fit-error and AWG sweeps, one review table, the accuracy the 1e-3 / 5e-4 budgets need. Medium.
6. Later: measured line response in `qpuforge.characterisation.v0`; gate-to-gate tail history (circuit level).

Multi-qubit spectator pulse simulation, the other half of the roadmap item, is separate and larger. **Decided:** it gets its
own issue ("Spectator multi-qubit simulation"), filed by Forge PM.

**Risks.** No measured line response exists, so every number is an assumption until hardware, and a precise-looking table
would mislead. `residual_ctrl_err` already covers distortion, drift, TLS and crosstalk, so a new term can double-count. The path
stays deterministic (no Monte Carlo). Calibration runtime grows, and `CZModel::kVersion` must change if default behaviour does.

**First step (decided):** item 4 in its fast form (effective offset from the stored curvature, so no `simulate_cz`
change), together with item 1's sensitivity table as its test; then items 2–3 if the fast path shows the full path is
worth its cost.

## Decisions (Finnian Bond, 2026-11-28)

Decided by Finnian Bond, 2026-11-28 (relayed by Forge PM). These replace the open questions of the proposal.

1. **q17 solver:** support **both Elmer and Palace** from the start. Both write to the one `qpuforge.fem_crosscheck.v0`
   format through an importer per solver; neither is a build dependency.
2. **q17 tolerance:** stage-one agreement under matching assumptions remains the bar before any claim about the real
   stack. The FEM-vs-model deviation on C, E_C and g is **reported first** and tightened into an enforced tolerance later.
3. **q18 term:** **additive** on top of `residual_ctrl_err`, **off by default**. No carve-out, so no committed result is
   regenerated; turning it on is a sensitivity study that double-counts by design.
4. **Spectators:** multi-qubit spectator simulation gets **its own issue**.

## Implementation plan

Each item is a branch off `main` (or stacked where noted), opened as a PR by Forge PM after QA's CI. Every branch
builds with 0 warnings under `-Werror` and passes `ctest -LE smoke` on g++ and clang. Nothing changes a committed number.

| # | Branch | Scope | Acceptance |
|---|---|---|---|
| 1 | `feat/q18-flux-term-optin` (off `main`) | `FluxLineSpec` (`tunable_coupler.flux_line`: static offset + exponential tails), `flux_line_eps_eff`, `flux_distortion_err`, `CZCalibration::err_distortion` added to `err_ctrl` → `err_2q_ctrl` after `residual_ctrl_err`; `qpuforge cz` reports it when on; `docs/design/flux-line-distortion.md` | off: calibration JSON and `err_ctrl` bit-identical, no key written, cache key and `CZModel::kVersion` unchanged; on: additive, reaches the emulators' per-pair error; the sensitivity table above reproduced from `results/cz_cache.json` |
| 2 | `feat/q17-fem-reference` (stacked on 1) | `fem.hpp`: `qpuforge.fem_crosscheck.v0` reader, importers for Elmer `CapacitanceMatrix.dat` (lumped → Maxwell) and Palace `terminal-C.csv`, `compare_cells()` reporting ΔC, ΔE_C, Δg per cell; small fixtures in `tests/data/fem/`; `DATA_CONTRACT.md` section | both importers give the same Maxwell matrix for the same fixture; deviations reported against the stated tolerance, never failing on it |
| 3 | (later) known-answer driver | `tools/fem_crosscheck.py` (gmsh → Elmer or Palace → JSON), CTest skip code without a solver; CPW and coplanar-capacitor known answers | stage-one agreement under matching assumptions |
| 4 | (later) committed references | the three cells, both stages, `PROVENANCE.md` rows, README measured deltas; then an enforced tolerance | Finnian's sign-off on the tolerance |
| 5 | (later) q18 full path | `FluxLine` IIR + inverse filter, waveform-driven `simulate_cz` (identity filter reproduces it), predistortion study | identity filter bit-for-bit; one review table |
| — | separate issue | spectator multi-qubit simulation | see the issue |
