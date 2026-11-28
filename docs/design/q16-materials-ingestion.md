# Design note: the atomforge → qpuforge materials path (q16)

Status: **assessment and decided plan.** Phase 1 is implemented on `feat/materials-flag-remaining`. Decisions D1–D5 were
**decided by Finnian Bond, 2026-11-28 (relayed by Forge PM)**; section 6 records them and section 5 is now the chosen
build order. Phases 2–5 are planned, not yet merged. Assessed against qpuforge `main` at `5495bce` (PRs #20, #22, #25, #26 merged) and atomforge
`main` at `b6f6d35` (PR #17 merged), 28 November 2026.

## 1. Summary

The *format* path is closed. A sidecar written by atomforge's own `search --materials` is read, checked and applied by
qpuforge without any manual step (section 3). The *value* path is empty: atomforge can only write placeholders (value
`null`), so nothing a real run produces changes a qpuforge number today. Closing that needs a producer of real
quantities (the DFT bridge, a07/a08, and a loss model), not more qpuforge reader code.

What was still missing on the qpuforge side, and what this design covers:

| # | Gap | Where | Phase |
|---|---|---|---|
| G1 | `--materials` was read by `analyze` only; `run`, `budget`, `fab`, `qec`, `cz`, `synth-chip`, `hwgate` ignored it **silently** (checked: `fab chip --materials s.json` ran as if the flag were absent). The values never reached the emulators or the logical-error runs | CLI | **1 (done)** |
| G2 | Applying a `junction_gap_ueV` moves every qubit frequency (a chip stores `Rn_kOhm`, not `E_J`); the report line did not say so | reader report | **1 (done: note on the line)**; D1 decided: keep fixed-Rn |
| G3 | A sidecar can only be applied in memory; there is no way to save the result, so a reviewed number cannot be traced to a sidecar | CLI / chip JSON | 2 |
| G4 | No cross-repo test: the qpuforge fixture is a byte copy of atomforge output at `ee734f9`; atomforge's CI never runs qpuforge's reader | CI | 3 |
| G5 | Stale documentation: atomforge's README links the contract at branch commit `f64238f`; qpuforge's `DATA_CONTRACT.md` quotes atomforge `ee734f9` | docs | 3 |
| G6 | Validation holes: `provenance` is free text (an EAM-computed sidecar could carry a number labelled anything); `conditions` (T, f) are ignored; `system` is not compared with anything; `contains_placeholders` is not cross-checked | reader | 4 |
| G7 | Only one chip-wide `Q_dielectric`; no per-element or participation-weighted loss. This is the contract's open question 1 and needs the FEM work (q17) | schema v1 | 5 |

## 2. What PR #22 gives qpuforge

| Aspect | State on `main` |
|---|---|
| Schema and version | `forge.materials.v0` required exactly; any other schema, a missing schema or a non-string schema is an error that names what was found. Additive members stay on v0 and are ignored: `structures_ref`, `atomforge_computed`, `not_in_contract`, `open_stack`, `notes`, and a quantity's `reason`, `relation`, `conditions`. An unknown *quantity* is an error |
| Quantities | `Q_dielectric` (`dimensionless`) → `ChipDesign::Q_dielectric`; `junction_gap_ueV` (`ueV`) → `ChipDesign::gap_ueV`; `loss_tangent` (`dimensionless`) → `Q = 1/tanδ`, opt-in (`--loss-tangent`), only when `Q_dielectric` is null/absent |
| Units | Exact unit strings (`dimensionless`, `ueV`; not `µeV`); no conversion, a different unit is an error. `T_eff_mK`, `Rn_spread_pct` and the other environment fields are rejected under `quantities` |
| Values | null (skipped, reported) or finite and > 0; a placeholder carrying a number is rejected; `qpuforge_default_for_reference` is never read; `maps_to` must match the contract |
| Errors | `std::runtime_error` before the design is touched, prefixed with the file name |
| Report | `MaterialsReport`, one line per quantity (applied/skipped, before/after, why) |
| CLI | `analyze chip.json --materials sidecar.json [--loss-tangent]`; `--loss-tangent` alone is an error |
| Tests | `materials` group (25 checks); five CLI guards; fixtures in `tests/data/`, one of them atomforge's sidecar byte for byte |
| Docs | `DATA_CONTRACT.md` "Implementation status", README, CHANGELOG |
| CI | the materials tests run in the normal matrix (gcc, clang, `-Werror`, sanitizers) because they are in `ctest -LE smoke` |

### Consumers of the two fields

| Field | Used in | Reaches |
|---|---|---|
| `Q_dielectric` | `analyze()`: `G_diel = 2π f01 / Q` → `T1`, `T1_diel`, `T2` (`device.hpp`, one line) | `noise_model(D, R)` → `run`, `budget`, `qec`, `cz` (decoherence part), `synth-chip`, `hwgate`; `fab_monte_carlo` (via `analyze`) |
| `gap_ueV` | `EJ_from_Rn_GHz(Rn_kOhm, gap)` in `analyze()` (f01, anharmonicity, collisions) and in the pair analysis of tunable qubits (CZ interaction point, `device.hpp`); `Rn_from_EJ_kOhm` at *design time* in `inverse_design()` (so `qpuforge new` and presets bake the default 180 µeV into every stored `Rn_kOhm`) | everything downstream of `analyze()` |

`Rn_spread_pct` stays a chip-JSON/process input, as the contract says.

## 3. End-to-end check with a real atomforge file

atomforge `main` `b6f6d35` (gcc, `-Werror`, 0 warnings, 32/32 ctest) ran `atomforge search Cu Au --quick --materials`;
qpuforge `main` `7380b58` (before #20; the reader and `analyze` are unchanged by it) read the result:

* The file has the schema, `produced_by` (commit `b6f6d35…`, `Zhou04-EAM`, ISO-8601 with offset), `system`,
  `structures_ref` (`report.json`, relative), three placeholder quantities (value `null`, `status` and `provenance`
  `placeholder`), `placeholders`, `atomforge_computed`, `not_in_contract`.
* `qpuforge analyze alpha1_v11.json --materials materials-contract-CuAu.json` prints the provenance line and three
  `skipped … placeholder, value null` lines, exit 0; the device report below it is byte-identical to the report without
  the flag.
* Not a bug but the point of G1/G7: nothing was applied, because nothing real was written.
* Before phase 1, `fab … --materials …` (and `run`, `qec`, …) ran without a word of complaint and without the file being
  read, so even a non-placeholder sidecar would have been dropped.

## 4. Remaining design

### Phase 1 (implemented): one loader for every design-consuming command

`load_design_cli(argc, argv, path)` loads the chip and applies `--materials [--loss-tangent]` in memory, printing the
per-quantity lines first. `analyze`, `run`, `budget`, `fab`, `qec`, `cz`, `synth-chip` and `hwgate` use it. `layout` and
`gds` refuse the flag with an error: both draw from the chip JSON, and `gds` synthesises cells to frequency targets that a
changed gap would move; silently ignoring was the failure, so they say no. Without `--materials` nothing changes: no
committed output moves. The `junction_gap_ueV` report line now says the frequencies scale with the gap.

### Phase 2: save an applied sidecar (G3), decided (D2: build it)

`qpuforge materials apply chip.json sidecar.json --out chip2.json [--loss-tangent]` writes the chip with the accepted values
and a `materials_source` block (`tool`, `commit`, `schema`, SHA-256 of the sidecar, per-quantity provenance) so a reviewed
chip names its origin. Alternatives:

| Option | For | Against |
|---|---|---|
| A. Keep in-memory only (today) | no new output format | not reproducible from the chip file; each command needs the sidecar path |
| B. `materials apply --out` (**chosen**, D2) | explicit, diffable, one new file; chips stay valid JSON for old readers if `materials_source` is ignored | `design_from_json` must accept or skip a new block; touches `device.hpp` |
| C. Writer in atomforge that emits a chip JSON | one step for the user | couples the repos (atomforge would need `ChipDesign`) — rejected |
| D. Library dependency (FetchContent) | one code path | contract says the repos are companions, not a code dependency — rejected |

### Phase 3: cross-repo test and doc sync (G4, G5), decided (D3: add the atomforge CI job)

* qpuforge: a `tools/refresh_materials_fixture.sh` that regenerates `tests/data/materials-contract-CuAu.json` from a given
  atomforge commit and records the commit in the fixture's test header; the existing byte-for-byte check then catches an
  atomforge writer change when the fixture is refreshed.
* atomforge CI job (non-blocking at first, `continue-on-error`; qpuforge ref is a workflow input defaulting to `main`, so a
  broken qpuforge `main` cannot block an atomforge merge): check out qpuforge, build `qpuforge`, generate the synthetic sidecar the
  test already builds, run `qpuforge analyze … --materials`, require exit 0 and the three `skipped` lines. This is the only
  test that fails when either side changes alone.
* atomforge README link: point at qpuforge `docs/DATA_CONTRACT.md` on `main` instead of commit `f64238f` (an atomforge PR).
* Alternative considered: a shared JSON Schema file in both repos. Cheaper to check, but a second source of truth next to
  the reader; rejected for now.

### Phase 4: tighter validation (G6), decided (D4: enumerate provenance now, in v0); implemented on `feat/materials-provenance-enum`

* `provenance` is an enumerated set on `forge.materials.v0` itself (no v0.1, no v1 wait): the four values of the contract's
  §3.2 example (`assumed`, `measured`, `dft-loss-tangent`, `participation-weighted`), `placeholder` (what atomforge
  writes), and `literature` (contract open question 5 names a literature value for the gap). `assumed` is the label for a
  hand estimate. Any other string is an error; a number without a provenance is an error; `literature` and `measured`
  name their source in a `reference` string. atomforge writes only `placeholder`, so no writer breaks; hand-written
  sidecars without a provenance on a number do.
* Values outside a physical range are rejected, not only non-positive ones (bounds in `DATA_CONTRACT.md`).
* Warn (not error) when `conditions.f_GHz` differs from the qubit's `f01` by more than a tolerance, or `T_K` from
  `T_eff_mK`, because `Q = 1/tanδ` is only meaningful under the uniform-loss model.
* `contains_placeholders` must agree with the `placeholders` list and the per-quantity status.
* The reader keeps ignoring `system` for application (a Cu–Au file is harmless), but prints it, as now.

### Phase 5: real producers and per-element loss (G7), depends on other work

* D5 decided **both**: first a hand-written Alpha-1 sidecar (`provenance: "assumed"`, every number explained in its `reason`
  and none presented as a sourced value), so the whole path is exercised now; later the DFT bridge writes sidecars.
* atomforge side: a07/a08 (DFT `ExternalPotential` bridge) gives energies in eV; it still gives no loss tangent. A
  tanδ needs a separate model (TLS density × dipole moment, or measured), so the first non-null quantity will more likely be
  `provenance: "assumed"` / `"measured"` from a person than `dft-loss-tangent` from the code.
* qpuforge side: participation ratios need the FEM cross-check (q17). Per-element Q means a `ChipDesign` change
  (per-qubit/per-resonator Q), hence **forge.materials.v1** and a version policy: the reader keeps v0 and reads v1
  explicitly; a v1 file is never silently down-converted.

## 5. Build order (decided 2026-11-28)

Each item is a new branch stacked on the one before it (qpuforge unless noted); none is merged without Finnian.

| # | Branch | Repo | Content | Risk |
|---|---|---|---|---|
| 1 | `feat/materials-flag-remaining` | qpuforge | `--materials` on every design-consuming command, `layout`/`gds` refuse; gap note; CLI tests; docs (done) | low: no output changes without the flag |
| 2 | `feat/materials-provenance-enum` | qpuforge | provenance enum in v0, `reference` for literature/measured, range checks, `placeholders`/`contains_placeholders` consistency (D4) | low–medium: hand-written sidecars without provenance now fail |
| 3 | `feat/materials-apply-out` | qpuforge | `materials apply --out` + `materials_source` block, fixed-Rn kept (D1, D2) | medium: `design_from_json` / `design_to_json` |
| 4 | `data/alpha1-sidecar-handwritten` | qpuforge | hand-written Alpha-1 sidecar, `assumed`, load-and-validate test (D5, first half) | low |
| 5 | `ci/qpuforge-reader-crosscheck` | atomforge | non-blocking CI job: build qpuforge at a pinned ref, read the sidecar atomforge emits (D3) | low (CI only) |
| 6 | (later) | qpuforge | q18 opt-in flux-noise term, then the q17 FEM reference | see the q17/q18 note |
| 7 | (later) | both | DFT bridge writes sidecars (D5, second half); v1 schema and per-element Q after q17 and a07/a08 | high; separate design note |

Not scheduled yet: the `conditions` warnings and the fixture refresh script (phase 3a).

## 6. Decisions (decided by Finnian Bond, 2026-11-28, relayed by Forge PM)

| # | Question | Decision |
|---|---|---|
| D1 | Gap semantics: fixed `Rn_kOhm` (f01 moves) or retune `Rn_kOhm` (f01 fixed)? | **Keep fixed-Rn.** A sidecar gap changes `E_J` at fixed `Rn_kOhm`; the report line says the frequencies move. `materials apply --out` keeps every `Rn_kOhm` too |
| D2 | Persist an applied sidecar into a chip JSON? | **Yes:** build `materials apply --out` (phase 2) |
| D3 | atomforge CI job that runs qpuforge's reader? | **Yes** (phase 3), non-blocking, qpuforge ref pinned by a parameter |
| D4 | Enumerate `provenance` in v0 or only in v1? | **Now, in v0** (phase 4) |
| D5 | First real quantity: hand-written sidecar or DFT bridge? | **Both:** a hand-written Alpha-1 sidecar first, the DFT bridge writing sidecars later |

The options as they were put (kept for the record):

* **D1, gap semantics.** A sidecar gap currently changes `E_J` at fixed `Rn_kOhm`, so every f01 moves (Alpha-1 v1.1, 180 → 190
  µeV: qubit 0 5.81 → 5.975 GHz, the design leaves its frequency plan). The alternative is to keep f01 fixed and retune
  `Rn_kOhm`. Recommendation: keep the physical fixed-Rn behaviour, which says what a different gap would do to the same junction,
  and surface it (done); a retune would hide a real effect.
* **D2, persistence.** Is saving an applied sidecar into a chip JSON (phase 2) wanted, or should the sidecar stay an
  in-memory overlay?
* **D3, cross-repo CI.** Run qpuforge's reader from an atomforge CI job (phase 3), accepting that atomforge CI then depends on
  qpuforge `main` being buildable?
* **D4, strictness.** Make `provenance` enumerated in v0 (stricter than the contract text), or only in v1?
* **D5, first real quantity.** Who supplies the first non-null value: a hand-written `assumed`/`measured` sidecar for Alpha-1
  (cheap, exercises the whole path now), or wait for the DFT bridge?
