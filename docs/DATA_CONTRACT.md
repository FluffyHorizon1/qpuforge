# Forge materials data contract (atomforge → qpuforge)

**Status:** Proposal for review. This document only; it is not an implemented API and this change adds no code.  
**Drafted:** Wednesday 25 November 2026, 22:40 Europe/London (UTC+0)  
**Filed:** Thursday 26 November 2026, Europe/London (UTC+0), as `docs/DATA_CONTRACT.md` on branch `docs/materials-data-contract` (based on qpuforge `main`).  
**Grounding:** atomforge `main` (`ac56e4b`); qpuforge `release/v0.3` (`cee64fa`). Citations to qpuforge sources are to `release/v0.3` (open PR #1). qpuforge `main` (`600e908`) does not yet contain the v0.3 device model or `materials` JSON keys described here.

The two repositories are programme companions, not a code dependency. There is no CMake `add_subdirectory`, FetchContent, or shared source. This note describes what each side already writes or reads, then proposes a versioned interchange file that maps onto those existing fields. It does not invent a C++ loader, a new qpuforge constructor, or a private atomforge RPC.

---

## Implementation status

*Added Saturday 28 November 2026 on branch `scs-bot/materials-reader` (qpuforge 0.4 integration, `25139c7`) and merged into `main` in PR #22. The writer it reads is atomforge's `include/atomforge/materials.hpp` (`ee734f9`). Sections 1–5 are the proposal as filed; where they say that no reader exists, this section replaces them.*

qpuforge now reads the sidecar. The reader is `include/qpuforge/materials.hpp`, included by `qpuforge.hpp`; `device.hpp` is unchanged (§3.3). `load_materials()` / `materials_from_json()` read and check a file, and `apply_materials()` writes the accepted values into a `ChipDesign` and returns one report line per quantity: applied or skipped, and why. `qpuforge analyze chip.json --materials sidecar.json [--loss-tangent]` applies a sidecar to the design in memory, prints those lines, then prints the device report. The chip JSON file is not changed.

| Quantity | Unit string | Written to | Accepted value |
|---|---|---|---|
| `Q_dielectric` | `dimensionless` | `ChipDesign::Q_dielectric` (`materials.Q_dielectric`) | null, or finite and > 0 |
| `junction_gap_ueV` | `ueV` | `ChipDesign::gap_ueV` (`materials.junction_gap_ueV`) | null, or finite and > 0 |
| `loss_tangent` | `dimensionless` | `ChipDesign::Q_dielectric` = 1 / `loss_tangent` | null, or finite and > 0 |

- `loss_tangent` is used only on request (`--loss-tangent`, `MaterialsApplyOptions::use_loss_tangent`) and only when `Q_dielectric` is null or absent. The conversion is the uniform-loss relation of §2.2, and the report line says so. A `Q_dielectric` with a value always wins.
- A null value changes nothing and is reported as skipped. This covers atomforge's placeholders (value null, `status` and `provenance` "placeholder", listed in `placeholders`). A placeholder that carries a number is an error, and `qpuforge_default_for_reference` is never read. atomforge's placeholder sidecar leaves Alpha-1 v1.1 unchanged.
- Errors (`std::runtime_error`, raised before the design is touched): a schema other than `forge.materials.v0`; an unknown quantity; `Rn_spread_pct` or an environment field under `quantities` (§2.1: they stay in the chip JSON); a missing or different unit (exactly `dimensionless` and `ueV`, as in the §3.2 example and atomforge's writer, not `µeV`); a value that is neither null nor a positive finite number; a `maps_to` naming another field. `loss_tangent` keeps `maps_to: null`: it is converted, never copied.
- Ignored, as additive fields on v0: `structures_ref`, `atomforge_computed`, `not_in_contract`, `open_stack`, `notes` and any other top-level block, and a quantity's other members (`reason`, `relation`, `conditions`). An unknown *quantity* is an error, not ignored, so that a number qpuforge cannot map is never dropped silently.

Not implemented: the `conditions` of a value are not interpreted, so Q and tanδ are not scaled to the qubit's `f01` or `T_eff_mK` (open question 3); there is no participation weighting, and a sidecar sets the one chip-wide `Q_dielectric` (question 1); only `analyze` takes `--materials`, so the other sub-commands still need the values copied into the chip JSON by hand (question 6).

Tests: the `materials` group (label `unit`, 25 checks) and six CLI guards (`cli_materials*`, `cli_usage_materials`). The fixtures are in `tests/data/`. `materials-contract-CuAu.json` is byte for byte what atomforge's `materials_contract_json()` writes for the synthetic Cu–Au search result and fixed options of its `tests/test_materials_contract.cpp`.

---

## 1. What atomforge outputs today

All of this is produced by `include/atomforge/io.hpp` and `src/main.cpp` (`atomforge search <A> <B> [--out DIR]`). The energy model on `main` is `ZhouEAM` (`include/atomforge/potential.hpp`), named `"Zhou04-EAM"` in the report.

### 1.1 Structures

Every candidate is written twice, as the README states:

- **POSCAR** (`write_poscar`) — VASP / Quantum ESPRESSO via ASE. Lattice in Å, `Direct` fractional coordinates, species grouped. Comment line is `Structure::label`.
- **Extended XYZ** (`write_extxyz`) — OVITO / ASE / VESTA. Lattice in Å, Cartesian positions in Å, header field `energy=` is the **total** energy in eV (`e_per_atom * natoms`). Also writes `pbc="T T T"` and `label=`.

There is no trajectory output, no charge file, and no phonon output. Dynamical stability, when run, is an MD heat-and-quench stored only in the report's `screen` block, not as a structure series.

### 1.2 Search report (`report.json`)

`write_report_json` emits a single JSON object. A committed example is `results/CuAu/report.json` (and the same schema for CuAg and NiAl). Fields:

| Field | Meaning | Units |
|---|---|---|
| `system` | `[A, B]` element symbols | — |
| `potential` | energy-model name (today `"Zhou04-EAM"`) | — |
| `reference.<el>.phase` | elemental reference structure label | — |
| `reference.<el>.mu` | elemental chemical potential | eV/atom |
| `hull[]` | convex-hull vertices `{x, ef, id}` | `x` is composition fraction of B (dimensionless); `ef` is formation energy in eV/atom |
| `candidates[]` | all relaxed candidates (see below) | |

Each candidate:

| Field | Meaning | Units |
|---|---|---|
| `id`, `label`, `source`, `parent` | identity and provenance | — |
| `natoms` | number of atoms | — |
| `x` | composition fraction of B | dimensionless |
| `e_per_atom` | total energy / N | eV/atom |
| `ef` | formation energy vs elemental refs | eV/atom |
| `e_hull` | energy above the convex hull | eV/atom |
| `vol_per_atom` | cell volume / N | Å³/atom |
| `converged` | FIRE relax success | bool |
| `screen` (optional, if `screened`) | stability block | see below |

`screen` (from `analysis.hpp`): `verdict` (e.g. `STABLE (on hull, mechanically + thermally stable)`), `bulk_GPa`, `shear_GPa`, `min_C_eig_GPa`, `C11`, `C12`, `C44` (GPa), `md_T` (K), `md_rms_disp` (Å), `dE_quench` (eV/atom), `max_site_shift` (Å).

Other units used in the code (`include/atomforge/math.hpp` `units`): `eV_per_A3_to_GPa = 160.21766208`, `kB = 8.617333262e-5` eV/K. Forces on the `Potential` interface are eV/Å; virial is eV (stress = virial / volume).

### 1.3 Planned, but not produced, DFT quantities

`quantum/pyscf_dft_probe.py` is a standalone fcc-Al PBE/GTH probe. It is **not** wired into the search. When run it prints:

- total energy in **Hartree** (`E=... Ha`);
- a force / gradient tensor (PySCF nuclear gradient; typical units Ha/Bohr, not converted);
- a stress tensor from `pyscf.pbc.grad.krks_stress` (PySCF convention, not converted to GPa or eV/Å³).

The README's planned `ExternalPotential` would stream structures to a persistent Python worker and must return the `EvalResult` the relaxer already expects: `energy` in eV, `forces` in eV/Å, `virial` in eV. That is still structure-and-energy data, not dielectric loss.

**atomforge does not currently output** a loss tangent, a dielectric quality factor, a T1, a participation ratio, a superconducting gap, or any per-layer stack description. Those would have to be computed by a later DFT (or FEM) step that does not exist in either repo.

---

## 2. What qpuforge needs today

Cited from `include/qpuforge/device.hpp` on `release/v0.3`. Frequencies inside the device model are GHz (E/h) unless a name says otherwise (file header comment).

### 2.1 The only materials hook: chip JSON `materials`

`ChipDesign` stores three materials fields, defaulted in the struct and read/written by `design_from_json` / `design_to_json`:

| JSON key (`materials.*`) | `ChipDesign` field | Default | Used for |
|---|---|---|---|
| `junction_gap_ueV` | `gap_ueV` | 180 | Ambegaokar–Baratoff \(E_J\) from normal-state resistance (`EJ_from_Rn_GHz`) |
| `Rn_spread_pct` | `Rn_spread_pct` | 2.0 | Fabrication-yield Monte Carlo junction spread. Alpha designs use 2; Beta-1 uses 0.5 (`presets.hpp` `beta1()`, `designs/beta1/beta1.json`) |
| `Q_dielectric` | `Q_dielectric` | `3e6` | Dielectric-loss contribution to T1. Alpha-1 / v1.1 JSON uses 3e6; Beta-1 uses 6e6 |

These are **chip-wide scalars**. There is no per-qubit Q, no per-material stack, no tanδ field, and no participation-ratio field.

Related environment fields (same struct, under JSON `environment`, not `materials`) also affect coherence but are not materials outputs from atomforge: `T_eff_mK` (default 45; Beta-1 35), `flux_noise_uPhi0` (2.0), `charge_noise_e` (1e-4), `Tphi_other_us` (1e4).

### 2.2 How Q becomes T1

In `analyze()` (`device.hpp`), with `f = f01 * 1e9` (Hz):

```
G_diel    = 2 * π * f / Q_dielectric
G_purcell = κ * (g / Δ)² / max(purcell_filter, 1)
G1        = G_diel + G_purcell
T1_us     = 1e6 / G1          // also T1_diel_us, T1_purcell_us
T2_us     = 1e6 / (0.5 * G1 + Gphi)
```

`κ` comes from `readout.kappa_MHz`, `g` from `readout.g_MHz`, `Δ` from `(f_res - f01)` in Hz, and `purcell_filter` from the readout design (a power-suppression factor, not a material parameter). If `G_purcell > G_diel` the report warns that T1 is Purcell-limited.

Under this uniform-loss model, \(Q = 1 / \tan\delta\) would be consistent with \(G_\mathrm{diel} = 2\pi f\,\tan\delta\). qpuforge never reads a tanδ; it only reads `Q_dielectric`. A materials loop that produced loss tangents would have to convert (and decide whether a single effective Q is enough) before those numbers can enter a chip JSON unchanged.

### 2.3 What the docs say the loop should eventually provide

`README.md` Roadmap tooling item 5 (and the same idea in the baseline): "Materials loop: atomforge → DFT for junction barriers, oxides and TLS defects. Loss tangents → Q → T1."

Item 1 is a related but distinct qpuforge-side step: "Full-wave/FEM cross-check of the cells; participation ratios → TLS-limited T1 per element." The Limitations section states that participation ratios are **not** extracted by the current MoM layout (`layout.hpp`). Until that exists, even a perfect material tanδ cannot be participation-weighted inside qpuforge; the only place to put a number is the single `Q_dielectric`.

`ROADMAP.md` Phase 0 describes atomforge's programme use as "materials screening for junction and base-layer metals." The named stack in the atomforge README is Al/AlOx junction barriers, Nb and Ta surface oxides, TiN, and the Si substrate. EAM cannot represent O, N or Si.

`Rn_spread_pct` is a process-control input (junction trimming), not something an atomistic search produces. `junction_gap_ueV` is the superconducting gap Δ used in Ambegaokar–Baratoff (design docs: Al/AlOx junctions with Δ = 180 µeV). atomforge does not compute a gap today. Both can be reserved in a file schema; only `Q_dielectric` (and later tanδ / participation) is a realistic near-term atomforge→DFT product.

---

## 3. Proposed interchange format

Propose a **sidecar JSON file**, not a change to `report.json` and not a new C++ API. Chip designs continue to be valid if the sidecar is absent: qpuforge already defaults the three materials keys. A human, or a later small loader, would copy accepted values into `materials` in the design JSON.

### 3.1 Filename and versioning

- Filename pattern: `materials-contract.json` (or `materials-contract-<system>.json` next to a search `report.json`).
- Required field: `schema` = `"forge.materials.v0"`.
- Breaking changes bump the integer (`v1`, `v2`). Additive optional fields stay on `v0`.
- `produced_by` records the tool, git commit, potential/DFT level of theory, and timestamp (ISO-8601, include offset).

### 3.2 Suggested `forge.materials.v0` schema

```json
{
  "schema": "forge.materials.v0",
  "produced_by": {
    "tool": "atomforge",
    "commit": "ac56e4b69b87e6876c0fd80e54a93c6f65d310ae",
    "potential_or_method": "Zhou04-EAM | PBE/GTH-PySCF | ...",
    "timestamp": "2026-11-25T22:40:00+00:00"
  },
  "system": ["Al", "O"],
  "structures_ref": {
    "report": "results/AlO/report.json",
    "format_notes": "POSCAR in Å, extxyz energy in eV (total), report energies in eV/atom"
  },
  "quantities": {
    "Q_dielectric": {
      "value": 3000000,
      "unit": "dimensionless",
      "maps_to": "ChipDesign.materials.Q_dielectric",
      "provenance": "assumed | measured | dft-loss-tangent | participation-weighted",
      "conditions": {"T_K": 0.01, "f_GHz": 5.0}
    },
    "loss_tangent": {
      "value": null,
      "unit": "dimensionless",
      "relation": "If the uniform-loss model in device.hpp holds, Q_dielectric = 1 / loss_tangent.",
      "maps_to": null
    },
    "junction_gap_ueV": {
      "value": 180,
      "unit": "ueV",
      "maps_to": "ChipDesign.materials.junction_gap_ueV",
      "note": "Reserved. atomforge does not compute this today. Default in device.hpp is 180."
    }
  },
  "not_in_contract": {
    "Rn_spread_pct": "Process spread used by fab_monte_carlo(); not an atomistic output.",
    "environment": "T_eff_mK, flux_noise_uPhi0, charge_noise_e, Tphi_other_us stay in chip JSON."
  },
  "open_stack": [],
  "notes": []
}
```

Units in this file must be explicit on every quantity. Recommended conventions, matching the two codebases:

| Quantity | Unit | Notes |
|---|---|---|
| Structure lattice / positions | Å | Already true of POSCAR and extxyz |
| atomforge energies | eV or eV/atom | Never Hartree in the interchange file; convert PySCF output first |
| Stress / elastic | GPa | Convert from eV/Å³ or Ha/Bohr³ before writing |
| `Q_dielectric` | dimensionless | Same number that goes into chip JSON |
| `loss_tangent` | dimensionless | tanδ |
| `junction_gap_ueV` | µeV | Same as `materials.junction_gap_ueV` |
| Temperature, if quoted | K | qpuforge's `T_eff_mK` is millikelvin; convert at the boundary |
| Frequency, if quoted with Q or tanδ | GHz | Matches qpuforge device units |

`loss_tangent` is allowed on `v0` even though qpuforge cannot ingest it yet, so a DFT workflow has somewhere to put a number. Until a participation-weighted path exists, the only automatic mapping is `Q_dielectric = 1 / loss_tangent` under the uniform-loss assumption, and that mapping should be spelled out in `provenance` rather than hidden.

### 3.3 What this file is not

- Not a replacement for `report.json`, POSCAR or extxyz. Those remain the structure/energy products.
- Not a proposal to add `load_materials_contract()` or any other new symbol to `device.hpp`.
- Not a demand that every chip JSON be generated from atomforge. Alpha-1 and Beta-1 already set `Q_dielectric` by assumption (3e6 / 6e6).
- Not a hardware characterisation format. Measured T1/T2 would be a later, separate mapping into the same chip JSON fields.

---

## 4. Open questions

1. **Is a single chip-wide `Q_dielectric` enough for v0 of the loop?** The code can only consume that today. A participation-weighted, per-element T1 needs the FEM work in qpuforge README tooling item 1, which is not built.
2. **Who converts tanδ → Q?** atomforge (or its DFT worker), a small standalone script, or a human editor of chip JSON? The conversion is only unique under the uniform-loss model implied by `G_diel = 2πf / Q`.
3. **At what temperature and frequency is tanδ quoted?** qpuforge T1 uses the qubit `f01` (several GHz) and `T_eff_mK` (~35–45 mK). DFT loss estimates are often 0 K or 300 K. The `conditions` object is required for a reason; the method is not chosen.
4. **Which material is first?** The READMEs name Al/AlOx barriers, Nb/Ta surface oxides, TiN, Si, and TLS defects. EAM cannot start those systems; the DFT bridge is a prerequisite. Ni–Al is the documented EAM failure, but it is not itself a QPU-stack dielectric.
5. **Can `junction_gap_ueV` ever come from atomforge?** It is a superconducting gap, not a ground-state DFT total energy. If it stays a literature / measured input, the schema should keep it reserved or drop it from v0.
6. **Should accepted contract values be copied into design JSON by hand, or is a loader in scope later?** This draft assumes hand copy until someone designs a reader. Do not add an API in advance of that decision.
7. **How should failed or qualitative EAM structures be marked** so they are not silently used as DFT inputs (Ni–Al B2 is the example in `results/NiAl/report.json`)?
8. **File location and ownership:** this proposal lives at `docs/DATA_CONTRACT.md` in qpuforge. Still open: should a copy or link also live in atomforge, next to atomforge `results/`, or next to a qpuforge design folder?

---

## 5. Recommended next software step (not part of this proposal's implementation)

Once the DFT `ExternalPotential` can relax and evaluate a named QPU-stack system, write one `forge.materials.v0` file by hand from those results, set `Q_dielectric` (or tanδ plus an explicit conversion), and patch a copy of `designs/alpha1-v11/alpha1_v11.json` to see the T1 change through existing `analyze()`. That exercise validates the contract against the code that already exists. It does not require a new API.
