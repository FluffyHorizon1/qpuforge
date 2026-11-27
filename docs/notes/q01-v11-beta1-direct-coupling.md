# Alpha-1 v1.1 and Beta-1 direct coupling (inherited 1 MHz): intent and impact

Issue file: `q01-v11-g.md`. On `main` 78b9f8b. Local runs on 27 Nov 2026, not committed.

## What PR #2 did

`d7ab24d` changed `alpha1()` to `s.g = 0.001` (1 MHz bare direct coupling, 0.035 fF per link; it was 4 MHz,
0.14 fF) so that the new three-transmon coupler model leaves −4 to −14 kHz idle ZZ. `alpha1_v11()` starts from
`alpha1(d)` and `beta1()` from `alpha1_v11(d)`, so both silently moved to 1 MHz as well. PR #2 sets
`tc_model = false` for v1.1 ("v1.1 pulse coupler retains its existing idle model"), which suggests the v1.1 / Beta-1
models were meant to stay as they were; the coupling change looks unintended. **Finnian needs to confirm.**

## Impact, measured

`analyze()` on the presets (designed pairs):

| | before PR #2 (`1cfa2b1`) | `main` |
|---|---|---|
| `alpha1_v11(3)` direct g | 4.000 MHz | **1.000 MHz** |
| `alpha1_v11(3)` max static (pair-diagonalised) ZZ | 61.0 kHz | 3.8 kHz |
| `beta1(5)` direct g | 4.000 MHz | **1.000 MHz** |
| collisions | 0 | 0 |

Committed outputs (provenance `bff9d10`, generated at 4 MHz):
- **Design files change.** `save_design(alpha1_v11(3))` on `main` differs from the committed
  `designs/alpha1-v11/alpha1_v11.json` in 75 lines, and `beta1(5)` from `designs/beta1/beta1.json` in 227 lines: every
  link `C_fF` 0.14 → 0.035 fF, and every qubit's shunt `C_fF` (e.g. 95.66 → 95.87 fF) and `Rn_kOhm`, because the
  inverse design re-balances the capacitance matrix. `qpuforge analyze designs/alpha1-v11/alpha1_v11.json` (4 MHz)
  and `analyze(alpha1_v11(3))` (1 MHz) no longer describe the same chip.
- **v1.1 review physics barely moves.** `alpha1_v11_review --smoke` built from `main` and from the proposed fix, same
  seeds and the committed CZ cache (16 of 16 calibrations cached in both): every QEC, leakage, budget, sensitivity and
  distance-5 number is bit-identical; CZ calibration and post-layout values differ only at about 1e-10 relative (the
  inverse design converges to the same qubit E_C and E_J either way, and `pulse.hpp` sets g12 per pair itself).
  The full 100,000-shot review and the Beta-1 review were not rerun; the GDS files were not regenerated.
- So the visible damage is to design-file and `analyze()` reproducibility, and no test pinned it (186 checks pass).

## Proposed fix (branch `fix/alpha1-v11-keep-4mhz-coupling`, code only)

- `presets.hpp`: new `alpha1_with_coupling(d, g_GHz)` holds the old `alpha1()` body. `alpha1()` calls it with 1 MHz
  (unchanged behaviour); `alpha1_v11()` calls it with 4 MHz, as before PR #2. `beta1()` follows v1.1.
- Test: `alpha1_v11(3)` and `beta1(5)` designed pairs have g = 4 MHz (2 checks; 188 in total, 21/21 CTest on g++ 14.2
  and clang 19.1.7, `-Werror`).
- With the fix, `save_design(alpha1_v11(3))` and `save_design(beta1(5))` reproduce the committed `alpha1_v11.json` and
  `beta1.json` exactly, apart from the two keys PR #2 added (`coupler_EC_GHz`, `model_coupler: false`). No committed
  output needs regenerating.

If 1 MHz is intended for v1.1 and Beta-1 instead, keep `main`, then regenerate `designs/alpha1-v11/*`,
`designs/beta1/*` (reviews, design JSON, GDS, `BETA-1.md` via `tools/fill_beta_spec.py`) and `results/PROVENANCE.md`,
and add a CHANGELOG line.
