# Design note: CZ leakage above the Alpha limit (CZ model version 4)

Status: adopted in CZ model version 4 (`CZModel::kVersion`). `results/cz_cache.json` was regenerated with it.

## Problem

The Alpha simulation gate (ROADMAP.md) requires, for every pair, CZ leakage of at most 2e-4 per gate and a coherent
error (1 − F_avg) of at most 5e-4. Alpha-1 v1.1, Beta-1 rev C and Gamma-1 share 16 frequency combinations (pair
classes). Under CZ model version 2 one class, mover 5.81 GHz / target 5.29 GHz (Q8–Q11 on AF-17), leaked 2.33e-4 per
gate with a coherent error of 4.06e-4. The mean leakage over the 16 classes was 6.54e-5 (5.7e-5 over the 24 AF-17
pairs). The next class, 5.87/5.29 GHz, leaked 1.20e-4: within the limit.

## Cause

`optimize_cz` minimised `1 - F_avg` alone. In the average gate fidelity a leaked population counts no more than a
conditional-phase error of similar size, so the optimiser has no particular reason to avoid leakage. For the two
classes with target 5.29 GHz and movers 5.81 and 5.87 GHz, the calibration landscape has two branches: a flat top of
about 22 ns and one of about 26 ns (one full |11⟩ → |20⟩ → |11⟩ cycle at a slightly different coupler "on" point). The
coarse scan that seeds Nelder–Mead lands on the 22 ns branch, which leaks several times more:

| class (GHz) | branch | 1 − F_avg | leakage | 1 − F_avg − leakage | t_cz | total CZ error |
|---|---|---|---|---|---|---|
| 5.81 / 5.29 | short (plain objective) | 4.06e-4 | 2.33e-4 | 1.74e-4 | 36.2 ns | 2.031e-3 |
| 5.81 / 5.29 | long (leakage-weighted repeat) | 3.55e-4 | 1.60e-5 | 3.39e-4 | 39.7 ns | 2.042e-3 |
| 5.87 / 5.29 | short (plain objective) | 2.08e-4 | 1.20e-4 | 8.8e-5 | 35.6 ns | 1.841e-3 |
| 5.87 / 5.29 | long (leakage-weighted repeat) | 2.68e-4 | 1.94e-5 | 2.49e-4 | 39.8 ns | 1.977e-3 |

"Total CZ error" includes the residual control allowance and decoherence during the gate (`err_total`). The long branch
trades leakage for non-leakage coherent error and a longer gate. For 5.81/5.29 that trade is needed to meet the 2e-4
leakage limit, and 1 − F_avg improves; the total error is about the same (+0.5 %). For 5.87/5.29 the trade is not needed
(the short branch is within the limit) and it is worse on 1 − F_avg and on total error, so it is not taken.

The motivation is the gate criterion, not the logical error rate: leakage is only a few per cent of the Alpha-1 v1.1
error budget (§4 of `designs/alpha1-v11/ALPHA-1-v1.1.md`), and this change is not expected to move the memory
experiments beyond their Monte Carlo error.

## Options considered

1. **Change the frequency plan for the class** (move the target or the mover). Changes the chips' frequency
   allocation, the collision analysis and every downstream number for one calibration artefact. Rejected as the most
   invasive.
2. **Change the pulse shape** (longer sin² ramps, or another ramp). Changes every pair's gate time and all cache keys,
   and is not needed: a solution within the limit already exists with the present shape.
3. **Penalise leakage for every pair**: `1 - F_avg + w * leakage`. With `w = 10` the worst class reaches 1.6e-5, but
   the penalty also trades phase error for leakage where none is needed: class 5.81/5.17 went from 1 − F = 3.6e-5 to
   4.3e-4. Rejected.
4. **Hinge penalty** `w * max(0, leakage - 1e-4)` for every pair. Non-smooth at the cap; Nelder–Mead stalled on the
   kink for class 5.93/5.29 (1 − F = 2.0e-3, leakage pinned at 1.000e-4). Rejected.
5. **Leakage fallback on leakage alone** (CZ model version 3, never released): repeat the calibration with
   `1 - F_avg + 10 * leakage` whenever leakage exceeds 1e-4, and keep the repeat if it leaks less. Rejected in review:
   it re-tuned 5.87/5.29, which already met the limit, to a pulse with worse 1 − F_avg (2.08e-4 → 2.68e-4), a longer gate
   and a higher total error, with no guard on any of them.
6. **Guarded leakage fallback** (adopted, version 4). Calibrate as before. Only if the result fails the leakage limit
   (`leakage > leak_limit = 2e-4`) calibrate again with `1 - F_avg + 10 * leakage`, and keep the repeat only if it meets
   the limit, its 1 − F_avg is no worse than the plain one's, and it stays within `coh_limit = 5e-4`. Otherwise the plain
   calibration stands and the gate reports the failure.

Option 6 is the least invasive: 15 of the 16 classes keep their version-2 calibration bit for bit, and only the class
that failed the gate changes. It is also what one would do on hardware: re-tune a pair that fails a leakage check, and
keep the new tune-up only if it is not worse otherwise.

## Implementation

- `optimize_cz(S, fluxc_on, t_rise, dt, max_eval, leak_weight = 0, leak_limit = 1, coh_limit = 1)`: the fallback
  above; `leak_weight = 0` gives the version-2 behaviour.
- `CZModel`: `kVersion = 4`, `leak_limit = 2e-4`, `coh_limit = 5e-4`, `leak_weight = 10`. `cz_cache_key` encodes all
  three, so a cache from another model version is rejected.
- `calibrate_cz_pair` passes the `CZModel` settings.
- Test: `pulse_cz` calibrates 5.81/5.29 and 5.87/5.29 with and without the fallback. 5.81/5.29 must start above the
  limit and end within it with 1 − F no worse; 5.87/5.29 must keep its plain calibration exactly.
- `alpha1_v11_review` computes the gate from the per-pair maxima (`exit_criteria`).

## Result

Over the 16 classes: largest leakage 1.20e-4 (5.87/5.29), mean 5.18e-5; largest coherent error 3.55e-4 (5.81/5.29).
The per-pair values are in the three reviews and in `results/cz_cache.json`.
