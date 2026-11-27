# Design note: CZ leakage above the Alpha limit (CZ model version 3)

Status: adopted in CZ model version 3 (`CZModel::kVersion`). `results/cz_cache.json` was regenerated with it.

## Problem

The Alpha simulation gate (ROADMAP.md) requires CZ leakage of at most 2e-4 per gate. Under CZ model version 2 the
worst of the 16 frequency combinations shared by Alpha-1 v1.1, Beta-1 rev C and Gamma-1, mover 5.81 GHz / target
5.29 GHz (Q8–Q11 on AF-17), leaked 2.33e-4 per gate (coherent error 4.06e-4). The mean over the 16 classes was 5.7e-5.
A second class, 5.87/5.29 GHz, leaked 1.20e-4: inside the limit, but with little margin.

## Cause

`optimize_cz` minimised `1 - F_avg` alone. In the average gate fidelity a leaked population counts no more than a
conditional-phase error of similar size, so the optimiser has no particular reason to avoid leakage. The calibration
landscape has two branches for these pairs: a flat top of about 22 ns and one of about 26 ns (one full |11⟩ → |20⟩ →
|11⟩ cycle at a slightly different coupler "on" point). For the two classes with the smallest mover–target detuning,
the coarse scan that seeds Nelder–Mead lands on the 22 ns branch, which leaks several times more. The 26 ns branch
exists for these pairs and is better on both counts:

| class (GHz) | branch | 1 − F_avg | leakage | flat top |
|---|---|---|---|---|
| 5.81 / 5.29 | v2 (plain objective) | 4.06e-4 | 2.33e-4 | 22.2 ns |
| 5.81 / 5.29 | v3 (leakage fallback) | 3.55e-4 | 1.60e-5 | 25.7 ns |
| 5.87 / 5.29 | v2 (plain objective) | 2.08e-4 | 1.20e-4 | 21.6 ns |
| 5.87 / 5.29 | v3 (leakage fallback) | 2.68e-4 | 1.94e-5 | 25.8 ns |

So this is a calibration problem, not a limit of the pulse shape or of the frequency plan.

## Options considered

1. **Change the frequency plan for the class** (move the target or the mover). Changes the chips' frequency
   allocation, the collision analysis and every downstream number for one calibration artefact. Rejected as the most
   invasive.
2. **Change the pulse shape** (longer sin² ramps, or another ramp). Changes every pair's gate time and the cache keys
   of all classes, and is not needed: a low-leakage solution already exists with the present shape.
3. **Penalise leakage in the objective for every pair**: `1 - F_avg + w * leakage`. With `w = 10` the worst class
   reaches 1.6e-5, but the penalty also trades phase error for leakage where none is needed: class 5.81/5.17 went from
   1 − F = 3.6e-5 to 4.3e-4. Rejected.
4. **Hinge penalty** `w * max(0, leakage - 1e-4)` for every pair. Non-smooth at the cap; Nelder–Mead stalled on the
   kink for class 5.93/5.29 (1 − F = 2.0e-3, leakage pinned at 1.000e-4). Rejected.
5. **Leakage fallback** (adopted). Calibrate as before. Only if the result leaks more than `leak_cap = 1e-4` (half the
   gate limit), calibrate again with the objective `1 - F_avg + 10 * leakage`, and keep that result if it leaks less.

Option 5 is the least invasive: the 14 classes under the cap get exactly the version-2 calibration (the plain branch
is bit-identical), and only the two classes above it change. It is also what one would do on hardware: re-tune the
pairs that fail a leakage check.

## Implementation

- `optimize_cz(S, fluxc_on, t_rise, dt, max_eval, leak_weight = 0, leak_cap = 1)`: the fallback above;
  `leak_weight = 0` gives the version-2 behaviour.
- `CZModel`: `kVersion = 3`, `leak_cap = 1e-4`, `leak_weight = 10`. `cz_cache_key` encodes both (`_lw10_lc0.0001`),
  so a version-2 cache is rejected and its entries cannot be reused by mistake.
- `calibrate_cz_pair` passes the `CZModel` settings.
- Test: `pulse_cz` calibrates class 5.81/5.29 both ways and checks that the plain calibration is on the leaky branch
  and that the fallback brings it under `leak_cap / 2` with 1 − F below 5e-4.

## Result

See the "CZ calibration" lines of the three reviews and `results/cz_cache.json` (16 classes, regenerated from an
empty cache).
