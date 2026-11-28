# Flux-line distortion term (q18, opt-in)

Status: **implemented, off by default.** Scope and decisions: `docs/design/q17-q18-feasibility.md` (q18), decided by
Finnian Bond on 2026-11-28: the term is **additive** to `residual_ctrl_err` and **off by default**, so no committed
result is regenerated.

## What it adds

The flux line between the AWG and the SQUID is a linear system. After predistortion, what is left of the step response
of each line (mover and coupler) is modelled as

    s(t) = 1 + static_offset + Σ_k a_k · exp(−t / τ_k)

`static_offset` is a gain or offset error that persists through the gate; each `(a_k, τ_k)` is an uncorrected exponential
tail, as a fraction of the step. This is the **fast path** of the feasibility note: no change to `simulate_cz`.

- Effective error over the hold of the CZ (`flux_line_eps_eff`, hold T = `t_flat_ns` of the calibrated pulse):
  `ε_eff = | static_offset + Σ_k a_k · (τ_k / T) · (1 − exp(−T / τ_k)) |`.
- Infidelity (`flux_distortion_err`), from the stored quasi-static curvature of 1 − F in flux
  (`curv = c_mover + c_coupler`, the `curv_flux` field of the CZ cache):
  `err_distortion = ½ · curv · (ε_eff · max(|ΔΦ_mover|, |ΔΦ_coupler|))²`.
  The cache stores only the sum of the two curvatures, so this is the upper bound of
  `½ (c_m (εΔΦ_m)² + c_c (εΔΦ_c)²)` (the "worst pair" column of the note's sensitivity table). Both idle fluxes are 0,
  so the steps are `flux1_int` and `fluxc_on`.
- `calibrate_chip_cz` adds it **after** the unchanged sum:
  `err_ctrl = (err_coherent − leakage) + err_flux_local + err_flux_quasistatic + residual_ctrl_err [+ err_distortion]`.
  `err_ctrl` is written to `PairReport::err_2q_ctrl`, which is what the QEC and circuit emulators read (not
  `analyze()`'s `e_2q`). `err_total` grows by the same amount. The leakage share is 0 in the fast path, so `leak_2q`
  is unchanged.

## Configuration

Design JSON, under `tunable_coupler`:

```json
"flux_line": {"enabled": true, "static_offset": 3e-4, "tails": [{"amplitude": 1e-3, "tau_ns": 200}]}
```

Absent, or `"enabled": false`, means off. A tail needs `tau_ns > 0` and a finite amplitude. The block is written back
only when enabled, so every design written without it is byte-for-byte what it was. `qpuforge cz <chip.json>` prints a
`flux-line distortion` line when the term is on, and `--out` then adds `flux_eps_eff` and `err_distortion` per pair.

## What "off" guarantees

Nothing is computed, no key is written, and `err_ctrl`, `err_total`, `err_2q_ctrl` and `leak_2q` are the same bits:
the extra term is a separate `if` after the unchanged expression. The CZ cache key and `CZModel::kVersion` are
unchanged, because the term is computed from cached quantities after the cache lookup, so existing cache entries stay
valid with the term on or off. The `flux_line` test group checks this on Alpha-1 v1.1 (configured-but-off and
on-with-zero-error both give identical calibration JSON / bits), and the thread-count determinism digest is unchanged.

## Limits

- No measured line response exists: every number is an assumption until hardware data arrives
  (`qpuforge.characterisation.v0`, later).
- `residual_ctrl_err` already budgets distortion among other things, so turning the term on double-counts by design.
  That is the decided behaviour: it is a sensitivity knob, not a re-split of the budget.
- Fast path only: ramp distortion (adiabaticity, |11>→|20> leakage, spectator Landau-Zener), coupler-strength change
  and gate-to-gate tail history are not modelled. A waveform-driven `simulate_cz` is the next step (note, q18 item 3).
