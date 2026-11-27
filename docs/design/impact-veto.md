# Design note: per-shot impact injection and an impact veto

Status: implemented on `scs-bot/gamma-decoding` (`ImpactSpec`, `VetoSpec`, `calibrate_veto`, `burst_floor_veto` in
`include/qpuforge/qec.hpp`).

## Background

In 0.4, high-energy impacts enter only through `burst_floor()`. It runs memory experiments with a chip-wide excess
decay or dephasing rate at a ladder of severities, then integrates the excess logical error over the burst profile
s(t) = exp(−t/τ). The decoder keeps its quiet-chip priors. Nothing in a shot "sees" an impact start, and nothing
reacts to one.

## Per-shot injection

`FrameRunOptions::impact` (`ImpactSpec`: a `BurstSpec` and `p_shot`). With probability `p_shot`, a shot contains an
impact:
- the onset time t0 is uniform over the shot, and the centre is a uniformly chosen qubit;
- from the onset, qubit q has severity s = exp(−(t − t0)/τ) · g(q), where g = exp(−r²/2σ²) of its lattice distance r
  from the centre (σ = `BurstSpec::radius` in data-qubit pitches; σ ≤ 0 means the whole chip);
- the excess rates s/T1_peak and s/Tφ_peak are applied in each time slice after the quiet noise. They are
  state-aware, as in `build_frame_program`: a qubit known to be in |1⟩ decays as a bit flip, one in |0⟩ is
  untouched, and one in superposition gets the Pauli-twirled channel. For this the frame program now records the
  known state of every qubit in every slice (`FMoment::known`).

Every draw of the impact comes from its own generator, keyed on (seed, shot) by `impact_shot_seed` (as
`frame_shot_seed`, with a different constant). A shot without an impact is therefore bit for bit the same as in a run
without injection. The test checks this with impacts that carry no excess noise in half of the shots: identical
failure counts and detection rates.

## Veto

`FrameRunOptions::veto` (`VetoSpec`). The shot's burst statistic is the largest number of detection events in a
space-time window. A window is the detectors whose stabilisers lie within `radius` lattice units (Chebyshev) of a
centre stabiliser, over `window_rounds` consecutive rounds; every centre and start round is a candidate. It is
computed from the fired detectors only: each fired detector increments the windows that contain it. A shot whose
statistic reaches `threshold` is discarded (post-selection). It is still decoded, so the result reports failures
with and without the discard. `MemoryResult` gains the impact and discard counts, `p_logical_kept`, `eps_kept` and
the statistic's histogram.

`calibrate_veto(M, shots, seed, target, spec)` runs burst-free shots, records the histogram, and returns the
smallest threshold whose tail fraction is at most `target`.

Window choice, on the d = 5 test chip (10 rounds, 5% of shots with a chip-wide T1 impact, T1_peak 1 µs, false-veto
target 1e-3): 3 rounds × radius 2 caught 77% of impacts. 2 × 1 caught 57%, 1 × 1 caught 29%, and the whole chip over
3 rounds caught 77%. The whole chip over the whole shot caught 71% and is better only for weak, long, chip-wide
severities: 26% against 10% of shots vetoed at s = 0.1. The default is 3 × 2. The impacts it misses mostly start late
in a shot, with too few rounds left to show.

## Floor over the burst profile

`burst_floor_veto` runs `burst_floor`'s ladder of chip-wide severities with the veto on. Each point gives the
discarded fraction v(s) and the logical error of the kept shots. Two figures follow:
- **shot-by-shot discard:** rounds lost = ∫ (1 − v(s)) · max(0, ε_kept(s) − ε_kept(0)) dt / t_round, and the time
  discarded = ∫ (v(s) − v(0)) dt;
- **latched veto:** once the veto fires at the onset, the chip is held for a time h. Only the tail after h costs
  logical errors: rounds lost = ∫_{t>h} max(0, ε(s) − ε0) dt / t_round, from the no-veto curve interpolated in
  log s (below the ladder the excess is taken as linear in s, as in `burst_floor`). The dead-time fraction is
  rate × h. Failures before the veto first fires are what the per-shot injection measures.

The probability that an impact flips the logical, p_fail = ½(1 − e^(−2·rounds lost)), and the floor
rate · t_round · p_fail are formed as in `burst_floor`.

What it shows: with shot-by-shot discard the floor does not come down. The burst lasts thousands of rounds, and the
weak part of it, which the veto cannot tell from noise in one shot, still randomises the logical qubit. On the d = 5
test chip (T1 burst, τ = 25 ms) the floor is 4.3e-8 per round with and without the discard. A latched hold of 10 τ
brings it to 9.8e-9 at 2.5% dead time; 15 τ gives 7.5e-11 at 3.75%. These points lie below the measured ladder and
use its linear-tail extrapolation. Gamma-1 numbers: `results/decoder_study.txt`, section 4.

## Validation

`impact_veto` (tests/test_qpuforge.cpp), d = 5, 8 rounds:
- impacts without excess noise leave every count unchanged, and the impact rate matches `p_shot`;
- a threshold calibrated for 2e-3 holds on fresh burst-free shots (3.75e-4 measured);
- chip-wide T1 impacts in 5% of shots: 69% caught. The impact excess in ε falls from 0.00413 − 0.00147 to
  0.00198 − 0.00147 (impacted-shot failures 422 → 84);
- a local impact (σ = 1) does less damage than a chip-wide one (0.14 against 0.43 failures per impact);
- d = 3 dephasing burst: the veto floor does not exceed the plain floor, the latched floor falls monotonically with
  the hold, and at 20 τ it is below 1% of the plain floor.

CLI: `qpuforge qec chip.json --bursts --veto [--veto-window R] [--veto-radius L] [--veto-false P] [--impact-p P]
[--burst-radius L]`.
