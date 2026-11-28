# Design note: flux-line response, predistortion, and what the CZ needs from them

Status: model, tests and a study (`include/qpuforge/fluxline.hpp`, test group `flux_predistortion`,
`examples/flux_distortion_study.cpp` → `results/flux_distortion_study.{txt,json}`). Nothing in `calibrate_chip_cz` uses
the model, so no calibration, cache or review output changes. In the plan for flux-pulse distortion
(`docs/design/q17-q18-feasibility.md`) this is the full path: line model, inverse filter, waveform-driven CZ and
predistortion study. The plan's opt-in error term for the calibration is separate work; this note supplies its
coefficients (`cz_offset_response`) and says which of two assumptions it has to make (see "Two ways to count").

## Why

`simulate_cz` (pulse.hpp) assumes that the flux at each SQUID follows the programmed pulse exactly. Whatever the flux
lines do to the pulse is inside one number: the residual control-error allowance of the CZ (`residual_ctrl_err`, 1e-3
for Alpha-1 v1.1 and 5e-4 for Beta-1 and Gamma-1), which also covers drift, defects and crosstalk. The reviews show
that this allowance sets most of the CZ error and that the logical error follows it closely (Alpha-1 v1.1: 1e-3 → 2e-3
raises the error per round by 37 %). This note turns the flux-line part of the allowance into requirements that can be
measured.

Measured lines are not flat. Rol et al. [1] report that distortion at short times causes leakage, that distortion at
long times makes a flux pulse act differently depending on the pulses before it, and that iterated digital correction
brought their step responses to within about 0.1 % of ideal. One published coupler line reaches 97 % of its final value
after 60 ns and needs more than 650 ns for 99 % before correction [2].

## Model

**The line** (`FluxLine`). First-order sections (`FluxSection`) in cascade, each with the step response 1 + a·exp(−t/τ), then an
optional single-pole low-pass (given by its 10–90 % rise time) and a delay. a > 0 is an overshoot that decays (a
high-pass element such as a bias tee); a < 0 is slow settling (skin effect, filters). The static gain is 1: it belongs to
the flux-per-volt calibration. The cascade is expanded in partial fractions into one low-pass state per pole, so the
response to a drive that is held constant over each interval has a closed form (`FluxLineFilter`). There is no
time-stepping error, and an idle stretch of any length costs one step.

**The AWG** (`FluxAwg`). The programmed waveform is sampled at the middle of each sample and held; an optional DAC
rounds it to a step of full scale / 2^bits. With no sample rate given the waveform is held on a 0.05 ns grid.

**Predistortion** (`Predistorter`). One digital first-order inverse per fitted section, at the AWG sample rate. For one
section and an exact fit, the line output averaged over each sample equals the programmed sample exactly. The filter
inverts a *fit* (`FluxChannel::correction`), which may differ from the line: the difference is what the study varies.
The low-pass and the delay are not corrected.

**The CZ on the delivered flux** (`CZDrive`, `simulate_cz_drive`). The same propagator, rotating frame and fidelity as
`simulate_cz`, for any piecewise-constant mover and coupler flux. `cz_ideal_drive` reproduces `simulate_cz` (fidelity to
3e-16). The time grid of a gate is 0.1 ns on the ramps, 0.5 ns on the flat top, and a settling time after the pulse
when a line has a low-pass or a delay.

**Pulse history** (`cz_train_drive`). Each line carries a train of pulses; a gate is simulated on the flux its two
lines deliver in its own time window, which includes what earlier pulses left behind.

**Calibration on a distorted line.** A hardware calibration tunes the gate on the line as it is, so a distortion that
is the same every time is not an error. The model has two levels of this:
- `cz_matched_command` scales each programmed step so that an isolated pulse delivers the calibrated flux on average
  over its flat top. It needs no optimisation and is used throughout the study.
- the study's local three-knob re-calibration (mover step, flat-top length, coupler step, Nelder–Mead from the matched
  start), used for the sample rate and bandwidth cases, where the pulse shape itself changes.

### Two ways to count

The fidelity of `simulate_cz` allows free single-qubit Z corrections after the gate. On hardware these virtual Z gates
are calibrated once, in one context. The study therefore reports two numbers for a gate that runs in another context
than the one it was calibrated in:
- **free**: the virtual Z gates are re-calibrated for the case shown. Only the conditional phase and leakage count.
- **fixed**: they stay as calibrated (`cz_fidelity_in_frame`). The shift of the single-qubit phases counts too.

The measure qubit moves 0.3–0.7 GHz to the interaction point, where its frequency changes by about 4–6 GHz per Φ0, so a
flux offset shifts its phase in first order: 1.18 rad per mΦ0 (median of the 16 pair classes). That makes the fixed
number about ten times the free one (10.1–10.4 times for a constant offset on the mover line). The curvature stored
with each calibration (`CZCalibration::err_flux_quasistatic`) is the free form.

`cz_offset_response` returns both forms for a gate: the added 1 − F to second order in constant offsets on the two
lines, and the mover's phase per flux. The study's JSON has it for every pair class.

## Checks (`flux_predistortion`, 32 checks)

- Partial fractions of two sections and a low-pass against a step-by-step integration of the sections in sequence
  (8e-7); a section followed by its analytic inverse is the identity (to rounding); the filter's integral against the
  closed form, in one step and in six.
- A pulse through a line: the delivered area equals the programmed area (unit gain), two pulses equal the sum of each
  (3e-16), a delay is a pure shift.
- Predistortion at 1 GS/s: exact for one section (3e-17 of a 0.2 Φ0 step), 1e-9 for two; a 10 % amplitude misfit leaves
  10 % of the tail. An 8-bit DAC stays on its code grid.
- The CZ: the drive against `simulate_cz` (3e-16); ideal lines through the line model (7e-7, the two discretisations).
  A 1 % tail far slower than the gate equals a pulse scaled by 1.01, and the amplitude calibration removes it (to 1e-7).
- History: the flux offset of a second pulse against an independent quadrature of the tail (to 2e-3 of the offset); the
  shift of the mover's phase against the frequency shift of the bare transmon (−0.057 against −0.059 rad).
- `cz_offset_response` against the curvature recipe of `calibrate_cz_pair`, against a direct simulation with both lines
  offset (free and fixed, to 1 %), and its phase coefficient against the bare transmon (1114 against 1160 rad/Φ0).

## Results

All numbers are from `results/flux_distortion_study.txt`: the 16 calibrated pair classes of the CZ cache (model
version 4), and the memory cycles of Alpha-1 v1.1 (24 CZ per round) and Beta-1 (80). "Added error" is the change of
1 − F against the same gate on ideal lines, or against the gate on its own.

**1. A static offset.** A flux offset of 1e-4 of the mover's step (about 2e-5 Φ0) adds 1.4e-5 free and 1.3e-4 fixed
(median; worst class 3.9e-5 and 3.7e-4). It grows as the square: 3e-4 adds 1.2e-3 fixed. Per fraction of its step, the
coupler line is 9 times (free) to 40 times (fixed) less sensitive.

**2. One tail, gate on its own.** A −0.1 % tail on both lines, programmed as calibrated, adds 5e-4 at τ = 30 ns and
1.6e-3 for slow tails: it is an amplitude error. With the amplitudes re-calibrated the rest is at most 2.3e-5 (mean,
at τ = 10–30 ns, nearly all leakage) and below 1e-6 from 300 ns on. A gate's own tail is a matter for its calibration.

**3. One tail in the memory cycle.** Each measure qubit's line carries up to four pulses back to back, round after
round. Gates calibrated on their own, −0.1 % tail on every line, mean added error over the CZ gates of a round (max in
brackets):

| τ | Alpha-1 v1.1, first round: fixed | free | steady state: fixed | free |
|---|---|---|---|---|
| 10 ns | 1.6e-4 (4.4e-4) | 3.9e-6 | 1.6e-4 | 3.9e-6 |
| 30 ns | 1.35e-3 (3.3e-3) | 1.2e-4 (2.9e-4) | 1.35e-3 | 1.2e-4 |
| 100 ns | 1.31e-3 (3.9e-3) | 1.3e-4 (3.8e-4) | 1.31e-3 | 1.3e-4 |
| 300 ns | 3.5e-4 | 3.4e-5 | 4.3e-4 | 4.2e-5 |
| 1 µs | 4.3e-5 | 4.3e-6 | 2.8e-4 | 2.6e-5 |
| 10 µs | 4e-7 | 5e-8 | 3.2e-4 (8.3e-4) | 3.0e-5 |

- Tails of 30–100 ns act within a round: the pulse before is still settling. Tails of 1 µs and more act across rounds:
  in steady state every flat top sits lower by the tail amplitude times the duty cycle of the line.
- +0.1 % gives the same numbers, and −0.3 % gives 9.0 times more: the added error grows as a².
- Beta-1 is 8–12 % higher at 30–100 ns and 1.8 times higher in steady state for slow tails (a shorter round, 0.72 µs
  against 0.92 µs, and more lines with four pulses).
- A quasi-static estimate (flux offset averaged over the flat top, times `cz_offset_response`) reproduces the simulated
  means to within 3 % (fixed) and 8 % (free) from τ = 30 ns on; at 10 ns it gives 0.76 of the fixed number.

**The tolerance.** For a mean added error below a tenth of the allowance, the corrected step response may deviate by
at most this much, at the worst time constant (30–100 ns):

| chip | allowance | fixed | free |
|---|---|---|---|
| Alpha-1 v1.1 | 1e-3 | 0.027 % | 0.089 % |
| Beta-1 | 5e-4 | 0.018 % | 0.059 % |

For slow tails (1 µs and more) the limits are 0.06 % and 0.19 % (Alpha-1 v1.1), 0.03 % and 0.10 % (Beta-1). Gamma-1 has
Beta-1's gates and allowance; its cycle was not run.

**4. Predistortion of an example line.** Sections (−5 %, 15 ns), (−2.5 %, 150 ns), (−1.2 %, 1.5 µs): an assumption, chosen
to settle roughly like the line in [2]. Uncorrected, the gates do not work (added error 0.3). An exact fit removes the
error (6e-10). A fit whose amplitudes are 3 % off leaves 3.7e-3 fixed and 3.0e-4 free; time constants 10 % off leave
6.2e-3 and 5.9e-4. Leaving out the 1.5 µs section costs 4.2e-2 fixed in steady state and 3e-3 in the first round. For a
tenth of the Alpha allowance the section amplitudes have to be fitted to within 0.5 % (fixed) or 1.7 % (free), and the
time constants to within 1.3 % or 4.1 %.

**5. AWG and timing** (three pair classes, re-calibrated on the line):
- 2.4 GS/s with a rise time of 0.5 or 1.0 ns costs nothing (−4e-5 to +2e-5).
- A rise time of 2.0 ns costs 2–3e-4 with the 7 ns ramps of the present pulse, with or without sampling at 2.4 GS/s.
- 1 GS/s sample-and-hold costs 3–7e-4 at 1.0 ns and up to 1.1e-3 at 2.0 ns, nearly all leakage: one sample is a seventh
  of the ramp, and a single-pole filter does not smooth the steps enough.
- DAC, full scale 1 Φ0: half a step on both lines costs at most 8.6e-5 at 14 bits and 5e-6 at 16 bits.
- A timing skew between the mover and the coupler line costs 1.9e-5 at 0.1 ns and 8e-5 at 0.2 ns (worst class); a tenth
  of the Alpha allowance allows 0.21 ns. The three-knob calibration cannot remove it.

## What this means for the chips

1. **Predistortion is required, and to 0.1 % it is not enough on its own.** A residual tail of 0.1 % with a time
   constant of 30–100 ns costs the CZ gates of a round 1.3e-3 each on average when they are used exactly as calibrated
   in isolation: more than the whole Alpha allowance.
2. **Calibrate the virtual Z gates of every CZ in the cycle it runs in.** That removes the first-order term: at 0.1 %
   the cost falls to 1.2e-4, about a tenth of the Alpha allowance and a quarter of Beta's. The characterisation plan
   should say so, and should measure the step response of every flux line (Cryoscope [1]) from 5 ns to 10 µs.
3. For Beta-1 and Gamma-1 the same criterion needs 0.06 %, or a pulse that cancels the offset: a bipolar (net-zero)
   pulse makes both the conditional phase and the mover's phase first-order insensitive to slow offsets [3]. The
   measure qubits idle at their sweet spot, so such a pulse is available in this design. It needs its own calibration
   and is the natural next step.
4. Slow tails matter through the duty cycle. An experiment that changes its pulse pattern (the first rounds of a
   memory, a lattice-surgery merge) shifts every gate by the steady-state numbers above until the lines settle.
5. Hardware: 2.4 GS/s with a rise time of 1 ns or less works and 1 GS/s sample-and-hold does not (or the ramps get
   longer); 14-bit DACs or better at 1 Φ0 full scale; the two pulses of a CZ aligned to 0.2 ns.

## Limits

- No measured line exists. The example line is an assumption; the tolerances per time constant do not depend on it.
- Sections are first order and real; reflections and ringing are not modelled. The low-pass is a single pole, which
  matters for the sample-and-hold result: a steeper reconstruction filter would do better at 1 GS/s.
- Each gate is a two-qubit simulation in its own window. The mover's other couplers are parked, and the lines do not
  couple to each other (no flux crosstalk). Lines with a low-pass get a settling time between pulses.
- The cycle is the compiled memory experiment; data qubits are not pulsed (no DQLR). Tails that remain while a qubit
  idles at its sweet spot shift its frequency only in second order and are left out (of order 1e-7 per round at 0.1 %,
  from the curvature of the transmon spectrum).
- The amplitude calibration is a rule, not an optimisation. A three-knob calibration in the cycle would also remove the
  steady-state part of the free numbers; the fixed and free numbers bracket what is left without it.
- The study does not change any logical-error number.

## Open points

1. An error term in the calibration (`err_distortion`): it has to choose between the free and the fixed coefficients,
   that is, state whether the virtual Z gates are calibrated in the cycle. The stored curvature is the free form. A
   bound from that curvature and the larger of the two flux steps (the plan's fast form) lies 2.6–4.8 times above the
   free number for equal fractional offsets on both lines, and 2.6–4.5 times below the fixed one. What enters as the
   offset matters as much: a gate's own tail is removed by its calibration (result 2); the part that is left comes from
   the pulses before it (result 3).
2. A net-zero CZ in `pulse.hpp` and its calibration, then this study again.
3. A measured step response per line in `qpuforge.characterisation.v0`, so that `qpuforge hwgate` can run the study's
   third section on the chip's own lines.
4. Outside this note: a local 1 − F polish of the 5.81/5.29 GHz calibration, the one class the leakage fallback of CZ
   model version 4 calibrated, gives 1 − F = 8.1e-5 with leakage 7.3e-5 (stored: 3.55e-4 and 1.6e-5), within both
   limits. A polish step after the fallback would lower the largest coherent error of all three chips from 3.6e-4 to
   2.1e-4; it changes the CZ cache and every review, so it is left for a model revision.

## References

1. M. A. Rol et al., "Time-domain characterization and correction of on-chip distortion of control pulses in a quantum
   processor", Appl. Phys. Lett. 116, 054001 (2020), arXiv:1907.04818.
2. "Sensitivity-Adapted Closed-Loop Optimization for High-Fidelity Controlled-Z Gates in Superconducting Qubits",
   arXiv:2412.17454 (its measured step response of a coupler flux line, before correction).
3. M. A. Rol et al., "Fast, High-Fidelity Conditional-Phase Gate Exploiting Leakage Interference in Weakly Anharmonic
   Superconducting Qubits", Phys. Rev. Lett. 123, 120502 (2019), arXiv:1903.02492.
