# AF-17 characterisation plan — from a cold chip to the Alpha hardware gate

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

| | |
|---|---|
| **Applies to** | Alpha-1 v1.1 (AF-17), 17 qubits, 24 tunable couplers. The same plan scales to Beta-1 with the Beta thresholds. |
| **Outcome** | One data file per cooldown (`qpuforge.characterisation.v0`) and an automatic verdict on the Alpha hardware gate of `ROADMAP.md` |
| **Tooling** | qpuforge v0.4: `qpuforge hwgate`, `qpuforge synth-chip` (`include/qpuforge/characterise.hpp`) |
| **Status** | Written before hardware exists. The tool chain has been rehearsed on a virtual chip (§7); no measurement has been taken. |

The Alpha phase closes only when measured numbers meet the hardware gate **and** the model, fed with those numbers,
predicts the measured logical error. This plan says what to measure, in which order, where each number goes in the
data file, and how the verdict is produced.

## 1. Order of work

Each stage depends on the one before it. A stage is done when its numbers are in the data file.

| Stage | Measurement | Result | Field in the data file |
|---|---|---|---|
| 0. Before cooldown | Room-temperature junction resistance on every junction and on the test structures | Spread against the design `Rn_kOhm`; predicted frequency shifts; trim decision | — (kept in the fabrication record) |
| A. Bring-up | Resonator spectroscopy; qubit spectroscopy against flux; π-pulse amplitude and DRAG | Resonator frequency and linewidth, sweet-spot frequency, working flux bias | `qubits[].f01_GHz` |
| B. Coherence | T1 (repeated over several hours), Hahn echo, Ramsey, residual excited-state population | Distribution of T1; echo T2; effective temperature | `T1_us`, `T2e_us`, `p_thermal` |
| C. Single-qubit gates | Randomised benchmarking on all qubits at once | Error per gate with neighbours active | `err_1q` |
| D. Readout and reset | Assignment matrix at the 500 ns integration time; active reset; leakage-removing reset | P(1 given 0), P(0 given 1); probability of \|1⟩ after reset | `ro_p01`, `ro_p10`, `reset_err`, `operation.leakage_removing_reset` |
| E. Couplers and CZ | Coupler spectroscopy; idle point; conditional Ramsey; CZ calibration; interleaved RB or XEB; leakage RB | Idle ZZ per pair; CZ duration, error and leakage per pair | `pairs[].zz_kHz`, `t_cz_ns`, `err_cz`, `leak_cz` |
| F. Error correction | Distance-3 memory in the Z and the X basis at 5, 15 and 25 rounds, at least 20,000 shots each | Logical failures; detection fraction per round | `memory[]` |

Notes on the stages that most often go wrong:

- **B. Use the median and report the worst.** T1 moves in time as two-level-system defects drift through the qubit
  frequency. Report the median of the repeated measurements for each qubit. The gate uses the median over qubits and
  prints the worst qubit beside it.
- **D. Measure readout with all 17 qubits read at once.** Readout crosstalk and measurement-induced dephasing of
  neighbouring data qubits do not appear in single-qubit readout tests. They do appear in the memory experiment, and
  they are the most likely cause of a model-versus-measurement failure (§6).
- **E. Calibrate the CZ with the three knobs the simulation uses:** the flux amplitude of the moving (measure) qubit,
  the flat-top duration, and the coupler-on flux. The third knob closes the off-resonant swap between \|01⟩ and \|10⟩.
  Calibrate the flux-pulse distortion first (step response), or the CZ error will carry a history-dependent phase.
- **E. Benchmark in the layers the code uses.** The surface-code round runs four layers of parallel CZs. Measure the CZ
  error in those layers, not pair by pair with the rest of the chip idle.
- **F. Record every shot's detection events,** not only the logical outcome. The detection fraction per round is the
  stability test, and it is the first place leakage shows up (a rising detection rate).

## 2. The data file

One JSON file per cooldown. Leave out anything that was not measured: the gate report then lists that criterion as
"not in the data" and the model keeps its own value for that quantity.

```json
{
  "schema": "qpuforge.characterisation.v0",
  "chip": "AF-17 wafer 3 die 12",
  "design": "designs/alpha1-v11/alpha1_v11.json",
  "date": "2027-03-14",
  "source": "measured",
  "operation": { "t_1q_ns": 25, "t_readout_ns": 500, "t_reset_ns": 200, "leakage_removing_reset": true },
  "qubits": [
    { "id": 0, "f01_GHz": 5.171, "T1_us": 71.0, "T2e_us": 64.0, "err_1q": 2.9e-4,
      "ro_p01": 0.006, "ro_p10": 0.011, "p_thermal": 0.004, "reset_err": 0.003 }
  ],
  "pairs": [
    { "a": 0, "b": 3, "t_cz_ns": 44.0, "err_cz": 1.9e-3, "leak_cz": 6e-5, "zz_kHz": -2.8 }
  ],
  "memory": [
    { "d": 3, "basis": "Z", "rounds": 25, "shots": 20000, "logical_errors": 2334,
      "detection_by_round": [0.031, 0.049, 0.050] }
  ],
  "notes": []
}
```

(The numbers above show the format only.)

| Field | Unit | Meaning |
|---|---|---|
| `qubits[].id`, `pairs[].a`, `pairs[].b` | — | Qubit ids of the design file |
| `T1_us`, `T2e_us` | µs | Energy relaxation time; Hahn-echo coherence time |
| `err_1q` | probability | Average error per single-qubit gate from simultaneous randomised benchmarking |
| `ro_p01`, `ro_p10` | probability | P(read 1 given prepared 0), P(read 0 given prepared 1) |
| `p_thermal` | probability | Excited-state population without reset |
| `reset_err` | probability | Population of \|1⟩ after the active reset |
| `err_cz` | probability | Average error per CZ from interleaved RB or XEB. It includes decoherence and leakage. |
| `leak_cz` | probability | Leakage per CZ from leakage RB |
| `zz_kHz` | kHz | Idle ZZ with the coupler parked |
| `t_cz_ns` | ns | Calibrated CZ duration, ramps included |
| `memory[].detection_by_round` | fraction | Detection fraction in each syndrome round, then one last entry for the final data readout |

## 3. How the model uses the numbers

`model_from_measurement` replaces the design's predicted values with the measured ones and leaves everything else
alone. Gate errors need one step of care: the emulator adds decoherence and leakage itself, so only the remainder of a
measured gate error is passed on as control error.

- Single qubit: control error = `err_1q` − (t/3)(1/T1 + 1/Tφ)
- CZ: control error = `err_cz` − (4/5)(t/2)·Σ(1/T1 + 1/Tφ) − `leak_cz`, the sum running over the two qubits

Both are floored at zero. Tφ comes from the echo time: 1/Tφ = 1/T2e − 1/(2 T1).

The single-qubit control error and the reset error enter the model as chip averages. Everything else is per qubit or
per pair.

## 4. Producing the verdict

```
qpuforge hwgate designs/alpha1-v11/alpha1_v11.json measured.json \
        --phase alpha --shots 50000 --cache results/cz_cache.json --out gate_report.json
```

The command exits with status 0 when the gate passes and 3 when it does not.

| Criterion | Threshold | Computed from |
|---|---|---|
| Median T1 | ≥ 50 µs | `T1_us` |
| Median T2 echo | ≥ 40 µs | `T2e_us` |
| Median single-qubit gate error | ≤ 1e-3 | `err_1q` |
| Median CZ error | ≤ 5e-3 | `err_cz` |
| Median CZ leakage per gate | ≤ 1e-3 | `leak_cz` |
| Median readout assignment error | ≤ 2 % | mean of `ro_p01` and `ro_p10` per qubit |
| Median idle ZZ | ≤ 20 kHz | absolute value of `zz_kHz` |
| Distance-3 logical error per round | ≤ 1.5 % | `memory[]`: fit over the round counts, mean of the bases present |
| Detection-rate drift over 25 rounds | < 10 % | longest distance-3 run: mean of syndrome rounds 2–6 against the last five |
| Model versus measurement | within a factor 2 | the same memory runs simulated with the measured parameters, per basis |

The gate passes only when every criterion was evaluated and met. For each median the report also prints the worst
qubit or pair: one bad element can pass the median test and still cap the logical error.

## 5. Reading a failure

Each failing criterion points at a line of the Alpha-1 v1.1 error budget (`ALPHA-1-v1.1.md` §4–5). The figures in
brackets are the change in the simulated logical error from that review.

| Failing criterion | Likely cause | First response |
|---|---|---|
| T1 | Dielectric loss, or Purcell decay through the readout | Surface treatment and base-layer material (doubling Q: −38 %); check the Purcell filter |
| T2 | Flux noise away from the sweet spot; residual photons in the resonators | Confirm sweet-spot biasing; improve readout-line filtering |
| CZ error | Flux-pulse distortion; spectator collisions; coupler-on point | Recalibrate predistortion; re-run the three-knob calibration (residual 2e-3 instead of 1e-3: +41 %) |
| CZ leakage | Ramp too fast, or the \|11⟩–\|20⟩ crossing is mis-set | Longer ramps; the leakage-removing reset tolerates about ten times the simulated rate |
| Readout | Low signal-to-noise or T1 decay during readout | TWPA bias; integration weights (300 ns readout at the same error: −25 %) |
| Idle ZZ | Coupler park point off the exchange null | Re-null by moving the coupler park frequency (post-layout value: 8.295 GHz) |
| Detection drift | Leakage accumulating | Verify the leakage-removing reset; consider data-qubit leakage removal |
| Model versus measurement | An error source that is not in the parameter list | See §6 |

## 6. When the model and the measurement disagree

A ratio above 2 means the memory experiment suffers from something the per-element measurements did not see. Check, in
this order:

1. **Measurement-induced dephasing and readout crosstalk** on data qubits while the measure qubits are read.
2. **Simultaneous-operation errors:** CZ and single-qubit errors measured in isolation, not in the code's layers.
3. **Leakage** above the benchmarked rate (compare the detection fraction per round with the model's).
4. **Slow drift** between calibration and the memory run (repeat the per-element measurements afterwards).
5. **Correlated events** (bursts of detection events in a few consecutive shots).

Add the source to the model (`NoiseModel`), re-run, and only then pass the gate. A ratio below 0.5 means the per-element
numbers are too pessimistic, usually because benchmarking included state-preparation and measurement error.

## 7. Rehearsal on a virtual chip

Until hardware exists the chain is tested on a virtual chip: the design with fabrication-like scatter (log-normal T1
spread of 25 %, 6 % of qubits on a two-level-system defect, scatter on dephasing, readout, CZ error and leakage), its
parameters "measured" with 3 % estimation error, and the memory experiment simulated on it. This is synthetic data; it
tests the software, not the chip.

```
qpuforge synth-chip designs/alpha1-v11/alpha1_v11.json virtual.json --seed 3 --shots 20000 --cache results/cz_cache.json
qpuforge hwgate     designs/alpha1-v11/alpha1_v11.json virtual.json --shots 20000 --cache results/cz_cache.json
```

Result of the rehearsal (virtual chip, seed 3; 20,000 shots per memory run):

- **All ten criteria were evaluated and met.** The virtual chip has one qubit on a defect, with a T1 of 22 µs against a
  median of 67 µs. The report shows it as the worst value beside the median.
- **Model versus measurement agrees.** The model with the "measured" parameters gives 0.523 % (Z) and 0.499 % (X) per
  round. The virtual chip gave 0.492 % and 0.498 %: ratios 0.94 and 1.00.
- **A hidden error source fails the gate, as it should.** With an 8 µs dephasing process acting during the memory runs
  and absent from the reported T2 (`--hidden-dephasing-us 8`), the measured error rises to 3.4 % and 4.2 %, the ratios to
  6.5 and 8.5, and two criteria fail: the logical error and model versus measurement.

The example files are `example_characterisation_synthetic.json` and `example_gate_report_synthetic.json` in this folder.
