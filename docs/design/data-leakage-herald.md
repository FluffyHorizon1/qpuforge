# Design note: heralding data-qubit leakage through DQLR

Status: implemented and merged into `main` in PR #35 (`LeakageSpec::dqlr_transfer`, `NoiseModel::dqlr_transfer`,
`FrameProgram::dqlr_transfer`, `qpuforge qec --leakage-removal lrr+dqlr --herald --dqlr-herald`). Off by default.
The measured result is negative on every case tried: see "Results".

## The problem

The herald of `heralded-leakage.md` is the three-state readout of the measure qubits. Data qubits are not read out
during the memory, so their leakage is invisible. On Gamma-1 the measure qubit is the moving (leaking) partner of every
CZ. A data qubit becomes leaked only when a leaked measure qubit hands its leakage over (transport, 5% per CZ). With
the leakage-removing reset alone it then stays leaked until transport hands the leakage back or it seeps out, which
takes several rounds, and in each of them the four neighbouring measurements are randomised.

## Which herald is physically available

Data-qubit leakage removal (DQLR) is already in the leakage model: after the measure qubit's reset, a pulse couples the
data qubit's |2⟩ to that measure qubit. In the 0.4 model the excitation is then dumped: the pulse sits inside the reset
slot, the reset takes it away, the data qubit returns in a random computational state with probability `dqlr_eff`
(0.98), and nothing records that it happened.

The readout that can see a data qubit's leakage is therefore the measure qubit's, and only if the excitation is still
on the measure qubit when it is next read. The model here is:

- The removal pulse comes after the measure qubit's leakage-removing reset, as before. With probability `dqlr_eff` it
  moves the leakage from the data qubit onto the measure qubit: the data qubit returns randomised (X and Z, each with
  probability 1/2), and the measure qubit is now leaked. Physically this is a full exchange |2⟩_data |0⟩_meas →
  |0⟩_data |2⟩_meas, for example the leakage iSWAP |20⟩ → |11⟩ followed by the |11⟩ → |02⟩ exchange that the CZ uses.
  This is an assumption of the model; no pulse for it has been simulated.
- Nothing dumps the excitation. The measure qubit stays leaked through the next round: its CZs put a random Z on each
  data partner, it can hand the leakage on by transport, and its syndrome bit for that round is lost (it reads |1⟩).
- **Which readout, and when.** The herald is the three-state readout of that measure qubit at the end of the round
  after the pulse: the same readout that heralds the measure qubit's own leakage, with the same fidelity and
  false-herald rate (`HeraldSpec::fidelity`, `false_rate`). For a pulse after round r, the herald arrives with the
  readout of round r + 1 and says: one of this measure qubit's DQLR partners was leaked before round r + 1, or this
  measure qubit leaked during round r + 1. The leakage-removing reset after that readout clears the measure qubit.
- If the measure qubit is already leaked when the pulse arrives, the pulse does nothing and the data qubit stays
  leaked. This cannot happen with the leakage-removing reset, except when two data partners of one measure qubit are
  leaked in the same round.

The extension to the leakage model is this one flag. `FrameProgram::FMoment::dqlr_to` records which measure qubit each
pulse couples to (the pairing of `prepare_memory`, unchanged). With the flag off, the sampler makes the same random
draws as before, so every earlier result is unchanged. `--dqlr-herald` needs `--leakage-removal lrr+dqlr`: without the
leakage-removing reset the transferred excitation would stay on the measure qubit until it seeps out.

A herald with less damage would read the measure qubit a second time right after the pulse (the |20⟩ → |11⟩ step
alone, then a two-state readout and a reset). That adds a readout and a reset to every round. It needs a change to the
circuit schedule and is not modelled here.

## Decoder: the window the leakage can have lasted

The weights come from the leakage chain of `derived-herald-weights.md`, which follows the excitation through the
transfer. For a heralded readout of measure qubit m in round r + 1, the backward pass finds:

- m leaked at one of its own CZs in round r + 1, with the consequences of `derived-herald-weights.md`;
- a data partner d was leaked and transferred after round r. The pass then continues backwards on d: through round r,
  and through earlier rounds with probability 1 − `dqlr_eff` for each removal that must have failed, until it finds
  the CZ where the leakage started. Every CZ of d in that window puts Z on the measure qubit on the other side with
  probability 1/2, weighted by the probability that d was leaked there; the transfer itself randomises d.

So the data qubit's edges are re-weighted over the window its leakage can have lasted, and the measure qubit's edges
over the round it carried the excitation. With `dqlr_eff` = 0.98 almost all of the weight is on the round before the
herald. Nothing in the decoder is specific to DQLR: the same code produces the weights with the flag on or off.

A fixed-factor rule for this herald does not exist. The footprint rule re-weights only the heralded measure qubit's own
round, which is half of what the herald points at.

## Validation

`dqlr_herald` (tests/test_qpuforge.cpp), d = 5 test chip, 10 rounds, the data qubit as the leaking partner of every
CZ at 2e-3, leakage-removing reset and DQLR, 20,000 shots, 1.59 leakage events per shot:
- True heralds per shot: 1.3093 sampled with the transfer and 1.3189 from the model (+0.7%); 0.0887 and 0.0879
  without it (dump), where only leakage that reaches a measure qubit by transport is seen.
- Leaked fraction, data / measure qubits: 4.05e-3 / 2.63e-4 with the dump, 4.67e-3 / 5.90e-3 with the transfer.
- Failures: transfer 1541 without heralds, 988 with derived weights (−36%). Dump: 918 without heralds, 910 with.
- The look-back leaves less than 0.01 of a herald unexplained.
- The flag without DQLR changes nothing.

The thread-count determinism digest includes a run with the transfer.

## Results

`results/decoder_study2.txt`, section 2. Gamma-1, memory-Z, 14 rounds, exact matching, 100,000 shots per point,
derived weights. "Designed" is Gamma-1 as calibrated: the measure qubit leaks at the CZ. "Data leak" is a stress case:
the data qubit is the leaking partner of every CZ, at the same rates. Schemes are compared on different syndromes (the
circuits differ); within a scheme the two columns are the same syndromes.

Failures without heralds → with heralds (derived weights):

| d | case | leakage | LRR alone | LRR + DQLR (dump) | LRR + DQLR with herald |
|---|---|---|---|---|---|
| 5 | designed | ×1 | 230 → 224 | 267 → 257 | 272 → 257 |
| 5 | designed | ×10 | 594 → 418 | 632 → 405 | 676 → 434 |
| 5 | data leak | ×1 | 302 → 282 | 254 → 253 | 314 → 281 |
| 5 | data leak | ×10 | 1734 → 1239 | 458 → 458 | 1124 → 609 |
| 7 | designed | ×1 | 38 → 31 | 43 → 41 | 52 → 46 |
| 7 | designed | ×10 | 182 → 94 | 137 → 70 | 168 → 96 |
| 7 | data leak | ×1 | 83 → 70 | 45 → 45 | 52 → 39 |
| 7 | data leak | ×10 | 851 → 545 | 94 → 94 | 356 → 146 |

What the table shows:

- **The herald works as a herald.** With the transfer, about as many true heralds arrive per shot as the model predicts
  (d = 7, data leak, ×10: 1.0121 sampled, 1.0176 from the model), and decoding with them removes 36 to 59% of that
  scheme's failures at ×10 (d = 7, data leak: 356 → 146).
- **It does not pay.** In every row at ×10 the transfer with its herald is worse than DQLR that simply dumps the
  leakage and uses no data-qubit herald at all: 609 against 458 failures (d = 5, data leak), 146 against 94 (d = 7,
  data leak), 434 against 405 and 96 against 70 on the designed chip. Leaving the excitation on a measure qubit for a
  round costs that qubit's syndrome bit and puts random phases on its four data partners; the herald recovers about
  much of that damage, not all of it (without heralds the transfer has 2.5 to 3.8 times the dump's failures at ×10 in the
  data-leak case).
- **At the chip's own leakage nothing is resolved.** The ×1 rows differ by a few failures to a few tens, on different
  syndromes; the one row where the transfer is ahead (d = 7, data leak, ×1: 39 against 45) is within the counting error.
- **On Gamma-1 as designed, DQLR itself is marginal**, as the Beta-1 review found: at the chip's own leakage it costs more
  than it saves at both distances (the round is 30 ns longer and each pulse has an error of 5e-4); at ×10 it is level at
  d = 5 (405 against 418 failures with heralds) and helps at d = 7 (70 against 94). The data qubits here leak only by transport.
- **Where data qubits leak directly, DQLR matters far more than any herald**: at d = 7, ×10, the leakage-removing
  reset alone gives 851 failures (545 with heralds), the dump 94. The measure-qubit heralds still help the reset-only
  scheme (−36%), because a data qubit's leakage is seen once transport hands it to a measure qubit, and the chain
  then points back at the rounds it spent on the data qubit.

Conclusion: with this model of the transfer, heralding data-qubit leakage by leaving it on the measure qubit is not
worth doing. DQLR that dumps the leakage is the better use of the same pulse. The option stays in the code, off by
default, as the measured comparison.

## Limits

- The transfer is modelled as one pulse with the efficiency and the Pauli error of the existing DQLR (`dqlr_eff`,
  `dqlr_err`). The extra |11⟩ → |02⟩ step would have its own error and leakage; neither is included, so the transfer
  is, if anything, flattered.
- The second-readout variant, which would avoid leaving a measure qubit leaked for a round, is not modelled.
- The comparison between schemes is between different sets of syndromes, so differences of a few percent are within
  Monte Carlo noise. The conclusions above rest on differences much larger than that.
- The limits of `derived-herald-weights.md` apply: first order in the leakage, each herald on its own.
