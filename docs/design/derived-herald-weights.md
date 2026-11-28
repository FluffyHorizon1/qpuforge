# Design note: herald weights derived from the leakage model

Status: implemented on `scs-bot/decoding-next` (`HeraldSpec::derived`, `HeraldModel`, `frame_detail::leak_chain`,
`build_herald_model`, `prepare_herald_model` in `include/qpuforge/qec.hpp`). The fixed-factor rules of
`heralded-leakage.md` are unchanged and remain the default.

## Purpose

A herald says that a measure qubit read out as |2⟩. The decoder of `heralded-leakage.md` answers with fixed factors: the
edges of that qubit's round are multiplied by 0.4 and its own measurement edge is erased. The factors were found by a
sweep on one chip and one leakage rate. They do not know how likely a herald is to be false, which CZ the leakage most
likely started at, or what a leaked qubit does to each neighbour.

The frame simulator has a definite leakage model, so those probabilities can be computed. This note derives the edge
weights from it.

## The leakage model, as one excitation moving through the circuit

`frame_detail::leak_chain` follows a single leakage excitation with exactly the transitions of
`frame_detail::sample_shot`:

| Event | What the sampler does | In the chain |
|---|---|---|
| CZ leaks | with the pair's probability the moving qubit leaves the subspace; its partner is randomised | an excitation starts on the mover; the partner gets X and Z, each with probability 1/2 |
| single-qubit gate leaks | probability `leak_1q` | an excitation starts; no other effect |
| CZ with one leaked qubit | random Z on the partner; with probability `transport` the leakage moves to the partner and the first qubit returns randomised | Z on the partner with probability 1/2; with probability `transport` the excitation moves and the old qubit gets X and Z, each 1/2 |
| measurement of a leaked qubit | reads \|1⟩ | the bit is flipped if its noiseless value is 0, unchanged if it is 1, and flipped with probability 1/2 if it is not deterministic |
| reset | removes leakage if the reset is leakage-removing | the excitation ends, with no error |
| DQLR pulse | removes a data qubit's leakage with probability `dqlr_eff`; the qubit returns randomised | the excitation ends (or moves to the measure qubit, see `data-leakage-herald.md`); X and Z on the data qubit, each 1/2 |
| seepage | ends the leakage at 2/T1; the qubit returns randomised | the same, accumulated between two operations of the qubit |

The chain is first order in the leakage probabilities: it never has two leaked qubits at once. With that, the position
of the excitation is a Markov chain over (qubit, time), and three passes give everything a herald needs.

1. **Forward.** One pass over the program gives a(q, t), the probability that qubit q is leaked just before each of its
   operations. At a measure-qubit measurement c this is the prior Z_c = P(leaked at readout c).
2. **Backward, per readout.** b(q, t) is the probability that an excitation on q at time t is on the heralded qubit at
   its readout. The pass runs backwards from the readout through the same transitions. Smoothing then gives the
   probability of everything that happened before the readout, given that the qubit is leaked there:
   - the leakage started at a given CZ: p_leak · b / Z_c. This is the answer to "when did it most likely start";
   - qubit q was leaked during a given CZ: a(q, t) · b(q, t) / Z_c;
   - the leakage moved from q to its partner at a given CZ: a(q, t) · transport · b(partner) / Z_c.

   The pass goes back at most `max_rounds` rounds of the heralded qubit (8 by default) and stops earlier once the starts
   it has found add up to 1 − `cutoff` (1e-3). What is left is the probability that the leakage began before the
   look-back window. It is reported as `HeraldModel::worst_unexplained`.
3. **Forward, per readout.** From the heralded qubit at its readout the excitation is followed forward until it has
   ended with probability 1 − `cutoff`, again for at most `max_rounds` rounds. With the leakage-removing reset this is
   one operation: the reset ends it. Without it the excitation goes on to the next rounds.

How long leakage persists follows from the removal scheme through these same transitions. With the leakage-removing
reset a measure qubit's leakage ends at its reset, so it lasts less than one round. A data qubit's leakage is ended only
by transport to a measure qubit (5% per CZ, four CZs per round), by seepage or by DQLR, so without DQLR the backward
pass has to look several rounds back whenever the herald comes from leakage that moved in from a data qubit.

Every consequence is a single Pauli placed after a time slice, or a flipped readout bit, with its probability given
the leak. `fault_signatures` (the linear-time fault enumeration of `linear-fault-enumeration.md`, opened to any list of
mechanisms) gives the detectors and the logical bit that each one flips. A consequence that flips one or two detectors
is a decoder-graph edge. Consequences on the same edge are combined as independent faults,
p ← p(1 − q) + q(1 − p). The result is, per readout bit c, the list `HeraldModel::cond[c]` of
(edge, P(fault on that edge | qubit leaked at c)).

Edges whose conditional probability is below `min_prob` (5e-3) are left out. They would barely change a weight, and
every re-weighted edge adds two terminals to the exact decoder's per-shot graph (`decode_scaled`). At d = 5 the pruning
removes 4% of the conditional probability and 93% of the entries, and leaves the failure counts unchanged to within one
failure in 400 (measured on 100,000 shots while choosing the limit; not part of the committed study).

`prepare_memory` does not build the model. `prepare_herald_model(M)` does, in SETUP_TIME at d = 7 and 14 rounds.

## From a herald to weights

For a heralded bit c, with herald fidelity f and false-herald rate r:

    π_c = Z_c f / (Z_c f + (1 − Z_c) r)

is the probability that the qubit really was leaked. Every edge e in `cond[c]` gains the fault probability π_c · q_e.
For each edge touched by any herald of the shot, the decoder combines its own prior p_e with these (as independent
faults), giving p'_e, and decodes with the log-likelihood ratio

    w'_e = log((1 − p'_e) / p'_e).

Two cases need care.

- **p'_e > 1/2.** The fault is more likely than not, so the weight would be negative. The edge is applied instead: its
  two detectors are toggled, its logical bit enters the prediction, and it keeps the weight of the complement,
  log(p'_e / (1 − p'_e)). This is the usual treatment of a negative edge in matching. It is what happens to the
  measurement edge of a Z-type measure qubit when a herald is probably true: a leaked qubit reads |1⟩, the noiseless
  value is 0, so the bit is wrong with probability π_c, not 1/2. The fixed rule erases that edge (probability 1/2).
- **Weights are only lowered.** `decode_scaled` is exact only when edges get cheaper, so p'_e is not allowed to fall
  below p_e. This matters only after an edge was applied with π_c above 1 − p_e (about 0.99); the edge then keeps its
  original weight.

Exact matching uses `decode_scaled` on the re-weighted graph, and union-find takes the new lengths, as for the fixed
rule. The streaming decoder now takes them too (`StreamOptions::scaled`): each window that holds a re-weighted edge
decodes with its new length.

## What the derived weights look like

HERALD_EXAMPLE

## Validation

`herald_weights` (tests/test_qpuforge.cpp), d = 5 test chip, 10 rounds, measure-qubit CZ leakage 2e-3, leakage-removing
reset, 20,000 shots:
- The model's leak probability at readout against the sampler: 1.4639 true heralds per shot from the sum of the priors
  times the fidelity, 1.4628 sampled (+0.1%; limit 5%).
- The 8-round look-back leaves at most 0.006 of a herald unexplained, and no conditional probability falls on a
  consequence without a graph edge.
- Same syndromes: exact matching 928 failures without heralds, 694 with the fixed factors, 685 with the derived
  weights; union-find 984, 763, 732.
- Heralds that never fire give the failure count of a run without heralds, and building the model does not change the
  fixed-factor result.
- Streaming decoder with heralds: one window over the whole shot gives exactly the batch union-find counts (732 derived,
  763 fixed). With a 6-round window committing 3, 1020 failures without heralds and 760 with derived weights.

`dqlr_herald` checks the chain through the DQLR transfer (`data-leakage-herald.md`). The thread-count determinism
digest includes the model and a derived-weight run, batch and streaming.

## Results on Gamma-1

`results/decoder_study2.txt`, section 1: memory-Z, 14 rounds, leakage-removing reset, herald fidelity 0.95 and
false-herald rate 1e-3, CZ leakage as calibrated (×1) and ×10. Each row is one set of syndromes decoded three ways.
Failures, with the change against no heralds:

RESULTS_TABLE

RESULTS_TEXT

## Limits

- **First order.** The chain has one leaked qubit at a time. Two excitations that meet (a CZ between two leaked
  qubits) are not modelled; at ten times Gamma-1's leakage a measure qubit is leaked in 0.1% of its rounds, so such
  meetings are rare.
- **Each herald on its own.** A herald's posterior uses that herald only. A readout that did *not* herald is not used
  as evidence, and two heralds caused by the same excitation are counted as two independent ones. With the
  leakage-removing reset this is nearly exact, because a measure qubit's leakage does not survive its reset. Without
  it, a qubit that stays leaked for several rounds is heralded several times and its faults are over-weighted.
- **No leakage prior without a herald.** The detector error model has no leakage faults. A shot without heralds is
  decoded exactly as before, so the leakage that the readout misses (5% at fidelity 0.95) and all data-qubit leakage
  before it reaches a measure qubit stay invisible to the decoder.
- **Independent faults on a graph.** The consequences of one leak are correlated (the partners of one leaked measure
  qubit all get their Z errors in the same round). The matching graph takes each edge's marginal probability and
  drops the correlation. Consequences that flip more than two detectors would be dropped; none occur in the memory
  circuits tested.
- **The herald parameters are taken as known.** The fidelity and the false-herald rate enter π_c directly. A wrong
  false-herald rate changes every weight; the sensitivity has not been measured.
- **Window.** Leakage that began more than `max_rounds` rounds before the herald is not attributed. On Gamma-1 without
  DQLR this is at most UNEXPLAINED of a herald.
