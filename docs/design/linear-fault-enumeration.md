# Design note: linear-time single-fault enumeration

Status: implemented on `scs-bot/gamma-decoding` (`enumerate_faults` in `include/qpuforge/qec.hpp`).

## Problem

The detector error model is built from the signature of every single fault: which detectors it flips and whether
it flips the logical. In 0.4, `enumerate_faults` propagated each fault through the rest of the program, following
only the qubits that carry a Pauli. A Pauli on a data qubit survives every later round, though, so each fault still
costs time proportional to the rounds after it. With a fixed number of faults per round, setup is quadratic in the
number of rounds. The Gamma-1 review shortened its long run in smoke mode for this reason.

## Method

Pauli-frame propagation is linear over GF(2). The bits a frame flips from some point onwards are the XOR of the bits
flipped by each of its single-qubit X and Z components. Detectors and the observable are XORs of bits, so the same
holds for them.

1. **Cuts.** The program is cut after every time slice that contains a measurement (roughly one or a few per round).
2. **Tails.** For each cut `c`, qubit `q` and Pauli `P` in {X, Z}, `tail[c][q][P]` is the set of detectors (and the
   observable bit) flipped by `P` on `q`, placed just after cut `c` and propagated to the end. Tails are built
   backwards. Propagate the single Pauli to the next cut, fold the bits it flipped into detector parities, then XOR
   in the tails of the frame left at that cut. Each step spans one segment, so building all tails costs
   O(cuts × qubits × segment). The detector sets stay small because a persistent data error flips the same
   stabiliser in every round, and consecutive-round detectors cancel in pairs.
3. **Faults.** Each fault is propagated only up to the first cut at or after it. The bits it flips are folded in, and
   the remaining frame is finished with tails. A fault placed right after a cut goes straight to tails.

Both phases reuse the sparse propagator of the 0.4 code (`frame_detail::SparseFrame`): only operations on qubits
that carry a Pauli are applied, and time slices without any are skipped. Workspaces are allocated once per thread
for all levels; per-level allocation of the O(rounds) bit and detector arrays would bring the quadratic cost back.
Tail levels run in sequence, with each level's qubits split across threads. Faults run in parallel and write to
their own index, so the output does not depend on the thread count.

The result is exact, fault by fault and in the same order as `enumerate_faults_reference`. The cut placement only
affects speed, not the result.

## Validation

`fault_enumeration` (tests/test_qpuforge.cpp) compares the reference, the 0.4 sparse version (now
`enumerate_faults_direct`) and the linear-time version fault by fault on ten experiments: d = 3 and 5, both bases,
1 to 8 rounds, with DQLR faults, and with unequal CZ durations (twice the time slices). Result: 188,932 faults,
0 differences. Time per fault on d = 5: 0.121 µs at 10 rounds and 0.180 µs at 40 rounds (best of three). The test
allows a factor of 2.5; the quadratic version takes about 4× longer per fault at 40 rounds.

Gamma-1 (d = 7, memory-Z; `results/decoder_study.txt`, section 1, four threads):

| rounds | faults | time slices | sparse (0.4) | linear-time | speed-up | signatures differing |
|---|---|---|---|---|---|---|
| 25 | 2,850,523 | 10,019 | 7.48 s | 0.66 s | 11× | 0 |
| 50 | 5,967,552 | 21,044 | 41.37 s | 1.30 s | 32× | 0 |

Doubling the rounds doubles the linear-time cost (0.66 → 1.30 s) and multiplies the sparse cost by 5.5. The whole
setup (enumeration, detector error model, matching decoder with its all-pairs shortest paths, herald map) takes
1.33 s at 25 rounds and 3.76 s at 50. Enumeration is no longer the larger part; the rest includes the decoder's
all-pairs shortest paths, which grow as the cube of the detector count.

## API

- `enumerate_faults(F, E)`: linear-time (used by `prepare_memory`).
- `enumerate_faults_direct(F, E)`: the 0.4 sparse version, quadratic in rounds; kept for timing comparisons.
- `enumerate_faults_reference(F, E)`: dense propagation; the test reference.
