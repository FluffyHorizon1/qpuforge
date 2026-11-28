# Design note: union-find decoder

Status: implemented on `scs-bot/gamma-decoding` (`include/qpuforge/unionfind.hpp`, `DecodeMethod::UnionFind`).

## Purpose

Phase Gamma needs real-time decoding at under 1 µs per round, first on an FPGA. The exact decoder (dynamic programme
or Edmonds' blossom on all-pairs shortest paths) stays the reference. Union-find is the usual hardware candidate: with
boundary lists it runs in almost-linear time, and it needs no shortest paths and does no matching. This implementation is the software model
of such a decoder, so its accuracy can be measured against the exact one on the same syndromes.

## Algorithm

The graph is the matching decoder's: detectors plus one boundary vertex. For each detector pair, the most likely
edge is kept, with weight w = log((1 − p)/p) and its logical-flip bit.

- **Weighted growth** (Huang, Newman and Brown 2020). Each edge has an integer length `round(16 w)`, at least 1.
  Every cluster with an odd number of defects that does not contain the boundary grows along all edges on its border:
  rate 1 from each odd cluster at an end, so 2 when both ends are growing. The step jumps straight to the next time
  an edge becomes fully grown, so the number of iterations does not depend on the length scale. Fully grown edges
  are merged in ascending edge order (union by size, path halving), and the parity and boundary flag of each merged
  cluster are updated. Growth stops when no odd cluster lacks the boundary.
- **Peeling** (Delfosse and Zémor). The fully grown edges form an erasure. A BFS spanning forest of them is built,
  rooted at the boundary first, then at any other reached vertex. It is peeled from the leaves inwards: a marked
  vertex adds its tree edge to the correction and moves its mark to its parent. The XOR of the correction's
  logical-flip bits is the prediction, and the correction itself is left in `UFScratch::correction`.
- **Reset.** Only the vertices and edges the previous call touched are reset, so a call costs time in proportion
  to the clusters it grows, not to the graph.

Because the boundary is a single vertex, every cluster that reaches it is neutral. This matches the matching
decoder's convention.

The same decoder accepts per-call length changes (`decode(defects, ws, &scaled)`). Heralded edges are shortened, and
an edge whose length rounds to 0 is an erasure: it is merged before growth starts. See `heralded-leakage.md`.

## Validation

`union_find` (tests/test_qpuforge.cpp), d = 3 and 5, both bases:
- 12,000 syndromes, sampled from the error model at 3× its probabilities or uniformly at random: every correction
  clears its syndrome.
- Every single-fault mechanism of the circuit that the graph represents (55,917): union-find decodes all of them
  correctly. Another 19 mechanisms are outweighed by a likelier path of the other parity; the exact decoder also
  gets them "wrong", and the test checks that their weight explains it.
- Logical error on the same syndromes (40,000 shots at d = 3, 6 rounds; 60,000 at d = 5, 5 rounds):
  ratio to exact 1.01 (d = 3, Z), 1.03 (d = 3, X), 1.08 (d = 5, Z), 1.08 (d = 5, X). The test requires < 1.5.

Gamma-1, 14 rounds, same syndromes (`results/decoder_study.txt`, section 2; 100,000 shots, 200,000 at d = 7):

| d | basis | exact ε per round (failures) | union-find ε (failures) | ratio | decode µs per round, exact / UF |
|---|---|---|---|---|---|
| 3 | Z | 8.77e-4 (1214) | 9.32e-4 (1289) | 1.06 | 0.11 / 0.06 |
| 3 | X | 8.11e-4 (1123) | 8.19e-4 (1135) | 1.01 | 0.03 / 0.06 |
| 5 | Z | 1.65e-4 (230) | 1.73e-4 (241) | 1.05 | 0.29 / 0.37 |
| 5 | X | 1.68e-4 (235) | 1.85e-4 (259) | 1.10 | 0.27 / 0.38 |
| 7 | Z | 3.65e-5 (102) | 3.57e-5 (100) | 0.98 | 0.87 / 0.78 |
| 7 | X | 3.50e-5 (98) | 3.65e-5 (102) | 1.04 | 0.88 / 0.75 |

At d = 7 the two are within Monte Carlo error of each other. Decode times are single-core software wall clock
(11.6 defects per shot at d = 7) and say little about hardware.

## Limits

- The software implementation is not optimised for speed. At a few defects per shot the exact dynamic programme is
  faster in software (test chip, d = 5: about 0.8 µs per shot against 1.6 µs); at Gamma-1 d = 7 union-find is
  already slightly faster (0.78 against 0.87 µs per round). The comparison that matters is the number of
  operations per defect on hardware, which this code does not model.
- There is no cluster-size bound and no sliding window: the decoder sees a whole shot. A streaming version is the
  next step towards the FPGA.
- Each growth step walks every member of each growing cluster, not a separate list of its boundary vertices. Below
  threshold clusters stay small and the cost is close to linear in the defects, but the worst case is quadratic in
  the cluster size. Keeping a boundary list per cluster gives the almost-linear bound of Delfosse and Nickerson.
