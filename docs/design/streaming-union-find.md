# Design note: streaming (sliding-window) union-find decoder

Status: implemented on `scs-bot/gamma-streaming` (`StreamingUnionFind` in `include/qpuforge/unionfind.hpp`).

## Purpose

A real-time decoder cannot wait for the end of an experiment. It sees the syndrome round by round, and its output for
the early rounds must be final while later rounds are still arriving. This decoder is the software reference for such
a scheme, built on the batch union-find decoder (`union-find-decoder.md`). It lets the cost of windowing be measured
on the same syndromes as the batch decoder.

## Scheme

Parameters: a window of W rounds and a commit region of C rounds, with 1 ≤ C < W. A detector's round is
`MemoryExperiment::detector_round`; the final data-readout detectors form the last round.

1. A window covers rounds [t, t + W). Its graph is the detector graph restricted to those rounds. An edge from the
   window into a later round goes to the boundary vertex: the future acts as a boundary, so a defect near the window's
   end can be absorbed there for now. Edges into rounds that are already committed are dropped.
2. The carried syndrome inside the window is decoded with weighted union-find.
3. Every correction edge that touches the first C rounds of the window is committed. Its logical-flip bit enters the
   prediction, and both of its ends are toggled in the carried syndrome. That clears the commit region, and an edge
   reaching into the buffer moves its defect there. The other correction edges are discarded.
4. The window advances by C rounds. The last window (the one that reaches the final round) commits everything.

Detector-graph edges span at most one round, so with W > C a committed edge never reaches beyond the window. A future
boundary edge always starts in the buffer, so it is never committed. When a defect in the buffer was merged with a
"real boundary" edge, the merge makes no difference, because that edge would only be committed from the commit
region. The union of the committed edges clears the whole syndrome: a commit-region vertex has an odd number of
incident correction edges in its window, and all of them are committed.

Each window has a fixed graph, built once when the decoder is constructed (about R/C windows). Each window also keeps
its own union-find scratch: the bulk windows have graphs of the same size but different edge lengths, and a shared
scratch would keep stale lengths. This was found by the thread-count determinism test. `UFScratch` now also records
which decoder it was set up for, and resets itself when it is used with another.

Latency is measured per window, as wall-clock time on one core. Dividing by C gives the time per committed round,
which is what has to stay below the round time.

## Validation

`streaming_uf` (tests/test_qpuforge.cpp):
- d = 3 and 5, 9 rounds, windows 2/1, 4/2 and 6/3: 6,000 syndromes, sampled from the error model at 3× its
  probabilities or uniformly at random. Every committed correction clears its syndrome. All 158,250 single faults that
  batch union-find decodes correctly are decoded correctly.
- A single window over the whole shot is exactly the batch decoder: 815 failures in both, on the same 30,000 shots
  (d = 5, 15 rounds).
- Logical error against batch union-find, same syndromes (d = 5, 15 rounds, 30,000 shots, 815 batch failures):

  | window / commit | windows | failures | × batch | µs per window (mean) | µs per committed round |
  |---|---|---|---|---|---|
  | 2 / 1 | 15 | 1835 | 2.25 | 0.70 | 0.70 |
  | 3 / 1 | 14 | 1102 | 1.35 | 1.06 | 1.06 |
  | 4 / 2 | 7 | 950 | 1.17 | 1.70 | 0.85 |
  | 6 / 3 | 5 | 855 | 1.05 | 2.40 | 0.80 |
  | 10 / 5 | 3 | 819 | 1.00 | 4.19 | 0.84 |

  The test requires less than 1.25× for windows of at least 4 rounds with C ≤ W/2 at this distance; at d = 7 the
  buffer must grow with the distance (below). A one-round buffer (2/1) is not
  enough: defects at the window's end escape to the future boundary instead of pairing across it.
- The thread-count determinism digest includes a streaming run.

Gamma-1, d = 7, 50 rounds, memory-Z, 50,000 shots, same syndromes (`results/decoder_study.txt`, section 5). Exact
matching gives 95 failures and batch union-find 94:

| window / commit | failures | × batch | µs per committed round |
|---|---|---|---|
| 2 / 1 | 2217 | 23.6 | 1.78 |
| 4 / 2 | 356 | 3.79 | 1.72 |
| 6 / 3 | 168 | 1.79 | 1.53 |
| 10 / 5 | 109 | 1.16 | 1.73 |
| 20 / 10 | 99 | 1.05 | 1.84 |

The window that d = 5 tolerates is far too short at d = 7. The rounds beyond the commit region (W − C) must be
comparable to the code distance: an error chain can stretch over about d rounds before its two ends meet. With
W − C ≥ d (10/5 and 20/10 here) the cost is 5–16%. On one core in software, the time per committed round (1.5–1.8 µs)
is about twice the 0.73 µs round time, so this implementation would not keep up as it stands.

## Limits

- Latency is software wall clock. The maximum per window is dominated by preemption and first-touch allocation
  (milliseconds against microseconds on average), so only the mean is used. A hardware estimate needs an operation
  count per window.
- Windows are decoded one after another. A pipelined or parallel-window scheme, in which independent windows are
  decoded concurrently and their seams are fixed afterwards, is not implemented.
- The decoder works on whole shots from the frame simulator. It does not take rounds incrementally from a stream,
  although the window order is the one a stream would deliver.
