# Design note: an operation-count model of the streaming union-find decoder for an FPGA

Status: a model, merged into `main` in PR #35. The software decoder is instrumented (`UFOpCounts`,
`UnionFindDecoder::decode_counted`, `StreamOptions::count`), its operations are counted on simulated syndromes
(`run_stream_profile`), and two stated hardware models turn the counts into clock cycles (`UFCostModel`).
**Nothing here was synthesised, placed, routed or run on hardware.** No clock frequency, resource count or timing in
this note is a measurement of an FPGA.

## Question

`streaming-union-find.md` ends with a software latency: 1.5 to 1.8 µs per committed round on one core at d = 7,
against Gamma-1's 0.73 µs round. Wall-clock time of a C++ program says little about an FPGA. What carries over is the
work the algorithm does: how many defects it loads, edges it touches, pointers it follows, clusters it merges. This
note counts that work, states a cost per operation, and asks whether the decoder keeps up with the rounds and, if not,
what would make it.

## Counters

`UnionFindDecoder::decode_counted` and `StreamingUnionFind::decode` with `StreamOptions::count` fill a `UFOpCounts` per
window. The counted decoder is the same template as the plain one, instantiated with a counter object; the plain
`decode()` gets a counter whose calls are empty, so it does no extra work. The test `stream_ops` checks that counted
and plain decoding give the same correction on 3,000 syndromes.

| Counter | What is counted |
|---|---|
| `loads` | defects in the window when its decode starts |
| `growth_iters` | growth steps: each jumps to the next time an edge becomes fully grown |
| `edge_visits` | (vertex, incident edge) pairs examined while collecting the border edges of the growing clusters |
| `edge_updates` | border edges whose growth is advanced, one per edge and growth step |
| `finds` | root look-ups |
| `hops` | parent-pointer dereferences made by those look-ups (path halving) |
| `unions` | merges of two different clusters |
| `peel_steps` | edges examined by the spanning-forest search, plus vertices peeled |
| `corrections` | correction edges produced by the window |
| `commits` | correction edges committed |
| `pre_pairs` | isolated defect pairs removed by the pre-decoder (below) |

The counts are those of this implementation. It keeps no border list: a growth step walks every vertex of each growing
cluster. `edge_visits` therefore grows with cluster size times the number of steps. On the measured syndromes
`edge_updates` is 89% of `edge_visits` (d = 7, window 10 / commit 5), so a border list would save little here; the cost is in the large
clusters themselves.

Everything is an integer count per window. `run_stream_profile` collects histograms of them, so its result does not
depend on the thread count (it is in the determinism digest).

## The two hardware models

Both are stated in `UFCostModel::serial_bram()` and `UFCostModel::parallel_growth()`, which the study uses. Common
assumptions:

- **Clock: 250 MHz** (4 ns). Gamma-1's round is 730.4 ns in the 50-round runs at d = 7, so a round is 182.6 cycles.
- A window with no defects costs 1 cycle (an OR over the window's syndrome register).
- The window graph is fixed logic or ROM. Edge lengths are the software's integers (at most about 120: 7 bits).
- Resetting state between windows is free (generation tags on the table entries, or a parallel clear of registers).

**S: serial, block RAM.** One finite-state machine does one thing at a time, as the software does.
- Memory: a vertex table (parent, size, parity, boundary flag, member-list links) and an edge table (growth, length) in
  true dual-port block RAM. One read and one write per table per cycle.
- Pipeline: reads are registered at the address and at the output, so data arrives 2 cycles after the address; the
  pipeline is 4 stages deep (address, RAM, RAM register, compute and write back).
- Independent accesses stream through the pipeline at one per cycle. A dependent access, where the next address is
  the data just read, waits the 2-cycle read.

| Operation | Cycles | Reason |
|---|---|---|
| defect loaded | 1 | a priority encoder over the syndrome register gives one defect per cycle; one table write |
| edge visit | 1 | one streamed read of the edge table |
| find | 1 | one streamed read of the vertex table |
| hop | 2 | a dependent read |
| edge update | 2 | one pass to find the step length (the minimum), one to add it |
| growth step | 8 | the pipeline drains twice: before the minimum is known and before the merges |
| union | 2 | two table writes (root record, member list) |
| peeling step | 1 | one streamed access |
| correction edge | 1 | commit test and syndrome toggle |
| fixed, per window with defects | 12 | three phases (load, growth, peeling), each filling the pipeline |

**P: parallel growth.** One cell per vertex (root register, parity, boundary bit) and one per edge (growth counter,
length), all in registers. The degree of parallelism is the window graph: at d = 7 a 10-round window has 240
vertex cells and 1274 edge cells, each edge cell with a 7-bit adder and comparator.

| Operation | Cycles | Reason |
|---|---|---|
| growth step | 5 | every edge cell flags whether it is on the border of a growing cluster (1); the smallest remaining length over the flagged edges, in a pipelined comparator tree (2); all flagged cells add (1); fully grown edges are found (1). The cost does not depend on the cluster sizes |
| union | 2 | the absorbed root is broadcast; every vertex cell compares and rewrites its root (2, pipelined). Roots stay flat, so a find is a register read |
| find, hop, edge visit, edge update, load, correction | 0 | done by all cells at once inside the steps above, or a register read |
| peeling step | 1 | peeling stays serial |
| fixed, per window with defects | 4 | parallel load, start, parallel commit |

Model P assumes that a comparator tree over about a thousand 7-bit values and a root broadcast to a few hundred cells
close timing at 250 MHz with the stated pipelining. That is plausible for a current FPGA and unverified here. It is the
assumption most likely to be wrong, and the one a first synthesis should test.

**Pre-decoder.** Before union-find runs on a window, an isolated defect pair is removed: two defects joined by an
edge, neither with another defect on any neighbouring vertex, and the edge no longer than the two defects' boundary
edges together. The pair's edge is its correction. Each vertex can decide this from its neighbours' syndrome bits, so
it is one layer of combinational logic over the window; the model charges 2 cycles when it clears the window, and
nothing extra otherwise. It is `StreamOptions::predecode` in software, so its effect on the logical error is measured,
not assumed. It is not what union-find would always do: an isolated pair near another cluster can be matched
differently.

**Latency.** Windows of a shot are decoded one after another. A window can start when its last round has arrived and
the previous window is done; its latency is the time from that arrival to the end of its decode. The shot starts with
an idle decoder. This is a single-server queue, so the mean cost has to stay below the round for the decoder to keep
up at all, and the tail decides how far it falls behind.

## Counts

`results/decoder_study2.txt`, section 3: Gamma-1, memory-Z, 50 rounds, 50,000 shots, the window / commit settings
of `decoder_study` section 5 and the same seed. Per window: mean / 99th percentile / 99.9th percentile / maximum.

d = 7, window 10 / commit 5 (500,000 windows):

| Counter | mean | 99th | 99.9th | max |
|---|---|---|---|---|
| defects loaded | 7.21 | 17 | 21 | 31 |
| growth steps | 3.63 | 12 | 20 | 61 |
| edge visits | 196.30 | 1052 | 3313 | 35356 |
| edge updates | 175.14 | 820 | 1959 | 11719 |
| finds | 57.54 | 271 | 1119 | 17455 |
| pointer hops | 14.02 | 152 | 919 | 16774 |
| unions | 5.68 | 20 | 36 | 118 |
| peeling steps | 21.15 | 67 | 117 | 388 |
| correction edges | 4.03 | 10 | 12 | 21 |
| commits | 2.14 | 6 | 8 | 14 |

Across distances and settings (per window; edge visits as mean / 99.9th / max):

| d | window / commit | vertices, edges | windows without defects | defects (mean) | growth steps (mean) | edge visits | unions (mean) |
|---|---|---|---|---|---|---|---|
| 3 | 2 / 1 | 8, 24 | 87.0% | 0.19 | 0.16 | 1.3 / 42 / 203 | 0.22 |
| 3 | 4 / 2 | 16, 58 | 76.0% | 0.42 | 0.33 | 4.3 / 103 / 329 | 0.43 |
| 3 | 6 / 3 | 24, 92 | 66.0% | 0.66 | 0.50 | 7.7 / 153 / 495 | 0.65 |
| 3 | 10 / 5 | 40, 160 | 51.7% | 1.08 | 0.80 | 14.2 / 213 / 729 | 1.05 |
| 3 | 20 / 10 | 80, 330 | 29.4% | 2.08 | 1.46 | 32.1 / 342 / 1087 | 2.00 |
| 5 | 2 / 1 | 24, 89 | 64.3% | 0.64 | 0.57 | 9.8 / 351 / 1403 | 0.76 |
| 5 | 4 / 2 | 48, 213 | 43.5% | 1.36 | 0.93 | 21.2 / 712 / 4391 | 1.24 |
| 5 | 6 / 3 | 72, 337 | 28.3% | 2.13 | 1.38 | 38.2 / 1058 / 6145 | 1.85 |
| 5 | 10 / 5 | 120, 585 | 13.9% | 3.49 | 2.12 | 73.3 / 1378 / 10860 | 2.95 |
| 5 | 20 / 10 | 240, 1205 | 3.2% | 6.74 | 3.63 | 175.1 / 2087 / 15214 | 5.58 |
| 7 | 2 / 1 | 48, 194 | 41.0% | 1.32 | 1.12 | 25.0 / 711 / 5482 | 1.55 |
| 7 | 4 / 2 | 96, 464 | 18.8% | 2.82 | 1.74 | 56.6 / 1712 / 16544 | 2.44 |
| 7 | 6 / 3 | 144, 734 | 7.9% | 4.41 | 2.48 | 101.5 / 2222 / 24641 | 3.59 |
| 7 | 10 / 5 | 240, 1274 | 2.2% | 7.21 | 3.63 | 196.3 / 3313 / 35356 | 5.68 |
| 7 | 20 / 10 | 480, 2624 | 0.2% | 13.91 | 5.81 | 465.5 / 5712 / 38625 | 10.68 |

The means are small and the tails are long. At d = 7 with a 10-round window, a window holds 7.2 defects on average
and its decode visits 196 edges; one window in a thousand visits more than 3313, and the largest seen visited
35356, 180 times the mean. The tail is a few large clusters: each growth step of a large cluster walks all of
it. Border edges are 89% of the visited edges on average (175 updates for 196 visits), which is why a border
list would not change the picture.

## Cycles per committed round against the round

A round is 182.6 cycles at d = 7 (182.0 at d = 5, 181.5 at d = 3). A window's budget is its committed rounds times that. Latency is for the
continuous stream: 1,024 shots (52,224 rounds) end to end, the queue carried from shot to shot.

**Model S (serial, block RAM):**

| d | window / commit | cycles per committed round: mean (share of the round) | 99th | 99.9th | max | windows over their budget | latency, rounds: mean | 99.9th | max |
|---|---|---|---|---|---|---|---|---|---|
| 3 | 2 / 1 | 10.1 (6%) | 98 | 210 | 748 | 0.220% | 0.06 | 1.2 | 4.4 |
| 3 | 4 / 2 | 12.1 (7%) | 105 | 202 | 570 | 0.200% | 0.14 | 2.3 | 8.7 |
| 3 | 6 / 3 | 12.9 (7%) | 106 | 192 | 573 | 0.127% | 0.22 | 3.3 | 11.5 |
| 3 | 10 / 5 | 13.6 (7%) | 87 | 158 | 493 | 0.050% | 0.39 | 4.6 | 21.5 |
| 3 | 20 / 10 | 14.1 (8%) | 74 | 124 | 370 | 0.010% | 0.85 | 7.7 | 23.9 |
| 5 | 2 / 1 | 46.7 (26%) | 743 | 1146 | 3668 | 4.365% | 2.79 | 125.3 | 230.3 |
| 5 | 4 / 2 | 46.9 (26%) | 338 | 1163 | 6864 | 3.558% | 0.88 | 35.2 | 103.0 |
| 5 | 6 / 3 | 51.7 (28%) | 343 | 1096 | 6726 | 3.522% | 1.38 | 55.2 | 146.0 |
| 5 | 10 / 5 | 57.3 (31%) | 349 | 880 | 7098 | 3.816% | 2.41 | 85.0 | 260.6 |
| 5 | 20 / 10 | 64.2 (35%) | 319 | 667 | 4919 | 4.037% | 5.17 | 125.0 | 535.3 |
| 7 | 2 / 1 | 106.8 (58%) | 1146 | 2019 | 8542 | 13.504% | 20.28 | 363.4 | 679.5 |
| 7 | 4 / 2 | 112.9 (62%) | 739 | 2714 | 26128 | 15.738% | 9.56 | 333.8 | 596.9 |
| 7 | 6 / 3 | 125.4 (69%) | 715 | 2289 | 26054 | 18.479% | 12.91 | 474.8 | 803.6 |
| 7 | 10 / 5 | 142.1 (78%) | 691 | 2083 | 22157 | 23.118% | 27.31 | 874.3 | 1204.6 |
| 7 | 20 / 10 | 161.2 (88%) | 634 | 1790 | 12163 | 30.396% | 61.95 | 1130.0 | 1337.5 |

**Model P (parallel growth):**

| d | window / commit | cycles per committed round: mean (share of the round) | 99th | 99.9th | max | windows over their budget | latency, rounds: mean | 99.9th | max |
|---|---|---|---|---|---|---|---|---|---|
| 3 | 2 / 1 | 3.4 (2%) | 25 | 40 | 85 | 0.000% | 0.02 | 0.2 | 0.5 |
| 3 | 4 / 2 | 2.9 (2%) | 20 | 28 | 50 | 0.000% | 0.03 | 0.3 | 0.6 |
| 3 | 6 / 3 | 2.7 (1%) | 15 | 22 | 41 | 0.000% | 0.05 | 0.4 | 0.7 |
| 3 | 10 / 5 | 2.5 (1%) | 11 | 16 | 31 | 0.000% | 0.07 | 0.4 | 0.9 |
| 3 | 20 / 10 | 2.2 (1%) | 8 | 10 | 23 | 0.000% | 0.12 | 0.6 | 1.2 |
| 5 | 2 / 1 | 9.2 (5%) | 80 | 100 | 179 | 0.000% | 0.05 | 0.6 | 1.2 |
| 5 | 4 / 2 | 7.3 (4%) | 33 | 67 | 153 | 0.000% | 0.08 | 0.7 | 1.7 |
| 5 | 6 / 3 | 6.8 (4%) | 28 | 51 | 117 | 0.000% | 0.11 | 0.9 | 1.9 |
| 5 | 10 / 5 | 6.2 (3%) | 23 | 37 | 99 | 0.000% | 0.17 | 1.0 | 2.7 |
| 5 | 20 / 10 | 5.3 (3%) | 16 | 25 | 64 | 0.000% | 0.30 | 1.4 | 4.1 |
| 7 | 2 / 1 | 16.9 (9%) | 96 | 130 | 245 | 0.003% | 0.09 | 0.8 | 2.4 |
| 7 | 4 / 2 | 13.0 (7%) | 52 | 102 | 303 | 0.017% | 0.14 | 1.1 | 3.5 |
| 7 | 6 / 3 | 11.9 (7%) | 42 | 76 | 240 | 0.002% | 0.20 | 1.3 | 3.9 |
| 7 | 10 / 5 | 10.8 (6%) | 34 | 58 | 173 | 0.000% | 0.30 | 1.6 | 4.7 |
| 7 | 20 / 10 | 9.3 (5%) | 24 | 40 | 103 | 0.000% | 0.52 | 2.4 | 7.6 |

**Serial model.** On the mean it fits: 142 cycles per committed round at d = 7 with a 10-round window, 78% of the
round, and 88% with a 20-round window. The mean is not what limits a real-time decoder. 23% of the 10-round windows take
longer than the five rounds they commit, one in a thousand needs 2083 cycles per committed round (11 rounds'
worth), and the largest needed 22157 (121 rounds' worth for every round it commits). With a mean load of 78% the queue
drains slowly after such a window: over the continuous stream the mean latency is 27 rounds, the 99.9th percentile
874 rounds and the maximum 1205 rounds (0.88 ms). The last quarter of the stream looks like the whole
(50 against 62 rounds of mean latency for 20 / 10, 28 against 27 for 10 / 5), so the queue is
stable, not growing without bound; it is just hundreds of rounds deep at its worst. At d = 5 the serial model has
a mean load of 26 to 35% of the round and still a 99.9th-percentile latency of 35 to 125 rounds (maximum 535). At d = 3
the mean load is 6 to 8%, at most 0.22% of windows exceed their budget, and the latency is 1 to 8 rounds at the 99.9th
percentile and 24 rounds at the maximum.

**Parallel-growth model.** 10.8 cycles per committed round on average at d = 7, 10 / 5 (6% of the round),
58 at the 99.9th percentile and 173 at the maximum, which is still below the 183 cycles of one round. No 10-round
or 20-round window exceeded its budget; of the shorter windows, 0.002 to 0.017% did (maximum 303 cycles per committed
round at 4 / 2). Latency over the continuous stream: mean 0.30 rounds, 99.9th percentile
1.6, maximum 4.7 rounds. Short windows are the worst case per committed round for this model (2 / 1: maximum
245 cycles), because a window's fixed and growth costs are spread over fewer rounds; long windows have the larger absolute
latency (7.6 rounds at most for 20 / 10), because a window waits for its last round before it starts.

## What closes the gap

Against the 182.6 cycles of a round at d = 7, the serial model's gap is in the tail, not in the mean. Three levers,
with their measured effect:

1. **Parallel cluster growth within a window (model P).** This closes it. It removes the term that carries the tail:
   the per-step walk over every edge of a large cluster. At d = 7, 10 / 5, the 99.9th percentile falls from
   2083 to 57 cycles per committed round and the maximum from 22157 to 173; the worst latency over the stream falls from
   1204 rounds to 4.7. The price is hardware in proportion to the window graph (240 vertex cells and 1,274 edge cells
   for a 10-round window at d = 7, 480 and 2,624 for 20 rounds) and the two timing assumptions stated above.

2. **A pre-decoder for isolated defect pairs.** It clears a useful fraction of windows and does not close the gap.
   At d = 7:

   | window / commit | windows without defects | cleared by the pre-decoder alone (all windows) | of the windows with defects | failures of 50,000: without → with pre-decoder | S: mean cycles per round, without → with | S: 99.9th, without → with |
   |---|---|---|---|---|---|---|
   | 2 / 1 | 41.0% | 23.4% | 39.7% | 2189 → 2180 | 107 → 72 | 2019 → 2038 |
   | 4 / 2 | 18.8% | 37.0% | 45.5% | 284 → 317 | 113 → 63 | 2714 → 3185 |
   | 6 / 3 | 7.9% | 37.8% | 41.1% | 125 → 159 | 125 → 67 | 2289 → 2603 |
   | 10 / 5 | 2.2% | 30.5% | 31.1% | 85 → 115 | 142 → 72 | 2083 → 2346 |
   | 20 / 10 | 0.2% | 15.2% | 15.2% | 80 → 103 | 161 → 77 | 1790 → 2489 |

   It halves the serial model's mean (142 → 72 cycles per committed round at 10 / 5). It leaves the tail where it was or
   worse (99.9th percentile 2083 → 2346), because the windows in the tail are the ones with large clusters, which no
   pair rule clears. And at d = 7 it costs logical errors: 85 → 115 failures at 10 / 5 and 80 → 103 at 20 / 10, about
   a third more, on the same syndromes. At d = 3 the failure counts are unchanged (2043 at 10 / 5 either way) and at d = 5
   they move by a few (400 → 416 at 10 / 5). The fraction of windows it clears falls as windows get longer, since a
   long window rarely holds nothing but isolated pairs: 54% of the windows with defects at d = 3 (10 / 5), 31% at
   d = 7 (10 / 5), 15% at d = 7 (20 / 10). As implemented here it is not worth its logical cost at d = 7. A stricter
   rule (for example, pairs only on the likeliest edges, or only away from the window's end) might keep the accuracy;
   that has not been tried.

3. **Window choice and clock.** Neither helps the serial model enough. Its mean per committed round rises slowly with
   the window (107 cycles at 2 / 1, 161 at 20 / 10) while the tail stays between 1790 and 2715 cycles at the 99.9th
   percentile, ten to fifteen rounds' worth. A clock twice as fast would halve every number and still leave the
   99.9th percentile five to seven rounds' worth per committed round.

So the answer of this model is: a serial, memory-bound union-find decoder at 250 MHz keeps up with Gamma-1 on average
and falls hundreds of rounds behind in its tail; growing all clusters of a window in parallel brings the worst
window seen to under one round's worth of cycles per committed round for windows of 10 rounds or more (at most 1.7
rounds' worth for the shorter ones), with a latency of at most 4.7 rounds at 10 / 5 and 7.6 at 20 / 10. The pre-decoder
is a throughput optimisation for a serial design, not a fix for the tail, and it is not free in logical error.

What would have to be checked before any of this is believed for hardware: the two timing assumptions of model P by
synthesis, the tail under bursts and stressed leakage, and the logical error of the window scheme itself, which
`streaming-union-find.md` measures (a 10-round window committing 5 gives 85 failures here, a 20-round window 80, on
50,000 shots).

## Limits

- **A model.** The cycle costs are assumptions about a design that does not exist. The operation counts are exact for
  this software implementation on these simulated syndromes; a hardware design would not do the same operations in the
  same order (model P already does not).
- **One clock, no resource estimate.** Model P's cell counts are given, but not LUTs, registers or routing. Whether its
  comparator tree and root broadcast meet 250 MHz is not checked.
- **Peeling is serial in both models**, and the window graph is taken as fixed. A parallel peeling scheme would lower
  model P's cost further.
- **No input or output.** Getting the syndrome from the readout chain into the decoder and the correction out to the
  control system is not in the budget. Neither is the readout and discrimination latency.
- **Windows one after another.** Decoding several windows at once, with their seams fixed afterwards, is not modelled
  (it is also not implemented in software).
- **Syndromes of a quiet chip.** The counts are for Gamma-1's noise without impacts and without stressed leakage. A
  burst produces many defects at once; the tail under bursts has not been measured.
- **Memory experiments only.** Lattice-surgery windows have a different graph.
- **Finite samples.** The 99.9th percentile rests on the top 0.1% of the windows of 50,000 shots; the maximum is
  the largest value seen, not a bound.
