# Design note: an operation-count model of the streaming union-find decoder for an FPGA

Status: a model, on `scs-bot/decoding-next`. The software decoder is instrumented (`UFOpCounts`,
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
`edge_updates` is about 90% of `edge_visits`, so a border list would save little here; the cost is in the large
clusters themselves.

Everything is an integer count per window. `run_stream_profile` collects histograms of them, so its result does not
depend on the thread count (it is in the determinism digest).

## The two hardware models

Both are stated in `UFCostModel::serial_bram()` and `UFCostModel::parallel_growth()`, which the study uses. Common
assumptions:

- **Clock: 250 MHz** (4 ns). Gamma-1's round is ROUND_NS ns in the 50-round runs, so a round is ROUND_CYCLES cycles.
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
length), all in registers. The degree of parallelism is the window graph: at d = 7 a 10-round window has P_VERTICES
vertex cells and P_EDGES edge cells, each edge cell with a 7-bit adder and comparator.

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

`results/decoder_study2.txt`, section 3: Gamma-1, memory-Z, 50 rounds, COUNT_SHOTS shots, the window / commit settings
of `decoder_study` section 5 and the same seed. Per window: mean / 99th percentile / 99.9th percentile / maximum.

COUNTS_TABLE

COUNTS_TEXT

## Cycles per committed round against the round

CYCLES_TABLE

CYCLES_TEXT

## What closes the gap

GAP_TEXT

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
- **Finite samples.** The 99.9th percentile rests on the top 0.1% of the windows of COUNT_SHOTS shots; the maximum is
  the largest value seen, not a bound.
