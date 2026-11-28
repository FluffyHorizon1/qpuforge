# Design note: lattice surgery

Status: implemented in `include/qpuforge/surgery.hpp` and merged into `main` in PR #21 (Phase Gamma, ROADMAP.md).
Joint and correlated decoding, routing strips and the chip account were added in PR #30.
Study: `examples/surgery_study.cpp`, results in `results/surgery_study.{txt,json}`.

## Goal

Phase Gamma needs logical operations between surface-code patches: a logical ZZ and XX measurement between two
patches, and from those a logical CNOT. qpuforge 0.4 could only run a single patch (`memory_experiment`). This note
describes the multi-patch layer: geometry, merge and split schedules, detectors and observables across the merge and
split boundaries, decoding, and how it is verified.

## Geometry

Coordinates are those of `surface.hpp`: data qubit (i, j) at (2i+1, 2j+1), face (A, B) at (2A, 2B). A *region* is a
w × h rectangle of data qubits with an even origin and odd w and h. Its faces follow the rotated-code rules of
`rotated_surface_code`: weight-2 X faces on the top and bottom boundaries, weight-2 Z faces on the left and right, X
where A + B is even, and the same CZ order (X faces NW NE SW SE, Z faces NW SW NE SE). Because the origin is even, the
type of a face depends only on its position, so a face keeps its type and its measure qubit whichever region it
belongs to. Only its support changes (weight 2 on a patch boundary, weight 4 inside a merged region).

A *patch* is a d × d region. Two patches one gap row apart (A at (0, 0), B at (0, d+1)) merge into a d × (2d+1) region.
The faces in the two seam rows that neither patch has are all Z type, and their product is Z on A's bottom row times Z
on B's top row, which is Z_A Z_B: every gap qubit lies in exactly two of them. The X faces of the seam rows are the
patches' own boundary faces, extended from weight 2 to weight 4 by two gap qubits. Two patches one gap column apart
merge into a (2d+1) × d region in the same way with the roles of X and Z exchanged. The lattice (`SurgeryLattice`) is
the union of every region the schedule uses; qubits are numbered row by row, as in `surface.hpp`, so a single patch
gets the numbering, the couplers and the chip of `rotated_surface_code`.

| Merge | Gap qubits prepared in | New seam faces | Product of their first outcomes | Gap read out at the split in |
|---|---|---|---|---|
| ZZ (patches stacked vertically) | \|+⟩ | Z | Z_A Z_B | X |
| XX (patches side by side) | \|0⟩ | X | X_A X_B | Z |

**Routing strips.** Two patches that are not adjacent are joined through a routing strip: a block of data qubits, d
wide and g long (g odd), between them. In the builders this is the `gap` argument; g = 1 is the single gap row or
column above. Patch B then sits at (0, d+g) for a ZZ measurement and at (d+g, 0) for XX, and the merged region is
d × (2d+g) or (2d+g) × d. Nothing else changes in the protocol: the strip qubits are prepared like gap qubits, the
strip is merged with both patches for the rounds of the measurement, and it is read out at the split.

- The new faces of the measured type are the (d+1)(g+1)/2 faces of that type in the strip and at its two ends. Every
  strip qubit lies in exactly two of them and every seam-row qubit of a patch in one, so their product is still
  Z_A Z_B (or X_A X_B).
- The strip also has (d−1)(g−1)/2 faces of the other type that touch no patch. Their data are prepared in their basis,
  so they are deterministic from the first merged round, and they end at the split (see Detectors).
- Each strip row adds d data qubits and d measure qubits (d = 5: a strip of 3 has 20 more qubits than adjacent
  patches). It adds no rounds: the measurement still takes d merged rounds, whatever the strip length.

Strips are straight. A route with a bend, or one that passes other patches, is not built (see Limitations).

## Schedules

A `SurgerySchedule` is a list of phases. Each phase cuts the lattice into regions (a patch on its own, or two patches
merged) and runs a number of syndrome rounds; the round circuit is the one of `memory_experiment` (H on the measure
qubits, four CZ steps with H on the data of X faces around each step, H, measurement, reset). Between phases, data
qubits that leave a region (gap qubits at a split) are read out in their basis and data qubits that join (gap qubits at
a merge) are reset and prepared. The first phase prepares each patch transversally in |0⟩ or |+⟩ (with an optional
logical X or Z for |1⟩ and |−⟩); after the last phase every data qubit is read out in its patch's basis. Ordering
details follow `memory_experiment` exactly, so a one-patch schedule produces the same circuit gate for gate.

The protocols:

- `surgery_parity(d, 'Z' | 'X', a, b, rounds, gap)`: one separate round, d merged rounds, one round after the split.
- `surgery_cnot(d, c, t, readout, rounds, require_observable, gap)`: control C at (0, 0), ancilla A at (0, d+g) in |+⟩,
  target T at (d+g, d+g), with g = `gap` (1 by default). One separate round, d rounds of the ZZ merge of C and A, d
  rounds of the XX merge of A and T, one separate round, then C and T are read out in the chosen basis and A in Z. The
  ZZ split and the XX merge happen at the same phase boundary. With g > 1 there is a strip between C and A and one
  between A and T, both of length g.

## Detectors

Detectors are derived from the schedule, not written out per protocol. When a face is measured, it is compared with its
previous measurement. Since then, data qubits may have left its support or joined it:

- a qubit that left must have been read out in the face's basis; its outcome joins the detector;
- a qubit that joined must have been prepared in the face's basis; it contributes nothing;
- otherwise the face's value is random and there is no detector.

A face measured for the first time gives a detector on its own if all its data were prepared in its basis (the first
round of a transversally prepared patch, and the faces of the other type inside a routing strip), and none otherwise
(the seam faces of a merge: their first outcome is random).
A face that ends at a split closes there if its whole support is read out in its basis at that split: its last
measurement and those outcomes form a detector. These are the faces of the other type inside a routing strip; without
them a readout error on a strip qubit could flip a logical that crosses the strip without firing any detector (tested).
A face that ends is forgotten, so a face that returns in a later phase starts afresh.
At the final readout each face whose type matches its patch's readout basis closes with its last measurement and its
data outcomes. Each detector has a type; a run decodes the detectors of one type, as the memory experiments do.

The noiseless value of every classical bit (`reference`) follows from the same rules: 0 for a face that is
deterministic from preparation, carried along by the comparisons, and random (−1) otherwise. It sets the asymmetric
readout error of each measurement, as in the memory experiments.

## Observables and the Pauli frame

Each observable is a set of classical bits whose XOR is deterministic in the noiseless circuit, with an expected value.

**Parity measurement.** With both inputs in the measured basis, the observables are the merge outcome m (the product of
the new seam faces in the first merged round) and each patch's logical next to the seam. With both inputs in the other
basis the merge must preserve the product of the other logicals: X_A X_B through a ZZ merge is column 0 of A, column 0
of B and the gap qubits between them (one per strip row), read out at the split.

**CNOT.** Write a for the ZZ outcome, b for the XX outcome, c_A for the ancilla's final Z, g for the A–T gap qubits in the
logical Z row and h for the C–A gap qubits in the logical X column (one qubit each for adjacent patches, one per strip
row or column otherwise), all read out at the splits, and f for the product of A's X faces in the first XX round. The logical X of A is used on its left column at the C–A seam and on its right column
at the A–T seam; the two differ by the product of all of A's X faces, which f supplies. Then

- Z_T(out) = Z_T(final, row 0) + Z_A(final, row 0) + g + a, so the correction is X_T^(a + c_A + g);
- X_C(out) = X_C(final, column 0) + h + f + b, so the correction is Z_C^(b + h + f);
- Z_C(out) = Z_C(final) and X_T(out) = X_T(final).

The protocol is Clifford and the corrections are Paulis, so they are applied to the outcomes in software rather than
fed forward. The rows and columns are chosen next to each seam so that they match the seam products: C's bottom row, A's
top row (both at the C–A seam and in the A–T merged row), T's row 0; C's and A's column 0 (C–A seam), A's last and T's
first column (A–T seam). Each observable is the output logical with the correction applied: for Z-basis inputs Z_C = z_C
and Z_T = z_C + z_T; for X-basis inputs X_C = x_C + x_T and X_T = x_T; for a Bell-pair input (control in X, target in
Z) Z_C Z_T = z_T or X_C X_T = x_C.

## Decoding

`prepare_surgery` builds the frame program with `build_frame_program` (the memory experiments' noise model: Pauli-
twirled T1/T_phi with state-aware decoherence, coherent-ZZ windows, control errors, readout and reset errors, leakage)
and enumerates every single fault once with `enumerate_faults` (`surgery_faults`). To get everything from one
enumeration, each observable is appended to the detector list as an extra detector and split off again, so a fault
carries the decoded detectors it flips and a mask of the observables it flips. A run decodes the detectors of one Pauli
type, the type of its observables. There are five decoders (`SurgeryDecoder`, `qpuforge surgery --decoder`):

- **Separate matching** (`Separate`, the default). One detector error model (`build_dem`) and one exact matching
  decoder (`MatchingDecoder`) per observable. The graphs have the same detectors; parallel edges merge differently
  depending on which observable they flip.
- **Joint matching** (`Joint`). One graph for all observables (`build_joint_dem`, `SurgerySetup::joint`). Each pair of
  detectors has one edge. Its probability is that of all mechanisms flipping exactly that pair, and it carries the
  observable mask of the most likely class among them. One exact matching returns the XOR of the masks along the
  correction, so all observables are read off the same correction and are consistent with each other. An edge whose
  mechanisms do not all flip the same observables is counted (`mixed_edges`, with `minority_p`, the probability of the
  mechanisms that disagree with the edge's mask): this is the error mass a single mask per edge cannot represent.
- **Union-find** (`SeparateUF`, `JointUF`). The weighted union-find decoder of `unionfind.hpp` on the same graphs.
  `JointUF` peels one correction and XORs the masks of its edges.
- **Correlated matching** (`Correlated`, needs `prepare_surgery(..., correlated = true)`). Two passes that use the
  detectors of both Pauli types. A mechanism that flips one or two detectors of each type (a Y error on a data qubit,
  most two-qubit errors of a CZ) is an edge e of the decoded graph and an edge o of the other type's graph at once.
  Pass 1 decodes the other type's detectors with union-find. For every edge o of that correction, the weight of each
  edge e that shares mechanisms with o is lowered to log((1 − p)/p) with p = p(e and o) / p(o), the probability of e
  given that o fired. Pass 2 is exact matching on the joint graph with those weights (`decode_scaled`, the routine
  that re-weights edges for heralds). Weights are only lowered, never raised, and the edge masks are not changed.

**Mechanisms that flip more than two decoded detectors** (hyperedges) are handled as before in every decoder: they are
left out of the graph and their total probability is reported (`hyperedge_p`). The correlated decoder does use
mechanisms that flip up to two detectors of each type, which the single-type graphs see only as two unrelated edges.

## Monte Carlo

`run_surgery_decoders` samples the frame program with the same per-shot function as `run_memory_frame`
(`frame_detail::sample_shot` in `qec.hpp`; the surgery runner uses neither heralds nor impacts) and keeps its determinism contract (docs/design/deterministic-mc.md): each shot reseeds from
`frame_shot_seed(seed, shot)`, shots are grouped in blocks of `kFrameBlock`, and block tallies are added in block order.
The result depends only on the setup, the shot count and the seed. For a single-patch schedule it reproduces
`run_memory_frame` shot for shot (tested), so a change to the shared sampler changes both runners together.

`run_surgery_decoders` decodes every shot with each decoder in a list, so decoders are compared on the same syndromes.
For each it returns the error of every observable, the shots with any observable wrong, and against the first decoder
of the list the shots whose prediction differs, the shots it fixes and the shots it breaks. It also returns the mean
decode time per shot, which is wall time and the only field that is not reproducible. `run_surgery_frame(M, shots,
seed, decoder)` is the same run with one decoder; its default is separate matching.

## Chips

`surgery_chip(lattice, preset)` lays out the lattice as `generate()` lays out a surface-code chip (frequency plan by
site colour, 45° placement, couplers in CZ order) and applies the physics of `alpha1_v11` or `beta1`. Both use Beta-1's
readout plan (at most 9 resonators per feedline): Alpha-1's three feedlines do not scale to 169 qubits. For one d = 3
patch with `beta1` the chip equals `beta1(3)` qubit for qubit (tested).

Strip qubits get their frequencies from the same colouring by site, so a lattice with strips uses the same eight
design frequencies and the same sixteen (measure, data) frequency pairs as one patch, and the CZ calibrations of
`results/cz_cache.json` cover it. The study's chip account (section 6, see Results) counts qubits, couplers, feedlines
and control lines for two and three d = 5 patches with and without strips and runs the collision rules and the
fabrication Monte Carlo on them. It is an account of the lattice only: there is no placement of control lines, no
cell layout and no mask for these chips.

## Verification

- `check_surgery_noiseless` runs the logical circuit on a stabiliser tableau (Aaronson–Gottesman) with random
  measurement outcomes: no detector of either type fires, every observable takes its expected value and every bit marked
  deterministic is deterministic with that value. The tests cover ZZ and XX on |00⟩ |01⟩ |10⟩ |11⟩ |++⟩ |+−⟩ |−+⟩ |−−⟩
  at d = 3 and 5, and the CNOT on the Z- and X-basis inputs at d = 3 and 5 and on Bell-pair inputs at d = 3, each
  with its own seed; at d = 3 also with the merge starting in the first round (`pre = 0`). Two negative controls (an
  X_A X_B observable without the split's gap qubit, and an inserted X error) must be caught.
- Process check with reference qubits (`surgery_cnot_reference`): two extra qubits start in Bell pairs with C and T
  (a controlled logical X, or a controlled logical Z followed by H on the reference), so every patch logical is random.
  After the CNOT, Z_RC Z_C, Z_RT Z_C Z_T, X_RC X_C X_T and X_RT X_T must all be +1, with the corrected output logicals,
  for Z- and X-basis patches at d = 3 and 5. This sees frame-assembly errors that fixed logical inputs hide: adding
  the control's round-0 faces of the readout type to its corrected logical is caught.
- A single patch through the surgery engine gives the same logical errors and detection rate as
  `memory_experiment`/`run_memory_frame` in both bases, after every prefix of 1 to 120 shots (so shot by shot) and over
  6000 shots, with plain noise and with leakage on (CZ leakage, a plain or a leakage-removing reset, DQLR).
- The d = 3 detector error models of ZZ, XX and the CNOT have no undetectable logical faults, and `prepare_surgery`
  refuses an observable that a single fault can flip without a decoded detector (a mis-typed observable). Distance
  check with uniform edge weights, built from every mechanism (parallel edges keep the observable flag with the larger
  probability) and decoding every mechanism, hyperedges included: on the Beta-1 chip every fault is corrected except a
  few coherent-ZZ twirl terms, 2 for the ZZ measurement's Z-type decoders (total p 6.7e-8), 3 for the CNOT's Z-type
  decoders (2.1e-7, one of them a hyperedge) and 1 for its X-type decoders (3.3e-8). They come from as-soon-as-possible
  scheduling: a seam measure qubit can start the next round while a neighbouring gap qubit is still in the CZ steps of
  the previous one, so a twirled conditional phase between them straddles two rounds and acts as a weight-2 error. The
  test allows twirl terms only, and bounds their total probability.
- The noiseless frame run of a ZZ measurement has no defects and no errors.
- Joint decoding (`surgery_joint` test group). A matching decoder built from a hand-made edge list returns the XOR of
  the masks along the matching, with matching and with union-find. For the d = 3 ZZ measurement and CNOT in both
  bases, the joint model has the same detector pairs as each per-observable model, and each edge's probability is
  that of the pair's two classes combined. With uniform weights the joint decoder corrects every single fault except
  the coherent-ZZ twirl terms listed above. On 10,000 shots, `run_surgery_decoders` reproduces `run_surgery_frame` for
  the separate and the joint decoder, the tallies of every decoder are consistent, union-find stays within 1.3× of
  matching, and the correlated decoder fixes more shots than it breaks. Without noise every decoder is silent.
- Routing strips (`surgery_routing` test group). Qubit and face counts of the strip at d = 3 and 5 for strips of 1, 3
  and 5, and one closing detector per strip face of the other type. A strip of length 1 gives the adjacent-patch
  circuit. The tableau check passes for ZZ and XX on the eight inputs (strips 3 and 5 at d = 3, strip 3 at d = 5; at
  d = 3 with a strip of 3 also merging from the first round), for the CNOT on Z- and X-basis inputs, on Bell-pair
  inputs (d = 3) and with reference qubits. Negative controls: dropping any strip qubit from X_A X_B, and an X error
  in the middle of the strip, are caught; with the strip's closing detectors removed, `prepare_surgery` refuses
  X_A X_B. The d = 3 detector error models with strips of 3 have no undetectable logical faults, and uniform-weight
  decoding leaves only coherent-ZZ twirl terms uncorrected (at most 4 mechanisms, total p 2.4e-7, for the CNOT's
  Z-type detectors; the same for separate and joint decoding). A strip of 5 raises the error of the ZZ outcome at
  d = 3 by more than 1.5×, at the same number of rounds.
- The thread-count determinism digest covers every decoder on the same shots and a routed XX measurement.

## Results

From `results/surgery_study.txt` (CZ model v4; 100,000 shots per point; one separate round, d rounds per merge, one round after the
split; ± is the 1σ Monte Carlo error). Logical error per operation:

| Operation | Observable | Alpha-1 v1.1, d = 3 | d = 5 | Beta-1, d = 3 | d = 5 |
|---|---|---|---|---|---|
| ZZ on \|00⟩ | m_ZZ | 2.35 % ± 0.05 | 0.82 % ± 0.03 | 0.74 % ± 0.03 | 0.14 % ± 0.01 |
| | any observable wrong | 4.11 % | 1.48 % | 1.26 % | 0.24 % |
| XX on \|++⟩ | m_XX | 2.29 % ± 0.05 | 0.80 % ± 0.03 | 0.76 % ± 0.03 | 0.16 % ± 0.01 |
| | any observable wrong | 4.05 % | 1.46 % | 1.25 % | 0.27 % |
| CNOT on \|00⟩ (Z readout) | Z_C / Z_T | 2.37 % / 7.81 % | 1.22 % / 3.71 % | 0.69 % / 2.24 % | 0.16 % / 0.58 % |
| | any observable wrong | 9.49 % | 4.69 % | 2.74 % | 0.71 % |
| CNOT on \|++⟩ (X readout) | X_C / X_T | 8.44 % / 2.47 % | 4.20 % / 1.25 % | 2.45 % / 0.67 % | 0.65 % / 0.15 % |
| | any observable wrong | 9.22 % | 4.47 % | 2.71 % | 0.70 % |

- Going from d = 3 to d = 5 lowers every error: by 1.9–2.9× with Alpha-1 v1.1 physics and by 3.8–6.4× with Beta-1
  physics.
- In the CNOT, the observable that carries the corrections (Z_T in the Z basis, X_C in the X basis) fails 3.0–4.3 times
  as often as the other: it collects the merge outcome, the ancilla's readout and a split byproduct, each
  of which is decoded through the merged region and its time boundaries.
- The CNOT takes 2d + 2 rounds (5.98 µs at d = 3, 9.06 µs at d = 5 on Beta-1) and 169 qubits at d = 5.
- The 52 noiseless tableau checks in the study (d = 3 and 5) all pass.

### Decoders on the same shots

Section 4 of the study decodes the shots of the table above again with every decoder. Any observable wrong, and the
mean decode time per shot on one thread (wall time, indicative only):

| Operation (d = 5) | Separate matching | Joint matching | Joint union-find | Correlated matching |
|---|---|---|---|---|
| Alpha-1 v1.1, ZZ on \|00⟩ | 1.48 % (35 µs) | 1.48 % (10 µs) | 1.59 % (9 µs) | 0.98 % (123 µs) |
| Alpha-1 v1.1, CNOT on \|00⟩ | 4.69 % (102 µs) | 4.69 % (42 µs) | 5.07 % (22 µs) | 2.89 % (604 µs) |
| Beta-1, ZZ on \|00⟩ | 0.24 % (12 µs) | 0.24 % (4 µs) | 0.25 % (4 µs) | 0.16 % (35 µs) |
| Beta-1, XX on \|++⟩ | 0.27 % (13 µs) | 0.27 % (4 µs) | 0.29 % (4 µs) | 0.17 % (35 µs) |
| Beta-1, CNOT on \|00⟩ | 0.71 % (38 µs) | 0.71 % (15 µs) | 0.75 % (13 µs) | 0.42 % (184 µs) |
| Beta-1, CNOT on \|++⟩ | 0.70 % (47 µs) | 0.70 % (18 µs) | 0.74 % (17 µs) | 0.39 % (224 µs) |

- **Joint matching gives the same prediction as the separate decoders in every one of the 1,600,000 shots** (16
  points, both chips, d = 3 and 5). The reason is in the model: at d = 5 no edge has mechanisms with different
  observable masks, and at d = 3 one or two edges do, with a total probability below 1e-7. The disagreement on ties
  that the PR #21 note expected does not occur at these sizes. What joint decoding buys here is consistency by
  construction and time: one matching instead of two or three, 1.9–3.4× less decode time.
- Joint and separate union-find also agree shot for shot. Union-find gets some observable wrong 1.02–1.09 times as
  often as matching.
- **Correlated matching lowers "any observable wrong" to 0.56–0.85 of joint matching**: 0.81–0.85 at d = 3 and
  0.56–0.69 at d = 5. With Beta-1 physics at d = 5 a CNOT on |00⟩ goes from 0.71 % to 0.42 % (389 shots fixed, 103
  broken, of 100,000). The gain grows with distance, as expected when the decoder uses more of the information.
- **Its cost: 7–16 times the decode time of joint matching** (d = 5 CNOT with Beta-1 physics: 184 µs per shot
  against 15 µs; with Alpha-1 v1.1 physics 604 µs against 42 µs). The second pass rebuilds shortest paths over the
  defects and the ends of every re-weighted edge for each shot. As it stands it is an offline reference, not a
  real-time candidate: a 12-round CNOT lasts 9 µs.

### Routing strips

Section 5 of the study, Beta-1 physics, joint matching, 100,000 shots per point (strip 1 is the adjacent-patch row of
the first table). Logical error per operation, and the footprint:

| d | Operation | Observable | Strip 1 | Strip 3 | Strip 5 | Strip 7 |
|---|---|---|---|---|---|---|
| 3 | ZZ on \|00⟩ | m_ZZ | 0.74 % | 1.31 % | 1.87 % | 2.17 % |
| | | any observable wrong | 1.26 % | 1.83 % | 2.42 % | 2.71 % |
| | | qubits | 41 | 53 | 65 | 77 |
| 3 | CNOT on \|00⟩ | Z_C / Z_T | 0.69 % / 2.24 % | 0.66 % / 3.50 % | 0.69 % / 4.87 % | 0.62 % / 6.02 % |
| | | any observable wrong | 2.74 % | 3.99 % | 5.39 % | 6.48 % |
| | | qubits | 65 | 89 | 113 | 137 |
| 5 | ZZ on \|00⟩ | m_ZZ | 0.14 % ± 0.01 | 0.24 % ± 0.02 | 0.35 % ± 0.02 | 0.45 % ± 0.02 |
| | | any observable wrong | 0.24 % | 0.37 % | 0.45 % | 0.55 % |
| | | qubits | 109 | 129 | 149 | 169 |
| 5 | XX on \|++⟩ | m_XX | 0.16 % ± 0.01 | 0.25 % ± 0.02 | 0.39 % ± 0.02 | 0.49 % ± 0.02 |
| | | any observable wrong | 0.27 % | 0.35 % | 0.48 % | 0.61 % |
| 5 | CNOT on \|00⟩ | Z_C / Z_T | 0.16 % / 0.58 % | 0.15 % / 0.92 % | 0.15 % / 1.32 % | 0.16 % / 1.69 % |
| | | any observable wrong | 0.71 % | 1.03 % | 1.44 % | 1.82 % |
| | | any wrong, correlated matching | 0.42 % | 0.66 % | 0.89 % | 1.11 % |
| | | qubits | 169 | 209 | 249 | 289 |
| 5 | CNOT on \|++⟩ | X_C / X_T | 0.65 % / 0.15 % | 0.88 % / 0.15 % | 1.20 % / 0.16 % | 1.48 % / 0.16 % |
| | | any observable wrong | 0.70 % | 0.93 % | 1.25 % | 1.53 % |

- **The strip costs error in the merge outcome only, and about linearly in its length.** At d = 5 the ZZ outcome
  goes from 0.14 % to 0.45 % between strips of 1 and 7, about 0.05 % per strip row. The patch logicals next to the
  seam (Z_A, Z_B, and Z_C or X_T in the CNOT) do not change within the Monte Carlo error. The merge outcome is the
  product of all seam faces of the first merged round, so each strip row adds faces whose measurement error can flip
  it, while the d merged rounds that protect it stay the same.
- Any observable wrong rises by 1.3–1.6× for a strip of 3 and 2.2–2.6× for a strip of 7, at d = 3 and at d = 5. Going
  from d = 3 to d = 5 lowers the error by about the same factor at every strip length (CNOT on |00⟩: 3.9× adjacent,
  3.6× with strips of 7).
- **Time: none.** The operations take the same number of rounds whatever the strip (5 and 7 rounds for a parity
  measurement at d = 3 and 5, 8 and 12 for the CNOT), and the same duration within 5 % (d = 5 CNOT: 9.06 µs adjacent,
  9.17 µs with strips of 7; the small differences come from the as-soon-as-possible schedule and the readout plan).
- **Footprint:** 2d qubits per strip row. A d = 5 ZZ or XX measurement needs 109 qubits with adjacent patches and
  169 with a strip of 7; the three-patch CNOT needs 169 and 289.
- The 108 noiseless tableau checks with strips of 3, 5 and 7 (d = 3 and 5) all pass.
- The correlated decoder keeps its gain with strips (d = 5 CNOT on |00⟩ with strips of 7: 1.82 % to 1.11 %).

### Chip account for two and three d = 5 patches

Section 6 of the study: `surgery_chip` with Beta-1 physics for each lattice. Lines are counted as in the design
reviews (one XY and one flux line per qubit, one flux line per coupler, two per feedline).

| Lattice (d = 5) | Qubits (data + measure) | Couplers | Feedlines | Lines |
|---|---|---|---|---|
| Beta-1, one patch (reference) | 49 (25 + 24) | 80 | 6 | 190 |
| Gamma-1, d = 7 memory (reference) | 97 (49 + 48) | 168 | 11 | 384 |
| 2 patches, adjacent (ZZ or XX) | 109 (55 + 54) | 188 | 13 | 432 |
| 2 patches, strip 3 | 129 (65 + 64) | 224 | 15 | 512 |
| 2 patches, strip 5 | 149 (75 + 74) | 260 | 17 | 592 |
| 2 patches, strip 7 | 169 (85 + 84) | 296 | 19 | 672 |
| 3 patches, CNOT, adjacent | 169 (85 + 84) | 296 | 19 | 672 |
| 3 patches, CNOT, strips 3 | 209 (105 + 104) | 368 | 24 | 834 |
| 3 patches, CNOT, strips 5 | 249 (125 + 124) | 440 | 28 | 994 |
| 3 patches, CNOT, strips 7 | 289 (145 + 144) | 512 | 33 | 1156 |

- **The 8-frequency plan fits.** Every lattice uses the same 8 design frequencies and the same 16 (measure, data)
  frequency pairs as Beta-1 and Gamma-1, so no new CZ calibration is needed; the collision rules flag no pair, and the
  largest residual ZZ is 5.0 kHz, as on Beta-1.
- Fabrication Monte Carlo (100 samples, static collision rules): with 0.5 % junction spread (after trimming) 99 of 100
  chips are collision-free for two adjacent patches, 98 with a strip of 3 and 99 for the three adjacent patches of the
  CNOT, against the 80 % of the Gamma gate. At 2 % spread (as fabricated) the numbers are 91, 86 and 82.
- **The smallest "2 × d5 + routing" chip is larger than Gamma-1.** Two adjacent patches (one ZZ or one XX
  measurement, no CNOT) need 109 qubits and 432 lines, against 97 and 384 for Gamma-1. A CNOT needs the third patch:
  169 qubits and 672 lines at least. Wiring, already the open item for Gamma-1, is the constraint here too.
- This is a count on the lattice. No design file, cell layout or mask is made for these chips, and one lattice serves
  one operation: a ZZ and an XX measurement between the same two patches need different lattices (stacked or side by
  side), or the L of the CNOT.

## Limitations

- The patches are placed on a fixed lattice per protocol (stacked for ZZ, side by side for XX, an L for the CNOT).
  Routing strips are straight and as wide as the patches, and the two strips of the CNOT have the same length; there
  is no route with a bend or around other patches, and no general router.
- Corrections are applied in post-processing; real-time feed-forward is not modelled.
- The merge still runs for d rounds whatever the strip length, so the error of the merge outcome grows with the strip
  (see Results). More merged rounds for long strips have not been studied.
- The joint decoder gives each edge one observable mask. Where mechanisms of one detector pair flip different
  observables it keeps the most likely class; at d = 3 and 5 this concerns at most 3 edges per graph and a
  probability below 1e-7.
- The correlated decoder's first pass is union-find, not matching; weights are lowered and never raised; the edge
  masks are not updated; and it is 7–16 times slower than joint matching (see Results). It has not been tried with a
  union-find second pass.
- Mechanisms that flip more than two decoded detectors stay out of every graph. The CNOT's coherent-ZZ twirls leave a
  few of them in the Z-basis model, of total probability 1.2e-6 at d = 3 and 2.4e-6 at d = 5 with adjacent patches,
  and 4.3e-6 at d = 5 with strips (`hyperedge_p` and `routing[].hyperedge_p` in the study's JSON).
- Rounds are scheduled as soon as possible, not in lockstep, so a few coherent-ZZ twirl terms straddle two rounds
  (see Verification; total probability at most 2.1e-7 on the Beta-1 chip at d = 3 with adjacent patches, 2.4e-7 with
  strips of 3).
- The chip account counts qubits, couplers and lines and checks frequencies. It has no layout, no control-line routing
  and no mask.
