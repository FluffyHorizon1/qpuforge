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

## Limitations

- The patches are placed on a fixed lattice per protocol (stacked for ZZ, side by side for XX, an L for the CNOT); there
  is no general router.
- Corrections are applied in post-processing; real-time feed-forward is not modelled.
- Independent decoders per observable (see Decoding).
- Rounds are scheduled as soon as possible, not in lockstep, so a few coherent-ZZ twirl terms straddle two rounds
  (see Verification; total probability at most 2.1e-7 on the Beta-1 chip at d = 3).
- The CNOT's coherent-ZZ twirls leave a few hyperedge mechanisms in the Z-basis model, of total probability 1.2e-6
  at d = 3 and 2.4e-6 at d = 5 (reported in the study), which the matching graph does not contain.
