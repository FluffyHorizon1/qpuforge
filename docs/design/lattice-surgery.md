# Design note: lattice surgery

Status: implemented in `include/qpuforge/surgery.hpp` (Phase Gamma, ROADMAP.md). Study: `examples/surgery_study.cpp`,
results in `results/surgery_study.{txt,json}`.

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

## Schedules

A `SurgerySchedule` is a list of phases. Each phase cuts the lattice into regions (a patch on its own, or two patches
merged) and runs a number of syndrome rounds; the round circuit is the one of `memory_experiment` (H on the measure
qubits, four CZ steps with H on the data of X faces around each step, H, measurement, reset). Between phases, data
qubits that leave a region (gap qubits at a split) are read out in their basis and data qubits that join (gap qubits at
a merge) are reset and prepared. The first phase prepares each patch transversally in |0⟩ or |+⟩ (with an optional
logical X or Z for |1⟩ and |−⟩); after the last phase every data qubit is read out in its patch's basis. Ordering
details follow `memory_experiment` exactly, so a one-patch schedule produces the same circuit gate for gate.

The protocols:

- `surgery_parity(d, 'Z' | 'X', a, b)`: one separate round, d merged rounds, one round after the split.
- `surgery_cnot(d, c, t)`: control C at (0, 0), ancilla A at (0, d+1) in |+⟩, target T at (d+1, d+1). One separate
  round, d rounds of the ZZ merge of C and A, d rounds of the XX merge of A and T, one separate round, then C and T are
  read out in the chosen basis and A in Z. The ZZ split and the XX merge happen at the same phase boundary.

## Detectors

Detectors are derived from the schedule, not written out per protocol. When a face is measured, it is compared with its
previous measurement. Since then, data qubits may have left its support or joined it:

- a qubit that left must have been read out in the face's basis; its outcome joins the detector;
- a qubit that joined must have been prepared in the face's basis; it contributes nothing;
- otherwise the face's value is random and there is no detector.

A face measured for the first time gives a detector on its own if all its data were prepared in its basis (the first
round of a transversally prepared patch), and none otherwise (the seam faces of a merge: their first outcome is random).
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
of B and the gap qubit between them, read out at the split.

**CNOT.** Write a for the ZZ outcome, b for the XX outcome, c_A for the ancilla's final Z, g for the A–T gap qubit in the
logical Z row and h for the C–A gap qubit in the logical X column, both read out at the splits, and f for the product of
A's X faces in the first XX round. The logical X of A is used on its left column at the C–A seam and on its right column
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
twirled T1/T_phi with state-aware decoherence, coherent-ZZ windows, control errors, readout and reset errors, leakage),
enumerates every single fault once with `enumerate_faults`, and builds one detector error model per decoded observable
with `build_dem`. To get all observables from one enumeration, each observable is appended to the detector list as an
extra detector and split off again. Each observable gets its own exact matching decoder (`MatchingDecoder`): the
matching graph is the same, but parallel edges merge differently depending on which observable they flip.

The decoders are independent, so "any observable wrong" combines matchings that can differ on ties. For the parity
measurement and the CNOT these differences are rare; a joint decoder would have to carry an observable mask per edge.

## Monte Carlo

`run_surgery_frame` samples the frame program with the same per-shot function as `run_memory_frame`
(`frame_detail::sample_shot` in `qec.hpp`; the surgery runner uses neither heralds nor impacts) and keeps its determinism contract (docs/design/deterministic-mc.md): each shot reseeds from
`frame_shot_seed(seed, shot)`, shots are grouped in blocks of `kFrameBlock`, and block tallies are added in block order.
The result depends only on the setup, the shot count and the seed. For a single-patch schedule it reproduces
`run_memory_frame` shot for shot (tested), so a change to the shared sampler changes both runners together.

## Chips

`surgery_chip(lattice, preset)` lays out the lattice as `generate()` lays out a surface-code chip (frequency plan by
site colour, 45° placement, couplers in CZ order) and applies the physics of `alpha1_v11` or `beta1`. Both use Beta-1's
readout plan (at most 9 resonators per feedline): Alpha-1's three feedlines do not scale to 169 qubits. For one d = 3
patch with `beta1` the chip equals `beta1(3)` qubit for qubit (tested).

## Verification

- `check_surgery_noiseless` runs the logical circuit on a stabiliser tableau (Aaronson–Gottesman) with random
  measurement outcomes: no detector of either type fires, every observable takes its expected value and every bit marked
  deterministic is deterministic with that value. The tests cover ZZ and XX on |00⟩ |01⟩ |10⟩ |11⟩ |++⟩ |+−⟩ |−+⟩ |−−⟩
  at d = 3 and 5, and the CNOT on the Z- and X-basis inputs at d = 3 and 5 and on Bell-pair inputs at d = 3. Two
  negative controls (an X_A X_B observable without the split's gap qubit, and an inserted X error) must be caught.
- A single patch through the surgery engine gives the same logical errors and detection rate as
  `memory_experiment`/`run_memory_frame` in both bases.
- The d = 3 detector error models of ZZ, XX and the CNOT have no undetectable logical faults. With uniform edge weights
  every single fault is corrected by every observable's decoder (circuit distance at least 3). Mechanisms below 1e-7
  are left out of that check: the coherent-ZZ twirl places a long idle window's phase at the window's end, and for
  measure qubits that idle until a merge a few such mechanisms (about 1e-8) end up parallel to a time-boundary edge,
  where uniform weights would tie with it.
- The noiseless frame run of a ZZ measurement has no defects and no errors.

## Results

From `results/surgery_study.txt` (100,000 shots per point; one separate round, d rounds per merge, one round after the
split; ± is the 1σ Monte Carlo error). Logical error per operation:

| Operation | Observable | Alpha-1 v1.1, d = 3 | d = 5 | Beta-1, d = 3 | d = 5 |
|---|---|---|---|---|---|
| ZZ on \|00⟩ | m_ZZ | 2.39 % ± 0.05 | 0.88 % ± 0.03 | 0.82 % ± 0.03 | 0.19 % ± 0.01 |
| | any observable wrong | 4.21 % | 1.52 % | 1.36 % | 0.30 % |
| XX on \|++⟩ | m_XX | 2.34 % ± 0.05 | 0.84 % ± 0.03 | 0.75 % ± 0.03 | 0.15 % ± 0.01 |
| | any observable wrong | 4.10 % | 1.47 % | 1.31 % | 0.25 % |
| CNOT on \|00⟩ (Z readout) | Z_C / Z_T | 2.45 % / 7.93 % | 1.23 % / 3.85 % | 0.69 % / 2.29 % | 0.15 % / 0.60 % |
| | any observable wrong | 9.68 % | 4.83 % | 2.79 % | 0.72 % |
| CNOT on \|++⟩ (X readout) | X_C / X_T | 8.48 % / 2.57 % | 4.22 % / 1.17 % | 2.42 % / 0.68 % | 0.65 % / 0.15 % |
| | any observable wrong | 9.30 % | 4.52 % | 2.67 % | 0.71 % |

- Going from d = 3 to d = 5 lowers every error: by 2.0–2.8× with Alpha-1 v1.1 physics and by 3.7–5.5× with Beta-1
  physics.
- In the CNOT, the observable that carries the corrections (Z_T in the Z basis, X_C in the X basis) fails 3.1–4.4 times
  as often as the other: it collects the merge outcome, the ancilla's readout and a split byproduct, each
  of which is decoded through the merged region and its time boundaries.
- The CNOT takes 2d + 2 rounds (5.98 µs at d = 3, 9.06 µs at d = 5 on Beta-1) and 169 qubits at d = 5.
- The 52 noiseless tableau checks in the study (d = 3 and 5) all pass.

## Limitations

- The patches are placed on a fixed lattice per protocol (stacked for ZZ, side by side for XX, an L for the CNOT); there
  is no general router.
- Corrections are applied in post-processing; real-time feed-forward is not modelled.
- Independent decoders per observable (see Decoding).
- The CNOT's coherent-ZZ twirls leave a few hyperedge mechanisms in the Z-basis model, of total probability 1.2e-6
  at d = 3 and 2.4e-6 at d = 5 (reported in the study), which the matching graph does not contain.
