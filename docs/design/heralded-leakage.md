# Design note: leakage-aware decoding from heralded leakage

Status: implemented on `scs-bot/gamma-decoding` (`HeraldSpec`, `MatchingDecoder::decode_scaled`,
`herald_footprints` in `include/qpuforge/qec.hpp`).

## Model

The measure qubits are read out with three states. In `run_memory_frame`, at each measure-qubit measurement:
- a leaked qubit is reported as |2⟩ with probability `fidelity` (default 0.95);
- an unleaked one is misreported as |2⟩ with probability `false_rate` (default 1e-3).

The bit passed to the detectors does not change: a leaked qubit still reads as |1⟩, as in 0.4. One extra uniform draw
per measure-qubit measurement comes from the shot's own stream, so runs with heralds on stay deterministic.
`use_in_decoder = false` makes the same draws and ignores them. That gives the same syndromes decoded with and
without heralds.

Data qubits are not heralded. Their final readout is the logical measurement.

## Decoder: which edges a herald points at

A heralded bit tells the decoder that a measure qubit was leaked at its measurement. The errors that leakage causes
are:
- the measurement of that round is wrong;
- the CZ at which the qubit leaked randomises its partner;
- every later CZ in that round applies a random phase to its data partner. For an X-type measure qubit, the H gates
  around its CZs turn those phases into X errors on the data.

Two rules are implemented.

- **Footprint (default).** The edges of every single fault owned by the heralded measure qubit, from after its
  previous measurement up to this one, are multiplied by `footprint_factor` (0.4). These are faults on its own
  operations, or on either qubit of its CZs. The bit's own measurement-error edge is multiplied by `erase_factor`
  (0, an erasure). The map from bits to edges is built once with the decoder (`herald_footprints`, using the
  `owner` of each fault from `frame_detail::list_faults`).
- **Detector.** The edge joining the two detectors built from the bit gets `erase_factor`; every other edge touching
  them gets `touch_factor`. This is the obvious rule, but it misses the main mechanism. In a Z memory, a leaked
  X-type measure qubit puts X errors on its data qubits, yet its bits feed no detector there, so this rule never sees
  them. Erasing everything around a heralded detector also opens free paths to the boundary.

Measured on the same syndromes (`results/decoder_sweeps.txt`, section 1; d = 5 test chip, 15 rounds, measure-qubit CZ leakage 2e-3 per CZ, transport 0.05,
leakage-removing reset, 40,000 shots, exact matching): logical error per round 0.00531 without heralds. With the
detector rule, erasing every touching edge gave 0.00869 and erasing only the measurement edge gave 0.00543. With the
footprint rule at factor 0.4, it is 0.00400. Factors from 0.3 to 0.5 are within noise of each other; 0.15 and 0 are
worse (0.00437 and 0.00530).

## Matching on a re-weighted graph

The exact decoder works on precomputed all-pairs shortest paths, and heralds change weights per shot.
`decode_scaled(defects, scaled)` builds a small graph whose terminals are the defects, the boundary and the ends of
the modified edges. Its edges are the static shortest distances between all terminal pairs and the modified edges at
their new weights. Floyd–Warshall on this graph gives the exact modified shortest paths between defects, with their
logical parities: any shortest path splits at its modified edges into segments of unmodified edges, and each segment
is at least the static distance between its ends. The usual DP/blossom matching then runs on those distances. The
union-find decoder takes the modified lengths directly.

The matching code was refactored onto distance and parity oracles to allow this. The original path reads the static
matrix through the same oracle, so results without heralds are unchanged: the `matching` and `memory_bases` tests,
and the union-find study counts, are identical before and after.

## Validation

`heralds` (tests/test_qpuforge.cpp):
- 300 random defect sets on randomly re-weighted d = 3 graphs: the matching weight equals a brute-force
  Floyd–Warshall plus DP on the whole modified graph (worst difference 8.9e-16). The union-find correction clears
  every syndrome.
- Herald sampling: perfect heralds are 100% leaked; zero fidelity and false rate gives no heralds; the default gives
  86% leaked.
- d = 5, 10 rounds, CZ leakage 2e-3, leakage-removing reset, 20,000 shots, same syndromes: exact matching
  928 → 694 failures with the footprint rule (1,035 with the detector rule); union-find 984 → 763.

Gamma-1, memory-Z, 14 rounds, leakage-removing reset, CZ leakage as calibrated (×1) and ×10, same syndromes with
and without heralds (`results/decoder_study.txt`, section 3):

| d | leakage | heralds per shot (leaked) | exact: ε without → with | union-find: ε without → with |
|---|---|---|---|---|
| 5 | ×1 | 0.40 (16%) | 1.68e-4 → 1.57e-4 (−6%) | 1.85e-4 → 1.74e-4 (−6%) |
| 5 | ×10 | 0.95 (65%) | 6.30e-4 → 3.95e-4 (−37%) | 6.71e-4 → 4.20e-4 (−37%) |
| 7 | ×1 | 0.81 (17%) | 3.25e-5 → 2.93e-5 (−10%) | 3.65e-5 → 3.14e-5 (−14%) |
| 7 | ×10 | 2.00 (66%) | 2.16e-4 → 1.17e-4 (−46%) | 2.27e-4 → 1.26e-4 (−45%) |

At the chip's own leakage most heralds are false (the false-herald rate of 1e-3 per readout dominates), and the
gains at ×1 rest on 9 to 15 failures out of about 90 to 260: they are within two standard errors. At ×10 the gain is
large and clear.
