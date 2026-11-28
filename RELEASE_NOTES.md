# qpuforge release notes

## 0.4.0 (unreleased)

Finishes the software side of Phases Alpha and Beta and starts Phase Gamma.

See `CHANGELOG.md` for the API changes, and the README "What's new in v0.4" table for a module-level summary.

### Highlights
- Exact minimum-weight matching at any number of defects (`matching.hpp`). The greedy fallback above 22 defects is gone.
- Memory experiments in both bases, and correlated error bursts from high-energy impacts (`qec.hpp`).
- Hardware-gate evaluation from a measured-data file (`characterise.hpp`, `qpuforge hwgate`), rehearsed on a virtual chip.
- Beta-1 rev C: readout plan and cells corrected. Gamma-1 (GF-97): a 97-qubit distance-7 memory.
- Error bars of fitted logical error rates corrected; they were 20–40 % too small.
- Thread-count-independent Monte Carlo: every result is bit-identical for any `OMP_NUM_THREADS`.
- Hardened JSON (RFC 8259, nesting limit) and GDSII (I/O errors) serialisation.
- Alpha-1 v1.1, Beta-1 and Gamma-1 keep the 4 MHz direct coupling their design files were generated with; the CZ
  flux-noise term in `e_2q` is opt-in (Alpha-1 v1.0).
- 321 validation checks (186 in 0.3.0) and three review smoke tests.

### Known issues
- CZ leakage: the worst frequency combination leaks 2.3e-04 per gate, above the 2e-4 limit of the Alpha
  simulation gate (mean 5.7e-05). Alpha-1, Beta-1 and Gamma-1 share the same 16 frequency combinations, so the
  pair class appears on all three chips (Q8–Q11 on AF-17).
- No layout passes the merged-etch width/space DRC (4 µm / 4 µm, projection metric): AF-17 has 34 width + 202 space violations,
  BF-49 rev C 98 width + 504 space violations, and GF-97 194 width + 1,016 space violations. The
  masks are not fabrication-ready.
- Setting up a memory experiment takes time quadratic in the number of rounds (single-fault enumeration).
- The hardware-gate chain has only been run on synthetic data from the virtual chip.

## 0.3.0

First tagged release. Adds pulse-level CZ calibration, a leakage-modelled QEC frame simulator, layout synthesis with
GDSII output, and the Alpha-1 v1.1 (AF-17) and Beta-1 (BF-49) designs.

See `CHANGELOG.md` for the API changes, and the README "What's new in v0.3" table for a module-level summary.

### Highlights
- Pulse-level tunable-coupler CZ (`pulse.hpp`), with a versioned on-disk calibration cache.
- Leakage-modelled QEC (`qec.hpp`): leakage is simulated; the decoder does not use it.
- Layout synthesis and GDSII (`layout.hpp`, `gds.hpp`), including placement/spacing checks and merged-etch width/space DRC.
- Alpha-1 v1.1 (AF-17) and Beta-1 (BF-49) design reviews, regenerated from code with provenance.
- 186 validation checks, plus labelled CTest entries and two review smoke tests.

### Known issues
- Alpha-1 v1.1 pair Q8–Q11 leaks 2.3e-4 per CZ, above the 2e-4 gate limit.
- The AF-17 and BF-49 first-cut layouts fail the merged-etch width/space DRC (4 µm / 4 µm, projection metric):
  AF-17 has 34 width + 202 space violations, and BF-49 has 98 width + 602 space violations. The masks are not
  fabrication-ready.
