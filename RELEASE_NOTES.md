# qpuforge 0.3.0 release notes

First tagged release. Adds pulse-level CZ calibration, a leakage-modelled QEC frame simulator, layout synthesis with
GDSII output, and the Alpha-1 v1.1 (AF-17) and Beta-1 (BF-49) designs.

See `CHANGELOG.md` for the API changes, and the README "What's new in v0.3" table for a module-level summary.

## Highlights
- Pulse-level tunable-coupler CZ (`pulse.hpp`), with a versioned on-disk calibration cache.
- Leakage-modelled QEC (`qec.hpp`): leakage is simulated; the decoder does not use it.
- Layout synthesis and GDSII (`layout.hpp`, `gds.hpp`), including placement/spacing checks and merged-etch width/space DRC.
- Alpha-1 v1.1 (AF-17) and Beta-1 (BF-49) design reviews, regenerated from code with provenance.
- 101 validation checks, plus labelled CTest entries and two review smoke tests.

## Known issues
- Alpha-1 v1.1 pair Q8–Q11 leaks 2.3e-4 per CZ, above the 2e-4 gate limit.
- The AF-17 and BF-49 first-cut layouts fail the merged-etch width/space DRC (4 µm / 4 µm, projection metric):
  AF-17 has 34 width + 202 space violations, and BF-49 has 98 width + 602 space violations. The masks are not
  fabrication-ready.
