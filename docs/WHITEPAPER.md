# Transmon anharmonicity and the Alpha-1 processor

Finnian Bond
Secured Cyber Solutions Ltd
First written 25 November 2026 (PR #2). Revised for qpuforge v0.4 on 28 November 2026; the revision is a draft for the author's review.

## Abstract

Alpha-1 (AF-17) is a 17-qubit distance-3 surface-code design. It has not been fabricated. This note records the junction physics the design uses, why the |11⟩–|20⟩ frequency has to be solved at the tuned Josephson energy, the Alpha-1 coupler-transmon and CZ flux-noise models added on 25 November, and the AF17.V11 metal layout. Every error rate quoted here is from the qpuforge emulator in this repository, not from a device.

## Status and known issues

Alpha-1 has not been fabricated. Every number in this note is a qpuforge simulation result, not a measurement. Two known issues of the development tree concern the AF-17 design this note describes (`CHANGELOG.md`, [Unreleased]):

1. **The AF-17 mask is DRC- and connectivity-clean in the model (modelled airbridges treated as ideal ground); airbridges, EM and inductance are not verified, and no foundry DRC deck has been run. A review found the SQUID ground and flux-line return sharing one airbridge, flux lines that cannot tune the SQUIDs, and a coupler-dominated XY drive (`docs/design/mask-connectivity.md`).** The AF17.V11 layout passes the placement and island-spacing checks, the merged-etch width/space DRC at 4 µm / 4 µm (projection metric: **0 + 0 violations**, 34 + 202 in 0.4), and the metal connectivity check with 122 modelled ground straps (airbridges); its control lines are within a 10× T1 budget. BF-49 rev C and GF-97 pass the DRC and connectivity checks; on their cells the flux line reaches 9.7× the design T1, just short of that budget, which is therefore reported, not gated. The straps' capacitance to the coupler and bus strips they cross, and the flux lines' inductive coupling to the SQUIDs, are not modelled.
2. **CZ leakage is within the gate for every pair.** CZ model version 4 (`docs/design/cz-leakage.md`) re-calibrates the one pair class that failed the 2e-4 limit: Q8–Q11 goes from 2.33e-4 to 1.6e-5 per CZ. The largest leakage is now 1.2e-4 (Q8–Q10), the 24-pair mean 4.8e-5.

Sources: `designs/alpha1-v11/alpha1_v11_review.json` (`layout.drc`, `layout.connectivity`, `layout.control_lines`, `cz_calibration`, `exit_criteria`), `designs/beta1/beta1_review.json` and `designs/gamma1/gamma1_review.json` (`layout`), `docs/design/mask-connectivity.md`.

## Two designs, both called AF-17

The note uses "v1.0" and "v1.1" for two models of the same 17-qubit chip, and gives their numbers in separate tables. They are never mixed.

| | Alpha-1 v1.0 | Alpha-1 v1.1 |
|---|---|---|
| Couplers | abstract: g_max = 40 MHz. The 24 Nov review used 4 MHz direct coupling and a typed 5 kHz residual ZZ; since PR #2 `alpha1()` uses 1 MHz and a coupler transmon at the exchange null (idle ZZ −4 to −14 kHz, computed) | physical coupler transmon per pair (f_max 8.5 GHz, α_c −130 MHz, g_qc 160 MHz) |
| CZ | 23.3 ns, 2e-3 assumed control error | pulse-level, 3-knob calibration, 35.6–49.2 ns |
| Leakage | folded into the control error | simulated in the CZ (pulse level) and in the QEC frame model |
| Layout | floor plan (SVG), no mask | MoM-synthesised cells, GDSII, DRC, connectivity check with ground straps, control-line budget (no foundry DRC yet) |
| Files | `designs/alpha1/` (spec `ALPHA-1.md`; review output from 24 Nov, pre-PR #2, historical: see `designs/alpha1/README.md`) | `designs/alpha1-v11/` (spec `ALPHA-1-v1.1.md`, review regenerated for v0.4, `results/PROVENANCE.md`) |

The abstract v1.0 model is kept because the transmon physics below and the v1.0 error-correction ladder are built on it. The v1.1 numbers are the ones that describe the design that has a mask.

## The junction

A superconducting tunnel junction is two condensates separated by a barrier thin enough for Cooper pairs to tunnel. On this design the barrier is aluminium oxide, with gap Δ = 180 µeV. The phase difference φ obeys

I = Ic sin φ,    dφ/dt = 2eV / ℏ.

The energy stored in the junction is EJ (1 − cos φ), with EJ = ℏ Ic / (2e). At zero temperature the Ambegaokar–Baratoff relation fixes Ic from the normal-state resistance Rn of the same barrier. For Rn = 10 kΩ and Δ = 180 µeV the code returns EJ/h = 14.043 GHz, matching the hand calculation of 14.04 GHz (`results/test_output.txt`). Junction spread is therefore Rn spread, which is what the fabrication Monte Carlo perturbs.

The island also has a charging energy EC = e² / (2C). The Hamiltonian, in frequency units, is

H = 4 EC n² − EJ cos φ.

Alpha-1 parks every transmon at the flux sweet spot of its SQUID, where the first derivative of EJ with respect to flux is zero.

## Anharmonicity is not −EC

Expanding the cosine gives a harmonic oscillator at the plasma frequency ωp = √(8 EJ EC), plus a negative φ⁴ perturbation. First order in that perturbation produces the textbook result

f01 ≈ ωp − EC,    α ≡ f12 − f01 ≈ −EC.

The φ⁶ term and the second-order action of φ⁴ do not cancel. The exact charge-basis spectrum is softer. `results/test_output.txt` shows it against the Koch et al. asymptotics: at EJ/EC = 80, α = −277.4 MHz exact against −EC = −250 MHz. For the four Alpha-1 idle frequencies, solving the charge basis for the f01 in the design and α = −220 MHz gives:

| f01 (GHz) | EJ/EC | EC (MHz) | α exact (MHz) | −EC (MHz) | f23 − f12 (MHz) | f01 charge dispersion |
|------:|------:|---------:|--------------:|----------:|----------------:|----------------------:|
| 5.11 | 88.8 | 200 | −220 | −200 | −238 | 22 Hz |
| 5.29 | 94.3 | 200 | −220 | −200 | −237 | 11 Hz |
| 5.81 | 111.1 | 202 | −220 | −202 | −235 | 1.3 Hz |
| 5.99 | 117.2 | 202 | −220 | −202 | −235 | 0.6 Hz |

EJ/EC, EC, α and f23 − f12 are from the standalone check at the end of this note. The charge dispersions are the `disp_kHz` column for the `alpha1()` preset (`qpuforge new build/alpha1.json --preset alpha1`), whose qubits sit at these f01. The −220 MHz in the design is this exact α: `alpha1()` analyses at −220.0 MHz on every qubit. Rn is what moves f01. Charge dispersion of f01 is hertz, not megahertz: these ratios are in the transmon regime. The committed `designs/alpha1/alpha1.json` differs slightly: it has the 1 MHz links but the shunt capacitances generated for the 4 MHz coupling, so its qubits sit 5.3–12.5 MHz higher, at α = −220.5 to −221.0 MHz.

The charge matrix element |⟨0|n|1⟩| at EJ/EC = 88.8 is 1.265. The large-EJ/EC asymptotic (EJ / 32 EC)^(1/4) is 1.291.

The next spacing is not the same α. f23 − f12 is about −235 to −238 MHz. A leakage-reduction pulse that assumes equal Kerr spacings is already wrong by that amount.

## The CZ point moves when the qubit is tuned

The gate brings |11⟩ into resonance with |20⟩ of the higher-frequency qubit:

f_tuned + α(EJ_tuned) = f_partner.

α depends on EJ, so the idle value is the wrong one to subtract: |α| grows as the measure qubit is tuned down. For EC/h = 200 MHz, an idle frequency of 5.90 GHz and a 5.20 GHz partner, idle α is −217.4 MHz and the idle-α estimate of the interaction is 5.4174 GHz. The solved point is 5.4192 GHz, where α is −219.2 MHz, so the idle estimate sits about 2 MHz low. (Standalone check below. The library does the same solve in `f01_for_11_20`, and `test_qpuforge transmon` prints this example.)

How qpuforge treats this:
- **v1.1** does not use the idle-α estimate. The CZ interaction flux is one of the three knobs the pulse-level calibration finds by simulating the three-body Hamiltonian. The calibrated interaction points are 5.32–5.50 GHz, which is an excursion of 309–665 MHz on the measure qubits (`alpha1_v11_review.json`, `cz_calibration.per_pair`).
- **v1.0** uses the solved point. Since v0.4 `analyze()` takes the interaction point from `f01_for_11_20` (`include/qpuforge/device.hpp`), at the anharmonicity the moving qubit has there; v0.3 used the idle-α estimate `lo.f01 - hi.alpha`. For the 24 `alpha1()` links the interaction points are 5.33–5.51 GHz and the excursions 298.9–657.5 MHz. The solved points lie 1.1–2.5 MHz above the idle-α estimate, which gives excursions of 300–660 MHz.

## Alpha-1 v1.0: the coupler has to be a circuit (added in PR #2)

The previous Alpha-1 file parked every coupler at a typed residual of 5 kHz. That number is not what a capacitive three-transmon circuit produces.

With a direct exchange of 4 MHz and a coupler–qubit coupling of 70 MHz, the bias that cancels the exchange (coupler near 4.23 GHz) still leaves −119 kHz of ZZ. Cutting the direct capacitance to 1 MHz, and setting g_qc = 40 MHz with EC/h = 200 MHz, moves the exchange null to 3.88 GHz on a 5.20/5.90 GHz pair and leaves −7 kHz. Across the 24 Alpha-1 links the nulls lie between 3.787 and 3.967 GHz and the idle ZZ lies between −4.3 and −13.8 kHz. Each coupler junction is sized so that null is its flux sweet spot. The 23 ns gate still takes g_max = 40 MHz from a separate on-state bias. This diagonalisation is three levels on each island. It is not an electromagnetic extract of the AF17.V11 polygons.

## Alpha-1 v1.0: flux noise during the 23 ns gate (added in PR #2)

Idle dephasing uses the sweet spot, where df/dΦ is zero and only the second derivative of frequency survives. The CZ does not happen there. On a 5.90 GHz qubit tuned onto a 5.20 GHz partner the gate sits at 0.177 Φ₀, where df/dΦ is −5.5 GHz/Φ₀.

The same 1/f spectrum as the idle model (A = 2 µΦ₀, S_Φ(f) = A²/|f|) then gives a short-time infidelity of 3.2×10⁻⁵ for that pair, and 2.0–4.6×10⁻⁵ across Alpha-1. Phase variance grows as t², so treating the dephasing rate as Markovian and multiplying by 23 ns would overstate the error by about two orders of magnitude. The reported CZ error is 2.43×10⁻³, still the 2.0×10⁻³ control allocation plus decoherence plus this term. The v1.0 logical-error table below was not recomputed with this term, the coupler model or the 1 MHz coupling. The term is opt-in: `ChipDesign::cz_flux_noise` (design JSON `gates.cz_flux_noise`; files without the key follow `tunable_coupler.model_coupler`) adds it to `e_2q`, and only `alpha1()` (v1.0) turns it on. `analyze()` still reports it for every design (the `flux_e` column) but leaves it out of `e_2q` elsewhere, so Gen-0 keeps its pre-PR #2 mean `e_2q` of 4.09e-3.

These two sections describe the v1.0 abstract model only. v1.1 does not use them: its idle ZZ and its flux-noise terms come from the pulse-level calibration below. Its design file keeps the 4 MHz bare direct coupling it was generated with (`alpha1_with_coupling`; Beta-1 and Gamma-1 inherit it), whereas `alpha1()` has used 1 MHz since PR #2. The pulse-level calibration sets the per-pair g12 (9.4–9.8 MHz) itself.

## Alpha-1 v1.1: the CZ gate

Alpha-1 is a rotated distance-3 surface code: 9 data qubits and 8 measure qubits, 24 tunable couplers, three readout feedlines. Data qubits idle at 5.11–5.29 GHz and measure qubits at 5.81–5.99 GHz, all at α ≈ −220 MHz. Readout is 500 ns with κ = 4 MHz, and the active reset is 200 ns (`ALPHA-1-v1.1.md` §2).

In v1.1 the measure qubit (the higher frequency) is the mover. Its flux pulse brings |11⟩ into resonance with |20⟩, and at the same time the coupler swings from 8.5 GHz down to 6.03–6.14 GHz. One full |11⟩ → |20⟩ → |11⟩ cycle gives a conditional phase of π. The simulation propagates the full qubit–coupler–qubit charge Hamiltonian, truncated to 4/4/4 levels, in piecewise-exact steps of 0.1 ns, including the single-qubit rotating frame so that idle ZZ counts as an error.

| Quantity over the 24 pairs | Value |
|---|---|
| Duration | 35.6–49.2 ns (mean 43.7 ns) |
| Coherent error (1 − F_avg, closed system) | 7.9e-5 mean, 2.1e-5 to 3.6e-4 per pair |
| **Leakage per CZ (moving qubit)** | **4.8e-5 mean, 1.2e-4 max (Q8–Q10)**, within the 2e-4 gate for every pair |
| 1/f flux noise at the excursion (local + quasi-static) | 3.2e-5 mean, 1.3e-5 to 5.4e-5 per pair |
| Residual idle ZZ at park | 3.4–4.5 kHz (mean 3.9 kHz) |
| Residual control allowance (assumed, not simulated) | 1.0e-3 |
| Total CZ error incl. T1/Tφ during the gate | 1.86e-3 mean |

Source: `alpha1_v11_review.json` (`cz_calibration`, produced by `alpha1_v11_review 100000`; see `results/PROVENANCE.md`). The 1.0e-3 control allowance dominates the total and is an assumption.

**Leakage is modelled in the CZ.** Since v0.3 qpuforge has a pulse-level, leakage-modelled CZ (`pulse.hpp`). The leakage figure above is the population that leaves the computational subspace during the calibrated pulse, averaged over the four computational input states (1 − Tr(M†M)/4 in `simulate_cz`, `pulse.hpp`); the error-correction model injects it on the moving qubit. Under CZ model version 2, Q8–Q11 leaked 2.33e-4 (the Lc = 4 recalibration; the earlier Lc = 3 truncation was not converged). Its plain calibration sits on a short-flat-top branch; version 4 re-calibrates a pair that fails the 2e-4 limit with leakage weighted more strongly and keeps the result only if it meets the limit without a worse 1 − F_avg. Q8–Q11 then leaks 1.6e-5, with 1 − F_avg 3.6e-4 and a 3.5 ns longer gate (`docs/design/cz-leakage.md`).

**Leakage in the error-correction model.** The Pauli-frame QEC simulator (`qec.hpp`) carries a per-qubit leakage flag. It models CZ leakage injection (with the per-pair values above), transport at 5 %, seepage at 2/T1, a leakage-removing reset of the measure qubits, and optional data-qubit leakage removal (DQLR). The decoder does not use the leakage: there is no 3-state readout and no erasure flag, so the model is leakage-modelled, not leakage-aware.

## Error correction

### v1.1 (pulse CZ, leakage, leakage-removing reset)

Distance-3 memory, 5–25 rounds, Pauli-frame Monte Carlo with exact matching, 100,000 shots per round count (200,000 for the v1.1 memory-Z headline) (`alpha1_v11_review.json`, `alpha1_v11_review.txt`):

| Model | Logical error per round | Round |
|---|---:|---:|
| v1.0 abstract model, run inside the v1.1 review (no CZ leakage), memory-Z | 0.454 % ± 0.006 % | 847 ns |
| **v1.1 (pulse CZ + leakage + leakage-removing reset), memory-Z** | **0.377 % ± 0.004 %** | **928 ns** |
| v1.1, memory-X (\|+⟩ logical state, phase flips) | 0.360 % ± 0.005 % | 928 ns |

The v1.0 figure here is not the 0.45 % of the v1.0 review below: they are separate runs with different settings, and this one uses the current `alpha1()` with the PR #2 coupler model. The two v1.1 bases differ by 0.017 %, 2.7 standard errors. The logical lifetime is 246 µs against a best physical T1 of 75.6 µs. At distance 5 with the same physics the review reports 0.123 % per round against 0.364 % at distance 3 in that comparison (a fit over 5–15 rounds), a suppression factor Λ = 2.95 ± 0.18.

Leakage handling, from separate Monte Carlo runs with their own seeds (hence 0.366 % here against the 0.377 % headline):

| Scheme | eps / round | P_L (50 rounds) | leaked data / measure (time-averaged over the 50 rounds) |
|---|---:|---:|---|
| Leakage not modelled | 0.366 % | 0.150 | — |
| No leakage removal | 0.395 % | 0.163 | 1.0e-3 / 1.4e-3 |
| **Leakage-removing reset (v1.1)** | **0.366 %** | **0.153** | 4.9e-5 / 9.7e-5 |
| LRR + data-qubit leakage removal | 0.402 % | 0.168 | 5.5e-6 / 8.8e-5 |

Leave-one-out error budget for v1.1 (change in eps): T1 −80 %, CZ/1q control −40 %, Tφ −32 %, readout −18 %, reset −13 %. Leakage (−5 %), thermal (−3 %) and ZZ (+1 %) are within a few per cent, at the level of the Monte Carlo noise. The full table and the sensitivity runs are in `alpha1_v11_review.txt`.

### v1.0 (abstract couplers)

The committed v1.0 review (`designs/alpha1/alpha1_review.txt`, 24 November, before PR #2: 4 MHz coupling, typed 5 kHz idle ZZ, no CZ flux term) is a Pauli-frame memory experiment on the abstract model. It is kept as a historical result and was not regenerated.

| Step | Logical error per round | Round |
|---|---:|---:|
| Fixed coupling, g = 3 MHz, 800 ns readout | 1.23 % | 1907 ns |
| Tunable couplers (idle ZZ 5 kHz) | 0.54 % | 1253 ns |
| Fast Purcell-filtered readout | 0.54 % | 953 ns |
| Active reset (Alpha-1 v1.0) | 0.45 % | 853 ns |

At distance 5 the same physics gives 0.18 % per round, Λ = 2.61. A Monte Carlo of junction resistance gives a collision-free yield of 98 % at 1 % and 2 % spread and 91 % at 4 %. The logical lifetime is about 187 µs against a best physical T1 of 76 µs. For the 24 Nov file, `qpuforge analyze` gave 24 couplers, 0 collisions, mean T1 70.8 µs, mean CZ error 2.40e-3 and a 23.3 ns CZ. On `main` (`designs/alpha1/alpha1.json` after PR #2) it gives mean T1 70.5 µs, mean CZ error 2.43e-3 and max |idle ZZ| 14.0 kHz, still with 0 collisions.

The committed v1.0 review was not regenerated for v0.3 or v0.4, but it is bit-reproducible under the conditions that produced it: the binary `examples/alpha1_review` as committed in `3fe4d23`, run with `OMP_NUM_THREADS=2`. That binary predates the thread-count-independent Monte Carlo: its `run_memory_frame` gave each OpenMP thread its own RNG stream, so its results depend on the thread count. It was also a g++ 13.3 `-march=native` build with FMA contraction, which slightly shifts the Gen-0 noise model. Reproduced that way, the JSON matches byte for byte (see `results/PROVENANCE.md`, v1.0 row, on `main` since PR #6; the recipe applies to the `3fe4d23` checkout only). A rebuild from v0.3 source at 8 threads still reproduces the headline within Monte Carlo and those two effects (Alpha-1 0.45–0.46 % per round at 853 ns, Λ = 2.60, yield 98/98/91 %), but not digit for digit: the earlier steps differ by up to 0.0003 per round (for example Gen-0 1.20 % against 1.23 %) and the lifetime is about 190 µs. Treat the v1.0 table as a sample, not as exact. Neither effect applies to the v1.1 numbers above: since v0.4 the frame Monte Carlo seeds every shot separately, so its results are bit-identical at any `OMP_NUM_THREADS`, and CMake builds every target with `-ffp-contract=off` (PR #6), so `-march` does not change them. Every other committed Monte Carlo output (the Alpha-1 v1.1, Beta-1 and Gamma-1 reviews, the design study and the test log) was regenerated that way; the CZ calibrations in `results/cz_cache.json` (CZ model version 4) and the `AF17_v11.gds` mask were regenerated with them, as `results/PROVENANCE.md` records.

### Limits that apply to both

These are emulator outputs. The reviews list what is missing: correlated and time-varying errors (quasiparticles, cosmic rays, TLS drift), a decoder that uses leakage information, and a real electromagnetic extract. On published hardware such effects commonly cost a factor of about 1.5–3 in logical error, so the numbers above are best cases of these models. The 1.0e-3 CZ control allowance is assumed.

## The AF17.V11 layout

The layout is generated for v1.1. The cells are synthesised against circuit targets with a quasi-static method-of-moments extraction, assembled into a 5.8 × 5.8 mm chip, and exported as GDSII (`designs/alpha1-v11/AF17_v11.gds`, cell `AF17_V11`). It is not a photograph and not a measured chip.

**The mask is DRC- and connectivity-clean in the model** (2D LVS-lite connectivity with modelled airbridges). Airbridges, EM and inductance are not verified, and no foundry DRC deck has been run.
- Placement and island-spacing checks: clean (0 problems).
- Merged-etch minimum width and space of 4 µm, projection metric: **0 + 0 violations** (34 + 202 in 0.4). BF-49 rev C and GF-97: 0 + 0.
- Metal connectivity (gated): clean. Each qubit SQUID bridges its stub to a ground lead cut out of the qubit's inner corner; each qubit flux line is shorted next to its SQUID; 122 ground straps (modelled airbridges) tie the enclosed ground (qubit inner corners, coupler–bus strips, plaquettes) to the plane, none over a qubit. No floating ground remains (143 fragments before). A strap's capacitance to the coupler or bus it crosses is not in the extraction.
- Control lines (reported, not gated; **do not hold**: the review's mask-window extraction finds the XY bump coupled ~1.4–1.6 fF to neighbouring coupler pads, and the flux line's mutual inductance to the SQUID is ~0.01 pH): XY 0.03 fF and flux line 0.03 fF effective per qubit; T1 limits ≥ 1.56 ms (XY) and 1.47 ms (flux), at least 22× the design T1. The flux lines' mutual inductance to the SQUIDs is not modelled. Design note: `docs/design/mask-connectivity.md`.

Post-layout circuit check (`alpha1_v11_review.txt`, section 7), from the extracted cells:
- α = −220 MHz on both qubit types and −130 MHz on the coupler; g_qc = 168.0 MHz; readout g = 154.9 MHz against a design value of 155.8 MHz.
- The direct coupling is g12 = 9.72 MHz against the 9.54 MHz needed for an exchange-free park at the coupler sweet spot, which would leave a 0.37 MHz gap. Parking each coupler at 8.398 GHz, 102 MHz below its sweet spot, nulls the exchange, and ZZ is then −4.4 kHz.

The extraction is quasi-static on a thick substrate. It is not a full-wave extract, and a full-wave cross-check is the next design step before tape-out.

### Figures

Previews rendered from `AF17_v11.gds` by the independent reader `tools/gds_preview.py` are committed with the design: [full chip](../designs/alpha1-v11/AF17_v11_layout.png) and [detail](../designs/alpha1-v11/AF17_v11_detail.png), with vector versions of the same two views ([full chip](../designs/alpha1-v11/AF17_v11_layout.svg), [detail](../designs/alpha1-v11/AF17_v11_detail.svg)). Vector redraws of the committed `AF17_v11.gds` over the three views of the earlier figures are in `docs/from-committed-layout/`: [full chip](from-committed-layout/af17-full.svg), [Q5, Q7, Q8 and Q10](from-committed-layout/af17-zoom.svg) and [Q5](from-committed-layout/af17-q5.svg). The ground straps are drawn in violet. Diagonal gaps are stair-stepped in all of them because the mask itself is drawn from rectangles (`AF17_v11.gds` has no non-Manhattan edges), not because of the rendering.

The three figures in the earlier version of this note (`docs/af17-full.png`, `af17-zoom.png`, `af17-q5.png` from PR #2, and their bitmap SVG wrappers `docs/af17-full.svg`, `af17-zoom.svg`, `af17-q5.svg`) are an earlier drawing and are not used here. They do not match `AF17_v11.gds`: the redraws above cover the same views, and the bump pad and junction positions and the cut-outs around Q5 differ, so they do not show the layout this note describes.

## Reproducing the numbers in this note

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
ctest --test-dir build -LE smoke                        # the test groups' checks are logged in results/test_output.txt
./build/test_qpuforge layout_presets                    # AF-17, BF-49, GF-97: DRC 0+0 and metal connectivity clean
./build/test_qpuforge designs_v04                       # Alpha-1 interaction points: 1.14 to 2.54 MHz above the idle-alpha estimate
./build/qpuforge new build/alpha1.json --preset alpha1  # writes and analyses alpha1(): charge dispersions, CZ points 5.33-5.51 GHz
./build/qpuforge analyze designs/alpha1/alpha1.json     # v1.0 (PR #2 model): 2.43e-3 CZ error
./build/qpuforge analyze designs/alpha1-v11/alpha1_v11.json
```

The v1.1 CZ, leakage and error-correction figures are in `designs/alpha1-v11/alpha1_v11_review.json`; `results/PROVENANCE.md` has the exact command (`alpha1_v11_review 100000`), and the JSON records the source commit and the compiler. The CZ calibrations come from the committed `results/cz_cache.json`. A rerun at any `OMP_NUM_THREADS` gives the same numbers bit for bit (`ctest -R determinism` checks the thread-count independence), and the `-ffp-contract=off` build makes them independent of `-march` as well (`ctest -R fp_contract`). The Beta-1 and Gamma-1 DRC and connectivity figures are in `designs/beta1/beta1_review.json` and `designs/gamma1/gamma1_review.json`; `beta1_review --smoke --out DIR` prints the Beta-1 line in a few minutes with the committed `results/cz_cache.json`.

The anharmonicity table and the 11–20 example come from this standalone program, which uses only header-only qpuforge (`g++ -std=c++17 -O2 -Iinclude check.cpp`):

```cpp
#include "qpuforge/device.hpp"
#include <cstdio>
using namespace qpuforge;
int main() {
    std::printf("E_J(Rn=10 kOhm, gap 180 ueV) = %.3f GHz\n", EJ_from_Rn_GHz(10.0, 180.0));
    for (double f01 : {5.11, 5.29, 5.81, 5.99}) {
        const auto [EJ, EC] = solve_EJ_EC(f01, -0.220);
        const auto T = transmon_levels(EJ, EC, 0.0, 6);
        std::printf("f01 %.2f  EJ/EC %.1f  EC %.0f MHz  alpha %.1f MHz  f23-f12 %.1f MHz  |n01| %.3f  (EJ/32EC)^(1/4) %.3f\n",
                    f01, EJ / EC, EC * 1e3, T.alpha() * 1e3, (T.E[3] - 2 * T.E[2] + T.E[1]) * 1e3,
                    std::abs(T.n(0, 1)), std::pow(EJ / (32 * EC), 0.25));
    }
    // |11>-|20> point: f_tuned + alpha(EJ_tuned) = f_partner, against the idle-alpha estimate
    const double EC = 0.200, f_idle = 5.90, f_partner = 5.20;
    const auto idle = transmon_fa(solve_EJ_for_f01(f_idle, EC), EC);
    double lo = std::log(0.2), hi = std::log(4000.0);
    for (int b = 0; b < 80; ++b) {
        const double m = 0.5 * (lo + hi);
        const auto fa = transmon_fa(std::exp(m), EC);
        ((fa.f01 + fa.alpha) < f_partner ? lo : hi) = m;
    }
    const auto tuned = transmon_fa(std::exp(0.5 * (lo + hi)), EC);
    std::printf("idle alpha %.1f MHz; idle-alpha point %.4f GHz; solved point %.4f GHz (alpha there %.1f MHz)\n",
                idle.alpha * 1e3, f_partner - idle.alpha, tuned.f01, tuned.alpha * 1e3);
}
```

## Licence

Software, design files and this note belong to Secured Cyber Solutions Ltd. All rights reserved. Third-party interatomic-potential parameters used by atomforge, if present alongside this repository, keep their own licence.
