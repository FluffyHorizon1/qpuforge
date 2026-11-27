# Transmon anharmonicity and the Alpha-1 processor

Finnian Bond
Secured Cyber Solutions Ltd
First written 25 November 2026 (the version on `main`, `docs/WHITEPAPER.md`, from PR #2). Revised for qpuforge v0.3 on 27 November 2026.

> **Merge proposal, not yet the whitepaper.** This file combines the PR #2 note on `main` with the v0.3 revision on
> `docs/whitepaper-v0.3` (`6671993`), checked against `main` at `78b9f8b`. Notes marked **[Cn]** point to the content
> choices in `docs/notes/whitepaper-merge-choices.md`. Once the author has chosen, this text replaces
> `docs/WHITEPAPER.md` and the [Cn] notes are deleted.

## Abstract

Alpha-1 (AF-17) is a 17-qubit distance-3 surface-code design. It has not been fabricated. This note records the junction physics the design uses, why the |11⟩–|20⟩ frequency has to be solved at the tuned Josephson energy, the Alpha-1 coupler-transmon and CZ flux-noise models added on 25 November, and the AF17.V11 metal layout. Every error rate quoted here is from the qpuforge emulator in this repository, not from a device. **[C1]**

## Status and known issues

Alpha-1 has not been fabricated. Every number in this note is a qpuforge simulation result, not a measurement. v0.3 ships with two known issues, and both concern the AF-17 design this note describes:

1. **The AF-17 mask fails design-rule checks.** The AF17.V11 first-cut layout passes the placement and island-spacing checks (0 problems) but **fails the merged-etch width/space DRC at 4 µm / 4 µm (projection metric) with 34 width + 202 space violations**. The Beta-1 layout (BF-49) fails it too, with **98 width + 602 space violations**. Neither mask is fabrication-ready.
2. **Worst-pair CZ leakage is over the gate.** In Alpha-1 v1.1 the pair Q8–Q11 leaks **2.33e-4 per CZ, above the 2e-4 simulation-gate limit** (`ROADMAP.md`). The 24-pair mean is 5.7e-5. The leakage was 1.3e-4 before the v0.3 recalibration at coupler truncation Lc = 4.

Sources: `RELEASE_NOTES.md` (Known issues), `designs/alpha1-v11/alpha1_v11_review.json` (`layout.drc`, `cz_calibration`), `designs/beta1/beta1_review.json` (`layout.drc`), `results/test_output.txt` (layout test).

## Two designs, both called AF-17

The note uses "v1.0" and "v1.1" for two models of the same 17-qubit chip, and gives their numbers in separate tables. They are never mixed.

| | Alpha-1 v1.0 | Alpha-1 v1.1 |
|---|---|---|
| Couplers | abstract: g_max = 40 MHz. The 24 Nov review used 4 MHz direct coupling and a typed 5 kHz residual ZZ; since PR #2 `alpha1()` uses 1 MHz and a coupler transmon at the exchange null (idle ZZ −4 to −14 kHz, computed) | physical coupler transmon per pair (f_max 8.5 GHz, α_c −130 MHz, g_qc 160 MHz) |
| CZ | 23.3 ns, 2e-3 assumed control error | pulse-level, 3-knob calibration, 35.6–49.2 ns |
| Leakage | folded into the control error | simulated in the CZ (pulse level) and in the QEC frame model |
| Layout | floor plan (SVG), no mask | MoM-synthesised cells, GDSII, DRC |
| Files | `designs/alpha1/` (spec `ALPHA-1.md`; review output from 24 Nov, pre-PR #2, historical: see `designs/alpha1/README.md`) | `designs/alpha1-v11/` (spec `ALPHA-1-v1.1.md`, review regenerated for v0.3, `results/PROVENANCE.md`) |

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

EJ/EC, EC, α and f23 − f12 are from the standalone check at the end of this note. The charge dispersions are the `disp_kHz` column of `qpuforge analyze designs/alpha1/alpha1.json`. The −220 MHz written into the design file is this exact α. Rn is what moves f01. Charge dispersion of f01 is hertz, not megahertz: these ratios are in the transmon regime.

The charge matrix element |⟨0|n|1⟩| at EJ/EC = 88.8 is 1.265. The large-EJ/EC asymptotic (EJ / 32 EC)^(1/4) is 1.291.

The next spacing is not the same α. f23 − f12 is about −235 to −238 MHz. A leakage-reduction pulse that assumes equal Kerr spacings is already wrong by that amount.

## The CZ point moves when the qubit is tuned

The gate brings |11⟩ into resonance with |20⟩ of the higher-frequency qubit:

f_tuned + α(EJ_tuned) = f_partner.

α depends on EJ, so the idle value is the wrong one to subtract: |α| grows as the measure qubit is tuned down. For EC/h = 200 MHz, an idle frequency of 5.90 GHz and a 5.20 GHz partner, idle α is −217.4 MHz and the idle-α estimate of the interaction is 5.4174 GHz. The solved point is 5.4192 GHz, where α is −219.2 MHz, so the idle estimate sits about 2 MHz low. (Standalone check below; this solve is not a function in the library.)

How v0.3 treats this:
- **v1.1** does not use the idle-α estimate. The CZ interaction flux is one of the three knobs the pulse-level calibration finds by simulating the three-body Hamiltonian. The calibrated interaction points are 5.32–5.50 GHz, which is an excursion of 309–665 MHz on the measure qubits (`alpha1_v11_review.json`, `cz_calibration.per_pair`).
- **v1.0** uses the idle-α estimate in `analyze()` (`include/qpuforge/device.hpp`, `cz_interaction_GHz = lo.f01 - hi.alpha`). For the 24 `alpha1()` links that gives interaction points of 5.33–5.51 GHz and excursions of 300–660 MHz. Those figures carry the ~2 MHz idle-α offset; the tuned-α solve gives 298.9–657.5 MHz. The PR #2 version of this note said `analyze()` already does the solve; it does not on `main`. **[C3]**

## Alpha-1 v1.0: the coupler has to be a circuit (added in PR #2)

The previous Alpha-1 file parked every coupler at a typed residual of 5 kHz. That number is not what a capacitive three-transmon circuit produces.

With a direct exchange of 4 MHz and a coupler–qubit coupling of 70 MHz, the bias that cancels the exchange (coupler near 4.23 GHz) still leaves −119 kHz of ZZ. Cutting the direct capacitance to 1 MHz, and setting g_qc = 40 MHz with EC/h = 200 MHz, moves the exchange null to 3.88 GHz on a 5.20/5.90 GHz pair and leaves −7 kHz. Across the 24 Alpha-1 links the nulls lie between 3.787 and 3.967 GHz and the idle ZZ lies between −4.3 and −13.8 kHz. Each coupler junction is sized so that null is its flux sweet spot. The 23 ns gate still takes g_max = 40 MHz from a separate on-state bias. This diagonalisation is three levels on each island. It is not an electromagnetic extract of the AF17.V11 polygons.

## Alpha-1 v1.0: flux noise during the 23 ns gate (added in PR #2)

Idle dephasing uses the sweet spot, where df/dΦ is zero and only the second derivative of frequency survives. The CZ does not happen there. On a 5.90 GHz qubit tuned onto a 5.20 GHz partner the gate sits at 0.177 Φ₀, where df/dΦ is −5.5 GHz/Φ₀.

The same 1/f spectrum as the idle model (A = 2 µΦ₀, S_Φ(f) = A²/|f|) then gives a short-time infidelity of 3.2×10⁻⁵ for that pair, and 2.0–4.6×10⁻⁵ across Alpha-1. Phase variance grows as t², so treating the dephasing rate as Markovian and multiplying by 23 ns would overstate the error by about two orders of magnitude. The reported CZ error is 2.43×10⁻³, still the 2.0×10⁻³ control allocation plus decoherence plus this term. The v1.0 logical-error table below was not recomputed with this term, the coupler model or the 1 MHz coupling. On `main` the term is also added to `e_2q` of every other design whose higher qubit is tunable (Gen-0 mean `e_2q` 4.09e-3 → 5.98e-3), not just Alpha-1. **[C5]**

These two sections describe the v1.0 abstract model only. v1.1 does not use them: its idle ZZ and its flux-noise terms come from the pulse-level calibration below, and it was generated with 4 MHz direct coupling. **[C6]**

## Alpha-1 v1.1: the CZ gate

Alpha-1 is a rotated distance-3 surface code: 9 data qubits and 8 measure qubits, 24 tunable couplers, three readout feedlines. Data qubits idle at 5.11–5.29 GHz and measure qubits at 5.81–5.99 GHz, all at α ≈ −220 MHz. Readout is 500 ns with κ = 4 MHz, and the active reset is 200 ns (`ALPHA-1-v1.1.md` §2).

In v1.1 the measure qubit (the higher frequency) is the mover. Its flux pulse brings |11⟩ into resonance with |20⟩, and at the same time the coupler swings from 8.5 GHz down to 6.03–6.13 GHz. One full |11⟩ → |20⟩ → |11⟩ cycle gives a conditional phase of π. The simulation propagates the full qubit–coupler–qubit charge Hamiltonian, truncated to 4/4/4 levels, in piecewise-exact steps of 0.1 ns, including the single-qubit rotating frame so that idle ZZ counts as an error.

| Quantity over the 24 pairs | Value |
|---|---|
| Duration | 35.6–49.2 ns (mean 43.6 ns) |
| Coherent error (1 − F_avg, closed system) | 8.1e-5 mean, 2.1e-5 to 4.1e-4 per pair |
| **Leakage per CZ (moving qubit)** | **5.7e-5 mean, 2.33e-4 max (Q8–Q11, over the 2e-4 gate)**; the next-worst pair is 1.2e-4 (Q8–Q10) |
| 1/f flux noise at the excursion (local + quasi-static) | 3.1e-5 mean, 1.1e-5 to 5.4e-5 per pair |
| Residual idle ZZ at park | 3.4–4.5 kHz (mean 3.9 kHz) |
| Residual control allowance (assumed, not simulated) | 1.0e-3 |
| Total CZ error incl. T1/Tφ during the gate | 1.86e-3 mean |

Source: `alpha1_v11_review.json` (`cz_calibration`, produced by `alpha1_v11_review 100000`; see `results/PROVENANCE.md`). The 1.0e-3 control allowance dominates the total and is an assumption.

**Leakage is modelled in the CZ.** v0.3 has a pulse-level, leakage-modelled CZ (`pulse.hpp`). The leakage figure above is the |2⟩ population of the moving qubit at the end of the calibrated pulse. The 2.33e-4 on Q8–Q11 is the Lc = 4 recalibration. The earlier Lc = 3 truncation was not converged (it left a 4e-3 rad conditional-phase error), and `CHANGELOG.md` records the change.

**Leakage in the error-correction model.** The Pauli-frame QEC simulator (`qec.hpp`) carries a per-qubit leakage flag. It models CZ leakage injection (with the per-pair values above), transport at 5 %, seepage at 2/T1, a leakage-removing reset of the measure qubits, and optional data-qubit leakage removal (DQLR). The decoder does not use the leakage: there is no 3-state readout and no erasure flag, so the model is leakage-modelled, not leakage-aware.

## Error correction

### v1.1 (pulse CZ, leakage, leakage-removing reset)

Distance-3 memory-Z, 5–25 rounds, Pauli-frame Monte Carlo with exact matching, 100,000 shots (`alpha1_v11_review.json`, `alpha1_v11_review.txt`):

| Model | Logical error per round | Round |
|---|---:|---:|
| v1.0 abstract model, run inside the v1.1 review (no leakage) | 0.440 % ± 0.005 % | 847 ns |
| **v1.1 (pulse CZ + leakage + leakage-removing reset)** | **0.368 % ± 0.003 %** | **928 ns** |

The v1.0 figure here differs from the 0.45 % of the v1.0 review below because they are separate runs with different settings. The logical lifetime is 252 µs against a best physical T1 of 75.6 µs. At distance 5 with the same physics the review reports 0.123 % per round against 0.383 % at distance 3 in that comparison, a suppression factor Λ = 3.12.

Leakage handling, from separate Monte Carlo runs with their own seeds (hence 0.387 % here against the 0.368 % headline):

| Scheme | eps / round | P_L (50 rounds) | leaked data / measure after 50 rounds |
|---|---:|---:|---|
| Leakage not modelled | 0.374 % | 0.152 | — |
| No leakage removal | 0.406 % | 0.171 | 1.2e-3 / 1.6e-3 |
| **Leakage-removing reset (v1.1)** | **0.387 %** | **0.159** | 8.1e-5 / 1.2e-4 |
| LRR + data-qubit leakage removal | 0.414 % | 0.170 | 7.6e-6 / 1.1e-4 |

Leave-one-out error budget for v1.1 (change in eps): T1 −78 %, CZ/1q control −39 %, Tφ −27 %, readout −19 %, reset −9 %, leakage −1 %. Thermal and ZZ (+2 %) are within Monte Carlo noise. The full table and the sensitivity runs are in `alpha1_v11_review.txt`.

### v1.0 (abstract couplers)

The committed v1.0 review (`designs/alpha1/alpha1_review.txt`, 24 November, before PR #2: 4 MHz coupling, typed 5 kHz idle ZZ, no CZ flux term) is a Pauli-frame memory experiment on the abstract model. It is kept as a historical result and was not regenerated. **[C2]**

| Step | Logical error per round | Round |
|---|---:|---:|
| Fixed coupling, g = 3 MHz, 800 ns readout | 1.23 % | 1907 ns |
| Tunable couplers (idle ZZ 5 kHz) | 0.54 % | 1253 ns |
| Fast Purcell-filtered readout | 0.54 % | 953 ns |
| Active reset (Alpha-1 v1.0) | 0.45 % | 853 ns |

At distance 5 the same physics gives 0.18 % per round, Λ = 2.61. A Monte Carlo of junction resistance gives a collision-free yield of 98 % at 1 % and 2 % spread and 91 % at 4 %. The logical lifetime is about 187 µs against a best physical T1 of 76 µs. For the 24 Nov file, `qpuforge analyze` gave 24 couplers, 0 collisions, mean T1 70.8 µs, mean CZ error 2.40e-3 and a 23.3 ns CZ. On `main` (`designs/alpha1/alpha1.json` after PR #2) it gives mean T1 70.5 µs, mean CZ error 2.43e-3 and max |idle ZZ| 14.0 kHz, still with 0 collisions.

The committed v1.0 review was not regenerated for v0.3, but it is bit-reproducible under the conditions that produced it: the binary `examples/alpha1_review` as committed in `3fe4d23`, run with `OMP_NUM_THREADS=2`. `run_memory_frame` gives each OpenMP thread its own RNG stream, so results depend on the thread count. That binary was also a g++ 13.3 `-march=native` build with FMA contraction, which slightly shifts the Gen-0 noise model. Reproduced that way, the JSON matches byte for byte (see `results/PROVENANCE.md`, v1.0 row, on `main` since PR #6; the recipe applies to the `3fe4d23` checkout only). A rebuild from v0.3 source at 8 threads still reproduces the headline within Monte Carlo and those two effects (Alpha-1 0.45–0.46 % per round at 853 ns, Λ = 2.60, yield 98/98/91 %), but not digit for digit: the earlier steps differ by up to 0.0003 per round (for example Gen-0 1.20 % against 1.23 %) and the lifetime is about 190 µs. Treat the v1.0 table as a sample, not as exact. Making the Monte Carlo independent of thread count is proposed on `feat/deterministic-mc-on-main` **[C7]**; the `-ffp-contract=off` build change is on `main` (PR #6).

### Limits that apply to both

These are emulator outputs. The reviews list what is missing: correlated and time-varying errors (quasiparticles, cosmic rays, TLS drift), a decoder that uses leakage information, and a real electromagnetic extract. On published hardware such effects commonly cost a factor of about 1.5–3 in logical error, so the numbers above are best cases of these models. The 1.0e-3 CZ control allowance is assumed.

## The AF17.V11 layout

The layout is generated for v1.1. The cells are synthesised against circuit targets with a quasi-static method-of-moments extraction, assembled into a 5.8 × 5.8 mm chip, and exported as GDSII (`designs/alpha1-v11/AF17_v11.gds`, cell `AF17_V11`). It is not a photograph and not a measured chip.

**The mask is a first cut and is not DRC clean.**
- Placement and island-spacing checks: clean (0 problems).
- Merged-etch minimum width and space of 4 µm, projection metric: **34 width + 202 space violations. Not clean.**
- The same check on the Beta-1 layout BF-49: 98 width + 602 space violations.
- The Alpha simulation gate checks placement and spacing only (`ALPHA-1-v1.1.md` section 6), so the gate reads as passed on means while the mask DRC is not clean. `ROADMAP.md` requires the mask width/space DRC to be clean before fabrication.

Post-layout circuit check (`alpha1_v11_review.txt`, section 7), from the extracted cells:
- α = −215 MHz on both qubit types and −129 MHz on the coupler; g_qc = 168.4 MHz; readout g = 157.6 MHz against a design value of 155.8 MHz.
- The direct coupling is g12 = 9.94 MHz against the 9.54 MHz needed for an exchange-free park at the coupler sweet spot, which would leave a 0.77 MHz gap. Parking each coupler at 8.295 GHz, 205 MHz below its sweet spot, nulls the exchange, and ZZ is then −5.0 kHz.

The extraction is quasi-static on a thick substrate. It is not a full-wave extract, and a full-wave cross-check is the next design step before tape-out.

### Figures

Previews rendered from `AF17_v11.gds` by the independent reader `tools/gds_preview.py` are committed with the design: [full chip](../designs/alpha1-v11/AF17_v11_layout.png) and [detail](../designs/alpha1-v11/AF17_v11_detail.png). They are raster previews (the gaps look stair-stepped) and they do not mark the DRC violations. **[C4]**

The three figures in the earlier version of this note (`docs/af17-full.png`, `af17-zoom.png`, `af17-q5.png` from PR #2, now on `main` under `docs/`) are not used here. Rendering the v0.3 `AF17_v11.gds` over the same windows gives different geometry from those figures (bump pad and junction positions, and the cut-outs around Q5 differ), so they do not show the layout this note describes.

## Reproducing the numbers in this note

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
ctest --test-dir build -LE smoke                        # 21 entries, 186 checks in total on main 78b9f8b
./build/test_qpuforge layout                            # prints the AF-17 DRC line: 34 width + 202 space
./build/qpuforge analyze designs/alpha1/alpha1.json     # v1.0 (PR #2 model): charge dispersions, 300-660 MHz excursions, 2.43e-3 CZ error
./build/qpuforge analyze designs/alpha1-v11/alpha1_v11.json
```

The v1.1 CZ, leakage and error-correction figures are in `designs/alpha1-v11/alpha1_v11_review.json`; `results/PROVENANCE.md` has the exact command (`alpha1_v11_review 100000`, 8 OpenMP threads). The Beta-1 DRC figures are in `designs/beta1/beta1_review.json`; `beta1_review --smoke --out DIR` prints them in a few minutes with the committed `results/cz_cache.json`.

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
