# Design note: thread-count-independent Monte Carlo

Status: **adopted in 0.4.0.** The prototype below (on `feat/deterministic-mc`, rebased as `stack/deterministic-mc-on-q13`)
was ported onto the 0.4 `run_memory_frame`, which adds `FrameRunOptions`, the decode-method choice and the syndrome and
decoder statistics. They follow the same rule: the integer counts (errors, detections per round, defect counts) are
exact in any order, and the leaked-time sums go through the fixed block order. `decode_us_per_shot` is a
wall-clock measurement and stays machine-dependent. Every committed review output was regenerated with it
(`results/PROVENANCE.md`), so `openmp_threads` in their provenance is informational. The rest of this note is the
design record as written for the prototype; its check counts refer to that branch.

## Problem

`run_memory_frame` (`include/qpuforge/qec.hpp`), the Pauli-frame Monte Carlo behind every logical-error number in
the Alpha-1 v1.0, Alpha-1 v1.1 and Beta-1 reviews, gave each OpenMP thread its own xoshiro256 stream seeded from
`(seed, thread id)` and dealt shots out with `s = tid; s += nthreads`. The sample is therefore a function of
`OMP_NUM_THREADS`. That is why the committed Alpha-1 v1.0 review reproduced only with 2 threads (plus the FMA build;
see `results/PROVENANCE.md`), and why every committed review records its thread count.

Measured spread of the old code (Alpha-1 v1.0 ladder, Gen-0 step, 100000 shots per round count, seeds 11 + rounds):

| build | 1 | 2 | 3 | 4 | 6 | 8 | 12 | 16 threads |
|---|---|---|---|---|---|---|---|---|
| committed v1.0 binary (FMA build) | 0.012151 | **0.012323** | | 0.012306 | | 0.012278 | | 0.012108 |
| source, g++ 14.2 `-O3 -ffp-contract=off` | 0.012002 | 0.012281 | 0.012081 | 0.012044 | 0.011844 | 0.012002 | 0.011904 | 0.012106 |

The source-build column has mean 0.012033 and standard deviation 0.000133; the reported `eps_err` is 0.000123.

A second, smaller source of non-determinism: the trajectory emulator's sampled XEB (`run()` in `emulator.hpp`)
added doubles into per-thread accumulators under `schedule(dynamic, 1)`. That sum depends on which thread ran which
trajectory, so it could differ in the last bits between thread counts and even between runs with the same count.

## Goal

The result of every Monte Carlo routine is a function of its inputs and seed only: bit-identical for any
`OMP_NUM_THREADS`, any OpenMP schedule and any run. The cost should be negligible, and API signatures stay the same.

## Options

### A. Counter-based per-trial seeding (chosen)

Each shot gets its own generator state derived from `(seed, shot index)`. Which thread runs it, and in what order,
cannot matter.

- Candidates for the derivation: Philox-4x32/4x64 (Random123) used directly as the generator, or SplitMix64 to seed
  the existing xoshiro256 per shot.
- Philox is a true counter-based generator: no state, and any draw is addressable as `(key, counter)`. Adopting it
  would mean replacing the generator at every draw site, and it needs a portable 64x64->128 multiply (`__int128` or a
  manual version). It gives nothing extra here: a shot draws at most a few thousand numbers, and we only need shots
  to be independent and addressable.
- **Chosen: SplitMix64-keyed xoshiro256 per shot.** `frame_shot_seed(seed, shot) = mix(seed ^ c) + shot * 4γ`,
  where `mix` is the SplitMix64 finaliser (a bijection) and γ = 0x9E3779B97F4A7C15. The existing `Xoshiro256`
  constructor consumes four SplitMix64 steps of γ each, so shot *k* uses outputs 4k+1 … 4k+4 of the SplitMix64
  sequence keyed by `mix(seed ^ c)`. Shots of one frame never share seed material. Adjacent frame seeds (the reviews
  use `seed + rounds`, so seeds 13, 15, 17 … are common) are spread to unrelated keys by `mix`. A clash between two
  frames' shot ranges needs two keys within `4 * shots` steps of each other: probability about `n_frames^2 * shots /
  2^62`, which is below 1e-10 for every review. SplitMix64 seeding of xoshiro is the initialisation that xoshiro's
  authors recommend.
- The draw sites are unchanged (`rng.uniform()`), so the physics code needs no edits.
- The "frame" index is already part of the seed: callers pass `seed + rounds`, and different experiments use
  different base seeds. A separate frame argument would change the API without making anything more reproducible.
- Debugging bonus: one shot can be replayed in isolation from `(seed, shot)`.

### B. Jump-ahead per trial block

Split shots into fixed blocks of B. Block *b* uses the base stream advanced by *b* `jump()` calls (2^128 steps
each).

- This is equally thread-count independent, if blocks are tied to block indices and never to threads.
- Getting block *b* costs *b* jumps (256 `next()` calls plus XORs each), so block seeds must be precomputed serially.
- The **block size becomes part of the result**: changing B (for example to tune load balance) changes every number.
- A shot's draws depend on how many draws earlier shots in its block consumed, so single shots cannot be replayed.
- Rejected in favour of A, which has neither drawback and costs about the same.

### C. Static scheduling with a fixed chunk size — insufficient

`schedule(static, C)` makes the shot-to-thread assignment fixed for a **given** team size: chunk *i* goes to thread
`i mod nthreads`. With per-thread streams, the shot-to-stream mapping still changes when `OMP_NUM_THREADS` changes,
and that is exactly the bug. The old `s += nthreads` loop was already a static cyclic schedule, and it reproduced
only at the same thread count. Static scheduling also cannot guarantee the team size:
`OMP_DYNAMIC`, `OMP_THREAD_LIMIT`, nested regions and library callers of `omp_set_num_threads` can all change it.
Static scheduling makes C a hidden result parameter too. Static scheduling only works if each chunk carries its own
stream, and then it is option B.

## Reduction-order determinism

Per-shot seeding fixes the sample. The summation must also be order-independent:

- **Integer tallies** (logical errors, fired-detector count, leakage events) are exact in any order. `fired_sum` and
  `leak_ev` were `double`s with integer increments, which was exact below 2^53. They are now `long`, to make that
  explicit.
- **Floating-point tallies** (leaked-time sums `leak_d`, `leak_m`: `duration_ns * n_leaked`) are not associative.
  `#pragma omp parallel reduction(+:...)` combines per-thread partials in an unspecified order, which can vary between
  runtimes and team sizes. The prototype uses **fixed blocks of `kFrameBlock` = 1024 shots**. Each block sums its own
  shots in shot order into a `Tally`. After the parallel loop, the block tallies are added serially in block order.
  The block size affects only last-bit rounding of the leakage fractions, never the sample. Blocks run under
  `schedule(dynamic, 1)` for load balance, because shot cost varies with the number of defects to match.
- **Kahan summation** is not needed for determinism, because fixed order is enough. For accuracy, 1e5 to 1e6
  positive terms summed naively in blocks carry a relative error of about 1e-13, many orders below the Monte Carlo
  standard error. Kahan can be added to the final block-order sum at no measurable cost if a future tally needs it.
- **Emulator XEB:** per-trajectory partial sums, added in trajectory order. Shot counts are integers, and the
  per-thread count maps (integers) are left as they are.

Other parallel loops were audited and already deterministic:

- Trajectory loops in `emulator.hpp` (`run`, `run_records`, `run_ghz_fidelity`) seed each trajectory from
  `(seed, t)` and write per-index.
- `calibrate_chip_cz` (per-key results).
- `extract_capacitance` (per-row fill; element-wise LU update).
- `fab_monte_carlo` is serial.

The `std::` distributions used by the trajectory emulator are implementation-defined (libstdc++ against libc++), and
that is outside this change. The frame Monte Carlo uses only its own generator.

## Prototype

- `qec.hpp`: `splitmix64_mix`, `frame_shot_seed`, and `run_memory_frame` rewritten as above. The signature is
  unchanged.
- `emulator.hpp`: XEB per-trajectory partial sums.
- `tests/test_determinism.cpp`, plus two CTest entries labelled `unit`. Like `fp_contract`, they are build/behaviour
  guards and **not** part of the 101 validation checks:
  - `determinism`: runs the workloads at 1, 2, 3, 4 and 8 threads in one process (`omp_set_num_threads`) and
    requires bit-identical results.
  - `determinism_env`: `tests/compare_thread_digests.cmake` runs `test_determinism --digest` as separate
    processes under `OMP_NUM_THREADS` = 1, 2, 4, 8 and requires identical output.
  - Workloads: distance-3 memory with leakage on (gate and CZ leakage, transport), 2 and 5 rounds, 5000 shots (four
    full blocks plus a partial one), checked to the last bit with `%a`. Also a 6-qubit XEB run (300 trajectories × 20
    shots).
  - Against the old headers, the same test fails at every thread count other than 1 (checked).
- `test_qpuforge`: **101 checks passed, 0 failed** with the new code (g++ 14.2). The docs' "101 checks" claim stays
  true. Three printed Monte Carlo lines in `results/test_output.txt` would change in the last digits (frame detection
  rate 0.0477 -> 0.0475, leaked fractions, detection rates). The file is not regenerated here.

## Performance

Same machine, g++ 14.2, `-O3 -ffp-contract=off`, OpenMP. Times cover the Monte Carlo only, for the Gen-0 point
(5 round counts × 100000 shots):

| threads | old | new |
|---|---|---|
| 1 | 11.1–11.2 s | 8.8–8.9 s |
| 8 | 1.93–2.05 s | 1.17–1.24 s |

Full Alpha-1 v1.0 review at 8 threads: old 124 s, new 111 s wall.

The extra work per shot is four SplitMix64 steps (a few ns) against microseconds of frame propagation and matching.
The new code came out **faster**. A likely cause is that block-local accumulators replace OpenMP reduction variables
captured by reference in the hot loop's lambda. This was not profiled further. The per-block `Tally` vector costs
40 bytes per 1024 shots.

## How much the numbers move

The new code draws a different (equally valid) sample, so each Monte Carlo number moves by roughly one Monte Carlo
standard error, in a random direction. That is the same size as changing the thread count with the old code. All runs
below are scratch runs (`--out` or a scratch copy of the tree). Nothing was committed.

**Gen-0 ensemble** (100000 shots per round count):

- New code, 30 base seeds (11, 1011, …, 29011): mean 0.012101, standard deviation 0.000135, range
  0.011823–0.012309.
- Old code, 8 thread counts (table above): mean 0.012033, standard deviation 0.000133.
- The two means differ by 0.000068 ± 0.000053 (1.3σ), so there is no bias.
- The empirical seed-to-seed standard deviation is 1.1× the reported `eps_err`, so `eps_err` is a fair noise scale.

**Full Alpha-1 v1.0 review**, new code, 8 threads, against the committed `results/alpha1_review.json`. The column
"old @8" is the old source rebuilt at 8 threads, which gives the same comparison for a thread-count or build change
alone. z = Δ / √(err₁² + err₂²).

| quantity | committed | old @8 | new | new vs committed | new vs old @8 |
|---|---|---|---|---|---|
| Gen-0 eps | 0.012323 | 0.012002 | 0.011958 | −3.0 % (z −2.09) | −0.4 % (z −0.25) |
| + tunable couplers | 0.005377 | 0.005590 | 0.005425 | +0.9 % (z +0.42) | −3.0 % (z −1.43) |
| + fast readout | 0.005379 | 0.005610 | 0.005546 | +3.1 % (z +1.48) | −1.1 % (z −0.57) |
| + active reset (ladder) | 0.004476 | 0.004632 | 0.004524 | +1.1 % (z +0.49) | −2.3 % (z −1.10) |
| **Alpha-1 headline eps** (2× shots) | 0.004568 | 0.004487 | 0.004467 | −2.2 % (z −1.46) | −0.4 % (z −0.28) |
| distance-5 eps | 0.001751 | 0.001724 | 0.001706 | −2.6 % (z −0.56) | −1.1 % (z −0.22) |
| Λ = eps3/eps5 | 2.608 | 2.602 | 2.619 | +0.4 % | +0.6 % |
| logical lifetime (µs) | 186.7 | 190.1 | 191.0 | +2.3 % | +0.4 % |

- **Over the six quantities with standard errors:**
  - new against committed: rms z 1.25, max |z| 2.09 (Gen-0);
  - new against old @8: rms z 0.79, max 1.43;
  - old @8 against committed, a pure thread-count and build change: rms z 1.58, max 2.06.

  **The deterministic RNG moves the numbers no more than a thread-count change already did.**
- **Budget and sensitivity tables** (15 values, no `eps_err` in the JSON): the standard error is approximated as
  `err_ref * sqrt(eps/eps_ref)` from the ladder's Alpha-1 point, which uses the same shots and rounds.
  - New against committed: rms z 1.00, max 2.04.
  - Old @8 against committed: rms 0.84, max 1.83.
  - These rows share seeds, so they are correlated and their mean z is not an independent statistic.
- **The headline conclusions are unchanged.** Alpha-1 rounds to 0.45 % per round (0.447 %, against 0.457 %
  committed). Λ is about 2.6. The error-budget ranking is unchanged (T1, then control, Tphi, readout, reset, ZZ,
  thermal). In the sensitivity table the top two are unchanged (CZ control error, then dielectric Q). Readout
  500 -> 300 ns and T_eff 45 -> 35 mK swap places: they are 1.3–1.6 % apart, which is within one standard error.

**Alpha-1 v1.1 review** (`alpha1_v11_review 100000 --out <scratch>`, 8 threads), new code against the committed
`designs/alpha1-v11/alpha1_v11_review.json` (old code, 8 threads, `bff9d10`):

| quantity | committed | new | change |
|---|---|---|---|
| headline eps, v1.0 model | 0.4401 % | 0.4501 % | +2.3 % (z +1.51) |
| **headline eps, v1.1 model** | **0.3679 %** | **0.3738 %** | +1.6 % (z +1.40) |
| leakage off | 0.3740 % | 0.3676 % | −1.7 % (z −1.08) |
| leakage, no removal | 0.4062 % | 0.4118 % | +1.4 % (z +0.90) |
| leakage-removing reset | 0.3868 % | 0.3767 % | −2.6 % (z −1.67) |
| LRR + DQLR | 0.4136 % | 0.4087 % | −1.2 % (z −0.78) |
| distance-5 eps | 0.1227 % | 0.1129 % | −8.0 % (≈ z −1.7, est.) |
| Λ | 3.12 | 3.32 | +6.5 % (≈ 1.3σ, est.) |
| logical lifetime (µs) | 252.2 | 248.2 | −1.6 % |
| runtime (8 threads) | 549 s | 336 s | −39 % |

- **Rms z over the six points with `eps_err`: 1.27, max |z| 1.67.** That is what an independent resample gives.
- The distance-5 run stores no `eps_err`. The estimate assumes 50000 shots at 15 rounds, so eps_err is about 4e-5
  per point, and Λ carries about 3.5 % relative error per run.
- The 50-round P_L values move by 0.2–1.5 %.
- Leave-one-out budget and sensitivity entries are *relative* changes, differences of two noisy estimates, so their
  small entries swing a lot in relative terms:
  - thermal +0.015 -> −0.006, ZZ +0.015 -> +0.021, flux ×3 +0.110 -> +0.068.
  - The largest absolute change is 0.058 (residual 2e-3).
  - The ranking of the large entries is unchanged (T1, control, Tphi, readout, reset).
- The leaked fractions for the v1.0 model (about 2–3e-5) rest on very few leakage events and move by up to 40 %.
- The CZ calibration block differs only in the 12th digit, because this run read `cz_cache.json`, which stores
  `%.12g`, and the committed run calibrated afresh. That is not Monte Carlo.
- **Conclusion unchanged:** v1.1 at about 0.37 % per round beats v1.0 at about 0.44–0.45 %; leakage-removing reset
  helps; DQLR does not pay at the default rates; Λ is above 3.

**Beta-1 review** (`beta1_review 60000 --out <scratch>`, 8 threads), new code against the committed
`designs/beta1/beta1_review.json`:

| quantity | committed | new | change |
|---|---|---|---|
| eps d=3 (5–20 rounds) | 0.0827 % ± 0.0019 | 0.0817 % ± 0.0019 | −1.2 % (z −0.36) |
| eps d=5 (5–20 rounds) | 0.0155 % ± 0.0012 | 0.0170 % ± 0.0012 | +9.5 % (z +0.88) |
| eps d=5, 9 rounds | 0.0180 % ± 0.0020 | 0.0209 % ± 0.0022 | +16 % (z +0.99) |
| eps d=7, 9 rounds (100000 shots) | 0.00456 % ± 0.00071 | 0.00289 % ± 0.00057 | −37 % (z −1.83) |
| Λ(3→5) | 5.33 ± 0.42 | 4.81 | −10 % (≈ −0.9σ) |
| Λ(5→7) | 3.96 ± 0.76 | 7.24 ± 1.6 | +83 % (≈ +1.8σ) |
| Alpha-1 v1.1 reference Λ | 2.40 | 2.44 | +1.6 % |
| 50-round P_L, d=5 LRR / LRR+DQLR | 4.68 % / 5.73 % | 4.78 % / 5.60 % | +2.0 % / −2.2 % |
| runtime (8 threads) | 872 s | 726 s | −17 % |

- Rms z over the four distance points: 1.15.
- **d = 7 is the statistics-limited point:** 41 logical errors in 100000 shots committed, 26 new. Λ(5→7) therefore
  moves a lot in absolute terms, but within its own error bar. BETA-1.md already quotes it with ± 0.76; any rerun
  should keep quoting the error bar. This is the largest shift of about 40 compared quantities, as expected for
  the most extreme of that many.
- **The d=3 leave-one-out budget** (60000 shots, relative changes against a 120000-shot reference) moves by up to
  0.10 in relative terms (leakage −0.18 -> −0.09, readout −0.11 -> −0.18, Tphi −0.16 -> −0.10). A relative change
  between two runs has a standard error of about 0.057, so these shifts are 1–1.7σ. The dominant terms (T1 −0.74,
  control −0.55 / −0.49) are unchanged in rank. The small ZZ and thermal entries change sign.
- **Conclusions unchanged:** below threshold, with Λ(3→5) about 5, and d = 5 at about 0.016 % per round.

**Summary of the shift estimate:** adopting per-shot seeding moves every Monte Carlo number by a fresh-sample amount.
- Over 16 quantities with a reported standard error across the three reviews, the rms z is 1.23 (mean z −0.11; 1.0 expected for an independent resample, χ² = 24 on 16 degrees of freedom). The largest
  is 2.1, Gen-0 against the committed FMA-build value.
- That is the same size as the old code's own spread across thread counts (rms z 1.58, old 8-thread source build
  against the committed v1.0).
- **No headline conclusion changes.** Numbers that are statistics-limited (d = 7, Λ(5→7), leave-one-out entries near
  zero, leaked fractions around 1e-5) move visibly, inside their error bars.


## Migration impact if adopted

- **Every frame-Monte-Carlo number changes once** (a new sample; see above):
  - `results/alpha1_review.{json,txt}` (= `designs/alpha1/`);
  - `designs/alpha1-v11/alpha1_v11_review.{json,txt}`;
  - `designs/beta1/beta1_review.{json,txt}` and the generated `BETA-1.md`;
  - three printed lines of `results/test_output.txt`;
  - the whitepaper and design-doc tables quoting them.

  The pulse calibrations and `results/cz_cache.json` are not Monte Carlo and do not change. Nor do the GDS files or
  the layout outputs.
- **If this merges without regenerating:** the committed outputs stay reproducible from their recorded source commits
  (already listed in `results/PROVENANCE.md`), but no longer from the branch head at any thread count. Recommendation:
  adopt and regenerate in the same release, then rewrite the PROVENANCE "Thread count matters" bullet as history.
  `openmp_threads` in each JSON's provenance becomes informational.
- **Runtime** for regeneration at 8 threads, measured in scratch with the new code: v1.0 111 s, v1.1 336 s (old 549 s), Beta-1 726 s (old 872 s).
- **API:** none. `run_memory_frame(M, shots, seed)` keeps its signature. Callers that relied on a specific
  per-thread sample (none in the repo) would see a different, thread-independent one.
- **Smoke tests:** unaffected (no golden values).
- **Choice of target:** this branch is stacked on `fix/alpha1-v1.0-repro` because it needs `-ffp-contract=off`.
  Without that flag, an `-march` choice still changes results through FMA contraction, independently of threads. It
  can be retargeted to `release/v0.3` (or `main` for the next release) once that PR is merged, depending on whether
  the regeneration goes into v0.3.
