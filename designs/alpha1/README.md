# Alpha-1 v1.0 files: which model they describe

PR #2 (`design/alpha1`, merged 27 Nov 2026) changed the Alpha-1 model on `main`: direct coupling 4 MHz → 1 MHz
(link capacitance 0.14 → 0.035 fF), a coupler transmon parked at the exchange null (idle ZZ computed, −4 to −14 kHz,
instead of a typed 5 kHz), and a CZ flux-noise term in `e_2q`. Not every file here was regenerated for it.

| File | Describes | Regenerated after PR #2? |
|---|---|---|
| `alpha1.json` (and `examples/alpha1.json`) | the current `alpha1()` preset: 1 MHz coupling, `model_coupler: true` | yes (PR #2) |
| `ALPHA-1.md` summary table and section 2 | the current preset | yes (PR #2) |
| `ALPHA-1.md` sections 3–7 | **Alpha-1 v1.0 review, 24 Nov 2026** (4 MHz, typed 5 kHz, no flux term) | **no: historical** |
| `alpha1_review.txt`, `alpha1_review.json` (= `results/alpha1_review.*`) | **Alpha-1 v1.0 review, 24 Nov 2026** | **no: historical** |
| `alpha1_layout.svg`, `alpha1_layout.png` | v1.0 floor plan | no |

The historical review files are generated output. They are left byte-for-byte as committed so that they stay
reproducible: check out `3fe4d23` and run the binary committed there, `OMP_NUM_THREADS=2 ./examples/alpha1_review`
(full recipe and the reasons, thread count and FMA contraction, in `results/PROVENANCE.md`). Do not edit them by hand.

Building `examples/alpha1_review` from `main` and running it gives a different review, because the Alpha-1 rows use the
PR #2 model. Replacing the historical files with such a run is a separate decision (issue file
`q04-stale-alpha1-v10.md`); it would also need `results/PROVENANCE.md`, `ALPHA-1.md` sections 3–7 and the whitepaper's
v1.0 table updating together.

The Alpha-1 design that has a mask and a v0.3 review is Alpha-1 v1.1, in `designs/alpha1-v11/`.
