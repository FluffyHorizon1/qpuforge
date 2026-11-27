# Whitepaper merge: choices for the author

Issue file: `q06-whitepaper-conflict.md`. Two texts exist for `docs/WHITEPAPER.md`:

- **A**, on `main`: the PR #2 note (115 lines, 25 Nov 2026). It adds the coupler-transmon and CZ flux-noise sections and
  the `docs/af17-*.png` figures. It also says `analyze()` solves the tuned-α |11⟩–|20⟩ point and quotes 299–658 MHz;
  `main` does not do that (issue file `q02-whitepaper-tuned-alpha.md`).
- **B**, on `docs/whitepaper-v0.3` (`6671993`, based on `release/v0.3`): the 208-line v0.3 revision. It adds Status and
  known issues, keeps v1.0 and v1.1 apart, states the leakage model, and drops the unreproducible 299–658 MHz claim. It
  was written before PR #2 merged, so it lacks the coupler and flux-noise sections.

Merging B into `main` is an add/add conflict on the same path, so this branch does not touch `docs/WHITEPAPER.md`. It adds
`docs/WHITEPAPER-merge-proposal.md`, which is B with A's abstract and A's two new sections ported in, every figure and
reference checked against `main` 78b9f8b, and [Cn] markers where you have to choose. Nothing in it is a new measurement
except the `main` figures quoted below (local runs, not committed).

## What the proposal already does (mechanical, no choice needed)

- Check count 101 → 186 and 20 → 21 CTest entries (`main` 78b9f8b).
- Provenance and `-ffp-contract=off` described as on `main` (PR #6), not "not yet on release/v0.3".
- `device.hpp:467` line reference replaced by the expression, because the line moved (548 on `main`).
- `analyze designs/alpha1/alpha1.json` on `main`: mean T1 70.5 µs, mean CZ error 2.43e-3, max |idle ZZ| 14.0 kHz,
  0 collisions (the 24 Nov file gave 70.8 µs, 2.40e-3, 5 kHz).
- The PR #2 figures are referred to as on `main` under `docs/`.

## Choices

| # | Question | Options | Proposal's default |
|---|---|---|---|
| C1 | Which text is the base? | (a) B with A's new sections ported in (this proposal); (b) A with B's corrections applied; (c) keep both files under different names | (a) |
| C2 | Alpha-1 v1.0 logical-error table (1.23 / 0.54 / 0.54 / 0.45 %) | (a) keep as the historical 24 Nov review, labelled pre-PR #2 (needs nothing else); (b) regenerate the v1.0 review with the PR #2 model and replace the table (also needs `designs/alpha1/alpha1_review.*`, `results/alpha1_review.*`, PROVENANCE, ALPHA-1.md §3–7) | (a) |
| C3 | Tuned-α |11⟩–|20⟩ point | (a) text says `analyze()` uses the idle-α estimate, 300–660 MHz (matches `main` today, and branch `docs/whitepaper-tuned-alpha-matches-code`); (b) merge `fix/cz-tuned-alpha-solve` and say 298.9–657.5 MHz (changes analyze outputs; regenerate) | (a) |
| C4 | Figures | (a) v0.3 GDS previews from `designs/alpha1-v11/` only, and say the PR #2 PNGs do not match `AF17_v11.gds`; (b) also show the PR #2 PNGs, labelled as an earlier drawing; (c) delete `docs/af17-*.png` from `main` | (a); the PNGs stay on `main` until you say otherwise |
| C5 | CZ flux-noise term scope | (a) describe it as on `main`: added for every design with a tunable higher qubit (moves Gen-0 and fixed-coupling numbers); (b) merge `fix/cz-flux-term-opt-in` and describe it as Alpha-1 only | text states (a) and flags it |
| C6 | Alpha-1 v1.1 / Beta-1 direct coupling | (a) v1.1 and Beta-1 stay at 4 MHz as generated (merge `fix/alpha1-v11-keep-4mhz-coupling`); (b) 1 MHz is intended for them too: regenerate v1.1 and Beta-1 outputs, then update the v1.1 tables | text assumes (a) |
| C7 | Thread-count-independent Monte Carlo | (a) mention as proposed (`feat/deterministic-mc-on-main`); (b) adopt it, regenerate all MC tables, drop the "sample, not exact" caveats | (a) |
| C8 | Author line and status | keep "draft for the author's review" until you approve, then date the revision | keep |
| C9 | Order of sections | Status and known issues first (B) or Abstract first (A) | Abstract, then Status |

Choosing (b) in C2, C3, C5, C6 or C7 means regenerating committed outputs; the related note or branch lists which.
