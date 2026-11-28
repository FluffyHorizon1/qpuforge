# Contributing to qpuforge

Thank you for helping with qpuforge. This repository is **header-only C++17** (the implementation lives under `include/qpuforge/`, with a small CLI in `src/` and tests under `tests/`).

## Licence

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

qpuforge is proprietary software, licensed under the [Secured Cyber Solutions Proprietary Software Licence v1.0](LICENSE). It is **not open source**: no use, copying, modification or distribution is permitted without written permission from SCS. Third-party components are listed in [NOTICE](NOTICE).

By contributing, you confirm that you have the right to submit the change and that SCS may use it under that licence.

## Branch and pull-request workflow

- Do your work on a **new branch**. Open a pull request for review.
- **Never push to `main` or to any `release/*` branch.** Do not merge, close, tag, or change repository settings unless the owner has approved that specific action.
- Prefer a short, descriptive branch name (`fix/…`, `docs/…`, `chore/…`, `feat/…`, `ci/…`, `test/…`).
- Keep pull requests focused. Describe what changed and how you tested it.

## Build and test

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
./build/test_qpuforge                 # 500 validation checks; or:
ctest --test-dir build -LE smoke      # one CTest entry per group, excluding end-to-end smokes
```

`ctest` labels are `unit`, `physics`, `pulse`, `qec`, `layout` and `smoke`. The four `smoke` entries run the Alpha-1 v1.1, Beta-1 and Gamma-1 review programs and the lattice-surgery study (`surgery_study`) with 200 shots; they write into the build tree and leave committed outputs untouched.

Please run `./build/test_qpuforge` or `ctest --test-dir build -LE smoke` before opening a pull request. Mention any tests you could not run, and why.

## Code style

- C++17, no compiler extensions required of contributors. Keep the library **header-only**.
- Prefer small, focused headers under `include/qpuforge/`. Match the existing indentation, naming, and include style.
- Do not add a build-time dependency without a clear reason in the pull request.
- If a change would alter committed review numbers or caches, say so in the pull request and leave provenance notes accurate (`results/PROVENANCE.md`).

Licensing enquiries: https://securedcybersolutions.co.uk
