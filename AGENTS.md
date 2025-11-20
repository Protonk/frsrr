# Project guide for contributors

This repo contains a small didactic R/C++ package for exploring fast reciprocal square root routines (Quake-style magic constants, Newton refinements, bin searches, and related diagnostics). Use this file as your single source of truth when touching any part of the project.

## What matters most
- Keep the package reproducible and well-documented; every behavioural change should be explained via roxygen comments in the relevant `.R` file.
- Preserve the teaching value: prefer clear names, short functions, and comments that connect R interfaces to their C++ helpers.
- Avoid surprise state: seed the RNG in examples/tests and note parameters (`magic`, `NRmax`, `tol`, `threads`).

## Repo layout
- `R/`: user-facing R functions and roxygen docs (edit here, not in `man/`).
- `src/`: C++ implementations called via `.Call`; each entry point wraps work in `Rcpp::RNGScope` and may use `RcppParallel`.
- `man/`: generated `.Rd` files. Never edit by hand.
- `tests/`: `testthat` coverage for exported functions and bridge behaviour.
- `README.md`: runnable examples and orientation.

## Build & test
- Fast feedback: `R -q -e "pkgload::load_all('.'); testthat::test_dir('tests/testthat', reporter='summary')"`.
- Full check (only when needed): `R CMD build .` then `R CMD check frsrr_<version>.tar.gz`.
- Keep temporary build artefacts (`frsrr.Rcheck/`, tarballs) out of git.

## Coding expectations
- R code: validate inputs early, keep vectorised operations where practical, and align argument names/behaviour across functions (`frsr`, `frsr_sample`, `frsr_bin`, `frsr_phase`, `frsr_NR`).
- C++ code: prefer modern C++20, keep functions small, and leave brief comments near tricky bit manipulations or Newton steps. Avoid try/catch around imports.
- Threading: respect user-provided `threads` or `options(frsrr.threads)`; avoid oversubscribing shared hardware.
- Performance notes belong in comments or docstrings, not hidden magic numbers.

## Documentation rules
- Update roxygen comments alongside code changes; regenerate docs with `devtools::document(roclets = c('rd','collate','namespace'))` when necessary.
- Examples should run quickly; use small inputs but include seeds so outputs are reproducible.

## Contribution checklist
1. Understand the data contract for the function you are editing (inputs, outputs, and any diagnostic columns).
2. Add or adjust tests that cover new behaviour or bug fixes.
3. Run the relevant tests/commands above and note them in your summary.
4. Keep commits focused and messages descriptive.

## Contact points inside the code
- Core entry points: `R/frsr.R`, `R/sample.R`, `R/bin.R`, `R/phase.R`, `R/nr.R` with matching C++ in `src/`.
- Thread config helper: `R/threads.R`.

Follow this guide to keep the project coherent and approachable for anyone exploring fast reciprocal square root techniques.
