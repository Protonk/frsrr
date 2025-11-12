

This package is an R/C++ lab for diagnosing Fast Reciprocal Square Root (FRSR) behavior

## Catalog

Agents map to the main workflows exposed in `R/`, `src/`, and `tests/testthat/`.

- ComputeAgent → Use `frsr()` when you need fast, vectorized reciprocal-square-root diagnostics or to replay the classic Quake-style bit hack under new parameters.
- SamplerAgent → Call `frsr_sample()` to build reproducible corpora of inputs/magic constants before evaluating error envelopes or downstream tuning.
- BinSearchAgent → Run `frsr_bin()` when you want to segment an input interval and locate the best restoring constant per bin using explicit objective metrics.
- PhaseAgent → Use `frsr_phase()` to inspect phase/exponent bias and pick the most stable magic constant over a log2 grid.
- NRFormulaAgent → Reach for `frsr_NR()` when experimenting with custom Newton-style iteration formulas expressed in R syntax.
- PackageQAAgent → Exercise `devtools`/`testthat` and `R CMD check` to keep `DESCRIPTION`, documentation, and the C++ bridge coherent.

## Invocation

- ComputeAgent: `R -q -e "pkgload::load_all('.'); print(frsrr::frsr(c(1,4,9,16), detail = TRUE, keep_params = TRUE))"` (repo root).
- SamplerAgent: `R -q -e "pkgload::load_all('.'); set.seed(123); print(frsrr::frsr_sample(64, keep_params = TRUE))"`.
- BinSearchAgent: `R -q -e "pkgload::load_all('.'); print(frsrr::frsr_bin(n_bins = 4, float_samples = 1024, magic_samples = 2048, objective = 'max_relative_error', dependent = 'avg_relative_error'))"`.
- PhaseAgent: `R -q -e "pkgload::load_all('.'); str(frsrr::frsr_phase(phases = 32L, exponents = -8L:8L, per_cell = 16L, magics = as.integer(seq(0x5f3750df, 0x5f3765df, by = 256L)), q = 0.95, NRmax = 0L))"`.
- NRFormulaAgent: `R -q -e "pkgload::load_all('.'); custom <- quote(y * (1.6 - 0.6 * x * y^2)); print(frsrr::frsr_NR(c(1,4,9), formula = custom, NRmax = 4, tol = 1e-6))"`.
- PackageQAAgent:
  - Tests only: `R -q -e "pkgload::load_all('.'); testthat::test_dir('tests/testthat', reporter = 'summary')"` (covers every helper listed in README/man files).
  - Release pass: `R CMD build . && R CMD check frsrr_<version>.tar.gz` (creates `frsrr_<version>.tar.gz` and `frsrr.Rcheck/`).

## Tools

- ComputeAgent: `frsrr::frsr(x, magic, NRmax, A, B, tol, detail, keep_params, threads)` lives in `R/frsr.R` and offloads to `src/frsr.cpp`; optionally parallelized with `RcppParallel`.
- SamplerAgent: `frsr::frsr_sample()` (`R/sample.R`) draws floats via `.Call('_frsrr_sample_inputs', ...)` (`src/sample.cpp`) before piping through `frsr(detail = TRUE, ...)`; `method` switches between midpoint grids, irrational rotations, uniform draws, or log-stratified sampling.
- BinSearchAgent: `frsrr::frsr_bin()` (`R/bin.R`) constructs evenly spaced bins, samples floats/magic ranges, and calls `.Call('_frsrr_search_optimal_constant', ...)` to evaluate metrics.
- PhaseAgent: `frsrr::frsr_phase()` (`R/phase.R`) wraps `.Call('_frsrr_phase_orchestrator', ...)` from `src/phase.cpp`, returning `{magic, J, R, phase_tbl, heat}` summaries.
- NRFormulaAgent: `frsrr::frsr_NR()` (`R/nr.R`) runs the bit-hack start (via `frsr(..., NRmax = 0)`) and iterates an R-side formula; handy when C++ constraints are too rigid.
- PackageQAAgent: `pkgload::load_all`, `devtools::test`, `testthat::test_dir`, `R CMD build`, and `R CMD check`; these touch `DESCRIPTION`, `NAMESPACE`, `man/`, `src/`, and `tests/`.

## I/O

- ComputeAgent: inputs numeric `x > 0` plus optional per-element `magic`, `NRmax`, `A`, `B`, `tol`, `threads`; outputs either numeric vectors (`detail = FALSE`) or data frames with columns `input`, `initial`, `after_one`, `final`, `error`, `enre`, `diff`, `iters`, and optional parameter columns if `keep_params = TRUE`.
- SamplerAgent: inputs `n`, float/magic bounds (`x_min`, `x_max`, `magic_min`, `magic_max`), `method`, and Newton args via `...`; outputs the same diagnostic schema as `frsr(detail = TRUE)` with reproducible samples.
- BinSearchAgent: inputs scalar bounds, counts (`n_bins`, `float_samples`, `magic_samples`), metric names drawn from `{max_relative_error, avg_relative_error, rmse_relative_error}`, and optional `NRmax`, `threads`; outputs a data frame with `N_bins`, `Location`, `Range_Min`, `Range_Max`, `Magic`, `Objective`, `Dependent`.
- PhaseAgent: inputs integers (`phases`, `exponents`, `per_cell`, `NRmax`), candidate `magics`, and `q ∈ (0, 1]`; outputs a list containing the winning magic, scalar metrics `J` (max phase quantile) and `R` (phase roughness), a `phase_tbl` data frame with per-phase stats, and a `heat` matrix (rows = exponents, cols = phases).
- NRFormulaAgent: inputs numeric `x`, optional `magic`, quoted `formula` using `y` and `x`, `NRmax ≥ 1`, `tol ≥ 0`; outputs data frames with `input`, `initial`, `final`, `error`, `converged`, `conv_rate`, `iters`.
- PackageQAAgent: consumes repo metadata (`DESCRIPTION`, `NAMESPACE`, `R/`, `src/`, `man/`, `tests/`) and produces console logs, `frsrr_<version>.tar.gz`, `frsrr.Rcheck/`, and the standard `check.log`.

## Runtime config

- Reproducibility: README’s guarantees come from `Rcpp::RNGScope` in every `.Call` (`src/frsr.cpp`, `src/sample.cpp`, `src/search.cpp`, `src/phase.cpp`), so `set.seed()` in R governs all randomness; document the seed when sharing diagnostics.
- Threads: `frsrr_configure_threads()` (`R/threads.R`) respects `getOption("frsrr.threads")` and falls back to `RcppParallel::defaultNumThreads()`; the heavy helpers (`frsr`, `frsr_bin`) accept a `threads` argument.
- Toolchain: `SystemRequirements: C++20` (DESCRIPTION), OpenMP flags in `src/Makevars`, and optional Intel TBB include/lib paths; CRAN-style builds expect LAPACK/BLAS too.
- Dependencies: Imports `Rcpp`, `RcppParallel`; Suggests `testthat (>= 3.0.0)`; workflows in README assume `devtools`, `pkgload`, and optionally `knitr` for table output.
- Layout reminders: `R/` holds the public wrappers, `src/` the numerics, `tests/testthat/` the behavioral specs, and `man/` mirrors the roxygen docs.

## Policies 

- Documentation: never hand-edit `.Rd` files in `man/`; update the roxygen blocks and ask the user to regenerate with `devtools::document(roclets = c('rd', 'collate', 'namespace'))`. `devtools` may not be loaded in the Codex environment.
- Keep the diagnostic focus: prefer sample sizes and bin grids that expose error structure before scaling to enormous workloads.
- Preserve reproducibility by noting seeds, Newton parameters, and thread counts alongside any reported metrics.
- Treat `threads` as advisory—omit or limit the value when running on shared machines to avoid oversubscription.
- Temporary build artifacts (`frsrr_<version>.tar.gz`, `frsrr.Rcheck/`) belong in `.gitignore` or should be cleaned once inspected.
- When extending formulas or samplers, validate bounds (`x_min > 0`, `x_max > x_min`, exponents within [-126, 127]) in R before invoking `.Call` to keep errors informative.

## Tests/health

- Interpreter check: `R --version` (already verified on this machine) ensures the runtime is present.
- Unit/regression tests: `R -q -e "pkgload::load_all('.'); testthat::test_dir('tests/testthat', reporter = 'summary')"` exercises `frsr`, `frsr_sample`, `frsr_bin`, `frsr_phase`, `frsr_NR`, and the low-level samplers defined in README examples.
- Release validation: `R CMD build .` followed by `R CMD check frsrr_<version>.tar.gz` confirms documentation, compiled code, and tests, mirroring CRAN expectations; inspect `frsrr.Rcheck/tests/*.Rout` for failures.
- Optional smoke demos: snippets from README (custom parameters, sampling tables, bin search) double as quick sanity checks for manual runs.

## Done criteria

- ComputeAgent: returns either a numeric vector matching `1 / sqrt(x)` within the requested tolerance or a diagnostic data frame whose `error`/`enre` columns stay finite for all inputs.
- SamplerAgent: yields exactly `n` diagnostic rows with parameter columns present whenever `keep_params = TRUE`, and replicated seeds reproduce identical tables.
- BinSearchAgent: produces `n_bins` rows with `Magic` inside the provided bounds and finite `Objective`/`Dependent` metrics for each bin.
- PhaseAgent: outputs a list where `phase_tbl$n` equals `length(exponents) * per_cell` for all phases and the `heat` matrix dimensions match `length(exponents) × phases`.
- NRFormulaAgent: respects `NRmax`, reports meaningful `conv_rate`, and marks `converged = TRUE` whenever a positive tolerance is satisfied before iterations expire.
- PackageQAAgent: records a clean `testthat` run and, when executed

## Documentation

Documentation is generated in R projects by a build process. When the agent is asked to document a function or object, the user is always asking for structured documentation to be added to .R files. The agent is not responsible for the files in `/man`.

## Agent running notes

The below are generated by audits or passes where Agents make changes to this file. It is meant to serve as a human and machine readable ephemeral hints to intent or ignorance. If AGENTS.md grows too large in size, trim the below first. 

### Discovery notes

- README.md supplied the package overview, reproducibility guarantees, and usage snippets for `frsr`, `frsr_sample`, and `frsr_bin`.
- DESCRIPTION and NAMESPACE clarified dependencies, exported symbols, and `SystemRequirements`.
- R sources (`R/frsr.R`, `R/bin.R`, `R/sample.R`, `R/nr.R`, `R/phase.R`, `R/threads.R`, `R/RcppExports.R`) described arguments, return shapes, and thread handling.
- C++ sources (`src/frsr.cpp`, `src/sample.cpp`, `src/search.cpp`, `src/phase.cpp`, plus `frsr.h`, `search.h`, `sample.h`, `Makevars`) showed the parallel/RNG behavior backing each tool.
- tests/testthat/ supplied behavioral expectations covering every exported helper and the underlying samplers.

### Unknowns

- Availability/configuration of TBB or alternative threading backends referenced in `src/Makevars`.
- Preferred limits for long-running sweeps (large `n`, `n_bins`, `phases`, or `per_cell`) when running on shared hardware.
