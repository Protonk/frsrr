This package is an R/C++ lab for diagnosing Fast Reciprocal Square Root (FRSR) behaviour. Each agent name maps directly to a function, shell command, and success criteria so you can grab the right workflow quickly.

## Workflows

| Agent | Use it for | Entry point |
| --- | --- | --- |
| ComputeAgent | Fast, vectorised reciprocal-square-root diagnostics or Quake-style bit-hack replays with new parameters | `frsr()` (`R/frsr.R`, `src/frsr.cpp`) |
| SamplerAgent | Reproducible corpora of `{input, magic}` pairs prior to error-envelope studies or parameter sweeps | `frsr_sample()` (`R/sample.R`, `src/sample.cpp`) |
| BinSearchAgent | Segmenting `[x_min, x_max]` into bins to locate the best restoring constant per bin under a named metric | `frsr_bin()` (`R/bin.R`, `src/search.cpp`) |
| PhaseAgent | Inspecting phase/exponent bias and picking the most stable magic constant over a log2 grid | `frsr_phase()` (`R/phase.R`, `src/phase.cpp`) |
| NRFormulaAgent | Experimenting with Newton-style iteration formulas written in R syntax | `frsr_NR()` (`R/nr.R`) |
| PackageQAAgent | Keeping the R API, documentation, and C++ bridge in sync | `testthat`, `devtools`, `R CMD check` |

## Invocation cheatsheet

- ComputeAgent: `R -q -e "pkgload::load_all('.'); print(frsrr::frsr(c(1,4,9,16), detail = TRUE, keep_params = TRUE))"`
- SamplerAgent: `R -q -e "pkgload::load_all('.'); set.seed(123); print(frsrr::frsr_sample(64, keep_params = TRUE))"`
- BinSearchAgent: `R -q -e "pkgload::load_all('.'); print(frsrr::frsr_bin(n_bins = 4, float_samples = 1024, magic_samples = 2048, objective = 'max_relative_error', dependent = 'avg_relative_error'))"`
- PhaseAgent: `R -q -e "pkgload::load_all('.'); str(frsrr::frsr_phase(phases = 32L, exponents = -8L:8L, per_cell = 16L, magics = as.integer(seq(0x5f3750df, 0x5f3765df, by = 256L)), q = 0.95, NRmax = 0L))"`
- NRFormulaAgent: `R -q -e "pkgload::load_all('.'); custom <- quote(y * (1.6 - 0.6 * x * y^2)); print(frsrr::frsr_NR(c(1,4,9), formula = custom, NRmax = 4, tol = 1e-6))"`
- PackageQAAgent: VS Code task `R: Test` in `.vscode/tasks.json`; `R: Check` runs CRAN-style `R CMD check` when requested.

## Input/output contracts

- `frsr()` accepts numeric `x > 0` plus optional per-element `magic`, `NRmax`, `A`, `B`, `tol`, and `threads`. With `detail = FALSE` it returns a numeric vector of `length(x)`; `detail = TRUE` exposes all diagnostics and (when `keep_params = TRUE`) the parameter columns.
- `frsr_sample()` takes `n`, float bounds (`x_min`, `x_max`), magic bounds (`magic_min`, `magic_max`), a sampler `method`, and Newton arguments via `...`. Output mirrors `frsr(detail = TRUE)` and must have exactly `n` rows.
- `frsr_bin()` consumes `x_min`, `x_max`, `n_bins`, `NRmax`, metric names, `float_samples`, `magic_samples`, optional sampler args, and returns `{N_bins, Location, Range_Min, Range_Max, Magic, Objective, Dependent}`.
- `frsr_phase()` takes `{phases, exponents, per_cell, magics, q, NRmax}` and returns `{magic, J, R, phase_tbl, heat}`. `phase_tbl$n` equals `length(exponents) * per_cell`, and `heat` has `length(exponents)` rows by `phases` columns.
- `frsr_NR()` expects numeric `x`, optional `magic`, a quoted `formula` using `y` and `x`, `NRmax >= 1`, `tol >= 0`. It returns columns `input`, `initial`, `final`, `error`, `converged`, `conv_rate`, `iters`.
- Package QA flows read `DESCRIPTION`, `NAMESPACE`, `R/`, `src/`, `man/`, `tests/` and can emit `frsrr_<version>.tar.gz` plus `frsrr.Rcheck/`.

## Tooling facts

- Reproducibility: every `.Call` helper (`src/frsr.cpp`, `src/sample.cpp`, `src/search.cpp`, `src/phase.cpp`) is wrapped in `Rcpp::RNGScope`. A single `set.seed()` governs deterministic runs—record it whenever you share diagnostics.
- Threads: `frsrr_configure_threads()` in `R/threads.R` honours `getOption("frsrr.threads")`, falls back to `RcppParallel::defaultNumThreads()`, and propagates the final count to the C++ workers. `frsr()` and `frsr_bin()` expose a `threads` argument.
- Toolchain: `SystemRequirements: C++20`. OpenMP flags live in `src/Makevars`; optional Intel TBB include/lib paths are respected. LAPACK/BLAS must be present for CRAN builds.
- Dependencies: Imports `Rcpp`, `RcppParallel`. Suggested helper packages: `testthat (>= 3.0.0)`, `pkgload`, `devtools`, `knitr`.
- Layout reminder: `R/` exposes the R API, `src/` hosts numerics, `tests/testthat/` holds behavioural specs, and `man/` is roxygen output—never edit `.Rd` files by hand.

## Policies

1. Document changes inside `.R` files via roxygen comments and regenerate with `devtools::document(roclets = c('rd','collate','namespace'))`. Never edit `man/*.Rd`.
2. Prefer diagnostic-friendly workloads: modest `n`, `n_bins`, `phases`, or `per_cell` surface bugs faster and keep shared machines responsive.
3. Always note `{seed, Newton parameters, threads}` when sharing measurements so runs stay reproducible.
4. Treat `threads` as advisory; omit or cap it instead of oversubscribing shared hardware.
5. Keep temporary build artefacts (`frsrr_<version>.tar.gz`, `frsrr.Rcheck/`) out of version control and delete them when done.
6. Validate bounds (`x_min > 0`, `x_max > x_min`, `exponents` within `[-126, 127]`, etc.) in R before calling `.Call` so errors remain informative.

## Health checks

- Interpreter: `R --version` already succeeds on this machine.
- Tests: `R -q -e "pkgload::load_all('.'); testthat::test_dir('tests/testthat', reporter = 'summary')"` hits every exported workflow plus the low-level samplers from the README.
- Release: `R CMD build .` then `R CMD check frsrr_<version>.tar.gz` mirrors CRAN. Inspect `frsrr.Rcheck/tests/*.Rout`. Only run when asked.
- Smoke demos: README snippets (custom Newton parameters, sampling tables, bin searches) are handy sanity checks.

## Done criteria

- `frsr()` returns a finite numeric vector (or matching diagnostic frame) whose errors respect the requested tolerance.
- `frsr_sample()` emits exactly `n` diagnostic rows; `keep_params = TRUE` adds the documented parameter columns. Equal seeds reproduce equal tables.
- `frsr_bin()` returns `n_bins` rows with magics inside the supplied bounds and finite `Objective` / `Dependent` metrics.
- `frsr_phase()` guarantees `phase_tbl$n = length(exponents) * per_cell` and the `heat` matrix uses the promised dimensions.
- `frsr_NR()` honours `NRmax`, reports a meaningful `conv_rate`, and marks `converged = TRUE` once each element drops below `tol`.
- Package QA runs: `testthat` passes, and (when invoked) `R CMD build/check` complete cleanly.

## Documentation contract

“Documenting” a function or object always means: add or update roxygen comments inside the corresponding `.R` file. Never touch the generated `man/*.Rd` files.

## Agent running notes

Trim this section first if the file grows too large.

### Discovery notes

- _(none yet)_

### Unknowns

- Availability/configuration of TBB or other threading backends referenced in `src/Makevars`.
- Preferred limits for long-running sweeps on shared hardware.
- Whether `devtools` is pre-installed in every Codex environment.
