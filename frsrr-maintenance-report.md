# Maintenance completion report

Completed against `471dd4d` (the commit reviewed in the handoff). The original working tree had no tracked changes; the untracked handoff was preserved. No commits or remote changes were made during that maintenance pass.

## Reproductions and local evidence

| Concern | Baseline finding |
| --- | --- |
| Equal magic bounds | Reproduced in both exports. With seed 42 and both bounds `0x5f3759df`, `frsr_sample()` returned unrelated constants and `frsr_bin()` selected `1443151060`. |
| Log-stratified endpoint | `[1, 1 + 2^-24)` incorrectly failed as empty. The fix also handles adjacent double bounds at large exponents whose logarithms coincide. |
| Arithmetic/error disagreement | Confirmed the three source-level contracts. For `x = pi`, one-step evaluation and bin search reported about `0.0004121806` and `0.0004121269`, respectively. |
| `enre` | With zero refinements, inputs 1 and 2 had different ordinary errors but the same `enre`, confirming a separate experiment that discarded exponent parity. |
| Phase sampling | Reversing two candidates with seed 2024 retained the winning magic but changed its J from about `0.03652475` to `0.03657245`. |
| Reduction reproducibility | On 131,072 common inputs, one versus two threads changed average error from about `0.0009372473` to `0.0009372684`. |
| Worker access | No crash was reproduced. Installed Rcpp numeric indexing uses a cached pointer for valid indices, but its bounds-error path calls Rcpp warning machinery. Replaced worker-held Rcpp vectors with supported accessors rather than claiming ordinary indexing always crashed. |
| Build | Unset TBB paths produced a bare `-I` that swallowed the next include flag; unused Fortran linkage then failed. RcppParallel 5.1.9 and 5.1.10 also failed to compile their bundled TBB headers under this C++20/libc++ toolchain. |

The original test suite passed after temporary build workarounds despite the numerical defects. Locale warnings in that baseline were environmental. No targeted numerical defect was already fixed in this checkout. The handoff's worker concern remains a source-level contract issue, not evidence of a reproduced crash.

## Implemented contracts

- Removed `frsr_NR()` and its implementation, export, generated help, and dedicated tests. `frsr(..., NRmax = 0)` still returns the seed. Removed the `enre` output and its entire second calculation; updated schemas, README, and contributor inventory.
- Inputs and coefficients round to IEEE-754 float32. The shared kernel rounds each operation in `y * (A - ((B * x) * y) * y)` separately, using explicit stores to prevent contraction across operations. `std::bit_cast` replaces union punning, with compile-time representation checks. This favors explicit arithmetic over peak throughput.
- Reference values and errors use double precision and target the float32-rounded input. Original input/coefficient values remain in diagnostic columns. `diff` is a double subtraction. Positive tolerance is checked after each step; zero steps produce NA intermediate/difference diagnostics.
- Rounded inputs must be positive normal floats; original inputs may not exceed the largest finite float32. Nonfinite inputs, underflow to subnormal/zero, invalid iteration counts/tolerances, and coefficient overflow fail before worker execution. Finite unusual coefficients and signed magic integers remain available for exploration. Nonfinite approximations produce infinite evaluation error and disqualify search candidates; an all-invalid candidate set errors.
- Samplers use `[x_min, x_max)`. Log-stratified sampling calculates integer float indices directly from the original bounds, with no logarithm round trip. Magic sampling uses offsets, including equal and reversed bounds. Double-valued samplers prevent rescaling from including the upper endpoint.
- Phase candidates share one sampled grid. Rounding at the top of an exponent slab is clamped to its largest float; phase labels still describe the original log2 draws. Exact ties compare J, roughness, then the smallest signed magic integer. Bin objective ties select the smallest signed magic integer.
- Bin measurements accumulate in double precision over fixed sample blocks, joined in fixed order. Parallel tasks cover both blocks and candidates. Identical samples, build and floating-point environment yield identical results across thread counts; no historical-machine or cross-platform bitwise guarantee is made.
- Bin results retain a `settings` attribute identifying metrics, NRmax, sampler, sample/candidate counts, bounds and threads. Phase results retain a `settings` component. These identify the experiment and its best tested candidate, not a global optimum; CSV does not preserve the bin attribute.
- Replaced manual BLAS/LAPACK/Fortran/OpenMP/TBB flags with the [supported RcppParallel integration](https://rcppcore.github.io/RcppParallel/#r-packages). Declared RcppParallel >= 5.1.11-2 after verifying that release builds cleanly here. Prepared one macOS R-package-check workflow, based on the [r-lib example](https://github.com/r-lib/actions/blob/v2/examples/check-standard.yaml). The workflow is retained locally and excluded from the 1.1.0 release, as requested after GitHub rejected workflow updates with the current token permissions.
- Updated roxygen comments and quick seeded examples, enabled Markdown rendering, and regenerated all affected help/NAMESPACE files. The previously ignored `.Rbuildignore` is now visible to git so clean checkouts exclude task reports and CI configuration from package archives.

## Validation actually run

Environment: macOS Sonoma 14.8.3, Apple Silicon, R 4.4.1, Apple clang 15.0.0, MacOSX14.4 SDK, Rcpp 1.0.13.1, testthat 3.2.1.1, roxygen2 7.3.2. RcppParallel 5.1.11-2 and 6.2.1 were installed in separate temporary libraries; the existing system installation (5.1.9) was not changed.

1. Ran the prescribed `pkgload::load_all()` plus `testthat::test_dir()` baseline. Its first attempts failed on the sandbox's processx restriction, malformed include flags, old TBB headers and missing unused Fortran linkage. Reproductions used external TBB paths and a temporary Makevars enabling the old libc++ type traits and clearing FLIBS; package numerical sources were still unchanged.
2. Ran `devtools::document(roclets = c('rd', 'collate', 'namespace'))` and the expanded development suite. Passed with both the current dependency and the verified minimum dependency. Also ran the suite against a clean optimized install with 5.1.11-2.
3. Computed seed/refinement fixture bits independently using exact rational multiplication/subtraction with nearest-even binary32 rounding. Tests cover adjacent values around 1 and 2, ordinary values, normal extremes, rounded-input references, custom coefficient rounding, finite-error filtering, powers-of-four scaling on an interior domain, sampling boundaries, candidate permutations/insertions/ties, metadata persistence, and one/two-thread execution over 131,072 inputs.
4. Built and installed a clean source archive, then ran `R CMD check --no-manual --as-cran`. Final result: **0 errors, 0 warnings, 3 notes; 282 passing assertions, 0 failures, 0 warnings and 0 skips in tests**. All package examples passed. Extracted README usage snippets also ran successfully against the installed package. `git diff --check` passed.

Clean-build commands, run from `/private/tmp/frsrr-maintenance` with `R_MAKEVARS_USER=/dev/null`, `R_LIBS=/private/tmp/frsrr-maintenance/library`, `LANG=en_US.UTF-8`, `LC_CTYPE=en_US.UTF-8` and LC_ALL unset:

```sh
R CMD build /Users/achyland/Desktop/Math/etak/frsrr
R CMD INSTALL --preclean --library=/private/tmp/frsrr-maintenance/library frsrr_1.0.0.tar.gz
R CMD check --no-manual --as-cran frsrr_1.0.0.tar.gz
```

The check used `RCPP_PARALLEL_NUM_THREADS=2`. Build products and raw logs are outside the repository in `/private/tmp/frsrr-maintenance`; `check-final.log`, `frsrr.Rcheck/00check.log`, `frsrr.Rcheck/tests/testthat.Rout`, `install-final.log`, and `baseline.log` contain the detailed evidence. The first full check also reported redundant versioned-LinkingTo metadata; that introduced note was fixed, retaining version enforcement in Imports.

The three remaining notes are:

- CRAN incoming feasibility: this is a new submission, and two existing reference URLs return Semantic Scholar redirect/202 and Stack Overflow 403 responses. These are retained references, not introduced runtime failures.
- Future timestamps: R could not contact its current-time verification service.
- GNU make: the package explicitly declares the requirement used by the supported linker integration.

## Limits and follow-ups

- The system RcppParallel 5.1.9 is below the new minimum. After this maintenance pass, the authorized local update installed RcppParallel 6.2.1 and rebuilt frsrr in the personal R library. A fresh R session selected both personal-library copies and passed all 282 assertions; the system copies and global R configuration were preserved.
- PDF manual generation was omitted because `pdflatex` is absent. Rd parsing, cross-references, documentation/code agreement and executable examples passed.
- Linux, Windows, other architectures/toolchains, and the TinyThread backend were not validated. The optional CI workflow remains local and has not run remotely. Fast-math, altered rounding modes and flush-to-zero environments are unsupported. No performance benchmark was conducted.
- The deferred full worked example should explain original versus rounded-input error, exponent parity and powers-of-four scaling, common-sample candidate comparisons, and saving search settings. Routine interface documentation is complete. Separately review the two reference URLs flagged by CRAN's remote checks.
