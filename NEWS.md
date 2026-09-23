# frsrr 1.1.0

Version 1.1.0 gives evaluation and magic-constant searches a common numerical contract and updates the supported API.

- **Removed `frsr_NR()`.** Custom Newton coefficients remain available through `frsr(A = ..., B = ...)`; use `frsr(..., NRmax = 0)` for the initial bit-hack approximation.
- **Removed `enre`** from detailed `frsr()` and `frsr_sample()` results. The `error` column reports the approximation's absolute relative error.
- **Unified arithmetic and error measurement.** Inputs, coefficients and each refinement operation round to float32. References and errors use double precision and target the float32-rounded input. Numerical results can therefore differ from 1.0.0.
- **Made input and convergence rules explicit.** Inputs must convert to positive normal float32 values, and original inputs must not exceed the largest finite float32. Iteration limits must be non-negative integers; tolerances must be finite and non-negative. Nonfinite exploratory approximations have infinite error and are excluded from candidate searches.
- **Corrected sampling boundaries.** Sampling uses `[x_min, x_max)`, including intervals containing just one admissible float32 value. Equal magic bounds select that constant, and reversed magic bounds remain supported.
- **Made candidate comparisons consistent.** `frsr_phase()` evaluates every candidate on the same sampled grid. Exact ties compare J, roughness R, then the smallest magic integer. `frsr_bin()` selects the smallest magic integer on an exact objective tie and uses a fixed accumulation order for repeatable measurements across thread counts with the same inputs and build.
- **Retained search configuration.** `frsr_bin()` results have a `settings` attribute recording metrics, refinement count, sampler, sample/candidate counts, bounds and threads. `frsr_phase()` results include a `settings` component. Search winners are the best tested candidates on the sampled inputs.
- **Requires RcppParallel >= 5.1.11-2**, with the existing C++20 requirement.
