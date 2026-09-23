# Fast Reciprocal Square Root R Package, frsrr

## Overview

This R package implements a parameterized Fast Reciprocal Square Root (FRSR) algorithm, also known as Fast Inverse Square Root (FISR), written in C++. You can read more about the FRSR at [0x5f37642f.com](https://0x5f37642f.com/) or see reasons why you might find instrumenting the output fun at [FastInverseSqrt-Visualized on GitHub](https://github.com/hyland-uw/FastInverseSqrt-Visualized).

## Why?

I like R! R doesn't have a type for 32 bit floats, so I wanted a way to mess with the FISR--a 32-bit specific implementation detail (ignoring for now the various extensions to other bases)--in R. Now you can mess with it too!

## Features

- Customizable parameters for fine-tuning accuracy and performance (or the reverse)
- C++ parallel implementation for speed so you can get the wrong answer faster
- Fast sampler to ease sampling over parameter ranges
- Optional detailed output including initial approximation, intermediate steps, and error metrics
- Compare sampled magic constants within input bins or across log2 phases.

## FISR or FRSR?

When the FRSR became famous (and only then), it was referred to as an inverse, meaning "[multiplicative inverse](https://en.wikipedia.org/wiki/Multiplicative_inverse)". 

By contrast, in both the original source of Quake's FRSR tucked away in [a math library since 1986](https://www.netlib.org/fdlibm/e_sqrt.c) and one of the first software libraries ever written--due to Alan Turing, D.G. Prinz, and Cecily Poppelwell--1/sqrt(x) is the "[reciproot](https://0x5f37642f.com/documents/ManchesterRecipRoot.pdf)". Mike Day also argues for the name FRSR in his [2023 generalization of the FRSR](https://arxiv.org/abs/2307.15600) to support any rational power or precision of base. Fame has an interia all its own, so we shall see which name prevails.

## Installation

Building requires C++20 and RcppParallel >= 5.1.11-2; older bundled TBB headers
can fail to compile with C++20. Install from GitHub using `devtools`:

```R
# install.packages("devtools")
devtools::install_github("Protonk/frsrr")
```

## Usage

```R
library(frsrr)

# NRmax = 0 retains the initial bit-hack approximation.
frsr(c(1, 2, 4), magic = 0x5f3759df, NRmax = 0, tol = 0, threads = 1)

# Custom coefficients, with intermediate results and relative error.
frsr(c(1, 4, 9, 16), magic = 0x5f375a86, NRmax = 2,
     A = 1.6, B = 0.6, tol = 0, detail = TRUE, threads = 1)

set.seed(123)
samples <- frsr_sample(4, magic_min = 0x5f3759df, magic_max = 0x5f3759df,
                       NRmax = 1, tol = 0, threads = 1, keep_params = TRUE)
samples

set.seed(123)
bins <- frsr_bin(n_bins = 4, float_samples = 64, magic_samples = 32,
                  NRmax = 1, threads = 1)
bins
attr(bins, "settings")
```

## Numerical contract and reproducibility

Inputs must convert to positive normal IEEE-754 float32 values; original values
must not exceed the largest finite float32, `(2 - 2^-23) * 2^127`.
The input and coefficients round to float32, and each operation of
`y * (A - ((B * x) * y) * y)` rounds separately without fused multiply-subtract.
The reference and error measurements use double precision and target the
float32-rounded input, rather than the original R double. Ordinary rounding
with gradual underflow is assumed; fast-math and altered rounding modes are
unsupported. This package instruments an arithmetic experiment, not a hardware
throughput benchmark or a historical-machine emulator.

`detail = TRUE` returns `input`, `initial`, `after_one`, `final`, `error`, `diff`
and `iters`. Nonfinite approximations from exploratory parameters have infinite
error. The separate mantissa-only `enre` experiment and the custom-formula
`frsr_NR()` API have been removed; configurable A/B coefficients remain supported.

Sampling uses half-open intervals `[x_min, x_max)`. Log-stratified sampling
selects representable normal float32 values directly; equal magic bounds select
that constant, and reversed magic bounds remain supported.

Use `set.seed()` for sampling and `threads` (or `options(frsrr.threads)`) to
control parallel execution. Random inputs are generated on the main R thread.
Bin candidates share samples, and phase candidates share a grid independent of
candidate order. Bin aggregation uses double precision and fixed reduction order,
so the same samples and build give identical measurements across thread counts.
This does not promise bitwise equality across platforms or toolchains.

Search results identify the best tested candidate on the sampled inputs. Exact
bin-objective ties select the smallest magic integer; phase ties compare J,
roughness R, then the smallest magic integer. Candidates with nonfinite
approximations are excluded. Bin results retain their configuration in a
`settings` attribute (preserved by `saveRDS()`, not CSV); phase results include a
`settings` component. See the function help for the recorded fields.

## Our friends the robots

This project was built with the paid assistance of two AI agents, [Perplexity AI](https://www.perplexity.ai/) and [OpenAI's Codex](https://chatgpt.com/codex). Perplexity AI was used until version `0.8.9`. Codex assisted with later development, including a large refactoring guaranteed by testing. 

## License

This software is released under a variant of the MIT license. Copyright claims have been made on various versions of the FRSR, including [an attempt to sue Microsoft over Copilot's regurgitation of the code](https://www.saverilawfirm.com/our-cases/github-copilot-intellectual-property-litigation?utm_source=chatgpt.com). All are likely to fail or have failed. License text follows:

Permission to use, copy, modify, distribute, and sell this software and its
documentation for any purpose is hereby granted without fee, provided that
the above copyright notice appear in all copies and that both that
copyright notice and this permission notice appear in supporting
documentation.  No representations are made about the suitability of this
software for any purpose.  It is provided "as is" without express or
implied warranty.