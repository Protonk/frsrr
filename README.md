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
- Ability to run the frsr with a custom iteration formula, specified in R formula syntax
- Bin input range and compute optimal magic constants for each bin, if efficiency really isn't your thing.

## FISR or FRSR?

When the FRSR became famous (and only then), it was referred to as an inverse, meaning "[multiplicative inverse](https://en.wikipedia.org/wiki/Multiplicative_inverse)". 

By contrast, in both the original source of Quake's FRSR tucked away in [a math library since 1986](https://www.netlib.org/fdlibm/e_sqrt.c) and one of the first software libraries ever written--due to Alan Turing, D.G. Prinz, and Cecily Poppelwell--1/sqrt(x) is the "[reciproot](https://0x5f37642f.com/documents/ManchesterRecipRoot.pdf)". Mike Day also argues for the name FRSR in his [2023 generalization of the FRSR](https://arxiv.org/abs/2307.15600) to support any rational power or precision of base. Fame has an interia all its own, so we shall see which name prevails.

## Installation

Install the package from GitHub using the `devtools` package:

```R
# install.packages("devtools")
devtools::install_github("Protonk/frsrr")
```

## Usage

```R
library(frsrr)

# Custom parameters
result <- frsr(c(1, 4, 9, 16), magic = 0x5f375a86, NRmax = 2, A = 1.6, B = 0.6)
## result is a vector of length 4
print(result)
# [1] 0.9990148 0.4995074 0.3337626 0.2497537

# Optional detail 
result.df <- frsr(c(pi, 2^-31, 0.4, 6.02e23), detail = TRUE)
## result.df is a dataframe with 4 rows and 7 columns
print(result)
#          input      initial    after_one        final        error          diff iters
# 1 3.141593e+00 5.735160e-01 5.639570e-01 5.639570e-01 0.0004121269 -9.558976e-03     1
# 2 4.656613e-10 4.693787e+04 4.632937e+04 4.632937e+04 0.0002499308 -6.085039e+02     1
# 3 4.000000e-01 1.632430e+00 1.578616e+00 1.578616e+00 0.0015955754 -5.381417e-02     1
# 4 6.020000e+23 1.306493e-12 1.288484e-12 1.288484e-12 0.0002824810 -1.800936e-14     1

## Generate 4 samples using default parameters and random input and magic values
# Optionally, parameters can be returned with keep_params = TRUE
set.seed(123)
samples <- frsr_sample(4, keep_params = TRUE)
# knitr's kable() makes tables look better on github
library(knitr)
kable(samples, format = "simple")
#      input     initial   after_one      final       error        diff   iters        magic   NRmax     A     B   tol
# ----------  ----------  ----------  ---------  ----------  -----------  ------  -----------  ------  ----  ----  ----
#  0.3035076   1.8942177    1.809921   1.809921   0.0028867   -0.0842963       1   1598040167       1   1.5   0.5     0
#  0.8480096   1.1454246    1.080945   1.080945   0.0045856   -0.0644797       1   1597974746       1   1.5   0.5     0
#  0.7123331   1.1783870    1.184784   1.184784   0.0000444    0.0063969       1   1597113118       1   1.5   0.5     0
#  0.9425029   0.9680235    1.024561   1.024561   0.0053300    0.0565371       1   1597011026       1   1.5   0.5     0

## Find optimal constant for 4 bins betweeon 0.25 and 1.0
set.seed(123)
bins <-  frsr_bin(n_bins = 4)
kable(bins, format = "simple")
#  Location   Range_Min   Range_Max        Magic   Avg_Relative_Error   Max_Relative_Error    N
# ---------  ----------  ----------  -----------  -------------------  -------------------  ---
#         1      0.2500      0.4375   1597469260            0.0218475            0.0331088    4
#         2      0.4375      0.6250   1597202340            0.0052858            0.0090786    4
#         3      0.6250      0.8125   1597247742            0.0083959            0.0136274    4
#         4      0.8125      1.0000   1597610966            0.0158552            0.0252852    4
```

## Reproducibility

- All C++ entry points wrap their random draws in `Rcpp::RNGScope`, so `set.seed()` in R fully controls the stochastic components exposed via `.Call()` and the higher level helpers.
- Functions that fan out across threads (`frsr()` and `frsr_bin()`) accept a `threads` argument (and respect `getOption("frsrr.threads")`) that internally calls `RcppParallel::setThreadOptions()`, keeping the worker count explicit and stable across sessions.
- Sampling helpers never touch `unif_rand()` from worker threads; instead, the main R thread prepares any random inputs and the parallel code only performs pure numeric transforms. This keeps the results independent of thread scheduling.
- The test suite pins a seed, runs the stochastic helpers twice, and checks for bitwise identical results so regressions in RNG wiring are caught automatically.

## Our friends the robots

This project was built with the paid assistance of two AI agents, [Perplexity AI](https://www.perplexity.ai/) and [OpenAI's Codex](https://chatgpt.com/codex). Perplexity AI was used until version `0.8.9`. Codex assisted with later development, including a large refactoring guaranteed by testing. 
