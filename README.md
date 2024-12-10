# Fast Reciprocal Square Root R Package, frsrr

## Overview

This R package implements a parameterized Fast Reciprocal Square Root (FRSR) algorithm, also known as Fast Inverse Square Root (FISR), written in C++. You can read more about the FRSR at [0x5f37642f.com](https://0x5f37642f.com/) or see reasons why you might find instrumenting the output fun at [FastInverseSqrt-Visualized on GitHub](https://github.com/hyland-uw/FastInverseSqrt-Visualized).

## Why?

I like R! R doesn't have a type for 32 bit floats, so I wanted a way to mess with the FISR--a 32-bit specific implementation detail (ignoring for now the various extensions to other bases)--in R. Now you can mess with it too!

## Features

- Customizable parameters for fine-tuning accuracy and performance (or the reverse)
- Supports vectorized input
- C++ implementation for speed.
- Optional detailed output including initial approximation, intermediate steps, and error metrics

## Installation

Install the package from GitHub using the `devtools` package:

```R
# install.packages("devtools")
devtools::install_github("Protonk/frsrr")
```

### Usage

```R
library(frsrr)

# Custom parameters
result <- frsr(c(1, 4, 9, 16), magic = 0x5f375a86, NR = 2, A = 1.6, B = 0.6)
## result is a vector of length 4
print(result)
# [1] 0.9990148 0.4995074 0.3337626 0.2497537

# Optional detail
result.df <- frsr.detail(c(pi, 2^-31, 0.4, 6.02e23))
## result.df is a dataframe with 4 rows and 8 columns
print(result)
#          input      magic      initial    after_one        final        error          diff iters
# 1 3.141593e+00 1597463007 5.735160e-01 5.639570e-01 5.639570e-01 0.0004121269 -9.558976e-03     1
# 2 4.656613e-10 1597463007 4.693787e+04 4.632937e+04 4.632937e+04 0.0002499308 -6.085039e+02     1
# 3 4.000000e-01 1597463007 1.632430e+00 1.578616e+00 1.578616e+00 0.0015955754 -5.381417e-02     1
# 4 6.020000e+23 1597463007 1.306493e-12 1.288484e-12 1.288484e-12 0.0002824810 -1.800936e-14     1
```