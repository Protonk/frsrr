# Fast Inverse Square Root (FRSR) R Package

## Overview

This R package implements a parameterized Fast Reciprocal Square Root (FRSR) algorithm, also known as Fast Inverse Square Root (FISR), written in C++. You can read more about the FRSR at [0x5f37642f.com](https://0x5f37642f.com/) or see reasons why you might find instrumenting the output fun at [FastInverseSqrt-Visualized on GitHub](https://github.com/hyland-uw/FastInverseSqrt-Visualized).

## Features

- Customizable parameters for fine-tuning accuracy and performance
- Supports vectorized input
- Provides detailed output including initial approximation, intermediate steps, and error metrics

## Installation

Install the package from GitHub using the `devtools` package:

```R
# install.packages("devtools")
devtools::install_github("Protonk/frsrr")
```

### Usage

```R
library(frsr)

# Custom parameters
result <- frsr(c(1, 4, 9, 16), magic = 0x5f375a86, NR = 2, A = 1.6, B = 0.6)
print(result)
```