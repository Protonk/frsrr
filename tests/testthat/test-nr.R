# FILE: tests/testthat/test-nr.R

test_that("frsr_NR works with standard formula", {
    x <- c(1, 4, 9, 16)
    guess <- c(1, 0.5, 0.333, 0.25)
    formula <- quote(y * (1.5 - 0.5 * x * y^2))
    result <- frsr_NR(x, formula = formula, NRmax = 1)
    expect_equal(result$final, frsr(x), tolerance = 1e-2)
})

test_that("frsr_NR works with different formulae", {
    x <- c(1, 4, 9, 16)
    ## not likely to be accurate, but should return a result
    wrong <- frsr_NR(x, formula = quote(y * (1.4 - 0.4 * x - y^2)), NRmax = 1)
    ## just to make sure we're not bamboozling ourselves
    std_formula <- quote(y * (1.5 - 0.5 * x * y^2))
    std_result <- frsr_NR(x, formula = std_formula, NRmax = 1)
    expect_false(isTRUE(abs(wrong$final - std_result$final) < 1e-3))
})

test_that("frsr_NR works with default parameters", {
    x <- c(1, 4, 9, 16)
    formula <- quote(y * (1.5 - 0.5 * x * y^2))
    result <- frsr_NR(x, formula = formula)
    expect_equal(result$final, frsr(x), tolerance = 1e-2)
})

test_that("frsr_NR works with multiple iterations", {
    x <- c(1, 4, 9, 16)
    formula <- quote(y * (1.5 - 0.5 * x * y^2))
    result <- frsr_NR(x, formula = formula, NRmax = 5)
    expect_equal(result$final, frsr(x, NRmax = 5), tolerance = 1e-3)
})

test_that("frsr_NR produces expected output structure", {
    x <- c(1, 4, 9, 16)
    formula <- quote(y * (1.5 - 0.5 * x * y^2))
    result <- frsr_NR(x, formula = formula)
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), length(x))
})