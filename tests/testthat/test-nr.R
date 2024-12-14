test_that("ynplusone works with standard formula", {
    x <- c(1, 4, 9, 16)
    guess <- c(1, 0.5, 0.333, 0.25)
    formula <- quote(y * (1.5 - 0.5 * x * y^2))
    result <- ynplusone(x, guess, formula)
    expect_equal(result, frsr(x), tolerance = 1e-2)
})

test_that("ynplusone works with different formulae", {
    x <- c(1, 4, 9, 16)
    guess <- c(1, 0.5, 0.333, 0.25)
    ## not likely to be accurate, but should return a result
    wrong <- ynplusone(x, guess, quote(y * (1.4 - 0.4 * x - y^2)))
    ## just to make sure we're not bamboozling ourselves
    std_formula <- quote(y * (1.5 - 0.5 * x * y^2))
    std_result <- ynplusone(x, guess, std_formula)
    expect_false(isTRUE(abs(wrong - std_result) < 1e-3))
})

test_that("customNR works with default parameters", {
    x <- c(1, 4, 9, 16)
    formula <- quote(y * (1.5 - 0.5 * x * y^2))
    result <- customNR(x, formula = formula)
    expect_equal(result$final, frsr(x), tolerance = 1e-2)
})

test_that("customNR works with multiple iterations", {
    x <- c(1, 4, 9, 16)
    formula <- quote(y * (1.5 - 0.5 * x * y^2))
    result <- customNR(x, formula = formula, NR = 5)
    expect_equal(result$final, frsr(x, NR = 5), tolerance = 1e-3)
})

test_that("customNR produces expected output structure", {
    x <- c(1, 4, 9, 16)
    formula <- quote(y * (1.5 - 0.5 * x * y^2))
    result <- customNR(x, formula = formula)
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), length(x))
})

