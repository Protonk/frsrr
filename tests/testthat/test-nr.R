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

test_that("customIteration works with default parameters", {
    x <- c(1, 4, 9, 16)
    formula <- quote(y * (1.5 - 0.5 * x * y^2))
    result <- customIteration(x, formula = formula)
    expected <- frsr(x)
    expect_equal(result, expected, tolerance = 1e-2)
})

test_that("customIteration works with multiple iterations", {
    x <- c(1, 4, 9, 16)
    formula <- quote(y * (1.5 - 0.5 * x * y^2))
    result <- customIteration(x, formula = formula, NR = 5)
    expected <- frsr(x, NR = 5)
    expect_equal(result, expected, tolerance = 1e-3)
})

test_that("customIteration respects tolerance", {
    x <- c(1, 4, 9, 16)
    formula <- quote(y * (1.5 - 0.5 * x * y^2))
    result <- customIteration(x, formula = formula, tol = 1e-5)
    expected <- frsr(x)
    expect_equal(result, expected, tolerance = 1e-4)
})
