test_that("frsr_sample returns correct number of samples", {
    result <- frsr_sample(4)
    expect_equal(nrow(result), 4)
})

test_that("frsr_sample returns parameters when keep_params is TRUE", {
    result <- frsr_sample(4, keep_params = TRUE)
    detail_comp <- frsr(1:4, detail = TRUE)
    expect_true( all(names(detail_comp) %in% names(result)) )
})

test_that("frsr_sample handles NULL magic_min and magic_max correctly", {
    result <- frsr_sample(4, magic_min = NULL,keep_params = TRUE)
    expect_equal(length(unique(result$magic)), 1)

    result <- frsr_sample(4, magic_max = NULL, keep_params = TRUE)
    expect_equal(length(unique(result$magic)), 1)
})

test_that("frsr_sample handles NULL x_min and x_max correctly", {
    result <- frsr_sample(4, x_min = NULL)
    expect_equal(length(unique(result$input)), 1)

    result <- frsr_sample(4, x_max = NULL)
    expect_equal(length(unique(result$input)), 1)
})

test_that("frsr_sample returns correct structure", {
    result <- frsr_sample(4)
    expect_true(all(c("input", "initial", "after_one", "final", "error", "enre", "diff", "iters") %in% names(result)))
})

test_that("boundedStratifiedSample handles narrow exponent ranges", {
    low <- log2(0.75)
    high <- log2(1)

    samples <- .Call(
        "_frsrr_boundedStratifiedSample",
        PACKAGE = "frsrr",
        32L,
        low,
        high
    )

    expect_length(samples, 32)
    expect_true(all(is.finite(samples)))
    expect_true(all(samples >= 0.75))
    expect_true(all(samples < 1))
})

test_that("boundedStratifiedSample tolerates zero bit draws", {
    samples <- .Call(
        "_frsrr_boundedStratifiedSample",
        PACKAGE = "frsrr",
        1024L,
        log2(0.25),
        log2(0.5)
    )

    expect_length(samples, 1024)
    expect_true(all(is.finite(samples)))
})
