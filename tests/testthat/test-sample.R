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
    expect_true(all(c("input", "initial", "after_one", "final", "error", "diff", "iters") %in% names(result)))
})
