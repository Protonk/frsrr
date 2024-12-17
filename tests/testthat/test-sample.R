test_that("frsr_sample returns correct number of samples", {
    result <- frsr_sample(4)
    expect_equal(nrow(result), 4)
})

test_that("frsr_sample returns parameters when keep_params is TRUE", {
    result <- frsr_sample(4, keep_params = TRUE)
    detail_comp <- frsr(1:4, detail = TRUE)
    expect_true( all(names(detail_comp) %in% names(result)) )
})
