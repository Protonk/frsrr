test_that("boundedStratifiedSample throws error for low < -126", {
    expect_error(boundedStratifiedSample(10, 2^-127, 2^2), "Subnormal numbers are not supported. 'low' must be >= -126")
})

test_that("boundedStratifiedSample returns samples within bounds", {
    result <- boundedStratifiedSample(100, 2^-2, 2^2)
    expect_true(all(result >= 2^-2 & result < 2^2))
})

test_that("frsr_sample returns correct number of samples", {
    result <- frsr_sample(4)
    expect_equal(nrow(result), 4)
})

test_that("frsr_sample returns parameters when keep_params is TRUE", {
    result <- frsr_sample(4, keep_params = TRUE)
    detail_comp <- frsr(1:4, detail = TRUE)
    expect_true( all(names(detail_comp) %in% names(result)) )
})
