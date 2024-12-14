test_that("boundedStratifiedSample returns correct number of samples", {
    result <- boundedStratifiedSample(10, 2^-2, 2^2)
    expect_length(result, 10)
})

test_that("boundedStratifiedSample throws error for low < -126", {
    expect_error(boundedStratifiedSample(10, 2^-127, 2^2), "Subnormal numbers are not supported. 'low' must be >= -126")
})

test_that("boundedStratifiedSample returns samples within bounds", {
    result <- boundedStratifiedSample(100, 2^-2, 2^2)
    expect_true(all(result >= 2^-2 & result < 2^2))
})

test_that("sample_frsr returns correct number of samples", {
    result <- sample_frsr(4)
    expect_equal(nrow(result), 4)
})

test_that("sample_frsr returns correct parameters when keep_params is TRUE", {
    result <- sample_frsr(4, keep_params = TRUE)
    expect_true("A" %in% names(result) && "B" %in% names(result))
})
