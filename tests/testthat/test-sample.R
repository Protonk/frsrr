test_that("boundedStratifiedSample returns correct number of samples", {
  result <- boundedStratifiedSample(10, -2, 2)
  expect_length(result, 10)
})

test_that("boundedStratifiedSample throws error for non-numeric n", {
  expect_error(boundedStratifiedSample("ten", -2, 2), "'n' must be a positive integer")
})

test_that("boundedStratifiedSample throws error for non-positive n", {
  expect_error(boundedStratifiedSample(-10, -2, 2), "'n' must be a positive integer")
})

test_that("boundedStratifiedSample throws error for non-integer n", {
  expect_error(boundedStratifiedSample(10.5, -2, 2), "'n' must be a positive integer")
})

test_that("boundedStratifiedSample throws error for non-numeric low/high", {
  expect_error(boundedStratifiedSample(10, "low", 2), "'low' and 'high' must be numeric")
  expect_error(boundedStratifiedSample(10, -2, "high"), "'low' and 'high' must be numeric")
})

test_that("boundedStratifiedSample throws error for low < -126", {
  expect_error(boundedStratifiedSample(10, -127, 2), "Lower bound is too small. 'low' must be >= -126")
})

test_that("boundedStratifiedSample returns samples within bounds", {
  result <- boundedStratifiedSample(10, -2, 2)
  expect_true(all(result >= 2^(-2) & result < 2^2))
})
