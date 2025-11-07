test_that("frsr_bin returns a data frame", {
  result <- frsr_bin()
  expect_s3_class(result, "data.frame")
})

test_that("frsr_bin handles different number of bins", {
  result <- frsr_bin(n_bins = 2, float_samples = 10, magic_samples = 10)
  expect_equal(nrow(result), 2)
  result <- frsr_bin(n_bins = 8, float_samples = 10, magic_samples = 10)
  expect_equal(nrow(result), 8)
})

test_that("frsr_bin validates parameters", {
  expect_error(
    frsr_bin(n_bins = 0),
    "`n_bins` must be at least 1",
    fixed = TRUE
  )

  expect_error(
    frsr_bin(magic_min = 1596980100L, magic_max = 1596980000L),
    "`magic_min` must be less than or equal to `magic_max`",
    fixed = TRUE
  )
})

test_that("frsr_bin returns expected bin metadata and magic bounds", {
  x_min <- 0.25
  x_max <- 1
  n_bins <- 4
  magic_min <- 1596980000L
  magic_max <- 1596980100L

  result <- frsr_bin(
    x_min = x_min,
    x_max = x_max,
    n_bins = n_bins,
    float_samples = 32,
    magic_samples = 32,
    magic_min = magic_min,
    magic_max = magic_max
  )

  edges <- seq(x_min, x_max, length.out = n_bins + 1)
  expect_equal(result$Range_Min, head(edges, -1))
  expect_equal(result$Range_Max, tail(edges, -1))
  expect_true(all(result$Magic >= magic_min))
  expect_true(all(result$Magic <= magic_max))
})
