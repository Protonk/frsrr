test_that("frsr_bin returns documented columns", {
  result <- frsr_bin(float_samples = 8, magic_samples = 8)
  expected_cols <- c(
    "N_bins",
    "Location",
    "Range_Min",
    "Range_Max",
    "Magic",
    "Objective_Max",
    "Metric_Mean"
  )

  expect_s3_class(result, "data.frame")
  expect_identical(names(result), expected_cols)
  expect_true(all(vapply(result, is.numeric, logical(1))))
  expect_identical(unique(result$N_bins), 4L)
})

test_that("frsr_bin supports weighted sampling", {
  result <- frsr_bin(float_samples = 8, magic_samples = 8, weighted = TRUE)
  expect_s3_class(result, "data.frame")
  expect_true(all(is.finite(unlist(result))))
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

  expect_error(
    frsr_bin(metrics_to_report = character()),
    "`metrics_to_report` must contain at least one metric name",
    fixed = TRUE
  )

  expect_error(
    frsr_bin(objective_metric = ""),
    "`objective_metric` must not be empty",
    fixed = TRUE
  )
})

test_that("frsr_bin supports configurable metrics", {
  result <- frsr_bin(
    float_samples = 8,
    magic_samples = 8,
    objective_metric = "rms",
    metrics_to_report = c("sum", "mean", "sum")
  )

  expected_metrics <- c("Magic", "Objective_RootMeanSquare", "Metric_Sum", "Metric_Mean")
  expect_identical(names(result)[5:8], expected_metrics)
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
