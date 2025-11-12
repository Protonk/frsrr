describe("frsr_bin", {
  it("returns documented columns", {
    result <- frsr_bin(float_samples = 8, magic_samples = 8)
    expected_cols <- c(
      "N_bins",
      "Location",
      "Range_Min",
      "Range_Max",
      "Magic",
      "Objective",
      "Dependent"
    )

    expect_s3_class(result, "data.frame")
    expect_identical(names(result), expected_cols)
    expect_true(all(vapply(result, is.numeric, logical(1))))
    expect_identical(unique(result$N_bins), 4L)
  })

  it("supports weighted sampling", {
    result <- frsr_bin(float_samples = 8, magic_samples = 8, weighted = TRUE)
    expect_s3_class(result, "data.frame")
    expect_true(all(is.finite(unlist(result))))
  })

  it("allows metric overrides", {
    result <- frsr_bin(
      float_samples = 8,
      magic_samples = 8,
      objective = "rmse_relative_error",
      dependent = "avg_relative_error"
    )
    expect_s3_class(result, "data.frame")
    expect_true(all(is.finite(result$Objective)))
    expect_true(all(is.finite(result$Dependent)))
  })

  it("handles different number of bins", {
    result <- frsr_bin(n_bins = 2, float_samples = 10, magic_samples = 10)
    expect_equal(nrow(result), 2)
    result <- frsr_bin(n_bins = 8, float_samples = 10, magic_samples = 10)
    expect_equal(nrow(result), 8)
  })

  it("only validates arguments needed for C++ bridge", {
    empty <- frsr_bin(n_bins = 0)
    expect_equal(nrow(empty), 0)

    expect_error(
      frsr_bin(objective = "bogus"),
      "'arg' should be one of",
      fixed = TRUE
    )
  })

  it("allows reversed magic bounds", {
    result <- frsr_bin(
      float_samples = 8,
      magic_samples = 8,
      magic_min = 1596980100L,
      magic_max = 1596980000L
    )
    expect_s3_class(result, "data.frame")
  })

  it("returns expected bin metadata and magic bounds", {
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

  it("trims vector arguments to scalars before calling C++ helpers", {
    result <- frsr_bin(
      n_bins = c(3, 8),
      float_samples = c(16, 32),
      magic_samples = c(32, 64),
      magic_min = c(1596980000L, 1596980500L),
      magic_max = c(1596980100L, 1596980600L)
    )

    expect_equal(nrow(result), 3)
    expect_identical(unique(result$N_bins), 3L)
  })

  it("is reproducible with matching seeds", {
    args <- list(
      float_samples = 8,
      magic_samples = 8,
      magic_min = 1596980000L,
      magic_max = 1596980100L
    )
    set.seed(77)
    first <- do.call(frsr_bin, args)
    set.seed(77)
    second <- do.call(frsr_bin, args)
    expect_identical(first, second)
  })
})
