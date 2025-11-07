test_that("detail argument produces expected output structure", {
  result <- frsr(c(1, 4, 9), detail = TRUE)
  expected_cols <- c(
    "input",
    "initial",
    "after_one",
    "final",
    "error",
    "enre",
    "diff",
    "iters"
  )

  expect_s3_class(result, "data.frame")
  expect_identical(names(result), expected_cols)
  expect_equal(nrow(result), 3)
  expect_true(all(vapply(result, is.numeric, logical(1))))

  result.vec <- frsr(c(1, 4, 9))
  expect_equal(length(result.vec), 3)
})

test_that("frsr respects custom parameters", {
  result1 <- frsr(4, magic = 0x5f375a86, NRmax = 2)
  result2 <- frsr(4)
  expect_false(isTRUE(all.equal(result1, result2)))

  result3 <- frsr(4, A = 3.6, B = 0.2)
  expect_false(isTRUE(all.equal(result3, result2)))

  result3_detail <- frsr(4, A = 3.6, B = 0.2, detail = TRUE)
  expect_false(isTRUE(all.equal(result3_detail$final, result1)))
})

test_that("frsr handles vector input", {
  result.vec <- frsr(c(1, 4, 9, 16))
  expect_equal(length(result.vec), 4)
  expect_equal(result.vec, 1/sqrt(c(1, 4, 9, 16)), tolerance = 1e-2)
  # frsr
  result <- frsr(c(1, 4, 9, 16), detail = TRUE)
  expect_equal(nrow(result), 4)
  expect_equal(result$final, 1/sqrt(c(1, 4, 9, 16)), tolerance = 1e-2)
})

test_that("frsr respects tolerance parameter", {
  # iterates to NRmax with tol = 0
  result1 <- frsr(0.0125, tol = 0, NRmax = 20, detail = TRUE)
  result2 <- frsr(0.0125, tol = 1e-7, NRmax = 20, detail = TRUE)
  expect_true(result2$iters[1] < result1$iters[1])
})

test_that("frsr produces accurate results with tolerance parameter", {
  result_vec <- frsr(4, tol = 1e-6, NRmax = 6)
  expect_equal(result_vec, 1/sqrt(4), tolerance = 1e-6)
})

test_that("frsr default returns a numeric vector", {
  # Test with a single value
  result_single <- frsr(4)
  expect_type(result_single, "double")
  expect_length(result_single, 1)

  # Test with multiple values
  input_vector <- c(1, 4, 9, 16)
  result_multiple <- frsr(input_vector)
  expect_type(result_multiple, "double")
  expect_length(result_multiple, length(input_vector))

  # Check that the output is indeed a vector and not a list or data frame
  expect_true(is.vector(result_multiple))
  expect_false(is.list(result_multiple))
  expect_false(is.data.frame(result_multiple))
  # with detail param set to TRUE we expect a data frame
  result_single_detail <- frsr(4, detail = TRUE)
  expect_s3_class(result_single_detail, "data.frame")
  expect_equal(nrow(result_single_detail), 1)

  result_multiple_detail <- frsr(input_vector, detail = TRUE)
  expect_s3_class(result_multiple_detail, "data.frame")
  expect_equal(nrow(result_multiple_detail), length(input_vector))
})

test_that("frsr toggles parameter columns according to keep_params", {
  base_cols <- c(
    "input",
    "initial",
    "after_one",
    "final",
    "error",
    "enre",
    "diff",
    "iters"
  )
  param_cols <- c("magic", "NRmax", "A", "B", "tol")

  result_with <- frsr(4, keep_params = TRUE, detail = TRUE)
  expect_identical(names(result_with), c(base_cols, param_cols))

  result_without <- frsr(4, keep_params = FALSE, detail = TRUE)
  expect_identical(names(result_without), base_cols)
})
