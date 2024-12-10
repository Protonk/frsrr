test_that("frsr detail produces expected output structure", {
  result <- frsr.detail(c(1, 4, 9))
  expect_s3_class(result, "data.frame")
  expect_named(result, c("input", "magic", "initial", "after_one", "final", "error", "diff", "iters"))
  expect_equal(nrow(result), 3)
  # frsr
  result.vec <- frsr(c(1, 4, 9))
  expect_equal(length(result.vec), 3)
})

test_that("frsr respects custom parameters", {
  result1 <- frsr(4, magic = 0x5f375a86, NR = 2)
  result2 <- frsr(4)
  expect_false(isTRUE(all.equal(result1, result2)))

  result3 <- frsr(4, A = 3.6, B = 0.2)
  expect_false(isTRUE(all.equal(result3, result2)))
  # frsr.detail
  result1_detail <- frsr.detail(4, magic = 0x5f375a86, NR = 2)
  result2_detail <- frsr.detail(4)
  expect_false(isTRUE(all.equal(result1_detail$final, result2_detail$final)))

  result3_detail <- frsr.detail(4, A = 3.6, B = 0.2)
  expect_false(isTRUE(all.equal(result3_detail$final, result2_detail$final)))
})

test_that("frsr handles vector input", {
  result.vec <- frsr(c(1, 4, 9, 16))
  expect_equal(length(result.vec), 4)
  expect_equal(result.vec, 1/sqrt(c(1, 4, 9, 16)), tolerance = 1e-2)
  # frsr.detail
  result <- frsr.detail(c(1, 4, 9, 16))
  expect_equal(nrow(result), 4)
  expect_equal(result$final, 1/sqrt(c(1, 4, 9, 16)), tolerance = 1e-2)
})

test_that("frsr respects tolerance parameter", {
  result1 <- frsr.detail(0.0125, tol = 0, NR = 20)
  result2 <- frsr.detail(0.0125, tol = 1e-7, NR = 20)
  expect_true(result2$iters[1] < result1$iters[1])
  # frsr
  result1_vec <- frsr(0.0125, tol = 0, NR = 20)
  result2_vec <- frsr(0.0125, tol = 1e-3, NR = 20)
  ## The results should differ by more than 1e-4 due to the tolerance parameter
  expect_true(abs(result2_vec - result1_vec) > 1e-4)
})

test_that("frsr produces accurate results with tolerance parameter", {
  result <- frsr.detail(4, tol = 1e-6, NR = 6)
  expect_equal(result$final[1], 1/sqrt(4), tolerance = 1e-6)
  # frsr
  result_vec <- frsr(4, tol = 1e-6, NR = 6)
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
  # frsr.detail
  result_single_detail <- frsr.detail(4)
  expect_s3_class(result_single_detail, "data.frame")
  expect_equal(nrow(result_single_detail), 1)
  
  result_multiple_detail <- frsr.detail(input_vector)
  expect_s3_class(result_multiple_detail, "data.frame")
  expect_equal(nrow(result_multiple_detail), length(input_vector))
})
