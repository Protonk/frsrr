test_that("frsr detail produces expected output structure", {
  result <- frsr.detail(c(1, 4, 9))
  expect_s3_class(result, "data.frame")
  expect_named(result, c("input", "magic", "initial", "after_one", "final", "error", "diff", "iters"))
  expect_equal(nrow(result), 3)
})

test_that("frsr respects custom parameters", {
  result1 <- frsr(4, magic = 0x5f375a86, NR = 2)
  result2 <- frsr(4)
  expect_false(isTRUE(all.equal(result1, result2)))

  result3 <- frsr(4, A = 3.6, B = 0.2)
  expect_false(isTRUE(all.equal(result3, result2)))
})

test_that("frsr handles vector input", {
  result <- frsr.detail(c(1, 4, 9, 16))
  expect_equal(nrow(result), 4)
  expect_equal(result$final, 1/sqrt(c(1, 4, 9, 16)), tolerance = 1e-2)
  result.vec <- frsr(c(1, 4, 9, 16))
  expect_equal(length(result.vec), 4)
  expect_equal(result.vec, 1/sqrt(c(1, 4, 9, 16)), tolerance = 1e-2)
})

test_that("frsr.detail respects tolerance parameter", {
  result1 <- frsr.detail(0.0125, tol = 0, NR = 20)
  result2 <- frsr.detail(0.0125, tol = 1e-7, NR = 20)
  expect_true(result2$iters[1] < result1$iters[1])
})

test_that("frsr.detail produces accurate results with tolerance parameter", {
  result <- frsr.detail(4, tol = 1e-6, NR = 6)
  expect_equal(result$final[1], 1/sqrt(4), tolerance = 1e-6)
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
})


