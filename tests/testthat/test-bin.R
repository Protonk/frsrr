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

