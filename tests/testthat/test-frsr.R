test_that("frsr produces expected output structure", {
  result <- frsr(c(1, 4, 9))
  expect_s3_class(result, "data.frame")
  expect_named(result, c("input", "magic", "initial", "after_one", "final", "error", "diff", "iters"))
  expect_equal(nrow(result), 3)
})

test_that("frsr respects custom parameters", {
  result1 <- frsr(4, magic = 0x5f375a86, NR = 2)
  result2 <- frsr(4)
  expect_false(isTRUE(all.equal(result1$final, result2$final)))

  result3 <- frsr(4, A = 3.6, B = 0.2)
  expect_false(isTRUE(all.equal(result3$final, result2$final)))
})

test_that("frsr handles vector input", {
  result <- frsr(c(1, 4, 9, 16))
  expect_equal(nrow(result), 4)
  expect_equal(result$final, 1/sqrt(c(1, 4, 9, 16)), tolerance = 1e-2)
})

test_that("frsr respects tolerance parameter", {
  result1 <- frsr(0.0125, tol = 0, NR = 20)
  result2 <- frsr(0.0125, tol = 1e-7, NR = 20)
  expect_true(result2$iters[1] < result1$iters[1])
})

test_that("frsr produces accurate results with tolerance parameter", {
  result <- frsr(4, tol = 1e-6, NR = 6)
  expect_equal(result$final[1], 1/sqrt(4), tolerance = 1e-6)
})




