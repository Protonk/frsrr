
test_that("sample_frsr generates correct number of samples", {
  n <- 10
  result <- sample_frsr(n)
  expect_equal(nrow(result), n)
})

test_that("sample_frsr respects input range", {
  n <- 10
  x_min <- 0.25
  x_max <- 1.0
  result <- sample_frsr(n, x_min = x_min, x_max = x_max)
  expect_true(all(result$input >= x_min & result$input <= x_max))
})

test_that("sample_frsr includes parameters when keep_params is TRUE", {
  n <- 10
  result <- sample_frsr(n, keep_params = TRUE)
  expect_true(all(c("magic", "NR", "A", "B", "tol") %in% colnames(result)))
})
