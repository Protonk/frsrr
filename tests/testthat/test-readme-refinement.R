# Run the retained calculation once, from the checkout or the installed package.
# Sourcing happens inside the first test that asks, so errors name a test.
readme_experiment <- local({
  experiment <- NULL
  function() {
    if (is.null(experiment)) {
      path <- testthat::test_path("..", "..", "inst", "examples", "readme-refinement.R")
      if (!file.exists(path)) {
        path <- system.file("examples", "readme-refinement.R", package = "frsrr",
                            mustWork = TRUE)
      }
      env <- new.env(parent = globalenv())
      sys.source(path, envir = env)
      experiment <<- env
    }
    experiment
  }
})

# Fenced blocks in README.md: info string, contents, and enclosing "## " heading.
readme_blocks <- function(path) {
  lines <- readLines(path, warn = FALSE)
  fences <- grep("^```", lines)
  starts <- fences[c(TRUE, FALSE)]
  ends <- fences[c(FALSE, TRUE)]
  headings <- grep("^## ", lines)
  lapply(seq_along(starts), function(i) {
    list(lang = sub("^```", "", lines[starts[i]]),
         text = lines[starts[i] + seq_len(ends[i] - starts[i] - 1)],
         section = lines[max(headings[headings < starts[i]])])
  })
}

# Evaluate code as typed at the console and return what it prints. Plots go to
# a null device; digits and width are fixed so output does not depend on options.
run_console <- function(code, envir) {
  old <- options(digits = 7, width = 80)
  grDevices::pdf(NULL)
  on.exit({
    grDevices::dev.off()
    options(old)
  })
  utils::capture.output(for (expr in parse(text = code)) {
    result <- withVisible(eval(expr, envir))
    if (result$visible) print(result$value)
  })
}

test_that("README code runs and its output blocks match the printed results", {
  readme <- testthat::test_path("..", "..", "README.md")
  skip_if_not(file.exists(readme), "README.md is only present in the source checkout")

  # An untyped block is the output of the R block just before it.
  env <- new.env(parent = globalenv())
  printed <- NULL
  n_outputs <- 0L
  for (block in readme_blocks(readme)) {
    if (block$section == "## Installation") next
    if (block$lang == "R") {
      printed <- run_console(block$text, env)
    } else if (block$lang == "") {
      expect_false(is.null(printed), label = "an output block follows an R block")
      expect_identical(block$text, printed)
      printed <- NULL
      n_outputs <- n_outputs + 1L
    }
  }
  expect_gt(n_outputs, 0L)

  # The README code and the retained script must build the same experiment.
  script <- readme_experiment()
  for (name in c("x_scan", "ordinary_scan", "candidates", "scan_max", "A_adjusted",
                 "d", "x_eval", "steps", "ordinary", "adjusted", "comparison",
                 "predicted")) {
    expect_identical(env[[name]], script[[name]], label = name)
  }
  for (name in c("float32", "measure", "max_error", "plot_refinement")) {
    expect_identical(deparse(env[[name]]), deparse(script[[name]]), label = name)
  }
})

test_that("README inputs and diagnostics follow the numerical contract", {
  e <- readme_experiment()
  for (x in list(e$x_scan, e$x_eval)) {
    expect_identical(float32(x), x)
    expect_true(all(x >= 1 & x < 4))
  }
  expect_true(all(e$x_scan %in% e$x_eval))
  expect_gt(length(setdiff(e$x_eval, e$x_scan)), length(e$x_scan))
  expect_identical(float32(e$A_adjusted), e$A_adjusted)

  # Check supplied parameters and the independently calculated signed errors.
  check_fit <- function(fit, x, A, steps) {
    expect_identical(fit$input, x)
    expect_true(all(fit$magic == 0x5f3759df & fit$B == 0.5 & fit$tol == 0))
    expect_true(all(fit$A == A & fit$NRmax == steps & fit$iters == steps))
    reference <- 1 / sqrt(x)
    signed <- (fit$final - reference) / reference
    expect_equal(fit$signed_error, signed, tolerance = 1e-14)
    expect_lt(max(abs(abs(signed) - fit$error)), 1e-14)
  }
  check_fit(e$ordinary_scan, e$x_scan, 1.5, 1L)
  for (i in seq_along(e$steps)) {
    check_fit(e$ordinary[[i]], e$x_eval, 1.5, e$steps[i])
    check_fit(e$adjusted[[i]], e$x_eval, e$A_adjusted, e$steps[i])
  }
})

test_that("README one-step and later-step comparisons hold on the stated grids", {
  e <- readme_experiment()
  expect_true(1.5 %in% e$candidates)
  expect_gt(e$A_adjusted, 1.5)
  expect_lt(min(e$scan_max), max(e$ordinary_scan$error))
  expect_gt(mean(e$ordinary_scan$signed_error < 0), 0.998)
  expect_lt(min(e$adjusted[[1]]$signed_error), 0)
  expect_gt(max(e$adjusted[[1]]$signed_error), 0)
  expect_true(any(e$adjusted[[1]]$error > e$ordinary[[1]]$error))

  # The README's "26%" and "above about -d/2" description of worsened inputs.
  worse <- e$adjusted[[1]]$error > e$ordinary[[1]]$error
  signed <- e$ordinary[[1]]$signed_error
  expect_equal(round(100 * mean(worse)), 26)
  expect_gt(min(signed[worse]), -e$d / 2 - 2e-5)
  expect_lt(max(signed[!worse]), -e$d / 2 + 2e-5)
  expect_true(any(worse & abs(e$x_eval - 2) < 0.05))

  expect_lt(e$comparison$adjusted[1], e$comparison$ordinary[1])
  expect_true(all(e$comparison$ordinary[-1] < e$comparison$adjusted[-1] / 100))
  expect_true(all(diff(e$comparison$adjusted) < 0))
  # Allow several float32 epsilons, rather than pinning incidental last bits.
  expect_lt(e$comparison$ordinary[3], 3 * 2^-23)
  expect_true(all(e$adjusted[[3]]$signed_error > 0))
  expect_lt(max(abs(e$adjusted[[3]]$signed_error - e$predicted)), 3 * 2^-23)
})
