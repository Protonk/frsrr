describe("frsr_phase", {
  it("returns expected structure", {
    set.seed(11)
    result <- frsr_phase(
      phases = 8L,
      exponents = -2L:2L,
      per_cell = 4L,
      magics = c(1597413411L, 1597200000L),
      q = 0.9
    )

    expect_type(result$magic, "integer")
    expect_type(result$J, "double")
    expect_type(result$R, "double")
    expect_s3_class(result$phase_tbl, "data.frame")
    expect_true(is.matrix(result$heat))
    expect_equal(nrow(result$phase_tbl), 8)
    expect_equal(result$phase_tbl$phase_id, seq_len(8))
    expect_true(all(result$phase_tbl$n == length(-2L:2L) * 4L))
    expect_equal(nrow(result$heat), length(-2L:2L))
    expect_equal(ncol(result$heat), 8)
  })

  it("worsens objective for clearly bad magic constants", {
    set.seed(99)
    good <- frsr_phase(
      phases = 4L,
      exponents = -3L:3L,
      per_cell = 6L,
      magics = 0x5f3759df
    )
    set.seed(99)
    bad <- frsr_phase(
      phases = 4L,
      exponents = -3L:3L,
      per_cell = 6L,
      magics = 0x5f100000
    )

    expect_lt(good$J, bad$J)
  })

  it("is reproducible with a fixed seed", {
    args <- list(
      phases = 6L,
      exponents = -2L:2L,
      per_cell = 4L,
      magics = c(1597413411L, 1597200000L),
      q = 0.8
    )
    set.seed(2024)
    first <- do.call(frsr_phase, args)
    set.seed(2024)
    second <- do.call(frsr_phase, args)
    expect_identical(first, second)
  })

  it("validates arguments before invoking C++ helpers", {
    expect_error(
      frsr_phase(phases = 0L),
      "`phases` must be a positive integer",
      fixed = TRUE
    )
    expect_error(
      frsr_phase(exponents = integer()),
      "`exponents` must contain at least one value",
      fixed = TRUE
    )
    expect_error(
      frsr_phase(exponents = 200L),
      "`exponents` must stay within [-126, 127]",
      fixed = TRUE
    )
    expect_error(
      frsr_phase(per_cell = 0L),
      "`per_cell` must be a positive integer",
      fixed = TRUE
    )
    expect_error(
      frsr_phase(magics = integer()),
      "`magics` must supply at least one integer constant",
      fixed = TRUE
    )
    expect_error(
      frsr_phase(q = 0),
      "`q` must satisfy 0 < q <= 1",
      fixed = TRUE
    )
    expect_error(
      frsr_phase(q = 2),
      "`q` must satisfy 0 < q <= 1",
      fixed = TRUE
    )
    expect_error(
      frsr_phase(NRmax = -1L),
      "`NRmax` must be a non-negative integer",
      fixed = TRUE
    )
  })
})

test_that('phase comparisons reuse each candidate experiment across candidate changes', {
    args <- list(phases = 6, exponents = -2:2, per_cell = 8, NRmax = 1, q = 0.8)
    magics <- as.integer(c(0x5f3759df, 0x5f375a86, 0x5f100000))
    run <- function(m) {
        set.seed(2024)
        do.call(frsr_phase, c(args, list(magics = m)))
    }
    individual <- lapply(magics, run)
    best <- order(vapply(individual, `[[`, 0, 'J'),
                  vapply(individual, `[[`, 0, 'R'), magics)[1]
    fit <- run(magics)
    fields <- c('magic', 'J', 'R', 'phase_tbl', 'heat')
    expect_identical(fit[fields], individual[[best]][fields])
    expect_identical(fit, run(rev(magics)))
    expect_identical(fit[fields], run(c(0L, magics))[fields])
    expect_identical(fit[fields], run(c(magics, magics))[fields])
    expect_identical(fit$settings$NRmax, 1L)
    expect_identical(fit$settings$magics, sort(magics))
})

test_that('phase statistics agree with independently reconstructed draws and evaluation', {
    set.seed(41)
    x <- float32(2^runif(64))
    fit <- frsr(x, NRmax = 2, tol = 0, threads = 1, detail = TRUE)
    signed <- (fit$final - 1 / sqrt(x)) / (1 / sqrt(x))
    set.seed(41)
    phase <- frsr_phase(phases = 1, exponents = 0, per_cell = 64,
                        magics = 0x5f3759df, q = 0.9, NRmax = 2)
    expect_equal(phase$J, unname(quantile(abs(signed), 0.9)), tolerance = 1e-15)
    expect_equal(phase$phase_tbl$mean_signed, mean(signed), tolerance = 1e-15)
    expect_identical(phase$phase_tbl$median_signed, median(signed))
    expect_identical(unname(phase$heat[1, 1]), median(signed))
})

test_that('phase ties are exact and normal exponent extremes remain supported', {
    ties <- as.integer(c(0x20000001, 0x20000000))
    set.seed(9)
    fit <- frsr_phase(phases = 1, exponents = 0, per_cell = 4,
                      magics = ties, NRmax = 1)
    set.seed(9)
    other <- frsr_phase(phases = 1, exponents = 0, per_cell = 4,
                        magics = rev(ties), NRmax = 1)
    expect_identical(fit, other)
    expect_identical(fit$magic, min(ties))
    set.seed(9)
    edge <- frsr_phase(phases = 4, exponents = c(-126, 127), per_cell = 32,
                       magics = 0x5f3759df, NRmax = 2)
    expect_true(is.finite(edge$J))
    expect_true(all(is.finite(edge$heat)))
})
