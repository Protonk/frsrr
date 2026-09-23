test_that('seed and refinements match independently rounded binary32 fixtures', {
    # Seed bits = 0x5f3759df - floor(input_bits / 2). Refinements were checked
    # with exact rational products/subtraction, rounded nearest-even after each
    # operation. Expected outputs are stored as bit patterns to avoid decimal parsing error.
    # Neighbors around 1 and 2 distinguish rounding/evaluation-order choices.
    bits <- c(0x3f7fffff, 0x3f800000, 0x3f800001, 0x3fffffff,
              0x40000000, 0x40000001, 0x3dcccccd, 0x40490fdb,
              0x00800000, 0x7f000000, 0x7f7fffff)
    seed_bits <- c(0x3f7759e0, 0x3f7759df, 0x3f7759df, 0x3f3759e0,
                   0x3f3759df, 0x3f3759df, 0x4050f379, 0x3f12d1f2,
                   0x5ef759df, 0x1fb759df, 0x1f7759e0)
    one <- float_from_bits(c(0x3f7f9110, 0x3f7f910f, 0x3f7f910d, 0x3f34f95e,
                             0x3f34f95e, 0x3f34f95e, 0x404a1017, 0x3f105f7d,
                             0x5eff910f, 0x1fb4f95e, 0x1f7f9110))
    two <- float_from_bits(c(0x3f7fffb8, 0x3f7fffb7, 0x3f7fffb7, 0x3f3504f3,
                             0x3f3504f1, 0x3f3504f1, 0x404a628f, 0x3f106eb8,
                             0x5effffb7, 0x1fb504f1, 0x1f7fffb8))
    x <- float_from_bits(bits)
    seed <- frsr(x, NRmax = 0, tol = 0, threads = 1, detail = TRUE)
    fit <- frsr(x, NRmax = 2, tol = 0, threads = 1, detail = TRUE)
    expect_identical(seed$initial, float_from_bits(seed_bits))
    expect_identical(seed$final, seed$initial)
    expect_identical(seed$iters, rep(0L, length(x)))
    expect_true(all(is.na(seed$after_one) & is.na(seed$diff)))
    expect_identical(fit$after_one, one)
    expect_identical(fit$final, two)
    expect_identical(fit$diff, two - one)
})

test_that('the error reference targets the rounded input in double precision', {
    x <- c(1 + 2^-25, 0.1, pi)
    fit <- frsr(x, NRmax = 1, threads = 1, detail = TRUE)
    reference <- 1 / sqrt(float32(x))
    expect_identical(fit$input, x)
    expect_identical(fit$error, abs((fit$final - reference) / reference))
    original_error <- abs((fit$final - 1 / sqrt(x)) / (1 / sqrt(x)))
    expect_true(all(fit$error != original_error))
    expect_identical(fit$final, frsr(float32(x), NRmax = 1, threads = 1))
    # The first fixture distinguishes mixed-double refinement from float32 steps.
    mixed <- float32(fit$initial * (1.5 - 0.5 * x * fit$initial * fit$initial))
    expect_true(any(mixed != fit$final))
})

test_that('coefficients round to float32 and vector tolerance stops after a step', {
    x <- c(1, 2, pi)
    A <- 1.5 + 2^-26
    B <- 0.5 + 2^-28
    expect_identical(frsr(x, A = A, B = B, threads = 1),
                     frsr(x, A = float32(A), B = float32(B), threads = 1))
    fit <- frsr(x, NRmax = 3, tol = c(1, 0, 1), threads = 1, detail = TRUE)
    expect_identical(fit$iters, c(1L, 3L, 1L))
    expect_identical(fit$final[c(1, 3)], fit$after_one[c(1, 3)])
})

test_that('unsupported inputs fail before worker computations', {
    for (x in c(0, -1, NA_real_, NaN, Inf, 2^-149, 2^-150, 2^128)) {
        expect_error(frsr(x, threads = 1), 'positive normal float32')
        expect_error(search_values(x, 0x5f3759df), 'positive normal float32')
    }
    expect_error(frsr(1, NRmax = -1), 'non-negative')
    expect_error(frsr(1, NRmax = 0.5), 'integers')
    expect_error(frsr_bin(NRmax = 0.5), 'integer')
    expect_error(frsr_phase(NRmax = 0.5), 'integer')
    expect_error(frsr(1, magic = NA_integer_), 'non-missing')
    expect_error(frsr(1, A = Inf), 'finite')
    expect_error(frsr(1, B = 2^128), 'float32')
    expect_error(frsr(1, tol = NA_real_), 'finite')
    expect_error(frsr(1, tol = -1), 'non-negative')
    expect_length(frsr(numeric(), threads = 1), 0)
})

test_that('bad exploratory parameters remain observable but cannot win searches', {
    bad <- frsr(1, magic = 0L, NRmax = 2, detail = TRUE, threads = 1)
    expect_false(is.finite(bad$final))
    expect_identical(bad$error, Inf)
    expect_error(search_values(1, 0L, NRmax = 2), 'No candidate')
    good <- search_values(1, c(0L, 0x5f3759df), NRmax = 2)
    expect_identical(good$Magic, as.integer(0x5f3759df))
    expect_true(is.finite(frsr(1, A = 3.6, B = -0.2, threads = 1)))
    set.seed(21)
    expect_error(frsr_phase(phases = 2, exponents = 0, per_cell = 2,
                           magics = 0L, NRmax = 2), 'No candidate')
})

test_that('ordinary powers-of-four scaling preserves the seed and refinement', {
    # Interior normal inputs/intermediates, default magic and coefficients only;
    # the property is not assumed through underflow/overflow or arbitrary magics.
    x <- c(1, 1.25, 1.75, 2, 3)
    for (steps in 0:2) {
        base <- frsr(x, NRmax = steps, threads = 1)
        for (k in c(-20, -1, 1, 20)) {
            expect_identical(frsr(x * 4^k, NRmax = steps, threads = 1) * 2^k, base)
        }
    }
    fit <- frsr(c(1, 2), NRmax = 0, detail = TRUE, threads = 1)
    expect_false('enre' %in% names(fit))
    expect_true(fit$error[1] != fit$error[2])
})

test_that('evaluation and all search metrics describe the same approximation', {
    x <- c(0.1, 1 + 2^-25, pi, 1 - 2^-24, 2 + 2^-22, 2^-126, 2^127)
    magics <- c(0x5f3759df, 0x5f375a86)
    metrics <- c('max_relative_error', 'avg_relative_error', 'rmse_relative_error')
    for (steps in 0:2) {
        for (magic in magics) {
            error <- frsr(x, magic = magic, NRmax = steps, tol = 0,
                          detail = TRUE, threads = 1)$error
            expected <- c(max(error), mean(error), sqrt(mean(error^2)))
            for (i in seq_along(metrics)) {
                found <- search_values(x, magic, steps, metrics[i], metrics[i])
                expect_equal(found$Objective, expected[i], tolerance = 1e-15)
                expect_identical(found$Objective, found$Dependent)
            }
        }
    }
    # Force exact objective ties: both near-zero results have double error 1.
    ties <- as.integer(c(0x1fc00000, 0x1fc00001))
    a <- search_values(1, ties, NRmax = 1)
    b <- search_values(1, rev(ties), NRmax = 1)
    expect_identical(a, b)
    expect_identical(a$Magic, min(ties))
})

test_that('fixed reduction order is identical across thread counts', {
    set.seed(918)
    # 64 blocks and several candidates: substantially more work than a serial
    # single-block case, with non-exact sums and close candidate rankings.
    x <- exp(runif(131072, log(0.25), log(4)))
    magics <- as.integer(0x5f3759df + c(-1, 0, 1, 2, 167))
    for (metric in c('max_relative_error', 'avg_relative_error', 'rmse_relative_error')) {
        RcppParallel::setThreadOptions(numThreads = 1L)
        one <- search_values(x, magics, 1, metric, 'rmse_relative_error')
        RcppParallel::setThreadOptions(numThreads = 2L)
        two <- search_values(x, magics, 1, metric, 'rmse_relative_error')
        expect_identical(one, two)
        expect_identical(two, search_values(x, rev(magics), 1, metric, 'rmse_relative_error'))
    }
    one <- frsr(x, threads = 1, detail = TRUE)
    two <- frsr(x, threads = 2, detail = TRUE)
    expect_identical(one, two)
})

test_that('removed API is absent from the supported namespace', {
    expect_false('frsr_NR' %in% getNamespaceExports('frsrr'))
    expect_false(exists('frsr_NR', envir = asNamespace('frsrr'), inherits = FALSE))
})
