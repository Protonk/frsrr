set.seed(42)

describe("frsr_sample", {
    it("returns correct number of samples", {
        result <- frsr_sample(4)
        expect_equal(nrow(result), 4)
    })

    it("validates scalar arguments up front", {
        expect_error(
            frsr_sample(-1),
            "`n` must be a non-negative scalar",
            fixed = TRUE
        )
        expect_error(
            frsr_sample(NA_integer_),
            "`n` must be a non-negative scalar",
            fixed = TRUE
        )
        expect_error(
            frsr_sample(1, x_min = 0),
            "`x_min` must be greater than 0",
            fixed = TRUE
        )
        expect_error(
            frsr_sample(1, x_max = 0),
            "`x_max` must be greater than 0",
            fixed = TRUE
        )
        expect_error(
            frsr_sample(1, x_min = Inf),
            "`x_min` must be finite when provided",
            fixed = TRUE
        )
        expect_error(
            frsr_sample(1, x_min = 1, x_max = 0.5),
            "`x_min` must be less than `x_max` when both are supplied",
            fixed = TRUE
        )
    })

    it("forwards Newton arguments through to frsr", {
        result <- frsr_sample(5, NRmax = 2)
        expect_identical(unique(result$iters), 2L)
    })

    it("returns parameters when keep_params is TRUE", {
        base_cols <- c("input", "initial", "after_one", "final", "error", "diff", "iters")
        param_cols <- c("magic", "NRmax", "A", "B", "tol")

        result <- frsr_sample(4, keep_params = TRUE)
        expect_identical(names(result), c(base_cols, param_cols))
    })

    it("handles NULL magic_min and magic_max correctly", {
        result <- frsr_sample(4, magic_min = NULL, keep_params = TRUE)
        expect_equal(length(unique(result$magic)), 1)

        result <- frsr_sample(4, magic_max = NULL, keep_params = TRUE)
        expect_equal(length(unique(result$magic)), 1)
    })

    it("handles NULL x_min and x_max correctly", {
        result <- frsr_sample(4, x_min = NULL)
        expect_equal(length(unique(result$input)), 1)

        result <- frsr_sample(4, x_max = NULL)
        expect_equal(length(unique(result$input)), 1)
    })

    it("returns documented columns", {
        result <- frsr_sample(4)
        expected_cols <- c("input", "initial", "after_one", "final", "error", "diff", "iters")

        expect_identical(names(result), expected_cols)
        expect_true(all(vapply(result, is.numeric, logical(1))))
    })

    it("is reproducible when seed is set", {
        set.seed(123)
        first <- frsr_sample(6, keep_params = TRUE)
        set.seed(123)
        second <- frsr_sample(6, keep_params = TRUE)
        expect_identical(first, second)
    })

})

describe("sample_inputs", {
    sample_call <- function(n, x_min, x_max, method = "log_stratified") {
        .Call(
            "_frsrr_sample_inputs",
            PACKAGE = "frsrr",
            as.integer(n),
            x_min,
            x_max,
            method
        )
    }

    it("validates inputs", {
        expect_error(sample_call(-1, 0.25, 1), "`n` must be non-negative")
        expect_error(sample_call(4, 1, 1), "`x_min` must be less than `x_max`")
        expect_error(sample_call(4, 0, 1, method = "log_stratified"), "must be > 0")
        expect_error(sample_call(4, 1, 2, method = "unknown"), "Unknown sampler method")
    })

    it("keeps samples within range", {
        methods <- c("log_stratified", "irrational", "uniform")
        for (method in methods) {
            draws <- sample_call(64, 0.25, 1, method = method)
            expect_true(all(draws >= 0.25))
            expect_true(all(draws <= 1))
        }
    })

    it("is reproducible for stochastic samplers", {
        specs <- c("log_stratified", "irrational", "uniform")
        for (method in specs) {
            set.seed(42)
            first <- sample_call(32, 0.5, 2, method = method)
            set.seed(42)
            second <- sample_call(32, 0.5, 2, method = method)
            expect_identical(first, second)
        }
    })
})

test_that('equal and reversed magic bounds use only the requested constants', {
    set.seed(104)
    magic <- as.integer(0x5f3759df)
    fit <- frsr_sample(16, magic_min = magic, magic_max = magic,
                       NRmax = 0, keep_params = TRUE, threads = 1)
    expect_identical(fit$magic, rep(magic, 16))
    fit <- frsr_sample(64, magic_min = magic + 1L, magic_max = magic,
                       NRmax = 0, keep_params = TRUE, threads = 1)
    expect_setequal(fit$magic, c(magic, magic + 1L))
    expect_type(frsrr:::.frsrr_draw_magics(8, -.Machine$integer.max,
                                          .Machine$integer.max), 'integer')
})

test_that('log-stratified sampling honors original half-open float boundaries', {
    set.seed(15)
    draw <- function(low, high) frsrr:::.frsrr_draw_inputs(64L, low, high)
    # Large exponents also test intervals whose log2 endpoints become identical.
    for (e in c(-126, -1, 0, 30, 100, 127)) {
        low <- 2^e
        next_float <- low * (1 + 2^-23)
        expect_identical(draw(low, low * (1 + 2^-24)), rep(low, 64))
        expect_identical(draw(low, next_float), rep(low, 64))
        expect_identical(draw(low * (1 + 2^-24), low * (1 + 2^-22)),
                         rep(next_float, 64))
        expect_error(draw(low * (1 + 2^-52), low * (1 + 2^-24)), 'No admissible')
        expect_error(draw(low * (1 + 2^-24), next_float), 'No admissible')
        values <- draw(low, low * (1 + 2^-22))
        expect_setequal(values, c(low, next_float))
        expect_true(all(values >= low & values < low * (1 + 2^-22)))
    }
    expect_identical(draw(2 - 2^-23, 2), rep(2 - 2^-23, 64))
    expect_setequal(draw(2 - 2^-23, 2 + 2^-22), c(2 - 2^-23, 2))
    largest <- (2 - 2^-23) * 2^127
    expect_identical(draw(largest, 2^128), rep(largest, 64))
    expect_error(draw(2^-127, 1), 'bounds')
    # Exactly one R double in this half-open interval, for the other samplers.
    for (method in c('irrational', 'uniform')) {
        values <- frsrr:::.frsrr_draw_inputs(64L, 1, 1 + 2^-52, method)
        expect_identical(values, rep(1, 64))
    }
})
