describe("frsr_sample", {
    it("returns correct number of samples", {
        result <- frsr_sample(4)
        expect_equal(nrow(result), 4)

        weighted_result <- frsr_sample(4, weighted = TRUE)
        expect_equal(nrow(weighted_result), 4)
    })

    it("forwards Newton arguments through to frsr", {
        result <- frsr_sample(5, NRmax = 2)
        expect_identical(unique(result$iters), 2L)
    })

    it("returns parameters when keep_params is TRUE", {
        base_cols <- c("input", "initial", "after_one", "final", "error", "enre", "diff", "iters")
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
        expected_cols <- c("input", "initial", "after_one", "final", "error", "enre", "diff", "iters")

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

    it("rejects weighted samplers for non-log methods", {
        expect_error(
            frsr_sample(2, method = "uniform", weighted = TRUE),
            "`weighted` can only be TRUE"
        )
    })

})

describe("sample_inputs", {
    sample_call <- function(n, x_min, x_max, weighted = FALSE, method = "log_stratified") {
        .Call(
            "_frsrr_sample_inputs",
            PACKAGE = "frsrr",
            as.integer(n),
            x_min,
            x_max,
            weighted,
            method
        )
    }

    it("validates inputs", {
        expect_error(sample_call(-1, 0.25, 1), "`n` must be non-negative")
        expect_error(sample_call(4, 1, 1), "`x_min` must be less than `x_max`")
        expect_error(sample_call(4, 0, 1, method = "log_stratified"), "must be > 0")
        expect_error(sample_call(4, 1, 2, method = "unknown"), "Unknown sampler method")
        expect_error(sample_call(4, 1, 2, weighted = TRUE, method = "uniform"), "only supported")
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
        specs <- list(
            list(method = "log_stratified", weighted = FALSE),
            list(method = "log_stratified", weighted = TRUE),
            list(method = "irrational", weighted = FALSE),
            list(method = "uniform", weighted = FALSE)
        )
        for (spec in specs) {
            set.seed(42)
            first <- sample_call(32, 0.5, 2, weighted = spec$weighted, method = spec$method)
            set.seed(42)
            second <- sample_call(32, 0.5, 2, weighted = spec$weighted, method = spec$method)
            expect_identical(first, second)
        }
    })
})
