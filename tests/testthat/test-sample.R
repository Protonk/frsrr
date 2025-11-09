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
})

describe("bounded_stratified_sample", {
    it("handles narrow exponent ranges", {
        low <- log2(0.75)
        high <- log2(1)

        samples <- .Call(
            "_frsrr_bounded_stratified_sample",
            PACKAGE = "frsrr",
            32L,
            low,
            high,
            FALSE
        )

        expect_length(samples, 32)
        expect_true(all(is.finite(samples)))
        expect_true(all(samples >= 0.75))
        expect_true(all(samples < 1))
    })

    it("tolerates zero bit draws", {
        samples <- .Call(
            "_frsrr_bounded_stratified_sample",
            PACKAGE = "frsrr",
            1024L,
            log2(0.25),
            log2(0.5),
            FALSE
        )

        expect_length(samples, 1024)
        expect_true(all(is.finite(samples)))
    })

    it("supports weighted sampling", {
        samples <- .Call(
            "_frsrr_bounded_stratified_sample",
            PACKAGE = "frsrr",
            64L,
            log2(0.5),
            log2(2),
            TRUE
        )

        expect_length(samples, 64)
        expect_true(all(is.finite(samples)))
        expect_true(all(samples >= 0.5))
        expect_true(all(samples < 2))
    })
})
