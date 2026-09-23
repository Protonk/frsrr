#' @useDynLib frsrr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# Sampler identifiers understood by the C++ draw helper.
.frsrr_sampler_methods <- c("log_stratified", "irrational", "uniform")

#' Internal wrapper around `_frsrr_sample_inputs`.
#'
#' Draws `n` floating-point samples inside `[x_min, x_max)` using one of the
#' supported sampling strategies. Keeping this helper in R makes argument
#' validation explicit before control passes to the parallel C++ workers.
#'
#' @param n Non-negative integer number of draws.
#' @param x_min,x_max Lower-inclusive, upper-exclusive bounds (must satisfy `0 < x_min < x_max`).
#' @param method Character scalar naming the sampler (see `.frsrr_sampler_methods`).
#'
#' @return
#' A numeric vector of length `n` whose values lie within `[x_min, x_max)`.
#'
#' @keywords internal
#' @noRd
.frsrr_draw_inputs <- function(n, x_min, x_max, method = "log_stratified") {
    method <- match.arg(method, .frsrr_sampler_methods)

    .Call(
        '_frsrr_sample_inputs',
        PACKAGE = 'frsrr',
        n, x_min, x_max, method
    )
}

#' Sample FRSR
#'
#' Generate samples for the Fast Reciprocal Square Root (FRSR) algorithm.
#'
#' @param n Number of samples to generate.
#' @param magic_min Minimum value for the magic number range. Default is \code{1596980000L}
#' @param magic_max Maximum value for the magic number range. Default is \code{1598050000L}
#' @param x_min Minimum value for the input range (must be > 0). Default is \code{0.25}
#' @param x_max Maximum value for the input range (must exceed \code{x_min}). Default is \code{1.0}
#' @param method Character scalar selecting the sampler. Options are
#'   \code{"irrational"}, \code{"uniform"}, or
#'   \code{"log_stratified"} (the legacy default).
#' @param ... Additional arguments passed to \code{frsr}.
#'
#' @details
#'
#' The default range for the magic number was determined by experiment.
#' Values within this range are relatively good restoring constants, with
#' numbers much higher or lower requiring more iterations to converge or
#' not converging at all.
#'
#' Equal magic bounds always use that constant; reversed bounds are supported.
#' Log-stratified sampling returns positive normal float32 values in the exact
#' half-open interval and errors if none exist; its bounds must lie in
#' `[2^-126, 2^128]`. Other samplers return R doubles in the half-open interval;
#' their subsequent float32 conversion follows [frsr()]. A NULL input bound
#' selects the other bound as a fixed input, rather than sampling an interval.
#'
#' Three sampler modes explore different coverage patterns over
#' \code{[x_min, x_max)}:
#' \itemize{
#'   \item{\strong{Log-stratified}:} Sample floats uniformly across exponent strata
#'     This method is the default and covers the FP subset of the reals well.
#'   \item{\strong{Irrational rotation}:} Step through \code{[0, 1)} via the golden
#'     ratio increment and rescale. This method provides low-discrepancy coverage
#'     of the unit interval.
#'   \item{\strong{Uniform}:} Draw from the standard R uniform sampler and rescale.
#'    Takes longer to smooth out than irrational rotation.
#' }
#'
#' @return
#' A data frame with \code{n} rows. When \code{keep_params = FALSE} (the default), the
#' columns match \code{frsr(..., detail = TRUE)}:
#'     \item{input}{Sampled input values}
#'     \item{initial}{Initial approximation from integer operations}
#'     \item{after_one}{Result after one Newton-Raphson iteration}
#'     \item{final}{Result from the last iteration}
#'     \item{error}{Absolute relative error of the final result}
#'     \item{diff}{Difference between the final and penultimate approximations}
#'     \item{iters}{Number of iterations performed}
#'
#' If \code{keep_params = TRUE}, the data frame will also include columns:
#'     \item{magic}{Magic constant(s) used for each sample}
#'     \item{NRmax}{Maximum number of Newton-Raphson iterations}
#'     \item{A}{Newton-Raphson parameter A}
#'     \item{B}{Newton-Raphson parameter B}
#'     \item{tol}{Specified tolerance}
#'
#' @seealso
#'
#' \code{\link{frsr}}
#'
#' @references
#'
#' Walker, A. J. (1974) Fast generation of uniformly distributed pseudorandom numbers with floating-point representation. Electronics Letters, 10, 533-534, \url{https://api.semanticscholar.org/CorpusID:110056594}
#'
#' Pharr, M. (2022) Sampling in Floating Point (2/3): 1D Intervals. Matt Pharr's Blog, \url{https://pharr.org/matt/blog/2022/03/14/sampling-float-intervals}
#'
#' @examples
#' set.seed(42)
#' frsr_sample(4, magic_min = 0x5f3759df, magic_max = 0x5f3759df,
#'             NRmax = 1, tol = 0, threads = 1, keep_params = TRUE)
#'
#' @export
#' @name frsr_sample
NULL

#' @rdname frsr_sample
#' @export
frsr_sample <- function(n,
                        magic_min = 1596980000L, magic_max = 1598050000L,
                        x_min = 0.25, x_max = 1.0,
                        method = c("log_stratified", "irrational", "uniform"),
                        ...) {
    method <- match.arg(method)
    n <- as.integer(n)[1]
    if (is.na(n) || n < 0L) {
        stop("`n` must be a non-negative scalar")
    }

    # bounds check hopefully adds to readability
    normalize_bound <- function(value, label) {
        if (is.null(value)) {
            return(NULL)
        }
        scalar <- as.numeric(value)[1]
        if (!is.finite(scalar)) {
            stop("`", label, "` must be finite when provided")
        }
        if (scalar <= 0) {
            stop("`", label, "` must be greater than 0")
        }
        scalar
    }

    x_min <- normalize_bound(x_min, "x_min")
    x_max <- normalize_bound(x_max, "x_max")
    if (!is.null(x_min) && !is.null(x_max) && x_min >= x_max) {
        stop("`x_min` must be less than `x_max` when both are supplied")
    }

    if (is.null(magic_min) && is.null(magic_max)) {
        stop("At least one magic bound is required")
    }
    if (is.null(x_min) && is.null(x_max)) {
        stop("At least one input bound is required")
    }
    # Determine magic numbers based on whether magic_min or magic_max is NULL
    magic_numbers <- if (is.null(magic_min)) {
        rep(magic_max, n)  # Use magic_max if magic_min is NULL
    } else if (is.null(magic_max)) {
        rep(magic_min, n)  # Use magic_min if magic_max is NULL
    } else {
        # Sample with replacement so we explore the full range even when n
        # exceeds the integer interval size.
        .frsrr_draw_magics(n, magic_min, magic_max)
    }
    # Determine inputs based on whether x_min or x_max is NULL
    inputs <- if (is.null(x_min)) {
        rep(x_max, n)  # Use x_max if x_min is NULL
    } else if (is.null(x_max)) {
        rep(x_min, n)  # Use x_min if x_max is NULL
    } else {
        .frsrr_draw_inputs(
            n = n,
            x_min = x_min,
            x_max = x_max,
            method = method
        )
    }
    # Call frsr with generated inputs and parameters
    # detail = TRUE keeps diagnostics users typically want
    frsr(x = inputs, magic = magic_numbers, detail = TRUE, ...)
}

# Sample offsets, avoiding sample()'s special treatment of a singleton integer.
# Signed R integers (except NA) remain valid exploratory bit patterns.
.frsrr_draw_magics <- function(n, lower, upper) {
    bounds <- c(lower, upper)
    if (length(bounds) != 2L || any(!is.finite(bounds) |
        abs(bounds) > .Machine$integer.max | bounds != trunc(bounds))) {
        stop("Magic bounds must be non-missing R integers", call. = FALSE)
    }
    lower <- as.double(lower)
    upper <- as.double(upper)
    step <- if (upper >= lower) 1 else -1
    as.integer(lower + step * (sample.int(abs(upper - lower) + 1, n, replace = TRUE) - 1))
}
