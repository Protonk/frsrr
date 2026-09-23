#' @useDynLib frsrr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL

#' FRSR Bin
#'
#' Select sampled magic constants for the Fast Reciprocal Square Root algorithm over specified bins
#' by minimizing an objective metric.
#'
#' @param x_min Numeric lower bound (> 0). Default is 0.25.
#' @param x_max Numeric upper bound (> x_min). Default is 1.0.
#' @param n_bins Integer. The number of bins to divide the range into. Default is 4.
#' @param NRmax Integer. The maximum number of Newton-Raphson iterations (default: 0).
#' @param objective Character scalar naming the metric to optimize. Defaults to
#'   \code{"max_relative_error"}.
#' @param dependent Character scalar naming the metric reported as dependent on
#'   the chosen objective. Defaults to \code{"avg_relative_error"}.
#'   Both `objective` and `dependent` accept: \code{"max_relative_error"},
#'   \code{"avg_relative_error"}, or \code{"rmse_relative_error"}.
#' @param float_samples Integer. The number of floating-point samples generated per bin.
#' @param magic_samples Integer. The number of magic constant samples to generate per bin.
#' @param magic_min Integer. The minimum magic constant to test (default: 1596980000).
#' @param magic_max Integer. The maximum magic constant to test (default: 1598050000).
#' @param threads Positive integer forwarded to \code{RcppParallel::setThreadOptions()}
#'   before running the parallel candidate search. Defaults to
#'   \code{getOption("frsrr.threads", NA)} which falls back to the package's
#'   compiled default thread count when unset.
#' @param ... Optional sampler arguments forwarded to the float sampler.
#'   Currently supports \code{method} (one of \code{"log_stratified"},
#'   \code{"irrational"}, or \code{"uniform"}). When omitted,
#'   \code{method = "log_stratified"}.
#'
#' @return
#' A data frame with columns:
#'     \item{N_bins}{Total number of bins}
#'     \item{Location}{Bin number (1-indexed)}
#'     \item{Range_Min}{Minimum value of the bin range}
#'     \item{Range_Max}{Maximum value of the bin range}
#'     \item{Magic}{Best tested magic constant on the sampled inputs, as an integer}
#'     \item{Objective}{Metric minimized for that bin (e.g., maximum relative error)}
#'     \item{Dependent}{Secondary metric reported for the winning magic.}
#'
#' @details
#' Each bin is a half-open interval `[Range_Min, Range_Max)`. All candidates
#' in a bin share the same sampled inputs and the float32 arithmetic and
#' double-precision rounded-input error contract of [frsr()], using A = 1.5,
#' B = 0.5 and tol = 0. Equal magic bounds select that constant; reversed
#' bounds are supported. Nonfinite approximations disqualify a candidate;
#' an error is raised if none remain. Exact objective ties select the smallest
#' magic integer. This is the best tested candidate on these samples, not a
#' claim of global optimality.
#'
#' Double sums use fixed sample blocks joined in a fixed order. For the same
#' build, floating-point environment and samples, measurements and selection
#' are identical across thread counts. This is not a cross-platform guarantee.
#' Use `set.seed()` to reproduce sampling.
#'
#' The data frame's `settings` attribute records `objective`, `dependent`,
#' `NRmax`, `method`, `float_samples`, `magic_samples`, `magic_min`, `magic_max`
#' and `threads`. Preserve this attribute when saving results (e.g. with
#' `saveRDS()`); plain CSV does not retain it.
#'
#' @examples
#' set.seed(42)
#' result <- frsr_bin(n_bins = 2, float_samples = 32, magic_samples = 16,
#'                    NRmax = 1, threads = 1)
#' result
#' attr(result, "settings")
#' @name frsr_bin
NULL

#' @rdname frsr_bin
#' @export
frsr_bin <- function(x_min = 0.25, x_max = 1.0,
                     n_bins = 4, NRmax = 0,
                     objective = c("max_relative_error", "avg_relative_error", "rmse_relative_error"),
                     dependent = c("avg_relative_error", "max_relative_error", "rmse_relative_error"),
                     float_samples = 1024, magic_samples = 2048,
                     magic_min = 1596980000L,
                     magic_max = 1598050000L,
                     threads = getOption("frsrr.threads", NA_integer_),
                     ...) {
  objective <- match.arg(objective)
  dependent <- match.arg(dependent)

  x_min <- as.numeric(x_min)[1]
  x_max <- as.numeric(x_max)[1]
  n_bins <- as.integer(n_bins)[1]
  float_samples <- as.integer(float_samples)[1]
  magic_samples <- as.integer(magic_samples)[1]
  magic_min <- as.integer(magic_min)[1]
  magic_max <- as.integer(magic_max)[1]
  NRmax <- as.numeric(NRmax)[1]
  if (!is.finite(NRmax) || NRmax < 0 || NRmax > .Machine$integer.max || NRmax != trunc(NRmax)) {
    stop("`NRmax` must be a non-negative integer", call. = FALSE)
  }
  NRmax <- as.integer(NRmax)
  threads <- frsrr_configure_threads(threads)

  dots <- list(...)
  dot_names <- names(dots)
  if (length(dots) && (is.null(dot_names) || any(!nzchar(dot_names)))) {
    stop("All arguments passed through `...` must be named")
  }
  allowed_sampler_args <- "method"
  unused <- setdiff(dot_names, allowed_sampler_args)
  if (length(unused)) {
    stop("Unused arguments in `...`: ", paste(unused, collapse = ", "))
  }

  sampler_method <- dots$method
  if (is.null(sampler_method)) {
    sampler_method <- "log_stratified"
  } else {
    sampler_method <- as.character(sampler_method)[1]
  }

  if (is.na(float_samples) || float_samples < 1L) stop("`float_samples` must be positive")
  if (is.na(magic_samples) || magic_samples < 1L) stop("`magic_samples` must be positive")
  if (is.na(magic_min) || is.na(magic_max)) stop("Magic bounds cannot be NA")
  sampler_method <- match.arg(sampler_method, .frsrr_sampler_methods)
  settings <- list(objective = objective, dependent = dependent, NRmax = NRmax,
                   method = sampler_method, float_samples = float_samples,
                   magic_samples = magic_samples, magic_min = magic_min,
                   magic_max = magic_max, threads = threads)

  # Argument coercions above intentionally drop vector inputs to a single scalar;
  # the downstream C++ helpers only read the first element, so we keep behavior
  # predictable by trimming here instead of letting implicit recycling occur.
  if (!is.finite(x_min) || !is.finite(x_max)) {
    stop("`x_min` and `x_max` must be finite")
  }
  if (x_min <= 0 || x_max <= 0) {
    stop("`x_min` and `x_max` must both be > 0 to keep log2 well-defined")
  }
  if (x_min >= x_max) {
    stop("`x_min` must be less than `x_max`")
  }
  if (is.na(n_bins) || n_bins < 1L) {
    return(structure(data.frame(
      N_bins = integer(0),
      Location = integer(0),
      Range_Min = numeric(0),
      Range_Max = numeric(0),
      Magic = integer(0),
      Objective = numeric(0),
      Dependent = numeric(0)
    ), settings = settings))
  }

  # Divide [x_min, x_max] into evenly spaced bin boundaries
  bin_edges <- seq(x_min, x_max, length.out = n_bins + 1)

  # Generate results for each bin
  bins <- lapply(seq_len(n_bins), function(i) {
    bin_min <- bin_edges[i]
    bin_max <- bin_edges[i + 1]
    # floats are generated directly from the sampler
    # They are independent of the choice of magic, or how many are sampled.
    floats <- .frsrr_draw_inputs(
      n = float_samples,
      x_min = bin_min,
      x_max = bin_max,
      method = sampler_method
    )
    # Magic constants are explored via simple sampling; drawing with replacement
    # keeps the runtime flat even when the range is narrower than magic_samples.
    magics <- .frsrr_draw_magics(magic_samples, magic_min, magic_max)
    # Call the C++ function to compute optimal magic constant
    result <- .Call('_frsrr_search_optimal_constant',
                    PACKAGE = 'frsrr',
                    floats, magics, NRmax, objective, dependent)

    # Return results as a data frame
    output <- data.frame(
      Location = i,
      Range_Min = bin_min,
      Range_Max = bin_max
    )
    cbind(output, result)
  })

  # Combine results from all bins into a single data frame
  result <- do.call(rbind, bins)
  # Each row inherits the global bin count here so callers can reshape or merge
  # without having to carry around the per-call metadata separately.
  result$N_bins <- rep.int(n_bins, nrow(result))
  result <- result[c(
    "N_bins",
    "Location",
    "Range_Min",
    "Range_Max",
    "Magic",
    "Objective",
    "Dependent"
  )]
  attr(result, "settings") <- settings
  result
}
