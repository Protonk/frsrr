#' @useDynLib frsrr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL

#' FRSR Bin
#'
#' Generate optimal magic constants for the Fast Reciprocal Square Root algorithm over specified bins
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
#' @param weighted Logical; if \code{TRUE}, weight the float sampler by the
#'   number of representable significands in each exponent stratum. Default is
#'   \code{FALSE}.
#' @param threads Positive integer forwarded to \code{RcppParallel::setThreadOptions()}
#'   before running the parallel candidate search. Defaults to
#'   \code{getOption("frsrr.threads", NA)} which falls back to the package's
#'   compiled default thread count when unset.
#'
#' @return
#' A data frame with columns:
#'     \item{N_bins}{Total number of bins}
#'     \item{Location}{Bin number (1-indexed)}
#'     \item{Range_Min}{Minimum value of the bin range}
#'     \item{Range_Max}{Maximum value of the bin range}
#'     \item{Magic}{Optimal magic constant as an integer}
#'     \item{Objective}{Metric minimized for that bin (e.g., maximum relative error)}
#'     \item{Dependent}{Secondary metric reported for the winning magic.}
#'
#' @details
#' This function divides the range [x_min, x_max] into n_bins bins and generates float_samples
#' floating-point samples within each bin. It then tests magic_samples magic constants within
#' the range [magic_min, magic_max] to find the optimal magic constant for each bin that minimizes
#' the maximum relative error per bin. The achievable error and optimal constants will change with
#' bin size and number of bins, as well as the integer and float samples.
#'
#' The FRSR is periodic over 0.25 to 1.0 with a magic constant of 0x5f3759df. I don't know if
#' this is generally true. Interestingly, the linear approximation to the logarithm used is
#' periodic over integral powers of two, and so loops twice before the FISR does once.
#'
#' @examples
#' \donttest{
#' # Generate optimal magic constants for the range [0.25, 1.0] divided into 4 bins
#' result <- frsr_bin()
#' print(result)
#' #   N_bins Location Range_Min Range_Max      Magic   Objective   Dependent
#' # 1      4        1 0.2500000 0.3535534 1597413411 0.0004549199 0.000134129
#' # 2      4        2 0.3535534 0.5000000 1597115686 0.0000617551 0.000033011
#' # 3      4        3 0.5000000 0.7071068 1597150657 0.0000348620 0.000021047
#' # 4      4        4 0.7071068 1.0000000 1597488127 0.0008389181 0.000244221
#' # }
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
                     weighted = FALSE,
                     threads = getOption("frsrr.threads", NA_integer_)) {
  objective <- match.arg(objective)
  dependent <- match.arg(dependent)

  x_min <- as.numeric(x_min)[1]
  x_max <- as.numeric(x_max)[1]
  n_bins <- as.integer(n_bins)[1]
  float_samples <- as.integer(float_samples)[1]
  magic_samples <- as.integer(magic_samples)[1]
  magic_min <- as.integer(magic_min)[1]
  magic_max <- as.integer(magic_max)[1]
  NRmax <- as.integer(NRmax)[1]
  weighted <- isTRUE(weighted)
  threads <- frsrr_configure_threads(threads)

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
    return(data.frame(
      N_bins = integer(0),
      Location = integer(0),
      Range_Min = numeric(0),
      Range_Max = numeric(0),
      Magic = integer(0),
      Objective = numeric(0),
      Dependent = numeric(0)
    ))
  }

  # Divide [x_min, x_max] into evenly spaced bin boundaries
  bin_edges <- seq(x_min, x_max, length.out = n_bins + 1)

  # Generate results for each bin
  bins <- lapply(seq_len(n_bins), function(i) {
    bin_min <- bin_edges[i]
    bin_max <- bin_edges[i + 1]
    floats <- .Call('_frsrr_sample_inputs',
                    PACKAGE = 'frsrr',
                    float_samples, bin_min, bin_max, weighted, "log_stratified")
    # Magic constants are explored via simple sampling; drawing with replacement
    # keeps the runtime flat even when the range is narrower than magic_samples.
    magics <- sample(magic_min:magic_max,
                     size = magic_samples,
                     replace = TRUE)
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
  result[c(
    "N_bins",
    "Location",
    "Range_Min",
    "Range_Max",
    "Magic",
    "Objective",
    "Dependent"
  )]
}
