#' @useDynLib frsrr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL

#' FRSR Bin
#'
#' Generate optimal magic constants for the Fast Reciprocal Square Root algorithm over specified bins
#' by minimizing a user-selected objective metric.
#'
#' @param x_min Numeric. Default is 0.25.
#' @param x_max Numeric. Default is 1.0.
#' @param n_bins Integer. The number of bins to divide the range into. Default is 4.
#' @param float_samples Integer. The number of floating-point samples generated per bin.
#' @param magic_samples Integer. The number of magic constant samples to generate per bin.
#' @param magic_min Integer. The minimum magic constant to test (default: 1596980000).
#' @param magic_max Integer. The maximum magic constant to test (default: 1598050000).
#' @param weighted Logical; if \code{TRUE}, weight the float sampler by the
#'   number of representable significands in each exponent stratum. Default is
#'   \code{FALSE}.
#' @param NRmax Integer. The maximum number of Newton-Raphson iterations (default: 0).
#' @param objective_metric Character scalar naming the metric used to choose the optimal magic
#'   constant. Supported values are \code{"max"}, \code{"mean"}, \code{"sum"}, \code{"rms"}, and
#'   \code{"variance"}.
#' @param metrics_to_report Character vector naming additional metrics to summarise for the selected
#'   magic constant. Each entry must be one of the supported metric names. Defaults to reporting the
#'   mean.
#'
#' @return
#' A data frame with columns:
#'     \item{N_bins}{Total number of bins}
#'     \item{Location}{Bin number (1-indexed)}
#'     \item{Range_Min}{Minimum value of the bin range}
#'     \item{Range_Max}{Maximum value of the bin range}
#'     \item{Magic}{Optimal magic constant as an integer}
#'     \item{Objective_<Metric>}{Value of the selected objective metric for the optimal magic}
#'     \item{Metric_<Metric>}{Summary metrics requested via \code{metrics_to_report}. The suffix
#'       matches the metric name (for example, \code{Metric_Mean}).}
#'
#' @details
#' This function divides the range [x_min, x_max] into n_bins bins and generates float_samples
#' floating-point samples within each bin. It then tests magic_samples magic constants within
#' the range [magic_min, magic_max] to find the optimal magic constant for each bin that minimizes
#' the objective metric per bin under the chosen configuration. The achievable measurements and optimal
#' constants will change with bin size and number of bins, as well as the integer and float samples.
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
#' #   Location Range_Min Range_Max      Magic Objective_Max Metric_Mean N_bins
#' # 1        1 0.2500000 0.3535534 1597413411   0.000454920 0.000170461      4
#' # 2        2 0.3535534 0.5000000 1597115686   0.000160304 0.000077194      4
#' # 3        3 0.5000000 0.7071068 1597150657   0.000135623 0.000060779      4
#' # 4        4 0.7071068 1.0000000 1597488127   0.000516480 0.000209730      4
#' # }
#' @name frsr_bin
NULL

#' @rdname frsr_bin
#' @export
frsr_bin <- function(x_min = 0.25, x_max = 1.0,
                     n_bins = 4, NRmax = 0,
                     float_samples = 1024, magic_samples = 2048,
                     magic_min = 1596980000L,
                     magic_max = 1598050000L,
                     weighted = FALSE,
                     objective_metric = "max",
                     metrics_to_report = c("mean")) {
  if (!is.numeric(x_min) || length(x_min) != 1L || !is.finite(x_min)) {
    stop("`x_min` must be a finite numeric scalar", call. = FALSE)
  }
  if (!is.numeric(x_max) || length(x_max) != 1L || !is.finite(x_max)) {
    stop("`x_max` must be a finite numeric scalar", call. = FALSE)
  }
  if (x_min >= x_max) {
    stop("`x_min` must be less than `x_max`", call. = FALSE)
  }

  if (!is.numeric(n_bins) || length(n_bins) != 1L || is.na(n_bins)) {
    stop("`n_bins` must be at least 1", call. = FALSE)
  }
  n_bins <- as.integer(n_bins)
  if (is.na(n_bins) || n_bins < 1L) {
    stop("`n_bins` must be at least 1", call. = FALSE)
  }

  if (!is.numeric(float_samples) || length(float_samples) != 1L || is.na(float_samples)) {
    stop("`float_samples` must be at least 1", call. = FALSE)
  }
  float_samples <- as.integer(float_samples)
  if (is.na(float_samples) || float_samples < 1L) {
    stop("`float_samples` must be at least 1", call. = FALSE)
  }

  if (!is.numeric(magic_samples) || length(magic_samples) != 1L || is.na(magic_samples)) {
    stop("`magic_samples` must be at least 1", call. = FALSE)
  }
  magic_samples <- as.integer(magic_samples)
  if (is.na(magic_samples) || magic_samples < 1L) {
    stop("`magic_samples` must be at least 1", call. = FALSE)
  }

  if (!is.numeric(magic_min) || length(magic_min) != 1L || !is.finite(magic_min)) {
    stop("`magic_min` must be a finite integer scalar", call. = FALSE)
  }
  if (!is.numeric(magic_max) || length(magic_max) != 1L || !is.finite(magic_max)) {
    stop("`magic_max` must be a finite integer scalar", call. = FALSE)
  }
  magic_min <- as.integer(magic_min)
  magic_max <- as.integer(magic_max)
  if (is.na(magic_min) || is.na(magic_max)) {
    stop("`magic_min` and `magic_max` must be representable as 32-bit integers", call. = FALSE)
  }
  if (magic_min > magic_max) {
    stop("`magic_min` must be less than or equal to `magic_max`", call. = FALSE)
  }

  if (!is.numeric(NRmax) || length(NRmax) != 1L || !is.finite(NRmax)) {
    stop("`NRmax` must be a finite numeric scalar", call. = FALSE)
  }
  NRmax <- as.integer(NRmax)
  if (is.na(NRmax) || NRmax < 0L) {
    stop("`NRmax` must be a non-negative integer", call. = FALSE)
  }

  if (!is.logical(weighted) || length(weighted) != 1L || is.na(weighted)) {
    stop("`weighted` must be a non-missing logical scalar", call. = FALSE)
  }

  if ((!is.character(objective_metric) && !is.factor(objective_metric)) ||
      length(objective_metric) != 1L) {
    stop("`objective_metric` must be a non-missing character scalar", call. = FALSE)
  }
  objective_metric <- as.character(objective_metric)
  if (is.na(objective_metric)) {
    stop("`objective_metric` must be a non-missing character scalar", call. = FALSE)
  }
  if (!nzchar(objective_metric)) {
    stop("`objective_metric` must not be empty", call. = FALSE)
  }

  if (!is.character(metrics_to_report) && !is.factor(metrics_to_report)) {
    stop("`metrics_to_report` must be a character vector", call. = FALSE)
  }
  metrics_to_report <- as.character(metrics_to_report)
  if (length(metrics_to_report) < 1L) {
    stop("`metrics_to_report` must contain at least one metric name", call. = FALSE)
  }
  if (any(is.na(metrics_to_report))) {
    stop("`metrics_to_report` must not contain missing values", call. = FALSE)
  }
  if (any(!nzchar(metrics_to_report))) {
    stop("`metrics_to_report` must not contain empty strings", call. = FALSE)
  }
  metrics_to_report <- unique(metrics_to_report)

  # Divide [x_min, x_max] into evenly spaced bin boundaries
  bin_edges <- seq(x_min, x_max, length.out = n_bins + 1)
  
  # Generate results for each bin
  bins <- lapply(seq_len(n_bins), function(i) {
    bin_min <- bin_edges[i]
    bin_max <- bin_edges[i + 1]
    # Pass exponent bounds because the sampler stratifies floats by log2 exponent range
    floats <- .Call('_frsrr_bounded_stratified_sample',
                    PACKAGE = 'frsrr',
                    float_samples, log2(bin_min), log2(bin_max), weighted)
    magics <- sample(magic_min:magic_max,
                     size = magic_samples,
                     replace = TRUE)
    # Call the C++ function to compute optimal magic constant
    result <- .Call('_frsrr_search_optimal_constant',
                    PACKAGE = 'frsrr',
                    floats, magics, NRmax,
                    objective_metric,
                    metrics_to_report)
    
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
  result$N_bins <- rep.int(n_bins, nrow(result))
  ordered_cols <- c(
    "N_bins",
    "Location",
    "Range_Min",
    "Range_Max",
    names(result)[!(names(result) %in% c("N_bins", "Location", "Range_Min", "Range_Max"))]
  )
  result[ordered_cols]
}
