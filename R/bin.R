#' @importFrom stats aggregate
NULL

#' FRSR Bin
#'
#' Generate optimal magic constants for the Fast Reciprocal Square Root algorithm over specified bins.
#'
#' @param x_min Numeric. Default is 0.25.
#' @param x_max Numeric. Default is 1.0.
#' @param n_bins Integer. The number of bins to divide the range into. Default is 4.
#' @param float_samples Integer. The number of floating-point samples generated per bin.
#' @param magic_samples Integer. The number of magic constant samples to generate per bin.
#' @param magic_min Integer. The minimum magic constant to test (default: 1596980000).
#' @param magic_max Integer. The maximum magic constant to test (default: 1598050000).
#' @param NRmax Integer. The maximum number of Newton-Raphson iterations (default: 0).
#'
#' @return
#' A data frame with columns:
#'     \item{N_bins}{Total number of bins}
#'     \item{Location}{Bin number (1-indexed)}
#'     \item{Range_Min}{Minimum value of the bin range}
#'     \item{Range_Max}{Maximum value of the bin range}
#'     \item{Magic}{Optimal magic constant as an integer}
#'     \item{Max_Relative_Error}{Maximum relative error for each bin}
#'
#' @details
#' This function divides the range [x_min, x_max] into n_bins bins and generates float_samples
#' floating-point samples within each bin. It then tests magic_samples magic constants within
#' the range [magic_min, magic_max] to find the optimal magic constant for each bin that minimizes
#' the relative error in the Fast Reciprocal Square Root algorithm.
#' 
#' The FRSR is periodic over 0.25 to 1.0 WITH the magic
#' constant of 0x5f3759df. I don't know if this is
#' generally true.
#' 
#' @examples
#' \donttest{
#' # Generate optimal magic constants for the range [0.25, 1.0] divided into 4 bins
#' result <- frsr_bin()
#' print(result)
#' # }
#' @name frsr_bin
NULL

#' @rdname frsr_bin
#' @export
frsr_bin <- function(x_min = 0.25, x_max = 1.0,
                     n_bins = 4, NRmax = 0,
                     float_samples = 1024, magic_samples = 2048,
                     magic_min = 1596980000L,
                     magic_max = 1598050000L) {

  # Calculate bin edges
  bin_edges <- n_bins + 1
  for (i in 0:n_bins) {
    bin_edges[i + 1] <- 2^((log2(x_min) + i * (log2(x_max) - log2(x_min)) / n_bins))
  }

  # Generate stratified samples for each bin and compute results
  bins <- lapply(seq_len(n_bins), function(i) {
    samples <- .Call('_frsrr_boundedStratifiedSample', PACKAGE = 'frsrr',
                     float_samples,
                     log2(bin_edges[i]), log2(bin_edges[i + 1]))
    magic <- sample(magic_min:magic_max, magic_samples, replace = TRUE)
    ## this could become large
    input_df <- expand.grid(samples, magic)
    names(input_df) <- c("float", "magic")
    # frsr is parallel C++, so this is spread out as well as is reasonable.
    input_df$error <- frsr(x = input_df[, "float"], magic = input_df[, "magic"],
                           NRmax = NRmax, tol = 0, detail = TRUE)$error
    input_df <- input_df[, c("magic", "error")]

    # Summarise to calculate the sum of 'error' for each group
    summed_errors <- aggregate(error ~ magic, data = input_df, sum)
    # Find the row with the minimum summed error
    result <- summed_errors[which.min(summed_errors$error), ]
    error_sum <- result[, "error"]
    best_magic <- result[, "magic"]

    data.frame(
      Location = i,
      Range_Min = bin_edges[i],
      Range_Max = bin_edges[i + 1],
      Magic = best_magic,
      # Don't bias smaller float samples
      Sum_Error = error_sum / float_samples
    )
  })
  # Combine all bins into a single data frame
  cbind(data.frame(Bins = n_bins), do.call(rbind, bins))
}
