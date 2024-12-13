#' Sample FRSR
#'
#' Generate samples for the Fast Reciprocal Square Root (FRSR) algorithm.
#'
#' @param n Number of samples to generate.
#' @param magic_min Minimum value for the magic number range. Default is \code{1596980000L}.
#' @param magic_max Maximum value for the magic number range. Default is \code{1598050000L}.
#' @param x_min Minimum value for the input range. Default is \code{0.25}.
#' @param x_max Maximum value for the input range. Default is \code{1.0}.
#'
#' The following parameters are passed to \code{\link{frsr.detail}}:
#' @param NR Number of iterations for the Newton-Raphson method. Default is \code{1}.
#' @param tol Tolerance level for stopping criterion. Default is \code{0}.
#' @param A Parameter A for the iteration formula. Default is \code{1.5}.
#' @param B Parameter B for the iteration formula. Default is \code{0.5}.
#' @param keep_params Logical indicating whether to keep parameters in the output. Default is \code{FALSE}.
#'
#' @return A data frame with the sampled values and details from the FRSR algorithm.
#'
#' @details
#' 
#' Uses R to sample floating point values with Log Stratified Sampling
#' and sample integers with replacement. Allows for flexible data collection
#' for adjustments to the FRSR algorithm.
#' 
#' @references
#' 
#' Walker, A. J. (1974) Fast generation of uniformly distributed pseudorandom numbers with floating-point representation. Electronics Letters, 10, 533-534, \url{https://api.semanticscholar.org/CorpusID:110056594}
#'
#' @examples
#' \donttest{
#' # Generate 10 samples using default parameters
#' samples <- sample_frsr(10)
#' print(samples)
#' }
#'
#' @name sample_frsr
NULL


# Define the logStratifiedSampler function
logStratifiedSampler <- function(min, max, n) {
  exp(runif(n, log(min), log(max)))
}

#' @export
sample_frsr <- function(n,
                        magic_min = 1596980000L, magic_max = 1598050000L,
                        x_min = 0.25, x_max = 1.0,
                        NR = 1, tol = 0, A = 1.5, B = 0.5,
                        keep_params = FALSE) {
  x <- logStratifiedSampler(x_min, x_max, n)
  magic <- sample(magic_min:magic_max, n, replace = TRUE)

  frsr.detail(x, magic = magic, NR = NR, A = A, B = B, tol = tol, keep_params = keep_params)
}
