#' @useDynLib frsrr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Sample FRSR
#'
#' Generate samples for the Fast Reciprocal Square Root (FRSR) algorithm.
#'
#' @param n Number of samples to generate.
#' @param magic_min Minimum value for the magic number range. Default is \code{1596980000L}
#' @param magic_max Maximum value for the magic number range. Default is \code{1598050000L}
#' @param x_min Minimum value for the input range. Default is \code{0.25}
#' @param x_max Maximum value for the input range. Default is \code{1.0}
#' @param ... Additional arguments passed to \code{frsr_compute}.
#'
#' @details 
#' 
#' The default range for the magic number was determined by experiment.
#' Values within this range are relatively good restoring constants, with
#' numbers much higher or lower requiring more iterations to converge or
#' not converging at all.
#'
#' Floating point values are searched by stratified sampling 
#' which samples uniformly within exponent ranges, as that is how
#' floating point numbers are distributed. The default range is 
#' chosen because the the error of the FISR is periodic 
#' from 0.25 to 1.0.
#'
#' @return
#' A data frame with \code{n} rows. When \code{keep_params = FALSE} (the default), the
#' columns match \code{frsr_compute(..., detail = TRUE)}:
#'     \item{input}{Sampled input values}
#'     \item{initial}{Initial approximation from integer operations}
#'     \item{after_one}{Result after one Newton-Raphson iteration}
#'     \item{final}{Result from the last iteration}
#'     \item{error}{Absolute relative error of the final result}
#'     \item{enre}{Exponent-normalized absolute relative error}
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
#' \code{\link{frsr_compute}}
#'
#' @references
#'
#' Walker, A. J. (1974) Fast generation of uniformly distributed pseudorandom numbers with floating-point representation. Electronics Letters, 10, 533-534, \url{https://api.semanticscholar.org/CorpusID:110056594}
#'
#' Pharr, M. (2022) Sampling in Floating Point (2/3): 1D Intervals. Matt Pharr's Blog, \url{https://pharr.org/matt/blog/2022/03/14/sampling-float-intervals}
#'
#' @examples 
#' 
#' #' \donttest{
#' # Generate 4 samples using default parameters
#' samples <- frsr_sample(4)
#' print(samples)
#' #       input  initial after_one    final        error         enre         diff iters
#' # 1 0.8100569 1.154508  1.108492 1.108492 0.0023223383 0.0023223383 -0.046015978     1
#' # 2 0.4073950 1.573407  1.566680 1.566680 0.0000273157 0.0000273157 -0.006727338     1
#' # 3 0.6304980 1.359966  1.247014 1.247014 0.0098225391 0.0098225391 -0.112952113     1
#' # 4 0.4316622 1.607666  1.514687 1.514687 0.0048355814 0.0048355814 -0.092979550     1
#' #}
#'
#' @export
#' @name frsr_sample
NULL

#' @rdname frsr_sample
#' @export
frsr_sample <- function(n, 
                        magic_min = 1596980000L, magic_max = 1598050000L, 
                        x_min = 0.25, x_max = 1.0, ...) {
    # Determine magic numbers based on whether magic_min or magic_max is NULL
    magic_numbers <- if (is.null(magic_min)) {
        rep(magic_max, n)  # Use magic_max if magic_min is NULL
    } else if (is.null(magic_max)) {
        rep(magic_min, n)  # Use magic_min if magic_max is NULL
    } else {
        sample(magic_min:magic_max, n, replace = TRUE)
    }
    # Determine inputs based on whether x_min or x_max is NULL
    inputs <- if (is.null(x_min)) {
        rep(x_max, n)  # Use x_max if x_min is NULL
    } else if (is.null(x_max)) {
        rep(x_min, n)  # Use x_min if x_max is NULL
    } else {
    # Forward log2 bounds because the C++ helper samples uniformly across exponent strata
    .Call('_frsrr_boundedStratifiedSample',
          PACKAGE = 'frsrr',
          n, log2(x_min), log2(x_max))
    }
    # Call frsr_compute with generated inputs and parameters
    frsr_compute(x = inputs, magic = magic_numbers, detail = TRUE, ...)
}
