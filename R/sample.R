#' Bounded stratified sampling of positive floating point numbers
#'
#' This function is not exported and is intended for internal use only.
#'
#' @param n Number of samples to generate
#' @param low Lower bound exponent
#' @param high Upper bound exponent
#' @return A numeric vector of length n with samples between 2^low and 2^high
#'
#' @details
#'
#' Samples uniformly within exponent ranges to draw unbiased samples
#' between \eqn{2^low} and \eqn{2^high}.
#'
#' @keywords internal
#' @name boundedStratifiedSample
NULL

boundedStratifiedSample <- function(n, low, high) {
    if (n <= 0 || is.integer(n)) {
        stop("'n' must be a positive integer")
    }
    ## Our C++ algorithm works with exponents,
    ## but expressing that for a user is a pain.
    low <- log2(low)
    high <- log2(high)

    tryCatch(
        .Call('_frsrr_boundedStratifiedSample', PACKAGE = 'frsrr', n, low, high),
        error = function(e) {
            stop(e$message, call. = FALSE)
        }
    )
}

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
#' @details The following parameters are passed to \code{\link{frsr.detail}}:
#'
#' \itemize{
#'  \item{"NR"}{Number of iterations for the Newton-Raphson method. Default is \code{1}.}
#'  \item{"A"}{Parameter A for the iteration formula. Default is \code{1.5}.}
#'  \item{"B"}{Parameter B for the iteration formula. Default is \code{0.5}.}
#'  \item{"tol"}{Tolerance level for stopping criterion. Default is \code{0}.}
#'  \item{"keep_params"}{Logical indicating whether to output parameters. Default is \code{FALSE}.}
#' }
#'
#' Floating point values are searched by stratified sampling 
#' which samples uniformly within exponent ranges, as that is how
#' floating point numbers are distributed. The default range is 
#' chosen because the the error of the FISR along \eqn{2^-2} and \eqn{2^0}
#' repeats over the whole range of the function.
#' 
#' The default range for the magic number was determined by experiment.
#' Values within this range are relatively good restoring constants, with
#' numbers much higher or lower requiring more iterations to converge or
#' not converging at all.
#'
#' @return A data frame with the sampled values and optional details from \code{frsr.detail}.
#' 
#' @seealso 
#' 
#' \code{\link{frsr.detail}}
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
#' samples <- sample_frsr(4)
#' print(samples)
#'       input  initial after_one    final       error        diff iters
#' 1 0.5472975 1.427882  1.345168 1.345168 0.004850830 -0.08271444     1
#' 2 0.4076866 1.659807  1.557596 1.557596 0.005469586 -0.10221040     1
#' 3 0.4417882 1.591603  1.496793 1.496793 0.005124446 -0.09481037     1
#' 4 0.7982250 1.176955  1.114741 1.114741 0.004051689 -0.06221318     1
#' }
#'
#' @export
#' @name sample_frsr
NULL

sample_frsr <- function(n, 
                         magic_min = 1596980000L, magic_max = 1598050000L, 
                         x_min = 0.25, x_max = 1.0,
                         NR = 1, A = 1.5, B = 0.5,
                         tol = 0,
                         keep_params = FALSE) {
     # Generate magic numbers using sample() without replacement
     # maybe in the future we can pass the sample option on
     # but practically speaking there's no need for that.
    magic_numbers <- sample(magic_min:magic_max, n, replace = FALSE)
    
    # Generate input values using bounded stratified sampling
    inputs <- boundedStratifiedSample(n, x_min, x_max)
  
    # Call frsr.detail with generated inputs and parameters
    frsr.detail(x = inputs, magic = magic_numbers, NR = NR, tol = tol, A = A, B = B, keep_params = keep_params)
}
