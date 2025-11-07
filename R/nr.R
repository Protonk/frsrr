#' Custom Newton-Raphson Formula Tools
#'
#' Set custom Newton-Raphson formula for the Fast Reciprocal Square Root.
#'
#' @param x A numeric vector of input values.
#' @param formula A custom formula as a function of \code{y} and \code{x}.
#' @param magic An optional magic number for the initial guess.
#' @param NRmax Number of iterations to perform. Default is \code{1}.
#' @param tol Tolerance level for stopping criterion. Default is \code{0}.
#'
#' @return
#' A data frame of \code{length(x)} rows with columns:
#'    \item{input}{The input values}
#'    \item{initial}{Initial approximation from integer operations}
#'    \item{final}{Result from final iteration}
#'    \item{error}{Absolute relative error of final versus standard library}
#'    \item{converged}{Whether the algorithm converged}
#'    \item{conv_rate}{Mean rate of convergence}
#'    \item{iters}{Number of iterations performed}
#'
#' @details
#' While \code{frsr} uses C++ internally, this is difficult to do when 
#' allowing the user to set a custom formula in an R-like idiom. This
#' function will be much slower than \code{frsr}.
#' 
#' When \eqn{(A - B * x * y_n^2)} is the formula, B = A - 1
#' must be true to ensure convergence. If you supply a different iteration step,
#' you must determine what constraint on A and B is necessary. See
#' the references for more information.
#' 
#' @references
#' 
#' Gec, (2024) Answer to "Why does this modified Newton's method fail to converge for N > 1 iterations?" Stack Overflow. \url{https://math.stackexchange.com/a/5007745/1025433}
#'
#' @examples
#' \donttest{
#' x <- c(1, 4, 9, 16)
#' ex_formula <- quote(y * (1.5 - 0.5 * x * y^2))
#' result <- frsr_NR(x, formula = ex_formula)
#' print(result$final)
#' # [1] 0.9990148 0.4995074 0.3337626 0.2497537
#' }
#' 
#' @name frsr_NR
NULL

#' @rdname frsr_NR
#' @export
frsr_NR <- function(x, magic = 0x5f3759df, formula, NRmax = 1, tol = 0) {
  y <- frsr(x, magic, NRmax = 0)
  reference <- 1 / sqrt(x)
  iter <- 0
  initial <- y
  error_initial <- abs(initial - reference) / reference
  ynplusone <- function(x, guess, formula) {
    y <- guess
    # Parse the formula
    f <- as.function(alist(y=, x=, eval(formula)))
    # Newton-Raphson iteration
    f(y, x)
  }

  # Iterate until NRmax is reached or error is below tol
  while (iter < NRmax) {
    new_y <- ynplusone(x, y, formula)
    error <- abs(new_y - reference) / reference
    y <- new_y
    iter <- iter + 1
    # Check the worst relative error so we only stop once every element
    # in the input vector has met the requested tolerance
    if (max(error) < tol && tol > 0) {
      break
    }
  }
  
  data.frame(
    input = x,
    initial = initial,
    final = y,
    error = error,
    converged = ifelse(tol > 0 && error < tol, TRUE, FALSE),
    conv_rate = (error_initial - error) / iter,
    iters = iter
  )
}