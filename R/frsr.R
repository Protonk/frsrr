#' @useDynLib frsrr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL

#' Fast Reciprocal Square Root (FRSR)
#'
#' A parameterized Fast Reciprocal Square Root algorithm written in C++.
#'
#' @param x A numeric vector of input values.
#' @param magic Integer restoring constant. Default is 0x5f3759df.
#' @param NRmax Integer specifying the maximum Newton-Raphson iterations. Default is 1.
#' @param A Newton-Raphson parameter where \eqn{(A - B * x * y_n^2)}. Default is \code{1.5}.
#' @param B Newton-Raphson parameter. Default is \code{0.5}.
#' @param tol The absolute relative error at which to stop early. Default is 0 (no early stopping).
#' @param detail Logical. If \code{TRUE}, a data frame with detailed results is returned. Default is \code{FALSE}.
#' @param keep_params Logical. If \code{TRUE}, generation parameters are included in detailed output. Default is \code{FALSE}.
#'
#' @return
#' \code{frsr_compute} returns a numeric vector of \code{length(x)}.
#'
#' If \code{detail = TRUE}, returns a data frame of \code{length(x)} rows with columns:
#'     \item{input}{The input values}
#'     \item{initial}{Initial approximation from integer operations}
#'     \item{after_one}{Result after one iteration of Newton-Raphson}
#'     \item{final}{Result from final iteration}
#'     \item{error}{Absolute relative error of final versus standard library}
#'     \item{enre}{Exponent-normalized absolute relative error}
#'     \item{diff}{Difference between final and penultimate approximations}
#'     \item{iters}{Number of iterations performed}
#'
#' If \code{keep_params = TRUE}, the data frame will also include columns:
#'    \item{magic}{The magic constant(s) used}
#'    \item{NRmax}{Maximum number of Newton-Raphson iterations}
#'    \item{A}{Newton-Raphson parameter A}
#'    \item{B}{Newton-Raphson parameter B}
#'    \item{tol}{Specified tolerance}
#'
#' @details
#'
#' The function supplies a Fast Reciprocal Square Root algorithm, which provides
#' an approximation of 1/sqrt(x). The user can specify their own parameters. The
#' default values are set to those used by the famous "fast inverse square
#' root" in Quake III Arena.
#'
#' The algorithm exploits the fact that the integer representation of a
#' floating-point number offers a piecewise linear approximation to the
#' logarithm function. Right-shifting the integer bits of a float is
#' equivalent to dividing the logarithm of the number by two.
#' By subtracting this from a carefully chosen constant, an approximation
#' of \eqn{-1/2 * log2(x)} can be obtained. Treating that result as a float
#' again by using integer bits in memory gives a good guess of
#' \eqn{exp(-1/2 * log2(x))}, which is \eqn{1/sqrt(x)}.
#'
#' The "magic" constant principally serves to restore the exponent bits lost
#' when the input float is right shifted. A restoring constant which does only
#' that is `0x5F400000`, given by Blinn 1997. Values of magic from roughly
#' `0x5f2ffb20` to `0x5f404ed0` will give acceptable levels of error.
#'
#' The Newton-Raphson step \eqn{y_{n+1} = y_n * (1.5 - 0.5 * x * y_n^2)}
#' is performed repeatedly until the specified maximum NRmax is reached. The
#' default is one. Grossly different values of magic from the default may
#' require many iterations to approach the correct output.
#'
#' Parameters in the Newton-Raphson step, \eqn{(A - B * x * y_n^2)} need
#' not be fixed at 1.5 and 0.5 and can be set by the user. Note that
#' if B =/= A - 1, the approximation may fail to converge.
#'
#' @references
#' J. F. Blinn, (July-Aug. 1997) "Floating-point tricks," in IEEE Computer Graphics and Applications, vol. 17, no. 4, pp. 80-84 \doi{10.1109/38.595279}
#'
#' J. T. Coonen, (1984) §2.3 "A Poor Man's Logarithm" in Contributions to a Proposed Standard for Binary Floating-Point Arithmetic. PhD Thesis, University of California Berkeley
#'
#' S. Summit, (2023) Answer to "Why does the integer representation of a floating point number offer a piecewise linear approximation to the logarithm?" Stack Overflow. \url{https://stackoverflow.com/a/75772363/1188479}
#'
#' @examples
#' \donttest{
#' # Custom Newton-Raphson parameters
#' result <- frsr_compute(c(1, 4, 9, 16), magic = 0x5f375a86, NRmax = 2, A = 1.6, B = 0.6)
#' ## result is a vector of length 4
#' print(result)
#' # [1] 0.9990148 0.4995074 0.3337626 0.2497537
#'
#' # Optional detail
#' result.df <- frsr_compute(c(pi, 2^-31, 0.4, 6.02e23), detail = TRUE)
#' ## result.df is a dataframe with 4 rows and 8 columns
#' print(result.df)
#' #         input      initial    after_one        final        error         enre          diff iters
#' # 1 3.141593e+00 5.735160e-01 5.639570e-01 5.639570e-01 0.0004122326 0.0004122326 -9.558976e-03     1
#' # 2 4.656613e-10 4.693787e+04 4.632937e+04 4.632937e+04 0.0002499308 0.0002499308 -6.085039e+02     1
#' # 3 4.000000e-01 1.632430e+00 1.578616e+00 1.578616e+00 0.0015955754 0.0015955754 -5.381417e-02     1
#' # 4 6.020000e+23 1.306493e-12 1.288484e-12 1.288484e-12 0.0002823969 0.0002823969 -1.800925e-14     1
#' # }
#' @name frsr_compute
NULL

#' @rdname frsr_compute
#' @export
frsr_compute <- function(x, magic = 0x5f3759df, NRmax = 1,
                         A = 1.5, B = 0.5, tol = 0,
                         detail = FALSE, keep_params = FALSE) {
  arg_df <- data.frame(x = x, magic = magic,
                       NRmax = NRmax, tol = tol,
                       A = A, B = B)
  .Call('_frsrr_frsr', PACKAGE = 'frsrr',
          arg_df, keep_params) -> result
  if (detail) {
    return(result)
  } else {
    return(result$final)
  }
}

#' @rdname frsr_compute
#' @export
frsr <- function(x, magic = 0x5f3759df, NRmax = 1,
                 A = 1.5, B = 0.5, tol = 0,
                 detail = FALSE, keep_params = FALSE) {
  frsr_compute(x = x, magic = magic, NRmax = NRmax,
               A = A, B = B, tol = tol,
               detail = detail, keep_params = keep_params)
}
