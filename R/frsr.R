#' @useDynLib frsrr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL

#' Fast Reciprocal Square Root (FRSR)
#'
#' A parameterized Fast Reciprocal Square Root algorithm written in C++.
#'
#' @param x Numeric inputs that convert to positive normal float32 values.
#'   Original values must not exceed `(2 - 2^-23) * 2^127`. Nonfinite values,
#'   zero and subnormal float32 inputs are unsupported.
#' @param magic Integer restoring constant. Default is 0x5f3759df.
#' @param NRmax Integer specifying the maximum Newton-Raphson iterations. Default is 1.
#' @param A Newton-Raphson parameter where \eqn{(A - B * x * y_n^2)}. Default is \code{1.5}.
#' @param B Newton-Raphson parameter. Default is \code{0.5}.
#' @param tol The absolute relative error at which to stop early. Default is 0 (no early stopping).
#' @param detail Logical. If \code{TRUE}, a data frame with detailed results is returned. Default is \code{FALSE}.
#' @param keep_params Logical. If \code{TRUE}, generation parameters are included in detailed output. Default is \code{FALSE}.
#' @param threads Positive integer passed to \code{RcppParallel::setThreadOptions()} before invoking
#'   the parallel C++ workers. Defaults to \code{getOption("frsrr.threads", NA)}, which falls back
#'   to \code{RcppParallel::defaultNumThreads()} when unset.
#'
#' @return
#' \code{frsr} returns a numeric vector of \code{length(x)}.
#'
#' If \code{detail = TRUE}, returns a data frame of \code{length(x)} rows with columns:
#'     \item{input}{The input values}
#'     \item{initial}{Initial approximation from integer operations}
#'     \item{after_one}{Result after one iteration of Newton-Raphson}
#'     \item{final}{Result from final iteration}
#'     \item{error}{Double-precision absolute relative error against the reciprocal
#'       square root of the float32-rounded input; Inf for a nonfinite approximation}
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
#' Inputs and coefficients are rounded to float32. The bit seed uses IEEE-754
#' binary32 representation; each operation in `y * (A - ((B * x) * y) * y)`
#' rounds to float32, with no fused multiply-subtract. This assumes the usual
#' round-to-nearest environment with gradual underflow; altered rounding modes,
#' flush-to-zero and fast-math builds are unsupported. Historical-machine or
#' cross-platform bitwise equivalence is not promised.
#'
#' The reference is `1 / sqrt(x_float32)` evaluated in double precision,
#' so `error` measures the algorithm rather than input-conversion error. The
#' `input`, `A` and `B` columns retain the supplied values. `diff` is a double
#' subtraction of the last two float32 approximations. `NRmax = 0` returns the
#' seed, zero iterations, and NA for `after_one` and `diff`. A positive `tol`
#' is checked after each refinement, so at least one step runs when `NRmax > 0`.
#' Finite exploratory coefficients and magic constants may yield poor or
#' nonfinite approximations; these are returned, with infinite error for the latter.
#' The former `enre` column (a separate mantissa-only experiment) has been removed.
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
#' frsr(c(1, 2, 4), magic = 0x5f3759df, NRmax = 0, tol = 0, threads = 1)
#' frsr(c(pi, 0.4), magic = 0x5f3759df, NRmax = 2, tol = 0,
#'      detail = TRUE, threads = 1)
#' @name frsr
NULL

#' @rdname frsr
#' @export
frsr <- function(x, magic = 0x5f3759df, NRmax = 1,
                 A = 1.5, B = 0.5, tol = 0,
                 detail = FALSE, keep_params = FALSE,
                 threads = getOption("frsrr.threads", NA_integer_)) {
  if (!is.numeric(x) || is.complex(x) || !is.numeric(A) || is.complex(A) ||
      !is.numeric(B) || is.complex(B) || !is.numeric(tol) || is.complex(tol)) {
    stop("`x`, `A`, `B` and `tol` must be numeric", call. = FALSE)
  }
  if (!is.numeric(NRmax) || any(!is.finite(NRmax) | NRmax < 0 |
                              NRmax > .Machine$integer.max | NRmax != trunc(NRmax))) {
    stop("`NRmax` must contain non-negative integers", call. = FALSE)
  }
  if (!is.numeric(magic) || any(!is.finite(magic) | abs(magic) > .Machine$integer.max |
                              magic != trunc(magic))) {
    stop("`magic` must contain non-missing R integers", call. = FALSE)
  }
  threads <- frsrr_configure_threads(threads)
  if (!length(x)) {
    magic <- NRmax <- integer()
    A <- B <- tol <- numeric()
  }
  arg_df <- data.frame(x = x, magic = as.integer(magic),
                       NRmax = as.integer(NRmax), tol = tol,
                       A = A, B = B)
  # Consolidate user inputs into a rectangular data frame so the C++ path can
  # run one parallel sweep; recycling happens here in R where the rules are
  # explicit instead of inside the worker threads.
  .Call('_frsrr_frsr', PACKAGE = 'frsrr',
          arg_df, keep_params) -> result
  if (detail) {
    # Returning the full data frame preserves intermediate diagnostics that are
    # expensive to recompute (e.g., per-iteration errors) and mirrors the
    # documented detail schema.
    return(result)
  } else {
    # The fast path only surfaces the final approximation vector so this
    # function can drop in as a near-drop-in replacement for 1 / sqrt(x).
    return(result$final)
  }
}
