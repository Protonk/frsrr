#' @importFrom RcppParallel defaultNumThreads setThreadOptions
NULL

#' Normalize thread counts before touching the parallel backends.
#'
#' Ensures user-provided `threads` values collapse to a single positive integer,
#' applies the package option fallback, and updates `RcppParallel` in one place
#' so every exported function inherits the same safety checks.
#'
#' @param threads Optional scalar overriding `getOption("frsrr.threads")`.
#'
#' @return
#' The sanitized integer thread count used for downstream calls.
#'
#' @keywords internal
#' @noRd
frsrr_configure_threads <- function(threads = getOption("frsrr.threads", NA_integer_)) {
  if (is.null(threads) || (length(threads) == 0L) || is.na(threads[1])) {
    threads <- RcppParallel::defaultNumThreads()
  }
  threads <- as.integer(threads)[1]
  if (!is.finite(threads) || threads < 1L) {
    stop("`threads` must be a positive integer", call. = FALSE)
  }
  RcppParallel::setThreadOptions(numThreads = threads)
  threads
}
