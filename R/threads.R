#' @importFrom RcppParallel defaultNumThreads setThreadOptions
NULL

# Internal helper to clamp requested thread counts to a safe value.
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
