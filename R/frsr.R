#' @useDynLib frsrr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
#' @rdname frsr
#' @export
frsr.detail <- function(x, magic = 0x5f3759df, NR = 1, A = 1.5, B = 0.5, tol = 0) {
  .Call('_frsrr_frsr_detail', PACKAGE = 'frsrr', x, magic, NR, A, B, tol)
}
#' @export
frsr <- function(x, magic = 0x5f3759df, NR = 1, A = 1.5, B = 0.5, tol = 0) {
  .Call('_frsrr_frsr', PACKAGE = 'frsrr', x, magic, NR, A, B, tol)
}