#' Phase-Oriented Magic Search
#'
#' Stratifies floating-point samples by exponent and phase
#' (fractional log2) bins, evaluates every candidate magic constant, and returns
#' the magic that minimizes the worst phase-wise tail error. Statistics
#' are returned to help diagnose phase bias.
#'
#' @details
#' log2 phases are a useful way to this function as it contains a
#' linear approximation to log2 whose performance varies greatly
#' depending on proximity to powers of two.
#'
#' @param phases Number of equally sized phase bins spanning `[0, 1)`. Must be
#'   a positive integer.
#' @param exponents Integer vector of unbiased exponents. Values must keep the
#'   generated floats in the normalized range (i.e. `[-126, 127]`).
#' @param per_cell Samples to generate for each `(phase, exponent)` cell.
#' @param magics Integer vector of candidate magic constants to evaluate.
#' @param q Quantile applied to the absolute relative errors inside each phase.
#'   Must satisfy `0 < q <= 1`.
#' @param NRmax Newton steps applied through [`frsr0()`]. Reuses existing
#'   Newton logic; the default sticks with the raw restoring constant path.
#'
#' @return
#' `frsr_phase()` returns a list with components:
#' \describe{
#'   \item{magic}{Chosen magic constant (integer).}
#'   \item{J}{Objective value, i.e. the maximum quantile across phases.}
#'   \item{R}{Mean absolute difference between consecutive phase-wise mean errors.}
#'   \item{phase_tbl}{Data frame with the per-phase quantile, mean, median, and
#'     sample count.}
#'   \item{heat}{A matrix indexed by exponents x phases storing per-cell median
#'     signed errors, suitable for plotting.}
#' }
#'
#' @examples
#' \donttest{
#' magics <- as.integer(seq.int(0x5f3750df, 0x5f3765df, by = 512L))
#' phase_fit <- frsr_phase(
#'   phases = 32L,
#'   exponents = -8L:8L,
#'   per_cell = 8L,
#'   magics = magics
#' )
#' phase_fit$magic
#' }
#'
#' @export
frsr_phase <- function(phases = 128L,
                       exponents = -62L:62L,
                       per_cell = 64L,
                       magics = c(1596980000L, 1598050000L),
                       q = 0.95,
                       NRmax = 0L) {
  phases <- as.integer(phases)[1]
  if (is.na(phases) || phases < 1L) {
    stop("`phases` must be a positive integer", call. = FALSE)
  }

  exponents <- as.integer(exponents)
  if (!length(exponents)) {
    stop("`exponents` must contain at least one value", call. = FALSE)
  }
  if (any(!is.finite(exponents))) {
    stop("`exponents` must be finite integers", call. = FALSE)
  }
  if (any(exponents < -126L | exponents > 127L)) {
    stop("`exponents` must stay within [-126, 127]", call. = FALSE)
  }

  per_cell <- as.integer(per_cell)[1]
  if (is.na(per_cell) || per_cell < 1L) {
    stop("`per_cell` must be a positive integer", call. = FALSE)
  }

  magics <- as.integer(magics)
  if (!length(magics) || any(is.na(magics))) {
    stop("`magics` must supply at least one integer constant", call. = FALSE)
  }

  q <- as.numeric(q)[1]
  if (!is.finite(q) || q <= 0 || q > 1) {
    stop("`q` must satisfy 0 < q <= 1", call. = FALSE)
  }

  NRmax <- as.integer(NRmax)[1]
  if (is.na(NRmax) || NRmax < 0L) {
    stop("`NRmax` must be a non-negative integer", call. = FALSE)
  }

  # Keep conversion/validation in R so the hot C++ path can assume scalars,
  # which avoids repeatedly checking lengths inside the tight sampling loops.
  .Call(
    '_frsrr_phase_orchestrator',
    PACKAGE = 'frsrr',
    phases,
    exponents,
    per_cell,
    magics,
    q,
    NRmax
  )
}
