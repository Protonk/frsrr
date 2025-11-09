#' Phase-Oriented Magic Search
#'
#' `frsr_phase()` stratifies floating-point samples by exponent and phase
#' (fractional log2) bins, evaluates every candidate magic constant, and returns
#' the magic that minimizes the worst phase-wise tail error. Statistics
#' are returned to help diagnose phase bias.
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
#' @param seed Optional numeric scalar. When supplied the sampler becomes
#'   deterministic; otherwise a fresh seed is drawn.
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
#'   magics = magics,
#'   seed = 1
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
                       NRmax = 0L,
                       seed = NULL) {
  phases <- as.integer(phases)[1]
  exponents <- as.integer(exponents)
  per_cell <- as.integer(per_cell)[1]
  magics <- as.integer(magics)
  q <- as.numeric(q)[1]
  NRmax <- as.integer(NRmax)[1]
  seed_arg <- if (is.null(seed)) NULL else as.numeric(seed)[1]

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
    NRmax,
    seed_arg
  )
}
