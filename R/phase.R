#' Phase-Oriented Magic Search
#'
#' Stratifies floating-point samples by exponent and phase
#' (fractional log2) bins, evaluates every candidate magic constant, and returns
#' the magic that minimizes the worst phase-wise tail error. Statistics
#' are returned to help diagnose phase bias.
#'
#' @details
#' Every candidate uses the same randomly sampled grid and the arithmetic/error
#' contract of [frsr()] with A = 1.5, B = 0.5 and tol = 0. Phase labels describe
#' the unrounded log2 draws; rounding can cross a phase boundary. Values are
#' clamped to the largest float in their exponent slab to avoid crossing into
#' the next exponent or overflowing at exponent 127.
#'
#' Candidates are ordered by J, then roughness R, then the smallest magic
#' integer, using exact ties. Permuting or inserting candidates does not change
#' the inputs assigned to existing candidates with the same seed. A candidate
#' with any nonfinite approximation is excluded; the call errors if none remain.
#' The winner is best among the tested candidates on this sampled grid.
#' The `settings` component retains NRmax, q, per_cell, exponents and the distinct
#' candidate magics. Use `set.seed()` to reproduce the grid.
#'
#' @param phases Number of equally sized phase bins spanning `[0, 1)`. Must be
#'   a positive integer.
#' @param exponents Integer vector of unbiased exponents. Values must keep the
#'   generated floats in the normalized range (i.e. `[-126, 127]`).
#' @param per_cell Samples to generate for each `(phase, exponent)` cell.
#' @param magics Integer vector of candidate magic constants to evaluate.
#' @param q Quantile applied to the absolute relative errors inside each phase.
#'   Must satisfy `0 < q <= 1`.
#' @param NRmax Number of float32 Newton steps. The default returns the raw
#'   restoring-constant approximation without refinement.
#'
#' @return
#' `frsr_phase()` returns a list with components:
#' \describe{
#'   \item{magic}{Chosen magic constant (integer).}
#'   \item{J}{Objective value, i.e. the maximum quantile across phases.}
#'   \item{R}{Mean absolute difference between consecutive phase-wise mean errors.}
#'   \item{phase_tbl}{Data frame with the per-phase quantile, mean, median, and
#'     sample count.}
#'   \item{settings}{Configuration needed to interpret the sampled search.}
#'   \item{heat}{A matrix indexed by exponents x phases storing per-cell median
#'     signed errors, suitable for plotting.}
#' }
#'
#' @examples
#' set.seed(42)
#' frsr_phase(phases = 8, exponents = -2:2, per_cell = 4,
#'            magics = c(0x5f3759df, 0x5f375a86), NRmax = 1)
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

  NRmax <- as.numeric(NRmax)[1]
  if (!is.finite(NRmax) || NRmax < 0 || NRmax > .Machine$integer.max || NRmax != trunc(NRmax)) {
    stop("`NRmax` must be a non-negative integer", call. = FALSE)
  }
  NRmax <- as.integer(NRmax)

  # Keep conversion/validation in R so the hot C++ path can assume scalars,
  # which avoids repeatedly checking lengths inside the tight sampling loops.
  result <- .Call(
    '_frsrr_phase_orchestrator',
    PACKAGE = 'frsrr',
    phases,
    exponents,
    per_cell,
    magics,
    q,
    NRmax
  )
  result$settings <- list(NRmax = NRmax, q = q, per_cell = per_cell,
                          exponents = exponents, magics = sort(unique(magics)))
  result
}
