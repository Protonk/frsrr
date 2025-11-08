#' Phase-Oriented Magic Search
#'
#' `frsr_phase()` stratifies floating-point samples by exponent and phase
#' (fractional log2) bins, evaluates every candidate magic constant, and returns
#' the magic that minimizes the worst phase-wise tail error. Phase statistics
#' and a coarse heatmap are returned to help diagnose phase bias.
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
#' grid <- default_magic_grid()
#' phase_fit <- frsr_phase(
#'   phases = 64,
#'   exponents = -10:10,
#'   per_cell = 32,
#'   magics = grid[seq(1, length(grid), by = 8)],
#'   seed = 1
#' )
#' print(phase_fit$magic)
#' }
#'
#' @export
frsr_phase <- function(phases = 128L,
                       exponents = -62L:62L,
                       per_cell = 64L,
                       magics = default_magic_grid(),
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

#' @rdname frsr_phase
#' @param lower,upper Inclusive bounds for the magic range (integers). Defaults
#'   match the bin search helper.
#' @param step Step size used to traverse the range. Must be strictly positive.
#'
#' @export
default_magic_grid <- function(lower = 1596980000L,
                               upper = 1598050000L,
                               step = 512L) {
  lower <- as.integer(lower)[1]
  upper <- as.integer(upper)[1]
  step <- as.integer(step)[1]
  stride <- abs(step)
  if (is.na(stride) || stride == 0L) {
    stride <- 1L
  }
  if (lower <= upper) {
    seq.int(lower, upper, by = stride)
  } else {
    seq.int(lower, upper, by = -stride)
  }
}

#' @rdname frsr_phase
#' @param heat Matrix returned by `frsr_phase()` (or compatible numeric matrix).
#' @param palette Character vector of colors used to render the heatmap. The
#'   default is a diverging ramp chosen for "error carpet" plots.
#' @param xlab,ylab Axes labels passed to [graphics::image()].
#' @param ... Additional arguments forwarded to [graphics::image()].
#'
#' @export
frsr_phase_heatmap <- function(heat,
                               palette = grDevices::colorRampPalette(
                                 c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026")
                               )(256),
                               xlab = "Phase Bin",
                               ylab = "Exponent",
                               ...) {
  heat <- as.matrix(heat)
  storage.mode(heat) <- "double"
  x <- seq_len(ncol(heat))
  row_names <- rownames(heat)
  y <- seq_len(nrow(heat))
  if (!is.null(row_names)) {
    suppressWarnings({
      y_numeric <- as.numeric(row_names)
    })
    if (!anyNA(y_numeric)) {
      y <- y_numeric
    }
  }
  graphics::image(
    x = x,
    y = y,
    z = t(heat),
    col = palette,
    xlab = xlab,
    ylab = ylab,
    useRaster = TRUE,
    ...
  )
  invisible(heat)
}
