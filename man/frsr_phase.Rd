\name{frsr_phase}
\alias{frsr_phase}
\alias{default_magic_grid}
\alias{frsr_phase_heatmap}
\title{Phase-Oriented Magic Search}
\description{
Optimize a single magic constant by probing worst-case errors across phase bins
and summarize the per-phase behavior for downstream analysis.
}
\usage{
frsr_phase(phases = 128L, exponents = -62L:62L, per_cell = 64L,
  magics = default_magic_grid(), q = 0.95, NRmax = 0L, seed = NULL)

default_magic_grid(lower = 1596980000L, upper = 1598050000L,
  step = 512L)

frsr_phase_heatmap(heat, palette = grDevices::colorRampPalette(c("#313695",
  "#74add1", "#ffffbf", "#f46d43", "#a50026"))(256), xlab = "Phase Bin",
  ylab = "Exponent", ...)
}
\arguments{
\item{phases}{Number of equally sized phase bins covering the fractional log2
range \verb{[0, 1)}.}

\item{exponents}{Integer vector of unbiased exponents. Values outside
\verb{[-126, 127]} are rejected because they produce subnormals or infinities.}

\item{per_cell}{Number of samples to draw for each combination of phase and
exponent.}

\item{magics}{Integer vector of candidate magic constants.}

\item{q}{Quantile applied to the absolute relative errors within each phase.
Must satisfy \verb{0 < q <= 1}.}

\item{NRmax}{Maximum number of Newton iterations passed through to
\code{frsr0()}.}

\item{seed}{Optional numeric scalar; when present the sampler is deterministic.}

\item{lower, upper}{Inclusive bounds for the helper grid.}

\item{step}{Stride used to traverse the helper grid. Must be strictly positive.}

\item{heat}{Numeric matrix produced by \code{frsr_phase()} (or any compatible
matrix where rows correspond to exponents and columns correspond to phase bins).}

\item{palette}{Character vector of colors used by the heatmap helper.}

\item{xlab, ylab}{Axis labels for the plotting helper.}

\item{...}{Additional arguments forwarded to \code{graphics::image()}.}
}
\value{
\code{frsr_phase()} returns a list containing the chosen magic, the objective
value (\eqn{J}), the roughness metric (\eqn{R}), the per-phase summary table,
and a matrix of exponent/phase median signed errors.

\code{default_magic_grid()} returns an integer vector of candidate magic
constants spaced over the specified interval.

\code{frsr_phase_heatmap()} draws a heatmap (via \code{graphics::image()}) and
returns the original matrix invisibly.
}
\examples{
\donttest{
grid <- default_magic_grid()
phase_fit <- frsr_phase(
  phases = 32,
  exponents = -8:8,
  per_cell = 16,
  magics = grid[seq(1, length(grid), by = 16)],
  seed = 42
)
print(phase_fit$magic)
frsr_phase_heatmap(phase_fit$heat)
}
}
