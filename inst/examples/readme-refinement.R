# Worked README experiment; requires frsrr >= 1.1.0.
# Source this file, then call plot_refinement() to draw the one-step comparison.
# All inputs are deterministic: no random sampling or RNG seed is needed.
library(frsrr)

float32 <- function(x) readBin(writeBin(as.double(x), raw(), size = 4),
                             "double", n = length(x), size = 4)
x_scan <- 1 + (0:3071) / 1024
measure <- function(x, A, steps = 1L) {
  fit <- frsr(x, magic = 0x5f3759df, A = A, B = 0.5, NRmax = steps,
              tol = 0, threads = 1, detail = TRUE, keep_params = TRUE)
  reference <- 1 / sqrt(x)
  fit$signed_error <- (fit$final - reference) / reference
  fit
}
ordinary_scan <- measure(x_scan, 1.5)

# Compare 121 candidates on the same grid, including the ordinary coefficient.
# Round before selection so A_adjusted stores the coefficient used by C++.
candidates <- unique(float32(1.5 + (0:120) * 1e-5))
scan_max <- vapply(candidates, function(A) max(measure(x_scan, A)$error),
                   numeric(1))
A_adjusted <- candidates[which.min(scan_max)]
d <- A_adjusted - 1.5

# Eight times the density, including every scan input and 21,504 new locations.
# Freeze A_adjusted for all evaluation runs; tol = 0 forces every requested step.
x_eval <- 1 + (0:24575) / 8192
steps <- c(1L, 2L, 4L)
ordinary <- lapply(steps, function(n) measure(x_eval, 1.5, n))
adjusted <- lapply(steps, function(n) measure(x_eval, A_adjusted, n))
max_error <- function(fit) max(fit$error)
comparison <- data.frame(steps,
                          ordinary = vapply(ordinary, max_error, numeric(1)),
                          adjusted = vapply(adjusted, max_error, numeric(1)))
predicted <- sqrt(1 + 2 * d) - 1

# Draw only one-step signed errors; the table uses these same evaluation runs.
plot_refinement <- function() {
  errors <- cbind(ordinary[[1]]$signed_error, adjusted[[1]]$signed_error)
  colours <- c("#0072B2", "#D55E00")
  matplot(x_eval, 1000 * errors, type = "l", lty = c(1, 2), lwd = 2,
          col = colours, xlim = c(1, 4), ylim = c(-1.8, 1.3), xlab = "Input x",
          ylab = "Signed relative error (x 1,000)")
  abline(h = 0, col = "grey40", lty = 3)
  legend("top", c("Ordinary: A = 1.5", "Adjusted: selected A"),
         col = colours, lty = c(1, 2), lwd = 2, bty = "n", horiz = TRUE)
}

# From the repository root, after loading the local package:
# source("inst/examples/readme-refinement.R")
# png("man/figures/readme-refinement.png", width = 1200, height = 720, res = 150,
#     type = "cairo")  # the macOS quartz device drops dashes on dense lines
# plot_refinement()
# dev.off()
# sprintf("A = %.17g; d = %.17g", A_adjusted, d)
# c(ordinary = max(ordinary_scan$error), adjusted = min(scan_max))
# print(comparison, digits = 6, row.names = FALSE)
# four_step <- adjusted[[3]]$signed_error
# c(predicted = predicted, min = min(four_step), max = max(four_step))
# The README's output blocks must match these printed results exactly; the
# README test runs its code and reports any block that needs pasting in again.
