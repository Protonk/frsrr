#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "frsr.h"

using namespace Rcpp;

namespace phase_detail {

// ---- Statistics helpers -----------------------------------------------------
// These operate on pre-sorted windows so we can keep allocation costs flat.
inline double QuantileFromSorted(const double* sorted, std::size_t count, double q) {
    if (count == 0) {
        return NA_REAL;
    }
    if (count == 1 || q <= 0.0) {
        return sorted[0];
    }
    if (q >= 1.0) {
        return sorted[count - 1];
    }
    const double position = q * static_cast<double>(count - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const double fraction = position - static_cast<double>(lower);
    double value = sorted[lower];
    if (fraction > 0.0 && lower + 1u < count) {
        value += fraction * (sorted[lower + 1] - sorted[lower]);
    }
    return value;
}

inline double MedianSorted(double* sorted, std::size_t count) {
    if (count == 0) {
        return NA_REAL;
    }
    if (count % 2 == 1) {
        return sorted[count / 2];
    }
    const double mid = sorted[count / 2 - 1];
    const double mid_next = sorted[count / 2];
    return 0.5 * (mid + mid_next);
}

// ---- Input validation helpers ----------------------------------------------
inline std::vector<int> ValidateExponents(const IntegerVector& exponents) {
    std::vector<int> values(exponents.size());
    for (int i = 0; i < exponents.size(); ++i) {
        if (exponents[i] == NA_INTEGER) {
            stop("`exponents` cannot contain NA");
        }
        if (exponents[i] < -126 || exponents[i] > 127) {
            stop("`exponents` must be within [-126, 127] to avoid subnormals/infinities");
        }
        values[i] = exponents[i];
    }
    return values;
}

inline std::vector<int> ValidateMagics(const IntegerVector& magics) {
    std::vector<int> values(magics.size());
    for (int i = 0; i < magics.size(); ++i) {
        if (magics[i] == NA_INTEGER) {
            stop("`magics` cannot contain NA");
        }
        values[i] = magics[i];
    }
    return values;
}

// ---- RNG helpers ------------------------------------------------------------
inline std::mt19937_64::result_type SeedFromRandomDevice() {
    std::random_device rd;
    // Two draws to give the 64-bit engine a wider seed surface.
    uint64_t high = static_cast<uint64_t>(rd()) << 32;
    uint64_t low = static_cast<uint64_t>(rd());
    return high ^ low;
}

inline std::mt19937_64 InitializeRng(const Nullable<NumericVector>& seed_nullable) {
    std::mt19937_64 rng;
    if (seed_nullable.isNotNull()) {
        NumericVector seed_vec(seed_nullable);
        if (seed_vec.size() != 1) {
            stop("`seed` must be NULL or a length-1 numeric value");
        }
        double seed_value = seed_vec[0];
        if (!std::isfinite(seed_value)) {
            stop("`seed` must be finite when provided");
        }
        int64_t rounded = static_cast<int64_t>(std::llround(seed_value));
        rng.seed(static_cast<std::mt19937_64::result_type>(rounded));
        return rng;
    }
    rng.seed(SeedFromRandomDevice());
    return rng;
}

}  // namespace phase_detail

// [[Rcpp::export]]
List phase_orchestrator(int phases,
                        IntegerVector exponents,
                        int per_cell,
                        IntegerVector magics,
                        double q,
                        int NRmax,
                        Nullable<NumericVector> seed_nullable) {
    if (phases <= 0) {
        stop("`phases` must be positive");
    }
    if (per_cell <= 0) {
        stop("`per_cell` must be positive");
    }
    if (exponents.size() == 0) {
        stop("`exponents` must contain at least one exponent");
    }
    if (magics.size() == 0) {
        stop("`magics` must contain at least one candidate");
    }
    if (!std::isfinite(q) || q <= 0.0 || q > 1.0) {
        stop("`q` must satisfy 0 < q <= 1");
    }
    if (NRmax < 0) {
        stop("`NRmax` must be non-negative");
    }

    // Copy R inputs into STL containers once so inner loops can stay pointer-friendly.
    std::vector<int> exponent_values = phase_detail::ValidateExponents(exponents);
    std::vector<int> magic_values = phase_detail::ValidateMagics(magics);

    // Keep sampling deterministic when the caller supplied a seed; otherwise rely on entropy.
    std::mt19937_64 rng = phase_detail::InitializeRng(seed_nullable);
    std::uniform_real_distribution<double> unit_dist(0.0, 1.0);

    const int num_exponents = static_cast<int>(exponent_values.size());
    const std::size_t samples_per_phase =
        static_cast<std::size_t>(num_exponents) * static_cast<std::size_t>(per_cell);
    if (samples_per_phase == 0) {
        stop("`per_cell` * length(`exponents`) must be non-zero");
    }
    if (samples_per_phase > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        stop("`per_cell` * length(`exponents`) must be representable as a 32-bit integer");
    }
    const std::size_t total_samples =
        samples_per_phase * static_cast<std::size_t>(phases);
    if (total_samples / samples_per_phase != static_cast<std::size_t>(phases)) {
        stop("Requested grid is too large for available memory");
    }

    // All scratch storage is allocated once to avoid quadratic reallocations while scanning magics.
    std::vector<double> phase_abs(total_samples);
    std::vector<double> phase_signed(total_samples);
    std::vector<double> phase_sum(static_cast<std::size_t>(phases));
    std::vector<std::size_t> offsets(static_cast<std::size_t>(phases));
    std::vector<double> candidate_heat(static_cast<std::size_t>(num_exponents) *
                                       static_cast<std::size_t>(phases));
    std::vector<double> cell_errors(static_cast<std::size_t>(per_cell));

    // When a candidate wins we copy its summaries into these vectors for the final report.
    std::vector<double> best_q_abs;
    std::vector<double> best_means;
    std::vector<double> best_medians;
    std::vector<double> best_heat;
    int best_magic = magic_values[0];
    double best_J = std::numeric_limits<double>::infinity();
    double best_R = std::numeric_limits<double>::infinity();
    bool have_best = false;

    const double phase_width = 1.0 / static_cast<double>(phases);

    for (int magic_idx = 0; magic_idx < static_cast<int>(magic_values.size()); ++magic_idx) {
        checkUserInterrupt();

        const int magic = magic_values[magic_idx];
        // Reset per-phase accumulators before reusing the buffers.
        std::fill(phase_sum.begin(), phase_sum.end(), 0.0);
        std::fill(offsets.begin(), offsets.end(), 0u);

        for (int exp_idx = 0; exp_idx < num_exponents; ++exp_idx) {
            const int exponent = exponent_values[exp_idx];

            for (int phase_idx = 0; phase_idx < phases; ++phase_idx) {
                const double phase_low = static_cast<double>(phase_idx) * phase_width;
                const double phase_high = phase_low + phase_width;

                for (int sample_idx = 0; sample_idx < per_cell; ++sample_idx) {
                    // Pick a mantissa inside the current phase slab to keep coverage uniform.
                    double u = unit_dist(rng);
                    double frac = phase_low + u * (phase_high - phase_low);
                    if (frac >= 1.0) {
                        frac = std::nextafter(1.0, 0.0);
                    }
                    // Mantissa rebuild keeps the sample inside the desired phase
                    // slice while std::ldexp reattaches the exponent under test.
                    const double mantissa = std::exp2(frac);
                    const double x_val = std::ldexp(mantissa, exponent);
                    const float xf = static_cast<float>(x_val);
                    const float approx = frsr0(xf, static_cast<uint32_t>(magic), NRmax);
                    const double exact = 1.0 / std::sqrt(static_cast<double>(xf));
                    const double epsilon = (static_cast<double>(approx) - exact) / exact;
                    const double abs_eps = std::abs(epsilon);

                    const std::size_t base = static_cast<std::size_t>(phase_idx) * samples_per_phase;
                    std::size_t offset = offsets[phase_idx];
                    if (offset >= samples_per_phase) {
                        stop("Internal error: exceeded allocated samples for a phase");
                    }
                    const std::size_t location = base + offset;
                    phase_abs[location] = abs_eps;
                    phase_signed[location] = epsilon;
                    offsets[phase_idx] = offset + 1u;
                    phase_sum[phase_idx] += epsilon;
                    cell_errors[static_cast<std::size_t>(sample_idx)] = epsilon;
                }

                // Median signed error for the exponent/phase cell becomes the heat-map entry.
                std::sort(cell_errors.begin(), cell_errors.end());
                const double cell_median =
                    phase_detail::MedianSorted(cell_errors.data(), cell_errors.size());
                candidate_heat[static_cast<std::size_t>(exp_idx) * static_cast<std::size_t>(phases) +
                               static_cast<std::size_t>(phase_idx)] = cell_median;
            }
        }

        std::vector<double> q_abs(static_cast<std::size_t>(phases));
        std::vector<double> means(static_cast<std::size_t>(phases));
        std::vector<double> medians(static_cast<std::size_t>(phases));

        for (int phase_idx = 0; phase_idx < phases; ++phase_idx) {
            const std::size_t base = static_cast<std::size_t>(phase_idx) * samples_per_phase;
            double* abs_begin = phase_abs.data() + base;
            double* signed_begin = phase_signed.data() + base;

            std::sort(abs_begin, abs_begin + samples_per_phase);
            std::sort(signed_begin, signed_begin + samples_per_phase);

            q_abs[phase_idx] = phase_detail::QuantileFromSorted(abs_begin, samples_per_phase, q);
            medians[phase_idx] = phase_detail::MedianSorted(signed_begin, samples_per_phase);
            means[phase_idx] = phase_sum[phase_idx] / static_cast<double>(samples_per_phase);
        }

        // J: the worst absolute error quantile across phases (lower is better).
        double J = 0.0;
        for (double val : q_abs) {
            if (val > J) {
                J = val;
            }
        }

        // Roughness: how jagged the signed means are as we walk around the circle of phases.
        double roughness = 0.0;
        for (int phase_idx = 0; phase_idx < phases; ++phase_idx) {
            const int prev = (phase_idx == 0) ? (phases - 1) : (phase_idx - 1);
            roughness += std::abs(means[phase_idx] - means[prev]);
        }
        roughness /= static_cast<double>(phases);

        // Prefer smaller J; break ties with the smoothest phase-to-phase transitions.
        if (!have_best || J < best_J || (std::abs(J - best_J) <= std::numeric_limits<double>::epsilon() && roughness < best_R)) {
            have_best = true;
            best_J = J;
            best_R = roughness;
            best_magic = magic;
            best_q_abs = q_abs;
            best_means = means;
            best_medians = medians;
            best_heat = candidate_heat;
        }
    }

    if (!have_best) {
        stop("No feasible magic constant found");
    }

    IntegerVector phase_ids(phases);
    for (int i = 0; i < phases; ++i) {
        phase_ids[i] = i + 1;
    }
    const int samples_per_phase_int = static_cast<int>(samples_per_phase);
    IntegerVector counts(phases, samples_per_phase_int);

    NumericVector q_abs_vec(best_q_abs.begin(), best_q_abs.end());
    NumericVector mean_vec(best_means.begin(), best_means.end());
    NumericVector median_vec(best_medians.begin(), best_medians.end());

    // Phase-level summary table for the winning magic.
    DataFrame phase_tbl = DataFrame::create(
        Named("phase_id") = phase_ids,
        Named("q_abs") = q_abs_vec,
        Named("mean_signed") = mean_vec,
        Named("median_signed") = median_vec,
        Named("n") = counts
    );

    // Heat map captures how each exponent behaved inside each phase slice.
    NumericMatrix heat(num_exponents, phases);
    CharacterVector row_names(num_exponents);
    CharacterVector col_names(phases);
    for (int r = 0; r < num_exponents; ++r) {
        row_names[r] = std::to_string(exponent_values[r]);
        for (int c = 0; c < phases; ++c) {
            if (r == 0) {
                col_names[c] = std::string("phase_") + std::to_string(c + 1);
            }
            const std::size_t idx =
                static_cast<std::size_t>(r) * static_cast<std::size_t>(phases) +
                static_cast<std::size_t>(c);
            heat(r, c) = best_heat[idx];
        }
    }
    heat.attr("dimnames") = List::create(row_names, col_names);

    return List::create(
        Named("magic") = best_magic,
        Named("J") = best_J,
        Named("R") = best_R,
        Named("phase_tbl") = phase_tbl,
        Named("heat") = heat
    );
}
