#include "search.h"

#include <RcppParallel.h>
#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "frsr.h"

using namespace Rcpp;
using namespace RcppParallel;

namespace {

// Convert an IEEE-754 binary32 payload encoded as an integer into a float.
inline float FromBits(uint32_t bits) {
    return std::bit_cast<float>(bits);
}

// Draw a random offset in [0, range) using R's RNG. The helper guards against
// the corner case where the random draw is exactly equal to the exclusive upper
// bound by clamping the result back inside the interval.
inline uint32_t sample_offset(uint32_t range) {
    if (range == 0u) {
        return 0u;
    }
    double u = unif_rand();
    double scaled = u * static_cast<double>(range);
    uint32_t draw = static_cast<uint32_t>(scaled);
    if (draw >= range) {
        draw = range - 1u;
    }
    return draw;
}

// Evaluate the absolute relative errors for the selected magic constant across
// every sampled float and accumulate only the scalars requested by the caller.
struct RequestedMetrics {
    bool want_sum;
    bool want_mean;
    bool want_rms;
    bool want_variance;
};

struct MetricAccumulators {
    double sum = 0.0;
    double sum_sq = 0.0;
};

MetricAccumulators accumulate_errors(const NumericVector& floats,
                                     uint32_t magic,
                                     int NRmax,
                                     const RequestedMetrics& request) {
    MetricAccumulators acc;
    if (!request.want_sum && !request.want_rms && !request.want_variance) {
        return acc;
    }

    for (int j = 0; j < floats.length(); ++j) {
        float x = floats[j];
        float approx = frsr0(x, magic, NRmax);
        float actual = 1.0f / std::sqrt(x);
        float rel_error = std::abs((approx - actual) / actual);

        if (request.want_sum || request.want_mean || request.want_variance) {
            acc.sum += static_cast<double>(rel_error);
        }
        if (request.want_rms || request.want_variance) {
            acc.sum_sq += static_cast<double>(rel_error) * static_cast<double>(rel_error);
        }
    }

    return acc;
}

}  // namespace

// [[Rcpp::export]]
// Enumerate float32 values with exponents inside [low, high] and sample from that
// discrete population. The sampler walks the exponent strata explicitly so that
// rare exponents are reachable even when the caller draws only a handful of
// samples.
NumericVector bounded_stratified_sample(int n, double low, double high, bool weighted) {
    if (n < 0) {
        throw std::invalid_argument("`n` must be non-negative");
    }
    if (!std::isfinite(low) || !std::isfinite(high)) {
        throw std::invalid_argument("`low` and `high` must be finite");
    }
    if (high <= low) {
        throw std::invalid_argument("`high` must be greater than `low`");
    }
    if (low < -126 || high > 128) {
        throw std::invalid_argument("Bounds must satisfy -126 <= low < high <= 128");
    }

    const long double lower_bound = static_cast<long double>(std::exp2(low));
    const long double upper_bound = static_cast<long double>(std::exp2(high));

    struct Stratum {
        int exponent;
        uint32_t smin;
        uint32_t smax;
        uint32_t count;
    };

    constexpr uint32_t SigCount = 1u << 23;
    const long double SigCountLD = static_cast<long double>(SigCount);
    const int emin = std::max(-126, static_cast<int>(std::floor(low)));
    const int emax = std::min(127, static_cast<int>(std::ceil(high)) - 1);

    std::vector<Stratum> strata;
    strata.reserve(std::max(0, emax - emin + 1));

    for (int e = emin; e <= emax; ++e) {
        long double twoe = std::ldexp(1.0L, e);

        uint32_t smin = 0u;
        if (e == emin) {
            long double ratio = lower_bound / twoe - 1.0L;
            if (ratio > 0.0L) {
                long double val = std::ceil(ratio * SigCountLD);
                if (val >= SigCountLD) {
                    continue;
                }
                if (val < 0.0L) {
                    val = 0.0L;
                }
                smin = static_cast<uint32_t>(val);
            }
        }

        uint32_t smax = SigCount - 1u;
        if (e == emax) {
            long double ratio = upper_bound / twoe - 1.0L;
            long double val = std::floor(ratio * SigCountLD) - 1.0L;
            if (val < 0.0L) {
                continue;
            }
            if (val > SigCountLD - 1.0L) {
                val = SigCountLD - 1.0L;
            }
            smax = static_cast<uint32_t>(val);
        }

        if (smax < smin) {
            continue;
        }

        strata.push_back(Stratum{e, smin, smax, static_cast<uint32_t>(smax - smin + 1)});
    }

    if (strata.empty()) {
        if (n == 0) {
            return NumericVector(0);
        }
        throw std::runtime_error("No admissible float32 strata within the requested bounds");
    }

    RNGScope scope;

    NumericVector result(n);

    if (!weighted) {
        // Deterministic stratified walk: rotate through the strata while taking a
        // random starting point within the cycle so that the workload distributes
        // evenly without clustering around the lower exponents.
        const std::size_t k = strata.size();
        const std::size_t start = sample_offset(static_cast<uint32_t>(k));
        for (int i = 0; i < n; ++i) {
            const Stratum& st = strata[(start + static_cast<std::size_t>(i)) % k];
            uint32_t range = st.count;
            uint32_t step = sample_offset(range);
            uint32_t significand = st.smin + step;
            uint32_t bits = (static_cast<uint32_t>(st.exponent + 127) << 23) | significand;
            float sample = FromBits(bits);
            result[i] = sample;
        }
    } else {
        // Weighted sampling: strata receive probability mass proportional to the
        // number of distinct significands they contain so that frequent exponents
        // are represented in proportion to their IEEE-754 multiplicity.
        std::vector<uint64_t> prefix(strata.size());
        uint64_t total = 0;
        for (std::size_t i = 0; i < strata.size(); ++i) {
            total += strata[i].count;
            prefix[i] = total;
        }

        for (int i = 0; i < n; ++i) {
            double u = unif_rand();
            double scaled = u * static_cast<double>(total);
            uint64_t target = static_cast<uint64_t>(scaled);
            if (target >= total) {
                target = total - 1;
            }

            auto it = std::lower_bound(prefix.begin(), prefix.end(), target + 1);
            std::size_t idx = static_cast<std::size_t>(it - prefix.begin());
            const Stratum& st = strata[idx];

            uint32_t step = sample_offset(st.count);
            uint32_t significand = st.smin + step;
            uint32_t bits = (static_cast<uint32_t>(st.exponent + 127) << 23) | significand;
            float sample = FromBits(bits);
            result[i] = sample;
        }
    }

    return result;
}

// Parallel reducer that scans the candidate magic constants and tracks the
// constant that minimises the maximum absolute relative error. Each worker
// keeps the best candidate for its partition; the final join retains the
// overall best.
struct MagicReducer : public Worker {
    const RVector<double> floats;
    const RVector<int> magics;
    const int NRmax;

    uint32_t best_magic;
    float min_max_error;
    float best_avg_error;

    MagicReducer(const NumericVector floats,
                 const IntegerVector magics,
                 int NRmax)
        : floats(floats),
          magics(magics),
          NRmax(NRmax),
          best_magic(static_cast<uint32_t>(magics[0])),
          min_max_error(std::numeric_limits<float>::max()),
          best_avg_error(0.0f) {}

    MagicReducer(const MagicReducer& reducer, Split)
        : floats(reducer.floats),
          magics(reducer.magics),
          NRmax(reducer.NRmax),
          best_magic(static_cast<uint32_t>(reducer.magics[0])),
          min_max_error(std::numeric_limits<float>::max()),
          best_avg_error(0.0f) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            uint32_t magic = static_cast<uint32_t>(magics[i]);
            float total_error = 0.0f;
            float max_error = 0.0f;

            // Evaluate the current candidate across all sampled floats. The
            // reducer keeps track of both the supremum of the error (used as the
            // optimisation objective) and the accumulated totals needed for
            // optional summaries.
            for (int j = 0; j < floats.length(); ++j) {
                float x = floats[j];
                float approx = frsr0(x, magic, NRmax);
                float actual = 1.0f / std::sqrt(x);
                float rel_error = std::abs((approx - actual) / actual);
                total_error += rel_error;
                if (rel_error > max_error) {
                    max_error = rel_error;
                }
            }

            float avg_error = total_error / floats.length();

            if (max_error < min_max_error) {
                min_max_error = max_error;
                best_magic = magic;
                best_avg_error = avg_error;
            }
        }
    }

    void join(const MagicReducer& rhs) {
        if (rhs.min_max_error < min_max_error) {
            min_max_error = rhs.min_max_error;
            best_magic = rhs.best_magic;
            best_avg_error = rhs.best_avg_error;
        }
    }
};

// [[Rcpp::export]]
DataFrame search_optimal_constant(NumericVector floats,
                                 IntegerVector magics,
                                 int NRmax,
                                 bool return_mean_error /*= true*/,
                                 bool return_sum_error /*= false*/,
                                 bool return_rms_error /*= false*/,
                                 bool return_variance_error /*= false*/) {
    if (magics.length() == 0) {
        stop("`magics` must contain at least one candidate");
    }
    if (floats.length() == 0) {
        stop("`floats` must contain at least one sample");
    }

    MagicReducer reducer(floats, magics, NRmax);
    parallelReduce(0, magics.length(), reducer);

    // Decide which expensive aggregates to compute for the winning magic. Each
    // additional scalar requires a full sweep of the float samples, so we honor
    // only the pieces of information the caller explicitly requested.
    RequestedMetrics request{
        /*want_sum=*/return_sum_error,
        /*want_mean=*/return_mean_error,
        /*want_rms=*/return_rms_error,
        /*want_variance=*/return_variance_error};

    MetricAccumulators acc = accumulate_errors(floats, reducer.best_magic, NRmax, request);

    List output;
    output.push_back(static_cast<double>(reducer.best_magic), "Magic");
    output.push_back(static_cast<double>(reducer.min_max_error), "Error_Objective");

    if (return_mean_error) {
        output.push_back(static_cast<double>(reducer.best_avg_error), "Error_Mean");
    }
    if (return_sum_error) {
        output.push_back(acc.sum, "Error_Sum");
    }
    if (return_rms_error) {
        double count = static_cast<double>(floats.length());
        double mean_square = acc.sum_sq / count;
        output.push_back(std::sqrt(mean_square), "Error_RootMeanSquare");
    }
    if (return_variance_error) {
        double count = static_cast<double>(floats.length());
        double mean_error = return_mean_error ? static_cast<double>(reducer.best_avg_error)
                                              : acc.sum / count;
        double mean_square = acc.sum_sq / count;
        double variance = mean_square - mean_error * mean_error;
        if (variance < 0.0 && std::abs(variance) < 1e-15) {
            variance = 0.0;  // Guard minor floating point drift below zero.
        }
        output.push_back(variance, "Error_Variance");
    }

    return DataFrame(output);
}
