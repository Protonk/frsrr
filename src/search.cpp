#include "search.h"

#include <RcppParallel.h>
#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "frsr.h"

using namespace Rcpp;
using namespace RcppParallel;

namespace detail {

// ---- Float Sampling Helpers -------------------------------------------------
// We turn the open-coded bit tricks into helpers so the exported functions read like intent.

inline float FromBits(uint32_t bits) {
    return std::bit_cast<float>(bits);
}

inline uint32_t SampleOffset(uint32_t range) {
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

// Each stratum represents all float32 values that share an exponent and live
// between two significand indices.
struct Stratum {
    int exponent;
    uint32_t smin;
    uint32_t smax;
    uint32_t count;
};

// Convert the caller's log2 bounds into a list of usable strata.
inline std::vector<Stratum> BuildStrata(double low, double high) {
    const long double lower_bound = static_cast<long double>(std::exp2(low));
    const long double upper_bound = static_cast<long double>(std::exp2(high));
    // Long-double keeps the mantissa bounds accurate when low/high have fractions.

    constexpr uint32_t kSigCount = 1u << 23;
    const long double SigCountLD = static_cast<long double>(kSigCount);
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

        uint32_t smax = kSigCount - 1u;
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

    return strata;
}

// Draw a single random float from a stratum by picking a significand.
inline float DrawSample(const Stratum& st) {
    uint32_t step = SampleOffset(st.count);
    uint32_t significand = st.smin + step;
    uint32_t bits = (static_cast<uint32_t>(st.exponent + 127) << 23) | significand;
    return FromBits(bits);
}

// ---- Metric Helpers ---------------------------------------------------------
// These helpers keep the metric bookkeeping in one spot so the reducer logic stays short.

constexpr const char* kMetricMaxRelativeError = "max_relative_error";
constexpr const char* kMetricAvgRelativeError = "avg_relative_error";
constexpr const char* kMetricRmseRelativeError = "rmse_relative_error";

struct MetricResult {
    float max_relative_error;
    float avg_relative_error;
    float rmse_relative_error;
};

struct ErrorAccumulator {
    float total_error = 0.0f;
    float squared_error = 0.0f;
    float max_error = 0.0f;
};

inline std::string MetricOptions() {
    return std::string(kMetricMaxRelativeError) + ", " + kMetricAvgRelativeError + ", " +
           kMetricRmseRelativeError;
}

inline void ValidateMetricName(const std::string& metric_name, const char* label) {
    if (metric_name == kMetricMaxRelativeError ||
        metric_name == kMetricAvgRelativeError ||
        metric_name == kMetricRmseRelativeError) {
        return;
    }
    std::string message = "`";
    message += label;
    message += "` must be one of: ";
    message += MetricOptions();
    throw std::invalid_argument(message);
}

inline float MetricValue(const MetricResult& metrics, const std::string& metric_name) {
    if (metric_name == kMetricMaxRelativeError) {
        return metrics.max_relative_error;
    }
    if (metric_name == kMetricAvgRelativeError) {
        return metrics.avg_relative_error;
    }
    if (metric_name == kMetricRmseRelativeError) {
        return metrics.rmse_relative_error;
    }
    throw std::logic_error("Unexpected metric requested");
}

}  // namespace detail

// -----------------------------------------------------------------------------
// bounded_stratified_sample
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
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

    // Step 1: carve the float32 line into exponent "strata" that respect the bounds.
    std::vector<detail::Stratum> strata = detail::BuildStrata(low, high);
    if (strata.empty()) {
        if (n == 0) {
            return NumericVector(0);
        }
        throw std::runtime_error("No admissible float32 strata within the requested bounds");
    }

    RNGScope scope;  // Make sure R's RNG state is saved/restored while we draw.
    NumericVector result(n);

    if (!weighted) {
        // Choose a random offset so repeated calls do not always start at exponent emin.
        const std::size_t k = strata.size();
        const std::size_t start = detail::SampleOffset(static_cast<uint32_t>(k));
        for (int i = 0; i < n; ++i) {
            const detail::Stratum& st = strata[(start + static_cast<std::size_t>(i)) % k];
            result[i] = detail::DrawSample(st);
        }
        return result;
    }

    // Weighted sampling mimics float density: wide exponent ranges deserve more draws.
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
        result[i] = detail::DrawSample(strata[idx]);
    }
    return result;
}

// -----------------------------------------------------------------------------
// MagicReducer
// -----------------------------------------------------------------------------

// Each worker thread aggregates error statistics for every candidate over a slice of the samples.
struct MagicReducer : public Worker {
    const float* samples;
    const float* rsqrt_values;
    const uint32_t* magics;
    const std::size_t magic_count;
    const int NRmax;

    std::vector<detail::ErrorAccumulator> accumulators;

    MagicReducer(const float* samples,
                 const float* rsqrt_values,
                 const uint32_t* magics,
                 std::size_t magic_count,
                 int NRmax)
        : samples(samples),
          rsqrt_values(rsqrt_values),
          magics(magics),
          magic_count(magic_count),
          NRmax(NRmax),
          accumulators(magic_count) {}

    MagicReducer(const MagicReducer& reducer, Split)
        : samples(reducer.samples),
          rsqrt_values(reducer.rsqrt_values),
          magics(reducer.magics),
          magic_count(reducer.magic_count),
          NRmax(reducer.NRmax),
          accumulators(magic_count) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t sample_idx = begin; sample_idx < end; ++sample_idx) {
            const float sample = samples[sample_idx];
            const float actual = rsqrt_values[sample_idx];

            for (std::size_t magic_idx = 0; magic_idx < magic_count; ++magic_idx) {
                const uint32_t magic = magics[magic_idx];
                const float approx = frsr0(sample, magic, NRmax);
                const float rel_error = std::abs((approx - actual) / actual);
                auto& acc = accumulators[magic_idx];
                acc.total_error += rel_error;
                acc.squared_error += rel_error * rel_error;
                if (rel_error > acc.max_error) {
                    acc.max_error = rel_error;
                }
            }
        }
    }

    void join(const MagicReducer& rhs) {
        for (std::size_t magic_idx = 0; magic_idx < magic_count; ++magic_idx) {
            accumulators[magic_idx].total_error += rhs.accumulators[magic_idx].total_error;
            accumulators[magic_idx].squared_error += rhs.accumulators[magic_idx].squared_error;
            if (rhs.accumulators[magic_idx].max_error > accumulators[magic_idx].max_error) {
                accumulators[magic_idx].max_error = rhs.accumulators[magic_idx].max_error;
            }
        }
    }
};

// [[Rcpp::export]]
DataFrame search_optimal_constant(NumericVector floats,
                                  IntegerVector magics,
                                  int NRmax,
                                  std::string objective_metric,
                                  std::string dependent_metric) {
    if (magics.length() == 0) {
        stop("`magics` must contain at least one candidate");
    }
    if (floats.length() == 0) {
        stop("`floats` must contain at least one value");
    }

    detail::ValidateMetricName(objective_metric, "objective");
    detail::ValidateMetricName(dependent_metric, "dependent");

    const std::size_t sample_count = static_cast<std::size_t>(floats.length());
    std::vector<float> sample_values(sample_count);
    std::vector<float> rsqrt_values(sample_count);
    for (std::size_t i = 0; i < sample_count; ++i) {
        const float value = static_cast<float>(floats[static_cast<R_xlen_t>(i)]);
        sample_values[i] = value;
        rsqrt_values[i] = 1.0f / std::sqrt(value);
    }

    const std::size_t magic_count = static_cast<std::size_t>(magics.length());
    std::vector<uint32_t> magic_values(magic_count);
    for (std::size_t i = 0; i < magic_count; ++i) {
        magic_values[i] = static_cast<uint32_t>(magics[static_cast<R_xlen_t>(i)]);
    }

    MagicReducer reducer(sample_values.data(),
                         rsqrt_values.data(),
                         magic_values.data(),
                         magic_count,
                         NRmax);
    // parallelReduce shards the float workload so even a small candidate set can saturate threads.
    parallelReduce(static_cast<std::size_t>(0), sample_count, reducer);

    uint32_t best_magic = magic_values[0];
    float best_objective_value = std::numeric_limits<float>::max();
    float best_dependent_value = 0.0f;
    const float inv_sample_count = 1.0f / static_cast<float>(sample_count);

    for (std::size_t i = 0; i < magic_count; ++i) {
        const detail::ErrorAccumulator& acc = reducer.accumulators[i];
        const detail::MetricResult metrics{
            acc.max_error,
            acc.total_error * inv_sample_count,
            std::sqrt(acc.squared_error * inv_sample_count)
        };
        const float objective_value = detail::MetricValue(metrics, objective_metric);
        if (objective_value < best_objective_value) {
            best_objective_value = objective_value;
            best_dependent_value = detail::MetricValue(metrics, dependent_metric);
            best_magic = magic_values[i];
        }
    }

    return DataFrame::create(
        Named("Magic") = best_magic,
        Named("Dependent") = best_dependent_value,
        Named("Objective") = best_objective_value
    );
}
