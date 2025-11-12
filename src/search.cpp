#include "search.h"

#include <RcppParallel.h>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "frsr.h"

using namespace Rcpp;
using namespace RcppParallel;

namespace detail {
// Metric helpers keep bookkeeping concentrated so the reducer logic stays short.

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

    RNGScope scope;  // Keep R's RNG state scoped even though this path is deterministic today.

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
