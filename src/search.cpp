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

namespace {
inline float FromBits(uint32_t bits) {
    return std::bit_cast<float>(bits);
}

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
}

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

    const long double lower_bound = static_cast<long double>(std::exp2(low));
    const long double upper_bound = static_cast<long double>(std::exp2(high));
    // Use long double so the derived significand bounds hold their precision
    // even when callers pass fractional log2 endpoints.

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
        const std::size_t k = strata.size();
        const std::size_t start = sample_offset(static_cast<uint32_t>(k));
        // Walking the strata with a rotating start index keeps the marginal
        // exponent distribution uniform without shuffling the container each time.
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
        std::vector<uint64_t> prefix(strata.size());
        uint64_t total = 0;
        for (std::size_t i = 0; i < strata.size(); ++i) {
            total += strata[i].count;
            prefix[i] = total;
        }
        // Weighted sampling approximates the natural float density where wider
        // exponent ranges contain more representable significands.

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

namespace {

constexpr const char* kMetricMaxRelativeError = "max_relative_error";
constexpr const char* kMetricAvgRelativeError = "avg_relative_error";
constexpr const char* kMetricRmseRelativeError = "rmse_relative_error";

struct MetricResult {
    float max_relative_error;
    float avg_relative_error;
    float rmse_relative_error;
};

inline std::string metric_options() {
    return std::string(kMetricMaxRelativeError) + ", " + kMetricAvgRelativeError + ", " +
           kMetricRmseRelativeError;
}

inline void validate_metric_name(const std::string& metric_name, const char* parameter_label) {
    if (metric_name == kMetricMaxRelativeError ||
        metric_name == kMetricAvgRelativeError ||
        metric_name == kMetricRmseRelativeError) {
        return;
    }
    std::string message = "`";
    message += parameter_label;
    message += "` must be one of: ";
    message += metric_options();
    throw std::invalid_argument(message);
}

inline float metric_value(const MetricResult& metrics, const std::string& metric_name) {
    if (metric_name == kMetricMaxRelativeError) {
        return metrics.max_relative_error;
    }
    if (metric_name == kMetricAvgRelativeError) {
        return metrics.avg_relative_error;
    }
    if (metric_name == kMetricRmseRelativeError) {
        return metrics.rmse_relative_error;
    }
    // Should never reach here because validate_metric_name should guard inputs.
    throw std::logic_error("Unexpected metric requested");
}

}

struct MagicReducer : public Worker {
    const RVector<double> floats;
    const RVector<int> magics;
    const int NRmax;
    const std::string objective_metric;
    const std::string dependent_metric;

    uint32_t best_magic;
    float best_objective_value;
    float best_dependent_value;

    MagicReducer(const NumericVector floats,
                 const IntegerVector magics,
                 int NRmax,
                 std::string objective_metric,
                 std::string dependent_metric)
        : floats(floats),
          magics(magics),
          NRmax(NRmax),
          objective_metric(std::move(objective_metric)),
          dependent_metric(std::move(dependent_metric)),
          best_magic(static_cast<uint32_t>(magics[0])),
          best_objective_value(std::numeric_limits<float>::max()),
          best_dependent_value(0.0f) {}

    MagicReducer(const MagicReducer& reducer, Split)
        : floats(reducer.floats),
          magics(reducer.magics),
          NRmax(reducer.NRmax),
          objective_metric(reducer.objective_metric),
          dependent_metric(reducer.dependent_metric),
          best_magic(static_cast<uint32_t>(reducer.magics[0])),
          best_objective_value(std::numeric_limits<float>::max()),
          best_dependent_value(0.0f) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            uint32_t magic = static_cast<uint32_t>(magics[i]);
            float total_error = 0.0f;
            float max_error = 0.0f;
            float squared_error = 0.0f;

            for (int j = 0; j < floats.length(); ++j) {
                // floats is an RVector wrapper, so reads here are zero-copy even
                // though this loop might be running on multiple worker threads.
                float x = floats[j];
                float approx = frsr0(x, magic, NRmax);
                float actual = 1.0f / std::sqrt(x);
                float rel_error = std::abs((approx - actual) / actual);
                total_error += rel_error;
                squared_error += rel_error * rel_error;
                if (rel_error > max_error) {
                    max_error = rel_error;
                }
            }

            MetricResult metrics{
                max_error,
                total_error / floats.length(),
                std::sqrt(squared_error / floats.length())
            };
            float objective_value = metric_value(metrics, objective_metric);
            if (objective_value < best_objective_value) {
                best_objective_value = objective_value;
                best_dependent_value = metric_value(metrics, dependent_metric);
                best_magic = magic;
            }
        }
    }

    void join(const MagicReducer& rhs) {
        if (rhs.best_objective_value < best_objective_value) {
            best_objective_value = rhs.best_objective_value;
            best_magic = rhs.best_magic;
            best_dependent_value = rhs.best_dependent_value;
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

    validate_metric_name(objective_metric, "objective");
    validate_metric_name(dependent_metric, "dependent");

    MagicReducer reducer(floats, magics, NRmax, objective_metric, dependent_metric);
    parallelReduce(0, magics.length(), reducer);

    return DataFrame::create(
        Named("Magic") = reducer.best_magic,
        Named("Dependent") = reducer.best_dependent_value,
        Named("Objective") = reducer.best_objective_value
    );
}
