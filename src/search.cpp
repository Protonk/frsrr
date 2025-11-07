#include "search.h"

#include <RcppParallel.h>
#include <algorithm>
#include <bit>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
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

enum class MetricId {
    kMax,
    kMean,
    kSum,
    kRms,
    kVariance,
};

std::string normalize_metric_name(const std::string& raw) {
    std::string normalized;
    normalized.reserve(raw.size());
    for (char ch : raw) {
        unsigned char uch = static_cast<unsigned char>(ch);
        if (std::isalnum(uch)) {
            normalized.push_back(static_cast<char>(std::tolower(uch)));
        }
    }
    return normalized;
}

MetricId parse_metric_name(const std::string& raw) {
    std::string normalized = normalize_metric_name(raw);
    if (normalized == "max" || normalized == "maximum") {
        return MetricId::kMax;
    }
    if (normalized == "mean" || normalized == "average") {
        return MetricId::kMean;
    }
    if (normalized == "sum" || normalized == "total") {
        return MetricId::kSum;
    }
    if (normalized == "rms" || normalized == "rootmeansquare") {
        return MetricId::kRms;
    }
    if (normalized == "variance") {
        return MetricId::kVariance;
    }

    stop("Unsupported metric: `" + raw + "`");
}

std::string metric_label(MetricId metric) {
    switch (metric) {
    case MetricId::kMax:
        return "Max";
    case MetricId::kMean:
        return "Mean";
    case MetricId::kSum:
        return "Sum";
    case MetricId::kRms:
        return "RootMeanSquare";
    case MetricId::kVariance:
        return "Variance";
    }
    return "";
}

std::string column_name(MetricId metric, bool objective) {
    std::string prefix = objective ? "Objective_" : "Metric_";
    return prefix + metric_label(metric);
}

struct CandidateStats {
    double sum = 0.0;
    double sum_sq = 0.0;
    double max = 0.0;
};

double evaluate_metric(MetricId metric, const CandidateStats& stats, double sample_count) {
    switch (metric) {
    case MetricId::kMax:
        return stats.max;
    case MetricId::kMean:
        return stats.sum / sample_count;
    case MetricId::kSum:
        return stats.sum;
    case MetricId::kRms: {
        double mean_square = stats.sum_sq / sample_count;
        return std::sqrt(mean_square);
    }
    case MetricId::kVariance: {
        double mean = stats.sum / sample_count;
        double mean_square = stats.sum_sq / sample_count;
        double variance = mean_square - mean * mean;
        if (variance < 0.0 && std::abs(variance) < 1e-15) {
            variance = 0.0;
        }
        return variance;
    }
    }
    stop("Unsupported metric id");
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
// constant that minimises the caller-selected objective metric. Each worker
// keeps the best candidate for its partition; the final join retains the
// overall best summary statistics.
struct MagicReducer : public Worker {
    const RVector<double> floats;
    const RVector<int> magics;
    const int NRmax;
    const MetricId objective;

    uint32_t best_magic;
    CandidateStats best_stats;
    double best_objective;
    bool has_best;

    MagicReducer(const NumericVector floats,
                 const IntegerVector magics,
                 int NRmax,
                 MetricId objective_metric)
        : floats(floats),
          magics(magics),
          NRmax(NRmax),
          objective(objective_metric),
          best_magic(0u),
          best_stats(),
          best_objective(std::numeric_limits<double>::infinity()),
          has_best(false) {}

    MagicReducer(const MagicReducer& reducer, Split)
        : floats(reducer.floats),
          magics(reducer.magics),
          NRmax(reducer.NRmax),
          objective(reducer.objective),
          best_magic(0u),
          best_stats(),
          best_objective(std::numeric_limits<double>::infinity()),
          has_best(false) {}

    void operator()(std::size_t begin, std::size_t end) {
        double sample_count = static_cast<double>(floats.length());
        for (std::size_t i = begin; i < end; ++i) {
            uint32_t magic = static_cast<uint32_t>(magics[i]);
            CandidateStats candidate_stats;

            // Evaluate the current candidate across all sampled floats and
            // accumulate the sufficient statistics required to evaluate any of
            // the scalar metrics supported by the reducer.
            for (int j = 0; j < floats.length(); ++j) {
                float x = floats[j];
                float approx = frsr0(x, magic, NRmax);
                float actual = 1.0f / std::sqrt(x);
                // The per-sample metric captures the absolute relative difference between the
                // approximation and the true reciprocal square root. Downstream aggregations derive
                // various scalar summaries (e.g., RMS, sum) from these measurements.
                double sample_measure = std::abs((approx - actual) / actual);
                candidate_stats.sum += sample_measure;
                candidate_stats.sum_sq += sample_measure * sample_measure;
                if (sample_measure > candidate_stats.max) {
                    candidate_stats.max = sample_measure;
                }
            }

            double objective_value = evaluate_metric(objective, candidate_stats, sample_count);

            if (!has_best || objective_value < best_objective) {
                has_best = true;
                best_objective = objective_value;
                best_magic = magic;
                best_stats = candidate_stats;
            }
        }
    }

    void join(const MagicReducer& rhs) {
        if (!rhs.has_best) {
            return;
        }
        if (!has_best || rhs.best_objective < best_objective) {
            has_best = true;
            best_objective = rhs.best_objective;
            best_magic = rhs.best_magic;
            best_stats = rhs.best_stats;
        }
    }
};

// [[Rcpp::export]]
DataFrame search_optimal_constant(NumericVector floats,
                                 IntegerVector magics,
                                 int NRmax,
                                 std::string objective_metric,
                                 CharacterVector metrics_to_report) {
    if (magics.length() == 0) {
        stop("`magics` must contain at least one candidate");
    }
    if (floats.length() == 0) {
        stop("`floats` must contain at least one sample");
    }

    MetricId objective = parse_metric_name(objective_metric);

    R_xlen_t n_metrics = metrics_to_report.length();
    std::vector<MetricId> metrics;
    metrics.reserve(static_cast<std::size_t>(n_metrics));
    for (R_xlen_t i = 0; i < n_metrics; ++i) {
        if (metrics_to_report[i] == NA_STRING) {
            stop("`metrics_to_report` must not contain missing values");
        }
        std::string raw = as<std::string>(metrics_to_report[i]);
        MetricId parsed = parse_metric_name(raw);
        if (std::find(metrics.begin(), metrics.end(), parsed) == metrics.end()) {
            metrics.push_back(parsed);
        }
    }

    MagicReducer reducer(floats, magics, NRmax, objective);
    parallelReduce(0, magics.length(), reducer);

    if (!reducer.has_best) {
        stop("No admissible magic constant was evaluated");
    }

    List output;
    output.push_back(static_cast<double>(reducer.best_magic), "Magic");
    output.push_back(reducer.best_objective, column_name(objective, /*objective=*/true));

    double sample_count = static_cast<double>(floats.length());
    for (MetricId metric : metrics) {
        double value = evaluate_metric(metric, reducer.best_stats, sample_count);
        output.push_back(value, column_name(metric, /*objective=*/false));
    }

    return DataFrame(output);
}
