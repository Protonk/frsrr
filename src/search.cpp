#include "search.h"

#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "frsr.h"

using namespace Rcpp;
using namespace RcppParallel;

namespace search_detail {
enum class Metric { maximum, average, rmse };

Metric ParseMetric(const std::string& name) {
    if (name == "max_relative_error") return Metric::maximum;
    if (name == "avg_relative_error") return Metric::average;
    if (name == "rmse_relative_error") return Metric::rmse;
    stop("Metric must be one of: max_relative_error, avg_relative_error, rmse_relative_error");
}

struct Accumulator {
    double sum = 0.0, squares = 0.0, maximum = 0.0;
    void add(double error) {
        sum += error;
        squares += error * error;
        maximum = std::max(maximum, error);
    }
    void join(const Accumulator& rhs) {
        sum += rhs.sum;
        squares += rhs.squares;
        maximum = std::max(maximum, rhs.maximum);
    }
    double value(Metric metric, std::size_t n) const {
        if (metric == Metric::maximum) return maximum;
        if (metric == Metric::average) return sum / static_cast<double>(n);
        return std::sqrt(squares / static_cast<double>(n));
    }
};

struct Blocks : public Worker {
    const std::vector<float>& samples;
    const std::vector<double>& references;
    const std::vector<int>& magics;
    const int NRmax;
    const std::size_t block_size;
    std::vector<Accumulator>& partials;

    Blocks(const std::vector<float>& samples, const std::vector<double>& references,
           const std::vector<int>& magics, int NRmax, std::size_t block_size,
           std::vector<Accumulator>& partials)
        : samples(samples), references(references), magics(magics), NRmax(NRmax),
          block_size(block_size), partials(partials) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t task = begin; task < end; ++task) {
            const std::size_t block = task / magics.size();
            const std::size_t m = task % magics.size();
            auto& acc = partials[task];
            const std::size_t stop = std::min(samples.size(), (block + 1) * block_size);
            for (std::size_t i = block * block_size; i < stop; ++i) {
                const float y = frsr0(samples[i], static_cast<std::uint32_t>(magics[m]), NRmax);
                acc.add(frsr_detail::RelativeError(y, references[i]));
            }
        }
    }
};
}  // namespace search_detail

// [[Rcpp::export]]
DataFrame search_optimal_constant(NumericVector floats,
                                  IntegerVector magics,
                                  int NRmax,
                                  std::string objective_metric,
                                  std::string dependent_metric) {
    RNGScope scope;
    if (magics.size() == 0) stop("`magics` must contain at least one candidate");
    if (floats.size() == 0) stop("`floats` must contain at least one value");
    if (NRmax == NA_INTEGER || NRmax < 0) stop("`NRmax` must be non-negative");
    const auto objective = search_detail::ParseMetric(objective_metric);
    const auto dependent = search_detail::ParseMetric(dependent_metric);
    const std::size_t n = floats.size(), m = magics.size();
    std::vector<float> samples(n);
    std::vector<double> references(n);
    for (std::size_t i = 0; i < n; ++i) {
        samples[i] = frsr_detail::CheckedInput(floats[i]);
        references[i] = frsr_detail::Reference(samples[i]);
    }
    std::vector<int> candidates(magics.begin(), magics.end());
    for (int magic : candidates) if (magic == NA_INTEGER) stop("`magics` cannot contain NA");

    // Fixed input-dependent blocks, unrelated to scheduler/thread count. A target
    // of 1024 samples amortizes scheduling; a cap of 64 bounds scratch
    // storage to 64 accumulators per candidate. Join strictly in block order.
    const std::size_t blocks = std::min<std::size_t>(64, (n + 1023) / 1024);
    const std::size_t block_size = (n + blocks - 1) / blocks;
    if (m > std::vector<search_detail::Accumulator>().max_size() / blocks) {
        stop("Too many candidate accumulators");
    }
    std::vector<search_detail::Accumulator> partials(blocks * m);
    search_detail::Blocks worker(samples, references, candidates, NRmax, block_size, partials);
    // Independent (block, candidate) tasks also parallelize small sample sets
    // with many candidates, without changing any candidate's summation order.
    parallelFor(0, partials.size(), worker);

    bool have_best = false;
    int best_magic = 0;
    double best_objective = std::numeric_limits<double>::infinity();
    double best_dependent = best_objective;
    for (std::size_t k = 0; k < m; ++k) {
        search_detail::Accumulator acc;
        for (std::size_t b = 0; b < blocks; ++b) acc.join(partials[b * m + k]);
        // Nonfinite approximations invalidate this candidate, regardless of metric.
        if (!std::isfinite(acc.maximum)) continue;
        const double value = acc.value(objective, n);
        if (!have_best || value < best_objective ||
            (value == best_objective && candidates[k] < best_magic)) {
            have_best = true;
            best_magic = candidates[k];
            best_objective = value;
            best_dependent = acc.value(dependent, n);
        }
    }
    if (!have_best) stop("No candidate has finite errors on every sampled input");
    return DataFrame::create(Named("Magic") = best_magic,
                             Named("Dependent") = best_dependent,
                             Named("Objective") = best_objective);
}
