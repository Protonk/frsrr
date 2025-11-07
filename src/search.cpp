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

struct MagicReducer : public Worker {
    const RVector<double> floats;
    const RVector<int> magics;
    const int NRmax;

    uint32_t best_magic;
    float min_max_error;
    float best_avg_error;

    MagicReducer(const NumericVector floats, const IntegerVector magics, int NRmax)
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
DataFrame search_optimal_constant(NumericVector floats, IntegerVector magics, int NRmax) {
    if (magics.length() == 0) {
        stop("`magics` must contain at least one candidate");
    }

    MagicReducer reducer(floats, magics, NRmax);
    parallelReduce(0, magics.length(), reducer);

    return DataFrame::create(
        Named("Magic") = reducer.best_magic,
        Named("Avg_Relative_Error") = reducer.best_avg_error,
        Named("Max_Relative_Error") = reducer.min_max_error
    );
}
