#include <Rcpp.h>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// Utility functions
inline float FromBits(uint32_t b) { return std::bit_cast<float>(b); }

// Accept log2 exponent bounds even though they arrive as doubles. Allowing
// fractional exponents gives callers partially open intervals when mapping back
// to float space.
// [[Rcpp::export]]
NumericVector boundedStratifiedSample(int n, double low, double high, bool weighted = false) {
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

    auto sample_offset = [&](uint32_t range) -> uint32_t {
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
    };

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