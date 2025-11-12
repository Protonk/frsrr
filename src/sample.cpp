#include "sample.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Rcpp;

namespace sample_detail {

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

struct Stratum {
    int exponent;
    uint32_t smin;
    uint32_t smax;
    uint32_t count;
};

inline std::vector<Stratum> BuildStrata(double low_log2, double high_log2) {
    const long double lower_bound = static_cast<long double>(std::exp2(low_log2));
    const long double upper_bound = static_cast<long double>(std::exp2(high_log2));

    constexpr uint32_t kSigCount = 1u << 23;
    const long double SigCountLD = static_cast<long double>(kSigCount);
    const int emin = std::max(-126, static_cast<int>(std::floor(low_log2)));
    const int emax = std::min(127, static_cast<int>(std::ceil(high_log2)) - 1);

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

inline float DrawSample(const Stratum& st) {
    uint32_t step = SampleOffset(st.count);
    uint32_t significand = st.smin + step;
    uint32_t bits = (static_cast<uint32_t>(st.exponent + 127) << 23) | significand;
    return FromBits(bits);
}

NumericVector LogStratified(int n, double x_min, double x_max) {
    if (n == 0) {
        return NumericVector(0);
    }
    if (x_min <= 0.0 || x_max <= 0.0) {
        throw std::invalid_argument("`x_min` and `x_max` must be > 0 for log-stratified sampling");
    }
    const double low = std::log2(x_min);
    const double high = std::log2(x_max);
    if (!std::isfinite(low) || !std::isfinite(high)) {
        throw std::invalid_argument("Log2 bounds must be finite");
    }
    if (high <= low) {
        throw std::invalid_argument("`x_max` must be greater than `x_min`");
    }
    if (low < -126 || high > 128) {
        throw std::invalid_argument("Bounds must satisfy -126 <= log2(x_min) < log2(x_max) <= 128");
    }

    std::vector<Stratum> strata = BuildStrata(low, high);
    if (strata.empty()) {
        throw std::runtime_error("No admissible float32 strata within the requested bounds");
    }

    RNGScope scope;
    NumericVector result(n);

    const std::size_t k = strata.size();
    const std::size_t start = SampleOffset(static_cast<uint32_t>(k));
    for (int i = 0; i < n; ++i) {
        const Stratum& st = strata[(start + static_cast<std::size_t>(i)) % k];
        result[i] = DrawSample(st);
    }
    return result;
}

inline double FractionalPart(double value) {
    double frac = value - std::floor(value);
    if (frac < 0.0) {
        frac += 1.0;
    }
    return frac;
}

NumericVector IrrationalRotation(int n, double x_min, double x_max) {
    NumericVector result(n);
    if (n == 0) {
        return result;
    }
    // (sqrt(5) - 1) / 2 precomputed to keep the constant constexpr without
    // relying on C++23's constexpr sqrt.
    constexpr double kAlpha = 0.6180339887498948482;
    RNGScope scope;
    double start = unif_rand();
    const double span = x_max - x_min;
    double current = start;
    for (int i = 0; i < n; ++i) {
        double frac = FractionalPart(current);
        result[i] = x_min + span * frac;
        current += kAlpha;
    }
    return result;
}

NumericVector Uniform(int n, double x_min, double x_max) {
    RNGScope scope;
    return Rcpp::runif(n, x_min, x_max);
}

}  // namespace sample_detail

// [[Rcpp::export]]
NumericVector sample_inputs(int n,
                            double x_min,
                            double x_max,
                            const std::string& method) {
    if (n < 0) {
        throw std::invalid_argument("`n` must be non-negative");
    }
    if (!std::isfinite(x_min) || !std::isfinite(x_max)) {
        throw std::invalid_argument("`x_min` and `x_max` must be finite");
    }
    if (x_min >= x_max) {
        throw std::invalid_argument("`x_min` must be less than `x_max`");
    }

    if (method == "log_stratified") {
        return sample_detail::LogStratified(n, x_min, x_max);
    }
    if (method == "irrational") {
        return sample_detail::IrrationalRotation(n, x_min, x_max);
    }
    if (method == "uniform") {
        return sample_detail::Uniform(n, x_min, x_max);
    }

    std::string message = "Unknown sampler method: ";
    message += method;
    throw std::invalid_argument(message);
}
