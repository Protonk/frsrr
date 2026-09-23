#include "sample.h"
#include "frsr.h"

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

inline std::vector<Stratum> BuildStrata(double lower_bound, double upper_bound) {
    constexpr uint32_t kSigCount = 1u << 23;
    std::vector<Stratum> strata;
    // There are only 254 normal exponent strata. Work directly from the original
    // bounds: a log2/exp2 round trip can move an endpoint across an adjacent float.
    for (int e = -126; e <= 127; ++e) {
        const double twoe = std::ldexp(1.0, e);
        const double low = std::max(lower_bound, twoe);
        const double high = std::min(upper_bound, 2.0 * twoe);
        if (low >= high) continue;
        // Float values are integer multiples of 2^(e-23). For [low, high),
        // ceil(high / spacing) - 1 is the last admissible integer, even when
        // high itself is not representable in float32. Power-of-two scaling is exact.
        const auto smin = static_cast<uint32_t>(std::ceil(std::ldexp(low, 23 - e)) - kSigCount);
        const double last = std::ceil(std::ldexp(high, 23 - e)) - kSigCount - 1.0;
        if (last < smin) continue;
        const auto smax = static_cast<uint32_t>(last);
        strata.push_back(Stratum{e, smin, smax, smax - smin + 1u});
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
    if (x_min < std::ldexp(1.0, -126) || x_max > std::ldexp(1.0, 128)) {
        throw std::invalid_argument("Log-stratified bounds must satisfy 2^-126 <= x_min < x_max <= 2^128");
    }
    std::vector<Stratum> strata = BuildStrata(x_min, x_max);
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
        result[i] = std::min(x_min + span * frac, std::nextafter(x_max, x_min));
        current += kAlpha;
    }
    return result;
}

NumericVector Uniform(int n, double x_min, double x_max) {
    RNGScope scope;
    NumericVector result = Rcpp::runif(n, x_min, x_max);
    for (int i = 0; i < n; ++i) {
        result[i] = std::min(result[i], std::nextafter(x_max, x_min));
    }
    return result;
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
