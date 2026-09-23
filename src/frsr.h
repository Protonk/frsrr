#ifndef FRSR_H
#define FRSR_H

#include <Rcpp.h>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>

static_assert(sizeof(float) == sizeof(std::uint32_t) &&
              std::numeric_limits<float>::is_iec559 &&
              std::numeric_limits<float>::radix == 2 &&
              std::numeric_limits<float>::digits == 24 &&
              std::numeric_limits<float>::min_exponent == -125 &&
              std::numeric_limits<float>::max_exponent == 128,
              "frsrr requires IEEE-754 binary32 floats");
static_assert(std::bit_cast<std::uint32_t>(1.0f) == 0x3f800000u,
              "frsrr requires matching float and integer bit ordering");

namespace frsr_detail {
// Called only on the main thread, before workers receive any inputs.
inline float CheckedInput(double x) {
    if (!std::isfinite(x) || x <= 0.0 || x > std::numeric_limits<float>::max()) {
        Rcpp::stop("`x` must convert to a positive normal float32 without overflow");
    }
    const float rounded = static_cast<float>(x);
    if (!std::isnormal(rounded)) {
        Rcpp::stop("`x` must convert to a positive normal float32 without underflow");
    }
    return rounded;
}

inline float Seed(float x, std::uint32_t magic) {
    // Unsigned subtraction intentionally wraps: unusual magics remain explorable.
    const auto bits = magic - (std::bit_cast<std::uint32_t>(x) >> 1);
    return std::bit_cast<float>(bits);
}

inline float Step(float x, float y, float A = 1.5f, float B = 0.5f) {
    // Each store specifies a binary32 rounding point and prevents fused multiply-
    // subtract across stages. Keep left-to-right B*x*y*y: y*y can overflow first.
    // These explicit stores favor a teachable arithmetic contract over throughput.
    volatile float bx = B * x;
    volatile float bxy = bx * y;
    volatile float bxyy = bxy * y;
    volatile float correction = A - bxyy;
    volatile float next = y * correction;
    return next;
}

inline double Reference(float x) {
    return 1.0 / std::sqrt(static_cast<double>(x));
}

inline double RelativeError(float y, double reference) {
    // Invalid exploratory approximations must never look like perfect candidates.
    if (!std::isfinite(y)) return std::numeric_limits<double>::infinity();
    return std::abs((static_cast<double>(y) - reference) / reference);
}
}  // namespace frsr_detail

Rcpp::DataFrame frsr(Rcpp::DataFrame input, bool keep_params);
float frsr0(float x, std::uint32_t magic, int NRmax);

#endif
