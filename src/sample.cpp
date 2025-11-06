#include <Rcpp.h>
#include <cmath>
#include <random>
#include <stdexcept>
#include <algorithm>
#include <bit>
#include <limits>

using namespace Rcpp;

// Utility functions
inline uint32_t ToBits(float f) { return std::bit_cast<uint32_t>(f); }
inline float FromBits(uint32_t b) { return std::bit_cast<float>(b); }
inline int Exponent(float f) { return ((ToBits(f) >> 23) & 0xff) - 127; }
constexpr int SignificandMask = (1 << 23) - 1;

// Accept log2 exponent bounds even though they arrive as doubles. Allowing
// fractional exponents gives callers partially open intervals when mapping back
// to float space.
// [[Rcpp::export]]
NumericVector boundedStratifiedSample(int n, double low, double high) {

    if (low < -126) {
        throw std::invalid_argument("Subnormal numbers are not supported. 'low' must be >= -126");
    }
    
    if (!(low < high)) {
        throw std::invalid_argument("'low' must be strictly less than 'high'");
    }

    NumericVector result(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> exponent_draw(1, std::numeric_limits<uint32_t>::max());
    std::uniform_int_distribution<uint32_t> significand_draw(0, SignificandMask);

    int emin = std::ceil(low);
    int emax = std::ceil(high);
    int exponent_span = std::max(1, emax - emin);
    
    for (int i = 0; i < n; ++i) {
        float sample;
        // Draw candidate bits and use the leading-zero count to rotate evenly
        // through the exponent range before composing the IEEE-754 float.
        do {
            int e = emax - 1;
            uint32_t bits = exponent_draw(gen);
            int lz = std::countl_zero(bits);
            e -= (lz % exponent_span);

            uint32_t significand = significand_draw(gen);
            sample = FromBits((uint32_t(e + 127) << 23) | significand);
        } while (sample < std::pow(2.0f, low) || sample >= std::pow(2.0f, high));

        result[i] = sample;
        // NumericVector stores doubles; assigning the float promotes it back to
        // R's native type while letting us keep the inner loop in 32-bit math.
    }

    return result;
}