#include <Rcpp.h>
#include <cmath>
#include <random>
#include <stdexcept>
#include <algorithm>

using namespace Rcpp;

// Utility functions
inline uint32_t ToBits(float f) { return std::bit_cast<uint32_t>(f); }
inline float FromBits(uint32_t b) { return std::bit_cast<float>(b); }
inline int Exponent(float f) { return ((ToBits(f) >> 23) & 0xff) - 127; }
constexpr int SignificandMask = (1 << 23) - 1;

// Don't be fooled by the argument types, these need to be exponents
// The double type is chosen because they could be fractional.
// [[Rcpp::export]]
NumericVector boundedStratifiedSample(int n, double low, double high) {

    if (low < -126) {
        throw std::invalid_argument("Subnormal numbers are not supported. 'low' must be >= -126");
    }
    
    NumericVector result(n);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> dis(0, std::numeric_limits<uint32_t>::max());
    
    int emin = std::ceil(low);
    int emax = std::ceil(high);
    
    for (int i = 0; i < n; ++i) {
        float sample;
        do {
            int e = emax - 1;
            uint32_t bits = dis(gen);
            int lz = __builtin_clz(bits);
            e -= (lz % (emax - emin));
            
            uint32_t significand = dis(gen) & SignificandMask;
            sample = FromBits((uint32_t(e + 127) << 23) | significand);
        } while (sample < std::pow(2.0f, low) || sample >= std::pow(2.0f, high));
        
        result[i] = sample;
    }
    
    return result;
}