#include <cmath>
#include "frsr0.h"

// The fast reciprocal square root documented without
// a certain word. An expurgated version.
float frsr0(float x, uint32_t magic, int NRmax) {
    // A union here allows us to represent the same memory bits
    // as both a float and an unsigned integer, giving us access
    // to the "evil floating point bit level hacking" of 
    // the https://en.wikipedia.org/wiki/Fast_inverse_square_root
    // without invoking undefined behavior.
    // " [T]]he bit pattern of a floating point number, 
    //   interpreted as an integer, gives a piecewise linear 
    //   approximation to the logarithm function"
    //  From Jim Blinn's "Floating-point tricks" (1997)
    //  The first step to get 1/sqrt(x) is to get -1/2 * log2(x),
    //  this is our "poor man's" logarithm.
    union {
        float f;
        uint32_t u;
    } y = {static_cast<float>(x)};
    // The magic constant is largely a restoring constant,
    // restoring the exponent bits lost when the float is right shifted.
    // But a specially chosen constant can give a better first guess.
    // This divides by -2, so now we have a guess at -1/2 * log2(x).
    y.u = magic - (y.u >> 1);
    // Blink and you'll miss it. Moving from y.u to y.f takes us
    // out of the logarithmic domain and approximates an exponential.
    // exp(-1/2 * log2(x)) = 1/sqrt(x), so that's our approximation
    // to feed into Newton's method.
    float guess = y.f;
    for (int i = 0; i <= NRmax; i++) {
        guess = guess * (1.5f - 0.5f * x * guess * guess);
    }
    return guess;
}