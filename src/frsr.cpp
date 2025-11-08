#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <limits>
#include "frsr.h"

using namespace Rcpp;
using namespace RcppParallel;

struct FRSRResult {
    NumericVector initial;
    NumericVector after_one;
    NumericVector final;
    NumericVector error;
    NumericVector enre;
    NumericVector diff;
    IntegerVector iters;

    FRSRResult(int n)
        : initial(n), after_one(n), final(n), error(n), enre(n), diff(n), iters(n) {}
};

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
    // R numeric vectors hand us doubles; narrowing to float here ensures the
    // bit reinterpretation below sees the 32-bit layout the algorithm expects
    // instead of trying to treat a double's 64-bit payload as a float.
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
    for (int i = 0; i < NRmax; i++) {
        guess = guess * (1.5f - 0.5f * x * guess * guess);
    }
    return guess;
}

struct FRSRWorker : public Worker {
    const NumericVector& x;
    const IntegerVector& magic;
    const IntegerVector& NRmax;
    const NumericVector& A;
    const NumericVector& B;
    const NumericVector& tol;
    
    FRSRResult& res;

    FRSRWorker(const NumericVector& x, const IntegerVector& magic, const IntegerVector& NRmax,
               const NumericVector& A, const NumericVector& B, const NumericVector& tol,
               FRSRResult& res)
        : x(x), magic(magic), NRmax(NRmax), A(A), B(B), tol(tol), res(res) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t j = begin; j < end; ++j) {
            double x_val = x[j];
            double reference = 1.0 / std::sqrt(x_val);
            float rel_error = std::numeric_limits<float>::max();
            int actual_iters = 0;

            int exponent = 0;
            double mantissa = std::frexp(x_val, &exponent);
            double canonical_input = mantissa * 2.0;
            int k = exponent - 1;

            // Evaluate the mantissa-only form so we can report exponent-normalized
            // errors later; this isolates the error the restoring constant causes
            // independent of the exponent bucket the caller sampled from.
            float canonical_y = frsr0(static_cast<float>(canonical_input), magic[j], 0);
            double canonical_reference = 1.0 / std::sqrt(canonical_input);
            float canonical_rel_error = std::abs(canonical_y - canonical_reference) / canonical_reference;

            // Initial guess
            float y = frsr0(static_cast<float>(x_val), magic[j], 0);
            res.initial[j] = y;
            rel_error = std::abs(y - reference) / reference;

            if (NRmax[j] == 0) {
                // When NRmax is 0, skip Newton-Raphson iterations
                // but write conforming results
                res.after_one[j] = NA_REAL;
                res.final[j] = y;
                res.diff[j] = NA_REAL;
                res.iters[j] = 0;
                res.error[j] = rel_error;
            } else {
                // Newton-Raphson iterations for the input value
                for (int i = 1; i <= NRmax[j]; i++) {
                    actual_iters++;
                    float prev_y = y;
                    y = y * (A[j] - B[j] * x_val * y * y);
                    // most implementations use 1 iteration only
                    if (actual_iters == 1) res.after_one[j] = y;
                    // For iters > 2, a difference is given so
                    // convergence can be checked on first
                    // and last iterations, rather than
                    // storing all of them.
                    res.diff[j] = y - prev_y;
                    rel_error = std::abs(y - reference) / reference;
                    // exit early if we are within tolerance &
                    // tolerance argument is set
                    if (tol[j] > 0 && rel_error <= tol[j]) break;
                }

                // Newton-Raphson iterations for the canonical input
                for (int i = 1; i <= NRmax[j]; ++i) {
                    canonical_y = canonical_y * (A[j] - B[j] * canonical_input * canonical_y * canonical_y);
                    canonical_rel_error = std::abs(canonical_y - canonical_reference) / canonical_reference;
                    if (tol[j] > 0 && canonical_rel_error <= tol[j]) break;
                }
            }
            // Final results
            res.final[j] = y;
            res.iters[j] = actual_iters;
            res.error[j] = rel_error;
            double scale = std::pow(2.0, 0.5 * static_cast<double>(k));
            if (!std::isfinite(scale) || scale == 0.0) {
                res.enre[j] = NA_REAL;
            } else {
                // Dividing out the power-of-two contribution lets consumers compare
                // errors across exponents, mirroring the canonical computation above.
                double normalized = static_cast<double>(canonical_y) / scale;
                res.enre[j] = std::abs(normalized - reference) / reference;
            }
        }
    }
};

// [[Rcpp::export]]
DataFrame frsr(DataFrame input, bool keep_params) {
    
    // Extract input parameters
    NumericVector x = input["x"];
    IntegerVector magic = input["magic"];
    IntegerVector NRmax = input["NRmax"];
    NumericVector A = input["A"];
    NumericVector B = input["B"];
    NumericVector tol = input["tol"];
    
    // faster than input.nrow(), and we explicitly
    // tell the user the rows will equal length(x)
    int n = x.size();
    FRSRResult res(n);

    FRSRWorker worker(x, magic, NRmax, A, B, tol, res);
    // The parallel loop is embarrassingly parallel and has no cross-element
    // sharing, so parallelFor can safely stride across inputs without locks.
    parallelFor(0, n, worker);

    if (keep_params) {
        return DataFrame::create(
            _["input"] = x,
            _["initial"] = res.initial,
            _["after_one"] = res.after_one,
            _["final"] = res.final,
            _["error"] = res.error,
            _["enre"] = res.enre,
            _["diff"] = res.diff,
            _["iters"] = res.iters,
            _["magic"] = magic,
            _["NRmax"] = NRmax,
            _["A"] = A,
            _["B"] = B,
            _["tol"] = tol
        );
    } else {
        return DataFrame::create(
            _["input"] = x,
            _["initial"] = res.initial,
            _["after_one"] = res.after_one,
            _["final"] = res.final,
            _["error"] = res.error,
            _["enre"] = res.enre,
            _["diff"] = res.diff,
            _["iters"] = res.iters
        );
    }
}




