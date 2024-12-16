#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <limits>

using namespace Rcpp;
using namespace RcppParallel;

struct FRSRResult {
    NumericVector initial;
    NumericVector after_one;
    NumericVector final;
    NumericVector error;
    NumericVector diff;
    IntegerVector iters;

    FRSRResult(int n) : initial(n), after_one(n), final(n), error(n), diff(n), iters(n) {}
};

float frsr0(float x, uint32_t magic) {
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
            float reference = 1.0 / std::sqrt(x[j]);
            float rel_error = std::numeric_limits<float>::max();
            int actual_iters = 0;

            // Initial guess
            float y = frsr0(x[j], magic[j]);
            res.initial[j] = y;

            if (NRmax[j] == 0) {
                // When NRmax is 0, skip Newton-Raphson iterations
                res.after_one[j] = NA_REAL;
                res.final[j] = y;
                res.diff[j] = NA_REAL;
                res.iters[j] = 0;
                res.error[j] = std::abs(y - reference) / reference;
            } else {
                // Newton-Raphson iterations
                for (int i = 1; i <= NRmax[j]; i++) {
                    actual_iters++;
                    float prev_y = y;
                    y = y * (A[j] - B[j] * x[j] * y * y);
                    // most implementations use 1 iteration only
                    if (i == 1) res.after_one[j] = y;
                    // For iters > 2, a difference is given so
                    // convergence can be checked on first
                    // and last iterations, rather than
                    // storing all of them.
                    res.diff[j] = y - prev_y;
                    rel_error = std::abs(y - reference) / reference;
                    if (tol[j] > 0 && rel_error <= tol[j]) break;
                }
            }

            // Final results
            res.final[j] = y;
            res.iters[j] = actual_iters;
            res.error[j] = rel_error;
        }
    }
};

// [[Rcpp::export]]
DataFrame frsr(DataFrame input, bool keep_params) {
    NumericVector x = input["x"];
    IntegerVector magic = input["magic"];
    IntegerVector NRmax = input["NRmax"];
    NumericVector A = input["A"];
    NumericVector B = input["B"];
    NumericVector tol = input["tol"];
    
    int n = x.size();
    FRSRResult res(n);
    
    FRSRWorker worker(x, magic, NRmax, A, B, tol, res);
    
    parallelFor(0, n, worker);


    if (keep_params) {
        return DataFrame::create(
            _["input"] = x,
            _["initial"] = res.initial,
            _["after_one"] = res.after_one,
            _["final"] = res.final,
            _["error"] = res.error,
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
            _["diff"] = res.diff,
            _["iters"] = res.iters
        );
    }
}





