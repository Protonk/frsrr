#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <limits>

using namespace Rcpp;
using namespace RcppParallel;

struct FRSRResult {
    NumericVector result;
    NumericVector initial;
    NumericVector after_one;
    NumericVector final;
    NumericVector error;
    NumericVector diff;
    IntegerVector iters_vec;
    
    FRSRResult(int n) : result(n), initial(n), after_one(n), final(n), error(n), diff(n), iters_vec(n) {}
};

// [[Rcpp::export]]
NumericVector frsr0(NumericVector x, uint32_t magic = 0x5f3759df) {
    NumericVector initial_guess(x.size());
    for (std::size_t j = 0; j < x.size(); ++j) {
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
        } y = {static_cast<float>(x[j])};
        // The magic constant is largely a restoring constant,
        // restoring the exponent bits lost when the float is right shifted.
        // But a specially chosen constant can give a better first guess.
        // This divides by -2, so now we have a guess at -1/2 * log2(x).
        y.u = magic - (y.u >> 1);
        // Blink and you'll miss it. Moving from y.u to y.f takes us
        // out of the logarithmic domain and approximates an exponential.
        // exp(-1/2 * log2(x)) = 1/sqrt(x), so that's our approximation
        // to feed into Newton's method.
        initial_guess[j] = y.f;
    }
    return initial_guess;
}

struct FRSRWorker : public Worker {
    const RVector<double> x;
    const uint32_t magic;
    const int NR;
    const float A, B, tol;
    FRSRResult& res;
    NumericVector initial_guess;

    FRSRWorker(const NumericVector x, uint32_t magic, int NR, float A, float B, float tol, FRSRResult& res)
        : x(x), magic(magic), NR(NR), A(A), B(B), tol(tol), res(res), initial_guess(frsr0(x, magic)) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t j = begin; j < end; ++j) {
            float x_val = x[j];
            float reference = 1 / std::sqrt(x_val);
            float rel_error = std::numeric_limits<float>::max();
            int actual_iters = 0;
            float y = initial_guess[j];
            float prev_y = y;
            res.initial[j] = y;

            for (int i = 1; i <= NR; i++) {
                actual_iters++;
                y = y * (A - B * x_val * y * y);
                rel_error = std::abs(y - reference) / reference;
                if (i == 1) res.after_one[j] = y;
                if (tol > 0 && rel_error <= tol) break;
            }

            res.diff[j] = y - prev_y;
            prev_y = y;
            res.final[j] = y;
            res.iters_vec[j] = actual_iters;
            res.error[j] = rel_error;
            res.result[j] = y;
        }
    }
};

FRSRResult core_frsr(NumericVector x, IntegerVector magic, IntegerVector NR, NumericVector A, NumericVector B, NumericVector tol) {
    int n = x.size();
    FRSRResult res(n);

    IntegerVector magic_vec = (magic.size() == 1) ? IntegerVector(n, magic[0]) : magic;
    IntegerVector NR_vec = (NR.size() == 1) ? IntegerVector(n, NR[0]) : NR;
    NumericVector A_vec = (A.size() == 1) ? NumericVector(n, A[0]) : A;
    NumericVector B_vec = (B.size() == 1) ? NumericVector(n, B[0]) : B;
    NumericVector tol_vec = (tol.size() == 1) ? NumericVector(n, tol[0]) : tol;

    FRSRWorker frsrWorker(x, magic_vec[0], NR_vec[0], A_vec[0], B_vec[0], tol_vec[0], res);
    parallelFor(0, n, frsrWorker);

    return res;
}

// [[Rcpp::export]]
NumericVector frsr(NumericVector x, IntegerVector magic, IntegerVector NR, NumericVector A, NumericVector B, NumericVector tol) {
    FRSRResult res = core_frsr(x, magic, NR, A, B, tol);
    return res.result;
}

// [[Rcpp::export]]
DataFrame frsr_detail(NumericVector x, IntegerVector magic, IntegerVector NR, NumericVector A, NumericVector B, NumericVector tol, bool keep_params = false) {
    FRSRResult res = core_frsr(x, magic, NR, A, B, tol);
    if (keep_params) {
        return DataFrame::create(
            _["input"] = x,
            _["initial"] = res.initial,
            _["after_one"] = res.after_one,
            _["final"] = res.final,
            _["error"] = res.error,
            _["diff"] = res.diff,
            _["iters"] = res.iters_vec,
            _["magic"] = magic,
            _["NR"] = NR,
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
            _["iters"] = res.iters_vec
        );
    }
}
