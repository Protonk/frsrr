#include <Rcpp.h>
#include <cmath>
#include <cstdint>
#include <RcppParallel.h>
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
        y.u = static_cast<uint32_t>(magic) - (y.u >> 1);
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
    : x(x), magic(magic), NR(NR), A(A), B(B), tol(tol), res(res), initial_guess(frsr0(NumericVector(x), magic)) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; ++j) {
      float x_val = x[j];
      float reference = 1 / std::sqrt(x_val);
      float rel_error = std::numeric_limits<float>::max();
      int actual_iters = 0;
      // Initial guess comes from frsr0()
      float y = initial_guess[j];
      res.initial[j] = y;
      // so we always have a diff available (may be 0 if NR == 0)
      float prev_y = y;

      for (int i = 1; i <= NR; i++) {
        actual_iters++;
        y = y * (A - B * x_val * y * y);
        rel_error = std::abs(y - reference) / reference;
        if (i == 1) res.after_one[j] = y;
        // If we're within tolerance, we can stop early
        if (tol > 0 && rel_error <= tol) break;
      }
      // Rather than record all iterations, we get the difference between the
      // final and the previous iteration. 
      res.diff[j] = y - prev_y;
      res.final[j] = y;
      res.iters_vec[j] = actual_iters;
      // Error is more useful to return than the reference value
      res.error[j] = rel_error;
      res.result[j] = y;
    }
  }
};

FRSRResult core_frsr(NumericVector x, uint32_t magic, int NR, float A, float B, float tol) {
  int n = x.size();
  FRSRResult res(n);
  FRSRWorker frsrWorker(x, magic, NR, A, B, tol, res);
  parallelFor(0, n, frsrWorker);
  return res;
}

// [[Rcpp::export]]
NumericVector frsr(NumericVector x, uint32_t magic = 0x5f3759df, int NR = 1, float A = 1.5, float B = 0.5, float tol = 0) {
  FRSRResult res = core_frsr(x, magic, NR, A, B, tol);
  return res.result;
}

// [[Rcpp::export]]
DataFrame frsr_detail(NumericVector x, uint32_t magic = 0x5f3759df, int NR = 1, float A = 1.5, float B = 0.5, float tol = 0) {
  FRSRResult res = core_frsr(x, magic, NR, A, B, tol);
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