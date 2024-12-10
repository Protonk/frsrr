#include <Rcpp.h>
#include <cmath>
#include <cstdint>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

struct FRSRWorker : public Worker {
  const RVector<double> x;
  const uint32_t magic;
  const int NR;
  const float A, B, tol;
  RVector<double> result;

  FRSRWorker(const NumericVector x, uint32_t magic, int NR, float A, float B, float tol, NumericVector result)
    : x(x), magic(magic), NR(NR), A(A), B(B), tol(tol), result(result) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; ++j) {
      float x_val = x[j];
      float reference = 1 / std::sqrt(x_val);
      float rel_error = std::numeric_limits<float>::max();

      // A union here allows us to represent the same memory bits
      // as both a float and an unsigned integer, giving us access
      // to the "evil floating point bit level hacking" of 
      // the https://en.wikipedia.org/wiki/Fast_inverse_square_root
      // without invoking undefined behavior.
      // "[T]]he bit pattern of a floating point number, 
      //  interpreted as an integer, gives a piecewise linear 
      //  approximation to the logarithm function"
      //  From Jim Blinn's "Floating-point tricks" (1997)
      //  The first step to get 1/sqrt(x) is to get -1/2 * log2(x),
      //  this is our "poor man's" logarithm.
      union {
        float f;
        uint32_t u;
      } y = {x_val};

      // The magic constant is largely a restoring constant,
      // restoring the exponent bits lost when the float is right shifted.
      // But a specially chosen constant can give a better first guess.
      // This divides by -2, so now we have a guess at -1/2 * log2(x).
      y.u = static_cast<uint32_t>(magic) - (y.u >> 1);

      // set NR to 0 to skip iteration entirely
      for (int i = 1; i <= NR; i++) {
        // Blink and you'll miss it. Moving from y.u to y.f takes us
        // out of the logarithmic domain and approximates an exponential.
        // exp(-1/2 * log2(x)) = 1/sqrt(x), so that's our approximation
        // to feed into Newton's method.
        y.f = y.f * (A - B * x_val * y.f * y.f);
        
        // Tolerance of 0 and below are ignored
        if (tol > 0) {
          rel_error = std::abs(y.f - reference) / reference;
          if (rel_error <= tol) {
            break;
          }
        }
      }

      result[j] = y.f;
    }
  }
};

struct FRSRDetailWorker : public Worker {
  const RVector<double> x;
  const uint32_t magic;
  const int NR;
  const float A, B, tol;
  RVector<double> initial, after_one, final, error, diff;
  RVector<int> iters_vec;

  FRSRDetailWorker(const NumericVector x, uint32_t magic, int NR, float A, float B, float tol,
                   NumericVector initial, NumericVector after_one, NumericVector final,
                   NumericVector error, NumericVector diff, IntegerVector iters_vec)
    : x(x), magic(magic), NR(NR), A(A), B(B), tol(tol),
      initial(initial), after_one(after_one), final(final),
      error(error), diff(diff), iters_vec(iters_vec) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; ++j) {
      float x_val = x[j];
      // Tracking rate of convergence (for some values of the magic constant,
      // this could be high)
      int actual_iters = 0;
      float reference = 1 / std::sqrt(x_val);
      float rel_error = std::numeric_limits<float>::max();

      union {
        float f;
        uint32_t u;
      } y = {x_val};

      y.u = static_cast<uint32_t>(magic) - (y.u >> 1);
      initial[j] = y.f;
      // so we always have a diff available (may be 0 if NR == 0)
      float prev_y = y.f;

      for (int i = 1; i <= NR; i++) {
        actual_iters++;
        y.f = y.f * (A - B * x_val * y.f * y.f);
        // Compute relative error here so we can track it
        // vis after the tolerance check for the main method
        rel_error = std::abs(y.f - reference) / reference;

        if (i == 1) {
          // If you take more than 1 iteration to converge and
          // you are modifying NR parameters, the first result
          // can be instructive.
          after_one[j] = y.f;
        }

        if (tol > 0 && rel_error <= tol) {
          break;
        }
      }
      // Rather than record all iterations, we get the final, first,
      // and the difference bertween the final and the output
      // of the last iteration. 
      diff[j] = y.f - prev_y;
      final[j] = y.f;
      // Error is more useful to return than the reference value
      error[j] = rel_error;
      iters_vec[j] = actual_iters;
    }
  }
};
// [[Rcpp::export]]
NumericVector frsr(NumericVector x, uint32_t magic = 0x5f3759df, int NR = 1, float A = 1.5, float B = 0.5, float tol = 0) {
  NumericVector result(x.size());
  FRSRWorker frsrWorker(x, magic, NR, A, B, tol, result);
  parallelFor(0, x.size(), frsrWorker);
  return result;
}
// [[Rcpp::export]]
DataFrame frsr_detail(NumericVector x, uint32_t magic = 0x5f3759df, int NR = 1, float A = 1.5, float B = 0.5, float tol = 0) {
  int n = x.size();
  NumericVector initial(n), after_one(n), final(n), error(n), diff(n);
  IntegerVector magic_vec(n, magic), iters_vec(n);

  FRSRDetailWorker frsrDetailWorker(x, magic, NR, A, B, tol, initial, after_one, final, error, diff, iters_vec);
  parallelFor(0, n, frsrDetailWorker);

  return DataFrame::create(
    _["input"] = x,
    _["magic"] = magic_vec,
    _["initial"] = initial,
    _["after_one"] = after_one,
    _["final"] = final,
    _["error"] = error,
    _["diff"] = diff,
    _["iters"] = iters_vec
  );
}