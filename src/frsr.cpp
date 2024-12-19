#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <limits>
#include "frsr0.h"

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
            float y = frsr0(x[j], magic[j], 0);
            res.initial[j] = y;
            rel_error = std::abs(y - reference) / reference;

            if (NRmax[j] == 0) {
                // When NRmax is 0, skip Newton-Raphson iterations
                res.after_one[j] = NA_REAL;
                res.final[j] = y;
                res.diff[j] = NA_REAL;
                res.iters[j] = 0;
                res.error[j] = rel_error;
            } else {
                // Newton-Raphson iterations
                for (int i = 1; i <= NRmax[j]; i++) {
                    actual_iters++;
                    float prev_y = y;
                    y = y * (A[j] - B[j] * x[j] * y * y);
                    // most implementations use 1 iteration only
                    if (actual_iters == 1) res.after_one[j] = y;
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





