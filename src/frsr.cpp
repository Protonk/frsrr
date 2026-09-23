#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <limits>
#include <vector>
#include "frsr.h"

using namespace Rcpp;
using namespace RcppParallel;

float frsr0(float x, std::uint32_t magic, int NRmax) {
    float y = frsr_detail::Seed(x, magic);
    for (int i = 0; i < NRmax; ++i) y = frsr_detail::Step(x, y);
    return y;
}

struct FRSRWorker : public Worker {
    const std::vector<float>& x;
    const RVector<int> magic, NRmax;
    const RVector<double> A, B, tol;
    RVector<double> initial, after_one, final, error, diff;
    RVector<int> iters;

    FRSRWorker(const std::vector<float>& x, IntegerVector magic, IntegerVector NRmax,
               NumericVector A, NumericVector B, NumericVector tol,
               NumericVector initial, NumericVector after_one, NumericVector final,
               NumericVector error, NumericVector diff, IntegerVector iters)
        : x(x), magic(magic), NRmax(NRmax), A(A), B(B), tol(tol),
          initial(initial), after_one(after_one), final(final), error(error),
          diff(diff), iters(iters) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t j = begin; j < end; ++j) {
            const double reference = frsr_detail::Reference(x[j]);
            float y = frsr_detail::Seed(x[j], static_cast<std::uint32_t>(magic[j]));
            initial[j] = y;
            double rel_error = frsr_detail::RelativeError(y, reference);
            for (int i = 0; i < NRmax[j]; ++i) {
                const float previous = y;
                y = frsr_detail::Step(x[j], y, static_cast<float>(A[j]), static_cast<float>(B[j]));
                if (i == 0) after_one[j] = y;
                // Diagnostic subtraction is double precision, like the error.
                diff[j] = static_cast<double>(y) - previous;
                iters[j] = i + 1;
                rel_error = frsr_detail::RelativeError(y, reference);
                if (tol[j] > 0.0 && rel_error <= tol[j]) break;
            }
            final[j] = y;
            error[j] = rel_error;
        }
    }
};

// [[Rcpp::export]]
DataFrame frsr(DataFrame input, bool keep_params) {
    RNGScope scope;
    NumericVector x = input["x"], A = input["A"], B = input["B"], tol = input["tol"];
    IntegerVector magic = input["magic"], NRmax = input["NRmax"];
    const R_xlen_t n = x.size();
    if (A.size() != n || B.size() != n || tol.size() != n ||
        magic.size() != n || NRmax.size() != n) stop("Parameter columns must have equal lengths");

    std::vector<float> rounded(n);
    for (R_xlen_t j = 0; j < n; ++j) {
        rounded[j] = frsr_detail::CheckedInput(x[j]);
        if (magic[j] == NA_INTEGER) stop("`magic` cannot contain NA");
        if (NRmax[j] == NA_INTEGER || NRmax[j] < 0) stop("`NRmax` must be non-negative");
        if (!std::isfinite(tol[j]) || tol[j] < 0.0) stop("`tol` must be finite and non-negative");
        if (!std::isfinite(A[j]) || !std::isfinite(B[j]) ||
            std::abs(A[j]) > std::numeric_limits<float>::max() ||
            std::abs(B[j]) > std::numeric_limits<float>::max()) {
            stop("`A` and `B` must be finite and within the float32 range");
        }
    }
    // R objects, NA initialization and accessor materialization stay on this thread.
    NumericVector initial(n), after_one(n, NA_REAL), final(n), error(n), diff(n, NA_REAL);
    IntegerVector iters(n);
    FRSRWorker worker(rounded, magic, NRmax, A, B, tol,
                      initial, after_one, final, error, diff, iters);
    parallelFor(0, static_cast<std::size_t>(n), worker);

    DataFrame result = DataFrame::create(
        _["input"] = x, _["initial"] = initial, _["after_one"] = after_one,
        _["final"] = final, _["error"] = error, _["diff"] = diff, _["iters"] = iters);
    if (keep_params) {
        result["magic"] = magic;
        result["NRmax"] = NRmax;
        result["A"] = A;
        result["B"] = B;
        result["tol"] = tol;
    }
    return result;
}
