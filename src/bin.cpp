#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "frsr.h"

using namespace Rcpp;
using namespace RcppParallel;

// Worker for parallel error computation
struct ErrorCalculator : public Worker {
    const RVector<double> floats;
    const RVector<int> magics;
    std::vector<double>& errors;
    const int NRmax; 

    ErrorCalculator(const NumericVector& floats, const IntegerVector& magics, std::vector<double>& errors, const int NRmax)
        : floats(floats), magics(magics), errors(errors), NRmax(NRmax) {}

    void operator()(std::size_t begin, std::size_t end) override {
        for (std::size_t i = begin; i < end; ++i) {
            double total_error = 0.0;
            for (std::size_t j = 0; j < floats.size(); ++j) {
                float approx = frsr0(floats[j], magics[i], NRmax);
                float actual = 1.0 / std::sqrt(floats[j]);
                total_error += std::abs((approx - actual) / actual);
            }
            errors[i] = total_error;
        }
    }
};

// Main function to compute errors and find optimal magic constant
// [[Rcpp::export]]
DataFrame computeOptimalMagic(const NumericVector& floats,
                              const IntegerVector& magics,
                              int NRmax) {
    std::vector<double> errors(magics.size(), 0.0);

    // Parallel computation of errors
    ErrorCalculator calculator(floats, magics, errors, NRmax);
    parallelFor(0, magics.size(), calculator);

    // Find the magic constant with the minimum error
    auto min_it = std::min_element(errors.begin(), errors.end());
    int best_magic_index = std::distance(errors.begin(), min_it);
    uint32_t best_magic = magics[best_magic_index];
    double min_error = *min_it;

    return DataFrame::create(
        Named("magic") = best_magic,
        Named("error") = min_error
    );
}
