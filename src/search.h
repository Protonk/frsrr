#ifndef FRSR_SEARCH_H
#define FRSR_SEARCH_H

#include <Rcpp.h>
#include <string>

using namespace Rcpp;

NumericVector bounded_stratified_sample(int n, double low, double high, bool weighted = false);
DataFrame search_optimal_constant(NumericVector floats,
                                 IntegerVector magics,
                                 int NRmax = 0,
                                 std::string objective_metric = "max",
                                 CharacterVector metrics_to_report = CharacterVector::create("mean"));

#endif // FRSR_SEARCH_H
