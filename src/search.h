#ifndef FRSR_SEARCH_H
#define FRSR_SEARCH_H

#include <Rcpp.h>

using namespace Rcpp;

NumericVector bounded_stratified_sample(int n, double low, double high, bool weighted = false);
DataFrame search_optimal_constant(NumericVector floats, IntegerVector magics, int NRmax = 0);

#endif // FRSR_SEARCH_H
