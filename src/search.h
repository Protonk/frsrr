#ifndef FRSR_SEARCH_H
#define FRSR_SEARCH_H

#include <Rcpp.h>
#include <string>

Rcpp::NumericVector bounded_stratified_sample(int n, double low, double high, bool weighted = false);
Rcpp::DataFrame search_optimal_constant(Rcpp::NumericVector floats,
                                        Rcpp::IntegerVector magics,
                                        int NRmax = 0,
                                        std::string objective_metric = "max_relative_error",
                                        std::string dependent_metric = "avg_relative_error");

#endif  // FRSR_SEARCH_H
