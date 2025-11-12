#ifndef FRSR_SAMPLE_H
#define FRSR_SAMPLE_H

#include <Rcpp.h>
#include <string>

Rcpp::NumericVector sample_inputs(int n,
                                  double x_min,
                                  double x_max,
                                  const std::string& method);

#endif  // FRSR_SAMPLE_H
