#ifndef FRSR_H
#define FRSR_H

#include <Rcpp.h>
#include <cstdint>

Rcpp::DataFrame frsr(Rcpp::DataFrame input, bool keep_params);
float frsr0(float x, std::uint32_t magic, int NRmax);

#endif
