#ifndef FRSR_H
#define FRSR_H

#include <cstdint>

using namespace Rcpp;

DataFrame frsr(DataFrame input, bool keep_params);
float frsr0(float x, uint32_t magic, int NRmax);

#endif