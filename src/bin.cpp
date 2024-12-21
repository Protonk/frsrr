#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <limits>
#include "frsr.h"

using namespace Rcpp;
using namespace RcppParallel;

struct MagicReducer : public Worker
{
   // input vectors
   const RVector<double> floats;
   const RVector<int> magics;
   const int NRmax;
   
   // accumulated values
   uint32_t best_magic;
   float min_max_error;
   float best_avg_error;
   
   // constructors
   MagicReducer(const NumericVector floats, const IntegerVector magics, int NRmax) 
     : floats(floats), magics(magics), NRmax(NRmax),
       best_magic(magics[0]), min_max_error(std::numeric_limits<float>::max()), best_avg_error(0.0f) {}
   
   MagicReducer(const MagicReducer& reducer, Split) 
     : floats(reducer.floats), magics(reducer.magics), NRmax(reducer.NRmax),
       best_magic(reducer.magics[0]), min_max_error(std::numeric_limits<float>::max()), best_avg_error(0.0f) {}
   
   // process a range of magics
   void operator()(std::size_t begin, std::size_t end) {
     for (std::size_t i = begin; i < end; ++i) {
       uint32_t magic = static_cast<uint32_t>(magics[i]);
       float total_error = 0.0f;
       float max_error = 0.0f;
       
       for (int j = 0; j < floats.length(); ++j) {
         float x = floats[j];
         float approx = frsr0(x, magic, NRmax);
         float actual = 1.0f / std::sqrt(x);
         float rel_error = std::abs((approx - actual) / actual);
         total_error += rel_error;
         if (rel_error > max_error) {
           max_error = rel_error;
         }
       }
       
       float avg_error = total_error / floats.length();
       
       if (max_error < min_max_error) {
         min_max_error = max_error;
         best_magic = magic;
         best_avg_error = avg_error;
       }
     }
   }
   
   // join my results with that of another MagicReducer
   void join(const MagicReducer& rhs) { 
     if (rhs.min_max_error < min_max_error) {
       min_max_error = rhs.min_max_error;
       best_magic = rhs.best_magic;
       best_avg_error = rhs.best_avg_error;
     }
   }
};

// [[Rcpp::export]]
DataFrame optimal_constant_search(NumericVector floats, IntegerVector magics, int NRmax = 0) {
   MagicReducer reducer(floats, magics, NRmax);
   parallelReduce(0, magics.length(), reducer);
   
   return DataFrame::create(
     Named("Magic") = reducer.best_magic,
     Named("Avg_Relative_Error") = reducer.best_avg_error,
     Named("Max_Relative_Error") = reducer.min_max_error
   );
}

