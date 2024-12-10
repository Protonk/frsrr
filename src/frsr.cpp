#include <Rcpp.h>
#include <cmath>
#include <cstdint>
using namespace Rcpp;
//' Fast Reciprocal Square Root (FRSR)
//'
//' A parameterized Fast Reciprocal Square Root algorithm written in C++.
//'
//' @param x A numeric vector of input values.
//' @param magic An unsigned 32-bit integer restoring constant. Default is 0x5f3759df.
//' @param NR An integer specifying the maximum number of Newton-Raphson iterations. Default is 1.
//' @param A Parameter for the Newton-Raphson iteration. Default is \code{1.5}.
//' @param B Parameter for the Newton-Raphson iteration. Default is \code{0.5}.
//' @param tol The absolute relative error at which to stop early. Default is 0 (no early stopping).
//'
//' @return 
//' \code{frsr} returns a numeric vector of the same length as \code{x},
//' containing the \code{final} approximations of 1/sqrt(x).
//'
//' \code{frsr.detail} returns a data frame with columns:
//'     \item{input}{The input values}
//'     \item{magic}{The integer restoring constant chosen}
//'     \item{initial}{Initial approximation from integer operations}
//'     \item{after_one}{Result after one iteration of Newton-Raphson}
//'     \item{final}{Result from final iteration}
//'     \item{error}{Absolute relative error of final versus standard library}
//'     \item{diff}{Difference between final and penultimate approximations}
//'     \item{iters}{Number of iterations performed}
//'
//' @details
//' The function supplies a Fast Reciprocal Square Root algorithm, which provides
//' an approximation of 1/sqrt(x). The user can specify their own parameters. The
//' default values are set to those used by the famous "fast inverse square
//' root" in Quake III Arena.
//'
//' The algorithm exploits the fact that the integer representation of a
//' floating-point number offers a piecewise linear approximation to the
//' logarithm function. Right-shifting the integer bits of a float is
//' equivalent to dividing the logarithm of the number by two.
//' By subtracting this from a carefully chosen constant, an approximation
//' of \eqn{-1/2 * log2(x)} can be obtained. Treating that result as a float
//' again by using integer bits in memory gives a good guess of
//' \eqn{exp(-1/2 * log2(x))}, which is \eqn{1/sqrt(x)}.
//'
//' The "magic" constant principally serves to restore the exponent bits lost
//' when the input float is right shifted. A restoring constant which does only
//' that is `0x5F400000`, given by Blinn 1997. Values of magic from roughly
//' `0x5f2ffb20` to `0x5f404ed0` will give acceptable levels of error.
//'
//' The Newton-Raphson step \eqn{y_{n+1} = y_n * (1.5 - 0.5 * x * y_n^2)}
//' is performed repeatedly until the specified maximum NR is reached. The
//' default is one. Grossly different values of magic from the default may
//' require many iterations to approach the correct output.
//'
//' Parameters in the Newton-Raphson step, \eqn{(A - B * x * y_n^2)} need
//' not be fixed at 1.5 and 0.5 and can be set by the user. Note that
//' if B =/= A - 1, the approximation may fail to converge.
//'
//' @references
//' J. F. Blinn, (July-Aug. 1997) "Floating-point tricks," in IEEE Computer Graphics and Applications, vol. 17, no. 4, pp. 80-84 \doi{10.1109/38.595279}
//'
//' J. T. Coonen, (1984) §2.3 "A Poor Man's Logarithm" in Contributions to a Proposed Standard for Binary Floating-Point Arithmetic. PhD Thesis, University of California Berkeley
//'
//' S. Summit, (2023). Answer to "Why does the integer representation of a floating point number offer a piecewise linear approximation to the logarithm?" Stack Overflow. \url{https://stackoverflow.com/a/75772363/1188479}
//'
//' @examples
//' \dontrun{
//' # Custom parameters
//' result <- frsr(c(1, 4, 9, 16), magic = 0x5f375a86, NR = 2, A = 1.6, B = 0.6)
//' ## result is a vector of length 4
//' print(result)
//' # [1] 0.9990148 0.4995074 0.3337626 0.2497537
//'
//' # Optional detail
//' result.df <- frsr.detail(c(pi, 2^-31, 0.4, 6.02e23))
//' ## result.df is a dataframe with 3 rows and 8 columns
//' print(result.df)
//' #          input      magic      initial    after_one        final        error          diff iters
//' # 1 3.141593e+00 1597463007 5.735160e-01 5.639570e-01 5.639570e-01 0.0004121269 -9.558976e-03     1
//' # 2 4.656613e-10 1597463007 4.693787e+04 4.632937e+04 4.632937e+04 0.0002499308 -6.085039e+02     1
//' # 3 4.000000e-01 1597463007 1.632430e+00 1.578616e+00 1.578616e+00 0.0015955754 -5.381417e-02     1
//' # 4 6.020000e+23 1597463007 1.306493e-12 1.288484e-12 1.288484e-12 0.0002824810 -1.800936e-14     1
//' }
//'
// [[Rcpp::export]]
NumericVector frsr(NumericVector x, uint32_t magic = 0x5f3759df, int NR = 1, float A = 1.5, float B = 0.5, float tol = 0) {
    int n = x.size();
    NumericVector result(n);

    for (int j = 0; j < n; ++j) {
        float x_val = x[j];
        float reference = 1 / std::sqrt(x_val);
        float rel_error = std::numeric_limits<float>::max();

        union {
            float f;
            uint32_t u;
        } y = {x_val};

        y.u = static_cast<uint32_t>(magic) - (y.u >> 1);

        for (int i = 1; i <= NR; i++) {
            y.f = y.f * (A - B * x_val * y.f * y.f);
            rel_error = std::abs(y.f - reference) / reference;

            if (tol > 0 && rel_error <= tol && i >= 1) {
                break;
            }
        }
        
        result[j] = y.f;
    }

    return result;
}
// [[Rcpp::export]]
DataFrame frsr_detail(NumericVector x, uint32_t magic = 0x5f3759df, int NR = 1, float A = 1.5, float B = 0.5, float tol = 0) {
  int n = x.size();
  NumericVector input(n), initial(n), after_one(n), final(n), error(n), diff(n);
  IntegerVector magic_vec(n, magic), iters_vec(n);

  for (int j = 0; j < n; ++j) {
    float x_val = x[j];
    int actual_iters = 0;
    float reference = 1 / std::sqrt(x_val);
    float rel_error = std::numeric_limits<float>::max();
    // Access bits of the float and treat them as bits in memory of an unsigned int
    // A union is used to reinterpret the float bits as an integer,
    // allowing for bit manipulation without violating strict aliasing rules.
    union {
      float f;
      uint32_t u;
    } y = {x_val};
    // Right shift the integer, which is the equivalent of log2(x)/2
    // See https://stackoverflow.com/q/75772363/1188479
    // subtract from a restoring constant to return exponent
    // bits lost in the shift & make that -log2(x)/2
    y.u = static_cast<uint32_t>(magic) - (y.u >> 1);
    // back as a float, this is exp(-log2(x)/2) ~ 1/sqrt(x)
    initial[j] = y.f;

    // this is to track lagged differences
    float prev_y = y.f;

    for (int i = 1; i <= NR; i++) {
      // tracking this rather than i lets us not worry about
      // where our break statements are
      actual_iters++;
      y.f = y.f * (A - B * x_val * y.f * y.f);
      rel_error = std::abs(y.f - reference) / reference;
      // capturing the first and last iteration can be handy
      // without needing to store all of them
      if (i == 1) {
        after_one[j] = y.f;
      }

      // A user may want to iterate to a tolerance
      // or a fixed number of iterations
      if (tol > 0 && rel_error <= tol && i >= 1) {
        break;
      }
    }
    diff[j] = y.f - prev_y;
    final[j] = y.f;
    error[j] = rel_error;
    iters_vec[j] = actual_iters;
  }

  return DataFrame::create(
    _["input"] = x,
    _["magic"] = magic_vec,
    _["initial"] = initial,
    _["after_one"] = after_one,
    _["final"] = final,
    _["error"] = error,
    _["diff"] = diff,
    _["iters"] = iters_vec
  );
}