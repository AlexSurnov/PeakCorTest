#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
// [[Rcpp::export]]
std::vector<double> smooth_n_cpp(const std::vector<double>& x, int& n) {
  std::vector<double> y(x.size());
  std::vector<double> z(x.size());
  z = x;
  if (n == 0) {
    return z;
  } else {
    for(int k=1; k <= n; k++){
      y = z;
      for(int j=1; j < (z.size()-1); j++){
        y[j] = (z[j-1]+z[j]+z[j+1])/3;
      }
      y[0] = (z[1]+z[0])/2;
      y[z.size()-1] = (z[z.size()-2]+z[z.size()-1])/2;
      z = y;
    }
    return y;
  }
}