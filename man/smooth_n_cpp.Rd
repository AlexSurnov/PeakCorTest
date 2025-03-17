\name{smooth_n_cpp}
\alias{smooth_n_cpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Iterative smoothing
}
\description{
A C++ code that averages adjacent elements of the input vector and runs this procedure iteratively the indicated number of times.
Here is the code.

#include <Rcpp.h>

#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]

std::vector<double> smooth_n_cpp(const std::vector<double>& x, int& n) \{

  std::vector<double> y(x.size());
  
  std::vector<double> z(x.size());
  
  z = x;
  
  if (n == 0) \{
  
    return z;
    
  \} else \{
  
    for(int k=1; k <= n; k++)\{
    
      y = z;
      
      for(int j=1; j < (z.size()-1); j++)\{
      
        y[j] = (z[j-1]+z[j]+z[j+1])/3;
        
      \}
      
      y[0] = (z[1]+z[0])/2;
      
      y[z.size()-1] = (z[z.size()-2]+z[z.size()-1])/2;
      
      z = y;
      
    \}
    
    return y;
    
  \}
  
\}

}

\usage{
smooth_n_cpp(x, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric: a vector
}
  \item{n}{
integer: degree of smoothing (nonnegative)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
numeric: if `n=0`, the output is `x`. If `n=1`, the output is smoothed `S1(x[n])=mean(x[(n-1):(n+1)])` for `n in 2:(length(x)-1)`, and `S1(x[n])` equals to the mean of the two adjacent marginal elements for `n in c(1,length(x))`. If `n>1`, the return is the output of `S1` iterated `n` times.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
smooth_n_cpp(x = c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1), n = 3)
lapply(X = as.list(0:21),
FUN = function(n) {
smooth_n_cpp(x = c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1), n = n)
})
}