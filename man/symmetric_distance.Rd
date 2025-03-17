\name{symmetric_distance}
\alias{symmetric_distance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Getting the p-value for the Peak Correlation Test
}
\description{
This function generates an upper estimate of the probability that the the two input vectors belong to the neighborhood of each other by chance.
}
\usage{
symmetric_distance(a_x, a_y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a_x}{
numeric: vector of nonnegative elements
}
  \item{a_y}{
numeric: vector of nonnegative elements, whose length (p) and the sum of elements (s) are equal to those of a_x
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
numeric: the mean of upper estimates of two asymmetric distances. The asymmetric distance between a_x and a_y is the volume of the intersection of a `p-1`-dimensional ball of radius `||a_x-a_y||` and the simplex `sum(x) = s` divided by the volume of that simplex.
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
symmetric_distance(a_x = c(1,2,3), a_y = c(3,2,1))
symmetric_distance(a_x = c(4,5,7,3,2), a_y = c(3,6,8,2,2))

## The function is currently defined as
function (a_x, a_y) 
{
    if (!(((is.numeric(a_x)) | (is.integer(a_x))) & ((is.numeric(a_y)) | 
        (is.integer(a_y))))) {
        stop("'a_x' and 'a_y' must be vectors of class 'numeric' or 'integer")
    }
    if (!(length(a_x) == length(a_y))) {
        stop("'a_x' and 'a_y' must have same length")
    }
    if (!(sum(a_x) == sum(a_y))) {
        stop("`sum(a_x)` must be equal to `sum(a_y)`")
    }
    if (!(all(a_x >= 0) & all(a_y >= 0))) {
        stop("all elements of both 'a_x' and 'a_y' bust be nonnegative")
    }
    Eucl_dist <- sqrt(sum((a_x - a_y)^2))
    p <- length(a_x)
    s <- sum(a_x)
    d <- pi^((p - 1)/2) * Eucl_dist^(p - 1) * factorial(p - 1)/gamma((p - 
        1)/2 + 1)/s^(p - 1)/sqrt(p)
    return(d)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
