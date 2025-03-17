\name{find_constant_segments}
\alias{find_constant_segments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Splitting a sequence into segments with constant value
}
\description{
This function identifies the substrings of indices with constant values of the input vector.
}
\usage{
find_constant_segments(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a vector (numeric, logical, or character)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list labeled with unique values of `x`. Each element is in turn a list labeled with the index of the substring (i.e. segment), and its elements are the indices of the input vector's elements that form the corresponding segment.
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
find_constant_segments(c("b","b","b","a","a","a","c","c","a","a","c","b","b","b","b"))

## The function is currently defined as
function (x) 
{
    y <- which(!(x[1:(length(x) - 1)] == x[2:length(x)]))
    y <- as.matrix(data.frame(c(1, (y + 1)), c((y), length(x))))
    x <- data.frame(x)
    x$index <- 1:nrow(x)
    x$segment <- NA
    for (n in 1:nrow(y)) {
        x$segment[y[n, 1]:y[n, 2]] <- n
    }
    x <- split(x, x$x)
    x <- lapply(X = x, FUN = function(y) {
        lapply(X = split(x = y, f = y$segment), FUN = function(z) {
            z$index
        })
    })
    return(x)
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
