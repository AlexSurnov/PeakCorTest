\name{peak_cor_test}
\alias{peak_cor_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Running the Peak Correlation Test
}
\description{
This function runs `loc_max_vector()` (and `add_all_possible_zeros()` if necessary) on `x` and `y` and computes `symmetric_distance()`. If `add_all_possible_zeros()` was run, the function returns the minimal value across all variants.
}
\usage{
peak_cor_test(x, y, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric: vector
}
  \item{y}{
numeric: vector
}
  \item{verbose}{
logical: if `TRUE`, the function prints a message in case of failure
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
numeric: an upper estimate of the probability that the the local maxima of the two input vectors belong to the neighborhood of each other by chance
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
peak_cor_test(x = c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1),
y = c(1,2,3,2,1,1,0,1,2,1,0,-1,0,1,2,1,2,3,2,2,1))
peak_cor_test(x = c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1),
y = -c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1))
peak_cor_test(x = c(1,2,1,1,1), y = c(1,2,1,2,1), verbose = TRUE)
peak_cor_test(x = c(1:5), y = c(1,2,1,2,1), verbose = TRUE)
peak_cor_test(x = c(1:5), y = c(1,2,1,2,1), verbose = FALSE)

## The function is currently defined as
function (x, y, verbose = FALSE) 
{
    a <- list(x = x, y = y)
    a <- lapply(X = a, FUN = function(x) {
        try(loc_max_vector(x = x)[["a"]], silent = !verbose)
    })
    if (length(a$x) == length(a$y)) {
        return(symmetric_distance(a_x = a$x, a_y = a$y))
    }
    else if (length(a$x) < length(a$y)) {
        a <- a[2:1]
        names(a) <- c("x", "y")
    }
    y_variants <- try(add_all_possible_zeros(a_y = a$y, desired_length = length(a$x)), 
        silent = !verbose)
    distance_variants <- try(apply(X = y_variants, MARGIN = 2, 
        FUN = function(a_y) {
            symmetric_distance(a_x = a$x, a_y = a_y)
        }), silent = !verbose)
    if (paste(class(distance_variants), collapse = ", ") == "try-error") {
        return(NA)
    }
    return(min(c(1, distance_variants)))
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
