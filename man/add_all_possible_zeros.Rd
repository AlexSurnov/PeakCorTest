\name{add_all_possible_zeros}
\alias{add_all_possible_zeros}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Insert zeros into a sequence to reach the desired length
}
\description{
This function generates all possible insertions of zeros before the nonzero elements of an input vector to reach the specified length.
}
\usage{
add_all_possible_zeros(a_y, desired_length)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a_y}{
numeric: vector of nonnegative elements
}
  \item{desired_length}{
integer: the desired length of the output vectors
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a matrix whose columns differ from `a_y` by insertion of zeros before the nonzero elements of `a_y`.
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
add_all_possible_zeros(a_y = 1:5, desired_length = 8)

## The function is currently defined as
function (a_y, desired_length) 
{
    if (!((is.numeric(a_y)) | (is.integer(a_y)))) {
        stop("'a_y' must be a vector of class 'numeric' or 'integer'")
    }
    if (!all(a_y >= 0)) {
        stop("all elements 'a_y' bust be nonnegative")
    }
    if (!((is.numeric(desired_length)) | (is.integer(desired_length)))) {
        stop("'desired_length' must be a vector of class 'numeric' or 'integer'
        whose length is equal to 1")
    }
    if (!(length(desired_length) == 1)) {
        stop("'desired_length' must be a vector of class 'numeric' or 'integer'
        whose length is equal to 1")
    }
    if (desired_length <= length(a_y)) {
        stop("'desired_length' cannot be smaller than `length(a_y)`")
    }
    y <- add_one_zero(input_vector = a_y)
    if (desired_length == (length(a_y) + 1)) {
        y <- as.matrix(as.data.frame(y))
        rownames(y) <- NULL
        colnames(y) <- NULL
        return(y)
    }
    for (n in 2:(desired_length - length(a_y))) {
        y <- lapply(X = y, FUN = add_one_zero)
        y <- unlist(y, recursive = FALSE)
        y <- as.data.frame(t(as.data.frame(y)))
        rownames(y) <- NULL
        y <- dplyr::distinct(y)
        y <- as.list(as.data.frame(t(y)))
    }
    y <- as.matrix(as.data.frame(y))
    rownames(y) <- NULL
    colnames(y) <- NULL
    return(y)
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
