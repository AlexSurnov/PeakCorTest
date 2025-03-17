\name{add_one_zero}
\alias{add_one_zero}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Insert one zero into a sequence
}
\description{
This function generates all possible insertions of a single zero between nonzero elements of an input vector.
}
\usage{
add_one_zero(input_vector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_vector}{
numeric: vector of nonnegative elements
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of numeric vectors, where each element differs from the input by the insertion of one zero before a nonzero element of the input.
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
add_one_zero(c(1,2,3,4,5))
add_one_zero(c(1,0,0,2,3,0,4,5))

## The function is currently defined as
function (input_vector) 
{
    nonzero_elements_index <- which(!(input_vector == 0))
    y <- list(c(0, input_vector))
    k <- 1
    for (n in nonzero_elements_index[2:length(nonzero_elements_index)]) {
        k <- k + 1
        y[[k]] <- c(input_vector[1:nonzero_elements_index[k] - 
            1], 0, input_vector[nonzero_elements_index[k]:length(input_vector)])
    }
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
