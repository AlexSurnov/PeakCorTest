\name{loc_max_vector}
\alias{loc_max_vector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Getting local maxima of a sequence
}
\description{
This function generates the sequence of local maxima of the input vector.
}
\usage{
loc_max_vector(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric: vector
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
list: `[[1]]` a vector whose first element is the index of the first local maximum, followed by the distances between adjacent local maxima, and concluded by the distance between the last local maximum and the last element of `x`. If a maximum is flat, the function considers the mean of its indices. `[[2]]` a data frame with the position of local maxima and the peak height.
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
loc_max_vector(x = c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1))
loc_max_vector(x = c(1,0,1,2,1,0,1))
loc_max_vector(x = c(-1,0,1,2,1,0,-1))
loc_max_vector(x = c(1:5))

## The function is currently defined as
function (x) 
{
    loc_extr <- x[2:length(x)] - x[1:(length(x) - 1)]
    loc_extr <- which((loc_extr[2:length(loc_extr)] * loc_extr[1:(length(loc_extr) - 
        1)]) <= 0) + 1
    if (length(loc_extr) == 0) {
        return(list(a = loc_extr, peak_height = data.frame(loc_max_coord = vector(length = 0), 
            peak_height = vector(length = 0))))
    }
    if (length(loc_extr) == 1) {
        peak_height <- x[loc_extr] - x[1]
        if (peak_height > 0) {
            return(list(a = c(loc_extr, length(x) - loc_extr), 
                peak_height = data.frame(loc_max_coord = loc_extr, 
                  peak_height = peak_height)))
        }
        else {
            return(list(a = c(1, length(x) - 1), peak_height = data.frame(loc_max_coord = c(1, 
                length(x)), peak_height = c(-peak_height, x[length(x)] - 
                x[loc_extr]))))
        }
    }
    if (any(loc_extr[2:length(loc_extr)] == (loc_extr[1:(length(loc_extr) - 
        1)] + 1))) {
        x_loc_extr <- x[loc_extr]
        x_loc_extr_const_seg <- find_constant_segments(x_loc_extr)
        x_loc_extr_const_seg_new <- x_loc_extr_const_seg
        for (n in 1:length(x_loc_extr_const_seg)) {
            for (k in 1:length(x_loc_extr_const_seg[[n]])) {
                xx <- loc_extr[x_loc_extr_const_seg[[n]][[k]]]
                xx <- x[c(xx[1] - 1, xx[1], xx[length(xx)], xx[length(xx)] + 
                  1)]
                xx <- c(xx[2] - xx[1], xx[4] - xx[3])
                if (sign(xx[1]) == sign(xx[2])) {
                  x_loc_extr_const_seg_new[[n]][[k]] <- NULL
                }
            }
        }
        y <- lapply(X = x_loc_extr_const_seg_new, FUN = function(z) {
            lapply(X = z, FUN = function(r) {
                round(x = mean(r), digits = 0)
            })
        })
        for (n in names(y)) {
            y[[n]] <- data.frame(extremum_index = unlist(y[[n]]))
            rownames(y[[n]]) <- NULL
        }
        y <- as.data.frame(data.table::rbindlist(y))
        loc_extr <- loc_extr[y[order(y$extremum_index), ]]
    }
    if (length(loc_extr) == 1) {
        peak_height <- x[loc_extr] > x[1]
        if (peak_height > 0) {
            return(list(a = c(loc_extr, length(x) - loc_extr), 
                peak_height = data.frame(loc_max_coord = loc_extr, 
                  peak_height = peak_height)))
        }
        else {
            return(list(a = c(1, length(x) - 1), peak_height = data.frame(loc_max_coord = c(1, 
                length(x)), peak_height = c(-peak_height, x[length(x)] - 
                x[loc_extr]))))
        }
    }
    if (x[loc_extr[1]] < x[loc_extr[1] + 1]) {
        loc_extr <- list(MIN = loc_extr[seq(from = 1, to = length(loc_extr), 
            by = 2)], MAX = loc_extr[seq(from = 2, to = length(loc_extr), 
            by = 2)])
    }
    else {
        loc_extr <- list(MIN = loc_extr[seq(from = 2, to = length(loc_extr), 
            by = 2)], MAX = loc_extr[seq(from = 1, to = length(loc_extr), 
            by = 2)])
    }
    if (length(loc_extr$MAX) == 0) {
        return(loc_extr$MAX)
    }
    if (loc_extr$MIN[1] < loc_extr$MAX[1]) {
        if (loc_extr$MIN[length(loc_extr$MIN)] < loc_extr$MAX[length(loc_extr$MAX)]) {
            peak_height <- data.frame(x[loc_extr$MAX] - x[loc_extr$MIN], 
                c(x[loc_extr$MAX[1:(length(loc_extr$MAX) - 1)]] - 
                  x[loc_extr$MIN[2:length(loc_extr$MIN)]], x[loc_extr$MAX[length(loc_extr$MAX)]] - 
                  x[length(x)]))
        }
        else {
            peak_height <- data.frame(x[loc_extr$MAX] - x[loc_extr$MIN[1:(length(loc_extr$MIN) - 
                1)]], x[loc_extr$MAX] - x[loc_extr$MIN[2:length(loc_extr$MIN)]])
        }
    }
    else {
        if (loc_extr$MIN[length(loc_extr$MIN)] < loc_extr$MAX[length(loc_extr$MAX)]) {
            peak_height <- data.frame(c(x[loc_extr$MAX[1]] - 
                x[1], x[loc_extr$MAX[2:length(loc_extr$MAX)]] - 
                x[loc_extr$MIN]), c(x[loc_extr$MAX[1:(length(loc_extr$MAX) - 
                1)]] - x[loc_extr$MIN], x[loc_extr$MAX[length(loc_extr$MAX)]] - 
                x[length(x)]))
        }
        else {
            peak_height <- data.frame(c(x[loc_extr$MAX[1]] - 
                x[1], x[loc_extr$MAX[2:length(loc_extr$MAX)]] - 
                x[loc_extr$MIN[1:(length(loc_extr$MIN) - 1)]]), 
                x[loc_extr$MAX] - x[loc_extr$MIN])
        }
    }
    peak_height <- apply(X = peak_height, MARGIN = 1, FUN = max)
    peak_height <- data.frame(loc_max_coord = loc_extr$MAX, peak_height = peak_height)
    loc_extr <- c(loc_extr$MAX, length(x))
    return(list(a = c(loc_extr[1], loc_extr[2:length(loc_extr)] - 
        loc_extr[1:(length(loc_extr) - 1)]), peak_height = peak_height))
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
