\name{peak_cor_test_with_smoothing}
\alias{peak_cor_test_with_smoothing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Running the Peak Correlation Test preceded by smoothing
}
\description{
#' This function runs `loc_max_vector()`, `smooth_n_cpp()` with `n in 1:length(input)` and `add_all_possible_zeros()` (if necessary) on `x` and `y` and computes `symmetric_distance()`. The function returns the minimal value across all variants.
}
\usage{
peak_cor_test_with_smoothing(x, y,
  max_peak_count_difference = 5,
  max_smoothing_degree = length(x),
  verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric: vector
}
  \item{y}{
numeric: vector
}
  \item{max_smoothing_degree}{
positive integer: maximal smoothing degree
}
  \item{max_peak_count_difference}{
numeric: the maximum difference in peak count, above which the test will be omitted. This parameter is included in order to avoid wasting efforts on computing the distance between the curves with obviously different number of peaks.
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
peak_cor_test_with_smoothing(x = c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1),
y = c(1,2,3,2,1,1,0,1,2,1,0,-1,0,1,2,1,2,3,2,2,1))
peak_cor_test_with_smoothing(x = c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1),
y = -c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1))
peak_cor_test_with_smoothing(x = c(1,2,1,1,1), y = c(1,2,1,2,1),
verbose = TRUE)
peak_cor_test_with_smoothing(x = c(1,2,1,1,1), y = c(1,2,1,2,1),
verbose = FALSE)
peak_cor_test_with_smoothing(x = c(1:5), y = c(1,2,1,2,1), verbose = TRUE)
peak_cor_test_with_smoothing(x = c(1:5), y = c(1,2,1,2,1), verbose = FALSE)

## The function is currently defined as
function (x, y, max_peak_count_difference = 5, verbose = FALSE) 
{
  L<-list(x = x, y = y)
  L<-lapply(X = L,
            FUN = function(x) {
              y<-1:max_smoothing_degree
              apply(X = as.matrix(y),
                    MARGIN = 1,
                    FUN = function(n) {
                      smooth_n_cpp(x = x, n = n)
                    })
            })
  peak_count<-lapply(X = L,
                     FUN = function(x) {
                       x<-apply(X = x,
                                MARGIN = 2,
                                FUN = function(y) {
                                  length(loc_max_vector(y)[["a"]])
                                })
                       x<-data.frame(smoothing_degree = 1:max_smoothing_degree,
                                     peak_count = x)
                       x<-split(x = x, f = x$peak_count)
                       x<-lapply(X = x,
                                 FUN = function(y) {
                                   y[1,]
                                 })
                       x<-as.data.frame(data.table::rbindlist(x))
                       x<-x[order(x$smoothing_degree),]
                     })
  peak_count_difference<-abs(
    matrix(
      data = rep(x = peak_count$x$peak_count,
                 times = nrow(peak_count$y)),
      nrow = nrow(peak_count$x)
    ) - t(matrix(
      data = rep(x = peak_count$y$peak_count,
                 times = nrow(peak_count$x)),
      nrow = nrow(peak_count$y)
    ))
  )
  peak_count_difference<-peak_count_difference<=max_peak_count_difference
  peak_count_difference<-Matrix::Matrix(peak_count_difference, sparse = TRUE)
  if (all(dim(peak_count_difference)==c(1,1)))
  {
    if (peak_count_difference[1,1])
    {
      peak_count_difference<-data.frame(x_smooting_degree = 1,
                                        y_smooting_degree = 1,
                                        symmetric_distance = NA)
    } else {
      return(data.frame(x_smooting_degree = NA,
                        y_smooting_degree = NA,
                        symmetric_distance = NA))
    }
  } else {
    peak_count_difference<-Matrix::summary(peak_count_difference)
    peak_count_difference<-peak_count_difference[,1:2]
    colnames(peak_count_difference)<-paste(c("x","y"),"smooting_degree",
    sep = "_")
  }
  for (n in c("x","y"))
  {
    nn<-paste(n,"smooting_degree", sep = "_")
    peak_count_difference[,nn]<-peak_count[[n]]$smoothing_degree[
      peak_count_difference[,nn]
    ]
  }
  peak_count_difference$symmetric_distance<-NA
  for (n in 1:nrow(peak_count_difference))
  {
    y<-peak_cor_test(x = L$x[,peak_count_difference$x_smooting_degree[n]],
                     y = L$y[,peak_count_difference$y_smooting_degree[n]],
                     verbose = verbose)
    if (paste(class(y), collapse = ", ")=="numeric")
    {
      peak_count_difference$symmetric_distance[n]<-y
    }
  }
  if ((nrow(peak_count_difference)==0) |
      (all(is.na(peak_count_difference$symmetric_distance))))
  {
    return(data.frame(x_smooting_degree = NA,
                      y_smooting_degree = NA,
                      symmetric_distance = NA))
  }
  n<-which(peak_count_difference$symmetric_distance==min(
    peak_count_difference$symmetric_distance, na.rm = TRUE
  ))
  return(peak_count_difference[n,])
}
}
