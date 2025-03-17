#' The function identifies the substrings of indices with constant values of the input vector.
#' @name find_constant_segments
#'
#' @param x a vector (numeric, logical or character)
#'
#' @returns A list tiled with unique values of `x`. Each element is in turn a list titled with the index of the substring (i.e. segment), and its elements are the indices of the input vector's elements that form the corresponding segment.
#' @examples
#' # example code
#'  find_constant_segments(c("b","b","b","a","a","a","c","c","a","a","c","b","b","b","b"))
#' @export
find_constant_segments<-function(x)
{
  ### y is a vector (numeric, logical or character)
  ### The function identifies the substrings of indices with constant values of y
  ### Example
  ###     y: c("b","b","b","a","a","a","c","c","a","a","c","b","b","b","b")
  ###     find_constant_segments(y):
  ###
  ###       $a
  ###       $a$`2`
  ###       [1] 4 5 6
  ###
  ###       $a$`4`
  ###       [1]  9 10
  ###
  ###
  ###       $b
  ###       $b$`1`
  ###       [1] 1 2 3
  ###
  ###       $b$`6`
  ###       [1] 12 13 14 15
  ###
  ###
  ###       $c
  ###       $c$`3`
  ###       [1] 7 8
  ###
  ###       $c$`5`
  ###       [1] 11
  
  y<-which(!(x[1:(length(x)-1)]==x[2:length(x)]))
  y<-as.matrix(data.frame(c(1,(y+1)),c((y),length(x))))
  x<-data.frame(x)
  x$index<-1:nrow(x)
  x$segment<-NA
  for (n in 1:nrow(y))
  {
    x$segment[y[n,1]:y[n,2]]<-n
  }
  x<-split(x, x$x)
  x<-lapply(X = x, FUN = function(y) {lapply(X = split(x = y, f = y$segment), FUN = function(z) {z$index})})
  return(x)
}

#' This function generates the sequence of local maxima of the input vector.
#' @name loc_max_vector
#'
#' @param x numeric: vector
#'
#' @returns list: `[[1]]` a vector whose first element is the index of the first local maximum, followed by the distances between adjacent local maxima, and concluded by the distance between the last local maximum and the last element of `x`. If a maximum is flat, the function consider the mean of its indices. `[[2]]` a data frame with the position of local maxima and the peak height.
#' @export
#' @examples
#' # example code
#'  loc_max_vector(x = c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1))
#'  loc_max_vector(x = c(1,0,1,2,1,0,1))
#'  loc_max_vector(x = c(-1,0,1,2,1,0,-1))
#'  loc_max_vector(x = c(1:5))
#' @importFrom data.table rbindlist
loc_max_vector<-function(x)
{
  ### identifying peaks and wells
  loc_extr<-x[2:length(x)]-x[1:(length(x)-1)]
  loc_extr<-which((loc_extr[2:length(loc_extr)]*loc_extr[1:(length(loc_extr)-1)])<=0)+1
  if (length(loc_extr)==0)
  {
    return(list(a = loc_extr,
                peak_height = data.frame(loc_max_coord = vector(length = 0),
                                         peak_height = vector(length = 0))))
  }
  if (length(loc_extr)==1)
  {
    peak_height<-x[loc_extr]-x[1]
    if (peak_height>0)
    {
      return(list(a = c(loc_extr,length(x)-loc_extr),
                  peak_height = data.frame(loc_max_coord = loc_extr,
                                           peak_height = peak_height)))
    } else {
      return(list(a = c(1,length(x)-1),
                  peak_height = data.frame(loc_max_coord = c(1,length(x)),
                                           peak_height = c(-peak_height,x[length(x)]-x[loc_extr]))))
    }
  }
  if (any(loc_extr[2:length(loc_extr)]==(loc_extr[1:(length(loc_extr)-1)]+1)))
  {
    x_loc_extr<-x[loc_extr]
    x_loc_extr_const_seg<-find_constant_segments(x_loc_extr)
    x_loc_extr_const_seg_new<-x_loc_extr_const_seg
    for (n in 1:length(x_loc_extr_const_seg))
    {
      for (k in 1:length(x_loc_extr_const_seg[[n]]))
      {
        xx<-loc_extr[x_loc_extr_const_seg[[n]][[k]]]
        xx<-x[c(xx[1]-1,xx[1],xx[length(xx)],xx[length(xx)]+1)]
        xx<-c(xx[2]-xx[1],xx[4]-xx[3])
        if (sign(xx[1])==sign(xx[2]))
        {
          x_loc_extr_const_seg_new[[n]][[k]]<-NULL
        }
      }
    }
    y<-lapply(X = x_loc_extr_const_seg_new,
              FUN = function(z) {
                lapply(X = z,
                       FUN = function(r) {
                         round(x = mean(r), digits = 0)
                       })
              })
    for (n in names(y))
    {
      y[[n]]<-data.frame(extremum_index = unlist(y[[n]]))
      rownames(y[[n]])<-NULL
    }
    y<-as.data.frame(data.table::rbindlist(y))
    loc_extr<-loc_extr[y[order(y$extremum_index),]]
  }
  if (length(loc_extr)==1)
  {
    peak_height<-x[loc_extr]>x[1]
    if (peak_height>0)
    {
      return(list(a = c(loc_extr,length(x)-loc_extr),
                  peak_height = data.frame(loc_max_coord = loc_extr,
                                           peak_height = peak_height)))
    } else {
      return(list(a = c(1,length(x)-1),
                  peak_height = data.frame(loc_max_coord = c(1,length(x)),
                                           peak_height = c(-peak_height,x[length(x)]-x[loc_extr]))))
    }
  }
  if (x[loc_extr[1]]<x[loc_extr[1]+1])
  {
    loc_extr<-list(MIN = loc_extr[seq(from = 1, to = length(loc_extr), by = 2)],
                   MAX = loc_extr[seq(from = 2, to = length(loc_extr), by = 2)])
  } else {
    loc_extr<-list(MIN = loc_extr[seq(from = 2, to = length(loc_extr), by = 2)],
                   MAX = loc_extr[seq(from = 1, to = length(loc_extr), by = 2)])
  }
  if (length(loc_extr$MAX)==0)
  {
    return(loc_extr$MAX)
  }
  ### computing peak height
  if (loc_extr$MIN[1]<loc_extr$MAX[1])
  {
    if (loc_extr$MIN[length(loc_extr$MIN)]<loc_extr$MAX[length(loc_extr$MAX)])
    {
      peak_height<-data.frame(x[loc_extr$MAX]-x[loc_extr$MIN],
                              c(x[loc_extr$MAX[1:(length(loc_extr$MAX)-1)]]-x[loc_extr$MIN[2:length(loc_extr$MIN)]],
                                x[loc_extr$MAX[length(loc_extr$MAX)]]-x[length(x)]))
    } else {
      peak_height<-data.frame(x[loc_extr$MAX]-x[loc_extr$MIN[1:(length(loc_extr$MIN)-1)]],
                              x[loc_extr$MAX]-x[loc_extr$MIN[2:length(loc_extr$MIN)]])
    }
  } else {
    if (loc_extr$MIN[length(loc_extr$MIN)]<loc_extr$MAX[length(loc_extr$MAX)])
    {
      peak_height<-data.frame(c(x[loc_extr$MAX[1]]-x[1],
                                x[loc_extr$MAX[2:length(loc_extr$MAX)]]-x[loc_extr$MIN]),
                              c(x[loc_extr$MAX[1:(length(loc_extr$MAX)-1)]]-x[loc_extr$MIN],
                                x[loc_extr$MAX[length(loc_extr$MAX)]]-x[length(x)]))
    } else {
      peak_height<-data.frame(c(x[loc_extr$MAX[1]]-x[1],
                                x[loc_extr$MAX[2:length(loc_extr$MAX)]]-x[loc_extr$MIN[1:(length(loc_extr$MIN)-1)]]),
                              x[loc_extr$MAX]-x[loc_extr$MIN])
    }
  }
  peak_height<-apply(X = peak_height, MARGIN = 1, FUN = max)
  peak_height<-data.frame(loc_max_coord = loc_extr$MAX, peak_height = peak_height)
  loc_extr<-c(loc_extr$MAX,length(x))
  return(list(a = c(loc_extr[1],loc_extr[2:length(loc_extr)]-loc_extr[1:(length(loc_extr)-1)]),
              peak_height = peak_height))
}

#' This function generates an upper estimate of the probability that the the two input vectors belong to the neighborhood of each other by chance.
#' @name symmetric_distance
#'
#' @param a_x numeric: vector of nonnegative elements
#' @param a_y numeric: vector of nonnegative elements, whose length (p) and the sum of elements (s) are equal to those of a_x
#'
#' @returns numeric: the mean of upper estimates of two asymmetric distances. The asymmetric distance between a_x and a_y is the volume of the intersection of a `p-1`-dimensional ball of radius `||a_x-a_y||` and the simplex `sum(x) = s` divided by the volume of that simplex.
#' @export
#' @examples
#' # example code
#'  symmetric_distance(a_x = c(1,2,3), a_y = c(3,2,1))
#'  symmetric_distance(a_x = c(4,5,7,3,2), a_y = c(3,6,8,2,2))
symmetric_distance<-function(a_x, a_y)
{
  ### 'a_x' and 'a_y' are the numeric vectors of the same length and of the same sum.
  if (!(
    (
      (is.numeric(a_x)) | (is.integer(a_x))
    ) & (
      (is.numeric(a_y)) | (is.integer(a_y))
    )
  ))
  {
    stop("'a_x' and 'a_y' must be vectors of class 'numeric' or 'integer")
  }
  if (!(length(a_x)==length(a_y)))
  {
    stop("'a_x' and 'a_y' must have same length")
  }
  if (!(sum(a_x)==sum(a_y)))
  {
    stop("`sum(a_x)` must be equal to `sum(a_y)`")
  }
  if (!(all(a_x>=0) & all(a_y>=0)))
  {
    stop("all elements of both 'a_x' and 'a_y' bust be nonnegative")
  }
  
  Eucl_dist<-sqrt(sum((a_x-a_y)^2))
  p<-length(a_x)
  s<-sum(a_x)
  d<-pi^((p-1)/2) * Eucl_dist^(p-1) * factorial(p-1) / gamma((p-1)/2+1) / s^(p-1) / sqrt(p)
  return(d)
}

#' This function generates all possible insertions of a single zero between nonzero elements of an input vector.
#' @name add_one_zero
#'
#' @param input_vector numeric: vector of nonnegative elements
#'
#' @returns a list of numeric vectors. Each its element differs from the input by the insertion of one zero before a nonzero element of the input.
#' @export
#' @examples
#' # example code
#'  add_one_zero(c(1,2,3,4,5))
#'  add_one_zero(c(1,0,0,2,3,0,4,5))
add_one_zero<-function(input_vector)
{
  nonzero_elements_index<-which(!(input_vector==0))
  y<-list(c(0,input_vector))
  k<-1
  for (n in nonzero_elements_index[2:length(nonzero_elements_index)])
  {
    k<-k+1
    y[[k]]<-c(input_vector[1:nonzero_elements_index[k]-1],0,input_vector[nonzero_elements_index[k]:length(input_vector)])
  }
  return(y)
}

#' This function generates all possible insertions of zeros between nonzero elements of an input vector to reach the specified length.
#' @name add_all_possible_zeros
#'
#' @param a_y numeric: vector of nonnegative elements
#' @param desired_length integer: the desired length of the output vectors
#'
#' @returns a matrix whose columns differ from `a_y` by insertion of zeros before the nonzero elements of `a_y`.
#' @export
#' @examples
#' # example code
#'  add_all_possible_zeros(a_y = 1:5, desired_length = 8)
#' @importFrom dplyr distinct
add_all_possible_zeros<-function(a_y, desired_length)
{
  if (!(
    (is.numeric(a_y)) | (is.integer(a_y))
  ))
  {
    stop("'a_y' must be a vector of class 'numeric' or 'integer'")
  }
  if (!all(a_y>=0))
  {
    stop("all elements 'a_y' bust be nonnegative")
  }
  
  if (!(
    (is.numeric(desired_length)) | (is.integer(desired_length))
  ))
  {
    stop("'desired_length' must be a vector of class 'numeric' or 'integer' whose length is equal to 1")
  }
  if (!(length(desired_length)==1))
  {
    stop("'desired_length' must be a vector of class 'numeric' or 'integer' whose length is equal to 1")
  }
  if (desired_length<=length(a_y))
  {
    stop("'desired_length' cannot be smaller than `length(a_y)`")
  }
  
  y<-add_one_zero(input_vector = a_y)
  if (desired_length==(length(a_y)+1))
  {
    y<-as.matrix(as.data.frame(y))
    rownames(y)<-NULL
    colnames(y)<-NULL
    return(y)
  }
  
  for (n in 2:(desired_length-length(a_y)))
  {
    y<-lapply(X = y,
              FUN = add_one_zero)
    y<-unlist(y, recursive = FALSE)
    ### removing duplicates
    y<-as.data.frame(t(as.data.frame(y)))
    rownames(y)<-NULL
    y<-dplyr::distinct(y)
    y<-as.list(as.data.frame(t(y)))
  }
  y<-as.matrix(as.data.frame(y))
  rownames(y)<-NULL
  colnames(y)<-NULL
  return(y)
}

#' This function runs `loc_max_vector` (and `add_all_possible_zeros` if necessary) on `x` and `y` and computes `symmetric_distance`. If `add_all_possible_zeros` was run, the function returns the minimal value across all variants.
#' @name peak_cor_test
#'
#' @param x numeric: vector
#' @param y numeric: vector
#' @param verbose logical: if `TRUE`, the function prints a message in case of failure
#'
#' @returns numeric: an upper estimate of the probability that the the local maxima of the two input vectors belong to the neighborhood of each other by chance
#' @export
#' @examples
#' # example code
#'  peak_cor_test(x = c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1), y = c(1,2,3,2,1,1,0,1,2,1,0,-1,0,1,2,1,2,3,2,2,1))
#'  peak_cor_test(x = c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1), y = -c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1))
#'  peak_cor_test(x = c(1,2,1,1,1), y = c(1,2,1,2,1), verbose = TRUE)
#'  peak_cor_test(x = c(1:5), y = c(1,2,1,2,1), verbose = TRUE)
#'  peak_cor_test(x = c(1:5), y = c(1,2,1,2,1), verbose = FALSE)
peak_cor_test<-function(x, y, verbose = FALSE)
{
  a<-list(x = x, y = y)
  a<-lapply(X = a, FUN = function(x) {try(loc_max_vector(x = x)[["a"]],
                                          silent = !verbose)})
  if (length(a$x)==length(a$y))
  {
    return(symmetric_distance(a_x = a$x, a_y = a$y))
  } else if (length(a$x)<length(a$y)) {
    a<-a[2:1]
    names(a)<-c("x","y")
  }
  y_variants<-try(add_all_possible_zeros(a_y = a$y, desired_length = length(a$x)),
                  silent = !verbose)
  distance_variants<-try(apply(X = y_variants,
                               MARGIN = 2,
                               FUN = function(a_y) {
                                 symmetric_distance(a_x = a$x, a_y = a_y)
                               }),
                         silent = !verbose)
  if (paste(class(distance_variants), collapse = ", ")=="try-error")
  {
    return(NA)
  }
  return(min(c(1,distance_variants)))
}

#' This function runs `loc_max_vector`, `smooth_n_cpp` with `n in 1:length(input)` and `add_all_possible_zeros` (if necessary) on `x` and `y` and computes `symmetric_distance`. The function returns the minimal value across all variants.
#' @name peak_cor_test_with_smoothing
#'
#' @param x numeric: vector
#' @param y numeric: vector
#' @param verbose logical: if `TRUE`, the function prints a message in case of failure
#' @param max_peak_count_difference numeric: the maximum difference in peak count, above which the test will be omitted. This parameter is included in order to avoid wasting efforts on computing the distance between the curves with obviously different number of peaks.
#' @param max_smoothing_degree positive integer: maximal smoothing degree
#'
#' @returns numeric: an upper estimate of the probability that the the local maxima of the two input vectors belong to the neighborhood of each other by chance
#' @export
#' @examples
#' # example code
#'  peak_cor_test_with_smoothing(x = c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1), y = c(1,2,3,2,1,1,0,1,2,1,0,-1,0,1,2,1,2,3,2,2,1))
#'  peak_cor_test_with_smoothing(x = c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1), y = -c(1,2,3,3,2,1,0,1,2,1,0,-1,0,1,1,2,1,2,3,2,1))
#'  peak_cor_test_with_smoothing(x = c(1,2,1,1,1), y = c(1,2,1,2,1), verbose = TRUE)
#'  peak_cor_test_with_smoothing(x = c(1,2,1,1,1), y = c(1,2,1,2,1), verbose = FALSE)
#'  peak_cor_test_with_smoothing(x = c(1:5), y = c(1,2,1,2,1), verbose = TRUE)
#'  peak_cor_test_with_smoothing(x = c(1:5), y = c(1,2,1,2,1), verbose = FALSE)
#' @import Matrix
#' @importFrom data.table rbindlist
#' @importClassesFrom Matrix `sparseMatrix-class`
peak_cor_test_with_smoothing<-function(x, y, max_peak_count_difference = 5, max_smoothing_degree = length(x), verbose = FALSE)
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
    colnames(peak_count_difference)<-paste(c("x","y"),"smooting_degree", sep = "_")
  }
  for (n in c("x","y"))
  {
    nn<-paste(n,"smooting_degree", sep = "_")
    peak_count_difference[,nn]<-peak_count[[n]]$smoothing_degree[peak_count_difference[,nn]]
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
  n<-which(peak_count_difference$symmetric_distance==min(peak_count_difference$symmetric_distance, na.rm = TRUE))
  return(peak_count_difference[n,])
}






