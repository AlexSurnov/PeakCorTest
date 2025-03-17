\name{PeakCorTest-package}
\alias{PeakCorTest-package}
\alias{PeakCorTest}
\docType{package}
\title{
  A statistical test for the local maxima correlation
}
\description{
  The method is invented by A. T. Fomenko.
  Two mathematical functions are defined on the same segment, and one aims to
  identify statistically, whether the locations of local maxima of the first function
  are correlated to those of the second function.
  
  Given that both functions are defined on the same finite segment, the locations
  of their local maxima are nothing else than the vectors that belong to a certain
  compact domain in a higher-dimensional space. (Technically speaking, the
  dimensionality of the second vector may differ from that of the first vector,
  however, one can correct this by inserting zero-valued elements into the shorter
  vector.) By adoptig the null-hypthesis that those vectors are distributed uniformily
  over that compact domain, one can estimate the probability that the distance
  between them is smaller than the measured value. This estimate represents the
  desired p-value.
}
\details{
  The two wrapper functions are peak_cor_test(x, y, verbose = FALSE) and
  peak_cor_test_with_smoothing(x, y, max_peak_count_difference = 5, verbose = FALSE)
}
\author{
Alex Surnov, alexsurnov88@gmail.com

Maintainer: Alex Surnov <alexsurnov88@gmail.com>
}
\references{
  bibentry(
  bibtype  = "Book",
  title    = "Empirico-Statistical Analysis of Narrative Material and its Applications to Historical Dating: Volume I: The Development of the Statistical Tools",
  author   = "Fomenko, A. T.",
  year     = "1994",
  chapter   = "4.4",
  publisher = "Springer"
)

}
\keyword{ package }
\seealso{
  The Peak Correlation Test description.pdf
}
