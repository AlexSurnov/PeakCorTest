# A statistical test for the local maxima correlation (co-localization)

## Aim

Two mathematical functions are defined on the same segment, and one aims to
identify statistically, whether the locations of local maxima of the first function
are correlated to those of the second function.

## Idea

Given that both functions are defined on the same finite segment, the locations
of their local maxima are nothing else than the vectors that belong to a certain
compact domain in a higher-dimensional space. (Technically speaking, the
dimensionality of the second vector may differ from that of the first vector,
however, one can correct this by inserting zero-valued elements into the shorter
vector.) By adoptig the null-hypthesis that those vectors are distributed uniformily
over that compact domain, one can estimate the probability that the distance
between them is smaller than the measured value. This estimate represents the
desired p-value.

## Main Functions

The two wrapper functions are `peak_cor_test(x, y)` and
`peak_cor_test_with_smoothing(x, y)`.

## Author

The method is invented by A. T. Fomenko. The original publication is
A. T. Fomenko “Empirico-Statistical Analysis of
Narrative Material and its Applications to Historical Dating: Volume I: The
Development of the Statistical Tools.” 1994 edition. Publisher: Springer; 1993.
234 p., chapter 4.4.

## Creator

The package is made by Alex Surnov. See the short description in `man/Fomenko peak correlation test.pdf`.
