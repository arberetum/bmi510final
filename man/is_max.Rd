% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmi510final.R
\name{is_max}
\alias{is_max}
\title{Return a logical vector that is TRUE at indices that have the maximum value
in the given vector.}
\usage{
is_max(x, na.rm = T)
}
\arguments{
\item{x}{An atomic vector.}

\item{na.rm}{A logical indicating whether to remove missing values}
}
\value{
A logical vector of the same length as x that is TRUE at indices that
have the maximum value in the given vector x
}
\description{
Return a logical vector that is TRUE at indices that have the maximum value
in the given vector.
}
\examples{
x <- c(8, 4, 8, 3, 5)
is_max(x)  # returns c(TRUE, FALSE, TRUE, FALSE, FALSE)
y <- c(3, NaN, 4, 6, 1)
is_max(y)  # returns c(FALSE, NA, FALSE, TRUE, FALSE)
is_max(y, FALSE)  # returns c(NA, NA, NA, NA, NA)
}
