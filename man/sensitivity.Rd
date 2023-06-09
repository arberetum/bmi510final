% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmi510final.R
\name{sensitivity}
\alias{sensitivity}
\title{Calculate sensitivity as the ratio of true positives to total ground truth
positives (true positives + false negatives).}
\usage{
sensitivity(pred, truth)
}
\arguments{
\item{pred}{A vector of predictions that are either logical or numeric, with
0 being false and 1 being true.}

\item{truth}{A vector of ground truths that are either logical or numeric, with
0 being false and 1 being true.}
}
\value{
The sensitivity as a single number.
}
\description{
Calculate sensitivity as the ratio of true positives to total ground truth
positives (true positives + false negatives).
}
\examples{
pred <- c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)
truth <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
sensitivity(pred, truth)  # returns 0.6
}
