% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmi510final.R
\name{specificity}
\alias{specificity}
\title{Calculate specificity as the ratio of true negatives to total ground truth
negatives (true positives + false negatives).}
\usage{
specificity(pred, truth)
}
\arguments{
\item{pred}{A vector of predictions that are either logical or numeric, with
0 being false and 1 being true.}

\item{truth}{A vector of ground truths that are either logical or numeric, with
0 being false and 1 being true.}
}
\value{
The specificity as a single number.
}
\description{
Calculate specificity as the ratio of true negatives to total ground truth
negatives (true positives + false negatives).
}
\examples{
pred <- c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)
truth <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
specificity(pred, truth)  # returns 0.5
}
