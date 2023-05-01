#' Randomly sample items from a vector or rows from a dataframe-like object.
#'
#' @param x A vector or dataframe-like object.
#' @param n An integer number of samples.
#' @param replace Whether to sample with replacement.
#' @return Items or rows sampled from the input object.
#' @export
#' @examples
#' x <- 1:12
#' # get a random number from [1,12]
#' rando(x)
#' # get 5 random numbers from [1,12] without replacement
#' rando(x, 12, FALSE)
#' a <- data.frame("a" = c(1, 2, 3), "b" = c(4, 6, 8), "c" = c(10, 44, 32))
#' # sample 2 rows with replacement
#' rando(a, 2)
rando = function(x,n=1,replace=T) {
  if (is.atomic(x)) {
    return(sample(x, n, replace))
  } else if (is.data.frame(x)) {
    rows = seq(1, nrow(x))
    rows_rand = sample(rows, n, replace)
    return(x[rows_rand, ])
  }
}


#' Return a logical vector that is TRUE at indices that have the minimum value
#'   in the given vector.
#'
#' @param x An atomic vector.
#' @param na.rm A logical indicating whether to remove missing values.
#' @return A logical vector of the same length as x that is TRUE at indices that
#' have the minimum value in the given vector x.
#' @export
#' @examples
#' x <- c(3, 4, 8, 3, 5)
#' is_min(x)  # returns c(TRUE, FALSE, FALSE, TRUE, FALSE)
#' y <- c(3, NaN, 4, 6, 1)
#' is_min(y)  # returns c(FALSE, NA, FALSE, FALSE, TRUE)
#' is_min(y, FALSE)  # returns c(NA, NA, NA, NA, NA)
is_min = function(x,na.rm=T) {
  out = x == min(x, na.rm=na.rm)
}


#' Return a logical vector that is TRUE at indices that have the maximum value
#' in the given vector.
#'
#' @param x An atomic vector.
#' @param na.rm A logical indicating whether to remove missing values
#' @return A logical vector of the same length as x that is TRUE at indices that
#' have the maximum value in the given vector x
#' @export
#' @examples
#' x <- c(8, 4, 8, 3, 5)
#' is_max(x)  # returns c(TRUE, FALSE, TRUE, FALSE, FALSE)
#' y <- c(3, NaN, 4, 6, 1)
#' is_max(y)  # returns c(FALSE, NA, FALSE, TRUE, FALSE)
#' is_max(y, FALSE)  # returns c(NA, NA, NA, NA, NA)
is_max = function(x,na.rm=T) {
  out = x == max(x, na.rm=na.rm)
}


#' Replicate the given matrix or dataframe a given number of times along rows
#'   and/or columns.
#'
#' @param x A matrix or dataframe to replicate.
#' @param M The integer number of times to replicate the input along rows.
#' @param N The integer number of times to replicate the input along columns
#' @return A matrix with the replicated input matrix/dataframe
#' @export
#' @examples
#' x <- rbind(c(1,2), c(3,4))
#' rep_mat(x, 3, 1)  # replicates the matrix x 3 times row-wise
rep_mat = function(x, M=1, N=1) {
  nrow_x = nrow(x)
  ncol_x = ncol(x)
  out = matrix(NA, nrow=M*nrow_x, ncol=N*ncol_x)
  for (i in seq(1,M)) {
    for (j in seq(1,N)) {
      out[( (i-1)*nrow_x + 1 ) : ( (i-1)*nrow_x + nrow_x),
          ( (j-1)*ncol_x + 1 ) : ( (j-1)*ncol_x + ncol_x) ] = x
    }
  }
  return(out)
}


#' Return a character vector containing the classes of each variable in the
#'   given tibble.
#'
#' @param x A tibble.
#' @return A vector containing the classes of each variable in the tibble.
#' @export
#' @examples
#' a <- tibble::tibble("a" = c(TRUE, FALSE), "b" = c(3, 4), "c" = c("Hi", "Hello"))
#' classes(a)  # returns   "logical"   "numeric" "character"
classes = function(x) {
  sapply(x, class)
}


#' Center and/or scale numeric columns of the given tibble.
#'
#' @param x A tibble.
#' @param center A logical indicating whether to center the numeric columns by
#' subtracting the mean.
#' @param scale A logical indicating whether to scale the numeric columns by
#' dividing by their standard deviations.
#' @return A centered and/or scaled version of the input tibble.
#' @export
#' @examples
#' a <- tibble::tibble("a" = c(TRUE, FALSE, FALSE), "b" = c(1, 4, 5), "c" = c("Hi", "Hello", "What's up"))
#' df_scale(a)  # scales the "b" column so it becomes c(-1.12, 0.32, 0.801)
#' df_scale(a, scale = FALSE)  # centers the "b" column but doesn't center it
#' df_scale(a, center = FALSE)  # scales the "b" column but doesn't scale it
#' df_scale(a, scale = FALSE, center = FALSE)  # returns a unchanged
df_scale = function(x, center = T, scale = T) {
  if (center && scale) {
    f = function(col) (col - mean(col)) / stats::sd(col)
  } else if (center) {
    f = function(col) col - mean(col)
  } else if (scale) {
    f = function(col) col / stats::sd(col)
  } else {
    f = function(col) col
  }
  x[classes(x) == "numeric"] = sapply(x[classes(x) == "numeric"], f)
  return(x)
}


#' Return the log-likelihood of the given sample under a normal distribution
#'   with the given mean and standard deviation.
#'
#' @param x A numeric vector of samples.
#' @param mean The mean of the normal distribution.
#' @param sd The standard deviation of the normal distribution.
#' @return The log-likelihood of the given sample under a normal distribution
#' with the given mean and standard deviation.
#' @export
#' @examples
#' a <- stats::rnorm(100, 1, 3)  # generate 100 samples from a normal distribution with mean=1, sd=3
#' log_likelihood_norm(a, 1, 3)
log_likelihood_norm = function(x, mean, sd) {
  sum(log(stats::dnorm(x, mean=mean, sd=sd)))
}


#' Return the log-likelihood of the given sample under a uniform distribution
#'   with the given minimum and maximum.
#'
#' @param x A numeric vector of samples.
#' @param min The minimum of the uniform distribution.
#' @param max The maximum of the uniform distribution.
#' @return The log-likelihood of the given sample under a uniform distribution
#' with the given minimum and maximum.
#' @export
#' @examples
#' a <- stats::runif(100, 1, 1000)  # generate 100 samples from a uniform distribution on the interval [1,1000]
#' log_likelihood_unif(a, 1, 1000)
log_likelihood_unif = function(x, min, max) {
  sum(log(stats::dunif(x, min=min, max=max)))
}


#' Return the log-likelihood of the given sample under a chi-squared distribution
#' with the given degrees of freedom.
#'
#' @param x A numeric vector of samples.
#' @param df Degrees of freedom (non-negative, but can be non-integer).
#' @return The log-likelihood of the given sample under a chi-squared distribution
#' with the given degrees of freedom.
#' @export
#' @examples
#' a <- stats::rchisq(100, 33)  # generate 100 samples from a chi-squared distribution with df=33
#' log_likelihood_chisq(a, 33)
log_likelihood_chisq = function(x, df) {
  sum(log(stats::dchisq(x, df=df)))
}


#' Return the log-likelihood of the given sample under an F distribution
#'   with the given degrees of freedom df1 and df2.
#'
#' @param x A numeric vector of samples.
#' @param df1 Degrees of freedom. Inf is allowed.
#' @param df2 Degrees of freedom. Inf is allowed.
#' @return The log-likelihood of the given sample under an F distribution
#' with the given degrees of freedom df1 and df2.
#' @export
#' @examples
#' a <- stats::rf(100, 33, 40)  # generate 100 samples from an F distribution with df1=33 and df2=40
#' log_likelihood_f(a, 33, 40)
log_likelihood_f = function(x, df1, df2) {
  sum(log(stats::df(x, df1=df1, df2=df2)))
}


#' Return the log-likelihood of the given sample under a t distribution
#'   with the given degrees of freedom.
#'
#' @param x A numeric vector of samples.
#' @param df Degrees of freedom (non-negative, but can be non-integer, Inf allowed).
#' @return The log-likelihood of the given sample under a t distribution
#' with the given degrees of freedom.
#' @export
#' @examples
#' a <- stats::rt(100, 33)  # generate 100 samples from a t distribution with df=33
#' log_likelihood_t(a, 33)
log_likelihood_t = function(x, df) {
  sum(log(stats::dt(x, df=df)))
}


#' Calculate sensitivity as the ratio of true positives to total ground truth
#'   positives (true positives + false negatives).
#'
#' @param pred A vector of predictions that are either logical or numeric, with
#'   0 being false and 1 being true.
#' @param truth A vector of ground truths that are either logical or numeric, with
#'   0 being false and 1 being true.
#' @return The sensitivity as a single number.
#' @export
#' @examples
#' pred <- c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)
#' truth <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
#' sensitivity(pred, truth)  # returns 0.6
sensitivity = function(pred,truth) {
  if (class(pred) != class(truth)) {
    stop("Prediction and truth vectors must have the same class")
  }
  if (class(pred) != "logical" && class(pred) != "numeric") {
    stop("Prediction and truth vectors must be either logical or numeric")
  }
  true_pos = sum((pred == TRUE | pred == 1) & pred == truth)
  total_pos = sum((truth == TRUE | truth == 1))
  true_pos / total_pos
}


#' Calculate specificity as the ratio of true negatives to total ground truth
#'   negatives (true positives + false negatives).
#'
#' @param pred A vector of predictions that are either logical or numeric, with
#'   0 being false and 1 being true.
#' @param truth A vector of ground truths that are either logical or numeric, with
#'   0 being false and 1 being true.
#' @return The specificity as a single number.
#' @export
#' @examples
#' pred <- c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)
#' truth <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
#' specificity(pred, truth)  # returns 0.5
specificity = function(pred,truth) {
  if (class(pred) != class(truth)) {
    stop("Prediction and truth vectors must have the same class")
  }
  if (class(pred) != "logical" && class(pred) != "numeric") {
    stop("Prediction and truth vectors must be either logical or numeric")
  }
  true_neg = sum((pred == FALSE | pred == 0) & pred == truth)
  total_neg = sum((truth == FALSE | truth == 0))
  true_neg / total_neg
}


#' Calculate precision as the ratio of true positives to total positive predictions; also
#'   referred to as Positive Predictive Value (PPV).
#'
#' @param pred A vector of predictions that are either logical or numeric, with
#'   0 being false and 1 being true.
#' @param truth A vector of ground truths that are either logical or numeric, with
#'   0 being false and 1 being true.
#' @return The precision as a single number.
#' @export
#' @examples
#' pred1 <- c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)
#' truth1 <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
#' precision(pred1, truth1)  # returns 0.75
precision = function(pred,truth) {
  if (class(pred) != class(truth)) {
    stop("Prediction and truth vectors must have the same class")
  }
  if (class(pred) != "logical" && class(pred) != "numeric") {
    stop("Prediction and truth vectors must be either logical or numeric")
  }
  true_pos = sum((pred == TRUE | pred == 1) & pred == truth)
  total_pos_pred = sum((pred == TRUE | pred == 1))
  true_pos / total_pos_pred
}


#' Calculate recall as the ratio of true positives to total ground truth
#'   positives (true positives + false negatives).
#'
#' @param pred A vector of predictions that are either logical or numeric, with
#'   0 being false and 1 being true.
#' @param truth A vector of ground truths that are either logical or numeric, with
#'   0 being false and 1 being true.
#' @return The recall as a single number.
#' @export
#' @examples
#' pred1 <- c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)
#' truth1 <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
#' recall(pred1, truth1)  # returns 0.6
recall = function(pred,truth) {
  sensitivity(pred, truth)
}


#' Return the accuracy as the ratio of correct predictions to total predictions
#'
#' @param pred A vector of predictions that are either logical or numeric, with
#'   0 being false and 1 being true.
#' @param truth A vector of ground truths that are either logical or numeric, with
#'   0 being false and 1 being true.
#' @return The accuracy as a single number.
#' @export
#' @examples
#' pred1 <- c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)
#' truth1 <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
#' accuracy(pred1, truth1)  # returns 0.57
accuracy = function(pred,truth) {
  if (class(pred) != class(truth)) {
    stop("Prediction and truth vectors must have the same class")
  }
  if (class(pred) != "logical" && class(pred) != "numeric") {
    stop("Prediction and truth vectors must be either logical or numeric")
  }
  true_pos = sum((pred == TRUE | pred == 1) & pred == truth)
  true_neg = sum((pred == FALSE | pred == 0) & pred == truth)
  (true_pos + true_neg) / length(pred)
}


#' Calculate the F1 score as the harmonic mean of precision and recall
#'
#' @param pred A vector of predictions that are either logical or numeric, with
#'   0 being false and 1 being true.
#' @param truth A vector of ground truths that are either logical or numeric, with
#'   0 being false and 1 being true.
#' @return The F1 score as a single number.
#' @export
#' @examples
#' pred1 <- c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)
#' truth1 <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
#' f1(pred1, truth1)  # returns 0.67
f1 = function(pred,truth) {
  precision = precision(pred, truth)
  recall = recall(pred, truth)
  2 * precision * recall / (precision + recall)
}


#' Return the minimum sample size needed for a two-sample t-test with the given
#'   Cohen's d value and power.
#'
#' @param d true difference in means (Cohen's d)
#' @param power power of test (1 minus Type II error probability)
#' @return The minimum sample size needed for a two-sample t-test with the given
#'   Cohen's d value and power level.
#' @export
#' @examples
#' minimum_n_per_group(1)  # returns 17
#' minimum_n_per_group(d = 0.01, power = 0.9)  # returns 210150
minimum_n_per_group = function(d,power = 0.8) {
  ceiling(stats::power.t.test(delta = d, power = power)$n)
}


#' Calculated the r-squared statistic between predicted and ground truth values.
#'
#' @param pred The predicted values (on a continuous scale).
#' @param truth The ground truth values (on a continuous scale).
#' @return The r-squared statistic between predicted and ground truth values.
#' @export
#' @examples
#' a <- c(6, 2, 7, 3, 9)
#' b <- c(3, 1, 7, 3, 6)
#' r2(pred = a, truth = b)  # returns 0.72
r2 = function(pred,truth) {
  model = stats::lm(pred ~ truth)
  summary(model)$r.squared
}


#' Calculate the adjusted r-squared statistic between predicted and ground truth
#'  values with a given number of parameterss n_p.
#'
#' @param pred The predicted values (on a continuous scale).
#' @param truth The ground truth values (on a continuous scale).
#' @param n_p The number of model parameters.
#' @return The adjusted r-squared statistic.
#' @export
#' @examples
#' a <- c(6, 2, 7, 3, 9)
#' b <- c(3, 1, 7, 3, 6)
#' adj_R2(pred = a, truth = b, 3)  # returns -0.11
adj_R2 = function(pred,truth,n_p) {
  n_obs = length(pred)
  1 - (1 - r2(pred=pred, truth=truth)) * (n_obs - 1) / (n_obs - n_p - 1)
}
