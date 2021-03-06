% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fns.R
\name{cv.GC_scout}
\alias{cv.GC_scout}
\title{Perform cross-validation for the GC-Scout regression}
\usage{
cv.GC_scout(
  x,
  y,
  d2,
  lam1s = seq(0.01, 0.25, length.out = 10),
  lam2s = seq(0.01, 0.25, length.out = 10),
  K = 5,
  gc_scout_01 = TRUE,
  gc_scout_21 = TRUE,
  gc_scout_11 = TRUE,
  MAPE = FALSE
)
}
\arguments{
\item{x}{A numeric matrix. Used for training.}

\item{y}{A numeric vector. Used for training.}

\item{d2}{A numeric vector. The indices of discrete variables in \code{x}.}

\item{lam1s}{A numeric number or vector. The tuning parameters for
regularization of the covariance matrix. If
\code{gc_scout_01 = TRUE}, then no covariance regularization
is taking place and \code{lam1s} will be ignored. If
\code{gc_scout_11 = TRUE},  as graphical lasso can be uncomfortably slow.}

\item{lam2s}{A numeric number or vector. The tuning parameters
for the \eqn{L_1} regularization of the regression coefficients,
using the regularized covariance matrix.}

\item{K}{An integer. Number of cross-validation folds to be performed; default is 5.}

\item{gc_scout_01}{A logical value. If TRUE, GC_Scout(0,1) (aka GC Lasso) will be performed. Default is TRUE.}

\item{gc_scout_21}{A logical value. If TRUE, GC_Scout(2,1) will be performed. Default is TRUE.}

\item{gc_scout_11}{A logical value. If TRUE, GC_Scout(1,1) will be performed. Default is TRUE.}

\item{MAPE}{Logical. If TRUE, MAPE (mean absolute percentage error) will be
used to measure cross-validation error; if FALSE, MSE (Mean Squared Error)
will be used instead. Default is FALSE.}
}
\value{
A list with the following elements:
\itemize{
  \item mse_gc_scout_01 - the beta estimates for GC_Scout(0,1) model (only if \code{gc_scout_01 = TRUE}; NULL otherwise)
  \item mse_gc_scout_11 - the beta estimates for GC_Scout(2,1) model (only if \code{gc_scout_21 = TRUE}; NULL otherwise)
  \item mse_gc_scout_21 - the beta estimates for GC_Scout(1,1) model (only if \code{gc_scout_11 = TRUE}; NULL otherwise)
  \item bestlam1_11 - the best \eqn{\lamdba_1} associated with lowest cv error (only if \code{gc_scout_11 = TRUE}; NULL otherwise)
  \item bestlam1_21 - the best \eqn{\lamdba_1} associated with lowest cv error (only if \code{gc_scout_21 = TRUE}; NULL otherwise)
  \item bestlam2_01 - the best \eqn{\lamdba_2} associated with lowest cv error (only if \code{gc_scout_01 = TRUE}; NULL otherwise)
  \item bestlam2_11 - the best \eqn{\lamdba_2} associated with lowest cv error (only if \code{gc_scout_11 = TRUE}; NULL otherwise)
  \item bestlam2_21 - the best \eqn{\lamdba_2} associated with lowest cv error (only if \code{gc_scout_21 = TRUE}; NULL otherwise)
}
}
\description{
This function returns cross-validation errors for a range of tuning parameter
values (lam1 and lam2). Error measures can be specified as either MAPE (mean absolute
percentage error) or MSE (Mean Squared Error). Note that this function allows GC_Scout(0,1),
GC_Scout(1,1), and GC_Scout(2,1) models to be either fit individually or altogether in one run.
}
\examples{
data(Sigma_witten)
n <- 40; d <- ncol(Sigma_witten); d3 <- 4;
delta = c(qnorm(0.5), qnorm(0.75))
dat <- mvtnorm::rmvnorm(n, sigma = Sigma_witten)
y_tilde = dat[,d]
y = y_tilde; y = exp(y_tilde);
x_tilde = dat[,-d]; x = x_tilde^5
x[,1:d3] = as.numeric(cut(x[,1:d3], breaks = c(-Inf,(delta)^5,Inf)))-1
cv.GC_scout(x, y, lam1s = seq(0.01,0.2,length.out = 10), K=3,
 lam2s = seq(0.01,0.2,length.out = 10), d2 = seq(1:4),
 gc_scout_01 = FALSE, gc_scout_21 = TRUE, gc_scout_11 = FALSE)
}
\seealso{
\code{\link{GC_scout}}
}
