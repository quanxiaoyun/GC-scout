% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Sigma_witten}
\alias{Sigma_witten}
\title{Data from Experiment 1}
\format{
A 9x9 correlation matrix.
}
\usage{
Sigma_witten
}
\description{
This is a correlation matrix that was first introduced in a numeric study
Witten & Tibshirani (2009) for similar covariance-regularized regression
problems, see Simulation 6 in Witten & Tibshirani (2009). After appending to
the column the correlation between \code{x} and
\code{y} followed by an extra step of normalization, we made this matrix
into a 9x9 correlation matrix. This is exactly the same matrix we used for
Simulation 5 in our paper (Quan et al (2020)).
}
\keyword{datasets}
