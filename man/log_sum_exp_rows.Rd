% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{log_sum_exp_rows}
\alias{log_sum_exp_rows}
\title{Log-sum-exp trick for rows}
\usage{
log_sum_exp_rows(log_K, log_v)
}
\arguments{
\item{log_K}{Log-kernel matrix (n x m)}

\item{log_v}{Log potential vector (m x 1)}
}
\value{
Log-sum-exp for each row (n x 1)
}
\description{
Computes log(sum(exp(log_v + log_K[i,]))) for each row i
using the log-sum-exp trick for numerical stability
}
