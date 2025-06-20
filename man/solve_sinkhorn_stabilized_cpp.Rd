% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{solve_sinkhorn_stabilized_cpp}
\alias{solve_sinkhorn_stabilized_cpp}
\title{Stabilized Sinkhorn algorithm in C++}
\usage{
solve_sinkhorn_stabilized_cpp(C, tau, max_iter, tol)
}
\arguments{
\item{C}{Cost matrix (n x m)}

\item{tau}{Entropy regularization parameter}

\item{max_iter}{Maximum number of iterations}

\item{tol}{Convergence tolerance}
}
\value{
Transport plan matrix S
}
\description{
Implements the log-domain stabilized Sinkhorn algorithm for optimal transport
}
