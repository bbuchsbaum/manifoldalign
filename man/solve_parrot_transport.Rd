% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parrot.R
\name{solve_parrot_transport}
\alias{solve_parrot_transport}
\title{Solve PARROT Transport Problem with Proximal Point Method}
\usage{
solve_parrot_transport(
  networks,
  rwr_features,
  anchor_info,
  lambda_e,
  lambda_n,
  lambda_p,
  tau,
  alpha,
  sigma,
  gamma,
  solver,
  max_iter,
  max_outer = 10,
  tol
)
}
\arguments{
\item{networks}{List of network structures}

\item{rwr_features}{List of RWR descriptor matrices}

\item{anchor_info}{List with anchor information for both networks}

\item{lambda_e}{Edge consistency weight}

\item{lambda_n}{Neighborhood consistency weight}

\item{lambda_p}{Anchor prior weight}

\item{tau}{Entropy regularization parameter}

\item{solver}{Transport solver method}

\item{max_iter}{Maximum iterations for inner Sinkhorn}

\item{max_outer}{Maximum outer proximal iterations}

\item{tol}{Convergence tolerance}
}
\value{
List with transport plan and related matrices
}
\description{
Implements the proximal point loop that iteratively refines the cost matrix
and solves the optimal transport problem. This is critical for convergence.
}
\keyword{internal}
