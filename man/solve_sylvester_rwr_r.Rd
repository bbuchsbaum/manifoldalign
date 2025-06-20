% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parrot.R
\name{solve_sylvester_rwr_r}
\alias{solve_sylvester_rwr_r}
\title{Solve Sylvester Equation for Cross-Graph RWR Cost (R implementation)}
\usage{
solve_sylvester_rwr_r(
  W1,
  W2T,
  Cnode,
  beta = 0.15,
  gamma = 0.1,
  tol = 1e-06,
  max_iter = 50
)
}
\arguments{
\item{W1}{First network's transition matrix (row-stochastic)}

\item{W2T}{Transpose of second network's transition matrix}

\item{Cnode}{Node feature cost matrix}

\item{beta}{RWR restart probability (corresponds to sigma in code)}

\item{gamma}{Cross-graph discount factor}

\item{tol}{Convergence tolerance}

\item{max_iter}{Maximum iterations}
}
\value{
Cross-graph RWR cost matrix C_rwr
}
\description{
Implements the Sylvester iteration from Eq. 4 in the paper:
C_rwr = (1+β)C_node + (1-β)γ * W1 * C_rwr * W2^T
}
\keyword{internal}
