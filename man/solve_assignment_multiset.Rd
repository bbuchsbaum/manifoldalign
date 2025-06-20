% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grasp_multiset.R
\name{solve_assignment_multiset}
\alias{solve_assignment_multiset}
\title{Solve Assignment for Multi-Graph GRASP}
\usage{
solve_assignment_multiset(E1, E2, distance = "cosine", solver = "auction")
}
\arguments{
\item{E1}{First embedding matrix}

\item{E2}{Second embedding matrix}

\item{distance}{Distance metric: "cosine" or "euclidean"}

\item{solver}{Assignment solver: "auction" or "linear"}
}
\value{
Assignment vector
}
\description{
Computes node assignment between two embedded graphs using
cosine or Euclidean distance.
}
\keyword{internal}
