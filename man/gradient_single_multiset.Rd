% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grasp_multiset.R
\name{gradient_single_multiset}
\alias{gradient_single_multiset}
\title{Gradient for Single Rotation in Multi-Graph Setting}
\usage{
gradient_single_multiset(basis_s, FtPhi_s, FtPhi_all, Ms, s, lambda)
}
\arguments{
\item{basis_s}{Spectral basis for graph s}

\item{FtPhi_s}{Descriptor-basis product for graph s}

\item{FtPhi_all}{List of all descriptor-basis products}

\item{Ms}{List of all rotation matrices}

\item{s}{Index of current graph}

\item{lambda}{Regularization parameter}
}
\value{
Gradient matrix
}
\description{
Computes the gradient of the objective with respect to one rotation matrix
in the context of multiple graphs.
}
\keyword{internal}
