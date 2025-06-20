% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/cone_align.R
\name{cone_align}
\alias{cone_align}
\title{CONE-Align: Consensus Optimization for Node Embedding Alignment}
\usage{
cone_align(data, ...)

cone_align(data, ...)
}
\arguments{
\item{data}{Input data object containing graph domains}

\item{...}{Additional arguments passed to specific methods. See 
\code{\link{cone_align.hyperdesign}} for details on method-specific 
parameters such as \code{preproc}, \code{ncomp}, \code{sigma}, 
\code{lambda}, \code{use_laplacian}, \code{solver}, \code{max_iter}, 
and \code{tol}}
}
\value{
A \code{multiblock_biprojector} object containing:
\itemize{
  \item \code{s}: Aligned embeddings for all nodes
  \item \code{v}: Primal vectors for out-of-sample projection
  \item \code{assignment}: Node-to-node correspondence assignments
  \item \code{rotation}: Orthogonal transformation matrices
  \item \code{sdev}: Standard deviations of aligned components
  \item Additional metadata for reconstruction and validation
}

A \code{multiblock_biprojector} object containing:
\itemize{
  \item \code{s}: Aligned embeddings for all nodes
  \item \code{v}: Primal vectors for out-of-sample projection
  \item \code{assignment}: Node-to-node correspondence assignments
  \item \code{rotation}: Orthogonal transformation matrices
  \item \code{sdev}: Standard deviations of aligned components
  \item Additional metadata for reconstruction and validation
}
}
\description{
Performs CONE-Align on hyperdesign data structures. Aligns graph embeddings
by alternating between orthogonal transformations and node assignments.

Performs CONE-Align on hyperdesign data structures. Aligns graph embeddings
by alternating between orthogonal transformations and node assignments.
}
\details{
CONE-Align tackles graph alignment by iteratively refining both orthogonal
transformations and node-to-node correspondences. The algorithm alternates
between two steps: (1) finding optimal orthogonal rotations given current
assignments via Procrustes analysis, and (2) updating assignments given
current transformations via linear assignment.


CONE-Align operates through the following algorithmic blocks:
\itemize{
  \item \strong{Spectral Embedding}: Compute Laplacian eigenmaps for each graph
  \item \strong{Orthogonal Alignment}: Find rotation matrices via Procrustes analysis
  \item \strong{Assignment Update}: Solve linear assignment problem for correspondences
  \item \strong{Convergence Check}: Iterate until assignments stabilize
}

The algorithm minimizes the objective:
\deqn{\sum_{i,j} ||Q_i^T Z_i - P_{ij} Q_j^T Z_j||_F^2}

where \eqn{Z_i} are the embeddings, \eqn{Q_i} are orthogonal transforms,
and \eqn{P_{ij}} are permutation matrices.

Key parameters:
\itemize{
  \item \code{ncomp}: Dimension of spectral embeddings
  \item \code{sigma}: Bandwidth for embedding computation
  \item \code{lambda}: Regularization for numerical stability
  \item \code{solver}: Assignment algorithm ("linear" or "auction")
}

CONE-Align tackles graph alignment by iteratively refining both orthogonal
transformations and node-to-node correspondences. The algorithm alternates
between two steps: (1) finding optimal orthogonal rotations given current
assignments via Procrustes analysis, and (2) updating assignments given
current transformations via linear assignment.


CONE-Align operates through the following algorithmic blocks:
\itemize{
  \item \strong{Spectral Embedding}: Compute Laplacian eigenmaps for each graph
  \item \strong{Orthogonal Alignment}: Find rotation matrices via Procrustes analysis
  \item \strong{Assignment Update}: Solve linear assignment problem for correspondences
  \item \strong{Convergence Check}: Iterate until assignments stabilize
}

The algorithm minimizes the objective:
\deqn{\sum_{i,j} ||Q_i^T Z_i - P_{ij} Q_j^T Z_j||_F^2}

where \eqn{Z_i} are the embeddings, \eqn{Q_i} are orthogonal transforms,
and \eqn{P_{ij}} are permutation matrices.

Currently supports *exactly two* domains. For multi-graph alignment with
three or more domains, see \code{\link{cone_align_multiple}}.

Key parameters:
\itemize{
  \item \code{ncomp}: Dimension of spectral embeddings
  \item \code{sigma}: Bandwidth for embedding computation
  \item \code{lambda}: Regularization for numerical stability
  \item \code{solver}: Assignment algorithm ("linear" or "auction")
  \item \code{knn}: Number of nearest neighbors for graph construction
}
}
\examples{
\donttest{
# Example with hyperdesign graph data
library(multidesign)

# Create synthetic graph domains
set.seed(123)
domain1 <- list(
  x = matrix(rnorm(100), 50, 2),  # Node features
  design = data.frame(node_id = 1:50)
)
domain2 <- list(
  x = matrix(rnorm(100), 50, 2),
  design = data.frame(node_id = 1:50)  
)

# Create hyperdesign
hd <- hyperdesign(list(domain1 = domain1, domain2 = domain2))

# Run CONE-Align with default parameters
result <- cone_align(hd, ncomp = 10)

# Access alignment results
node_assignment <- result$assignment
aligned_embeddings <- result$s

# Use exact solver for optimal results
result_exact <- cone_align(hd, ncomp = 10, solver = "linear")
}

\donttest{
# Example with hyperdesign graph data
library(multidesign)
library(tibble)

# Create synthetic graph domains (node coordinates for graph construction)
set.seed(123)
X1 <- matrix(rnorm(100), 50, 2)  # Node coordinates for domain 1
X2 <- matrix(rnorm(100), 50, 2)  # Node coordinates for domain 2

# Create design data frames with node indices (CONE-Align doesn't use features)
design1 <- tibble(node_id = 1:50)
design2 <- tibble(node_id = 1:50)

# Create multidesign objects
md1 <- multidesign(X1, design1)
md2 <- multidesign(X2, design2)

# Create hyperdesign from multidesign objects
hd <- hyperdesign(list(domain1 = md1, domain2 = md2))

# Run CONE-Align (purely spectral, uses graph structure from X)
result <- cone_align(hd, ncomp = 10)

# Access alignment results
node_assignment <- result$assignment
aligned_embeddings <- result$s

# Use exact solver for optimal results
result_exact <- cone_align(hd, ncomp = 10, solver = "linear")
}

}
\references{
Heimann, M., Shen, H., Safavi, T., & Koutra, D. (2018). REGAL: Representation 
learning-based graph alignment. In Proceedings of the 27th ACM International 
Conference on Information and Knowledge Management (pp. 117-126).

Heimann, M., Shen, H., Safavi, T., & Koutra, D. (2018). REGAL: Representation 
learning-based graph alignment. In Proceedings of the 27th ACM International 
Conference on Information and Knowledge Management (pp. 117-126).
}
\seealso{
\code{\link{cone_align.hyperdesign}}

\code{\link{cone_align.hyperdesign}}
}
