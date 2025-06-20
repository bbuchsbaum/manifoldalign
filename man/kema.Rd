% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R
\name{kema}
\alias{kema}
\title{Kernel Manifold Alignment (KEMA)}
\usage{
kema(data, y, ...)
}
\arguments{
\item{data}{Input data object. Can be a hyperdesign object (for 
\code{kema.hyperdesign}) or a multidesign object (for 
\code{kema.multidesign})}

\item{y}{Name of the label variable to use for alignment. Can contain NA 
values for unlabeled samples in semi-supervised learning scenarios}

\item{...}{Additional arguments passed to specific methods. See 
\code{\link{kema.hyperdesign}} and \code{\link{kema.multidesign}} for 
details on method-specific parameters such as \code{preproc}, 
\code{ncomp}, \code{knn}, \code{sigma}, \code{u}, \code{kernel}, 
\code{sample_frac}, \code{use_laplacian}, \code{solver}, \code{dweight}, 
\code{rweight}, \code{simfun}, \code{disfun}, \code{lambda}, and for 
\code{kema.multidesign} the additional \code{subject} parameter}
}
\value{
A \code{multiblock_biprojector} object containing:
\itemize{
  \item \code{s}: Scores (embedded coordinates) for all samples
  \item \code{v}: Primal vectors (feature weights) for out-of-sample 
    projection
  \item \code{sdev}: Standard deviations of the components
  \item \code{alpha}: Dual coefficients in kernel space
  \item Additional metadata for reconstruction and validation
}
}
\description{
Performs Kernel Manifold Alignment for supervised/semi-supervised domain 
adaptation. Projects data from multiple domains into a shared latent space.
}
\details{
KEMA is designed for multi-domain data where you want to find a common 
representation that preserves both the intrinsic geometry of each domain and 
the class structure across domains. It supports semi-supervised learning 
with missing labels (NA values).

The algorithm offers two solver methods:
- "regression": Fast approximation using spectral regression (default). This 
  method first solves the eigenvalue problem on graph Laplacians, then uses 
  ridge regression to find kernel coefficients. It's much faster but may be 
  less accurate for non-linear kernels.
- "exact": Precise solution using the correct generalized eigenvalue 
  formulation. This method solves the mathematically correct KEMA 
  optimization problem but is more computationally intensive, especially for 
  large datasets.


KEMA solves the following optimization problem:
\deqn{\max_{\alpha} \frac{\alpha^T K A K^T \alpha}{\alpha^T K B K^T \alpha}}

where:
\itemize{
  \item \code{A = u*L + (1-u)*Ls} captures manifold structure (L) and 
    same-class alignment (Ls)
  \item \code{B = rweight*Lr + dweight*Ld + lambda*I} captures class 
    separation and regularization
  \item \code{K} is the block-diagonal kernel matrix across domains
}

The trade-off parameter \code{u} controls the balance between preserving 
manifold geometry and enforcing class alignment. The solver parameter 
determines the computational approach:
\itemize{
  \item "regression": Fast two-step approximation (eigendecomposition + 
    ridge regression)
  \item "exact": Direct solution of the generalized eigenvalue problem
}

For large datasets, use \code{sample_frac < 1} to enable REKEMA, which uses 
landmark points to reduce computational complexity from O(n^2) to O(r^2) where 
r is the number of landmarks.
}
\examples{
\donttest{
# Example with hyperdesign data
library(multivarious)

# Create synthetic multi-domain data
domain1 <- list(x = matrix(rnorm(100), 50, 2), 
               labels = sample(c("A", "B"), 50, TRUE))
domain2 <- list(x = matrix(rnorm(100), 50, 2), 
               labels = sample(c("A", "B"), 50, TRUE))
hd <- list(domain1 = domain1, domain2 = domain2)

# Run KEMA with default settings
result <- kema(hd, y = labels, ncomp = 2)

# Semi-supervised learning with missing labels
domain1$labels[1:10] <- NA  # Mark some samples as unlabeled
result_semi <- kema(hd, y = labels, ncomp = 2)

# Use exact solver for highest accuracy
result_exact <- kema(hd, y = labels, solver = "exact", ncomp = 2)

# Use REKEMA for large datasets
result_rekema <- kema(hd, y = labels, sample_frac = 0.5, ncomp = 2)
}

}
\references{
Tuia, D., & Camps-Valls, G. (2016). Kernel manifold alignment for domain 
adaptation. PLoS ONE, 11(2), e0148655.
}
\seealso{
\code{\link{kema.hyperdesign}}, \code{\link{kema.multidesign}}
}
