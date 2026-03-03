# Kernel Manifold Alignment (KEMA)

Performs Kernel Manifold Alignment for supervised/semi-supervised domain
adaptation. Projects data from multiple domains into a shared latent
space.

Performs Kernel Manifold Alignment on multidesign data structures.
Automatically splits data by subject variable and aligns domains.

Performs Kernel Manifold Alignment on hyperdesign data structures.
Projects data from multiple domains into a shared latent space while
preserving manifold structure and aligning same-class samples.

## Usage

``` r
kema(data, y, ...)

# S3 method for class 'multidesign'
kema(
  data,
  y,
  subject,
  preproc = center(),
  ncomp = 2,
  knn = 5,
  sigma = 0.73,
  u = 0.5,
  kernel = coskern(),
  sample_frac = 1,
  use_laplacian = TRUE,
  solver = "regression",
  dweight = 0.1,
  rweight = 0,
  simfun = neighborweights::binary_label_matrix,
  disfun = NULL,
  lambda = 1e-04,
  centre_kernel = FALSE,
  ...
)

# S3 method for class 'hyperdesign'
kema(
  data,
  y,
  preproc = center(),
  ncomp = 2,
  knn = 5,
  sigma = NULL,
  u = 0.5,
  kernel = NULL,
  sample_frac = 1,
  use_laplacian = TRUE,
  solver = "regression",
  dweight = 0.1,
  rweight = 0,
  simfun = neighborweights::binary_label_matrix,
  disfun = NULL,
  lambda = 1e-04,
  centre_kernel = FALSE,
  ...
)

# Default S3 method
kema(data, ...)
```

## Arguments

- data:

  A hyperdesign object containing multiple data domains

- y:

  Name of the label variable to use for alignment (can contain NA for
  unlabeled samples)

- ...:

  Additional arguments (currently unused)

- subject:

  Name of the subject variable that defines the domains/strata

- preproc:

  Preprocessing function to apply to the data (default: center())

- ncomp:

  Number of components to extract (default: 2)

- knn:

  Number of nearest neighbors for graph construction (default: 5)

- sigma:

  Kernel bandwidth parameter (default: 0.73)

- u:

  Trade-off parameter between data geometry and class alignment (0-1,
  default: 0.5)

- kernel:

  Kernel function to use (default: coskern())

- sample_frac:

  Fraction of samples to use for kernel approximation (default: 1)

- use_laplacian:

  Whether to use Laplacian normalization (default: TRUE)

- solver:

  Solver method: "regression" for fast approximation (default) or
  "exact" for precise solution

- dweight:

  Weight for dissimilarity/repulsion terms (default: 0.1)

- rweight:

  Weight for repulsion graph (default: 0)

- simfun:

  Function to compute similarity between labels

- disfun:

  Function to compute dissimilarity between labels (optional)

- lambda:

  Regularization parameter for matrix conditioning (default: 0.0001)

- centre_kernel:

  Whether to center kernel matrices (default: FALSE).
  \*\*\[EXTENSION\]\*\* The original paper uses uncentered kernels. Set
  TRUE for centered variant.

## Value

A `multiblock_biprojector` object containing:

- `s`: Scores (embedded coordinates) for all samples

- `v`: Primal vectors (feature weights) for out-of-sample projection

- `sdev`: Standard deviations of the components

- `alpha`: Dual coefficients in kernel space

- Additional metadata for reconstruction and validation

A multiblock_biprojector object containing the KEMA alignment

A multiblock_biprojector object containing the KEMA alignment

## Details

KEMA is designed for multi-domain data where you want to find a common
representation that preserves both the intrinsic geometry of each domain
and the class structure across domains. It supports semi-supervised
learning with missing labels (NA values).

The algorithm offers two solver methods: - "regression": Fast
approximation using spectral regression (default). This method first
solves the eigenvalue problem on graph Laplacians, then uses ridge
regression to find kernel coefficients. It's much faster but may be less
accurate for non-linear kernels. - "exact": Precise solution using the
correct generalized eigenvalue formulation. This method solves the
mathematically correct KEMA optimization problem but is more
computationally intensive, especially for large datasets.

KEMA solves the following optimization problem: \$\$\max\_{\alpha}
\frac{\alpha^T K A K^T \alpha}{\alpha^T K B K^T \alpha}\$\$

where:

- `A = u*L + (1-u)*Ls` captures manifold structure (L) and same-class
  alignment (Ls)

- `B = rweight*Lr + dweight*Ld + lambda*I` captures class separation and
  regularization

- `K` is the block-diagonal kernel matrix across domains

The trade-off parameter `u` controls the balance between preserving
manifold geometry and enforcing class alignment. The solver parameter
determines the computational approach:

- "regression": Fast two-step approximation (eigendecomposition + ridge
  regression)

- "exact": Direct solution of the generalized eigenvalue problem

For large datasets, use `sample_frac < 1` to enable REKEMA, which uses
landmark points to reduce computational complexity from O(n^2) to O(r^2)
where r is the number of landmarks.

This implementation follows the Tuia & Camps-Valls (2016) paper with
extensions:

\*\*Core KEMA (from paper):\*\* - Generalized eigenvalue problem:
Phi(L+mu\*Ls)Phi^T\*v = lambda \* Phi\*Ld\*Phi^T\*v (Eq. 4) -
Kernelization: K(L+mu\*Ls)K\*Lambda = lambda \* K\*Ld\*K \* Lambda (Eq.
6) - Reduced-rank KEMA (REKEMA) for computational efficiency -
Matrix-free eigensolver with Jacobi preconditioning

\*\*Extensions (not in original paper):\*\* - `solver="regression"`:
Fast spectral regression approximation (default) - `rweight`: Additional
repulsion graph Lr for within-domain separation - Semi-supervised
support: Handles NA labels for unlabeled samples - Enhanced numerical
stability and error handling

The algorithm offers two solver methods: - "regression":
\*\*\[EXTENSION\]\*\* Fast approximation using spectral regression
(default). This method first solves the eigenvalue problem on graph
Laplacians, then uses ridge regression to find kernel coefficients. Much
faster but may be less accurate for non-linear kernels. - "exact":
Precise solution using the correct generalized eigenvalue formulation
from the paper. This method solves the mathematically correct KEMA
optimization problem but is more computationally intensive, especially
for large datasets.

## References

Tuia, D., & Camps-Valls, G. (2016). Kernel manifold alignment for domain
adaptation. PLoS ONE, 11(2), e0148655.

Tuia, D., & Camps-Valls, G. (2016). Kernel manifold alignment for domain
adaptation. PLoS ONE, 11(2), e0148655.

## See also

`kema.hyperdesign`, `kema.multidesign`

## Examples

``` r
# \donttest{
# Example with hyperdesign data
library(multivarious)
#> 
#> Attaching package: ‘multivarious’
#> The following objects are masked from ‘package:manifoldalign’:
#> 
#>     apply_transform, block_indices
#> The following objects are masked from ‘package:stats’:
#> 
#>     residuals, screeplot
#> The following objects are masked from ‘package:base’:
#> 
#>     transform, truncate
library(multidesign)
library(tibble)

# Create synthetic multi-domain data
set.seed(123)
X1 <- matrix(rnorm(40), 20, 2)
X2 <- matrix(rnorm(40), 20, 2)
labels <- sample(c("A", "B"), 20, TRUE)

# Create design data frames
design1 <- data.frame(labels = labels)
design2 <- data.frame(labels = labels)

# Create multidesign objects
md1 <- multidesign(X1, design1)
md2 <- multidesign(X2, design2)

# Create hyperdesign
hd <- hyperdesign(list(domain1 = md1, domain2 = md2))

# Run KEMA with default settings
result <- kema(hd, y = labels, ncomp = 2, knn = 3)
#> Semi-supervised KEMA: 40 labeled samples, 0 unlabeled samples
#> Auto-selected sigma = 1.1158 using median distance heuristic
#> Using RBF kernel with auto-tuned sigma
#> normalize_laplacian(): 40 isolated nodes detected - treating as disconnected components.
#> Warning: Regression solver produced poor results (subspace angle: 19.6 deg, best match: 0.801). Automatically retrying with solver='exact' for higher fidelity.
#> Retry with exact solver completed successfully.

# Semi-supervised learning with missing labels
design1$labels[1:4] <- NA  # Mark a few samples as unlabeled
md1_semi <- multidesign(X1, design1)
hd_semi <- hyperdesign(list(domain1 = md1_semi, domain2 = md2))
result_semi <- kema(hd_semi, y = labels, ncomp = 2)
#> Semi-supervised KEMA: 36 labeled samples, 4 unlabeled samples
#> Auto-selected sigma = 1.1417 using median distance heuristic
#> Using RBF kernel with auto-tuned sigma
#> Warning: class_graph failed for stratum with labels containing NAs. Creating empty class graph. Error: Cannot create a graph object because the adjacency matrix contains NAs.
#> normalize_laplacian(): 4 isolated nodes detected - treating as disconnected components.
#> normalize_laplacian(): 40 isolated nodes detected - treating as disconnected components.
#> Warning: Regression solver produced poor results (subspace angle: 14.9 deg, best match: 0.95). Automatically retrying with solver='exact' for higher fidelity.
#> Retry with exact solver completed successfully.

# Use exact solver for highest accuracy
result_exact <- kema(hd, y = labels, solver = "exact", ncomp = 2)
#> Semi-supervised KEMA: 40 labeled samples, 0 unlabeled samples
#> Auto-selected sigma = 1.1558 using median distance heuristic
#> Using RBF kernel with auto-tuned sigma
#> normalize_laplacian(): 40 isolated nodes detected - treating as disconnected components.

# Use REKEMA for large datasets
result_rekema <- kema(hd, y = labels, sample_frac = 0.5, ncomp = 2)
#> Semi-supervised KEMA: 40 labeled samples, 0 unlabeled samples
#> Auto-selected sigma = 1.1287 using median distance heuristic
#> Using RBF kernel with auto-tuned sigma
#> REKEMA block 1: 20 x 10 kernel matrix
#> REKEMA block 2: 20 x 10 kernel matrix
#> normalize_laplacian(): 40 isolated nodes detected - treating as disconnected components.
#> Warning: REKEMA regression approximation fidelity is below target:
#>   Best component match: 0.479
#>   Subspace angle (deg): 40.4
#>   Samples: 40, Kernel rank: 20
#> Consider using solver='exact' or increasing sample_frac.
#> Warning: Regression solver produced poor results (subspace angle: 40.4 deg, best match: 0.479). Automatically retrying with solver='exact' for higher fidelity.
#> REKEMA: Solving reduced 20 x 20 eigenvalue problem instead of 40 x 40
#> REKEMA: Scaling lambda from 1e-04 to 2e-04 (factor: 2) to maintain regularization energy
#> REKEMA: Added regularization to ill-conditioned B matrix
#> Retry with exact solver completed successfully.
# }

# \donttest{
# Example with multidesign data
library(multidesign)

# Create synthetic multi-subject data
set.seed(123)
data_design <- expand.grid(
  subject = factor(1:4),
  condition = factor(c("A", "B")),
  trial = 1:10
)

# Generate synthetic data matrix
n_obs <- nrow(data_design)
n_features <- 20
X <- matrix(rnorm(n_obs * n_features), n_obs, n_features)

# Create multidesign object
md <- multidesign(X, data_design)

# Run KEMA alignment across subjects
result <- kema(md, y = condition, subject = subject, ncomp = 2)
#> Semi-supervised KEMA: 80 labeled samples, 0 unlabeled samples
#> normalize_laplacian(): 80 isolated nodes detected - treating as disconnected components.
#> Warning: Regression solver produced poor results (subspace angle: 80.2 deg, best match: 0.602). Automatically retrying with solver='exact' for higher fidelity.
#> Retry with exact solver completed successfully.

# Semi-supervised learning with missing labels
data_design$condition[sample(nrow(data_design), 20)] <- NA
md_semi <- multidesign(X, data_design)
result_semi <- kema(md_semi, y = condition, subject = subject, ncomp = 2)
#> Semi-supervised KEMA: 60 labeled samples, 20 unlabeled samples
#> Warning: class_graph failed for stratum with labels containing NAs. Creating empty class graph. Error: Cannot create a graph object because the adjacency matrix contains NAs.
#> Warning: class_graph failed for stratum with labels containing NAs. Creating empty class graph. Error: Cannot create a graph object because the adjacency matrix contains NAs.
#> Warning: class_graph failed for stratum with labels containing NAs. Creating empty class graph. Error: Cannot create a graph object because the adjacency matrix contains NAs.
#> Warning: class_graph failed for stratum with labels containing NAs. Creating empty class graph. Error: Cannot create a graph object because the adjacency matrix contains NAs.
#> normalize_laplacian(): 20 isolated nodes detected - treating as disconnected components.
#> normalize_laplacian(): 80 isolated nodes detected - treating as disconnected components.
#> Warning: Regression solver produced poor results (subspace angle: 71 deg, best match: 0.875). Automatically retrying with solver='exact' for higher fidelity.
#> Retry with exact solver completed successfully.
# }

# \donttest{
# Example with hyperdesign data
# Create synthetic multi-domain data
set.seed(123)
domain1 <- list(
  x = matrix(rnorm(100), 50, 2),
  design = data.frame(labels = sample(c("A", "B"), 50, TRUE))
)
domain2 <- list(
  x = matrix(rnorm(100), 50, 2),
  design = data.frame(labels = sample(c("A", "B"), 50, TRUE))
)
hd <- structure(list(domain1 = domain1, domain2 = domain2), class = "hyperdesign")

# Run KEMA with default settings
result <- kema(hd, y = labels, ncomp = 2)
#> Semi-supervised KEMA: 100 labeled samples, 0 unlabeled samples
#> Auto-selected sigma = 1.1446 using median distance heuristic
#> Using RBF kernel with auto-tuned sigma
#> normalize_laplacian(): 100 isolated nodes detected - treating as disconnected components.
#> Warning: Regression solver produced poor results (subspace angle: 29 deg, best match: 0.797). Automatically retrying with solver='exact' for higher fidelity.
#> Retry with exact solver completed successfully.

# Semi-supervised learning with missing labels
hd_semi <- hd
hd_semi$domain1$design$labels[1:10] <- NA  # Mark some samples as unlabeled
result_semi <- kema(hd_semi, y = labels, ncomp = 2)
#> Semi-supervised KEMA: 90 labeled samples, 10 unlabeled samples
#> Auto-selected sigma = 1.0682 using median distance heuristic
#> Using RBF kernel with auto-tuned sigma
#> Warning: class_graph failed for stratum with labels containing NAs. Creating empty class graph. Error: Cannot create a graph object because the adjacency matrix contains NAs.
#> normalize_laplacian(): 10 isolated nodes detected - treating as disconnected components.
#> normalize_laplacian(): 100 isolated nodes detected - treating as disconnected components.
#> Warning: Regression solver produced poor results (subspace angle: 31.3 deg, best match: 0.863). Automatically retrying with solver='exact' for higher fidelity.
#> Retry with exact solver completed successfully.

# Use exact solver for highest accuracy
result_exact <- kema(hd, y = labels, solver = "exact", ncomp = 2)
#> Semi-supervised KEMA: 100 labeled samples, 0 unlabeled samples
#> Auto-selected sigma = 1.236 using median distance heuristic
#> Using RBF kernel with auto-tuned sigma
#> normalize_laplacian(): 100 isolated nodes detected - treating as disconnected components.

# Use REKEMA for large datasets
result_rekema <- kema(hd, y = labels, sample_frac = 0.5, ncomp = 2)
#> Semi-supervised KEMA: 100 labeled samples, 0 unlabeled samples
#> Auto-selected sigma = 1.1201 using median distance heuristic
#> Using RBF kernel with auto-tuned sigma
#> REKEMA block 1: 50 x 25 kernel matrix
#> REKEMA block 2: 50 x 25 kernel matrix
#> normalize_laplacian(): 100 isolated nodes detected - treating as disconnected components.
#> Warning: REKEMA regression approximation fidelity is below target:
#>   Best component match: 0.512
#>   Subspace angle (deg): 37.9
#>   Samples: 100, Kernel rank: 50
#> Consider using solver='exact' or increasing sample_frac.
#> Warning: Regression solver produced poor results (subspace angle: 37.9 deg, best match: 0.512). Automatically retrying with solver='exact' for higher fidelity.
#> REKEMA: Solving reduced 50 x 50 eigenvalue problem instead of 100 x 100
#> REKEMA: Scaling lambda from 1e-04 to 2e-04 (factor: 2) to maintain regularization energy
#> REKEMA: Added regularization to ill-conditioned B matrix
#> Retry with exact solver completed successfully.
# }
```
