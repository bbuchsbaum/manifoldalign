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
  backend = "auto",
  backend_control = NULL,
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
  backend = "auto",
  backend_control = NULL,
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

  Deprecated compatibility argument; ignored.

- solver:

  Deprecated compatibility argument; accepted values are
  \`"regression"\` and \`"exact"\`, but both currently route to the
  original KEMA solver.

- backend:

  Backend for the original eigensolver. One of \`"auto"\`,
  \`"full_exact"\`, \`"reduced_exact"\`, or \`"operator_exact"\`.

- backend_control:

  Optional list controlling auto backend thresholds and fidelity checks
  (passed through to \`kema_orig()\`).

- dweight:

  Deprecated compatibility argument; ignored.

- rweight:

  Deprecated compatibility argument; ignored.

- simfun:

  Deprecated compatibility argument; ignored.

- disfun:

  Deprecated compatibility argument; ignored.

- lambda:

  Regularization parameter for matrix conditioning (default: 0.0001)

- centre_kernel:

  Deprecated compatibility argument; ignored.

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

Current behavior routes \`kema()\` to a paper-faithful implementation
(\`kema_orig\`) of the original Tuia & Camps-Valls generalized
eigenproblems. Legacy extension arguments are still accepted for
compatibility.

KEMA solves the original paper objective: \$\$K(L+\mu L_s)K\Lambda =
\lambda K L_d K\Lambda\$\$ and its reduced-rank REKEMA counterpart when
\`sample_frac \< 1\`.

\`kema()\` now delegates to the paper-faithful \`kema_orig()\` backend
and solves the original generalized eigenproblems from Tuia &
Camps-Valls (2016), including the reduced-rank REKEMA form when
\`sample_frac \< 1\`.

Legacy extension arguments remain in the API for backward compatibility
but are ignored by the current implementation.

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
#> Warning: KEMA fidelity checks failed for backend 'full_exact': max_rel_residual=1, max_B_orth_offdiag=5.63e-11

# Semi-supervised learning with missing labels
design1$labels[1:4] <- NA  # Mark a few samples as unlabeled
md1_semi <- multidesign(X1, design1)
hd_semi <- hyperdesign(list(domain1 = md1_semi, domain2 = md2))
result_semi <- kema(hd_semi, y = labels, ncomp = 2)
#> Warning: KEMA fidelity checks failed for backend 'full_exact': max_rel_residual=1, max_B_orth_offdiag=6.7e-11

# Use exact solver for highest accuracy
result_exact <- kema(hd, y = labels, solver = "exact", ncomp = 2)
#> Warning: KEMA fidelity checks failed for backend 'full_exact': max_rel_residual=0.0252, max_B_orth_offdiag=1.15e-10

# Use REKEMA for large datasets
result_rekema <- kema(hd, y = labels, sample_frac = 0.5, ncomp = 2)
#> REKEMA block 1: 20 x 10 kernel matrix
#> REKEMA block 2: 20 x 10 kernel matrix
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
#> Warning: KEMA fidelity checks failed for backend 'full_exact': max_rel_residual=0.00307, max_B_orth_offdiag=1.01e-12

# Semi-supervised learning with missing labels
data_design$condition[sample(nrow(data_design), 20)] <- NA
md_semi <- multidesign(X, data_design)
result_semi <- kema(md_semi, y = condition, subject = subject, ncomp = 2)
#> Warning: KEMA fidelity checks failed for backend 'full_exact': max_rel_residual=0.0353, max_B_orth_offdiag=3.61e-13
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
#> Warning: KEMA fidelity checks failed for backend 'full_exact': max_rel_residual=0.986, max_B_orth_offdiag=3.47e-10

# Semi-supervised learning with missing labels
hd_semi <- hd
hd_semi$domain1$design$labels[1:10] <- NA  # Mark some samples as unlabeled
result_semi <- kema(hd_semi, y = labels, ncomp = 2)
#> Warning: KEMA fidelity checks failed for backend 'full_exact': max_rel_residual=0.907, max_B_orth_offdiag=5.75e-10

# Use exact solver for highest accuracy
result_exact <- kema(hd, y = labels, solver = "exact", ncomp = 2)
#> Warning: KEMA fidelity checks failed for backend 'full_exact': max_rel_residual=0.986, max_B_orth_offdiag=3.47e-10

# Use REKEMA for large datasets
result_rekema <- kema(hd, y = labels, sample_frac = 0.5, ncomp = 2)
#> REKEMA block 1: 50 x 25 kernel matrix
#> REKEMA block 2: 50 x 25 kernel matrix
# }
```
