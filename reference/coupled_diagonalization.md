# Coupled Diagonalization for Multi-Modal Alignment

Performs Coupled Diagonalization alignment on multi-modal data. Finds
coupled eigenvectors across multiple modalities that jointly diagonalize
graph Laplacians while respecting cross-modal correspondences.

Performs Coupled Diagonalization alignment on hyperdesign data
structures containing multiple modalities. Finds coupled eigenvectors
that jointly diagonalize graph Laplacians while respecting cross-modal
correspondences.

## Usage

``` r
coupled_diagonalization(data, ...)

# S3 method for class 'hyperdesign'
coupled_diagonalization(
  data,
  correspondence = NULL,
  preproc = center(),
  ncomp = 10,
  ncomp_per_domain = 20,
  mu_coupling = 1,
  alpha_match = 0,
  lambda_target = NULL,
  knn = 10,
  sigma = NULL,
  spectral_cache = NULL,
  max_iter = 200,
  step_size = 0.3,
  tol = 1e-06,
  verbose = FALSE,
  ...
)

# Default S3 method
coupled_diagonalization(data, ...)
```

## Arguments

- data:

  A hyperdesign object containing multiple data domains

- ...:

  Additional arguments (currently unused)

- correspondence:

  Optional list of correspondence matrices between modalities. If NULL
  (default), assumes full correspondence based on sample order. Each
  matrix should be n_i x n_j sparse matrix indicating correspondences
  between modality i and j.

- preproc:

  Preprocessing function to apply to the data (default: center())

- ncomp:

  Number of coupled eigenvectors to extract (default: 10)

- ncomp_per_domain:

  Number of eigenvectors per domain to compute before coupling (default:
  20). Must be \>= ncomp.

- mu_coupling:

  Coupling weight controlling alignment vs diagonalization trade-off
  (default: 1). Higher values emphasize cross-modal alignment.

- alpha_match:

  Optional weight for matching the diagonal of `A_i^T Lambda_i A_i` to
  target eigenvalues (default: 0).

- lambda_target:

  Optional list of length `m` providing target eigenvalue vectors for
  each modality (only the first `ncomp` entries are used).

- knn:

  Number of nearest neighbors for graph construction (default: 10)

- sigma:

  Gaussian kernel bandwidth for graph weights. If NULL (default), uses
  median of nearest neighbor distances.

- spectral_cache:

  Optional precomputed spectral cache to reuse Laplacian subspaces
  across repeated runs (for example
  `previous_fit$diagnostics$spectral_cache`). When supplied,
  graph/Laplacian/eigendecomposition steps are skipped.

- max_iter:

  Maximum iterations for Stiefel manifold optimization (default: 200)

- step_size:

  Base step size for gradient descent (default: 0.3)

- tol:

  Convergence tolerance for relative cost change (default: 1e-6)

- verbose:

  Whether to print convergence information (default: FALSE)

## Value

A `multiblock_biprojector` object containing:

- `s`: Coupled eigenvectors (scores) for all samples

- `v`: Projection matrices for out-of-sample extension

- `coupled_bases`: List of coupled bases for each modality

- `sdev`: Standard deviations of the coupled components

- Additional metadata for reconstruction and diagnostics

A multiblock_biprojector object containing:

- `s`: Coupled eigenvectors (scores) concatenated across domains

- `v`: Projection matrices for out-of-sample extension

- `coupled_bases`: List of coupled bases V_i for each modality

- `sdev`: Standard deviations of coupled components

- `eigenvalues`: List of eigenvalues for each domain

- `diagnostics`: List with per-iteration cost diagnostics

- `converged`: Logical indicating convergence

- `iterations`: Number of iterations performed

- `final_cost`: Final objective value

## Details

Coupled Diagonalization (Eynard et al. 2015) extends spectral methods to
multi-modal data by finding bases that simultaneously diagonalize
individual Laplacians while being aligned through correspondence
constraints. The algorithm uses Stiefel manifold optimization to
maintain orthogonality constraints.

The algorithm minimizes the objective function: \$\$F(A) = sum_i
\|\|A_i^T L_i A_i - diag(A_i^T L_i A_i)\|\|\_F^2 + mu_c sum\_{i\<j}
\|\|F_i^T U_i A_i - F_j^T U_j A_j\|\|\_F^2\$\$

where:

- `A_i` are the coupling matrices on the Stiefel manifold

- `L_i` are the eigenvalues of the i-th Laplacian

- `U_i` are the eigenvectors of the i-th Laplacian

- `F_i` are the correspondence matrices between modalities

- `mu_c` controls the coupling strength

The first term encourages diagonalization of individual Laplacians,
while the second term enforces alignment through correspondences. The
algorithm uses projected gradient descent on the Stiefel manifold to
maintain orthogonality.

Key features:

- Handles partial correspondences between modalities

- Preserves manifold structure through graph Laplacians

- Maintains orthogonality via Stiefel manifold optimization

- Supports semi-supervised learning with missing correspondences

The algorithm minimizes the objective

    F = sum_i ||A_i^T Lambda_i A_i - diag(A_i^T Lambda_i A_i)||_F^2 +
        mu_c sum_{i<j} ||F_i^T U_i A_i - F_j^T U_j A_j||_F^2

It uses projected gradient descent on the Stiefel manifold to maintain
orthogonality of the coupling matrices A_i. When `alpha_match > 0`, the
optional match term nudges the diagonal entries of `A_i^T Lambda_i A_i`
toward user-specified target eigenvalues (without adding extra
off-diagonal penalty).

## References

Eynard, D., Kovnatsky, A., Bronstein, M. M., Glashoff, K., & Bronstein,
A. M. (2015). Multimodal manifold analysis by simultaneous
diagonalization of Laplacians. IEEE Transactions on Pattern Analysis and
Machine Intelligence, 37(12), 2505-2517.

Eynard et al. (2015). Multimodal manifold analysis by simultaneous
diagonalization of Laplacians. IEEE TPAMI, 37(12), 2505-2517.

## See also

`coupled_diagonalization.hyperdesign`

## Examples

``` r
# \donttest{
# Example with hyperdesign data
library(multidesign)

# Create synthetic multi-modal data
set.seed(123)
# Domain 1: 3D data
X1 <- matrix(rnorm(150), 50, 3)
# Domain 2: 4D data  
X2 <- matrix(rnorm(200), 50, 4)

# Create design with full correspondence
design1 <- data.frame(sample_id = 1:50)
design2 <- data.frame(sample_id = 1:50)

# Create multidesign objects
md1 <- multidesign(X1, design1)
md2 <- multidesign(X2, design2)

# Create hyperdesign
hd <- hyperdesign(list(domain1 = md1, domain2 = md2))

# Run coupled diagonalization with default settings
result <- coupled_diagonalization(hd, ncomp = 3)
#> sigma is 1.53794640409146
#> sigma is 1.85511070143682

# Access coupled bases
coupled_scores <- result$s
coupled_bases <- result$coupled_bases

# With partial correspondences
# Specify correspondence matrix (e.g., only first 30 samples correspond)
F_partial <- Matrix::sparseMatrix(i = 1:30, j = 1:30, x = 1, dims = c(50, 50))
result_partial <- coupled_diagonalization(hd, 
                                        correspondence = list(F_partial, F_partial),
                                        ncomp = 3)
#> sigma is 1.53794640409146
#> sigma is 1.85511070143682
# }

# \donttest{
library(multidesign)
# Create synthetic multi-modal data
X1 <- matrix(rnorm(150), 50, 3)
X2 <- matrix(rnorm(200), 50, 4)

# Create hyperdesign
hd <- hyperdesign(list(
  domain1 = multidesign(X1, data.frame(id = 1:50)),
  domain2 = multidesign(X2, data.frame(id = 1:50))
))

# Run coupled diagonalization
result <- coupled_diagonalization(hd, ncomp = 3)
#> sigma is 1.48550141814904
#> sigma is 1.96203797806432
# }
```
