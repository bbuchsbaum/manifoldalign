# Position-Aware Random Transport (PARROT) Network Alignment

Performs PARROT alignment on hyperdesign data structures. Aligns
networks using regularized optimal transport with position-aware
features and consistency constraints.

Performs PARROT alignment on hyperdesign data structures. Aligns
networks using regularized optimal transport with position-aware
features and consistency constraints.

Performs PARROT alignment on hyperdesign data structures containing
network domains.

## Usage

``` r
parrot(data, anchors, ...)

parrot(data, anchors, ...)

# S3 method for class 'hyperdesign'
parrot(
  data,
  anchors,
  preproc = center(),
  ncomp = NULL,
  sigma = 0.15,
  lambda = 0.1,
  lambda_e = NULL,
  lambda_n = NULL,
  lambda_p = NULL,
  tau = 0.05,
  alpha = 0.2,
  gamma = 0.1,
  solver = c("sinkhorn", "exact"),
  max_iter = 100,
  tol = 1e-06,
  use_cpp = FALSE,
  verbose = FALSE,
  ...
)

# Default S3 method
parrot(data, anchors, ...)
```

## Arguments

- data:

  A hyperdesign object containing multiple network domains

- anchors:

  Name of the anchor/correspondence variable for semi-supervised
  alignment

- ...:

  Additional arguments (currently unused)

- preproc:

  Preprocessing function to apply to the data (default: center())

- ncomp:

  Number of latent dimensions for barycenter (default: NULL,
  auto-determine)

- sigma:

  RWR restart probability and diffusion parameter (default: 0.15)

- lambda:

  Overall consistency regularization weight (default: 0.1). Can also
  specify lambda_e, lambda_n, lambda_p separately for fine control

- lambda_e:

  Edge consistency weight (default: lambda \* 0.5)

- lambda_n:

  Neighborhood consistency weight (default: lambda \* 0.5)

- lambda_p:

  Anchor prior weight (default: lambda \* 0.1)

- tau:

  Entropy regularization parameter for Sinkhorn (default: 0.01)

- alpha:

  Weight for feature vs RWR cost in Sylvester equation (default: 0.5)

- gamma:

  Cross-graph mixing parameter for Sylvester iteration (default: 0.1)

- solver:

  Transport solver: "sinkhorn" for entropic (default), "exact" for
  unregularized

- max_iter:

  Maximum number of iterations (default: 100)

- tol:

  Convergence tolerance (default: 1e-6)

- use_cpp:

  Whether to use C++ optimizations (default: FALSE)

- verbose:

  Logical flag to print progress information during fitting

## Value

A `multiblock_biprojector` object containing:

- `s`: Aligned network embeddings

- `v`: Primal vectors for out-of-sample projection

- `alignment_matrix`: Soft alignment/transport plan between networks

- `transport_plan`: Dense transport matrix S

- `sdev`: Standard deviations of aligned components

- Additional metadata for reconstruction and validation

A `multiblock_biprojector` object containing:

- `s`: Aligned network embeddings

- `v`: Primal vectors for out-of-sample projection

- `alignment_matrix`: Soft alignment/transport plan between networks

- `transport_plan`: Dense transport matrix S

- `sdev`: Standard deviations of aligned components

- Additional metadata for reconstruction and validation

A multiblock_biprojector object containing the PARROT alignment results

## Details

PARROT tackles network alignment by formulating it as a regularized
optimal transport problem. The method incorporates position-aware
features through Random Walk with Restart (RWR) descriptors and enforces
structural consistency through neighborhood-preserving regularization
terms.

PARROT operates through the following algorithmic components:

- **Position-Aware Features**: Compute RWR descriptors capturing network
  position

- **Cross-Network Cost**: Build transport cost matrix between networks

- **Consistency Regularization**: Add structural similarity constraints

- **Optimal Transport**: Solve regularized transport problem via
  Sinkhorn

The algorithm minimizes the objective: \$\$L(S) = \langle C, S \rangle +
\lambda_1 \Omega_1(S) + \lambda_2 \Omega_2(S) + \tau H(S)\$\$

where \\C\\ is the position-aware cost matrix, \\\Omega_1, \Omega_2\\
are consistency regularizers, and \\H(S)\\ is the entropy regularization
term with parameter \\\tau\\.

Key parameters:

- `sigma`: RWR restart probability and diffusion parameter

- `lambda`: Consistency regularization weights

- `tau`: Entropy regularization parameter for Sinkhorn

- `solver`: Transport solver ("sinkhorn" or "exact")

PARROT tackles network alignment by formulating it as a regularized
optimal transport problem. The method incorporates position-aware
features through Random Walk with Restart (RWR) descriptors and enforces
structural consistency through neighborhood-preserving regularization
terms.

PARROT operates through the following algorithmic components:

- **Position-Aware Features**: Compute RWR descriptors capturing network
  position

- **Cross-Network Cost**: Build transport cost matrix between networks

- **Consistency Regularization**: Add structural similarity constraints

- **Optimal Transport**: Solve regularized transport problem via
  Sinkhorn

The algorithm minimizes the objective: \$\$L(S) = \langle C, S \rangle +
\lambda_1 \Omega_1(S) + \lambda_2 \Omega_2(S) + \tau H(S)\$\$

where \\C\\ is the position-aware cost matrix, \\\Omega_1, \Omega_2\\
are consistency regularizers, and \\H(S)\\ is the entropy regularization
term with parameter \\\tau\\.

Key parameters:

- `sigma`: RWR restart probability and diffusion parameter

- `lambda`: Consistency regularization weights

- `tau`: Entropy regularization parameter for Sinkhorn

- `solver`: Transport solver ("sinkhorn" or "exact")

## References

Wang, S., Chen, Z., Yu, X., Li, T., Yang, J., & Liu, X. (2022). PARROT:
Position-aware regularized optimal transport for network alignment. In
Proceedings of the 28th ACM SIGKDD Conference on Knowledge Discovery and
Data Mining (pp. 1896-1905).

Wang, S., Chen, Z., Yu, X., Li, T., Yang, J., & Liu, X. (2022). PARROT:
Position-aware regularized optimal transport for network alignment. In
Proceedings of the 28th ACM SIGKDD Conference on Knowledge Discovery and
Data Mining (pp. 1896-1905).

## See also

`parrot.hyperdesign`

`parrot.hyperdesign`

## Examples

``` r
# \donttest{
# Example with hyperdesign network data

# Create synthetic network domains
set.seed(123)
domain1 <- list(
  x = matrix(rnorm(200), 100, 2),  # Node features
  design = data.frame(
    node_id = 1:100,
    anchors = c(1:10, rep(NA, 90))  # First 10 nodes are anchors
  )
)
domain2 <- list(
  x = matrix(rnorm(200), 100, 2),
  design = data.frame(
    node_id = 1:100,
    anchors = c(1:10, rep(NA, 90))  # Corresponding anchors
  )
)

# Create hyperdesign
hd <- structure(list(domain1 = domain1, domain2 = domain2), class = "hyperdesign")

# Run PARROT alignment with default parameters
result <- parrot(hd, anchors = anchors)

# Access alignment results
transport_plan <- result$transport_plan
aligned_embeddings <- result$s

# Use different regularization settings
result_strong <- parrot(hd, anchors = anchors, lambda = 0.5, tau = 0.1)
# }

# \donttest{
# Example with hyperdesign network data
library(multidesign)
library(tibble)

# Create synthetic network domains (node features for PARROT)
set.seed(123)
X1 <- matrix(rnorm(200), 100, 2)  # Node features for domain 1
X2 <- matrix(rnorm(200), 100, 2)  # Node features for domain 2

# Create design data frames with anchor correspondences (PARROT uses anchors)
design1 <- tibble(
  node_id = 1:100,
  anchors = c(1:10, rep(NA, 90))  # First 10 nodes are anchors
)
design2 <- tibble(
  node_id = 1:100,
  anchors = c(1:10, rep(NA, 90))  # Corresponding anchors
)

# Create multidesign objects
md1 <- multidesign(X1, design1)
md2 <- multidesign(X2, design2)

# Create hyperdesign from multidesign objects
hd <- hyperdesign(list(domain1 = md1, domain2 = md2))

# Run PARROT alignment (uses anchors from design component)
result <- parrot(hd, anchors = anchors)

# Access alignment results
transport_plan <- result$transport_plan
aligned_embeddings <- result$s

# Use different regularization settings
result_strong <- parrot(hd, anchors = anchors, lambda = 0.5, tau = 0.1)
# }
```
