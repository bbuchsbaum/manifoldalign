# Generalized Orthogonal Procrustes Alignment

Performs Generalized Orthogonal Procrustes alignment to find orthogonal
transformations. Aligns data from multiple domains by minimizing squared
differences between corresponding task observations.

## Usage

``` r
generalized_procrustes(data, ...)

# S3 method for class 'hyperdesign'
generalized_procrustes(
  data,
  y,
  preproc = center(),
  max_iter = 100,
  tol = 1e-06,
  tol_type = c("relative", "absolute"),
  verbose = FALSE,
  svd_method = c("irlba", "base"),
  svd_opts = list(),
  ...
)
```

## Arguments

- data:

  A hyperdesign object containing multiple data domains

- ...:

  Additional arguments (currently unused)

- y:

  Name of the task/label variable to use for alignment

- preproc:

  Preprocessing function to apply to the data (default: `center()`)

- max_iter:

  Maximum number of iterations (default: 100)

- tol:

  Convergence tolerance (default: 1e-6)

- tol_type:

  Type of tolerance check (default: "relative")

- verbose:

  Whether to print progress messages (default: FALSE)

- svd_method:

  SVD method to use (default: "irlba")

- svd_opts:

  Options for SVD method (default: empty list)

## Value

The return value depends on the specific method:

- For hyperdesign objects: A list containing orthogonal transformation
  matrices, consensus matrix, convergence information, and domain
  metadata

- For direct matrix input: A list with transformation matrices and
  alignment results

## Details

This method extends classical Procrustes analysis to handle partial task
observations, where each domain may observe only a subset of a global
set of tasks. The algorithm uses an efficient Generalized Power Method
(GPM) with sparse matrix operations and robust initialization to find
optimal orthogonal transformations.

The Generalized Procrustes problem seeks to find orthogonal matrices
\\O_i\\ for each domain \\i\\ that minimize: \$\$\sum\_{i,j} \sum\_{k
\in T\_{ij}} \|\|O_i^T A_i(:,k) - O_j^T A_j(:,k)\|\|^2\$\$

where \\T\_{ij}\\ represents the set of tasks observed by both domains
\\i\\ and \\j\\.

The algorithm handles several key challenges:

- **Partial observations**: Each domain may observe different subsets of
  tasks

- **Sparse structure**: Uses efficient sparse matrix operations for
  scalability

- **Robust initialization**: SVD-based initialization with fallback to
  random orthogonal matrices

- **Convergence guarantees**: Monotonic improvement with configurable
  tolerance

Key features:

- Vectorized operations using sparse matrix algebra

- Efficient projection onto the orthogonal group O(d)

- Optional tightness certificate for global optimality validation

- Flexible preprocessing and multiple input formats

## References

Gower, J. C. (1975). Generalized procrustes analysis. Psychometrika,
40(1), 33-51.

Ten Berge, J. M. F. (1977). Orthogonal procrustes rotation for two or
more matrices. Psychometrika, 42(2), 267-276.

## See also

`generalized_procrustes.hyperdesign`

## Examples

``` r
# \donttest{
# Example with hyperdesign data
library(multidesign)

# Create example domains with partial task overlap
d1_data <- matrix(rnorm(50), 5, 10)  # 5 tasks x 10 features
d1_design <- data.frame(task = factor(c("A", "B", "C", "D", "E")))
d1 <- multidesign(d1_data, d1_design)

d2_data <- matrix(rnorm(40), 4, 10)  # 4 tasks x 10 features
d2_design <- data.frame(task = factor(c("A", "C", "D", "F")))
d2 <- multidesign(d2_data, d2_design)

# Create hyperdesign
hd <- hyperdesign(list(domain1 = d1, domain2 = d2))

# Perform alignment
result <- generalized_procrustes(hd, task)
#> Warning: irlba::irlba failed: max(nu, nv) must be strictly less than min(nrow(A), ncol(A)). Falling back to random orthogonal matrices.

# Check convergence and results
print(result$converged)
#> [1] FALSE
print(dim(result$A_est))  # Features x total tasks
#> [1] 10  6
# }

# \donttest{
# Create example hyperdesign data
library(multidesign)

# Domain 1: 5 tasks x 10 features
d1_data <- matrix(rnorm(50), 5, 10)
d1_design <- data.frame(task = factor(c("A", "B", "C", "D", "E")))
d1 <- multidesign(d1_data, d1_design)

# Domain 2: 4 tasks x 10 features (partial overlap)
d2_data <- matrix(rnorm(40), 4, 10) 
d2_design <- data.frame(task = factor(c("A", "C", "D", "F")))
d2 <- multidesign(d2_data, d2_design)

# Create hyperdesign
hd <- hyperdesign(list(domain1 = d1, domain2 = d2))

# Perform alignment
result <- generalized_procrustes(hd, task)
#> Warning: irlba::irlba failed: max(nu, nv) must be strictly less than min(nrow(A), ncol(A)). Falling back to random orthogonal matrices.

# Access results
print(result$converged)
#> [1] FALSE
print(dim(result$A_est))  # 10 features x 6 total tasks
#> [1] 10  6
# }
```
