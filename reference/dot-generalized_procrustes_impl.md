# Generalized Orthogonal Procrustes with Partial Tasks

Extends Generalized Orthogonal Procrustes to handle partial task
observations. Uses efficient sparse matrix operations and robust
initialization.

## Usage

``` r
.generalized_procrustes_impl(
  A_list,
  task_labels_list,
  L,
  max_iter = 100,
  tol = 1e-06,
  tol_type = c("relative", "absolute"),
  verbose = FALSE,
  svd_method = c("irlba", "base"),
  svd_opts = list()
)
```

## Arguments

- A_list:

  A list of length \\n\\, where each element \\A_i\\ is a \\d x li\\
  numeric matrix representing the observed data for subject \\i\\. Each
  subject might have different \\li \< L\\.

- task_labels_list:

  A list of length \\n\\, where \\task_labels_list\[\[i\]\]\\ is an
  integer (or numeric) vector of length \\li\\, specifying which global
  tasks each column of \\A_i\\ corresponds to (in \\1..L\\).

- L:

  Integer. The total number of distinct tasks across all subjects.

- max_iter:

  Integer. Maximum number of GPM iterations. Default 100.

- tol:

  Numeric tolerance for convergence. Default `1e-6`.

- tol_type:

  Character. Type of tolerance ("relative" or "absolute"). Note: The
  relative tolerance is calculated as \`sqrt(diff_sq) / norm0\`, where
  \`norm0\` is the Frobenius norm of the initial \`O_mats\` (plus
  epsilon for stability). The absolute tolerance uses \`sqrt(diff_sq) /
  (sqrt(n\*d) + 1)\`. Default "relative".

- verbose:

  Logical. If `TRUE`, prints iteration info. Default `FALSE`.

- svd_method:

  Backend for the initial SVD (one of `"irlba"`, `"base"`). `"irlba"` is
  recommended and used by default, with fallback to random
  initialization if needed. `"base"` forces the use of dense SVD, which
  can be slow/memory-intensive. Default `"irlba"`.

- svd_opts:

  List of extra arguments passed to
  [`irlba::irlba`](https://rdrr.io/pkg/irlba/man/irlba.html) or
  [`base::svd`](https://rdrr.io/r/base/svd.html).

## Value

A list with:

- O_mats:

  List of \\n\\ orthogonal \\d \times d\\ matrices.

- A_est:

  A \\d \times L\\ consensus matrix. Columns corresponding to tasks
  observed by no subjects will contain `NA`.

- iterations:

  Number of iterations performed.

- converged:

  Logical indicating whether the GPM converged.

- final_diff:

  The final difference value (\`sqrt(diff_sq)\`) used for the
  convergence check.

## Details

**Sparse Construction of D**: We form a large sparse matrix \\D \in
R^{(n d) \times L}\\ using
[`Matrix::sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html).
Triplet indices (row, column, value) are collected efficiently.

**Robust Initialization**: Initialization uses the top \\d\\ left
singular vectors of \\D\\, typically computed via
[`irlba::irlba`](https://rdrr.io/pkg/irlba/man/irlba.html). If `irlba`
fails or returns a rank deficient result (less than \\d\\ vectors), it
falls back to initializing with random orthogonal matrices to ensure
robustness.

**Vectorized GPM Iteration**: The core iteration avoids \\O(n^2)\\
complexity by computing the update term \\CO\\ using two efficient
sparse matrix multiplications:

1.  \\Y \<- D^T \* stack(O\_{old})\\

2.  \\CO\_{stack} \<- D \* Y\\

where \\stack(O)\\ reshapes the \\d \times d \times n\\ array of
rotation matrices into an \\(nd) \times d\\ matrix. The resulting \\(nd)
\times d\\ matrix \\CO\_{stack}\\ is then unstacked back into a \\d
\times d \times n\\ array.

**Efficient Projection**: Each \\d \times d\\ slice of the unstacked
update is projected onto the orthogonal group O(d). An efficient
projection method (`.projO`) is used, switching between an
eigenvalue-based method for small \\d\\ and an SVD-based method for
larger \\d\\.

**Convergence Check**: Convergence is assessed by comparing the
Frobenius norm of the difference between consecutive projected iterates
\\O\_{new}\\ and \\O\_{old}\\. For relative tolerance, this difference
is scaled by the Frobenius norm of the initial \\O\\ matrices
(calculated once).

**Final Consensus A_est**: After convergence, the consensus matrix
\\A\_{est}\\ is computed. For each task \\j\\, we average the
transformed data \\O_i^T A_i(\*, j_i)\\ only over subjects \\i\\ in
\\S_j\\ that observed task \\j\\. Task counts (number of subjects
observing each task \\j\\) are computed by directly iterating over
`task_labels_list`, ensuring correctness even if input matrices contain
explicit zeros.

## Examples

``` r
# \donttest{
# Simple example with 2 subjects and 3 tasks
A1 <- matrix(rnorm(12), 3, 4)  # Subject 1: 3 features x 4 tasks
A2 <- matrix(rnorm(9), 3, 3)   # Subject 2: 3 features x 3 tasks

# Task labels (partial overlap)
task_labels1 <- c(1, 2, 3, 4)  # Subject 1 observes tasks 1,2,3,4
task_labels2 <- c(1, 3, 5)     # Subject 2 observes tasks 1,3,5

# Run alignment
result <- generalized_procrustes(
  A_list = list(A1, A2),
  task_labels_list = list(task_labels1, task_labels2),
  L = 5,  # Total of 5 unique tasks
  max_iter = 50,
  verbose = FALSE
)

# Check results
print(result$converged)
#> [1] TRUE
print(dim(result$A_est))  # Should be 3x5
#> [1] 3 5
# }
```
