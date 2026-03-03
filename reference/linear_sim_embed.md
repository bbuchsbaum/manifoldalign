# Linear Similarity Embedding using Optimal Transport

Performs linear dimensionality reduction by optimizing a weight matrix W
(d x m) such that the pairwise similarities P computed from the
projected data (Y = X match a target similarity matrix T as closely as
possible, considering a mask M. The optimization balances similarity
preservation (Js) with an orthogonality constraint (Jp) on W, minimizing
(1-alpha_p)\*Js + alpha_p\*Jp.

## Usage

``` r
linear_sim_embed(
  X,
  T = NULL,
  M = NULL,
  sigma_P = "auto",
  ncomp = 2,
  alpha_p = 0.1,
  alpha_schedule = FALSE,
  maxit = 500,
  tol = 1e-06,
  batch_size = 1,
  use_cpp = FALSE,
  verbose = FALSE,
  lr = 0.005,
  formula = NULL,
  data = NULL,
  ...
)
```

## Arguments

- X:

  Input data matrix (n x d), where n is samples, d is features.

- T:

  Target similarity matrix (n x n). If NULL, will be computed from data.

- M:

  Mask matrix (n x n), 1 to consider pair (i,j), 0 to ignore. If NULL,
  uses all pairs.

- sigma_P:

  Numeric scalar \> 0 or "auto". Scale parameter for Gaussian kernel. If
  "auto", uses log-space grid search with histogram-spread heuristic
  (default: "auto").

- ncomp:

  Integer \> 0. Number of dimensions for the embedding (m) (default: 2).

- alpha_p:

  Numeric in \[0,1\]. Weight for the orthogonality regularizer. Uses
  convex combination: (1-alpha_p)\*Js + alpha_p\*Jp (default: 0.1).

- alpha_schedule:

  Logical. If TRUE, linearly decay alpha_p from 1 to specified value
  over first 50 iterations to avoid early orthogonality trapping
  (default: FALSE).

- maxit:

  Integer. Maximum number of iterations (default: 500).

- tol:

  Numeric. Convergence tolerance. For ADAM (R): change in objective
  function. For L-BFGS-B (C++): gradient norm tolerance (default: 1e-6).

- batch_size:

  Numeric in (0, 1\]. Fraction of data for stochastic updates in ADAM (R
  only) (default: 1, i.e., full batch).

- use_cpp:

  Logical. If TRUE and C++ backend is available, use L-BFGS-B from C++.
  Otherwise, use the R ADAM implementation (default: FALSE).

- verbose:

  Logical. Print optimization progress (default: FALSE).

- lr:

  Numeric \> 0. Learning rate for the ADAM optimizer (R only) (default:
  5e-3).

- formula:

  Optional formula interface for supervised targets (e.g., ~ label).

- data:

  Optional data.frame when using formula interface.

- ...:

  Extra arguments (currently ignored).

## Value

A `simembed` object (S3 class) containing: - `weights` (W): The
optimized projection matrix (d x m). - `scores` (Y): The projected data
(n x m). - `sdev`: Standard deviations of the scores. - `preproc`: The
preprocessing object used on X. - `center`: Centering vector used in
preprocessing. - `scale`: Scaling vector used in preprocessing. -
`sigma_P`: Final sigma_P value used (important if auto-selected). -
`alpha_p`: Final alpha_p value used. - `objective_trace`: Vector of
objective function values during optimization. - Metadata: `target_sim`,
`mask`, `optimizer`, `convergence`.

## Details

The algorithm follows the Linear Similarity Embedding Framework (SEF)
from Passalis & Tefas (2016). The optimization uses either ADAM (R
implementation) or L-BFGS-B (C++ implementation).

Similarities are computed using a Gaussian kernel: P_ij = exp(-\|\|Y_i -
Y_j\|\|^2 / sigma_P). The orthogonality penalty is Jp = \|\|W'W -
I\|\|^2_F / (2\*m^2).

Key algorithmic features:

- Automatic sigma_P selection using histogram-spread heuristic (Step 2,
  Fig. 2)

- PCA/KPCA initialization for faster convergence (Step 3, Fig. 2)

- Enforced similarity matrix symmetry for numerical stability

- Alpha_p scheduling option for improved convergence

## References

Passalis, N., & Tefas, A. (2016). Learning deep representations with
probabilistic knowledge transfer. In Proceedings of the European
Conference on Computer Vision (pp. 268-284).

## See also

[`predict.simembed`](https://bbuchsbaum.github.io/manifoldalign/reference/predict.simembed.md)

## Examples

``` r
# \donttest{
# Basic usage with automatic sigma_P selection
X <- matrix(rnorm(100 * 10), 100, 10)
result <- linear_sim_embed(X, ncomp = 3, verbose = TRUE)
#> Creating target similarity matrix from data distances...
#> Initializing with PCA...
#> Auto-selecting sigma_P using histogram-spread heuristic...
#>   sigma=1.00e-05, spread=0.0000
#>   sigma=3.16e-03, spread=0.0000
#>   sigma=1.00e+00, spread=0.0369
#>   Early stopping: spread decreasing
#> Selected sigma_P = 10.000000
#> it=0001  obj=2.400115e-02  Js=2.6668e-02  Jp=1.0750e-31
#> Warning: Optimization stalled at iteration 16 (objective unchanged for 5 iterations). This often indicates poorly scaled input. Check similarity matrix scaling.

# Use predict method for new data
X_new <- matrix(rnorm(20 * 10), 20, 10)
Y_new <- predict(result, X_new)

# Custom target similarity matrix
D_orig <- as.matrix(dist(X))
T_sim <- exp(-D_orig^2 / median(D_orig^2))
result2 <- linear_sim_embed(X, T = T_sim, sigma_P = 1.0)

# Formula interface for supervised embedding
df <- data.frame(X, label = sample(c("A", "B", "C"), 100, replace = TRUE))
result3 <- linear_sim_embed(~ label, data = df, ncomp = 2)

# Use C++ backend for speed
result_cpp <- linear_sim_embed(X, use_cpp = TRUE, verbose = TRUE)
#> Creating target similarity matrix from data distances...
#> Initializing with PCA...
#> Auto-selecting sigma_P using histogram-spread heuristic...
#>   sigma=1.00e-05, spread=0.0000
#>   sigma=3.16e-03, spread=0.0000
#>   sigma=1.00e+00, spread=0.2060
#>   Early stopping: spread decreasing
#> Selected sigma_P = 3.162278
#> Warning: When called from R, the RNG seed has to be set at the R level via set.seed()
#> Warning: C++ optimization (L-BFGS-B) did not converge. Message: ERROR: ABNORMAL_TERMINATION_IN_LNSRCH
# }
```
