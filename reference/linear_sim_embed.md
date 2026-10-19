# Linear Similarity Embedding using Optimal Transport

Performs linear dimensionality reduction by optimizing a weight matrix W
(d x m) such that the pairwise similarities P computed from the
projected data (\`Y = X match a target similarity matrix T as closely as
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

  Symmetric target similarity matrix (n x n). If NULL, will be computed
  from data.

- M:

  Symmetric non-negative mask matrix (n x n), with positive entries
  selecting or weighting pairs. If NULL, uses all off-diagonal pairs.

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
  A failed compiled fit is retried with exact-objective base-R L-BFGS-B.
  If FALSE, use the R ADAM implementation (default: FALSE).

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
`mask`, `optimizer`, and `convergence`; the latter records normalized
status, backend, objective, evaluation counts when available, and any
fallback source.

## Details

The algorithm follows the Linear Similarity Embedding Framework (SEF)
from Passalis & Tefas (2016). The optimization uses either ADAM (R
implementation) or L-BFGS-B (C++ implementation).

Compiled L-BFGS-B status codes are reported as \`converged\` (0),
\`iteration_limit\` (1), \`optimizer_warning\` (51), \`optimizer_error\`
(52), or \`unknown_failure\`. A nonzero compiled status triggers an
exact-objective base-R L-BFGS-B retry. The function returns a fit only
if one backend converges; otherwise it stops instead of returning
unchecked weights.

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
set.seed(42)
X <- matrix(rnorm(100 * 10), 100, 10)
result <- linear_sim_embed(X, ncomp = 3, verbose = TRUE)
#> Creating target similarity matrix from data distances...
#> Initializing with PCA...
#> Auto-selecting sigma_P using histogram-spread heuristic...
#>   sigma=1.00e-05, spread=0.0000
#>   sigma=3.16e-03, spread=0.0000
#>   sigma=1.00e+00, spread=0.0418
#>   Early stopping: spread decreasing
#> Selected sigma_P = 10.000000
#> it=0001  obj=2.571296e-02  Js=2.8570e-02  Jp=3.0249e-31
#> it=0021  obj=2.115490e-02  Js=2.1070e-02  Jp=2.1915e-02
#> it=0041  obj=2.003754e-02  Js=2.0820e-02  Jp=1.2995e-02
#> it=0061  obj=1.936702e-02  Js=1.9814e-02  Jp=1.5341e-02
#> it=0081  obj=1.891715e-02  Js=1.9358e-02  Jp=1.4949e-02
#> it=0101  obj=1.860067e-02  Js=1.8898e-02  Jp=1.5927e-02
#> it=0121  obj=1.842203e-02  Js=1.8665e-02  Jp=1.6234e-02
#> it=0141  obj=1.833636e-02  Js=1.8543e-02  Jp=1.6477e-02
#> it=0161  obj=1.828937e-02  Js=1.8474e-02  Jp=1.6624e-02
#> it=0181  obj=1.825617e-02  Js=1.8428e-02  Jp=1.6713e-02
#> it=0201  obj=1.822763e-02  Js=1.8389e-02  Jp=1.6772e-02
#> it=0221  obj=1.820082e-02  Js=1.8354e-02  Jp=1.6819e-02
#> it=0241  obj=1.817490e-02  Js=1.8321e-02  Jp=1.6862e-02
#> it=0261  obj=1.815018e-02  Js=1.8289e-02  Jp=1.6902e-02
#> it=0281  obj=1.812798e-02  Js=1.8260e-02  Jp=1.6936e-02
#> Converged at iteration 283 (obj change < 1.000e-06)

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
#>   sigma=1.00e+00, spread=0.2033
#>   Early stopping: spread decreasing
#> Selected sigma_P = 3.162278
# }
```
