# Translation-invariant unbalanced OT (KL) via TI-Sinkhorn

Solve entropic unbalanced OT with KL marginal penalties using the
translation-invariant (TI) Sinkhorn updates from Séjourné, Vialard,
Peyré (AISTATS 2022). Supports dense costs and sparse neighbourhood
graphs.

## Usage

``` r
uot_ti_sinkhorn_kl(
  cost,
  alpha,
  beta,
  epsilon,
  rho1,
  rho2,
  max_iter = 2000,
  tol = 1e-06
)
```

## Arguments

- cost:

  A dense numeric matrix (n x m) or a sparse cost list. Sparse costs
  must include CSR fields \`row_ptr\`, \`col_idx\`, \`cost\`,
  \`n_rows\`, \`n_cols\`. If CSC fields \`col_ptr\`, \`row_idx\`,
  \`cost_csc\` are present they are used for faster updates.

- alpha:

  Source masses (length n, nonnegative).

- beta:

  Target masses (length m, nonnegative).

- epsilon:

  Entropic regularization parameter (\> 0).

- rho1:

  KL penalty on the first marginal (\> 0).

- rho2:

  KL penalty on the second marginal (\> 0).

- max_iter:

  Maximum number of TI-Sinkhorn iterations.

- tol:

  Stopping tolerance on iterate change.

## Value

A list with translation-invariant potentials \`fbar\`, \`gbar\`,
translated dual potentials \`f\`, \`g\`, translation \`lambda\`,
\`iterations\`, \`converged\`, and \`residual\`.

## Examples

``` r
set.seed(1)
C <- matrix(runif(30), 5, 6)
a <- rep(1, 5)
b <- rep(1, 6)
fit <- uot_ti_sinkhorn_kl(C, a, b, epsilon = 0.5, rho1 = 1, rho2 = 1,
                          max_iter = 200, tol = 1e-8)
str(fit)
#> List of 8
#>  $ fbar      : num [1:5, 1] -0.319 -0.398 -0.2 -0.316 -0.405
#>  $ gbar      : num [1:6, 1] -0.0128 0.0685 -0.0292 0.127 -0.0443 ...
#>  $ lambda    : num 0.0765
#>  $ f         : num [1:5, 1] -0.242 -0.322 -0.124 -0.24 -0.328
#>  $ g         : num [1:6, 1] -0.08933 -0.00797 -0.10572 0.05055 -0.12077 ...
#>  $ iterations: int 7
#>  $ converged : logi TRUE
#>  $ residual  : num 8.22e-10
```
