# KL translation-invariant Sinkhorn (dense cost)

CPU implementation of the KL translation-invariant Sinkhorn updates for
entropic unbalanced OT with KL marginal penalties.

## Usage

``` r
uot_ti_sinkhorn_kl_dense_cpp(
  cost,
  alpha,
  beta,
  epsilon,
  rho1,
  rho2,
  max_iter = 2000L,
  tol = 1e-06
)
```

## Arguments

- cost:

  Dense cost matrix (n x m).

- alpha:

  Source masses (length n, nonnegative, not all zero).

- beta:

  Target masses (length m, nonnegative, not all zero).

- epsilon:

  Entropic regularization parameter (\> 0).

- rho1:

  KL penalty on the first marginal (\> 0).

- rho2:

  KL penalty on the second marginal (\> 0).

- max_iter:

  Maximum number of iterations.

- tol:

  Stopping tolerance on the infinity-norm iterate difference.

## Value

A list with translation-invariant potentials \`fbar\`, \`gbar\`,
translated dual potentials \`f\`, \`g\`, translation \`lambda\`,
iteration count, convergence flag, and last residual.
