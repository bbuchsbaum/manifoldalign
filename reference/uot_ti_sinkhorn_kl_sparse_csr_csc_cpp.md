# KL translation-invariant Sinkhorn (sparse CSR+CSC cost)

Sparse CPU implementation using both CSR and CSC representations. This
improves cache locality and enables per-column parallelism.

## Usage

``` r
uot_ti_sinkhorn_kl_sparse_csr_csc_cpp(
  row_ptr,
  col_idx,
  cost,
  col_ptr,
  row_idx,
  cost_csc,
  n_rows,
  n_cols,
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

- row_ptr:

  Integer vector of length n+1 with CSR row offsets (0-based or 1-based
  accepted; if 1-based, it must start at 0 or 1 and be nondecreasing).

- col_idx:

  Integer vector of length nnz giving 1-based column indices.

- cost:

  Numeric vector of length nnz giving cost values aligned with
  \`col_idx\`.

- col_ptr:

  Integer vector of length m+1 with CSC column offsets (0-based or
  1-based accepted).

- row_idx:

  Integer vector of length nnz giving 1-based row indices.

- cost_csc:

  Numeric vector of length nnz aligned with \`row_idx\`.

- n_rows:

  Number of source nodes (n).

- n_cols:

  Number of target nodes (m).

- alpha:

  Source masses (length n).

- beta:

  Target masses (length m).

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
iteration count, convergence flag, last residual, and backend tag.
