# Compute log second marginal and transported feature means (sparse CSR)

Same as \`uot_kl_logpi2_meanF_dense_cpp\` but for a sparse cost graph.

## Usage

``` r
uot_kl_logpi2_meanF_sparse_cpp(
  row_ptr,
  col_idx,
  cost,
  n_rows,
  n_cols,
  alpha,
  beta,
  fbar,
  gbar,
  epsilon,
  F = NULL
)
```

## Arguments

- row_ptr:

  Integer vector of length n+1 with CSR row offsets (0-based or 1-based
  accepted; if 1-based, it must start at 0 or 1 and be nondecreasing).

- col_idx:

  Integer vector of length nnz giving 1-based column indices.

- cost:

  Numeric vector of length nnz aligned with \`col_idx\`.

- n_rows:

  Number of source nodes.

- n_cols:

  Number of target nodes.

- alpha:

  Source masses (length n_rows).

- beta:

  Target masses (length n_cols).

- fbar:

  Translation-invariant source potential (length n_rows).

- gbar:

  Translation-invariant target potential (length n_cols).

- epsilon:

  Entropic regularization parameter (\> 0).

- F:

  Optional source features (n_rows x D).

## Value

A list with \`log_pi2\` and optionally \`mean_F\`.
