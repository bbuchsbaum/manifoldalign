# Apply UOT coupling as a barycentric map (sparse, vector signal)

Sparse version using a CSR-like row pointer representation.

## Usage

``` r
uot_apply_map_sparse_vec_cpp(
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
  signal,
  delta = 1e-08
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

- signal:

  Source signal vector (length n_rows).

- delta:

  Stabilizer added to the denominator.

## Value

A numeric vector (length n_cols) in template space.
