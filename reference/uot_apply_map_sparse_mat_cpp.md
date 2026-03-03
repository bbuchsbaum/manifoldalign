# Apply UOT coupling as a barycentric map (sparse, matrix signal via CSC)

Sparse matrix version for time series or multi-contrast signals.
Requires a CSC representation of the cost graph (column pointers + row
indices).

## Usage

``` r
uot_apply_map_sparse_mat_cpp(
  col_ptr,
  row_idx,
  cost_csc,
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

- col_ptr:

  Integer vector of length m+1 with CSC column offsets (0-based or
  1-based accepted; if 1-based, it must start at 0 or 1 and be
  nondecreasing).

- row_idx:

  Integer vector of length nnz giving 1-based row indices.

- cost_csc:

  Numeric vector of length nnz aligned with \`row_idx\`.

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

  Source signal matrix (T x n_rows).

- delta:

  Stabilizer added to the denominator.

## Value

A numeric matrix (T x n_cols) in template space.

## Details

This computes \\\hat S = S W / (\pi_2 + \delta)\\ without materializing
the full \\\pi\\. The implementation uses per-column log-stabilization
to avoid overflow in the exponential weights.
