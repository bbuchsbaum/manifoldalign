# Extract a sparse coupling / map operator from UOT dual potentials (CSC)

Given translation-invariant dual potentials (fbar,gbar), compute edge
weights on a sparse cost graph in CSC form and return a \`dgCMatrix\`.

## Usage

``` r
uot_extract_coupling_sparse_csc_cpp(
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
  weight_type = 0L,
  delta = 1e-08,
  topk_col = 0L,
  threshold = 0
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

- weight_type:

  See above (0,1,2).

- delta:

  Stabilizer for weight_type=0.

- topk_col:

  Optional top-k pruning per column (0 keeps all).

- threshold:

  Optional threshold pruning (\<= 0 keeps all).

## Value

A list with \`coupling\` (dgCMatrix) and \`log_pi2\` (length n_cols).

## Details

Supported weight types: - weight_type = 0: map operator weights, w_ij /
(pi2_j + delta) - weight_type = 1: conditional weights, w_ij / pi2_j
(column-stochastic) - weight_type = 2: raw coupling weights, w_ij

Pruning is column-wise: - keep only weights \>= threshold (after
weight_type transform) - optionally keep only the topk_col largest
weights per column
