# Extract a sparse coupling or mapping operator (dgCMatrix)

Given translation-invariant dual potentials from
\[uot_ti_sinkhorn_kl()\], this function computes edge weights on a
sparse cost graph and returns a \`Matrix::dgCMatrix\`. This is useful
for downstream pipelines that want to: (1) prune the coupling and/or (2)
reuse the same operator across multiple signals without recomputing
exponentials.

## Usage

``` r
uot_extract_coupling(
  cost,
  alpha,
  beta,
  fbar,
  gbar,
  epsilon,
  weights = c("map", "cond", "coupling"),
  delta = 1e-08,
  prune_topk_col = NULL,
  prune_threshold = NULL
)
```

## Arguments

- cost:

  A sparse cost list as returned by \[uot_build_cost()\] containing CSC
  fields \`col_ptr\`, \`row_idx\`, \`cost_csc\`, plus dimensions
  \`n_rows\`, \`n_cols\`.

- alpha:

  Source masses (length \`n_rows\`, nonnegative, not all zero).

- beta:

  Target masses (length \`n_cols\`, nonnegative, not all zero).

- fbar:

  Translation-invariant source potential (length \`n_rows\`).

- gbar:

  Translation-invariant target potential (length \`n_cols\`).

- epsilon:

  Entropic regularization parameter (\> 0).

- weights:

  Which weights to compute: - \`"map"\`: barycentric map weights
  \\w\_{ij} / (\pi\_{2,j} + \delta)\\ - \`"cond"\`: conditional weights
  \\w\_{ij} / \pi\_{2,j}\\ (column-stochastic) - \`"coupling"\`: raw
  coupling weights \\w\_{ij}\\ (may overflow for extreme parameter
  regimes).

- delta:

  Stabilizer for \`"map"\` weights.

- prune_topk_col:

  Optional integer top-k pruning per template column. \`NULL\` keeps all
  edges.

- prune_threshold:

  Optional nonnegative threshold pruning applied after the weight
  transform. \`NULL\` keeps all edges.

## Value

A list with: - \`coupling\`: a \`dgCMatrix\` of dimension \`n_rows x
n_cols\` - \`log_pi2\`: numeric vector (length \`n_cols\`) with
\`log(pi2)\` on the sparse graph
