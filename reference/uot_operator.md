# Extract a reusable map/coupling operator from a \`uot_pair_fit\`

Extract a reusable map/coupling operator from a \`uot_pair_fit\`

## Usage

``` r
uot_operator(
  fit,
  weights = c("map", "cond", "coupling"),
  delta = 1e-08,
  prune_topk_col = NULL,
  prune_threshold = NULL
)
```

## Arguments

- fit:

  An object returned by \[uot_fit_pair()\].

- weights:

  Passed to \[uot_extract_coupling()\].

- delta:

  Stabilizer for \`"map"\` weights.

- prune_topk_col:

  Optional top-k pruning per template column.

- prune_threshold:

  Optional threshold pruning.

## Value

A \`dgCMatrix\` of dimension \`n x m\`.
