# Cross-Validated Row-Wise Alignment Scoring

Runs explicit row-index cross-validation for multi-domain alignment
models by delegating fold construction and execution to
[`multidesign::cv_rows()`](https://rdrr.io/pkg/multidesign/man/cv_rows.html)
and
[`multidesign::cross_validate()`](https://rdrr.io/pkg/multidesign/man/cross_validate.html).
Held-out rows are scored by how similar their cross-block neighbours are
in an external feature space, using either latent
[`oos_predict`](https://bbuchsbaum.github.io/manifoldalign/reference/oos_predict.md)
projections or method-specific `predict(..., type = "weights")` support
when available.

## Usage

``` r
cv_alignment_rows(
  data,
  rows,
  fit_fn,
  features,
  k = 5L,
  feature_similarity = c("cosine", "correlation"),
  target_pool = c("analysis", "assessment", "both"),
  prediction_mode = c("auto", "embedding", "weights")
)
```

## Arguments

- data:

  A `hyperdesign`-compatible object. If `data` is not already a
  multidesign hyperdesign, it is coerced internally.

- rows:

  Explicit held-out row specification forwarded to
  [`multidesign::cv_rows()`](https://rdrr.io/pkg/multidesign/man/cv_rows.html).
  For hyperdesigns, each fold should be a named list mapping block names
  or positions to held-out row indices.

- fit_fn:

  Function taking the analysis split for a fold and returning a fitted
  alignment object.

- features:

  External feature matrices aligned to the original rows of each block.
  Supply either a named list of matrices, one per block, or a
  hyperdesign-compatible object whose `$x` matrices are treated as
  features.

- k:

  Positive integer number of latent nearest neighbours used when scoring
  held-out rows.

- feature_similarity:

  Similarity used for the external feature space. Either `"cosine"`
  (default) or `"correlation"`.

- target_pool:

  Which rows are available as retrieval targets for each held-out query
  block: `"analysis"` (default) uses training rows from the fitted fold,
  `"assessment"` uses only other held-out rows, and `"both"`
  concatenates the two pools when latent projection is available.

- prediction_mode:

  How neighbour rankings are produced. `"auto"` (default) uses latent
  [`oos_predict`](https://bbuchsbaum.github.io/manifoldalign/reference/oos_predict.md)
  projection for embedding-oriented fits and prefers
  `predict(..., type = "weights")` for transport-style fits when
  `target_pool = "analysis"`. `"embedding"` forces latent neighbour
  search and `"weights"` forces weight-based ranking against training
  targets.

## Value

A `multidesign::cv_result` with one row per fold. The score table
includes `mean_top1_similarity`, `mean_topk_similarity`,
`oracle_top1_similarity`, `oracle_topk_similarity`, `top1_gap`,
`topk_gap`, `n_queries`, and `n_pairs`.

## Details

The external feature matrices are used only for held-out evaluation. If
your fitting procedure uses correspondence tables keyed by original
rows, include a stable row identifier in each block's design before
calling this function so the training correspondences can be rebuilt
inside `fit_fn`.
