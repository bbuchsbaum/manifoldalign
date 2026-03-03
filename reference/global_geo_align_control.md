# Control settings for \`global_geo_align.hyperdesign()\`

Control settings for \`global_geo_align.hyperdesign()\`

## Usage

``` r
global_geo_align_control(
  knn = 15L,
  metric = c("euclidean", "cosine"),
  corr_edge_weight = 0,
  n_landmarks_total = 1000L,
  landmarks_per_dataset = NULL,
  landmark_corr_frac = 0.3,
  k_embed = NULL,
  scale_method = c("none", "pairwise_eta", "robust_median"),
  scale_graph_reg = 1e-04,
  scale_anchor_max = 128L,
  max_pairs_per_label = 100L,
  ridge_eps = 1e-04,
  ridge_relative = TRUE,
  ginv_method = c("auto", "woodbury", "direct"),
  eig_tol = 1e-08,
  embed_block_size = 2000L,
  seed = 1L,
  verbose = TRUE
)
```

## Arguments

- knn:

  Number of within-dataset nearest neighbors used to build each domain
  graph.

- metric:

  Distance metric for within-dataset graph edges (\`"euclidean"\` or
  \`"cosine"\`).

- corr_edge_weight:

  Weight assigned to correspondence edges.

- n_landmarks_total:

  Total number of landmarks across all datasets.

- landmarks_per_dataset:

  Optional fixed landmark count per dataset.

- landmark_corr_frac:

  Fraction of landmarks reserved for correspondence nodes when
  available.

- k_embed:

  Intermediate geometry rank used in landmark MDS. If \`NULL\`, defaults
  to \`max(3 \* ncomp, ncomp + 5)\` inside the solver.

- scale_method:

  Dataset scale harmonization strategy.

- scale_graph_reg:

  Ridge penalty for global log-scale fitting.

- scale_anchor_max:

  Maximum matched anchors per dataset pair used during scale estimation.

- max_pairs_per_label:

  Maximum correspondences sampled per shared label when \`y\` is used to
  generate correspondences.

- ridge_eps:

  Ridge regularization base value for linear map estimation.

- ridge_relative:

  Whether to scale ridge by per-dataset feature energy.

- ginv_method:

  Method for solving the regularized least squares problem.

- eig_tol:

  Minimum eigenvalue threshold retained in spectral steps.

- embed_block_size:

  Column block size used for landmark out-of-sample embedding.

- seed:

  Random seed for landmark/correspondence sampling.

- verbose:

  Logical toggle for progress messages.

## Value

A list of class \`global_geo_align_control\`.

## Examples

``` r
ctrl <- global_geo_align_control(knn = 10, verbose = FALSE)
ctrl$knn
#> [1] 10
```
