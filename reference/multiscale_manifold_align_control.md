# Control settings for multiscale manifold alignment backends

Control settings for multiscale manifold alignment backends

## Usage

``` r
multiscale_manifold_align_control(
  enabled = FALSE,
  backend = c("spectral", "randomized_dwt"),
  operator_mode = c("compressed_graph_core"),
  rank_per_domain = 128L,
  svd_tol = 1e-06,
  eigen_tol = 1e-08,
  max_levels = 8L,
  eps_rank = 0.001,
  eigen_solver = c("lobpcg", "lanczos"),
  eigen_k_max = 512L,
  knn = 12L,
  sigma = NULL,
  cross_edge_weight = 1,
  mnn_k = 5L,
  max_pairs_per_label = 100L,
  max_unsup_pairs_per_pair = 200L,
  randomized_sketch_size = 256L,
  randomized_oversample = 16L,
  randomized_power_iters = 1L,
  shift_invert_delta = 0.001,
  cg_tol = 1e-06,
  cg_max_iter = 500L,
  scale_selection = c("energy", "supervised", "hybrid"),
  supervision_weight = 0.5,
  min_scale_weight = 0.05,
  max_basis_per_scale = 32L,
  regularization = 1e-06,
  tune = FALSE,
  candidate_cross_edge_weight = NULL,
  candidate_rank_per_domain = NULL,
  candidate_max_levels = NULL,
  seed = 1L,
  verbose = FALSE
)
```

## Arguments

- enabled:

  Logical; enable multiscale mode.

- backend:

  Character backend name. One of \`"spectral"\` or \`"randomized_dwt"\`.

- operator_mode:

  Character operator mode. Currently \`"compressed_graph_core"\`.

- rank_per_domain:

  Integer per-domain truncation rank for compressed operators.

- svd_tol:

  Numeric positive tolerance for truncated SVD decisions.

- eigen_tol:

  Numeric positive tolerance for non-trivial eigenvalues.

- max_levels:

  Integer maximum hierarchy depth.

- eps_rank:

  Numeric relative truncation threshold in \`(0, 1\]\`.

- eigen_solver:

  Character eigensolver selector. One of \`"lobpcg"\` or \`"lanczos"\`.

- eigen_k_max:

  Integer cap on eigenvectors for spectral hierarchy.

- knn:

  Integer nearest-neighbor graph degree for within-domain graphs.

- sigma:

  Optional positive scalar graph bandwidth. If \`NULL\`, auto-tuned per
  domain.

- cross_edge_weight:

  Numeric non-negative weight for cross-domain correspondence edges.

- mnn_k:

  Integer nearest-neighbor count for unsupervised MNN correspondences.

- max_pairs_per_label:

  Integer cap per label/pair when building label correspondences.

- max_unsup_pairs_per_pair:

  Integer cap per domain pair in unsupervised MNN mode.

- randomized_sketch_size:

  Integer sketch size for randomized backend.

- randomized_oversample:

  Integer oversampling budget for randomized backend.

- randomized_power_iters:

  Integer number of subspace power iterations.

- shift_invert_delta:

  Numeric positive regularization for shift-invert style operators.

- cg_tol:

  Numeric positive tolerance for iterative solves.

- cg_max_iter:

  Integer maximum iterations for iterative solves.

- scale_selection:

  How to weight scales when assembling the multiscale dictionary. One of
  \`"energy"\`, \`"supervised"\`, or \`"hybrid"\`.

- supervision_weight:

  Weight in \`\[0, 1\]\` used when \`scale_selection = "hybrid"\`.

- min_scale_weight:

  Relative threshold in \`\[0, 1)\` for keeping a scale in the final
  dictionary.

- max_basis_per_scale:

  Integer cap on the number of basis vectors retained from each scale
  after compression.

- regularization:

  Positive ridge regularization for the final dictionary generalized
  eigenproblem.

- tune:

  Logical; whether to tune \`cross_edge_weight\`, \`rank_per_domain\`,
  and \`max_levels\` over small candidate grids before the final fit.

- candidate_cross_edge_weight:

  Optional numeric vector of candidate cross-edge weights used when
  \`tune = TRUE\`.

- candidate_rank_per_domain:

  Optional integer vector of candidate \`rank_per_domain\` values used
  when \`tune = TRUE\`.

- candidate_max_levels:

  Optional integer vector of candidate \`max_levels\` values used when
  \`tune = TRUE\`.

- seed:

  Integer random seed.

- verbose:

  Logical; print backend progress messages.

## Value

A list with class \`multiscale_manifold_align_control\`.
