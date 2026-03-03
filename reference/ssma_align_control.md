# Control settings for \`ssma_align()\`

Control settings for \`ssma_align()\`

## Usage

``` r
ssma_align_control(
  knn = 15L,
  sigma = NULL,
  corr_edge_weight = 1,
  max_pairs_per_label = 100L,
  mnn_k = 5L,
  max_unsup_pairs_per_pair = 200L,
  rank_per_domain = 128L,
  svd_tol = 1e-06,
  eig_tol = 1e-08,
  regularization = 1e-06,
  solver = c("auto", "reduced", "operator"),
  operator_tol = 1e-06,
  operator_maxiter = 2000L,
  operator_power_iter = 25L,
  seed = 1L,
  verbose = FALSE,
  serial_control = ssma_serial_control()
)
```

## Arguments

- knn:

  Integer k used to build within-domain kNN graphs.

- sigma:

  Optional positive heat-kernel bandwidth for within-domain graph
  weights. If \`NULL\`, per-domain auto-tuning via \[choose_sigma()\] is
  used.

- corr_edge_weight:

  Non-negative scalar correspondence edge weight.

- max_pairs_per_label:

  Integer cap per shared label/domain pair when generating
  correspondences from labels.

- mnn_k:

  Integer k used for unsupervised MNN correspondence discovery.

- max_unsup_pairs_per_pair:

  Integer cap per domain pair for unsupervised correspondences.

- rank_per_domain:

  Integer rank cap for per-domain compression basis.

- svd_tol:

  Relative singular-value tolerance used in compression.

- eig_tol:

  Non-negative threshold for filtering trivial generalized eigenvalues.

- regularization:

  Positive ridge regularizer added to reduced \`B\` matrix before
  whitening.

- solver:

  One of \`"auto"\`, \`"reduced"\`, or \`"operator"\`.

- operator_tol:

  Positive convergence tolerance for operator eigensolve.

- operator_maxiter:

  Maximum ARPACK iterations for operator eigensolve.

- operator_power_iter:

  Number of power-iterations used to estimate operator spectral scale.

- seed:

  Integer random seed used in correspondence sampling.

- verbose:

  Logical progress toggle.

- serial_control:

  Serial-correlation handling settings from \[ssma_serial_control()\].

## Value

A list of class \`ssma_align_control\`.
