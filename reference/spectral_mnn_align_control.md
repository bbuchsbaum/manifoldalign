# Control settings for \`spectral_mnn_align()\`

Control settings for \`spectral_mnn_align()\`

## Usage

``` r
spectral_mnn_align_control(
  knn = 15L,
  sigma = NULL,
  use_laplacian = TRUE,
  k_embed = NULL,
  q_descriptors = 16L,
  time_mode = c("auto", "fixed"),
  time_range = c(0.1, 50),
  time_scale = 1,
  mnn_k = 10L,
  ratio = 0.9,
  max_anchors = 512L,
  min_anchors = NULL,
  max_pairs_per_label = 100L,
  dist_quantile = 0.3,
  refine_rounds = 1L,
  force_SO = TRUE,
  store_basis = FALSE,
  store_descriptors = FALSE,
  seed = 1L,
  verbose = TRUE
)
```

## Arguments

- knn:

  Integer k for within-domain kNN graph construction.

- sigma:

  Numeric bandwidth for heat-kernel affinities. Use \`NULL\` to
  auto-tune per domain via \[choose_sigma()\].

- use_laplacian:

  Logical; if TRUE (default) use the normalized Laplacian.

- k_embed:

  Optional integer number of eigenpairs used to build heat kernel
  signatures. If \`NULL\`, defaults to \`max(3\*ncomp, ncomp + 10)\`
  inside the solver.

- q_descriptors:

  Integer number of heat kernel signature time steps.

- time_mode:

  Either \`"auto"\` (derive a log-spaced time range from the retained
  eigenvalues) or \`"fixed"\` (use \`time_range\`).

- time_range:

  Length-2 numeric vector giving min/max times when
  \`time_mode="fixed"\`.

- time_scale:

  Positive scalar multiplier applied to the time grid.

- mnn_k:

  Integer k for mutual nearest-neighbor (MNN) discovery in descriptor
  space.

- ratio:

  Ratio test threshold in (0,1\]. Used to discard ambiguous nearest
  neighbors when \`mnn_k \>= 2\`.

- max_anchors:

  Maximum number of anchor pairs used per domain alignment.

- min_anchors:

  Minimum number of anchors required to estimate a stable Procrustes
  rotation. If \`NULL\`, defaults to \`max(10, 3\*ncomp)\` inside the
  solver.

- max_pairs_per_label:

  Maximum correspondences sampled per shared label when \`y\` is
  provided and explicit \`correspondences\` are absent.

- dist_quantile:

  Quantile of (directed) NN distances kept as candidates prior to
  enforcing mutuality and one-to-one pairing. Lower values are more
  selective.

- refine_rounds:

  Integer number of ICP-style refinement rounds. Each round updates
  anchors in the \*aligned\* embedding space and recomputes the
  Procrustes rotation.

- force_SO:

  Logical; if TRUE (default), enforce det(R)=+1.

- store_basis:

  Logical; if TRUE, store eigenpairs in the returned object.

- store_descriptors:

  Logical; if TRUE, store HKS descriptors in the result.

- seed:

  Integer seed for sampling label-based correspondences.

- verbose:

  Logical toggle for progress messages.

## Value

A list of class \`spectral_mnn_align_control\`.

## Examples

``` r
ctrl <- spectral_mnn_align_control(knn = 10, k_embed = 20)
str(ctrl)
#> List of 20
#>  $ knn                : int 10
#>  $ sigma              : NULL
#>  $ use_laplacian      : logi TRUE
#>  $ k_embed            : int 20
#>  $ q_descriptors      : int 16
#>  $ time_mode          : chr "auto"
#>  $ time_range         : num [1:2] 0.1 50
#>  $ time_scale         : num 1
#>  $ mnn_k              : int 10
#>  $ ratio              : num 0.9
#>  $ max_anchors        : int 512
#>  $ min_anchors        : NULL
#>  $ max_pairs_per_label: int 100
#>  $ dist_quantile      : num 0.3
#>  $ refine_rounds      : int 1
#>  $ force_SO           : logi TRUE
#>  $ store_basis        : logi FALSE
#>  $ store_descriptors  : logi FALSE
#>  $ seed               : int 1
#>  $ verbose            : logi TRUE
#>  - attr(*, "class")= chr "spectral_mnn_align_control"
```
