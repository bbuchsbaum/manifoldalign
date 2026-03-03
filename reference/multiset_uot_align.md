# Multi-subject UOT alignment to a shared template

Iteratively aligns multiple subject measures to a fixed-support template
by solving subject-to-template UOT problems and updating template
weights (and optionally template features) from transported marginals.

## Usage

``` r
multiset_uot_align(
  datasets,
  template,
  omega = NULL,
  epsilon,
  rho1,
  rho2,
  lambda_anat = 1,
  lambda_feat = 0,
  neighbor_mode = c("auto", "dense", "knn", "radius", "hybrid"),
  k_neighbors = NULL,
  radius = NULL,
  maxk = 128L,
  min_neighbors = 1L,
  constraint_fields = NULL,
  prior_fields = NULL,
  lambda_prior = 0,
  prior_sigma = NULL,
  dense_max_bytes = 2.56e+08,
  ensure_cols = NULL,
  max_outer = 10,
  max_inner = 2000,
  tol_inner = 1e-06,
  tol_outer = 1e-04,
  rescale_pi2 = c("alpha", "unit"),
  target_mass = "mean_alpha",
  learn_template_features = FALSE,
  delta = 1e-08,
  parallel = FALSE,
  ncores = 1,
  verbose = FALSE
)
```

## Arguments

- datasets:

  A list of datasets. Each element must be a list with fields \`X\` (n x
  3 coordinates), \`alpha\` (length n masses), and optional \`F\` (n x D
  features).

- template:

  A list with fields \`Y\` (m x 3 coordinates), \`beta\` (length m
  masses), and optional \`G\` (m x D features).

- omega:

  Optional nonnegative weights over datasets (length K). If NULL, uses
  uniform weights.

- epsilon:

  Entropic regularization parameter (\> 0).

- rho1:

  KL penalty on the first marginal (\> 0).

- rho2:

  KL penalty on the second marginal (\> 0).

- lambda_anat:

  Weight for anatomical squared distance.

- lambda_feat:

  Weight for feature squared distance.

- neighbor_mode:

  Neighborhood mode passed to \[uot_build_cost()\].

- k_neighbors:

  Integer k for kNN neighborhoods (used by \`"knn"\` and \`"hybrid"\`
  modes).

- radius:

  Positive scalar radius for \`"radius"\` / \`"hybrid"\` modes.

- maxk:

  Maximum number of candidate neighbors requested in \`"radius"\` /
  \`"hybrid"\` mode (prevents quadratic blowups).

- min_neighbors:

  Minimum number of within-radius neighbors required per row in
  \`"hybrid"\` mode before falling back to kNN.

- constraint_fields:

  Optional character vector of field names present on each dataset and
  the template that define \*\*hard transport constraints\*\*. Edges are
  only allowed when all fields match (implemented via grouped neighbor
  search). Typical examples: \`c("hemi")\` or \`c("hemi","network")\`.

- prior_fields:

  Optional character vector of field names present on each dataset and
  the template that define a \*\*soft identity prior\*\*. When provided
  together with \`lambda_prior \> 0\`, each dataset builds a
  \`prior_map\` by matching these fields to the template and adds a soft
  bias term via \[uot_build_cost()\].

- lambda_prior:

  Nonnegative weight for the soft prior term.

- prior_sigma:

  Optional positive scale for the prior term (in the units of
  \`template\$Y\`).

- dense_max_bytes:

  Maximum dense cost size (in bytes) allowed when
  \`neighbor_mode="auto"\`.

- ensure_cols:

  Logical; if TRUE, adds reverse 1NN edges to ensure every template
  column has at least one incoming edge (kNN/hybrid modes only).

- max_outer:

  Maximum number of template update iterations.

- max_inner:

  Maximum number of TI-Sinkhorn iterations per subject.

- tol_inner:

  TI-Sinkhorn tolerance.

- tol_outer:

  Template weight tolerance.

- rescale_pi2:

  How to scale each subject's transported mass before aggregating
  (\`"alpha"\` or \`"unit"\`).

- target_mass:

  Target total mass for template weights. If \`"mean_alpha"\`, uses mean
  subject mass; if \`"none"\`, no rescaling; if numeric, uses that.

- learn_template_features:

  Logical; update template features \`G\`.

- delta:

  Small stabilizer for feature update denominators.

- parallel:

  Logical; if TRUE, uses \`parallel::mclapply\` on Unix.

- ncores:

  Number of cores for parallel subject solves.

- verbose:

  Logical; print outer-loop progress.

## Value

A list with updated \`template\` and per-subject \`fits\`.

## Examples

``` r
# \donttest{
set.seed(1)
K <- 3
datasets <- lapply(seq_len(K), function(k) {
  X <- matrix(rnorm(30), 10, 3)
  F <- matrix(rnorm(40), 10, 4)
  list(X = X, alpha = rep(1, 10), F = F)
})
template <- list(Y = matrix(rnorm(36), 12, 3),
                 beta = rep(1, 12),
                 G = matrix(0, 12, 4))
fit <- multiset_uot_align(datasets, template,
  epsilon = 0.5, rho1 = 1, rho2 = 1,
  lambda_anat = 1, lambda_feat = 0.1,
  k_neighbors = 6, max_outer = 3, max_inner = 200
)
str(fit$template$beta)
#>  num [1:12] 1.589 1.806 1.307 0.171 0.21 ...
# }
```
