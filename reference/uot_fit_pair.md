# Fit a pairwise UOT alignment to a fixed template support

Convenience wrapper that builds a dense/sparse cost graph via
\[uot_build_cost()\] and then solves entropic KL-UOT via
\[uot_ti_sinkhorn_kl()\].

This is the recommended entry point for \*\*pairwise\*\* workflows; it
bundles the pieces you need for mapping and coupling extraction into a
single object.

## Usage

``` r
uot_fit_pair(
  X,
  Y,
  alpha,
  beta,
  F = NULL,
  G = NULL,
  lambda_anat = 1,
  lambda_feat = 0,
  prior_map = NULL,
  lambda_prior = 0,
  prior_sigma = NULL,
  neighbor_mode = c("auto", "dense", "knn", "radius", "hybrid"),
  k_neighbors = NULL,
  radius = NULL,
  maxk = 128L,
  min_neighbors = 1L,
  group_X = NULL,
  group_Y = NULL,
  dense_max_bytes = 2.56e+08,
  ensure_cols = NULL,
  epsilon,
  rho1,
  rho2,
  max_iter = 2000,
  tol = 1e-06
)
```

## Arguments

- X:

  Source coordinates (n x 3) or any numeric matrix.

- Y:

  Template coordinates (m x 3) or any numeric matrix.

- alpha:

  Source masses (length n, nonnegative, not all zero).

- beta:

  Target masses (length m, nonnegative, not all zero).

- F:

  Optional source features (n x D).

- G:

  Optional template features (m x D).

- lambda_anat:

  Weight for anatomical squared distance.

- lambda_feat:

  Weight for feature squared distance.

- prior_map:

  Optional integer vector (length n) giving a preferred target index
  (1..m) for each source node. Passed to \[uot_build_cost()\].

- lambda_prior:

  Nonnegative weight for the soft prior term.

- prior_sigma:

  Optional positive scale for the prior term.

- neighbor_mode:

  Neighborhood mode passed to \[uot_build_cost()\].

- k_neighbors:

  Integer k for kNN neighborhoods (used by \`"knn"\` and \`"hybrid"\`).

- radius:

  Positive scalar radius for \`"radius"\` / \`"hybrid"\` modes.

- maxk:

  Maximum number of candidate neighbors requested in \`"radius"\` /
  \`"hybrid"\` mode (prevents quadratic blowups).

- min_neighbors:

  Minimum number of within-radius neighbors required per row in
  \`"hybrid"\` mode before falling back to kNN.

- group_X:

  Optional vector (length n) of hard constraint group labels for source
  nodes. When provided with \`group_Y\`, neighborhoods are computed
  within matching groups only.

- group_Y:

  Optional vector (length m) of hard constraint group labels for
  template nodes.

- dense_max_bytes:

  Maximum dense cost size (in bytes) allowed when
  \`neighbor_mode="auto"\`.

- ensure_cols:

  Logical; if TRUE, adds reverse 1NN edges to ensure every template
  column has at least one incoming edge (kNN/hybrid modes only).

- epsilon:

  Entropic regularization parameter (\> 0).

- rho1:

  KL penalty on the first marginal (\> 0).

- rho2:

  KL penalty on the second marginal (\> 0).

- max_iter:

  Maximum number of TI-Sinkhorn iterations.

- tol:

  TI-Sinkhorn tolerance.

## Value

An object of class \`uot_pair_fit\` with fields \`cost\`, \`sol\`,
\`alpha\`, \`beta\`, \`epsilon\`, \`rho1\`, \`rho2\`.
