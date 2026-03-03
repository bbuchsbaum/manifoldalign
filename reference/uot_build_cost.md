# Build dense or sparse costs for multi-set UOT alignment

Cost builder for the multi-set alignment loop. Supports dense costs for
small parcellations and sparse neighborhood costs (kNN, radius, or
hybrid) for large voxel/vertex problems.

## Usage

``` r
uot_build_cost(
  X,
  Y,
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
  ensure_cols = NULL
)
```

## Arguments

- X:

  Source coordinates (n x 3) or any numeric matrix.

- Y:

  Template coordinates (m x 3) or any numeric matrix.

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
  (1..m) for each source node. When provided with \`lambda_prior \> 0\`,
  adds a soft identity/anatomical bias term based on squared distances
  in template space between \`Y\[j,\]\` and \`Y\[prior_map\[i\],\]\`.

- lambda_prior:

  Nonnegative weight for the soft prior term.

- prior_sigma:

  Optional positive scale (in the units of \`Y\`). When set, the prior
  penalty uses \`\|\|Y\[j\]-Y\[pref\]\|\|^2 / prior_sigma^2\`.

- neighbor_mode:

  Neighborhood mode: \`"auto"\`, \`"dense"\`, \`"knn"\`, \`"radius"\`,
  or \`"hybrid"\`.

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
  template nodes. Must share the same label alphabet as \`group_X\`.

- dense_max_bytes:

  Maximum dense cost size (in bytes) allowed when
  \`neighbor_mode="auto"\`. Above this threshold, \`"auto"\` selects a
  sparse neighborhood.

- ensure_cols:

  Logical; if TRUE (default for sparse kNN/hybrid modes), adds a minimal
  set of reverse 1NN edges to ensure every template column has at least
  one incoming edge in the sparse graph.

## Value

A dense numeric matrix (n x m) when \`neighbor_mode="dense"\` (or when
\`neighbor_mode="auto"\` selects dense); otherwise a sparse cost list
with CSR+CSC fields.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(30), 10, 3)
Y <- matrix(rnorm(36), 12, 3)
C <- uot_build_cost(X, Y, neighbor_mode = "dense")
str(C)
#>  num [1:10, 1:12] 4.87 1.9 8.4 20.25 1.9 ...
```
