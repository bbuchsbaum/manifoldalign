# Build a sparse cost graph from neighbour indices/distances (CSR+CSC)

Helper for CPU pipelines: given neighbor indices and distances (e.g.
from \`RANN::nn2\`), compute a sparse bipartite cost graph between
source nodes and a fixed template.

## Usage

``` r
uot_build_sparse_cost_knn_cpp(
  nn_idx,
  nn_dists,
  n_cols,
  lambda_anat = 1,
  F = NULL,
  G = NULL,
  lambda_feat = 0,
  max_dist,
  extra_i = NULL,
  extra_j = NULL,
  extra_dists = NULL,
  Y_template = NULL,
  prior_map = NULL,
  lambda_prior = 0,
  prior_sigma = NA_real_
)
```

## Arguments

- nn_idx:

  Integer matrix (n x k) of 1-based neighbor indices into the template
  support.

- nn_dists:

  Numeric matrix (n x k) of Euclidean distances aligned with \`nn_idx\`.

- n_cols:

  Number of template nodes (m).

- lambda_anat:

  Weight for anatomical distance-squared term.

- F:

  Optional source features (n x D).

- G:

  Optional template features (m x D).

- lambda_feat:

  Weight for feature distance-squared term.

- max_dist:

  Optional maximum distance threshold; edges with \`nn_dists \>
  max_dist\` are dropped.

- extra_i:

  Optional 1-based row indices for additional edges.

- extra_j:

  Optional 1-based column indices for additional edges.

- extra_dists:

  Optional distances for additional edges (aligned to \`extra_i\`,
  \`extra_j\`).

- Y_template:

  Optional template coordinates/features (m x d). When provided together
  with \`prior_map\` and \`lambda_prior \> 0\`, adds a soft
  diagonal/identity bias term based on squared distances in template
  space.

- prior_map:

  Optional integer vector (length n) giving the preferred template index
  (1..m) for each source node. Use NA for “no prior” on a row.

- lambda_prior:

  Nonnegative weight for the soft prior term.

- prior_sigma:

  Optional positive scale (in the units of \`Y_template\`). When finite,
  the prior penalty uses \`\|\|y_j - y_pref\|\|^2 / prior_sigma^2\`.

## Value

A list with CSR fields \`row_ptr\`, \`col_idx\`, \`cost\` and CSC fields
\`col_ptr\`, \`row_idx\`, \`cost_csc\`, plus dimensions.
