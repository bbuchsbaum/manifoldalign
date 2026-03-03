# Fit Token-OT Graph alignment on a pair of domains

Fit Token-OT Graph alignment on a pair of domains

## Usage

``` r
# S3 method for class 'token_ot_graph_aligner'
fit_pair(
  algo,
  X_i,
  X_j,
  links = NULL,
  ncomp = 20L,
  control = token_ot_graph_align_control(),
  ...
)
```

## Arguments

- algo:

  A token_ot_graph_aligner object.

- X_i:

  First domain data matrix (nodes x features).

- X_j:

  Second domain data matrix (nodes x features).

- links:

  Optional anchor links. Accepts either list(vec1=, vec2=) or list(v1,
  v2), where shared values indicate anchors and NAs denote unknowns.

- ncomp:

  Latent dimension used for baseline embeddings (default: 20).

- control:

  Control settings from token_ot_graph_align_control().

- ...:

  Additional arguments (unused).

## Value

An object of class token_ot_graph_pair_fit with fields transport (sparse
coupling), assignment, objective, and dimensions.
