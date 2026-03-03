# Compute anchor prior gradient

Gradient of KL divergence from anchor prior

## Usage

``` r
compute_anchor_gradient_cpp(
  S_prev,
  anchor_idx1,
  anchor_idx2,
  anchor_vals1,
  anchor_vals2,
  eps = 1e-16
)
```

## Arguments

- S_prev:

  Previous transport plan

- anchor_idx1:

  Indices of anchors in network 1

- anchor_idx2:

  Indices of anchors in network 2

- anchor_vals1:

  Anchor values for network 1

- anchor_vals2:

  Anchor values for network 2

- eps:

  Small epsilon for numerical stability

## Value

Gradient matrix
