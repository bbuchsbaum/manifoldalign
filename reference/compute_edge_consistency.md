# Compute Sparse Edge Consistency Regularization

Implements sparse version of L_edge for memory efficiency Only computes
distances for existing edges

## Usage

``` r
compute_edge_consistency(networks, transport_plan, lambda_e)
```

## Arguments

- networks:

  List of network structures

- transport_plan:

  Current transport plan S

- lambda_e:

  Edge regularization weight

## Value

Edge consistency cost matrix (sparse computation)
