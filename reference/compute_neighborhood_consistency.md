# Compute Neighborhood Consistency Regularization

Implements neighborhood consistency via KL divergence between
transported neighborhoods

## Usage

``` r
compute_neighborhood_consistency(networks, transport_plan, lambda_n)
```

## Arguments

- networks:

  List of network structures

- transport_plan:

  Current transport plan S

- lambda_n:

  Neighborhood regularization weight

## Value

Neighborhood consistency cost matrix
