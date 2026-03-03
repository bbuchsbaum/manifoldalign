# Compute PARROT objective

Computes the full objective: \\\<C_rwr, S\> + lambda_e \* L_e(S) +
lambda_n \* L_n(S) + lambda_p \* L_p(S)\\.

## Usage

``` r
compute_parrot_objective(
  S,
  C_rwr,
  networks,
  anchor_info,
  lambda_e,
  lambda_n,
  lambda_p
)
```

## Arguments

- S:

  Transport plan matrix

- C_rwr:

  Position-aware cost matrix

- networks:

  List of network structures

- anchor_info:

  Anchor information

- lambda_e:

  Edge consistency weight

- lambda_n:

  Neighborhood consistency weight

- lambda_p:

  Anchor prior weight

## Value

Objective function value
