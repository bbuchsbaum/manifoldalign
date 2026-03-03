# Compute RWR Features for PARROT (R implementation)

Compute RWR Features for PARROT (R implementation)

## Usage

``` r
compute_parrot_rwr_r(networks, anchor_info, sigma, max_iter, tol)
```

## Arguments

- networks:

  List of network structures

- anchor_info:

  List with anchor information for both networks

- sigma:

  RWR restart probability (beta in paper)

- max_iter:

  Maximum iterations

- tol:

  Convergence tolerance

## Value

List of RWR descriptor matrices
