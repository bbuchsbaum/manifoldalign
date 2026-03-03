# Gradient for Single Rotation in Multi-Graph Setting

Computes the gradient of the objective with respect to one rotation
matrix in the context of multiple graphs.

## Usage

``` r
gradient_single_multiset(basis_s, A_list, Ms, s, lambda, Gss = NULL)
```

## Arguments

- basis_s:

  Spectral basis for graph s

- A_list:

  List of descriptor-basis products for all graphs

- Ms:

  List of all rotation matrices

- s:

  Index of current graph

- lambda:

  Regularization parameter

- Gss:

  Optional list of A_s^T A_s matrices for reuse

## Value

Gradient matrix
