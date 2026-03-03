# Iterative CONE-Align Algorithm

Iterative CONE-Align Algorithm

## Usage

``` r
cone_align_iterate(embeddings, solver, max_iter, tol, lambda, anchors = NULL)
```

## Arguments

- embeddings:

  List of embedding matrices

- solver:

  Assignment solver

- max_iter:

  Maximum iterations

- tol:

  Convergence tolerance

- lambda:

  Regularization parameter

## Value

List with final rotation, permutation, and iteration count
