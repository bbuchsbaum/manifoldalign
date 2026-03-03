# Solve Procrustes Problem for CONE-Align

Solve Procrustes Problem for CONE-Align

## Usage

``` r
solve_procrustes_cone(Z1, Z2, P, lambda = 0, anchors = NULL)
```

## Arguments

- Z1:

  Source embedding matrix

- Z2:

  Target embedding matrix

- P:

  Permutation matrix (or vector of indices)

- lambda:

  Regularization parameter for Procrustes problem

## Value

Orthogonal transformation matrix Q
