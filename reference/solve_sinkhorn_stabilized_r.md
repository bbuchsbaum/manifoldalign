# Log-domain Stabilized Sinkhorn Algorithm (R implementation)

Implements numerically stable Sinkhorn in log domain. For a square
matrix (n1=n2), it produces a doubly-stochastic matrix (row/col sums =
1). For a rectangular matrix, it produces a transport plan with uniform
probability marginals.

## Usage

``` r
solve_sinkhorn_stabilized_r(C, tau, max_iter, tol)
```

## Arguments

- C:

  Cost matrix

- tau:

  Entropy regularization parameter

- max_iter:

  Maximum iterations

- tol:

  Convergence tolerance

## Value

Transport plan matrix
