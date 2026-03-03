# Solve Sylvester equation for RWR cost

Iterative solution of: C_rwr = (1+β)C_node + (1-β)γ \* W1 \* C_rwr \*
W2^T

## Usage

``` r
solve_sylvester_rwr_cpp(
  W1,
  W2T,
  C_node,
  beta = 0.15,
  gamma = 0.1,
  tol = 1e-06,
  max_iter = 50L
)
```

## Arguments

- W1:

  Transition matrix of network 1

- W2T:

  Transpose of transition matrix of network 2

- C_node:

  Node-level cost matrix

- beta:

  RWR restart probability

- gamma:

  Cross-graph discount factor

- tol:

  Convergence tolerance

- max_iter:

  Maximum iterations

## Value

Solution matrix C_rwr
