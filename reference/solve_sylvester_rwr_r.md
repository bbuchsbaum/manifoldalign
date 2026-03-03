# Solve Sylvester Equation for Cross-Graph RWR Cost (R implementation)

Implements the Sylvester iteration from Eq. 4 in the paper: C_rwr =
(1+beta) \* C_node + (1-beta) \* gamma \* W1 \* C_rwr \* W2^T

## Usage

``` r
solve_sylvester_rwr_r(
  W1,
  W2T,
  Cnode,
  beta = 0.15,
  gamma = 0.1,
  tol = 1e-06,
  max_iter = 50
)
```

## Arguments

- W1:

  First network's transition matrix (row-stochastic)

- W2T:

  Transpose of second network's transition matrix

- Cnode:

  Node feature cost matrix

- beta:

  RWR restart probability (corresponds to sigma in code)

- gamma:

  Cross-graph discount factor

- tol:

  Convergence tolerance

- max_iter:

  Maximum iterations

## Value

Cross-graph RWR cost matrix C_rwr
