# Solve Sylvester Equation for Cross-Graph RWR Cost (Dispatcher)

Solve Sylvester Equation for Cross-Graph RWR Cost (Dispatcher)

## Usage

``` r
solve_sylvester_rwr(
  W1,
  W2T,
  Cnode,
  beta = 0.15,
  gamma = 0.1,
  tol = 1e-06,
  max_iter = 50,
  use_cpp = FALSE
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

- use_cpp:

  Whether to use C++ implementation

## Value

A numeric matrix solving the Sylvester equation for RWR computation
