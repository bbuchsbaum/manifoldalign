# Solve PARROT Transport Problem with Proximal Point Method

Implements the proximal point loop that iteratively refines the cost
matrix and solves the optimal transport problem. This is critical for
convergence.

## Usage

``` r
solve_parrot_transport(
  networks,
  rwr_features,
  anchor_info,
  lambda_e,
  lambda_n,
  lambda_p,
  tau,
  alpha,
  sigma,
  gamma,
  solver,
  max_iter,
  max_outer = 10,
  tol,
  use_cpp = FALSE
)
```

## Arguments

- networks:

  List of network structures

- rwr_features:

  List of RWR descriptor matrices

- anchor_info:

  List with anchor information for both networks

- lambda_e:

  Edge consistency weight

- lambda_n:

  Neighborhood consistency weight

- lambda_p:

  Anchor prior weight

- tau:

  Entropy regularization parameter

- solver:

  Transport solver method

- max_iter:

  Maximum iterations for inner Sinkhorn

- max_outer:

  Maximum outer proximal iterations

- tol:

  Convergence tolerance

## Value

List with transport plan and related matrices
