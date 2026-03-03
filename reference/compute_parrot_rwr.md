# Compute RWR Features for PARROT (Dispatcher)

Compute RWR Features for PARROT (Dispatcher)

## Usage

``` r
compute_parrot_rwr(
  networks,
  anchor_info,
  sigma,
  max_iter,
  tol,
  use_cpp = FALSE
)
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

- use_cpp:

  Whether to use C++ implementation

## Value

A list of RWR feature matrices, one per network
