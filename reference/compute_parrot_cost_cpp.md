# Compute PARROT cost matrix

Complete cost matrix computation including Sylvester equation

## Usage

``` r
compute_parrot_cost_cpp(
  X1,
  X2,
  R1,
  R2,
  W1,
  W2,
  alpha = 0.5,
  sigma = 0.15,
  gamma = 0.1
)
```

## Arguments

- X1:

  Features of network 1

- X2:

  Features of network 2

- R1:

  RWR descriptors of network 1

- R2:

  RWR descriptors of network 2

- W1:

  Transition matrix of network 1

- W2:

  Transition matrix of network 2

- alpha:

  Weight for attribute vs RWR cost

- sigma:

  RWR restart probability

- gamma:

  Cross-graph discount factor

## Value

Position-aware cost matrix
