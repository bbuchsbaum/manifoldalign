# Compute Position-Aware Cost Matrix (Dispatcher)

Compute Position-Aware Cost Matrix (Dispatcher)

## Usage

``` r
compute_parrot_cost(
  networks,
  rwr_features,
  anchor_info,
  alpha = 0.5,
  sigma = 0.15,
  gamma = 0.1,
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

- alpha:

  Weight for RWR vs attribute cost (default: 0.5)

- sigma:

  RWR restart probability (beta in paper)

- gamma:

  Cross-graph discount factor for Sylvester equation

- use_cpp:

  Whether to use C++ implementation

## Value

A numeric cost matrix between networks
