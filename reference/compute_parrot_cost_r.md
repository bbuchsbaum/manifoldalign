# Compute Position-Aware Cost Matrix (R implementation)

Follows the two-stage process from the paper: 1. Compute C_node = alpha
\* cost_RWR + (1-alpha) \* cost_attr (Eq. 3) 2. Solve Sylvester equation
for C_rwr (Eq. 4)

## Usage

``` r
compute_parrot_cost_r(
  networks,
  rwr_features,
  anchor_info,
  alpha = 0.5,
  sigma = 0.15,
  gamma = 0.1
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

## Value

Position-aware cost matrix C_rwr
