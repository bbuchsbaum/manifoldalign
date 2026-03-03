# Unified Sinkhorn algorithm returning dual potentials for warm-starting

This variant returns both the transport plan and the (log) scaling
vectors that can be reused to warm-start subsequent Sinkhorn solves.

## Usage

``` r
sinkhorn_unified_potentials(
  cost,
  a,
  b,
  epsilon,
  log_u0 = NULL,
  log_v0 = NULL,
  max_iter = 1000L,
  tol = 1e-09,
  stabilized = TRUE
)
```

## Arguments

- cost:

  Cost matrix (n x m)

- a:

  Source marginal (n x 1)

- b:

  Target marginal (m x 1)

- epsilon:

  Entropic regularization parameter

- log_u0:

  Optional warm-start log scaling for rows (length n)

- log_v0:

  Optional warm-start log scaling for columns (length m)

- max_iter:

  Maximum iterations

- tol:

  Convergence tolerance

- stabilized:

  Use log-domain stabilization for small epsilon

## Value

A list with elements pi (plan), log_u, log_v
