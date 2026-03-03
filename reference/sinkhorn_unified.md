# Unified Sinkhorn algorithm with optional log-domain stabilization

Unified Sinkhorn algorithm with optional log-domain stabilization

## Usage

``` r
sinkhorn_unified(
  cost,
  a,
  b,
  epsilon,
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

- max_iter:

  Maximum iterations

- tol:

  Convergence tolerance

- stabilized:

  Use log-domain stabilization for numerical stability

## Value

Transport plan matrix
