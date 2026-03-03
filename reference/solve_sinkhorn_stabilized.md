# Log-domain Stabilized Sinkhorn Algorithm (Dispatcher)

Implements numerically stable Sinkhorn in log domain. Routes to either R
or C++ implementation based on use_cpp flag.

## Usage

``` r
solve_sinkhorn_stabilized(C, tau, max_iter, tol, use_cpp = FALSE)
```

## Arguments

- C:

  Cost matrix

- tau:

  Entropy regularization parameter

- max_iter:

  Maximum iterations

- tol:

  Convergence tolerance

- use_cpp:

  Whether to use C++ implementation

## Value

Transport plan matrix
