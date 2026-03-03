# Enhanced ADAM optimizer with alpha scheduling

Enhanced ADAM optimizer with alpha scheduling

## Usage

``` r
.optimize_W_enhanced(
  X,
  T,
  M,
  W0,
  sigma = 1,
  alpha_p = 0.2,
  alpha_schedule = FALSE,
  lr = 0.005,
  batch = 1,
  maxit = 500,
  tol = 1e-06,
  verbose = FALSE
)
```

## Value

List with components W (optimized weight matrix), trace (objective
values), convergence (status code), message (convergence message), and
iterations (count)
