# Frank-Wolfe solver for FPGW

Frank-Wolfe solver for FPGW

## Usage

``` r
fw_fpgw(
  C,
  Cx,
  Cy,
  p,
  q,
  omega1,
  lambda,
  rho,
  P0,
  max_iter = 200,
  inner_max_iter = 50,
  tol = 1e-06,
  verbose = FALSE
)
```

## Value

A list with transport plan P, distance, convergence status, iterations,
and final gradient.

## Details

`lambda` is retained in this internal signature for compatibility and
must be zero. Supply `rho` for the fixed-mass partial feasible set.
