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

For TV penalty (lambda \> 0), the objective includes the smooth
quadratic term: lambda \* (\|mu\|^2 + \|nu\|^2 - 2\|gamma\|^2) where
\|gamma\| = sum(gamma) is the total transported mass. This is NOT an L1
penalty but a smooth penalty that encourages partial transport. The
Frank-Wolfe oracle remains standard OT, with the penalty handled in step
size.
