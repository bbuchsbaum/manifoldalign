# Apply a UOT coupling to map signals into template space

Applies the entropic UOT coupling implied by translation-invariant dual
potentials, without materializing the full transport plan.

## Usage

``` r
uot_apply_map(cost, alpha, beta, fbar, gbar, epsilon, signal, delta = 1e-08)
```

## Arguments

- cost:

  A dense numeric matrix (n x m) or a sparse cost list with CSR fields
  from
  [`uot_build_cost`](https://bbuchsbaum.github.io/manifoldalign/reference/uot_build_cost.md).

- alpha:

  Source masses (length n).

- beta:

  Target masses (length m).

- fbar:

  Translation-invariant source potential (length n).

- gbar:

  Translation-invariant target potential (length m).

- epsilon:

  Entropic regularization parameter (\> 0).

- signal:

  Source signal: numeric vector (length n) or matrix (T x n).

- delta:

  Stabilizer added to the denominator.

## Value

A numeric vector (length m) or matrix (T x m).

## Examples

``` r
set.seed(1)
C <- matrix(runif(30), 5, 6)
a <- rep(1, 5)
b <- rep(1, 6)
fit <- uot_ti_sinkhorn_kl(C, a, b, epsilon = 0.5, rho1 = 1, rho2 = 1,
                          max_iter = 200, tol = 1e-8)
s <- rnorm(5)
y <- uot_apply_map(C, a, b, fit$fbar, fit$gbar, epsilon = 0.5, signal = s)
str(y)
#>  num [1:6] 0.384 0.593 0.338 0.456 0.534 ...
```
