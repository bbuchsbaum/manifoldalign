# Fused-Partial Gromov-Wasserstein Distance

Computes the Fused-Partial Gromov-Wasserstein (FPGW) distance between
domains in a hyperdesign object. This method extends Gromov-Wasserstein
by combining feature and structural information. Partial transport is
supported through an explicit transported-mass budget, `rho`.

## Usage

``` r
fpgw(data, ...)

# S3 method for class 'hyperdesign'
fpgw(
  data,
  omega1 = 0.001,
  lambda = 0,
  rho = NULL,
  epsilon = 0.01,
  metric = "euclidean",
  max_iter = 200,
  tol = 1e-06,
  inner_max_iter = 50,
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  A hyperdesign object containing multiple data domains

- ...:

  Additional arguments passed to methods:

  - `omega1`: Weight of the feature term (0 ≤ omega1 ≤ 1). The
    structural term weight omega2 = 1 - omega1 is computed automatically
    (default: 0.001)

  - `lambda`: Retained for API compatibility and required to be zero.
    The former `lambda > 0` route is withdrawn because it used a
    balanced transport oracle and therefore could not produce a partial
    coupling. Use `rho` instead (default: 0).

  - `rho`: Mass budget in (0, min(\|μ\|,\|ν\|)\] for the
    mass-constrained variant. If `rho` is supplied, the solver
    transports exactly that mass; otherwise it uses classical balanced
    FGW.

  - `epsilon`: Initial entropic regularization for warm-start (default:
    0.01). Set to 0 to disable warm-start

  - `metric`: Distance metric for within-domain distances (default:
    "euclidean")

  - `max_iter`: Maximum iterations for Frank-Wolfe optimization
    (default: 200)

  - `tol`: Convergence tolerance for the Frank-Wolfe gap (default: 1e-6)

  - `inner_max_iter`: Maximum iterations for inner optimization
    (default: 50)

  - `verbose`: Print convergence information (default: FALSE)

- omega1:

  Weight of the feature term (0 ≤ omega1 ≤ 1). Default: 0.001

- lambda:

  Compatibility parameter that must be zero. Use `rho` for partial
  transport. Default: 0

- rho:

  Mass budget in (0, min(\|μ\|,\|ν\|)\] for mass-constrained variant

- epsilon:

  Initial entropic regularization for warm-start. Default: 0.01

- metric:

  Distance metric for within-domain distances. Default: "euclidean"

- max_iter:

  Maximum iterations for Frank-Wolfe. Default: 200

- tol:

  Convergence tolerance for Frank-Wolfe gap. Default: 1e-6

- inner_max_iter:

  Maximum iterations for inner optimization. Default: 50

- verbose:

  Print convergence information. Default: FALSE

## Value

An fpgw object (inheriting from multiblock_biprojector) containing:

- `transport_plans`: List of optimal transport plans between domains

- `distances`: Matrix of pairwise FPGW distances between domains

- `converged`: Convergence status for each optimization

- `n_samples`: Number of samples per domain

- `omega1`: Feature weight used

- `lambda`: Compatibility field; always zero for new fits

- `rho`: Mass budget (if mass-constrained variant)

- `domain_names`: Names of the domains

- `training_data`: Original domain matrices for prediction

- `metric`: Distance metric used

## Details

The verified fixed-mass FPGW objective combines feature and structural
information through: \$\$L(\gamma) = \omega_1 \langle C, \gamma
\rangle + \omega_2 \sum\_{i,j,i',j'} \|C\_{X,ii'} - C\_{Y,jj'}\|^2
\gamma\_{ij} \gamma\_{i'j'}\$\$

where C is the feature cost matrix and C_X, C_Y are within-domain
distance matrices. With `rho`, the coupling has row and column marginals
bounded by the input masses and satisfies \\\|\gamma\| = \rho\\. Without
`rho`, the coupling has the full balanced marginals.

The implementation supports two verified modes:

- Mass-constrained FPGW: transports exactly `rho` mass

- Classical FGW: transports the full balanced mass

The former penalized mode selected balanced Frank-Wolfe vertices, fixing
the transported mass at one and making `lambda` ineffective as a partial
transport control. It is rejected rather than silently returning a
balanced result. The Frank-Wolfe algorithm uses closed-form step sizes.

## References

Bai et al. (2025). Fused Partial Gromov-Wasserstein for Structured
Objects. arXiv:2502.09934.

## Examples

``` r
# \donttest{
library(multidesign)

# Example 1: Basic FPGW between two domains with different dimensions
set.seed(123)
n <- 30

# Domain 1: 3D data with two clusters
X1 <- matrix(rnorm(n * 3), n, 3)
X1[1:15, ] <- X1[1:15, ] + 2

# Domain 2: 5D data with similar structure
X2 <- matrix(rnorm(n * 5), n, 5)
X2[1:15, ] <- X2[1:15, ] + 2

# Create hyperdesign
design <- data.frame(id = 1:n, cluster = rep(1:2, each = 15))
hd <- hyperdesign(list(
  visual = multidesign(X1, design),
  semantic = multidesign(X2, design)
))

# Classical FGW with balanced feature/structure weight
result <- fpgw(hd, omega1 = 0.5)
print(result)
#> Fused-Partial Gromov-Wasserstein
#> ================================
#> Number of domains: 2 
#> Domain names: visual, semantic 
#> Feature weight (omega1): 0.5 
#> Mode: Classical Fused GW
#> 
#> Pairwise distances:
#>        [,1]   [,2]
#> [1,] 0.0000 0.4898
#> [2,] 0.4898 0.0000
#> 
#> Warning: Some optimizations did not converge
#> Non-converged pairs:
#>   (visual, semantic)

# Example 2: Mass-constrained for noisy data
# Add outliers to domain 2
X2_noisy <- rbind(X2, matrix(rnorm(10 * 5, sd = 5), 10, 5))
design_noisy <- data.frame(id = 1:(n + 10))

hd_noisy <- hyperdesign(list(
  clean = multidesign(X1, design[1:n,]),
  noisy = multidesign(X2_noisy, design_noisy)
))

# Transport only 75% of mass to avoid outliers
result_partial <- fpgw(hd_noisy, omega1 = 0.3, rho = 0.75)

# Check transported mass
P <- result_partial$transport_plans[[1]]
cat("Transported mass:", sum(P), "\n")
#> Transported mass: 0.75 

# Example 3: Multi-domain alignment
X3 <- matrix(rnorm(n * 4), n, 4)
X3[16:30, ] <- X3[16:30, ] + 1.5

hd_multi <- hyperdesign(list(
  modality1 = multidesign(X1, design),
  modality2 = multidesign(X2, design),
  modality3 = multidesign(X3, design)
))

# Compute all pairwise alignments
result_multi <- fpgw(hd_multi, omega1 = 0.2, max_iter = 50)

# Examine distance matrix
print(result_multi$distances)
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.6979227 0.5977135
#> [2,] 0.6979227 0.0000000 0.6298517
#> [3,] 0.5977135 0.6298517 0.0000000
# }
```
