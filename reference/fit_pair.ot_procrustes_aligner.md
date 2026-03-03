# Fit OT-Procrustes on a pair of domains

Estimates an orthogonal map between two domains without known
correspondences by alternating between an entropic OT coupling and a
Procrustes update.

## Usage

``` r
# S3 method for class 'ot_procrustes_aligner'
fit_pair(
  algo,
  X_i,
  X_j,
  links = NULL,
  i = NULL,
  j = NULL,
  a = NULL,
  b = NULL,
  epsilon0 = 1,
  epsilon_min = 0.05,
  decay = 0.95,
  n_init = 1L,
  seed = NULL,
  max_iter = 50,
  tol = 1e-06,
  sinkhorn_max_iter = 200,
  sinkhorn_tol = 1e-06,
  stabilized = (epsilon_min < 0.01),
  cost_scale = c("median", "max", "none"),
  init = c("identity", "signature", "random_orthogonal"),
  force_SO = FALSE,
  store_transport = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- algo:

  An ot_procrustes_aligner object.

- X_i:

  First domain data matrix (samples x features).

- X_j:

  Second domain data matrix (samples x features).

- links:

  Optional correspondence links (currently unused).

- i:

  Internal domain index used by \[align_many()\] (unused).

- j:

  Internal domain index used by \[align_many()\] (unused).

- a:

  Optional source marginal (length nrow(X_i)); defaults to uniform.

- b:

  Optional target marginal (length nrow(X_j)); defaults to uniform.

- epsilon0:

  Initial entropic regularization (\> 0).

- epsilon_min:

  Minimum entropic regularization (\> 0).

- decay:

  Multiplicative decay factor in (0,1\] applied each iteration.

- n_init:

  Number of random restarts. When \`init="identity"\` and \`n_init\>1\`,
  the method tries the identity, then a deterministic distance-signature
  initialization, and uses remaining runs for random orthogonal
  restarts. The best run (lowest loss) is returned.

- seed:

  Optional seed for reproducible initialization. When called via
  \[align_many()\], the internal \`i\` and \`j\` indices are used to
  derive a per-edge seed from this value.

- max_iter:

  Maximum number of outer (alternating) iterations.

- tol:

  Convergence tolerance on max entrywise change in the rotation.

- sinkhorn_max_iter:

  Maximum number of Sinkhorn iterations.

- sinkhorn_tol:

  Sinkhorn convergence tolerance.

- stabilized:

  Logical; use log-domain stabilization in Sinkhorn (slower but helpful
  when \`epsilon_min\` is very small).

- cost_scale:

  Cost normalization used before Sinkhorn: \`"median"\` divides by the
  median positive cost, \`"max"\` divides by the maximum cost, and
  \`"none"\` disables scaling.

- init:

  Initialization method: \`"identity"\`, \`"signature"\`, or
  \`"random_orthogonal"\`.

- force_SO:

  Logical; if TRUE enforce det(R)=+1 for the learned rotation.

- store_transport:

  Logical; if TRUE store the final transport plan.

- verbose:

  Logical; print basic progress information.

- ...:

  Additional arguments (unused).

## Value

An object of class \`ot_procrustes_pair_fit\` containing the learned
rotation (i-\>j), loss, convergence diagnostics, and optionally the
final transport plan.

## Examples

``` r
# \donttest{
set.seed(1)
X <- matrix(rnorm(150), 50, 3)
algo <- ot_procrustes_aligner()
fit <- fit_pair(algo, X, X)
fit$converged
#> [1] FALSE
# }
```
