# Control RAPID-MA optimization and bounded sparse components

Control RAPID-MA optimization and bounded sparse components

## Usage

``` r
rapid_ma_control(
  max_iter = 8L,
  min_iter = 1L,
  tol = 1e-04,
  objective_tolerance = 1e-08,
  max_backtracks = 4L,
  backtrack_factor = 0.5,
  map_type = c("ridge", "orthogonal", "identity"),
  map_ridge = 0.05,
  map_step = 0.75,
  prototype_step = 1,
  uot_weight = 1,
  structure_weight = 0.25,
  supervised_weight = 0.5,
  relation = NULL,
  diffusion = NULL,
  prototype = NULL,
  uot = NULL,
  relmatch = NULL,
  seed = 1L,
  verbose = FALSE
)
```

## Arguments

- max_iter:

  Maximum accepted alternating updates.

- min_iter:

  Minimum updates before convergence can be declared.

- tol:

  Relative objective improvement tolerance.

- objective_tolerance:

  Absolute tolerance for monotone acceptance.

- max_backtracks:

  Maximum damped retries after an objective increase.

- backtrack_factor:

  Multiplicative step reduction in \`(0, 1)\`.

- map_type:

  Shared-coordinate map: \`"ridge"\`, \`"orthogonal"\`, or
  \`"identity"\`.

- map_ridge:

  Positive ridge penalty for small latent map solves.

- map_step:

  Initial map update fraction in \`(0, 1\]\`.

- prototype_step:

  Initial prototype update fraction in \`(0, 1\]\`.

- uot_weight, structure_weight, supervised_weight:

  Nonnegative objective component weights.

- relation, diffusion, prototype, uot, relmatch:

  Named nested control lists.

- seed:

  Nonnegative deterministic seed.

- verbose:

  Emit iteration progress.

## Value

A \`rapid_ma_control\` object.
