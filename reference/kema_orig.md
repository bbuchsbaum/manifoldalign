# Original Kernel Manifold Alignment (KEMA)

Paper-faithful implementation of Tuia & Camps-Valls (2016) using the
generalized eigenproblems:

## Usage

``` r
kema_orig(data, y, ...)

# S3 method for class 'multidesign'
kema_orig(
  data,
  y,
  subject,
  preproc = center(),
  ncomp = 2,
  knn = 5,
  sigma = 0.73,
  mu = 0.5,
  kernel = coskern(),
  sample_frac = 1,
  lambda = 0,
  backend = "auto",
  backend_control = NULL,
  ...
)

# S3 method for class 'hyperdesign'
kema_orig(
  data,
  y,
  preproc = center(),
  ncomp = 2,
  knn = 5,
  sigma = 0.73,
  mu = 0.5,
  kernel = coskern(),
  sample_frac = 1,
  lambda = 0,
  backend = "auto",
  backend_control = NULL,
  ...
)

# Default S3 method
kema_orig(data, ...)
```

## Arguments

- data:

  Input object. A \`hyperdesign\` or \`multidesign\` object.

- y:

  Label column used for class graphs.

- ...:

  Method-specific arguments.

- subject:

  Subject/domain column for splitting a multidesign object.

- preproc:

  Preprocessing recipe (default \`center()\`).

- ncomp:

  Number of components.

- knn:

  Number of nearest neighbors for within-domain graph.

- sigma:

  Kernel graph scale for local geometry.

- mu:

  Weight for class-pull Laplacian in \\L + \mu L_s\\.

- kernel:

  Kernel function.

- sample_frac:

  Fraction of samples retained as landmarks (\<1 =\> REKEMA).

- lambda:

  Non-negative regularization added to RHS generalized matrix. Set to
  \`0\` for strict paper formulation.

- backend:

  Backend for the generalized eigensolver. One of \`"auto"\`,
  \`"full_exact"\`, \`"reduced_exact"\`, or \`"operator_exact"\`.

- backend_control:

  Optional list for backend auto-selection and fidelity diagnostics.
  Supported keys: \`full_exact_max_n\`, \`reduced_exact_max_r\`,
  \`fidelity_residual_tol\`, \`fidelity_orth_tol\`, \`fidelity_action\`,
  \`primme_tol\`, and \`primme_method\`.

## Value

A \`multiblock_biprojector\` object with KEMA embeddings.

## Details

Full KEMA (Eq. 6): \$\$K (L + \mu L_s) K \Lambda = \lambda K L_d K
\Lambda\$\$

REKEMA (Eq. 10, reduced-rank): \$\$K\_{rn} (L + \mu L_s) K\_{nr} \Lambda
= \lambda K\_{rn} L_d K\_{nr} \Lambda\$\$

This function intentionally omits extension layers (regression solver,
repulsion term \\L_r\\, fail-soft retries, and kernel centering
variants).

## References

Tuia, D., & Camps-Valls, G. (2016). Kernel manifold alignment for domain
adaptation. PLoS ONE, 11(2), e0148655.
