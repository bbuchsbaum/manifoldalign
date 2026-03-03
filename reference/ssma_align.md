# Semi-Supervised Manifold Alignment (SSMA)

Baseline SSMA implementation based on Wang et al. feature-level manifold
alignment (\`Z^T L Z f = Z^T D Z f\`) solved in a reduced-rank space.

## Usage

``` r
ssma_align(data, ...)

# S3 method for class 'hyperdesign'
ssma_align(
  data,
  y = NULL,
  correspondences = NULL,
  serial_index = NULL,
  preproc = multivarious::center(),
  ncomp = 20L,
  control = ssma_align_control(),
  ...
)

# S3 method for class 'list'
ssma_align(
  data,
  y = NULL,
  correspondences = NULL,
  serial_index = NULL,
  preproc = multivarious::center(),
  ncomp = 20L,
  control = ssma_align_control(),
  ...
)

# Default S3 method
ssma_align(data, ...)
```

## Arguments

- data:

  A \`hyperdesign\` object or list of matrices.

- ...:

  Additional arguments passed to methods.

- y:

  Optional unquoted label column used for label-based correspondences
  when explicit \`correspondences\` are not provided.

- correspondences:

  Optional data.frame with columns
  \`domain_i,index_i,domain_j,index_j\[,weight\]\`.

- serial_index:

  Optional serial/time index specification used for within-domain graph
  decontamination. Can be: 1) \`NULL\`, 2) a single design-column name,
  or 3) a list with one numeric index vector per domain.

- preproc:

  Optional preprocessing function/preprocessor applied per domain.
  Default: \`multivarious::center()\`.

- ncomp:

  Number of aligned components.

- control:

  Control settings from \[ssma_align_control()\].
