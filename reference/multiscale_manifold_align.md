# Multiscale manifold alignment (Wang-Mahadevan style)

Experimental multiset manifold alignment with a QR-free spectral
hierarchy. The implementation builds a joint graph Laplacian, compresses
each domain with truncated SVD, forms a reduced operator, and extracts
nested multiscale bases from its spectrum.

## Usage

``` r
multiscale_manifold_align(data, ...)

# S3 method for class 'hyperdesign'
multiscale_manifold_align(
  data,
  y = NULL,
  correspondences = NULL,
  preproc = multivarious::center(),
  ncomp = 20L,
  control = multiscale_manifold_align_control(),
  verbose = TRUE,
  ...
)

# S3 method for class 'list'
multiscale_manifold_align(
  data,
  y = NULL,
  correspondences = NULL,
  preproc = multivarious::center(),
  ncomp = 20L,
  control = multiscale_manifold_align_control(),
  verbose = TRUE,
  ...
)

# Default S3 method
multiscale_manifold_align(data, ...)
```

## Arguments

- data:

  A \`hyperdesign\` object or list of matrices.

- ...:

  Additional arguments passed to methods.

- y:

  Optional unquoted label column used to build cross-domain
  correspondences.

- correspondences:

  Optional data.frame with columns
  domain_i,index_i,domain_j,index_j\[,weight\].

- preproc:

  Optional preprocessing function/preprocessor applied per domain.

- ncomp:

  Number of target embedding components.

- control:

  Control settings from \[multiscale_manifold_align_control()\].

- verbose:

  Logical progress toggle.

## Value

A \`multiblock_biprojector\`-compatible object with subclass
\`multiscale_manifold_align\`.
