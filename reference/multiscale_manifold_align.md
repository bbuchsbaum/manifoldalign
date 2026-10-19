# Multiscale manifold alignment (Wang-Mahadevan style)

Deprecated experimental multiset manifold alignment with a QR-free
spectral hierarchy. This implementation is retained temporarily for
reproducibility, but it has unresolved solver/completion failures and is
not recommended. Use \[rapid_ma()\] for supplied sparse relations and
structured attributes, \[ssma_align()\] for semi-supervised feature
manifolds, or \[spectral_mnn_align()\] for an unsupervised spectral
baseline. Removal is scheduled for manifoldalign 1.0.0.

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

## Migration

Given a list \`domains\` of feature matrices, choose the replacement by
input contract: “\`r \# Supplied sparse relations, positions, or
structured attributes rapid_ma(domains, relations = relations, positions
= positions, attributes = attributes)

\# Semi-supervised feature manifolds ssma_align(domains, correspondences
= correspondences)

\# Unsupervised feature-derived spectral alignment
spectral_mnn_align(domains) “\`

## Bounded repair gate

A repair will be considered only with a coherent algorithmic owner and
specification, zero malformed-formula or convergence-code errors on all
mandatory fixtures, at least 95 percent completion on applicable frozen
cases, passing numerical invariants and an independent oracle, and
evidence of a unique capability or non-dominated regime. Otherwise
removal remains scheduled for manifoldalign 1.0.0.
