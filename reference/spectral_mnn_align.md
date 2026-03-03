# Spectral MNN Alignment

A fast, correspondence-optional aligner that: 1) builds a within-domain
kNN graph per domain, 2) computes Laplacian eigenmaps, 3) constructs
intrinsic Heat Kernel Signatures (HKS), 4) discovers anchors via mutual
nearest neighbors (or optional supervision), 5) aligns domains to a
reference via Procrustes.

## Usage

``` r
spectral_mnn_align(data, ...)

# S3 method for class 'hyperdesign'
spectral_mnn_align(
  data,
  y = NULL,
  correspondences = NULL,
  preproc = multivarious::center(),
  ncomp = 20L,
  ref_idx = 1L,
  control = spectral_mnn_align_control(),
  ...
)

# S3 method for class 'list'
spectral_mnn_align(data, ...)

# Default S3 method
spectral_mnn_align(data, ...)
```

## Arguments

- data:

  A \`hyperdesign\` object or list of matrices.

- ...:

  Additional arguments passed to methods.

- y:

  Optional unquoted column name in design tables for label-based anchors

- correspondences:

  Optional data.frame with explicit anchor pairs

- preproc:

  Preprocessing function or list. Default: multivarious::center()

- ncomp:

  Number of components to extract. Default: 20

- ref_idx:

  Reference domain index. Default: 1

- control:

  Control settings from spectral_mnn_align_control()

## Value

A \`multiblock_biprojector\` with subclass \`spectral_mnn_align\`.

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(50), 10, 5)
X2 <- matrix(rnorm(50), 10, 5)
data_list <- list(domain1 = X1, domain2 = X2)
result <- spectral_mnn_align(data_list, ncomp = 3)
#> spectral_mnn_align: preprocessing domains
#> spectral_mnn_align: computing spectral bases (k_embed=13)
#> spectral_mnn_align: computing HKS descriptors (q=16)
#> spectral_mnn_align: domain 2 has insufficient anchors (2); using identity rotation
dim(result$scores)
#> NULL
# }
```
