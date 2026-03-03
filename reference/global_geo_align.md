# Global-Geometry Multi-Dataset Alignment

Aligns multiple domains by preserving global graph geometry induced by
within-domain kNN neighborhoods and cross-domain correspondence edges.

## Usage

``` r
global_geo_align(data, ...)

# S3 method for class 'hyperdesign'
global_geo_align(
  data,
  y = NULL,
  correspondences = NULL,
  feature_correspondence = NULL,
  preproc = multivarious::center(),
  ncomp = 20L,
  control = global_geo_align_control(),
  ...
)

# S3 method for class 'list'
global_geo_align(
  data,
  y = NULL,
  correspondences = NULL,
  feature_correspondence = NULL,
  preproc = multivarious::center(),
  ncomp = 20L,
  control = global_geo_align_control(),
  ...
)

# Default S3 method
global_geo_align(data, ...)
```

## Arguments

- data:

  A hyperdesign object or list of matrices.

- ...:

  Additional arguments passed to methods.

- y:

  Optional unquoted column name in design tables for label-based
  correspondences

- correspondences:

  Optional data.frame with columns domain_i/index_i/domain_j/index_j

- feature_correspondence:

  Optional list specifying spatial anchor basis settings

- preproc:

  Preprocessing function or list. Default: multivarious::center()

- ncomp:

  Number of components to extract. Default: 20

- control:

  Control settings from global_geo_align_control()

## Value

A \`multiblock_biprojector\` object with subclass \`global_geo_align\`.

## Examples

``` r
# \donttest{
set.seed(123)
X1 <- matrix(rnorm(50 * 10), 50, 10)
X2 <- matrix(rnorm(60 * 12), 60, 12)
y1 <- rep(1:5, each = 10)
y2 <- rep(1:5, each = 12)
domains <- list(
  domain1 = list(x = X1, design = data.frame(label = y1)),
  domain2 = list(x = X2, design = data.frame(label = y2))
)
class(domains) <- "hyperdesign"
result <- global_geo_align(domains, y = label, ncomp = 3)
#> global_geo_align: preprocessing 2 domains.
#> global_geo_align: running Dijkstra from 110 landmarks.
# }
```
