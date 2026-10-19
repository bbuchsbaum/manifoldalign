# Scalable structure-aware manifold alignment

Scalable structure-aware manifold alignment

## Usage

``` r
rapid_ma(data, ...)

# S3 method for class 'list'
rapid_ma(
  data,
  labels = NULL,
  relations = NULL,
  positions = NULL,
  attributes = NULL,
  ncomp = 16L,
  control = rapid_ma_control(),
  ...
)

# S3 method for class 'hyperdesign'
rapid_ma(
  data,
  y = NULL,
  relations = NULL,
  positions = NULL,
  attributes = NULL,
  ncomp = 16L,
  control = rapid_ma_control(),
  ...
)
```

## Arguments

- data:

  A list of numeric domain matrices or a \`hyperdesign\`.

- ...:

  Method-specific arguments.

- labels:

  Optional list of partially observed label vectors.

- relations:

  Optional list of custom named sparse relations per domain.

- positions:

  Optional list of spatial coordinate matrices per domain.

- attributes:

  Optional list of mixed attribute data.frames per domain.

- ncomp:

  Shared embedding dimension.

- control:

  Settings from \[rapid_ma_control()\].

- y:

  Optional unquoted or string label column in each domain design.

## Value

A \`rapid_ma\` embedding result with sparse couplings and diagnostics.
