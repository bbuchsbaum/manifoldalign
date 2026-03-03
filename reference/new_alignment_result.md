# Construct an alignment result with consistent S3 metadata

Construct an alignment result with consistent S3 metadata

## Usage

``` r
new_alignment_result(
  scores,
  loadings,
  preproc = NULL,
  feature_blocks,
  subclass,
  extras = list()
)
```

## Arguments

- scores:

  Matrix of concatenated scores (samples x components)

- loadings:

  Matrix of loadings/primal vectors

- preproc:

  Pre-processing object or NULL

- feature_blocks:

  Feature block indices (list or matrix)

- subclass:

  Character vector of additional classes to prepend

- extras:

  Named list of extra slots to attach to the result

## Value

Object inheriting from multiblock_biprojector and subclass
