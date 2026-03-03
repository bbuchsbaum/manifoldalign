# Control settings for \`gpca_align.hyperdesign()\`

Control settings for \`gpca_align.hyperdesign()\`

## Usage

``` r
gpca_align_control(
  knn = NA_integer_,
  knn_mode = c("symmetric", "mutual"),
  balance = c("none", "within", "all"),
  balance_power = 1,
  normalize = c("fro", "edges", "degree"),
  verbose = TRUE,
  max_dense_elems = 2e+08,
  memory_fraction = 0.5,
  memory_safety_bytes = 512 * 1024^2
)
```

## Arguments

- knn:

  Optional integer k for k-NN sparsification (applied to within and
  between graphs)

- knn_mode:

  Whether to use a "symmetric" (union) or "mutual" (intersection) k-NN
  graph

- balance:

  Whether to apply label balancing to "none", "within", or "all" edges

- balance_power:

  Exponent for inverse-frequency label scaling (\>= 0)

- normalize:

  Normalisation mode for the component metrics: "fro", "edges", or
  "degree"

- verbose:

  Logical toggle for informational messages

- max_dense_elems:

  Maximum dense elements permitted when forming the block matrix (set
  \`NA\` to auto-tune)

- memory_fraction:

  Fraction of available RAM to consider safe when auto-tuning
  densification

- memory_safety_bytes:

  Absolute byte cushion to subtract when auto-tuning densification

## Value

A list of class \`gpca_align_control\`

## Examples

``` r
ctrl <- gpca_align_control(knn = 15, normalize = "fro")
```
