# Predict method for FPGW

Transport new samples using learned FPGW alignment

## Usage

``` r
# S3 method for class 'fpgw'
predict(
  object,
  newdata,
  from,
  to,
  type = c("transport", "weights", "embedding"),
  k = NULL,
  ...
)
```

## Arguments

- object:

  An fpgw object

- newdata:

  New data from one domain to transport

- from:

  Source domain index or name

- to:

  Target domain index or name

- type:

  Type of prediction: "transport" (default), "weights", or "embedding"

- k:

  Number of nearest neighbors used for barycentric interpolation

- ...:

  Additional arguments

## Value

Transported samples, barycentric weights, or embeddings (when
implemented)

## Details

This method provides out-of-sample extension for FPGW by:

- Finding nearest neighbors in the source domain

- Using the transport plan to map to target domain

- Interpolating based on transport weights

The object must contain the training data used for fitting, otherwise an
error is thrown.

## Examples

``` r
if (FALSE) { # \dontrun{
# After computing FPGW alignment
result <- fpgw(hd, omega1 = 0.1)

# Transport new samples from domain 1 to domain 2
new_samples <- matrix(rnorm(10 * 5), 10, 5)
transported <- predict(result, new_samples, from = 1, to = 2)
head(transported)
} # }
```
