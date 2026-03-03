# Predict method for Gromov-Wasserstein

Transport new samples using the learned GW alignment.

## Usage

``` r
# S3 method for class 'gromov_wasserstein'
predict(
  object,
  newdata,
  from,
  to,
  type = c("weights", "transport"),
  k = 5,
  ...
)
```

## Arguments

- object:

  A gromov_wasserstein object

- newdata:

  New data from the source domain (samples x features)

- from:

  Source domain index or name

- to:

  Target domain index or name

- type:

  Prediction type. \`"weights"\` (default) returns barycentric weights
  over target-domain samples. \`"transport"\` returns barycentric
  combinations of target-domain features when available.

- k:

  Number of nearest neighbours used to build the barycentric combination
  (default: 5).

- ...:

  Reserved for future use.

## Value

Matrix of barycentric weights (\`type = "weights"\`) or transported
samples (\`type = "transport"\`).

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(30), 10, 3)
X2 <- matrix(rnorm(30), 10, 3)
library(multidesign)
md1 <- multidesign(X1, data.frame(id = 1:10))
md2 <- multidesign(X2, data.frame(id = 1:10))
hd <- hyperdesign(list(d1 = md1, d2 = md2))
fit <- gromov_wasserstein(hd, epsilon = 0.5, max_iter = 10)
newX <- matrix(rnorm(6), 2, 3)
pred <- predict(fit, newX, from = 1, to = 2)
# }
```
