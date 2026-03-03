# Out-of-sample projection/prediction

Out-of-sample projection/prediction

## Usage

``` r
oos_predict(fit_or_transform, newX, side, ...)
```

## Arguments

- fit_or_transform:

  A fit object or transform object

- newX:

  New data matrix to project (samples x features)

- side:

  Which domain side the new data belongs to ("i" or "j")

- ...:

  Additional arguments passed to methods

## Value

A matrix of predicted/projected coordinates for the new data in the
aligned space.

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(50), 25, 2)
X2 <- matrix(rnorm(50), 25, 2)
newX <- matrix(rnorm(20), 10, 2)
# Create a simple 2-block projector (loadings stacked by feature block)
v1 <- prcomp(X1)$rotation[, 1, drop = FALSE]
v2 <- prcomp(X2)$rotation[, 1, drop = FALSE]
v <- rbind(v1, v2)
s <- rbind(X1 %*% v1, X2 %*% v2)
fit <- structure(
  list(
    v = v,
    s = s,
    sdev = apply(s, 2, sd),
    preproc = NULL,
    block_indices = list(
      i = seq_len(ncol(X1)),
      j = (ncol(X1) + 1L):(ncol(X1) + ncol(X2))
    )
  ),
  class = "multiblock_biprojector"
)
pred <- oos_predict(fit, newX, side = "i")
# }
```
