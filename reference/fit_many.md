# Native multi-domain entry point (optional)

Native multi-domain entry point (optional)

## Usage

``` r
fit_many(algo, domains, ...)
```

## Arguments

- algo:

  An aligner object that supports multi-domain fitting

- domains:

  A hyperdesign object or list of domain data

- ...:

  Additional arguments passed to methods

## Value

The return value depends on the aligner method. Typically a
multiblock_biprojector or similar alignment result object.

## Examples

``` r
# \donttest{
set.seed(1)
library(multidesign)
X1 <- matrix(rnorm(60), 30, 2)
X2 <- matrix(rnorm(60), 30, 2)
labels <- sample(c("A", "B"), 30, TRUE)
design1 <- data.frame(labels = labels)
design2 <- data.frame(labels = labels)
md1 <- multidesign(X1, design1)
md2 <- multidesign(X2, design2)
hd <- hyperdesign(list(d1 = md1, d2 = md2))
algo <- kema_aligner()
result <- fit_many(algo, hd, y = labels, ncomp = 2)
#> Semi-supervised KEMA: 60 labeled samples, 0 unlabeled samples
#> Auto-selected sigma = 1.1419 using median distance heuristic
#> Using RBF kernel with auto-tuned sigma
#> normalize_laplacian(): 60 isolated nodes detected - treating as disconnected components.
#> Warning: Regression solver produced poor results (subspace angle: 35.4 deg, best match: 0.759). Automatically retrying with solver='exact' for higher fidelity.
#> Retry with exact solver completed successfully.
# }
```
