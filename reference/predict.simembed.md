# Predict Method for Similarity Embedding

Projects new data using trained similarity embedding model.

## Usage

``` r
# S3 method for class 'simembed'
predict(object, newdata, ...)
```

## Arguments

- object:

  A simembed object from linear_sim_embed()

- newdata:

  New data matrix (n_new x d) to project

- ...:

  Additional arguments (currently ignored)

## Value

Matrix of projected coordinates (n_new x ncomp)

## Examples

``` r
# \donttest{
X <- matrix(rnorm(100 * 5), 100, 5)
model <- linear_sim_embed(X, ncomp = 2)
X_new <- matrix(rnorm(20 * 5), 20, 5)
Y_new <- predict(model, X_new)
# }
```
