# Calculate Gradient of Similarity Term (Js) w.r.t W (Correct Implementation)

Calculate Gradient of Similarity Term (Js) w.r.t W (Correct
Implementation)

## Usage

``` r
grad_Js(X, W, T, M, sigma, sum_M_batch)
```

## Arguments

- X:

  Preprocessed data matrix (n x d) or batch (n_b x d)

- W:

  Current weight matrix (d x m)

- T:

  Target similarity matrix (n x n) or batch (n_b x n_b)

- M:

  Mask matrix (n x n) or batch (n_b x n_b)

- sigma:

  Scaling parameter (sigma_P)

- sum_M_batch:

  The sum of the \*batch\* mask matrix M (for consistent normalization)

## Value

Gradient matrix (d x m)
