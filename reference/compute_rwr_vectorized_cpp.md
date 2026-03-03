# Vectorized RWR computation

Compute Random Walk with Restart for multiple restart vectors
simultaneously

## Usage

``` r
compute_rwr_vectorized_cpp(W_transpose, E, sigma, max_iter, tol)
```

## Arguments

- W_transpose:

  Transpose of transition matrix

- E:

  Matrix of restart vectors (columns are different restart
  distributions)

- sigma:

  Restart probability

- max_iter:

  Maximum iterations

- tol:

  Convergence tolerance

## Value

RWR result matrix
