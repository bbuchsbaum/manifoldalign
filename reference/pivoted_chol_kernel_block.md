# Simple pivoted Cholesky kernel rank selection.

Simple pivoted Cholesky kernel rank selection.

## Usage

``` r
pivoted_chol_kernel_block(
  X,
  kernel,
  tol = 1e-06,
  max_rank = Inf,
  diag_fun = function(X) rep(1, nrow(X))
)
```

## Arguments

- X:

  Data matrix (samples x features).

- kernel:

  kernlab kernel object.

- tol:

  Residual tolerance.

- max_rank:

  Optional rank cap.

- diag_fun:

  Optional shortcut for diagonal entries of kernel matrix.

## Value

List with selected indices and factor matrix W (n x r).

## Examples

``` r
# \donttest{
X <- matrix(rnorm(50), 10, 5)
kernel <- kernlab::rbfdot(sigma = 0.1)
result <- pivoted_chol_kernel_block(X, kernel, tol = 1e-4, max_rank = 5)
str(result)
#> List of 2
#>  $ indices: int [1:5] 1 6 10 8 9
#>  $ W      : num [1:10, 1:5] 1 0.513 0.5 0.788 0.449 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:5] "new_col" "new_col" "new_col" "new_col" ...
# }
```
