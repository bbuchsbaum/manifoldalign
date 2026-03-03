# Test Out-of-Sample Reconstruction Accuracy

Validates KEMA's out-of-sample reconstruction capability against the L2
error reported in the paper appendix.

## Usage

``` r
validate_out_of_sample_reconstruction(
  expected_error = 0.14,
  tolerance = 0.05,
  test_fraction = 0.2
)
```

## Arguments

- expected_error:

  Expected L2 reconstruction error (default: 0.14)

- tolerance:

  Numerical tolerance for comparison

- test_fraction:

  Fraction of data to use for testing

## Value

List with reconstruction validation results

## Examples

``` r
# \donttest{
result <- validate_out_of_sample_reconstruction()
print(result$success)
#> [1] TRUE
# }
```
