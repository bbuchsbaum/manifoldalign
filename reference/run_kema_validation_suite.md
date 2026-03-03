# Comprehensive KEMA Numerical Validation

Runs all numerical validation tests for KEMA implementation.

## Usage

``` r
run_kema_validation_suite(verbose = TRUE)
```

## Arguments

- verbose:

  Whether to print detailed results

## Value

List with all validation results

## Examples

``` r
# \donttest{
results <- run_kema_validation_suite(verbose = FALSE)
#> Warning: Eigenvalue extraction failed: upper value must be greater than lower value
print(results$overall_success)
#> [1] FALSE
# }
```
